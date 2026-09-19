#! /usr/bin/env nextflow

def rna_ref_map = [:]
if (params.rna_ref_species){
    new File(params.rna_ref_species).eachLine { line ->
        def (species, ref) = line.split('\t')
        rna_ref_map[species] = ref
    }
}

// Load source-unit metadata: complete FASTQ basename -> STAR read-group fields.
// Collision-safe FASTQ names and this table are produced upstream.
def rg_meta = [:]
if (params.rg_metadata){
    def first = true
    new File(params.rg_metadata).eachLine { line ->
        if (first) { first = false; return }  // skip header
        def fields = line.split('\t')
        // fields: fnbase, library, bp_id, sample_idx, lane, flowcell, lane_num, pu, rg_id, rg_string
        if (fields.length >= 10){
            // For STAR, store the individual RG fields (not the full @RG\t string)
            // STAR --outSAMattrRGline wants: ID:xxx SM:xxx PL:xxx PU:xxx DS:xxx
            rg_meta[fields[0]] = "ID:${fields[8]} SM:${fields[1]} PL:Illumina PU:${fields[7]} DS:${fields[2]}"
        }
    }
}

/**
 * Strip optional __RunTag from a filename component before extracting
 * the library name.
 */
def strip_run_tag(s){
    return s.replaceFirst(/__Run\d+$/, '')
}

/**
 * Canonicalize both legacy names (..._19_S3_L001) and bcl-convert names
 * carrying pre-sample annotations (..._19_L1_L2_S3_L001).
 */
def canonical_library_name(s){
    def before_sample = strip_run_tag(s).replaceFirst(/_S\d+_L\d+$/, '')
    return before_sample.replaceFirst(/(?:_L\d+)+$/, '')
}

/**
 * Size RNA mapping memory from the compressed, trimmed FASTQs.
 *
 * The upper-envelope model was calibrated from 36 completed TET libraries:
 *     peak GiB <= 4 + 0.224 * trimmed FASTQ GiB
 *
 * The initial allocation adds 25% operational headroom and rounds upward to
 * an 8 GiB scheduling bucket. params.memgb remains an explicit floor, and an
 * OOM retry adds 32 GiB without changing the underlying mapping inputs.
 */
def rna_memory_request(reads1, reads2, minimum_gb, attempt){
    def bytes_per_gib = 1024d * 1024d * 1024d
    def trimmed_gib = (reads1 + reads2).collect{ read -> read.size() }.sum() / bytes_per_gib
    def upper_envelope_gb = 4d + (0.224d * trimmed_gib)
    def calibrated_gb = Math.ceil((1.25d * upper_envelope_gb) / 8d) * 8d
    def initial_gb = Math.max(minimum_gb.toString().toDouble(), calibrated_gb)
    def requested_gb = initial_gb + (32d * (attempt.toInteger() - 1))
    return "${requested_gb.toInteger()} GB"
}

def rna_geometry = params.rna_geometry ?: 'long-r2'
if (!(rna_geometry in ['long-r2', 'pe150'])){
    error("rna_geometry must be long-r2 or pe150; received: " + rna_geometry)
}

process map_rna{
    cpus params.threads
    time { 120.hour * task.attempt }
    memory { rna_memory_request(reads1, reads2, params.memgb, task.attempt) }
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input: 
    tuple val(lib),
        file(reads1),
        file(reads2),
        val(basenames),
        val(refname),
        file(ref),
        file(whitelist)

    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: { path -> path }

    output:
    tuple val(lib),
        path("gex.bam"),
        path("gex.bam.bai"),
        path("Barcodes.stats"),
        path("Features.stats"),
        path("Summary.csv"),
        path("UMIperCellSorted.txt"),
        path("CellReads.stats.gz"),
        path("STAR_Log.out"),
        path("STAR_Log.final.out"),
        path("STAR_SJ.out.tab.gz"),
        path("raw/*"),
        path("filtered/*"),
        path("*Unmapped.out.mate*", optional: true)

    script:
    def reftrunc = refname.split('/')[-1]    
    def r1 = reads1.join(',')
    def r2 = reads2.join(',')
    def lib2 = lib.replace('/', '_')
    def sortmem = (params.memgb.toInteger() - 1) * 1024 * 1024 * 1024
    def missing_rg = basenames.findAll{ bn -> !rg_meta.containsKey(bn) }
    if (missing_rg){
        error("Missing RNA read-group metadata for input unit(s): " + missing_rg.join(', '))
    }
    def rgline = basenames.collect{ bn -> rg_meta[bn] }.join(' , ')
    // Standard 10X uses a separate, final barcode read: cDNA R2 then barcode R1.
    // 5' PE150 embeds CB16+UMI12 in mate 1. Upstream trimming removes the internal
    // TSO, so STAR clips only the retained 28 bp barcode block and
    // aligns the downstream R1 cDNA together with R2.
    def readfiles = rna_geometry == 'pe150' ? "${r1} ${r2}" : "${r2} ${r1}"
    def geometry_args = (rna_geometry == 'pe150'
        ? '--soloBarcodeMate 1 --clip5pNbases 28 0'
        : '--soloBarcodeMate 0')
    """
    if [ ! -d ${reftrunc} ]; then
        mkdir ${reftrunc}
        mv ${ref} ${reftrunc}
    fi
    if [ \$( file -L --mime-type -b ${whitelist} | grep "gzip" | wc -l ) -gt 0 ]; then 
        zcat ${whitelist} > wl_unzip.txt
    else 
        cp ${whitelist} wl_unzip.txt
    fi

    STAR --genomeDir ${reftrunc} \
     --runThreadN ${params.threads} \
     --readFilesIn ${readfiles} \
     --readFilesCommand zcat \
     --clipAdapterType CellRanger4 \
     --outBAMsortingThreadN 1 \
     --limitBAMsortRAM ${sortmem} \
     --outFileNamePrefix ${lib2} \
     --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ${rgline} \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 \
     --soloCBlen 16 \
     --soloUMIstart 17 \
     --soloUMIlen 12 \
     ${geometry_args} \
     --soloCBwhitelist wl_unzip.txt \
     --outFilterScoreMinOverLread 0.33 \
     --outFilterMatchNminOverLread 0.33 \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
     --soloUMIfiltering MultiGeneUMI_CR \
     --soloUMIdedup 1MM_CR \
     --soloCellFilter EmptyDrops_CR \
     --soloCellReadStats Standard \
     --soloBarcodeReadLength 0 \
     --limitSjdbInsertNsj 5000000 \
     --soloFeatures GeneFull_Ex50pAS \
     --soloMultiMappers EM \
     --outReadsUnmapped Fastx

    mv ${lib2}Aligned.sortedByCoord.out.bam gex.bam
    samtools index gex.bam
    mv ${lib2}Solo.out/Barcodes.stats Barcodes.stats
    cp ${lib2}Solo.out/GeneFull_Ex50pAS/Features.stats .
    cp ${lib2}Solo.out/GeneFull_Ex50pAS/Summary.csv .
    cp ${lib2}Solo.out/GeneFull_Ex50pAS/UMIperCellSorted.txt .
    gzip -c ${lib2}Solo.out/GeneFull_Ex50pAS/CellReads.stats > CellReads.stats.gz
    mv ${lib2}Log.out STAR_Log.out
    mv ${lib2}Log.final.out STAR_Log.final.out
    gzip -c ${lib2}SJ.out.tab > STAR_SJ.out.tab.gz
    gzip ${lib2}Solo.out/GeneFull_Ex50pAS/filtered/*
    gzip ${lib2}Solo.out/GeneFull_Ex50pAS/raw/*
    mv ${lib2}Solo.out/GeneFull_Ex50pAS/filtered .
    mv ${lib2}Solo.out/GeneFull_Ex50pAS/raw .
    
    if [ -f ${lib2}Unmapped.out.mate1 ]; then
        mv ${lib2}Unmapped.out.mate1 Unmapped.out.mate1
    fi
    if [ -f ${lib2}Unmapped.out.mate2 ]; then
        mv ${lib2}Unmapped.out.mate2 Unmapped.out.mate2
    fi
    """
}

workflow align_rna_demux_species{
    main:
    
    if (!params.rna_ref_species){
        error("rna_ref_species is required")
    }
    if (!params.demux_species){
        error("demux_species output dir is required")
    }
    if (!params.rna_whitelist){
        error("rna_whitelist is required")
    }
    if (!params.rg_metadata){
        error("rg_metadata is required; source FASTQ provenance may not be omitted")
    }
    
    def read_pairs = Channel.fromFilePairs("${params.demux_species}/*/*/GEX_*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        def libn = reads[0].toString().split('/')[-3]
        def species = reads[0].toString().split('/')[-2]
        if (! rna_ref_map[species]){
            error("Species " + species + " does not have an RNA-seq reference specified")
        }
        [libn + "/" + species, id, reads, rna_ref_map[species], file(rna_ref_map[species] + "/*") ]
    }.groupTuple().map{ lib, ids, reads, refname, ref ->
        def r1s = []
        def r2s = []
        for (elt in reads){
            r1s.add(elt[0])
            r2s.add(elt[1])    
        }
        return [lib, r1s, r2s, ids, refname[0], ref[0] ]
    }.combine(Channel.fromPath(params.rna_whitelist))
    map_rna(read_pairs)
}

workflow align_rna{
    take:
    libs
    
    main:
    
    if (!params.rna_ref){
        error("RNA reference is required")
    }
    if (!params.rna_whitelist){
        error("RNA whitelist is required")
    }
    if (!params.rna_dir){
        error("RNA directory is required")
    }
    if (!params.rg_metadata){
        error("rg_metadata is required; source FASTQ provenance may not be omitted")
    }
    
    def rna_idx = Channel.fromPath("${params.rna_ref}/**").collect().map{ x -> 
        [params.rna_ref, x]}
    
    def read_pairs = Channel.fromFilePairs("${params.rna_dir}/*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        [canonical_library_name(id), id, reads]
    }.groupTuple().map{ lib, ids, reads ->
        def r1s = []
        def r2s = []
        for (elt in reads){
            r1s.add(elt[0])
            r2s.add(elt[1])    
        }
        return [lib, r1s, r2s, ids]
    }
    
    map_rna(libs.cross(read_pairs).map{ lib, tup -> tup }.combine(rna_idx).combine(Channel.fromPath(params.rna_whitelist)))
}
