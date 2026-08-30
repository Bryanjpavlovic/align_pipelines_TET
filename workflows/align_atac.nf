#! /usr/bin/env nextflow
import java.util.zip.GZIPInputStream
import java.nio.file.Files

def atac_map = [:]
if (!params.demux_species){
    if (params.atac_map){
        new File(params.atac_map).eachLine { line ->
            def (atac, lib) = line.split('\t')
            atac_map[atac] = lib
        }
    } else{
        if (!params.libs){
            error("libs file required.")
        }
        new File(params.libs).eachLine { line ->
            def linetrim = line.trim()
            atac_map[linetrim] = linetrim
        }
    }
}

def atac_ref_map = [:]
if (params.atac_ref_species){
    new File(params.atac_ref_species).eachLine { line ->
        def (species, ref) = line.split('\t')
        atac_ref_map[species] = ref
    }
}

// Load source-unit metadata: complete FASTQ basename -> minimap2 read-group line.
// Collision-safe FASTQ names and this table are produced upstream.
def rg_meta = [:]
if (params.rg_metadata){
    def first = true
    new File(params.rg_metadata).eachLine { line ->
        if (first) { first = false; return }  // skip header
        def fields = line.split('\t')
        // fields: fnbase, library, bp_id, sample_idx, lane, flowcell, lane_num, pu, rg_id, rg_string
        if (fields.length >= 10){
            rg_meta[fields[0]] = fields[9]
        }
    }
}

/**
 * Strip optional __RunTag from a filename component before extracting
 * the library name. Upstream consolidation appends the run tag as
 * __Run001 between _S##_L## and _R#_.
 *
 * After NF regex extracts match[2] (everything before _R3):
 *   Tet_2025_Multiome-ATAC_3_S3_L002__Run001 -> Tet_2025_Multiome-ATAC_3_S3_L002
 *   Tet_2025_Multiome-ATAC_3_S3_L002         -> Tet_2025_Multiome-ATAC_3_S3_L002
 *
 * fnbase keeps the full name (with run tag) as the exact metadata key.
 */
def strip_run_tag(s){
    return s.replaceFirst(/__Run\d+$/, '')
}

/** Canonicalize legacy and bcl-convert names with trailing _L# annotations. */
def canonical_library_name(s){
    def before_sample = strip_run_tag(s).replaceFirst(/_S\d+_L\d+$/, '')
    return before_sample.replaceFirst(/(?:_L\d+)+$/, '')
}

process preproc_atac_files{
    time '24h'
    
    input:
    tuple val(lib),
    val(basename),
    file(r1),
    file(r2),
    file(r3),
    file(idx),
    file(wl)
    
    output:
    tuple val(lib),
    val(basename),
    file("preproc/*_R1*.fastq.gz"),
    file("preproc/*_R2*.fastq.gz"),
    file(idx)

    script:
    """
    mkdir preproc
    ${baseDir}/atac_fq_preprocess -1 ${r1} -2 ${r2} -3 ${r3} -o preproc -w ${wl}
    """
}
process preproc_atac_files_multiome{
    time '24h'
    
    input:
    tuple val(lib),
    val(basename),
    file(r1),
    file(r2),
    file(r3),
    file(idx),
    file(wl_rna),
    file(wl_atac)
    
    output:
    tuple val(lib),
    val(basename),
    file("preproc/*_R1*.fastq.gz"),
    file("preproc/*_R2*.fastq.gz"),
    file(idx)

    script:
    """
    mkdir preproc
    ${baseDir}/atac_fq_preprocess -1 ${r1} -2 ${r2} -3 ${r3} -o preproc -w ${wl_rna} -W ${wl_atac}
    """
}

process align_atac_files{
    time { 36.hour * task.attempt }
    cpus params.threads
    memory params.memgb + ' GB'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    val(basename),
    val(rg_string),
    val(num),
    file(idx),
    file(r1),
    file(r2)
    
    output:
    tuple val(lib), file("${basename}_${num}_sorted.bam"), file("${basename}_${num}_sorted.bam.bai")

    script:
    """
    minimap2 -t ${params.threads} -y -a -x sr -R "${rg_string}" ${idx} ${r1} ${r2} \
| samtools sort -n - | samtools fixmate -m - - | samtools sort -o ${basename}_${num}_sorted.bam
    samtools index ${basename}_${num}_sorted.bam    
    """    
   
}

process cat_atac_bams{
    time { 12.hour * task.attempt }
    cpus params.threads
    memory params.memgb + ' GB'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    file(bams),
    file(bais)
    
    output:
    tuple val(lib), file("atac_merged.bam"), file("atac_merged.bam.bai")

    script:
    """
    samtools merge -@ ${params.threads} atac_merged.bam ${bams}
    samtools index atac_merged.bam
    """
}

process atac_mkdup{
    time '24h'
    cpus params.threads
    
    input:
    tuple val(lib),
    file(bam),
    file(bai)
    
    output:
    tuple val(lib),
    file("atac.bam"),
    file("atac.bam.bai")

    script:
    """
    samtools markdup -@ ${params.threads} -m t --barcode-tag CB ${bam} atac.bam
    samtools index atac.bam    
    """
}

process atac_namesort{
    time '12h'
    cpus params.threads

    input:
    tuple val(lib),
    file(bam),
    file(bai)
    
    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: {path -> path }
    
    output:
    tuple val(lib),
    file("atac.bam"),
    file("atac.bam.bai"),
    file("atac_namesort.bam")

    script:
    """
    samtools sort -@ ${params.threads} -n -o atac_namesort.bam ${bam}
    """
}

process atac_fragments{
    time { 12.hour * task.attempt }
    cpus params.threads
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    file(bam),
    file(bai),
    file(bam_namesort)
    
    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: {path -> path }
    
    output:
    tuple val(lib),
    file("atac_fragments.tsv.gz"),
    file("atac_fragments.tsv.gz.tbi")

    script:
    """
    sinto fragments -p ${params.threads} --collapse_within -m 30 -t CB \
--use_chrom "." -b ${bam} -f atac_fragments.tsv
    cat atac_fragments.tsv | sort -k1,1V -k2,2n -k3,3n | bgzip > atac_fragments.tsv.gz
    tabix -s 1 -b 2 -e 3 atac_fragments.tsv.gz 
    """
}

/**
 * Checks the first line of a gzipped FASTQ file to see if a cell barcode tag
 * has been appended.
 */
def peek_check_bc(filePath){
    def line = ""
    Files.newInputStream(filePath).withCloseable { stream ->
        new GZIPInputStream(stream).withReader { reader ->
            line = new BufferedReader(reader).readLine()        
        }
    }
    def lnsplit = line.split(' ')
    if ( lnsplit[-1] ==~ /^CB:Z:[ACGT]+$/){
        return true
    }
    else{
        return false
    }
}

workflow align_atac_demux_species{
    main:
    
    if (!params.atac_ref_species){
        error("ATAC reference is required")
    }
    if (!params.demux_species){
        error("demux_species output dir is required")
    }
    if (!params.rg_metadata){
        error("rg_metadata is required; source FASTQ provenance may not be omitted")
    }
    
    def atac_triples = Channel.fromPath("${params.demux_species}/*/*/ATAC*_R3*.fastq.gz").map{ fn -> 
        def r3 = fn.toString().trim()
        def libn = r3.split('/')[-3]
        def species = r3.split('/')[-2]
        
        def match = (r3 =~ /(.*)\/(.*)\_R3(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
        def dirn = ""
        if (match[1] != null && match[1] != ""){
            dirn += match[1] + '/'
        }
        def end = ""
        if (match[3] != null){
            end = match[3]
        }
        def gz = ""
        if (match[5] != null){
            gz = match[5]
        }
        def r1 = dirn + match[2] + "_R1" + end + '.' + match[4] + gz
        def r2 = dirn + match[2] + "_R2" + end + '.' + match[4] + gz
        def fnbase = match[2]
        if (! atac_ref_map[species]){
            error("Species " + species + " does not have an ATAC reference specified")
        }
        return [ libn + "/" + species, fnbase, file(r1), file(r2), file(r3), file(atac_ref_map[species]) ]
    }
    
    if (params.multiome){
        if (!params.rna_whitelist || !params.atac_whitelist){
            error("Both rna_whitelist and atac_whitelist are required if multiome data")
        }
        def wl_rna = Channel.fromPath(params.rna_whitelist)
        def wl_atac = Channel.fromPath(params.atac_whitelist)
        atac_preproc1 = preproc_atac_files_multiome(atac_triples.combine(wl_rna).combine(wl_atac))
    }
    else{
        if (!params.atac_whitelist){
            error("ATAC whitelist required")
        }
        def wl = Channel.fromPath(params.atac_whitelist)
        atac_preproc1 = preproc_atac_files(atac_triples.combine(wl))
    } 

    def atac_pairs = Channel.fromFilePairs("${params.demux_species}/*/*/ATAC*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        def libn = reads[0].toString().split('/')[-3]
        def species = reads[0].toString().split('/')[-2]
        if (! atac_ref_map[species]){
            error("Species " + species + " does not have an ATAC reference specified")
        }
        [ libn + "/" + species, id, reads[0], reads[1], file(atac_ref_map[species]) ]
    }.filter{ lib, fnbase, r1, r2, idx ->
        return peek_check_bc(r1) && peek_check_bc(r2)
    }
    
    atac_bams = atac_preproc1.concat(atac_pairs).map{ tup ->
        def lib = tup[0]
        def fsub = tup[1]
        if (!rg_meta.containsKey(fsub)){
            error("Missing ATAC read-group metadata for input unit: " + fsub)
        }
        def rg = rg_meta[fsub]
        [lib, fsub, rg, 0, tup[-1], tup[2], tup[3]]
    } | align_atac_files
    
    atac_bams.groupTuple() | cat_atac_bams | atac_mkdup | atac_namesort | atac_fragments
}

process split_reads_atac{
    input:
    tuple val(libname), val(fnbase), path(r1), path(r2), path(idx)

    output:
    tuple val(libname),
    val(fnbase),
    path(idx),
    path("*_R1_001.*.fastq.gz"),
    path("*_R2_001.*.fastq.gz")

    script:
    """
    ${baseDir}/split_read_files -1 ${r1} -2 ${r2} -o . -n ${params.num_chunks}
    """
}

workflow align_atac{
    take:
    libs
    
    main:
    
    if (!params.atac_ref){
        error("ATAC reference is required")
    }
    if (!params.atac_dir){
        error("ATAC reads directory is required")
    }
    if (!params.rg_metadata){
        error("rg_metadata is required; source FASTQ provenance may not be omitted")
    }
    if (params.num_chunks < 2){
        error("num_chunks must be at least 2")
    }

    def idx_atac = Channel.fromPath(params.atac_ref)
    
    // Get non-preprocessed files
    // strip_run_tag removes an optional __Run001 tag before library name extraction.
    // fnbase keeps the full name (with run tag) as the exact metadata key.
    def atac_triples = libs.cross(Channel.fromPath("${params.atac_dir}/*_R3*.fastq.gz").map{ fn ->
        def r3 = fn.toString().trim()
        def match = (r3 =~ /(.*)\/(.*)\_R3(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
        def libname_atac = canonical_library_name(match[2])
        def dirn = ""
        if (match[1] != null && match[1] != ""){
            dirn += match[1] + '/'
        }
        def end = ""
        if (match[3] != null){
            end = match[3]
        }
        def gz = ""
        if (match[5] != null){
            gz = match[5]
        }
        def r1 = dirn + match[2] + "_R1" + end + '.' + match[4] + gz
        def r2 = dirn + match[2] + "_R2" + end + '.' + match[4] + gz
        def fnbase = match[2]
        return [ atac_map[libname_atac], fnbase, file(r1), file(r2), file(r3)]
    }).map{ lib, tup -> tup }
    
    if (params.multiome){
        if (!params.atac_whitelist || !params.rna_whitelist){
            error("rna_whitelist and atac_whitelist are required for multiome data.")
        }
        def wl_rna = Channel.fromPath(params.rna_whitelist)
        def wl_atac = Channel.fromPath(params.atac_whitelist)    
        atac_preproc1 = preproc_atac_files_multiome(atac_triples.combine(idx_atac).combine(wl_rna).combine(wl_atac))
    }
    else{
        if (!params.atac_whitelist){
            error("atac_whitelist is required")
        }
        def wl = Channel.fromPath(params.atac_whitelist)
        atac_preproc1 = preproc_atac_files(atac_triples.combine(idx_atac).combine(wl))
    } 
    
    // Get pre-processed files
    // Same strip_run_tag logic for already-preprocessed files
    def atac_pairs = libs.cross(
        Channel.fromFilePairs("${params.atac_dir}/*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        [ atac_map[canonical_library_name(id)], id, reads[0], reads[1]]
    }).map{ lib, tup -> tup }.filter{ lib, fnbase, r1, r2 ->
        return peek_check_bc(r1) && peek_check_bc(r2)
    }
    atac_preproc2 = atac_pairs.combine(idx_atac)
    
    atac_bams = split_reads_atac(atac_preproc1.concat(atac_preproc2)).flatMap{ 
        ln, fsub, idx, files1, files2 ->
            if (!rg_meta.containsKey(fsub)){
                error("Missing ATAC read-group metadata for input unit: " + fsub)
            }
            def rg = rg_meta[fsub]
            files2.indices.collect{ i -> [ln, fsub, rg, i, idx, files1[i], files2[i]] }
    } | align_atac_files

    atac_bams.groupTuple() | cat_atac_bams | atac_mkdup | atac_namesort | atac_fragments

}
