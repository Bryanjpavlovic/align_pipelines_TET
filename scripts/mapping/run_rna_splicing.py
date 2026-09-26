#!/usr/bin/env python3
"""Count spliced/unspliced/ambiguous molecules from an existing STARsolo BAM.

No remapping. Uses gene-aware molecule grouping and velocyto 0.17.17's
Permissive10X splicing logic. Fragments must identify one gene; mates are
intersected before grouping by cell, gene and corrected UMI. Run per library.
The supplied GTF must be the annotation used for the original alignment.

Output: called-cell splicing.h5ad (cells x genes, X=spliced), summary, and
STAR-style splicing_filtered matrices. --raw-barcodes also exports splicing_raw
for that exact unfiltered roster. These are not native STARsolo counts.
Expression matrices are never edited. --strand follows STAR --soloStrand: forward means
the single cDNA read / first genomic mate is sense; the second mate is anti.
It is an alignment convention, not an inference from 3' versus 5' chemistry.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

RELEASE = "2026-09-26-splicing-v3"
BACKEND_VERSION = "0.17.17"
LAYERS = ("spliced", "unspliced", "ambiguous")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def cigar_segments(cigar, pos, patch_indels=3):
    """Return velocyto's inclusive segments; only M/D/N/=/X consume reference."""
    segments, merge_at = [], set()
    ref_skipped = False
    clip_left = clip_right = 0
    matches = (0, 7, 8)
    for i, (op, length) in enumerate(cigar):
        if op in matches:
            if i and cigar[i - 1][0] in matches:
                segments[-1] = (segments[-1][0], pos + length - 1)
            else:
                segments.append((pos, pos + length - 1))
            pos += length
        elif op == 3:
            ref_skipped = True
            pos += length
        elif op in (1, 2):
            if (length <= patch_indels and segments and i > 0
                    and i + 1 < len(cigar)
                    and cigar[i - 1][0] in matches
                    and cigar[i + 1][0] in matches):
                merge_at.add(len(segments) - 1)
            if op == 2:
                pos += length
        elif op == 4:
            if not segments:
                clip_left += length
            else:
                clip_right += length
        elif op not in (5, 6):
            raise ValueError(f"Unsupported CIGAR operation: {op}")
    for removed, index in enumerate(sorted(merge_at)):
        index -= removed
        start = segments.pop(index)[0]
        segments[index] = (start, segments[index][1])
    return segments, ref_skipped, clip_left, clip_right


def prepare_annotation(source, destination, references):
    """Alias contigs and number each transcript's retained exons in strand order.

    Gene/transcript IDs, names, and the annotated exonic bases remain unchanged.
    Missing or gapped exon numbers are bookkeeping, not new splice sites.
    Exact duplicate exon coordinates are written once. Immediately adjacent
    exons of the same transcript are joined because there is no intron between
    them; this prevents the backend from constructing a zero-length intron.
    Numbering statistics describe the retained, normalized exon rows.
    """
    aliases = {name: f"VCONTIG{i:08d}" for i, name in enumerate(references)}
    genes, transcripts, exons = {}, {}, {}
    attribute = re.compile(r'(\w+)\s+"([^\"]+)"')
    exon_number = re.compile(r'(^|;)\s*exon_number\s+(?:"([^\"]*)"|([^;\s]+))\s*(?:;|$)')
    stats = {"input_exon_rows": 0, "written_exons": 0,
             "duplicate_exon_rows_removed": 0, "adjacent_exon_rows_merged": 0,
             "exon_numbers_changed": 0,
             "exon_numbers_missing": 0, "transcripts_renumbered": 0}
    with open_text(source) as src:
        for line_number, line in enumerate(src, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {line_number}: expected nine tab-separated columns")
            if fields[2] != "exon":
                continue
            stats["input_exon_rows"] += 1
            chrom, strand = fields[0], fields[6]
            if chrom not in aliases:
                raise ValueError(f"GTF contig absent from BAM: {chrom}; use the matching GTF")
            if strand not in ("+", "-"):
                raise ValueError(f"GTF exon lacks a strand at line {line_number}")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"GTF line {line_number}: exon coordinates must be integers") from exc
            if start < 1 or end < start:
                raise ValueError(f"GTF line {line_number}: invalid exon coordinates {start}-{end}")
            tags = dict(attribute.findall(fields[8]))
            gene, transcript = tags.get("gene_id"), tags.get("transcript_id")
            if not gene or not transcript:
                raise ValueError(f"GTF line {line_number}: every exon needs gene_id and transcript_id")
            locus = (chrom, strand)
            if gene in genes and genes[gene] != locus:
                raise ValueError(f"gene_id {gene} spans contigs/strands; use distinct reference gene IDs")
            if transcript in transcripts and transcripts[transcript] != (gene, *locus):
                raise ValueError(f"transcript_id {transcript} is reused across genes/loci")
            genes[gene] = locus
            transcripts[transcript] = (gene, *locus)
            transcript_exons = exons.setdefault(transcript, {})
            if (start, end) in transcript_exons:
                stats["duplicate_exon_rows_removed"] += 1
                continue
            transcript_exons[start, end] = fields
    if not genes:
        raise ValueError("No transcript exons found in GTF")
    with open(destination, "w") as dst:
        for transcript, transcript_exons in exons.items():
            coordinates = sorted(transcript_exons)
            for left, right in zip(coordinates, coordinates[1:]):
                if right[0] <= left[1]:
                    raise ValueError(f"transcript_id {transcript} has overlapping nonidentical exons "
                                     f"{left[0]}-{left[1]} and {right[0]}-{right[1]}")
            joined = []
            for start, end in coordinates:
                if joined and start == joined[-1][1] + 1:
                    previous_start, previous_end = joined[-1]
                    fields = transcript_exons.pop((previous_start, previous_end))
                    fields[4] = str(end)
                    transcript_exons[previous_start, end] = fields
                    joined[-1] = (previous_start, end)
                    stats["adjacent_exon_rows_merged"] += 1
                else:
                    joined.append((start, end))
            coordinates = joined
            if transcripts[transcript][-1] == "-":
                coordinates.reverse()
            renumbered = False
            for number, coordinate in enumerate(coordinates, 1):
                fields = transcript_exons[coordinate]
                match = exon_number.search(fields[8])
                old_number = (match.group(2) or match.group(3)) if match else None
                if old_number is None:
                    stats["exon_numbers_missing"] += 1
                if old_number != str(number):
                    stats["exon_numbers_changed"] += 1
                    renumbered = True
                # Remove only the numbering attribute and preserve all others.
                attributes = exon_number.sub(r'\1', fields[8]).strip().rstrip(";").strip()
                fields[8] = f'{attributes}; exon_number "{number}";'
                fields[0] = aliases[fields[0]]
                dst.write("\t".join(fields) + "\n")
                stats["written_exons"] += 1
            stats["transcripts_renumbered"] += int(renumbered)
    annotated = {chrom for chrom, strand in genes.values()}
    stats["genes"] = len(genes)
    stats["transcripts"] = len(transcripts)
    stats["annotated_contigs"] = len(annotated)
    return {chrom: alias for chrom, alias in aliases.items() if chrom in annotated}, stats


def counter_class(vcy, pysam, sparse, aliases, strand_mode, threads):
    import numpy as np

    class FragmentRead(vcy.Read):
        __slots__ = ("fragment",)

    class StarBamCounter(vcy.ExInCounter):
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self.cellbarcode_str, self.umibarcode_str = "CB", "UB"
            self.passes = []
            self.molecule_statistics = Counter()
            self.count_cell_batch = self.count_gene_batch

        def count_gene_batch(self):
            """Resolve fragments to genes before collapsing corrected UMIs.

            A fragment compatible with several genes is excluded, rather than
            duplicated into each gene. Compatible transcript models within a
            gene are intersected across supporting fragments, as in velocyto.
            Mates are intersected first so discordant pairs cannot become two
            molecules merely because their ends overlap different genes.
            """
            molecules = {}
            pairs = {}
            stats = Counter()

            def add_fragment(barcode, umi, mappings):
                if not mappings:
                    stats["fragments_without_common_transcript"] += 1
                    return
                gene_ids = {model.geneid for model in mappings}
                if len(gene_ids) != 1:
                    stats["fragments_ambiguous_between_genes"] += 1
                    return
                gene = next(iter(gene_ids))
                key = (barcode, gene, umi)
                if key not in molecules:
                    molecules[key] = vcy.Molitem()
                molecules[key].add_mappings_record(mappings)
                stats["gene_assigned_fragments"] += 1

            self.reads_to_count.sort()
            for read in self.reads_to_count:
                mappings = self.feature_indexes[read.chrom + read.strand].find_overlapping_ivls(read)
                if not mappings:
                    stats["alignments_without_compatible_transcript"] += 1
                    continue
                if read.fragment is None:
                    add_fragment(read.bc, read.umi, mappings)
                    continue
                key = (read.bc, read.fragment)
                if key not in pairs:
                    pairs[key] = (read.umi, vcy.Molitem())
                umi, fragment = pairs[key]
                if umi != read.umi:
                    pairs[key] = (None, fragment)
                fragment.add_mappings_record(mappings)
            for (barcode, _), (umi, fragment) in pairs.items():
                if umi is None:
                    stats["fragments_with_inconsistent_corrected_UMI"] += 1
                else:
                    add_fragment(barcode, umi, fragment.mappings_record)

            barcodes = sorted(self.cell_batch)
            barcode_index = {barcode: i for i, barcode in enumerate(barcodes)}
            shape = (len(self.geneid2ix), len(barcodes))
            layers = {name: np.zeros(shape, dtype=np.uint32) for name in LAYERS}
            reuse = Counter((barcode, umi) for barcode, gene, umi in molecules)
            stats["cross_gene_UMI_groups"] += sum(n > 1 for n in reuse.values())
            for (barcode, gene, umi), molecule in molecules.items():
                result = self.logic.count(molecule, barcode_index[barcode], layers, self.geneid2ix)
                if result:
                    stats[f"unclassified_gene_molecules_code_{result}"] += 1
                else:
                    stats["counted_gene_molecules"] += 1
                    if reuse[barcode, umi] > 1:
                        stats["counted_molecules_in_cross_gene_UMI_groups"] += 1
            stats["candidate_gene_molecules"] += len(molecules)
            self.molecule_statistics.update(stats)
            # Retain only sparse count arrays once this cell batch completes.
            return {name: sparse.csr_matrix(matrix) for name, matrix in layers.items()}, barcodes

        def iter_alignments(self, bamfiles, unique=True, yield_line=False):
            for path in bamfiles:
                self._current_bamfile = str(path)
                stats = Counter()
                with pysam.AlignmentFile(path, "rb", threads=max(1, threads - 1)) as bam:
                    for record in bam:
                        stats["records"] += 1
                        if stats["records"] % 10000000 == 0:
                            logging.info("Read %d million alignments", stats["records"] // 1000000)
                        if record.flag & (4 | 256 | 512 | 2048):
                            stats["unmapped_secondary_supplementary_qcfail"] += 1
                            continue
                        if not record.has_tag("NH") or record.get_tag("NH") != 1:
                            stats["nonunique_or_missing_NH"] += 1
                            continue
                        cb = record.get_tag("CB") if record.has_tag("CB") else None
                        if cb not in self.valid_bcset:
                            stats["outside_called_cells_or_missing_CB"] += 1
                            continue
                        ub = record.get_tag("UB") if record.has_tag("UB") else ""
                        if not ub or any(base not in "ACGT" for base in ub):
                            stats["invalid_or_missing_corrected_UB"] += 1
                            continue
                        if record.reference_name not in aliases:
                            stats["contig_without_exon_annotation"] += 1
                            continue
                        if record.is_paired and record.is_read1 == record.is_read2:
                            raise ValueError("Paired BAM record must identify exactly one mate")
                        reverse = record.is_reverse ^ (record.is_paired and record.is_read2)
                        if strand_mode == "reverse":
                            reverse = not reverse
                        segments, skipped, left, right = cigar_segments(
                            record.cigartuples or [], record.reference_start + 1, vcy.PATCH_INDELS)
                        if not segments:
                            continue
                        obj = FragmentRead(cb, ub, aliases[record.reference_name],
                                       "-" if reverse else "+", record.reference_start + 1,
                                       segments, left, right, skipped)
                        obj.fragment = ((record.get_tag("RG") if record.has_tag("RG") else "",
                                         record.query_name) if record.is_paired else None)
                        stats["retained_records"] += 1
                        if record.is_paired:
                            stats["retained_paired_records"] += 1
                        yield (obj, record.to_string()) if yield_line else obj
                self.passes.append(dict(stats))
                logging.info("BAM pass: %s", json.dumps(dict(stats), sort_keys=True))
                yield (None, None) if yield_line else None

    return StarBamCounter


def annotation_classes(vcy):
    """Keep every supplied exon/intron and visit the index's final interval."""
    class FullTranscriptModel(vcy.TranscriptModel):
        def chop_if_long_intron(self, maxlen=None):
            # Legacy velocyto silently removes the 5' portion of models with
            # introns >1 Mb. Retain the actual annotation, including long genes.
            pass

    class CompleteFeatureIndex(vcy.FeatureIndex):
        def __init__(self, ivls=()):
            features = list(ivls)
            if features:
                # The legacy scan stops at index < len(ivls)-1. A nonmatching
                # sentinel makes its loop visit every real interval, including
                # a chromosome/strand containing only one annotated exon.
                features.append(vcy.Feature(start=2**63 - 1, end=2**63 - 1,
                                            kind=ord("e"), exin_no=0))
            super().__init__(features)

    return FullTranscriptModel, CompleteFeatureIndex


def sort_by_cell(args, barcodes_file, destination, samtools):
    # Filter in compiled samtools before sorting. No coordinate-BAM copy is
    # needed; the original is read directly for intron validation.
    expression = 'exists([NH]) && [NH]==1 && exists([UB]) && [UB]=~"^[ACGT]+$"'
    view = [samtools, "view", "-u", "-F", "2820", "-D", f"CB:{barcodes_file}",
            "-e", expression, str(args.bam)]
    sort = [samtools, "sort", "-@", str(max(0, args.threads - 2)),
            "-m", f"{args.sort_memory_mb}M", "-l", "1", "-t", "CB",
            "-T", str(destination.parent / "sort"), "-o", str(destination), "-"]
    logging.info("Sorting selected-barcode alignments; no remapping")
    reader = subprocess.Popen(view, stdout=subprocess.PIPE)
    try:
        result = subprocess.run(sort, stdin=reader.stdout, check=False)
        reader.stdout.close()
        rc = reader.wait()
        if rc or result.returncode:
            raise RuntimeError(f"samtools filtering/sort failed ({rc}, {result.returncode})")
    finally:
        if reader.poll() is None:
            reader.terminate()
            reader.wait()


def write_matrix_directory(destination, layers, gene_ids, gene_names, barcodes):
    """Write genes-by-cells matrices from cells-by-genes sparse count layers.

    The supplied gene and barcode order is retained, including zero-count genes
    and cells. Compatibility describes the file format, not STARsolo counting
    or cell-calling semantics. A raw export must receive its actual raw roster.
    """
    import numpy as np
    from scipy import sparse
    from scipy.io import mmwrite
    from threadpoolctl import threadpool_limits

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    gene_ids = list(gene_ids)
    gene_names = list(gene_names)
    barcodes = list(barcodes)
    shape = (len(barcodes), len(gene_ids))
    if len(gene_names) != shape[1]:
        raise ValueError("Gene names and gene IDs must have the same length")

    with gzip.open(destination / "features.tsv.gz", "wt", compresslevel=1,
                   encoding="utf-8", newline="") as stream:
        for gene_id, gene_name in zip(gene_ids, gene_names):
            stream.write(f"{gene_id}\t{gene_name or gene_id}\tGene Expression\n")
    with gzip.open(destination / "barcodes.tsv.gz", "wt", compresslevel=1,
                   encoding="utf-8", newline="") as stream:
        for barcode in barcodes:
            stream.write(f"{barcode}\n")

    summary = {"genes": shape[1], "barcodes": shape[0], "layers": {}}
    for name in ("spliced", "unspliced", "ambiguous"):
        source = layers[name]
        if not sparse.issparse(source) or source.shape != shape:
            raise ValueError(f"{name} must be a sparse cells-by-genes count matrix")
        if not np.issubdtype(source.dtype, np.integer) or np.any(source.data < 0):
            raise ValueError(f"{name} must contain nonnegative integer counts")
        matrix = source.astype(np.int64, copy=False).transpose().tocoo(copy=True)
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        with gzip.open(destination / f"{name}.mtx.gz", "wb", compresslevel=1) as stream:
            if matrix.nnz:
                # Avoid the fast Matrix Market writer starting one thread per
                # machine core inside an otherwise small Slurm allocation.
                with threadpool_limits(limits=1):
                    mmwrite(stream, matrix, field="integer", symmetry="general")
            else:
                # Some SciPy versions emit a real header for empty matrices
                # even with field="integer". Retain the integer file contract.
                stream.write(("%%MatrixMarket matrix coordinate integer general\n%\n"
                              f"{shape[1]} {shape[0]} 0\n").encode("ascii"))
        summary["layers"][name] = {
            "nonzero_entries": int(matrix.nnz),
            "molecules": int(matrix.data.sum(dtype=np.int64)),
        }
    return summary


def read_barcodes(path):
    with open_text(path) as handle:
        barcodes = [line.strip().split("\t")[0] for line in handle if line.strip()]
    if not barcodes or len(barcodes) != len(set(barcodes)):
        raise ValueError(f"Barcode roster must be nonempty and unique: {path}")
    return barcodes


def run(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import pysam
    from scipy import sparse
    import velocyto as vcy

    installed = importlib.metadata.version("velocyto")
    if installed != BACKEND_VERSION:
        raise ValueError(f"This adapter requires velocyto=={BACKEND_VERSION}; found {installed}")
    samtools = shutil.which("samtools")
    if not samtools:
        raise ValueError("samtools (>=1.13) must be on PATH")
    barcodes = read_barcodes(args.barcodes)
    count_barcodes = read_barcodes(args.raw_barcodes) if args.raw_barcodes else barcodes
    count_index = {barcode: i for i, barcode in enumerate(count_barcodes)}
    if not set(barcodes).issubset(count_index):
        raise ValueError("The raw barcode roster must include every called-cell barcode")
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        if bam.header.to_dict().get("HD", {}).get("SO") != "coordinate":
            raise ValueError("Input must be the coordinate-sorted genomic BAM")
        references = bam.references
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "splicing.h5ad"
    if output.exists():
        raise ValueError(f"Output already exists: {output}; choose a new output directory")
    matrix_root = args.matrix_output_dir or args.output_dir
    for name in ("splicing_filtered", "splicing_raw") if args.raw_barcodes else ("splicing_filtered",):
        destination = matrix_root / name
        if destination.exists() and any(destination.iterdir()):
            raise ValueError(f"Matrix output already exists: {destination}; choose a new --matrix-output-dir")
    scratch = args.work_dir or args.output_dir
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="splicing_work_", dir=scratch) as temp:
        work = Path(temp)
        annotation = work / "annotation.gtf"
        aliases, annotation_stats = prepare_annotation(args.gtf, annotation, references)
        roster = work / "barcodes.tsv"
        roster.write_text("\n".join(count_barcodes) + "\n")
        CounterClass = counter_class(vcy, pysam, sparse, aliases, args.strand, args.threads)
        counter = CounterClass(sampleid=args.sample_id, logic=vcy.Permissive10X,
                               valid_bcset=set(barcodes), loom_numeric_dtype="uint32",
                               outputfolder=str(work))
        original_classes = vcy.TranscriptModel, vcy.FeatureIndex
        vcy.TranscriptModel, vcy.FeatureIndex = annotation_classes(vcy)
        try:
            counter.read_transcriptmodels(str(annotation))
            logging.info("Validating annotated introns against existing alignments")
            counter.mark_up_introns((str(args.bam),), multimap=False)
            if not counter.passes or not counter.passes[0].get("retained_records", 0):
                raise ValueError("No usable NH=1, called-cell, corrected-UB alignments found")
            # Keep intron validation based on called cells in both modes, so
            # requesting additional raw columns never changes filtered counts.
            counter.valid_bcset = set(count_barcodes)
            sorted_bam = work / "cellsorted.bam"
            sort_by_cell(args, roster, sorted_bam, samtools)
            logging.info("Counting gene-aware molecules using Permissive10X classification")
            batches, observed_barcodes = counter.count((str(sorted_bam),), multimap=False,
                                                      cell_batch_size=args.cell_batch_size)
        finally:
            vcy.TranscriptModel, vcy.FeatureIndex = original_classes
        if len(observed_barcodes) != len(set(observed_barcodes)):
            raise RuntimeError("Barcode sorting produced repeated cell batches")
        reorder = sparse.csr_matrix(
            (np.ones(len(observed_barcodes), dtype=np.uint32),
             ([count_index[barcode] for barcode in observed_barcodes],
              np.arange(len(observed_barcodes)))),
            shape=(len(count_barcodes), len(observed_barcodes)))
        gene_ids = sorted(counter.geneid2ix, key=counter.geneid2ix.get)
        counted_layers = {}
        for key in LAYERS:
            observed = sparse.hstack(batches.pop(key), format="csr").T.tocsr()
            counted_layers[key] = (reorder @ observed).tocsr()
        selected = np.asarray([count_index[barcode] for barcode in barcodes], dtype=np.int64)
        layers = {key: matrix[selected].tocsr() for key, matrix in counted_layers.items()}
        reverse_alias = {alias: chrom for chrom, alias in aliases.items()}
        genes = [counter.genes[gene] for gene in gene_ids]
        cell_ids = [f"{args.sample_id}:{barcode}" for barcode in barcodes]
        obs = pd.DataFrame({"barcode": barcodes, "library": args.sample_id}, index=cell_ids)
        obs.index.name = "cell_id"
        var = pd.DataFrame({"gene_name": [gene.genename for gene in genes],
                            "contig": [reverse_alias[gene.chrom] for gene in genes],
                            "strand": [gene.strand for gene in genes]}, index=gene_ids)
        var.index.name = "gene_id"
        total = sum((matrix.astype(np.uint64) for matrix in layers.values()),
                    sparse.csr_matrix((len(barcodes), len(gene_ids)), dtype=np.uint64))
        summary = {
            "release": RELEASE, "backend": f"velocyto {installed}",
            "logic": "Permissive10X", "bam": str(args.bam), "gtf": str(args.gtf),
            "called_barcodes": str(args.barcodes), "sample_id": args.sample_id,
            "raw_barcodes": str(args.raw_barcodes) if args.raw_barcodes else "not requested",
            "counted_barcode_columns": len(count_barcodes),
            "strand_relative_to_first_genomic_mate": args.strand,
            "barcode_tag": "CB", "umi_tag": "UB", "alignment_filter": "primary NH=1; no QC-fail",
            "multimapper_EM": False, "remapped": False,
            "annotation_policy": "complete supplied models; adjacent exons joined without changing exonic bases; "
                                 "normalized exon numbers (statistics describe retained rows); no long-intron truncation",
            "annotation_normalization": annotation_stats,
            "molecule_key": "CB + uniquely assigned gene_id + corrected UB",
            "fragment_policy": "intersect annotated mates; discard multi-gene fragments; allow one usable mate",
            "intron_validation": "called-cell alignments only, including when raw output is requested",
            "X_definition": "spliced; original expression matrix is separate",
            "cells": len(barcodes), "genes": len(gene_ids),
            "cells_without_retained_alignments": len(set(barcodes) - set(observed_barcodes)),
            "cells_with_zero_splicing_counts": int(np.count_nonzero(np.asarray(total.sum(axis=1)).ravel() == 0)),
            "molecules": {key: int(matrix.sum(dtype=np.uint64)) for key, matrix in layers.items()},
            "input_alignment_statistics": counter.passes[0],
            "molecule_statistics": dict(counter.molecule_statistics),
            "molecule_statistics_scope": "all selected barcodes, including raw when requested",
        }
        if not total.nnz:
            raise ValueError("No called-cell molecules classified; check the matching GTF and --strand")
        gene_names = [gene.genename for gene in genes]
        summary["filtered_matrix_export"] = write_matrix_directory(
            matrix_root / "splicing_filtered", layers, gene_ids, gene_names, barcodes)
        if args.raw_barcodes:
            summary["raw_matrix_export"] = write_matrix_directory(
                matrix_root / "splicing_raw", counted_layers, gene_ids, gene_names, count_barcodes)
        data = ad.AnnData(X=layers["spliced"].copy(), obs=obs, var=var, layers=layers)
        data.uns["splicing"] = summary
        temporary_output = args.output_dir / "splicing.partial.h5ad"
        data.write_h5ad(temporary_output, compression="gzip")
        temporary_output.replace(output)
        (args.output_dir / "counts_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"COMPLETE: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for option in ("bam", "gtf", "barcodes", "output-dir"):
        parser.add_argument(f"--{option}", required=True, type=Path)
    parser.add_argument("--raw-barcodes", type=Path,
                        help="Unfiltered RNA barcode roster, including called cells; also write splicing_raw")
    parser.add_argument("--matrix-output-dir", type=Path,
                        help="Parent for splicing_filtered/raw; default --output-dir")
    parser.add_argument("--sample-id", required=True, help="Library identity; run independent libraries separately")
    parser.add_argument("--strand", choices=("forward", "reverse"), default="forward",
                        help="Match original STAR --soloStrand, relative to first genomic mate (default: forward)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--sort-memory-mb", type=int, default=2048, help="samtools sort memory per thread")
    parser.add_argument("--cell-batch-size", type=int, default=25, help="Cells retained as read objects per counting batch")
    parser.add_argument("--work-dir", type=Path, help="Temporary barcode-sort directory; default output directory (use BeeGFS)")
    parser.add_argument("--version", action="version", version=RELEASE)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        for name in ("bam", "gtf", "barcodes", "raw_barcodes", "output_dir", "matrix_output_dir", "work_dir"):
            path = getattr(args, name)
            if path is not None:
                setattr(args, name, path.expanduser().resolve())
        if args.threads < 2 or args.sort_memory_mb < 1 or args.cell_batch_size < 1:
            raise ValueError("threads must be >=2; sort memory and batch size must be positive")
        run(args)
    except (ValueError, RuntimeError, OSError, ImportError) as exc:
        logging.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
