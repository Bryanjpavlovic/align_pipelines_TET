# align_pipelines

Nextflow workflows for single-cell RNA-seq, single-cell ATAC-seq, and genomic
DNA alignment and variant calling without Cell Ranger.

The local workflow set adds source-FASTQ provenance, collision-safe handling of
consolidated FASTQ names, selectable RNA read geometry, and matching behavior
across RNA, ATAC, and DNA. See [CHANGELOG.md](CHANGELOG.md) for the complete
history relative to the supplied original repository copies.

## Pipeline layout

- `align_pipelines.nf`: RNA, ATAC, and DNA entrypoint
- `make_ref.nf`: STAR and minimap2 reference builder
- `nextflow.config`: shared parameters and task environment activation
- `nextflow.config.slurm`: optional SLURM overlay
- `workflows/align_rna.nf`: STARsolo RNA mapping
- `workflows/align_atac.nf`: barcode preprocessing and minimap2 ATAC mapping
- `workflows/align_dna.nf`: minimap2 DNA mapping and FreeBayes calling
- `scripts/mapping/`: 10X RNA/ATAC staging, mapping, QC, and plotting control layer
- `examples/`: parameter-file templates

## Cluster deployment

The deployed cluster copy is normally loaded with:

```bash
module purge
module load miniforge/3 nextflow/latest htslib/1.20 samtools/1.20 align_pipelines/bjp

PIPELINE_DIR="${ALIGN_PIPELINES_HOME:-/nvme/software/packages/align_pipelines/bjp}"
```

The deployment directory must contain the compiled helper programs
`atac_fq_preprocess`, `split_read_files`, `vcf_depth_filter`, and
`bin/rna_bam_evidence` in addition to
the scripts listed above. Install the Python mapping control layer under
`$PIPELINE_DIR/bin/mapping`; it is part of this repository and is no longer
installed by CellBouncer.

For a standalone source build, install the bioinformatics tools used by the
workflows, make the `htswrapper` dependency available at
`dependencies/htswrapper`, and run `make`. Existing helpers use C++11; the RNA
evidence profiler target uses C++17 and links directly against htslib and zlib.
The package-scoped deployment target is:

```bash
make rna_bam_evidence
make install-bjp BJP_PREFIX=/nvme/software/packages/align_pipelines/bjp
```

That target installs only into `align_pipelines/bjp`.

## 10X mapping orchestration and QC

The high-level 10X runner is:

```bash
python3 "$PIPELINE_DIR/bin/mapping/orchestrate_10x_mapping_qc.py" --help
```

Its mapping helpers are installed beside it. From the source checkout, the same
runner can be invoked from `scripts/mapping/orchestrate_10x_mapping_qc.py`.
Both layouts resolve helper scripts from the runner's directory and resolve
`align_rna.nf` and `align_atac.nf` from the repository/package `workflows/`
directory. The orchestrator retains explicit override options for unusual
layouts, but the standard source and installed layouts need none.

The old `mapping/NextflowConfigs` copy is not part of the finalized layout.
Top-level configuration files and `workflows/` are the authoritative copies.

### Opt-in one-pass RNA BAM evidence

The broad evidence profiler is intentionally opt-in. One orchestrator
invocation creates the complete harvesting DAG:

1. lightweight STAR diagnostic promotion;
2. one all-library BAM-reader array bounded by
   `--rna3-bam-evidence-max-concurrent` (default and absolute ceiling: 3);
3. BAM-free gather, trimming joins, plots, and final validation after every
   reader succeeds.

There are no pilot, benchmark, phase, continuation-audit, or manual-release
steps.

```bash
python3 "$PIPELINE_DIR/bin/mapping/orchestrate_10x_mapping_qc.py" \
  --run-name rna3_all40_full_saturation_20260830_v2 \
  --rna3-runs BP12952 BP12953 BP14051 BP14052 BP16462 BP16480 BP16481 BP16482 \
  --rna3-bam-evidence-from-bam \
  --rna3-bam-evidence-baseline-root \
    /mnt/beegfs/tetmultiome_rna_mapped/mapping_output \
  --rna3-bam-evidence-source-order \
    BP12952 BP12953 BP14051 BP14052 BP16462 BP16480 BP16481 BP16482 \
  --rna3-bam-evidence-class-manifest /path/validated_classes.tsv \
  --rna3-bam-evidence-cpus 4 \
  --rna3-bam-evidence-memory-gb 48 \
  --rna3-bam-evidence-max-concurrent 3 \
  --rna3-bam-evidence-hash-bins 100 \
  --resume \
  --submit
```

Use `--rna3-bam-evidence-no-biological-classification` instead of a manifest
only when unavailable mitochondrial, rRNA, species, and reference-class values
are intentional. Exactly one classification choice is required. A normalized
class TSV uses `entity_type` (`contig` or `feature`), `identifier`, and
optional `species`, `mitochondrial`, `rrna`, and `reference_class`
columns. Legacy `contig`, `feature_id`, or `GX` identifier columns are
also accepted.

To replace outputs from the rejected pre-2.3 deployment, first ensure every
recorded RNA-evidence and diagnostic job is terminal, then add
`--rna3-bam-evidence-reset-failed-run` to the command once. Recovery archives
the old scientific contract, scripts, per-library evidence directories, and
gathered output below `control/rna3_bam_evidence_failed_run_migrations/`.
Mapping data are never removed. The replacement invocation submits the complete
bounded array and its downstream dependency chain.

Dry-run generation remains the default. Add `--submit` to submit the complete
DAG. Existing outputs are reused only after full product and audit
revalidation. The deprecated `--rna3-cell-reads-backfill-from-bam` spelling
remains an alias to this compiled stage and never reaches the retired AWK
scanner.

The independent `rna3_star_diagnostics` job matches retained Nextflow work tasks
to published `Summary.csv` digests and atomically promotes `STAR_Log.out`,
`STAR_Log.final.out`, and `STAR_SJ.out.tab.gz`. Ambiguous matches are fatal. It
never opens a BAM and never deletes the work directory. Manual deletion remains
allowed only after complete validation. Reuse requires a manifest-bound JSON
marker and audit plus current hashes for all three outputs; arbitrary preexisting
diagnostic files or a completed Slurm ID alone are not accepted.

The generated manifest and compiled command are explicitly bound to the
producing STARsolo settings `GeneFull_Ex50pAS`, `MultiGeneUMI_CR`, `1MM_CR`,
and `EM`; a missing or different declaration is fatal before the BAM is opened.
The `EM` declaration records the producing run but does not make the standard
uppercase BAM tags an exact representation of STARsolo's separate multi-gene
EM matrix. The compiled pass counts primary mapped reads, validates required
tags, keeps all corrected barcodes rather than only called cells, and examines
all mapped alignment records while deduplicating ordinary matrix molecules by
`(CB,GX,UB)`. Here `UB` is STARsolo's
accepted, corrected UMI after `1MM_CR` correction and `MultiGeneUMI_CR`
filtering, not the raw `UR`. Each
molecule slot contains two 32-bit interned IDs,
a 64-bit packed/interned UMI token, separate 64-bit RG and physical-source
presence masks, and the minimum stable 64-bit QNAME hash. Multiple lane RGs
from one `bp_id` share its chronological source bit. The 40-byte slot at a
maximum 0.70 load factor is approximately 57.1 table bytes per molecule,
excluding other structures. A whole-process `RLIMIT_AS` ceiling is set to 90%
of the SLURM request. Molecules, all interners, barcode metrics, barcode-RG
metrics, correction aggregates, dense numeric contig counters, rehash peaks,
and sorting vectors share a conservative 70% admission budget inside that
ceiling. The remainder is reserved for htslib, BGZF threads, libc, allocator
overhead, and I/O. High-cardinality tables are compacted and sorted in place.
Each record's tags and CIGAR contribution are decoded once and reused for its
library, RG, barcode, barcode-RG, and contig accumulators.
The process fails before an unsafe allocation instead of exceeding the task
budget, and its audit records component peaks, cardinalities, elapsed time, and
complete process MaxRSS.

`fnv1a64_seeded_v1` starts from the 64-bit FNV-1a offset basis XOR the recorded
seed, then XORs each unsigned QNAME byte and multiplies by the 64-bit FNV prime
modulo 2^64. A nested bin is `floor(hash * bins / 2^64)`. The scientific hash
does not use `std::hash`; the seed and bin count are recorded in every audit.
Ordinary STARsolo counted-read evidence is one logical read per `(RG,QNAME)`
with a declared RG, valid `CB`, and one unambiguous `GX` in the ordinary feature
roster. The NH=1 path selects one fragment representative; the NH>1 path merges
primary and secondary alignments for the exact `(RG,QNAME)`, including a gene
tag exposed only on a secondary record. `NH` is unrestricted, and a valid `UB`
is not required because STARsolo records `countedU` before UMI
collapse/filtering. Ordinary
matrix molecules are distinct
`(CB,GX,accepted-corrected-UB)` tuples and are also unrestricted by `NH`; a `UB`
null sentinel from rejected `MultiGeneUMI_CR` evidence cannot form a molecule.
The C++ pass exactly merge-checks every raw and filtered MatrixMarket coordinate
and count; the wrapper also checks filtered per-barcode UMIs, called-cell count,
ordinary read total, STARsolo median, and native `CellReads.stats` when present.
Only a passing gate publishes `CellReads.countedU.from_bam.tsv.gz`.

The separately labelled `nh_gt1_unique_gene_countedU_reads` metric is the
`NH>1` subset of those ordinary singleton-gene counted reads; it is included in
the ordinary total and is not an EM membership test. Exact multi-gene evidence
for STARsolo's `UniqueAndMult-EM.mtx` is unavailable from the standard uppercase
`GX`/`UB` BAM tags and is not claimed or emitted. The source-derived fixture
includes an `NH>1` singleton-gene read in the ordinary outputs, a one-mismatch
raw `UR` collapsed to an existing corrected `UB`, and an equal-support
cross-gene UMI rejected with `UB:-`.

When no class manifest applies, class counts and fractions are blank and status
is `unavailable`; zero is never used to mean not measured. Feature/GX classes
take precedence over contig classes when both are supplied. No class name is
inferred heuristically.

The trimming aggregator accepts an optional
`barcode_correction_bridge` manifest column. When present it consumes
`raw_to_corrected_barcode_counts.tsv.gz`, reports aggregate conflicts as
unresolved only when they remain ambiguous within the matching RG, and does
not reopen the BAM. Cross-RG differences are resolved by the manifest's exact
RG/source key. Its historical BAM/read-name path is retained only for runs
without the bridge.

Focused tests live under `tests/rna_bam_evidence`. Pure Python and script
generation checks run everywhere. The complete tagged-BAM test requires an
htslib 1.20 profiler and samtools 1.20. The test compiles against real htslib
automatically when `pkg-config htslib` is available; it never uses fake headers:

```bash
module purge
module load miniforge/3 genomics-base/latest
module load htslib/1.20 samtools/1.20

make rna_bam_evidence
RNA_BAM_EVIDENCE_BIN="$PWD/rna_bam_evidence" \
  python3 -m unittest discover -s tests/rna_bam_evidence -v
```

The BAM-free gather revalidates every per-library source and every gathered
output against recorded schema, sort order, row count, byte size, and SHA-256,
including on resume. It stops at reusable evidence tables. Run downstream
questions independently, without reopening BAMs. Both large analyses accept a
repeatable `--library` option for independent shards:

```bash
MAPPING="$PIPELINE_DIR/bin/mapping"

python3 "$MAPPING/analyze_gained_cells.py" \
  --current-root /path/current/mapping_output \
  --baseline-root /path/historical/mapping_output \
  --output-dir /path/qc/rna3/gained_cells

python3 "$MAPPING/analyze_depth_history.py" \
  --current-root /path/current/mapping_output \
  --baseline-root /path/historical/mapping_output \
  --fastq-inventory /path/raw_fastq_inventory.tsv \
  --output-dir /path/qc/rna3/depth_history

python3 "$MAPPING/join_trim_cell_metrics.py" \
  --current-root /path/current/mapping_output \
  --trim-by-source-tsv /path/all_libraries_trim_by_barcode_by_source.tsv.gz \
  --output-dir /path/qc/rna3/trim_cell_links

python3 "$MAPPING/export_repooling_evidence.py" \
  --fastq-inventory /path/raw_fastq_inventory.tsv \
  --source-yield-tsv /path/qc/rna3/bam_evidence/all_libraries_source_yield.tsv \
  --output-dir /path/qc/rna3/repooling
```

`analyze_gained_cells.py` also accepts optional demultiplexing, trim, and
precomputed ambient-similarity tables keyed by library and barcode, so those
joins can be refreshed later without BAM access. It writes rows incrementally
per library through SQLite shards and reads raw matrices through a disk-backed
12-byte-per-column numeric mmap, not raw-barcode Python objects. All shared,
gained, and lost cells are
detailed. Background means neither filtered set plus nonzero current/historical
raw matrix evidence or current/historical BAM-observed corrected-barcode
evidence. Only a deterministic bounded background sample is detailed, while
exact fixed-bin distributions use every evidence-bearing background barcode.
Zero-only raw-whitelist entries are not exported.

`analyze_depth_history.py` uses fixed-roster numeric vectors, labels native
observed endpoints, chronological physical-source threshold proxies, and
final-BAM random-thinning threshold proxies separately, and exports median UMI
and read trajectories. Lane RGs sharing a `bp_id` are combined into one source
prefix. `export_repooling_evidence.py` preserves every original inventory column
and labels raw FASTQ representation as authoritative.

`plot_rna_evidence.py` accepts any available gained-cell, depth, source-yield,
and existing `processedstats.tsv` tables. It renders real fixed-bin distributions
for read, UMI, gene, saturation, mapping/gene-assignment, mitochondrial, rRNA,
source-dominance, and trimming measures by cell category. It also renders three
separately labelled cell-count trajectory panels, fixed-roster UMI/read panels,
and six library panels for mapping, unique mapping, gene assignment, valid
barcode, cell-associated reads, and saturation. Histogram bins, quantiles, and
all trajectory values are written to companion TSVs. BAM duplicate flags are
never described as UMI or PCR duplication, and final UMI/gene counts remain
authoritative to STARsolo matrices.

## Required source-unit metadata

Every RNA, ATAC, or DNA mapping run requires `rg_metadata`. It is a tab-delimited
table with this exact header and column order:

```text
fnbase  library  bp_id  sample_idx  lane  flowcell  lane_num  pu  rg_id  rg_string
```

`fnbase` is the complete consolidated FASTQ basename before `_R1`, `_R2`, or
`_R3`, including a collision-avoidance suffix such as `__Run001` when present.
The workflows use that exact key to assign a source-unit read group. They remove
the run tag and terminal sequencing suffix only when deriving the logical
library name.

`rg_string` is the serialized minimap2 read-group argument, normally beginning
with `@RG` and containing literal `\t` escapes between fields. RNA constructs
the corresponding STAR fields from `rg_id`, `library`, `pu`, and `bp_id`.

A missing metadata row is fatal. Retain this table with the run outputs because
it is the durable link from each BAM read group to its source FASTQ unit.

## Parameter files

Copy a supplied template and edit it outside the installed software directory:

```bash
cp "$PIPELINE_DIR/examples/example.yml" my_mapping_params.yml
```

For `align_pipelines.nf`, always set:

- `output_directory`
- `rg_metadata`
- `libs`, unless processing CellBouncer `demux_species` output
- the input directory and reference parameters for each modality being run

Remove unused modality sections from the copied YAML. If `rna_ref`, `atac_ref`,
and `dna_ref` are all present, the entrypoint will launch all three workflows.

### RNA geometry

`rna_geometry` accepts two values:

- `long-r2` (default): separate R1 barcode/UMI read and R2 cDNA; STAR receives
  cDNA R2 first and barcode R1 last.
- `pe150`: CB16+UMI12 at the start of mate 1 followed by retained R1 cDNA; STAR
  clips the first 28 bases and maps both mates.

### FASTQ naming and multiple sequencing runs

Upstream FASTQ consolidation may append `__RunNNN` when otherwise identical
basenames occur in multiple runs. Do not concatenate or rename those files after
the metadata table is created. The workflows group them into the same logical
library while preserving each full basename as a distinct read-group lookup
key.

## Running on SLURM

The project-level `nextflow.config` supplies defaults. Add the SLURM settings as
a soft configuration overlay:

```bash
nextflow -c "$PIPELINE_DIR/nextflow.config.slurm" \
    run "$PIPELINE_DIR/align_pipelines.nf" \
    -params-file my_mapping_params.yml
```

The overlay uses the `compute` partition, preserves the work cache, enables
resume, limits Nextflow to ten active tasks, and limits submission to ten jobs
per minute. These are Nextflow executor limits, not SLURM array limits.

To build STAR and minimap2 references:

```bash
cp "$PIPELINE_DIR/examples/example_make_ref.yml" my_ref_params.yml

nextflow -c "$PIPELINE_DIR/nextflow.config.slurm" \
    run "$PIPELINE_DIR/make_ref.nf" \
    -params-file my_ref_params.yml
```

For CellBouncer `demux_species` output, use
`examples/example_demux_species.yml` with `align_pipelines.nf`.

## Outputs and provenance

- RNA publishes `gex.bam`, its index, STARsolo summaries and matrices, and
  optional unmapped FASTQs under each logical library.
- ATAC publishes `atac.bam`, its index, name-sorted BAM, and fragment files
  under each logical library.
- DNA publishes library BAMs and combined, indexed variant-call outputs.

Final BAM filenames remain library-scoped. Source FASTQ identity is stored in
BAM read groups rather than being embedded in the merged BAM filename or read
QNAME.

## Resume and cleanup

The SLURM overlay sets `resume = true` and `cleanup = false`. Preserve the run's
Nextflow `work/` directory while the run may need to resume. After successful
completion and final validation, use `nextflow clean -n` to preview eligible
cache removal before deleting anything.

## ATAC preprocessing safety

`atac_fq_preprocess` writes R1 and R2 with their original basenames beneath
`--output_dir`. The program refuses an output directory that would resolve an
output to any input FASTQ, preventing gzip write mode from truncating a source
file.

