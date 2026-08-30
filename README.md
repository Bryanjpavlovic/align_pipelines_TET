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
- `examples/`: parameter-file templates

## Cluster deployment

The deployed cluster copy is normally loaded with:

```bash
module purge
module load miniforge/3
module load nextflow/latest
module load align_pipelines/latest

PIPELINE_DIR="${ALIGN_PIPELINES_HOME:-/nvme/software/packages/align_pipelines/bjp}"
```

The deployment directory must contain the compiled helper programs
`atac_fq_preprocess`, `split_read_files`, and `vcf_depth_filter` in addition to
the scripts listed above.

For a standalone source build, install the bioinformatics tools used by the
workflows, make the `htswrapper` dependency available at
`dependencies/htswrapper`, and run `make`. The Makefile uses C++11 and links
against zlib and htslib.

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
