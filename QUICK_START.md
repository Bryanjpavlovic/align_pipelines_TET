# align_pipelines quick start

## 1. Load the cluster software

```bash
module purge
module load miniforge/3
module load nextflow/latest
module load align_pipelines/latest

PIPELINE_DIR="${ALIGN_PIPELINES_HOME:-/nvme/software/packages/align_pipelines/bjp}"
```

`PIPELINE_DIR` must contain `align_pipelines.nf`, `make_ref.nf`,
`nextflow.config`, and the `workflows/` directory.

## 2. Build a reference

```bash
cp "$PIPELINE_DIR/examples/example_make_ref.yml" my_ref_params.yml
# Edit paths, genome_base, memgb, and threads in my_ref_params.yml.

nextflow -c "$PIPELINE_DIR/nextflow.config.slurm" \
    run "$PIPELINE_DIR/make_ref.nf" \
    -params-file my_ref_params.yml
```

## 3. Run RNA, ATAC, or DNA mapping

```bash
cp "$PIPELINE_DIR/examples/example.yml" my_mapping_params.yml
# Retain only the modality sections you intend to run and fill in every path.

nextflow -c "$PIPELINE_DIR/nextflow.config.slurm" \
    run "$PIPELINE_DIR/align_pipelines.nf" \
    -params-file my_mapping_params.yml
```

The mapping entrypoint requires `rg_metadata`, the ten-column source-unit table
generated during upstream FASTQ consolidation:

```text
fnbase  library  bp_id  sample_idx  lane  flowcell  lane_num  pu  rg_id  rg_string
```

The complete consolidated FASTQ basename, including any `__RunNNN` suffix, must
match `fnbase` exactly. Keep this table with the run outputs because it links BAM
read groups back to the source FASTQ unit.

For RNA, set `rna_geometry` to:

- `long-r2` for the standard separate R1 barcode/UMI read and R2 cDNA layout;
- `pe150` for CB16+UMI12 embedded in mate 1 followed by retained R1 cDNA.

## 4. Use the bundled 10X mapping orchestrator

For staged 10X RNA/ATAC trimming, mapping, QC, plotting, and validation, use:

```bash
python3 "$PIPELINE_DIR/bin/mapping/orchestrate_10x_mapping_qc.py" --help
```

The orchestrator and all seven runtime helpers now belong to this repository.
They are installed together below `$PIPELINE_DIR/bin/mapping` and resolve the
RNA/ATAC workflows from `$PIPELINE_DIR/workflows` automatically.

## 5. Submit the Nextflow controller through SLURM

If your cluster policy requires the Nextflow controller itself to run as a
SLURM job, create `submit_align.sbatch`:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=align_pipeline
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00

set -euo pipefail
module purge
module load miniforge/3
module load nextflow/latest
module load align_pipelines/latest

PIPELINE_DIR="${ALIGN_PIPELINES_HOME:-/nvme/software/packages/align_pipelines/bjp}"

nextflow -c "$PIPELINE_DIR/nextflow.config.slurm" \
    run "$PIPELINE_DIR/align_pipelines.nf" \
    -params-file my_mapping_params.yml
```

Submit it with:

```bash
sbatch submit_align.sbatch
```

The SLURM overlay limits Nextflow to ten active tasks and ten submissions per
minute. Those settings throttle Nextflow tasks; they do not control SLURM array
concurrency. The bundled mapping orchestrator provides
`--array-max-concurrent` for its arrays and `--max-cores` for its supported
RNA-only total-core ceiling.

## 6. Resume and clean up

Resume is enabled by the supplied SLURM overlay. You may also pass `-resume`
explicitly. Do not delete the run's `work/` directory until the run is complete
and you no longer need its resume cache.

## 7. Supplied parameter templates

- `examples/example.yml`: combined RNA, ATAC, and DNA template
- `examples/example_call_vars.yml`: DNA mapping and variant calling
- `examples/example_demux_species.yml`: CellBouncer `demux_species` output
- `examples/example_make_ref.yml`: STAR and minimap2 reference construction
