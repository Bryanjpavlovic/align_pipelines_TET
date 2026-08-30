# Alignment pipeline changelog

This changelog records local workflow, helper-program, configuration,
documentation, and packaging changes relative to the supplied original
`nkschaefer/align_pipelines` repository copies. It replaces the former RNA-only
changelog and the historical `align_rna.nf` backup files.

## 2026-08-30 — Final workflow consolidation

### Shared source-FASTQ provenance

- RNA, ATAC, and DNA mapping now require `rg_metadata`, a tab-delimited table
  with this header and column order:

  ```text
  fnbase  library  bp_id  sample_idx  lane  flowcell  lane_num  pu  rg_id  rg_string
  ```

- Upstream FASTQ consolidation may append `__RunNNN` to an otherwise repeated
  FASTQ basename. The mapping workflows do not create this suffix and do not
  rewrite read names.
- Each workflow removes `__RunNNN`, the terminal `_S#_L#` sequencing suffix,
  and supported pre-sample `_L#` annotations only when deriving the logical
  library name. The complete, unmodified basename remains the exact metadata
  lookup key.
- The source association is carried into BAM read groups. RNA constructs STAR
  read-group fields from `rg_id`, `library`, `pu`, and `bp_id`; ATAC and DNA pass
  the complete serialized `rg_string` to minimap2.
- A missing metadata row is fatal. The metadata TSV must be retained with the
  mapping run because it is part of the link between a BAM read group and its
  source FASTQ unit.
- Final BAM filenames remain library-scoped (`gex.bam`, `atac.bam`, or the DNA
  library name). Source provenance is stored in BAM read groups rather than in
  the merged BAM filename or read QNAME.

### RNA (`workflows/align_rna.nf`)

- Preserved the correct geometry-dependent STARsolo input order:
  - `long-r2` (default): cDNA R2 first, separate barcode/UMI R1 last, with
    `--soloBarcodeMate 0`.
  - `pe150`: R1 then R2, with CB16+UMI12 embedded in mate 1,
    `--soloBarcodeMate 1`, and `--clip5pNbases 28 0`.
- Added ordered STAR read groups for comma-separated FASTQ inputs, with
  `ID`, `SM`, `PL`, `PU`, and `DS` fields.
- Added explicit provenance preflight to both normal and species-demultiplexed
  entry points.
- Retained the current RNA mapping policies introduced by the historical
  versions: Cell Ranger 4 adapter clipping, normalized alignment thresholds of
  0.33, and separate optional unmapped FASTQ output.
- Modernized process output declarations from `file` to `path` and made
  unmapped-output publication conditional so absent optional mate files do not
  fail the process.

### ATAC (`workflows/align_atac.nf`)

- Replaced the original inline minimal read group with the complete metadata
  `rg_string` supplied to minimap2 `-R`.
- Applied one source-unit read group to every chunk derived from that FASTQ
  unit, preserving provenance through sort, merge, fixmate, and duplicate
  marking.
- Preserved minimap2 `-y`, which carries ATAC FASTQ header comments such as the
  cell-barcode tag; this behavior is ATAC-specific and is not copied to DNA.
- Added collision-safe library-name canonicalization for both raw R3 inputs and
  already-preprocessed R1/R2 inputs.
- Added explicit provenance preflight to the species-demultiplexed entry point.
- Retained the species-demultiplexed tuple repair that supplies a chunk number
  and read group to `align_atac_files`.
- Made the barcode-header parsing variable closure-local.

### DNA (`workflows/align_dna.nf`)

- Added the same collision-safe basename canonicalization and metadata-driven
  read-group provenance used by ATAC.
- Replaced the synthesized `ID`/`SM`/`PL` read group with the complete metadata
  `rg_string` supplied to minimap2 `-R`.
- Replaced the broken library filter, which compared a set of library-name
  strings with an entire FASTQ tuple, with the same keyed channel matching used
  by RNA and ATAC.
- Added fail-closed FASTQ discovery and missing-metadata errors.
- Made variant calling accept both gzip-compressed and plain FASTA references,
  matching the indexing workflow's stated input support.
- Restored the original 96-hour variant-calling time limit. The unexplained
  local 36-hour value was not paired with a timeout retry condition.

### Entrypoint, reference builder, and configuration

- Declared `rg_metadata` and `rna_geometry` in `nextflow.config`; the default RNA
  geometry is `long-r2`.
- Replaced noninteractive `conda init` with deterministic shell-hook activation
  of the `align_pipelines` environment.
- Made the main entrypoint require provenance metadata, fail on missing library
  and species files, ignore blank/comment-only library-list rows, and keep the
  species-file staging process explicit.
- Made reference input discovery fail closed, validated positive thread and
  memory settings, created the STAR index directory before use, and connected
  requested CPUs to STAR and minimap2 index construction.
- Kept `nextflow.config.slurm` as an overlay for the base configuration rather
  than a replacement.
- Removed SLURM `withName` selectors for nonexistent process names
  (`STAR_align`, `minimap2_align`, and `freebayes`). Process directives and the
  generic cluster defaults now govern real tasks.
- Removed the stale global `/mnt/beegfs/home/$USER/nextflow_work` setting so the
  launch context or orchestrator can select the run-specific work directory.
- Kept resume enabled and disabled automatic cleanup; enabling both previously
  deleted the cache required by resume.
- Retained the SLURM executor cap of 10 active tasks and submission limit of 10
  jobs per minute. These settings throttle Nextflow-submitted tasks, not SLURM
  array concurrency.

### ATAC preprocessing helper

- Restored protection against using an input FASTQ as an output destination in
  `src/atac_fq_preprocess.cpp`. Opening such a path with gzip write mode would
  truncate the source before barcode preprocessing began.
- Replaced the removed `std::filesystem::equivalent` implementation with a
  C++11-compatible canonical-path comparison, preserving the Makefile's C++11
  build target while also detecting symlinks and not-yet-created outputs.
- Added explicit output-directory and gzip-output-open validation so filesystem
  errors fail before read processing.

### Documentation and examples

- Replaced quick-start commands that referenced absent wrapper scripts, SLURM
  templates, and example files with direct commands for the files present in
  this repository.
- Documented the required ten-column `rg_metadata` contract, exact `fnbase`
  lookup behavior, `__RunNNN` collision suffixes, BAM read-group provenance,
  and both supported RNA geometries.
- Updated every mapping parameter template to include `rg_metadata`, added RNA
  geometry where relevant, and exposed reference-builder thread control.
- Documented `nextflow.config.slurm` as a soft overlay applied with
  `nextflow -c ... run`, preserving the base parameter configuration.

## 2026-08-23 — RNA/ATAC consolidation integration

- Added the ten-column source-unit metadata interface and fail-closed read-group
  lookup to RNA and ATAC.
- Added collision-safe `__RunNNN` awareness and canonical handling of legacy and
  bcl-convert filename annotations.
- Added RNA `long-r2` and `pe150` geometry selection.
- Propagated full FASTQ basenames through normal and species-demultiplexed RNA
  channels so STAR can associate ordered inputs with ordered read groups.
- Propagated full FASTQ basenames and read groups through ATAC preprocessing and
  chunking.

## 2025-10-21 — RNA filtering and unmapped-output policy

- Removed the explicit fixed `--outFilterScoreMin 30` setting.
- Explicitly set both `--outFilterScoreMinOverLread` and
  `--outFilterMatchNminOverLread` to 0.33. This is more permissive than STAR's
  normalized defaults; no unmeasured mapping-rate improvement is claimed here.
- Enabled `--outReadsUnmapped Fastx` and published stable
  `Unmapped.out.mate1`/`Unmapped.out.mate2` names.
- The first historical implementation declared names that did not match STAR's
  prefixed files; the subsequent version added renaming. The final workflow
  performs those renames conditionally and declares the outputs optional.

## 2025-10-13 — RNA adapter policy

- Replaced manual 3-prime poly-A clipping (`--clip3pAdapterSeq polyA` and
  `--clip3pAdapterMMp 0.1`) with STAR `--clipAdapterType CellRanger4`.
- Removed `--outSAMunmapped Within`; unmapped reads are no longer retained in
  `gex.bam`. Later changes added separate optional FASTQ export instead.

## Retired historical files

The following files are intentionally excluded from the finalized set because
their history is captured above and their intended behavior is incorporated in
the active workflow:

- `align_rna.nf.backup_20251013_144520`
- `align_rna_backup_2025_1021.nf`
- `align_rna.nf_backup2`
- `CHANGELOG_align_rna.md`
