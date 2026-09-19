# Alignment pipeline changelog

This changelog records local workflow, helper-program, configuration,
documentation, and packaging changes relative to the supplied original
`nkschaefer/align_pipelines` repository copies. It replaces the former RNA-only
changelog and the historical `align_rna.nf` backup files.

## 2026-09-05 - One-pass RNA BAM evidence profiler and analysis layer

- Removed the pilot, reader benchmark, continuation-audit, and phase-release
  machinery. One invocation now submits every selected library in a single
  array capped at three concurrent readers, followed automatically by BAM-free
  gather and the existing downstream dependency chain. Failed-run recovery
  first proves prior jobs terminal, then moves partial evidence and gathered
  outputs into a timestamped archive before regenerating the complete array.
- Repaired the first real-cluster reconciliation failure. The initial release
  incorrectly treated `NH=1` as the boundary of ordinary STARsolo evidence and
  therefore omitted `NH>1` reads whose annotation union resolves to one gene.
  Ordinary counted reads and matrix molecules are now NH-unrestricted. The
  separately reported `nh_gt1_unique_gene_countedU_reads` value is explicitly
  an ordinary singleton-gene subset; exact multi-gene EM evidence remains
  unavailable from the standard uppercase `GX`/`UB` BAM tags.
- A nonblocking per-library lock prevents concurrent publishers, while normal
  resume reuses only fully validated same-signature outputs.
- Corrected htslib 1.20 header access, coordinate-sorted the compiled test
  fixture before indexing, and replaced the unsupported `samtools idxstats -X`
  preflight with a `samtools view -c -X` region probe that opens the exact
  manifest-declared BAM index. Documented the required `genomics-base/latest`
  test environment and added direct custom-index regressions.
- Added the C++17 `rna_bam_evidence` target. It decodes each coordinate-sorted
  BAM exactly once through htslib, uses threaded BGZF decompression, validates
  `RG CB CR GX GN UB UR NH AS nM`, and emits barcode, RG, numeric contig,
  correction, and aggregate molecule evidence without per-read output.
- Split STARsolo semantics explicitly. Ordinary counted reads are one logical
  `(RG,QNAME)` with a declared RG, valid `CB`, and a single feature-roster `GX`;
  NH=1 uses one fragment representative, while NH>1 merges primary/secondary
  records so a secondary-only `GX` is retained. They do not require an accepted
  `UB`. Ordinary molecules are distinct `(CB,GX,accepted-corrected-UB)` tuples
  across all mapped records and are also NH-unrestricted. The
  manifest and compiled process are bound to `GeneFull_Ex50pAS`,
  `MultiGeneUMI_CR`, `1MM_CR`, and `EM` before the BAM is opened. The
  `nh_gt1_unique_gene_countedU_reads` field is an ordinary unique-gene subset,
  not an EM metric; exact multi-gene EM evidence is explicitly unavailable. The
  compiled process
  merge-checks every raw and filtered MatrixMarket coordinate and count, while
  the wrapper checks Summary and filtered-cell equalities before publication.
- Replaced string-tree hot-loop lookups with reserved numeric interners and
  open-addressed tables, reused tag decoding, and replaced contig maps with a
  dense RG-by-contig vector. Each record's numeric contribution, including its
  CIGAR walk, is computed once and reused across library, RG, barcode, and
  contig accumulators. Molecules, barcode metrics, barcode-RG metrics,
  correction aggregates, interners, rehash peaks, and output sorting share a
  conservative process admission budget. A whole-process `RLIMIT_AS` guard,
  runtime/allocator reserve, complete MaxRSS, component peaks, and relevant
  cardinalities are recorded. Hash-bin aggregation and high-cardinality sorts
  reuse compacted storage in place.
- Reworked `run_rna_bam_evidence.py` validation and BAM-free gather as sorted
  adjacent-key and merge-style streams. Resume now revalidates every source and
  gathered output for readability, schema, sort order, row count, byte size,
  and SHA-256 before reporting success. Only filtered-cell vectors and small
  RG/library summaries remain in memory.
- Bound every reader audit to the Slurm allocation and array task that actually
  produced it. Failed-run reset refuses active or unprovable prior reader jobs,
  and gather continues to validate every manifest row, audit, product schema,
  row count, and checksum before publication.
- Removed `align_pipelines/bjp` and Nextflow from BAM-reader job module loads;
  the jobs use only miniforge, htslib 1.20, and samtools 1.20 and invoke the
  frozen profiler by absolute path. This avoids the package module's unrelated
  conda `libcurl`/toolchain contamination.
- Restored STAR diagnostic promotion as a lightweight independent step that
  publishes `STAR_Log.out`, `STAR_Log.final.out`, and `STAR_SJ.out.tab.gz` from
  an unambiguous retained task matched by Summary SHA-256. It opens no BAM,
  invokes no retired AWK scanner, and never deletes the work directory. Reuse
  now independently verifies a manifest-bound audit/marker and all output
  hashes; unproven preexisting files or a completed Slurm ID are insufficient.
- Added explicit within-RG, cross-RG, and global raw-to-corrected barcode
  conflict fields. Exact RG trim joins reject only within-RG ambiguity and do
  not discard a valid mapping because another RG differs.
- Added validated contig and feature/GX classification. Missing classification
  produces blank mitochondrial/rRNA values and explicit `unavailable` status,
  never numeric zero. The orchestrator requires either a manifest or a clear
  no-classification acknowledgment before any BAM pass.
- Reworked gained-cell analysis into independent disk-backed library shards.
  It details all shared/gained/lost cells, samples only evidence-bearing
  background barcodes, and computes exact bounded histograms from the complete
  evidence-bearing population. Depth history uses fixed-roster numeric vectors
  and keeps observed endpoints, RG-prefix threshold proxies, and final-BAM
  random-thinning proxies separately labelled.
- Added true category histogram and quantile tables plus separate cell-count,
  fixed-roster UMI/read, source-yield, and six-panel library figures. All bins,
  quantiles, and trajectories are exported as TSV instead of substituting bars
  of category medians.
- Expanded tests for the exact focus set, a source-derived post-STARsolo SAM
  and MatrixMarket contract with NH>1 ordinary inclusion and EM unavailability,
  one-mismatch raw `UR`
  collapse to a corrected `UB`, and `MultiGeneUMI_CR` rejection, real-htslib
  compilation when available, one-thread versus multithread invariance, both conflict scopes,
  unavailable classes, corrupt and missing resume products, 120,000-row
  streaming validation, diagnostic promotion, and direct bounded full-array
  generation.
  Real compiled tests skip explicitly when this environment lacks htslib or
  samtools and never use fake declarations.

## 2026-09-04 — Validated parallel backfill for historical RNA mappings

- Clarified the retention boundary. RNA `gex.bam`, its index, raw and filtered
  matrices, `Barcodes.stats`, `Features.stats`, `Summary.csv`, and
  `UMIperCellSorted.txt` were already published to each library's final
  `mapping_output/<library>/` directory; they are not copied again.
- Added `--rna3-cell-reads-backfill-from-bam` for completed 3' RNA mappings
  made before native `CellReads.stats` was enabled. The mode streams each
  published `gex.bam` through `samtools` and `awk` without creating an
  intermediate SAM or BAM copy, and writes
  `CellReads.countedU.from_bam.tsv.gz` beside the BAM.
- BAM-derived counts use one primary alignment per read (`-F 2304`) and require
  a called `CB` plus nonmissing `GX` and `UB` tags. A result is published only
  if its barcode set, cell count, summed unique gene-assigned reads, and
  STAR-style median exactly reproduce the library's `filtered/barcodes.tsv.gz`
  and `Summary.csv`. An audit JSON and completion marker record the successful
  validation; the file is deliberately not named native `CellReads.stats`
  because its other diagnostic columns cannot be reconstructed from the BAM.
- Backfill is scheduled as a smallest-BAM pilot followed by a per-library
  SLURM array. The default task uses four CPUs, the array is capped at eight
  simultaneous libraries to bound shared-filesystem traffic, and the existing
  global array throttle can impose a stricter cap. One failed pilot therefore
  blocks the remaining scans instead of wasting work across every library.
- The same worker identifies the successful historical Nextflow task by an
  exact `Summary.csv` checksum and atomically promotes the previously stranded
  `STAR_Log.out`, `STAR_Log.final.out`, and `STAR_SJ.out.tab.gz` into that
  library's final mapping directory. Ambiguous or unverified work-task matches
  are fatal.
- The mapping plotter accepts either future native `CellReads.stats.gz` or the
  explicitly named and validated BAM-derived table. Mapping plots depend on
  completion of the backfill branch, and final validation requires its table,
  audit record, completion marker, and promoted STAR diagnostics.
- Resume planning now treats a newly introduced job label as a new upstream
  branch and automatically invalidates its prior downstream plot/validation
  submissions. Previously completed trimming and mapping job IDs remain
  reusable.

## 2026-09-04 — Durable mapping outputs and manual work-cache cleanup

- Established that a completed mapping run must remain scientifically usable
  after its Nextflow `work/` directory is manually deleted. The orchestrator
  never deletes a work directory automatically.
- Corrected the cell-read-statistics implementation. STAR defaults
  `--soloCellReadStats` to `None`, so merely declaring `CellReads.stats` as a
  Nextflow output could not create the file. RNA mapping now explicitly passes
  `--soloCellReadStats Standard`, compresses the resulting table, and publishes
  `CellReads.stats.gz` within every library directory.
- RNA mapping now also publishes the small STAR products that are useful for
  audit and reanalysis but were previously stranded in `work/`:
  `STAR_Log.out`, `STAR_Log.final.out`, and `STAR_SJ.out.tab.gz`.
- Expanded final validation to require the BAM and index, complete raw and
  filtered STARsolo matrices, aggregate statistics, per-cell read statistics,
  STAR logs, splice-junction table, source/read-group manifests, exact
  Nextflow parameters and configuration, execution report/trace/timeline, and
  generated mapping submission script. ATAC validation now also requires the
  published name-sorted BAM.
- Successful validation writes
  `validation/WORK_CLEANUP_READY.ok` and
  `validation/work_cleanup_targets.tsv`. These identify the exact work-cache
  directories that may then be deleted manually and explicitly warn that doing
  so removes Nextflow resume capability while leaving published results intact.
- Reference-index copies, input symlinks, decompressed whitelist copies,
  temporary sort data, Nextflow task wrappers, split FASTQs, per-source BAMs,
  pre-mark-duplicate BAMs, and other reconstructable intermediates remain
  work-cache products and are not copied into permanent output directories.
- Historical RNA mappings that did not request
  `--soloCellReadStats Standard` never created `CellReads.stats`; keeping their
  work directories cannot recover a file that STAR did not generate. Such runs
  require a validated BAM-derived backfill or mapping regeneration for an exact
  cell-level read distribution.

## 2026-09-04 — Cell-level RNA read distributions

- Added `reads_per_cell_distribution.png`, a cell-level view of the exact
  unique GeneFull_Ex50pAS-assigned read counts used by STARsolo for its
  `Median Reads per Cell` summary statistic. The figure combines a histogram
  and KDE across called cells with a 50,000-read reference line.
- Added a second panel showing the 10th percentile, interquartile range,
  median, and 90th percentile for every library. When a prior mapping root is
  supplied, this panel instead shows the distribution of read-count gains for
  the same library/barcode pairs in both snapshots.
- Added `reads_per_cell.tsv.gz`, preserving every plotted library, called-cell
  barcode, snapshot label, and read count for downstream analysis.
- Read distributions are taken from the `countedU` column of STARsolo's
  `CellReads.stats`, joined to `filtered/barcodes.tsv[.gz]`. Before plotting,
  every library must exactly reproduce the cell count, unique-read total, and
  median recorded in its `Summary.csv`; mismatches are fatal.
- Updated `workflows/align_rna.nf` to publish `CellReads.stats.gz` with each
  library for future mappings. The plotter may also recover a legacy copy from
  a retained Nextflow work directory, but only when that historical STAR run
  explicitly requested cell-read statistics and therefore created the file.
- The orchestrator now requests and verifies the cell-read TSV and PNG in the
  normal RNA mapping plot job. Added optional
  `--rna3-mapping-baseline-root` and
  `--rna3-mapping-baseline-work-root` controls for old-to-new cell-level
  comparisons.
- Histograms always use all called cells. KDE rendering uses a deterministic
  maximum of 100,000 cells per snapshot, and displayed axes stop at the 99.5th
  percentile to keep long tails from obscuring the main distribution.

## 2026-09-04 — Longitudinal RNA mapping deltas

- Added optional old-to-new comparison support to
  `scripts/mapping/plot_mapping_stats_V5.py`. A prior `processedstats.tsv` can
  now be supplied with `--baseline-stats`; `--current-stats` also permits a
  comparison to be regenerated directly from two saved statistics tables.
- Added `mapping_delta_dashboard.png`, showing baseline plus added median
  reads per cell, the change in estimated cell count, baseline versus current
  sequencing saturation, and changes in median UMIs and genes per cell.
- Added `mapping_stats_deltas.tsv`, an exact long-form table containing the
  baseline, current, absolute delta, and percent delta for every shared numeric
  metric and library.
- Added the orchestrator options `--rna3-mapping-baseline-stats`,
  `--rna3-mapping-baseline-label`, and `--rna3-mapping-current-label`. These
  pass the comparison into the normal RNA mapping plot job; no standalone or
  manually sequenced plotting step is required.
- Treat reporting options as warning-only runtime configuration during
  `--resume`, so a baseline can be added to an existing mapping run without
  invalidating completed trimming or mapping work.
- Require the mapping plot job to verify every requested PNG and TSV before it
  writes `MAPPING_PLOTS_COMPLETE.ok`. Plot-generation exceptions now cause a
  nonzero exit instead of being silently converted into a successful job.
- Corrected existing panels that were labelled as medians but read the
  corresponding mean columns. Cell-quality, complexity, and the reads-per-cell
  annotation in the saturation panel now use the actual median columns.
- Sequencing saturation remains a non-additive rate and is therefore compared
  between snapshots rather than divided into stacked run components. Source
  read-group contributions require a separate BAM-derived read-count product;
  cumulative cell-count convergence requires cell calling at each cumulative
  sequencing depth.

## 2026-09-01 — Calibrated per-library RNA memory

- Replaced the fixed RNA `map_rna` SLURM allocation with an input-sized request
  based on the combined compressed size of each library's trimmed R1 and R2
  FASTQs. The calibration used 36 completed TET RNA libraries and four
  independently observed 80 GiB OOM libraries.
- The completed-library fit was `peak RSS GiB = 1.459 + 0.22318 × trimmed GiB`
  (`R² = 0.9963`; residual standard error 0.77 GiB; largest leave-one-library-out
  underprediction 2.62 GiB). The deployed rule deliberately uses the higher
  empirical envelope `4 + 0.224 × trimmed GiB`, adds 25% operational headroom,
  and rounds upward to an 8 GiB scheduling bucket.
- `--memgb` remains a user-controlled minimum allocation. With the established
  80 GiB floor, the four OOM libraries are initially assigned 112 GiB (library
  38), 112 GiB (library 2), 120 GiB (library 36), and 136 GiB (library 1).
- An actual OOM retry adds 32 GiB to the library-specific request for each
  subsequent attempt. Non-memory failures still terminate immediately.
- Removed fixed memory declarations from the generated per-run
  `nextflow.config`; those declarations overrode the workflow's dynamic
  resource directive and would otherwise silently disable the estimator.
- Kept STAR's established `--limitBAMsortRAM` calculation based on `--memgb`.
  This preserves the mapping command and Nextflow cache identity for completed
  libraries, and it matches the exact STAR behavior used to calibrate the
  allocation model. Only the SLURM allocation changes.

## 2026-08-30 — Restore historical 3′ TSO matching

- Reverted the 3′ RNA front-adapter definition in
  `scripts/mapping/integrated_RNA_trim_map_pipeline_V10.py` from the newly
  anchored `TSO=^AAGCAGTGGTATCAACGCAGAGTACATGGG` form to the established V9
  form, `^TSO=AAGCAGTGGTATCAACGCAGAGTACATGGG`.
- In Cutadapt named-adapter syntax, the historical leading caret is part of the
  adapter name rather than an anchor on the sequence. This intentionally
  restores regular 5′ matching, including partial 5′ TSO occurrences permitted
  by the existing minimum-overlap setting.
- The anchored form was an unvalidated trimming-behavior change and was not
  required for barcode-linked trimming QC. The aggregator already canonicalizes
  the historical `^TSO` adapter name correctly.
- No 5′ RNA adapter definitions, other trimming adapters, mapping behavior, QC,
  plotting, orchestration, or workflow files changed in this correction.

## 2026-08-30 — Mapping control-layer migration

- Moved the complete eight-file 10X RNA/ATAC mapping control layer from the
  CellBouncer source/deployment layout into `scripts/mapping` in this
  repository. The installed layout is now `bin/mapping` below the deployed
  align_pipelines package.
- Kept the orchestrator and its seven runtime helpers together so sibling
  driver, aggregation, QC, and plotting resources continue to resolve without
  additional command-line configuration.
- Changed the orchestrator and RNA/ATAC drivers to derive `align_rna.nf` and
  `align_atac.nf` from the containing repository/package root. This supports
  both `scripts/mapping` in the source checkout and `bin/mapping` in the
  installed package while eliminating hard-coded ownership by the CellBouncer
  deployment tree.
- Changed standalone ATAC QC job generation to embed the collector's actual
  source or installed path instead of the former copied path in the production
  mapping-output directory.
- Retired the duplicate `mapping/NextflowConfigs` deployment model. The
  repository's top-level `.nf`/configuration files and `workflows/` directory
  are the only authoritative workflow copies.
- Updated the Nextflow manifest homepage to the
  `Bryanjpavlovic/align_pipelines_TET` fork and retained original and fork
  authorship.
- This migration changes code ownership and resource discovery; it does not
  alter the established trimming, alignment, barcode aggregation, QC, plotting,
  resume, node-placement, or array-throttling behavior of the migrated scripts.

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
