# STARsolo semantic fixture contract

This is a source-modelled post-STARsolo SAM fixture. It is not represented as
the output of a STAR invocation in this development environment. Its semantics
come from STAR 2.7.11b, the version used by the mapping workflow, and the exact
workflow settings `GeneFull_Ex50pAS`, `MultiGeneUMI_CR`, `1MM_CR`, and `EM`.
The fixture is consumed through a BAM produced by real `samtools`, and the
compiled tests require a profiler linked to real htslib.

## Ordinary unique-gene reads are not restricted to NH=1

STAR separates genomic mapping multiplicity from gene-assignment
multiplicity. In STAR 2.7.11b `source/SoloReadFeature_record.cpp`, `nTr > 1`
sets the genomic-multimapper flag, but a read enters the ordinary unique-gene
path whenever its read-wide feature set has exactly one gene. Consequently an
`NH:i:2` read can contribute to ordinary `matrix.mtx`.

Uppercase `GX` is alignment-specific. In
`source/ReadAlign_alignBAM.cpp`, STAR writes a valid uppercase `GX` only when
the read-wide feature set has one gene and that particular alignment overlaps
exactly one gene. The primary record is not guaranteed to be the record that
exposes the accepted gene. Corrected `CB` and `UB` are read-level values added
to every BAM record during coordinate sorting by
`source/SoloFeature_addBAMtags.cpp` and
`source/BAMbinSortByCoordinate.cpp`.

STAR's own `extras/scripts/soloCountMatrixFromBAM.awk` reconstructs the
ordinary matrix by examining all BAM records, retaining records with valid
`GX`, `CB`, and `UB`, and counting distinct UBs within each cell/gene
coordinate. It does not apply a primary-record or `NH==1` gate. Repeated
alignments and mates collapse through the `(CB,GX,UB)` molecule identity.

The fixture exercises those rules directly:

- `q08` is one `NH:i:2` read with complete `HI:i:1` and `HI:i:2` records. Its
  primary alignment has `GX:Z:-`; its secondary alignment exposes `GX:Z:G1`.
  It contributes one ordinary read and one G1 molecule.
- The supplementary `q01` record has the same read name and molecule tags as
  the primary `q01` record. It contributes neither an additional read nor an
  additional molecule.
- `q23` is a paired `NH:i:1` read with first- and second-mate records. Both
  records expose the same G1 molecule, so the pair contributes one read and
  one molecule.

## Read counts and molecule counts have different gates

STAR records `CellReads.stats` `countedU` when a corrected cell barcode is
assigned to a unique gene. This happens in
`source/SoloReadFeature_inputRecords.cpp` before UMI correction and
`MultiGeneUMI_CR` filtering. A unique-gene read can therefore contribute to
`countedU` even when its eventual `UB` is `-` and it contributes no ordinary
matrix molecule.

Records `q19` and `q20` model equal support for one raw UMI at G1 and G2.
Both remain unique-gene reads for `countedU`, but their rejected `UB:Z:-`
values exclude them from `matrix.mtx`. Record `q18` models `1MM_CR`: its raw
`UR` differs from `q01` by one base, while its corrected `UB` equals `q01`, so
it adds a read but not a molecule.

Records `q16`, `q17`, and `q21` are non-feature reads with `GX:Z:-`. They
exist only to exercise within-RG and cross-RG raw-to-corrected barcode
conflicts. Their valid read-level CB/UB tags do not make them gene-assigned
reads or matrix molecules.

## True multi-gene EM evidence is separate

STAR's ordinary `matrix.mtx` contains unique-gene UMI counts. A
`UniqueAndMult-EM.mtx` contains ordinary counts plus separately allocated
multi-gene UMI mass. This distinction is documented in `docs/STARsolo.md` and
implemented in `source/SoloFeature_collapseUMIall.cpp`.

`q22` is a true multi-gene, two-alignment read whose read-wide candidate gene
set is G2/G3. Its uppercase `GX` is `-` on both BAM records, as STAR emits for
a non-unique gene assignment, and it contributes nothing to ordinary
`matrix.mtx` or `countedU`. In barcode C, G2 and G3 both have zero ordinary
unique-gene molecules. STAR's EM initialization and update therefore retain
an exact equal allocation: 0.5 to G2 and 0.5 to G3. This adds one unit only to
the checked-in `UniqueAndMult-EM.mtx`.

## Exact checked-in totals

The ordinary matrix contains nine molecules in five nonzero cell/gene
coordinates:

- barcode A: four molecules across G1, G2, and G3;
- barcode C: five molecules across G1 and G4;
- total: nine molecules.

Unique-gene read counts are eleven for barcode A and five for barcode C,
including the two `UB:Z:-` reads and counting each multi-record read only once.
The checked-in Summary therefore contains sixteen unique reads, median reads
per cell eleven, median UMIs per cell four, and median genes per cell three.
STAR sorts UMI counts in descending order before selecting index
`n_cells / 2`, while its read and gene vectors are sorted in ascending order;
the intentionally even two-cell fixture records that exact behavior.
The EM matrix contains the nine ordinary molecules plus the one unit from
`q22`, for a total mass of ten across seven nonzero coordinates.

`STARsolo_semantic_expectations.tsv` has one row per SAM record. Its
`ordinary_candidate_read` field is a one-time contribution per input read,
placed on the record that establishes the assignment; it must sum to sixteen.
`ordinary_matrix_tag_record` records whether that BAM record has valid
matrix-bearing `CB/GX/UB` tags and may be duplicated across alignments or
mates. `ordinary_new_molecule` marks first observation in fixture order and
must sum to nine. The two EM mass columns sum to one.
