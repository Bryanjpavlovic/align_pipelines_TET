#!/usr/bin/env python3
"""Focused synthetic contract tests for the RNA BAM evidence pipeline."""

from __future__ import annotations

import csv
import gzip
import gc
import os
import shutil
import shlex
import subprocess
import sys
import tempfile
import tracemalloc
import unittest
import argparse
import json
from unittest import mock
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MAPPING = ROOT / "scripts" / "mapping"
sys.path.insert(0, str(MAPPING))

import run_rna_bam_evidence as runner
import trim_barcode_aggregator as trim_aggregator
import analyze_gained_cells
import analyze_depth_history
import export_repooling_evidence
import join_trim_cell_metrics
import plot_rna_evidence
import promote_rna_star_diagnostics
import rna_evidence_common


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stable_hash(qname: str, seed: int = runner.DEFAULT_HASH_SEED) -> int:
    value = 14695981039346656037 ^ seed
    for byte in qname.encode("utf-8"):
        value ^= byte
        value = (value * 1099511628211) & ((1 << 64) - 1)
    return value


def write_fixture_scientific_config(
    path: Path,
    fields: list[str],
    values: dict[str, object],
    profiler: Path,
) -> dict[str, object]:
    row = {field: str(values[field]) for field in fields}
    source_path = Path(row["source_order"])
    source_contents = source_path.read_text(encoding="utf-8")
    class_path = Path(row.get("class_manifest", ""))
    if class_path.is_file():
        classification = {
            "status": "manifest", "path": str(class_path.resolve()),
            "sha256": runner.sha256(class_path),
            "contents": class_path.read_text(encoding="utf-8"),
        }
    else:
        classification = {
            "status": "explicit_unavailable", "path": "",
            "sha256": "", "contents": "",
        }
    config: dict[str, object] = {
        "schema_version": 1,
        "manifest_fields": fields,
        "manifest_rows": [row],
        "source_order": {
            "path": str(source_path.resolve()), "provenance": "synthetic_fixture",
            "values": ["BP1", "BP2"], "sha256": runner.sha256(source_path),
            "contents": source_contents,
        },
        "biological_classification": classification,
        "hash": {
            "algorithm": runner.HASH_ALGORITHM,
            "seed": runner.DEFAULT_HASH_SEED, "bins": 100,
        },
        "profiler_executable": {
            "path": str(profiler.resolve()),
            "sha256": runner.sha256(profiler),
            "bytes": profiler.stat().st_size,
        },
        "starsolo": {
            "feature": runner.STARSOLO_FEATURE,
            "umi_filtering": runner.STARSOLO_UMI_FILTERING,
            "umi_dedup": runner.STARSOLO_UMI_DEDUP,
            "multimappers": runner.STARSOLO_MULTIMAPPERS,
            "ordinary_countedU_read_definition":
                runner.ORDINARY_COUNTEDU_READ_DEFINITION,
            "ordinary_molecule_definition": runner.ORDINARY_MOLECULE_DEFINITION,
            "multimapper_definition": runner.MULTIMAPPER_DEFINITION,
            "summary_unique_read_metric": runner.SUMMARY_UNIQUE_READ_METRIC,
            "nh_gt1_unique_gene_countedU_definition":
                runner.NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION,
            "starsolo_EM_evidence_availability":
                runner.STARSOLO_EM_EVIDENCE_AVAILABILITY,
        },
    }
    config["configuration_hash"] = runner.canonical_payload_hash(config)
    path.write_text(json.dumps(config, sort_keys=True) + "\n", encoding="utf-8")
    return config


class PythonContractTest(unittest.TestCase):
    def test_bam_preflight_uses_samtools_120_custom_index_interface(self) -> None:
        paths = {
            "bam": Path("/fixture/library/gex.bam"),
            "bam_index": Path("/fixture/indexes/gex.bam.bai"),
        }
        completed = (
            subprocess.CompletedProcess([], 0, "", ""),
            subprocess.CompletedProcess(
                [], 0,
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@RG\tID:RG1\tSM:Fixture\n",
                "",
            ),
            subprocess.CompletedProcess([], 0, "0\n", ""),
        )
        with mock.patch.object(
            runner, "run_checked", side_effect=completed
        ) as checked:
            result = runner.bam_preflight(paths)

        self.assertEqual(
            checked.call_args_list[2].args[0],
            [
                "samtools", "view", "-c", "-X",
                "/fixture/library/gex.bam",
                "/fixture/indexes/gex.bam.bai",
                "chr1:1-1",
            ],
        )
        self.assertEqual(result["indexed_contigs"], 1)
        self.assertEqual(result["index_probe_region"], "chr1:1-1")
        self.assertEqual(result["index_probe_records"], 0)

    def test_resume_signature_excludes_only_operational_cpu_and_memory(self) -> None:
        with tempfile.TemporaryDirectory(prefix="evidence_signature_") as temporary:
            executable = Path(temporary) / "profiler"
            executable.write_text("fixture\n", encoding="utf-8")
            row = {
                "library": "L1",
                "starsolo_feature": runner.STARSOLO_FEATURE,
                "starsolo_umi_filtering": runner.STARSOLO_UMI_FILTERING,
                "starsolo_umi_dedup": runner.STARSOLO_UMI_DEDUP,
                "starsolo_multimappers": runner.STARSOLO_MULTIMAPPERS,
            }
            config = {
                "configuration_hash": "fixture-config",
                "hash": {
                    "algorithm": runner.HASH_ALGORITHM,
                    "seed": runner.DEFAULT_HASH_SEED, "bins": 100,
                },
                "profiler_executable": {
                    "path": str(executable.resolve()),
                    "sha256": runner.sha256(executable),
                    "bytes": executable.stat().st_size,
                },
                "starsolo": {
                    "feature": runner.STARSOLO_FEATURE,
                    "umi_filtering": runner.STARSOLO_UMI_FILTERING,
                    "umi_dedup": runner.STARSOLO_UMI_DEDUP,
                    "multimappers": runner.STARSOLO_MULTIMAPPERS,
                    "ordinary_countedU_read_definition":
                        runner.ORDINARY_COUNTEDU_READ_DEFINITION,
                    "ordinary_molecule_definition": runner.ORDINARY_MOLECULE_DEFINITION,
                    "multimapper_definition": runner.MULTIMAPPER_DEFINITION,
                    "summary_unique_read_metric": runner.SUMMARY_UNIQUE_READ_METRIC,
                    "nh_gt1_unique_gene_countedU_definition":
                        runner.NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION,
                    "starsolo_EM_evidence_availability":
                        runner.STARSOLO_EM_EVIDENCE_AVAILABILITY,
                },
            }
            first = runner.build_signature(
                row, "fixture-v1", executable,
                argparse.Namespace(threads=1, max_memory_gb=32), config,
            )
            second = runner.build_signature(
                row, "fixture-v1", executable,
                argparse.Namespace(threads=8, max_memory_gb=96), config,
            )
            self.assertEqual(first, second)

    def test_matrix_contract(self) -> None:
        stats, umis, genes = runner.matrix_shape_and_counts(
            HERE / "matrix.mtx", 2, keep_per_barcode=True
        )
        self.assertEqual(stats, {"features": 4, "barcodes": 2, "nnz": 5, "molecules": 9})
        self.assertEqual(umis, [4, 5])
        self.assertEqual(genes, [3, 2])

    def test_starsolo_ordinary_and_multimapper_semantics_are_distinct(self) -> None:
        expectation_rows = read_tsv(
            HERE / "STARsolo_semantic_expectations.tsv"
        )
        feature_order = [
            line.split("\t", 1)[0]
            for line in (HERE / "features.tsv").read_text(encoding="utf-8").splitlines()
            if line
        ]
        barcode_order = [
            line.split("\t", 1)[0]
            for line in (HERE / "barcodes.tsv").read_text(encoding="utf-8").splitlines()
            if line
        ]
        feature_index = {value: index + 1 for index, value in enumerate(feature_order)}
        barcode_index = {value: index + 1 for index, value in enumerate(barcode_order)}
        rg_ids: set[str] = set()
        alignments: list[dict[str, object]] = []
        with (HERE / "fixture.sam").open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("@RG"):
                    rg_ids.add(next(
                        field[3:] for field in line.rstrip("\n").split("\t")
                        if field.startswith("ID:")
                    ))
                elif not line.startswith("@"):
                    fields = line.rstrip("\n").split("\t")
                    tags = {
                        field[:2]: field[5:]
                        for field in fields[11:]
                        if len(field) >= 6 and field[2] == ":"
                    }
                    alignments.append({
                        "qname": fields[0],
                        "flag": int(fields[1]),
                        "tags": tags,
                    })

        self.assertEqual(len(alignments), len(expectation_rows))
        records_by_id: dict[str, dict[str, object]] = {}
        for alignment, expected in zip(alignments, expectation_rows):
            tags = alignment["tags"]
            assert isinstance(tags, dict)
            self.assertEqual(alignment["qname"], expected["qname"])
            self.assertEqual(alignment["flag"], int(expected["flag"]))
            self.assertEqual(int(tags["HI"]), int(expected["HI"]))
            self.assertEqual(int(tags["NH"]), int(expected["NH"]))
            records_by_id[expected["record_id"]] = alignment

        ordinary_molecules: set[tuple[str, str, str]] = set()
        ordinary_entries: dict[tuple[int, int], int] = {}
        ordinary_reads_by_cb = {barcode: 0 for barcode in barcode_order}
        em_entries: dict[tuple[int, int], float] = {}
        for alignment, expected in zip(alignments, expectation_rows):
            flag = int(alignment["flag"])
            tags = alignment["tags"]
            assert isinstance(tags, dict)
            gx = tags.get("GX", "")
            cb = tags.get("CB", "")
            ub = tags.get("UB", "")
            valid_gene = (
                gx not in {"", "-", "0"}
                and ";" not in gx and "," not in gx
                and gx in feature_index
            )
            valid_cb = cb not in {"", "-", "0"} and cb in barcode_index
            valid_ub = ub not in {"", "-", "0"}
            mapped = not flag & 0x4
            matrix_tag_record = int(
                mapped and valid_cb and valid_gene and valid_ub
            )
            self.assertEqual(
                matrix_tag_record,
                int(expected["ordinary_matrix_tag_record"]),
                expected["record_id"],
            )
            new_molecule = 0
            if matrix_tag_record:
                molecule = (cb, gx, ub)
                if molecule not in ordinary_molecules:
                    ordinary_molecules.add(molecule)
                    coordinate = (feature_index[gx], barcode_index[cb])
                    ordinary_entries[coordinate] = ordinary_entries.get(coordinate, 0) + 1
                    new_molecule = 1
            self.assertEqual(
                new_molecule,
                int(expected["ordinary_new_molecule"]),
                expected["record_id"],
            )

        # STARsolo countedU is a logical-read count. NH==1 uses one fragment
        # representative, while NH>1 must consider all mapped alignments
        # because the sole accepted uppercase GX can occur on a secondary.
        grouped: dict[tuple[str, str], list[int]] = {}
        for index, alignment in enumerate(alignments):
            tags = alignment["tags"]
            assert isinstance(tags, dict)
            grouped.setdefault(
                (str(tags.get("RG", "")), str(alignment["qname"])), []
            ).append(index)
        counted_record_indexes: set[int] = set()
        nh_gt1_unique_gene_indexes: set[int] = set()
        for (rg, _qname), indexes in grouped.items():
            nh_values = {
                int(alignments[index]["tags"]["NH"])  # type: ignore[index]
                for index in indexes
            }
            self.assertEqual(len(nh_values), 1)
            nh = next(iter(nh_values))
            eligible = [
                index for index in indexes
                if not int(alignments[index]["flag"]) & (0x4 | 0x800)
            ]
            if nh == 1:
                representatives = []
                for index in eligible:
                    flag = int(alignments[index]["flag"])
                    if flag & 0x100:
                        continue
                    if not flag & 0x1 or flag & 0x40 or (
                        flag & 0x80 and flag & 0x8
                    ):
                        representatives.append(index)
                self.assertLessEqual(len(representatives), 1)
                if not representatives:
                    continue
                representative = representatives[0]
                tags = alignments[representative]["tags"]
                assert isinstance(tags, dict)
                gx = str(tags.get("GX", ""))
                cb = str(tags.get("CB", ""))
                if (
                    rg in rg_ids
                    and cb in barcode_index
                    and gx in feature_index
                ):
                    counted_record_indexes.add(representative)
                    ordinary_reads_by_cb[cb] += 1
                continue

            feature_records = [
                index for index in eligible
                if str(alignments[index]["tags"].get("GX", ""))  # type: ignore[union-attr]
                in feature_index
            ]
            feature_ids = {
                str(alignments[index]["tags"]["GX"])  # type: ignore[index]
                for index in feature_records
            }
            barcode_ids = {
                str(alignments[index]["tags"].get("CB", ""))  # type: ignore[union-attr]
                for index in eligible
                if str(alignments[index]["tags"].get("CB", ""))  # type: ignore[union-attr]
                in barcode_index
            }
            if rg not in rg_ids or len(feature_ids) != 1 or len(barcode_ids) != 1:
                continue
            representative = feature_records[0]
            counted_record_indexes.add(representative)
            nh_gt1_unique_gene_indexes.add(representative)
            ordinary_reads_by_cb[next(iter(barcode_ids))] += 1

        for index, expected in enumerate(expectation_rows):
            self.assertEqual(
                int(index in counted_record_indexes),
                int(expected["ordinary_candidate_read"]),
                expected["record_id"],
            )
            self.assertEqual(
                int(index in nh_gt1_unique_gene_indexes),
                int(expected["nh_gt1_unique_gene_ordinary_read"]),
                expected["record_id"],
            )

        q01 = records_by_id["q01_primary"]["tags"]
        q18 = records_by_id["q18_primary"]["tags"]
        assert isinstance(q01, dict) and isinstance(q18, dict)
        self.assertEqual(
            sum(
                left != right
                for left, right in zip(str(q01["UR"]), str(q18["UR"]))
            ),
            1,
        )
        self.assertNotEqual(q01["UR"], q18["UR"])
        self.assertEqual(q01["UB"], q18["UB"])

        q19 = records_by_id["q19_primary"]["tags"]
        q20 = records_by_id["q20_primary"]["tags"]
        assert isinstance(q19, dict) and isinstance(q20, dict)
        self.assertEqual(q19["UR"], q20["UR"])
        self.assertNotEqual(q19["GX"], q20["GX"])
        self.assertEqual(q19["UB"], "-")
        self.assertEqual(q20["UB"], "-")
        self.assertEqual(
            {
                row["record_id"]
                for row in expectation_rows
                if row["ordinary_candidate_read"] == "1"
                and row["ordinary_matrix_tag_record"] == "0"
            },
            {"q19_primary", "q20_primary"},
        )

        q08_primary = records_by_id["q08_primary"]
        q08_secondary = records_by_id["q08_secondary"]
        self.assertFalse(int(q08_primary["flag"]) & 0x100)
        self.assertTrue(int(q08_secondary["flag"]) & 0x100)
        self.assertEqual(q08_primary["tags"]["GX"], "-")  # type: ignore[index]
        self.assertEqual(q08_secondary["tags"]["GX"], "G1")  # type: ignore[index]

        q01_supplementary = records_by_id["q01_supplementary"]
        self.assertTrue(int(q01_supplementary["flag"]) & 0x800)
        self.assertEqual(
            next(
                row for row in expectation_rows
                if row["record_id"] == "q01_supplementary"
            )["ordinary_new_molecule"],
            "0",
        )
        q23_read1 = records_by_id["q23_read1"]
        q23_read2 = records_by_id["q23_read2"]
        self.assertTrue(int(q23_read1["flag"]) & 0x40)
        self.assertTrue(int(q23_read2["flag"]) & 0x80)
        self.assertEqual(q23_read1["tags"]["UB"], q23_read2["tags"]["UB"])  # type: ignore[index]

        def matrix_entries(path: Path) -> dict[tuple[int, int], float]:
            result: dict[tuple[int, int], float] = {}
            with path.open("r", encoding="utf-8") as handle:
                line = handle.readline()
                self.assertTrue(line.startswith("%%MatrixMarket matrix coordinate"))
                line = handle.readline()
                while line.startswith("%"):
                    line = handle.readline()
                _features, _barcodes, declared = map(int, line.split())
                for line in handle:
                    if not line.strip() or line.startswith("%"):
                        continue
                    feature, barcode, count = line.split()
                    result[(int(feature), int(barcode))] = float(count)
            self.assertEqual(len(result), declared)
            return result

        self.assertEqual(
            matrix_entries(HERE / "matrix.mtx"),
            {key: float(value) for key, value in ordinary_entries.items()},
        )
        em_entries.update(
            {key: float(value) for key, value in ordinary_entries.items()}
        )
        for alignment, expected in zip(alignments, expectation_rows):
            tags = alignment["tags"]
            assert isinstance(tags, dict)
            cb = str(tags.get("CB", ""))
            for feature, field in (
                ("G2", "em_G2_mass"), ("G3", "em_G3_mass")
            ):
                mass = float(expected[field])
                if not mass:
                    continue
                coordinate = (feature_index[feature], barcode_index[cb])
                em_entries[coordinate] = em_entries.get(coordinate, 0.0) + mass
        self.assertEqual(matrix_entries(HERE / "UniqueAndMult-EM.mtx"), em_entries)

        self.assertEqual(
            sum(int(row["ordinary_candidate_read"]) for row in expectation_rows),
            16,
        )
        self.assertEqual(
            sum(int(row["ordinary_new_molecule"]) for row in expectation_rows),
            9,
        )
        self.assertEqual(
            sum(
                int(row["nh_gt1_unique_gene_ordinary_read"])
                for row in expectation_rows
            ),
            1,
        )
        self.assertEqual(
            sum(
                int(row["true_multigene_em_candidate_read"])
                for row in expectation_rows
            ),
            1,
        )
        q08_expectations = [
            row for row in expectation_rows if row["qname"] == "q08"
        ]
        q22_expectations = [
            row for row in expectation_rows if row["qname"] == "q22"
        ]
        self.assertEqual(
            sum(int(row["ordinary_candidate_read"]) for row in q08_expectations),
            1,
        )
        self.assertEqual(
            sum(int(row["true_multigene_em_candidate_read"]) for row in q08_expectations),
            0,
        )
        self.assertEqual(
            sum(int(row["ordinary_candidate_read"]) for row in q22_expectations),
            0,
        )
        self.assertEqual(
            sum(int(row["true_multigene_em_candidate_read"]) for row in q22_expectations),
            1,
        )
        self.assertTrue(
            all(
                alignment["tags"].get("GX") == "-"  # type: ignore[union-attr]
                for alignment in alignments
                if alignment["qname"] == "q22"
            ),
            "standard uppercase BAM tags must not pretend to expose exact EM genes",
        )
        ordinary, _umis, _genes = runner.matrix_shape_and_counts(
            HERE / "matrix.mtx", 2, keep_per_barcode=False
        )
        self.assertEqual(ordinary["molecules"], 9)
        self.assertEqual(sum(em_entries.values()), 10.0)
        summary = runner.summary_values(HERE / "Summary.csv")
        per_barcode_umis = [
            sum(
                count for (feature, column), count in ordinary_entries.items()
                if column == barcode_index[barcode]
            )
            for barcode in barcode_order
        ]
        per_barcode_genes = [
            sum(
                1 for (feature, column) in ordinary_entries
                if column == barcode_index[barcode]
            )
            for barcode in barcode_order
        ]
        self.assertEqual(summary["Estimated Number of Cells"], len(barcode_order))
        self.assertEqual(
            summary["Unique Reads in Cells Mapped to GeneFull_Ex50pAS"],
            sum(ordinary_reads_by_cb.values()),
        )
        self.assertEqual(summary["Median Reads per Cell"], 11)
        self.assertEqual(
            summary["Median UMI per Cell"],
            sorted(per_barcode_umis, reverse=True)[len(per_barcode_umis) // 2],
        )
        self.assertEqual(
            summary["Median GeneFull_Ex50pAS per Cell"],
            sorted(per_barcode_genes)[1],
        )
        self.assertEqual(
            {row["CB"]: int(row["countedU"]) for row in read_tsv(HERE / "CellReads.stats")},
            ordinary_reads_by_cb,
        )

    def test_hash_definition(self) -> None:
        self.assertEqual(stable_hash("q01"), stable_hash("q01"))
        self.assertNotEqual(stable_hash("q01"), stable_hash("q02"))
        self.assertTrue(0 <= ((stable_hash("q01") * 100) >> 64) < 100)

    def test_slurm_execution_identity_is_exact_and_fail_closed(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(
                runner.slurm_execution_identity(), {"kind": "not_slurm"}
            )
        with mock.patch.dict(
            os.environ, {"SLURM_JOB_ID": "101"}, clear=True
        ):
            self.assertEqual(
                runner.slurm_execution_identity(),
                {
                    "kind": "standalone",
                    "SLURM_JOB_ID": "101",
                    "accounting_job_task_id": "101",
                },
            )
        with mock.patch.dict(
            os.environ,
            {
                "SLURM_JOB_ID": "205",
                "SLURM_ARRAY_JOB_ID": "202",
                "SLURM_ARRAY_TASK_ID": "2",
            },
            clear=True,
        ):
            self.assertEqual(
                runner.slurm_execution_identity(),
                {
                    "kind": "array_task",
                    "SLURM_JOB_ID": "205",
                    "SLURM_ARRAY_JOB_ID": "202",
                    "SLURM_ARRAY_TASK_ID": 2,
                    "accounting_job_task_id": "202_2",
                },
            )
        with mock.patch.dict(
            os.environ,
            {"SLURM_JOB_ID": "205", "SLURM_ARRAY_JOB_ID": "202"},
            clear=True,
        ):
            with self.assertRaisesRegex(
                runner.EvidenceError, "must be present together"
            ):
                runner.slurm_execution_identity()

    def test_library_output_lock_rejects_concurrent_writer(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rna_evidence_lock_") as temporary:
            output = Path(temporary) / "bam_evidence"
            with runner.exclusive_output_lock(output):
                with self.assertRaisesRegex(
                    runner.EvidenceError, "holds the library output lock"
                ):
                    with runner.exclusive_output_lock(output):
                        self.fail("a second writer acquired the output lock")

    def test_gather_rg_sort_key_keeps_missing_rg_last(self) -> None:
        declared = {
            "library": "L1", "RG": "RG1", "source_index": "0"
        }
        missing = {
            "library": "L1", "RG": "__MISSING_OR_UNDECLARED__",
            "source_index": "-1",
        }
        self.assertLess(
            runner.product_sort_key("all_libraries_rg_summary.tsv.gz", declared),
            runner.product_sort_key("all_libraries_rg_summary.tsv.gz", missing),
        )

    def test_trim_bridge_avoids_guessing_conflicts(self) -> None:
        with tempfile.TemporaryDirectory(prefix="trim_bridge_fixture_") as temporary:
            bridge = Path(temporary) / "bridge.tsv.gz"
            with gzip.open(bridge, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "CR\tRG\tCB\tread_count\twithin_rg_conflict\t"
                    "cross_rg_conflict\tglobal_conflict\n"
                )
                handle.write("AAAAAAAAAAAAAAAA\tRG1\tAAAAAAAAAAAAAAAA-1\t4\t0\t0\t0\n")
                handle.write("GGGGGGGGGGGGGGGG\tRG1\tAAAAAAAAAAAAAAAA-1\t1\t1\t0\t1\n")
                handle.write("GGGGGGGGGGGGGGGG\tRG1\tCCCCCCCCCCCCCCCC-1\t1\t1\t0\t1\n")
                handle.write("GGGGGGGGGGGGGGGG\tRG2\tAAAAAAAAAAAAAAAA-1\t1\t0\t0\t1\n")
                handle.write("TTTTTTTTTTTTTTTT\tRG1\tAAAAAAAAAAAAAAAA-1\t1\t0\t1\t1\n")
                handle.write("TTTTTTTTTTTTTTTT\tRG2\tCCCCCCCCCCCCCCCC-1\t1\t0\t1\t1\n")
            runner.validate_correction_conflicts(bridge)
            loaded = trim_aggregator.load_profiler_cb_bridge(
                bridge, {"AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"}
            )
            self.assertEqual(loaded.source_kind, "profiler_aggregate_bridge")
            self.assertEqual(loaded.assignments["AAAAAAAAAAAAAAAA"], "AAAAAAAAAAAAAAAA")
            self.assertIn("TTTTTTTTTTTTTTTT", loaded.conflicts)
            corrected, status = trim_aggregator.resolve_barcode(
                "read", "TTTTTTTTTTTTTTTT",
                {"AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"}, loaded,
            )
            self.assertIsNone(corrected)
            self.assertEqual(status, "starsolo_conflict_unresolved")
            rg1 = trim_aggregator.load_profiler_cb_bridge(
                bridge,
                {"AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"},
                rg_filter="RG1",
            )
            rg2 = trim_aggregator.load_profiler_cb_bridge(
                bridge,
                {"AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"},
                rg_filter="RG2",
            )
            self.assertEqual(rg1.assignments["TTTTTTTTTTTTTTTT"], "AAAAAAAAAAAAAAAA")
            self.assertEqual(rg2.assignments["TTTTTTTTTTTTTTTT"], "CCCCCCCCCCCCCCCC")
            self.assertIn("GGGGGGGGGGGGGGGG", rg1.conflicts)
            self.assertEqual(
                rg2.assignments["GGGGGGGGGGGGGGGG"], "AAAAAAAAAAAAAAAA"
            )

    def test_source_dominance_aggregates_lane_rgs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="source_dominance_fixture_") as temporary:
            table = Path(temporary) / "barcode_rg.tsv"
            runner.write_tsv(
                table,
                ["CB", "RG", "source_id", "candidate_countedU_reads"],
                [
                    {"CB": "A-1", "RG": "RG1", "source_id": "BP1", "candidate_countedU_reads": 3},
                    {"CB": "A-1", "RG": "RG2", "source_id": "BP1", "candidate_countedU_reads": 4},
                    {"CB": "A-1", "RG": "RG3", "source_id": "BP2", "candidate_countedU_reads": 5},
                ],
            )
            barcode, metrics = next(
                analyze_gained_cells.iter_source_metrics(table)
            )
            self.assertEqual(barcode, "A-1")
            self.assertEqual(metrics["sources"], 2)
            self.assertEqual(metrics["total"], 12)
            self.assertEqual(metrics["maximum"], 7)
            self.assertEqual(metrics["dominant"], "BP1")

    def test_bam_free_downstream_joins_and_physical_source_prefixes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rna_evidence_downstream_") as temporary:
            root = Path(temporary)
            current_root = root / "current"
            baseline_root = root / "baseline"
            library = "FixtureRNA"
            current = current_root / library
            baseline = baseline_root / library
            for directory in (current, baseline):
                (directory / "raw").mkdir(parents=True)
                (directory / "filtered").mkdir()
                shutil.copy2(HERE / "Summary.csv", directory / "Summary.csv")
                for name in ("barcodes.tsv", "matrix.mtx"):
                    with (HERE / name).open("rb") as source, gzip.open(
                        directory / "raw" / f"{name}.gz", "wb"
                    ) as target:
                        shutil.copyfileobj(source, target)
            with (HERE / "barcodes.tsv").open("rb") as source, gzip.open(
                current / "filtered" / "barcodes.tsv.gz", "wb"
            ) as target:
                shutil.copyfileobj(source, target)
            with (HERE / "historical_filtered_barcodes.tsv").open("rb") as source, gzip.open(
                baseline / "filtered" / "barcodes.tsv.gz", "wb"
            ) as target:
                shutil.copyfileobj(source, target)

            evidence = current / "bam_evidence"
            evidence.mkdir()
            runner.write_tsv(
                evidence / "barcode_read_metrics.tsv.gz",
                [
                    "CB", "current_filtered_member", "historical_filtered_member",
                    "all_records", "primary_mapped_reads",
                    "unique_gene_tagged_reads", "candidate_countedU_reads",
                    "nh_gt1_unique_gene_countedU_reads",
                    "candidate_matrix_molecules", "mitochondrial_reads", "rrna_reads",
                    "biological_classification_status",
                ],
                [
                    {
                        "CB": "AAAAAAAAAAAAAAAA-1", "current_filtered_member": 1,
                        "historical_filtered_member": 1, "all_records": 12,
                        "primary_mapped_reads": 10, "unique_gene_tagged_reads": 8,
                        "candidate_countedU_reads": 8, "candidate_matrix_molecules": 4,
                        "nh_gt1_unique_gene_countedU_reads": 0,
                        "mitochondrial_reads": 1, "rrna_reads": 0,
                        "biological_classification_status": "available",
                    },
                    {
                        "CB": "CCCCCCCCCCCCCCCC-1", "current_filtered_member": 1,
                        "historical_filtered_member": 0, "all_records": 6,
                        "primary_mapped_reads": 5, "unique_gene_tagged_reads": 3,
                        "candidate_countedU_reads": 3, "candidate_matrix_molecules": 3,
                        "nh_gt1_unique_gene_countedU_reads": 1,
                        "mitochondrial_reads": "", "rrna_reads": "",
                        "biological_classification_status": "unavailable",
                    },
                ],
                gzip_output=True,
            )
            barcode_rg_fields = [
                "CB", "RG", "source_id", "source_index", "primary_mapped_reads",
                "unique_gene_tagged_reads", "candidate_countedU_reads",
                "nh_gt1_unique_gene_countedU_reads", "candidate_matrix_molecules",
                "mitochondrial_reads", "rrna_reads",
                "biological_classification_status",
            ]
            runner.write_tsv(
                evidence / "barcode_rg_metrics.tsv.gz",
                barcode_rg_fields,
                [
                    {
                        "CB": "AAAAAAAAAAAAAAAA-1", "RG": "RG1", "source_id": "BP1",
                        "source_index": 0, "primary_mapped_reads": 5,
                        "unique_gene_tagged_reads": 5, "candidate_countedU_reads": 5,
                        "nh_gt1_unique_gene_countedU_reads": 0,
                        "candidate_matrix_molecules": 4, "mitochondrial_reads": 1,
                        "rrna_reads": 0, "biological_classification_status": "available",
                    },
                    {
                        "CB": "CCCCCCCCCCCCCCCC-1", "RG": "RG2", "source_id": "BP2",
                        "source_index": 1, "primary_mapped_reads": 4,
                        "unique_gene_tagged_reads": 3, "candidate_countedU_reads": 3,
                        "nh_gt1_unique_gene_countedU_reads": 0,
                        "candidate_matrix_molecules": 3, "mitochondrial_reads": "",
                        "rrna_reads": "", "biological_classification_status": "unavailable",
                    },
                    {
                        "CB": "CCCCCCCCCCCCCCCC-1", "RG": "RG3", "source_id": "BP2",
                        "source_index": 1, "primary_mapped_reads": 1,
                        "unique_gene_tagged_reads": 0, "candidate_countedU_reads": 0,
                        "nh_gt1_unique_gene_countedU_reads": 1,
                        "candidate_matrix_molecules": 1, "mitochondrial_reads": "",
                        "rrna_reads": "", "biological_classification_status": "unavailable",
                    },
                ],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "rg_summary.tsv",
                ["RG", "source_id", "source_index", "primary_mapped_reads"],
                [
                    {"RG": "RG1", "source_id": "BP1", "source_index": 0, "primary_mapped_reads": 5},
                    {"RG": "RG2", "source_id": "BP2", "source_index": 1, "primary_mapped_reads": 4},
                    {"RG": "RG3", "source_id": "BP2", "source_index": 1, "primary_mapped_reads": 1},
                ],
            )
            runner.write_tsv(
                evidence / "molecule_source_hash_bins.tsv.gz",
                ["CB", "source_mask_hex", "nested_min_hash_bin", "n_hash_bins", "candidate_matrix_molecules"],
                [
                    {"CB": "AAAAAAAAAAAAAAAA-1", "source_mask_hex": "0x3", "nested_min_hash_bin": 2, "n_hash_bins": 100, "candidate_matrix_molecules": 4},
                    {"CB": "CCCCCCCCCCCCCCCC-1", "source_mask_hex": "0x2", "nested_min_hash_bin": 8, "n_hash_bins": 100, "candidate_matrix_molecules": 3},
                ],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "raw_to_corrected_barcode_counts.tsv.gz",
                [
                    "CR", "RG", "CB", "read_count", "within_rg_conflict",
                    "cross_rg_conflict", "global_conflict",
                ],
                [
                    {
                        "CR": "AAAAAAAAAAAAAAAA", "RG": "RG1",
                        "CB": "AAAAAAAAAAAAAAAA-1", "read_count": 1,
                        "within_rg_conflict": 0, "cross_rg_conflict": 0,
                        "global_conflict": 0,
                    },
                    {
                        "CR": "CCCCCCCCCCCCCCCC", "RG": "RG2",
                        "CB": "CCCCCCCCCCCCCCCC-1", "read_count": 1,
                        "within_rg_conflict": 0, "cross_rg_conflict": 0,
                        "global_conflict": 0,
                    },
                ],
                gzip_output=True,
            )
            (evidence / "audit.json").write_text(
                json.dumps({"source_order": ["BP1", "BP2"]}) + "\n",
                encoding="utf-8",
            )

            gained_output = root / "gained"
            analyze_gained_cells.run(argparse.Namespace(
                current_root=str(current_root), baseline_root=str(baseline_root),
                output_dir=str(gained_output), demux_tsv=None, trim_tsv=None,
                ambient_similarity_tsv=None,
            ))
            gained_rows = read_tsv(gained_output / "cell_category_metrics.tsv.gz")
            self.assertEqual({row["cell_transition"] for row in gained_rows}, {"shared", "gained"})
            unavailable = next(
                row for row in gained_rows
                if row["CB"] == "CCCCCCCCCCCCCCCC-1"
            )
            self.assertEqual(unavailable["biological_classification_status"], "unavailable")
            self.assertEqual(unavailable["mitochondrial_primary_read_fraction"], "")
            self.assertEqual(unavailable["rrna_primary_read_fraction"], "")

            inventory = root / "inventory.tsv"
            runner.write_tsv(
                inventory, ["library", "bp_id", "raw_reads"],
                [
                    {"library": library, "bp_id": "BP1", "raw_reads": 100},
                    {"library": library, "bp_id": "BP2", "raw_reads": 200},
                ],
            )
            depth_output = root / "depth"
            analyze_depth_history.run(argparse.Namespace(
                current_root=str(current_root), baseline_root=str(baseline_root),
                fastq_inventory=str(inventory), output_dir=str(depth_output),
                thresholds=[1], random_steps=2,
            ))
            depth_rows = read_tsv(depth_output / "depth_trajectories.tsv")
            prefixes = [
                row for row in depth_rows
                if row["trajectory_type"] == "rg_prefix_threshold_proxy"
            ]
            self.assertEqual(len(prefixes), 2)
            self.assertEqual(prefixes[-1]["source_prefix"], "BP1|BP2")

            trim = root / "trim.tsv"
            runner.write_tsv(
                trim,
                [
                    "library", "cell_barcode", "source_id", "read_mate",
                    "total_reads", "reads_with_adapter", "total_bp_trimmed",
                    "frac_reads_with_adapter", "mean_bp_trimmed",
                ],
                [
                    {"library": library, "cell_barcode": "AAAAAAAAAAAAAAAA", "source_id": "RG1", "read_mate": "R1", "total_reads": 10, "reads_with_adapter": 1, "total_bp_trimmed": 10, "frac_reads_with_adapter": 0.1, "mean_bp_trimmed": 1.0},
                    {"library": library, "cell_barcode": "AAAAAAAAAAAAAAAA", "source_id": "RG1", "read_mate": "R2", "total_reads": 10, "reads_with_adapter": 2, "total_bp_trimmed": 30, "frac_reads_with_adapter": 0.2, "mean_bp_trimmed": 3.0},
                    {"library": library, "cell_barcode": "AAAAAAAAAAAAAAAA", "source_id": "NO_BAM_RG", "read_mate": "R2", "total_reads": 6, "reads_with_adapter": 6, "total_bp_trimmed": 120, "frac_reads_with_adapter": 1.0, "mean_bp_trimmed": 20.0},
                    {"library": library, "cell_barcode": "CCCCCCCCCCCCCCCC", "source_id": "RG2", "read_mate": "R2", "total_reads": 8, "reads_with_adapter": 2, "total_bp_trimmed": 16, "frac_reads_with_adapter": 0.25, "mean_bp_trimmed": 2.0},
                    {"library": library, "cell_barcode": "CCCCCCCCCCCCCCCC", "source_id": "RG3", "read_mate": "R1", "total_reads": 4, "reads_with_adapter": 3, "total_bp_trimmed": 20, "frac_reads_with_adapter": 0.75, "mean_bp_trimmed": 5.0},
                ],
            )
            trim_output = root / "trim_join"
            join_trim_cell_metrics.run(argparse.Namespace(
                current_root=str(current_root), trim_by_source_tsv=str(trim),
                output_dir=str(trim_output),
            ))
            joined = read_tsv(trim_output / "trim_cell_metrics.tsv.gz")
            self.assertEqual(len(joined), 5)
            unmatched_trim = next(row for row in joined if row["source_id"] == "NO_BAM_RG")
            self.assertEqual(unmatched_trim["bam_match_status"], "no_bam_evidence")
            self.assertEqual(unmatched_trim["frac_reads_with_adapter"], "1.0")
            self.assertEqual(unmatched_trim["bam_primary_mapped_reads"], "")
            self.assertEqual(
                unmatched_trim["bam_biological_classification_status"],
                "unavailable_no_bam_evidence",
            )
            matched_joined = [row for row in joined if row["bam_match_status"] == "matched"]
            self.assertEqual({row["matrix_umis"] for row in matched_joined}, {"4", "5"})
            self.assertTrue(all(row["matrix_value_authority"] == "STARsolo_raw_matrix" for row in matched_joined))
            unavailable = next(
                row for row in joined
                if row["cell_barcode"] == "CCCCCCCCCCCCCCCC" and row["source_id"] == "RG2"
            )
            self.assertEqual(
                unavailable["bam_biological_classification_status"],
                "unavailable",
            )
            self.assertEqual(unavailable["bam_mitochondrial_primary_fraction"], "")
            self.assertEqual(unavailable["bam_rrna_primary_fraction"], "")
            aggregate = read_tsv(
                trim_output / "trim_cell_metrics_by_barcode.tsv.gz"
            )
            self.assertEqual(len(aggregate), 2)
            self.assertEqual(
                len({(row["library"], row["cell_barcode"]) for row in aggregate}),
                2,
            )
            a_row = next(
                row for row in aggregate if row["cell_barcode"] == "AAAAAAAAAAAAAAAA"
            )
            self.assertEqual(a_row["bam_match_status"], "partial")
            self.assertEqual(a_row["trim_row_count"], "3")
            self.assertEqual(a_row["bam_primary_mapped_reads"], "10")
            self.assertEqual(a_row["read_mates"], "R1,R2")
            self.assertEqual(a_row["source_ids"], "NO_BAM_RG,RG1")
            self.assertAlmostEqual(float(a_row["frac_reads_with_adapter"]), 9 / 26)
            c_row = next(
                row for row in aggregate if row["cell_barcode"] == "CCCCCCCCCCCCCCCC"
            )
            self.assertEqual(c_row["bam_match_status"], "matched")
            self.assertEqual(c_row["RG"], "RG2,RG3")
            self.assertEqual(c_row["read_mates"], "R1,R2")
            self.assertEqual(c_row["bam_primary_mapped_reads"], "5")
            distributions = read_tsv(trim_output / "trim_only_distributions.tsv")
            self.assertTrue(any(
                row["bam_match_status"] == "no_bam_evidence"
                and row["metric"] == "frac_reads_with_adapter"
                for row in distributions
            ))

            # The one-row-per-barcode product is directly consumable by the
            # gained-cell optional join without duplicate-key failure.
            gained_with_trim = root / "gained_with_trim"
            analyze_gained_cells.run(argparse.Namespace(
                current_root=str(current_root), baseline_root=str(baseline_root),
                output_dir=str(gained_with_trim), demux_tsv=None,
                trim_tsv=str(trim_output / "trim_cell_metrics_by_barcode.tsv.gz"),
                ambient_similarity_tsv=None,
            ))
            gained_trim_rows = read_tsv(
                gained_with_trim / "cell_category_metrics.tsv.gz"
            )
            self.assertTrue(all("trim_bam_match_status" in row for row in gained_trim_rows))

            source_yield = root / "source_yield.tsv"
            runner.write_tsv(
                source_yield,
                [
                    "library", "source_id", *export_repooling_evidence.SUM_FIELDS,
                    "marginal_candidate_molecules_first_observed",
                ],
                [
                    {
                        "library": library, "source_id": "BP1",
                        **{field: 1 for field in export_repooling_evidence.SUM_FIELDS},
                        "marginal_candidate_molecules_first_observed": 1,
                    },
                    {
                        "library": library, "source_id": "BP2",
                        **{field: 2 for field in export_repooling_evidence.SUM_FIELDS},
                        "marginal_candidate_molecules_first_observed": 2,
                    },
                ],
            )
            repool_output = root / "repool"
            export_repooling_evidence.run(argparse.Namespace(
                fastq_inventory=str(inventory), source_yield_tsv=str(source_yield),
                output_dir=str(repool_output), library_column="library",
                source_column="bp_id", raw_reads_column="raw_reads",
            ))
            repool = read_tsv(repool_output / "repooling_evidence.tsv")
            self.assertTrue(all(row["bam_evidence_match"] == "1" for row in repool))

            plot_output = root / "plots"
            processedstats = root / "processedstats.tsv"
            runner.write_tsv(
                processedstats,
                [
                    "Library", "Reads Mapped to Genome: Unique+Multiple",
                    "Reads Mapped to Genome: Unique",
                    "Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS",
                    "Reads With Valid Barcodes", "Fraction of Unique Reads in Cells",
                    "Sequencing Saturation",
                ],
                [{
                    "Library": library,
                    "Reads Mapped to Genome: Unique+Multiple": 0.9,
                    "Reads Mapped to Genome: Unique": 0.8,
                    "Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS": 0.7,
                    "Reads With Valid Barcodes": 0.95,
                    "Fraction of Unique Reads in Cells": 0.6,
                    "Sequencing Saturation": 0.5,
                }],
            )
            plot_rna_evidence.run(argparse.Namespace(
                cell_category_tsv=str(gained_output / "cell_category_metrics.tsv.gz"),
                depth_trajectories_tsv=str(depth_output / "depth_trajectories.tsv"),
                source_yield_tsv=str(source_yield),
                processedstats_tsv=str(processedstats),
                output_dir=str(plot_output),
            ))
            self.assertTrue((plot_output / "cell_category_distributions.png").is_file())
            self.assertTrue((plot_output / "plot_values_depth_trajectories.tsv").is_file())
            self.assertTrue((plot_output / "depth_cell_count_curves.png").is_file())
            self.assertTrue((plot_output / "fixed_roster_trajectories.png").is_file())
            self.assertTrue((plot_output / "library_mapping_context.png").is_file())
            panel_rows = read_tsv(
                plot_output / "plot_values_library_mapping_context.tsv"
            )
            self.assertEqual(len(panel_rows), 6)
            histogram_rows = read_tsv(
                gained_output / "cell_category_histograms.tsv"
            )
            self.assertTrue(any(int(row["count"]) > 0 for row in histogram_rows))
            histogram_metrics = {row["metric"] for row in histogram_rows}
            self.assertTrue(
                {
                    "candidate_countedU_reads",
                    "current_raw_umis",
                    "current_raw_genes",
                    "bam_conditional_saturation",
                    "primary_mapping_record_fraction",
                    "unique_gene_tagged_primary_fraction",
                    "mitochondrial_primary_read_fraction",
                    "rrna_primary_read_fraction",
                    "source_dominance_fraction",
                }.issubset(histogram_metrics)
            )

    def test_bam_free_gather_streams_compact_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bam_free_gather_") as temporary:
            root = Path(temporary)
            evidence = root / "library" / "bam_evidence"
            evidence.mkdir(parents=True)
            metric_values = {field: 0 for field in runner.CORE_ADDITIVE_FIELDS}
            metric_values.update({"all_records": 2, "primary_mapped_reads": 2, "candidate_countedU_reads": 2})
            runner.write_tsv(
                evidence / "barcode_read_metrics.tsv.gz",
                [
                    "CB", "current_raw_member", "current_filtered_member",
                    "historical_raw_member", "historical_filtered_member", "barcode_category",
                    "candidate_matrix_molecules", "biological_classification_status",
                    *runner.CORE_ADDITIVE_FIELDS,
                ],
                [{
                    "CB": "A-1", "current_raw_member": 1, "current_filtered_member": 1,
                    "historical_raw_member": 0, "historical_filtered_member": 0,
                    "barcode_category": "current filtered cell",
                    "candidate_matrix_molecules": 1,
                    "biological_classification_status": "unavailable",
                    **metric_values,
                }],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "barcode_rg_metrics.tsv.gz",
                [
                    "CB", "RG", "source_id", "source_index",
                    "biological_classification_status",
                    *runner.CORE_ADDITIVE_FIELDS,
                ],
                [{
                    "CB": "A-1", "RG": "RG1", "source_id": "BP1",
                    "source_index": 0,
                    "biological_classification_status": "unavailable",
                    **metric_values,
                }],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "molecule_source_hash_bins.tsv.gz",
                ["CB", "source_mask_hex", "nested_min_hash_bin", "n_hash_bins", "candidate_matrix_molecules"],
                [{"CB": "A-1", "source_mask_hex": "0x1", "nested_min_hash_bin": 0, "n_hash_bins": 100, "candidate_matrix_molecules": 1}],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "rg_summary.tsv",
                ["library", "RG", "source_id", "source_index", "candidate_matrix_molecules", *runner.CORE_ADDITIVE_FIELDS],
                [{"library": "L1", "RG": "RG1", "source_id": "BP1", "source_index": 0, "candidate_matrix_molecules": 1, **metric_values}],
            )
            runner.write_tsv(
                evidence / "rg_contig_class_summary.tsv.gz",
                [
                    "library", "RG", "source_index", "contig_index", "contig",
                    "contig_class_available", "species", "mitochondrial", "rrna",
                    "reference_class", "biological_classification_status",
                ],
                [{
                    "library": "L1", "RG": "RG1", "source_index": 0,
                    "contig_index": 0, "contig": "chr1",
                    "contig_class_available": 0, "species": "",
                    "mitochondrial": "", "rrna": "", "reference_class": "",
                    "biological_classification_status": "unavailable",
                }],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "raw_to_corrected_barcode_counts.tsv.gz",
                [
                    "CR", "RG", "CB", "read_count", "within_rg_conflict",
                    "cross_rg_conflict", "global_conflict",
                ],
                [{
                    "CR": "A", "RG": "RG1", "CB": "A-1", "read_count": 2,
                    "within_rg_conflict": 0, "cross_rg_conflict": 0,
                    "global_conflict": 0,
                }],
                gzip_output=True,
            )
            runner.write_tsv(
                evidence / "CellReads.countedU.from_bam.tsv.gz",
                ["CB", "countedU"],
                [{"CB": "A-1", "countedU": 2}],
                gzip_output=True,
            )
            product_stats = runner.validate_published(evidence)
            compatibility_stats = runner.validate_compatibility(
                evidence / "CellReads.countedU.from_bam.tsv.gz"
            )
            (evidence / "audit.json").write_text(
                json.dumps({
                    "library": "L1", "reconciliation": {"status": "PASS"},
                    "resume_signature": {"fixture": True},
                    "profiler_summary": {
                        "total_records": "2", "primary_mapped_reads": "2",
                        "candidate_countedU_reads": "2",
                        "nh_gt1_unique_gene_countedU_reads": "0",
                        "candidate_matrix_molecules": "1",
                    },
                    "outputs": product_stats,
                    "compatibility_product": compatibility_stats,
                }) + "\n",
                encoding="utf-8",
            )
            (evidence / "BAM_EVIDENCE_COMPLETE.ok").write_text("ok\n", encoding="utf-8")
            manifest = root / "manifest.tsv"
            fields = sorted(runner.REQUIRED_MANIFEST_COLUMNS)
            values = {field: "fixture" for field in fields}
            source_order = root / "source_order.tsv"
            source_order.write_text("source_id\nBP1\n", encoding="utf-8")
            profiler = root / "rna_bam_evidence"
            profiler.write_bytes(b"synthetic-profiler\n")
            scientific_config_path = root / "scientific_config.json"
            values.update({
                "library": "L1", "library_dir": str(root / "library"),
                "output_dir": str(evidence),
                "source_order": str(source_order),
                "scientific_config": str(scientific_config_path),
                "starsolo_feature": runner.STARSOLO_FEATURE,
                "starsolo_umi_filtering": runner.STARSOLO_UMI_FILTERING,
                "starsolo_umi_dedup": runner.STARSOLO_UMI_DEDUP,
                "starsolo_multimappers": runner.STARSOLO_MULTIMAPPERS,
                "biological_classification_intent": "explicit_unavailable",
            })
            for input_column in (
                "bam", "bam_index", "summary", "raw_barcodes",
                "filtered_barcodes", "raw_features", "raw_matrix",
                "filtered_matrix", "rg_metadata",
            ):
                input_path = root / f"fixture_{input_column}"
                input_path.write_text(
                    f"synthetic {input_column}\n", encoding="utf-8"
                )
                values[input_column] = str(input_path.resolve())
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(values)
            source_contents = source_order.read_text(encoding="utf-8")
            scientific_config = {
                "schema_version": 1,
                "manifest_fields": fields,
                "manifest_rows": [values],
                "source_order": {
                    "path": str(source_order), "provenance": "fixture",
                    "values": ["BP1"],
                    "sha256": runner.sha256(source_order),
                    "contents": source_contents,
                },
                "biological_classification": {
                    "status": "explicit_unavailable", "path": "",
                    "sha256": "", "contents": "",
                },
                "hash": {
                    "algorithm": runner.HASH_ALGORITHM,
                    "seed": runner.DEFAULT_HASH_SEED, "bins": 100,
                },
                "profiler_executable": {
                    "path": str(profiler.resolve()),
                    "sha256": runner.sha256(profiler),
                    "bytes": profiler.stat().st_size,
                },
                "starsolo": {
                    "feature": runner.STARSOLO_FEATURE,
                    "umi_filtering": runner.STARSOLO_UMI_FILTERING,
                    "umi_dedup": runner.STARSOLO_UMI_DEDUP,
                    "multimappers": runner.STARSOLO_MULTIMAPPERS,
                    "ordinary_countedU_read_definition":
                        runner.ORDINARY_COUNTEDU_READ_DEFINITION,
                    "ordinary_molecule_definition": runner.ORDINARY_MOLECULE_DEFINITION,
                    "multimapper_definition": runner.MULTIMAPPER_DEFINITION,
                    "summary_unique_read_metric": runner.SUMMARY_UNIQUE_READ_METRIC,
                    "nh_gt1_unique_gene_countedU_definition":
                        runner.NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION,
                    "starsolo_EM_evidence_availability":
                        runner.STARSOLO_EM_EVIDENCE_AVAILABILITY,
                },
            }
            scientific_config["configuration_hash"] = runner.canonical_payload_hash(
                scientific_config
            )
            scientific_config_path.write_text(
                json.dumps(scientific_config, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            expected_parameters = runner.scientific_parameters(
                values, scientific_config
            )
            fixture_audit = json.loads(
                (evidence / "audit.json").read_text(encoding="utf-8")
            )
            fixture_profiler_version = "rna_bam_evidence 2.3.0"
            fixture_audit.update({
                "release": runner.RELEASE,
                "profiler_version": fixture_profiler_version,
                "manifest": str(manifest.resolve()),
                "manifest_row_index": 0,
                "scientific_configuration": {
                    "path": str(scientific_config_path.resolve()),
                    "configuration_hash": scientific_config["configuration_hash"],
                },
            })
            fixture_audit["resume_signature"] = {
                "release": runner.RELEASE,
                "library": "L1",
                "manifest_row": values,
                "inputs": {
                    column: runner.file_signature(Path(values[column]))
                    for column in runner.SIGNATURE_INPUT_COLUMNS
                    if values.get(column)
                },
                "profiler": fixture_profiler_version,
                "profiler_executable": runner.file_signature(profiler),
                "parameters": expected_parameters,
                "scientific_configuration_hash": scientific_config["configuration_hash"],
            }
            fixture_audit["parameters"] = expected_parameters
            (evidence / "audit.json").write_text(
                json.dumps(fixture_audit, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            output = root / "gather"
            runner.gather(argparse.Namespace(
                manifest=str(manifest), output_dir=str(output), replace_stale=False
            ))
            self.assertTrue((output / "BAM_EVIDENCE_GATHER_COMPLETE.ok").is_file())
            source_rows = read_tsv(output / "all_libraries_source_yield.tsv")
            self.assertEqual(source_rows[0]["candidate_matrix_molecules"], "1")
            self.assertEqual(source_rows[0]["marginal_candidate_molecules_first_observed"], "1")
            self.assertEqual(
                source_rows[0]["biological_classification_status"],
                "unavailable",
            )
            self.assertEqual(source_rows[0]["mitochondrial_reads"], "")
            self.assertEqual(source_rows[0]["rrna_reads"], "")

            # Product checksums alone do not make a stale pilot compatible.
            passing_audit = (evidence / "audit.json").read_text(encoding="utf-8")
            stale_audit = json.loads(passing_audit)
            stale_audit["resume_signature"]["manifest_row"] = {
                **values, "summary": "different-pilot-summary.csv",
            }
            (evidence / "audit.json").write_text(
                json.dumps(stale_audit, sort_keys=True) + "\n", encoding="utf-8"
            )
            with self.assertRaises(runner.EvidenceError):
                runner.gather(argparse.Namespace(
                    manifest=str(manifest), output_dir=str(output),
                    replace_stale=False,
                ))
            (evidence / "audit.json").write_text(passing_audit, encoding="utf-8")

            # A completion marker is not trusted on resume. Every gathered
            # product must still match its audited schema, row count, size,
            # and digest.
            (output / "all_libraries_source_yield.tsv").write_text(
                "library\tsource_id\nL1\tCORRUPTED\n", encoding="utf-8"
            )
            with self.assertRaises(runner.EvidenceError):
                runner.gather(argparse.Namespace(
                    manifest=str(manifest), output_dir=str(output),
                    replace_stale=False,
                ))
            runner.gather(argparse.Namespace(
                manifest=str(manifest), output_dir=str(output),
                replace_stale=True,
            ))
            (output / "all_libraries_barcode_summary.tsv.gz").unlink()
            with self.assertRaises(runner.EvidenceError):
                runner.gather(argparse.Namespace(
                    manifest=str(manifest), output_dir=str(output),
                    replace_stale=False,
                ))
            runner.gather(argparse.Namespace(
                manifest=str(manifest), output_dir=str(output),
                replace_stale=True,
            ))
            barcode_rg_path = evidence / "barcode_rg_metrics.tsv.gz"
            barcode_rg_bytes = barcode_rg_path.read_bytes()
            barcode_rg_path.unlink()
            with self.assertRaises(runner.EvidenceError):
                runner.gather(argparse.Namespace(
                    manifest=str(manifest), output_dir=str(output),
                    replace_stale=False,
                ))
            barcode_rg_path.write_bytes(barcode_rg_bytes)

            # Gather remains BAM-free but must reject changed scientific input
            # identities even when every evidence product and checksum survives.
            (root / "fixture_summary").write_text(
                "changed summary input\n", encoding="utf-8"
            )
            with self.assertRaises(runner.EvidenceError):
                runner.gather(argparse.Namespace(
                    manifest=str(manifest), output_dir=str(output),
                    replace_stale=False,
                ))

    def test_large_sorted_product_validation_has_bounded_python_memory(self) -> None:
        with tempfile.TemporaryDirectory(prefix="streaming_validation_") as temporary:
            path = Path(temporary) / "barcode_read_metrics.tsv.gz"
            fields = sorted(runner.SCHEMAS["barcode_read_metrics.tsv.gz"])

            def rows():
                for index in range(120_000):
                    row = {field: 0 for field in fields}
                    row.update({
                        "CB": f"CB{index:08d}-1",
                        "barcode_category": "raw nonfiltered droplet",
                        "biological_classification_status": "unavailable",
                        "mitochondrial_reads": "",
                        "rrna_reads": "",
                    })
                    yield row

            runner.write_tsv(path, fields, rows(), gzip_output=True)
            gc.collect()
            tracemalloc.start()
            stats = runner.validate_product(
                path, "barcode_read_metrics.tsv.gz"
            )
            _current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            self.assertEqual(stats["rows"], 120_000)
            self.assertLess(
                peak, 16 * 1024 * 1024,
                f"streaming validation peak was {peak} bytes",
            )

    def test_large_trim_left_join_uses_disk_backed_project_state(self) -> None:
        def measured_peak(row_count: int) -> int:
            with tempfile.TemporaryDirectory(prefix="large_trim_join_") as temporary:
                root = Path(temporary)
                current = root / "current"
                placeholder = current / "MappedButNotInTrim"
                placeholder.mkdir(parents=True)
                (placeholder / "Summary.csv").write_text(
                    "metric,value\n", encoding="utf-8"
                )
                trim = root / "trim.tsv"
                with trim.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=[
                            "library", "cell_barcode", "source_id", "read_mate",
                            "total_reads", "reads_with_adapter", "total_bp_trimmed",
                            "frac_reads_with_adapter", "mean_bp_trimmed",
                        ],
                        delimiter="\t", lineterminator="\n",
                    )
                    writer.writeheader()
                    for index in range(row_count):
                        writer.writerow({
                            "library": "UnmappedTrimLibrary",
                            "cell_barcode": f"CB{index:08d}",
                            "source_id": f"RG{index % 4}",
                            "read_mate": "R1" if index % 2 else "R2",
                            "total_reads": 10,
                            "reads_with_adapter": index % 11,
                            "total_bp_trimmed": index % 101,
                            "frac_reads_with_adapter": (index % 11) / 10,
                            "mean_bp_trimmed": (index % 101) / 10,
                        })
                gc.collect()
                tracemalloc.start()
                join_trim_cell_metrics.run(argparse.Namespace(
                    current_root=str(current), trim_by_source_tsv=str(trim),
                    output_dir=str(root / "output"),
                ))
                _current, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                with gzip.open(
                    root / "output" / "trim_cell_metrics.tsv.gz", "rt",
                    encoding="utf-8",
                ) as handle:
                    self.assertEqual(sum(1 for _ in handle) - 1, row_count)
                return peak

        small_peak = measured_peak(2_000)
        large_peak = measured_peak(30_000)
        self.assertLess(large_peak, 32 * 1024 * 1024)
        self.assertLess(large_peak, small_peak * 2.5)

    def test_large_raw_matrix_roster_uses_disk_backed_numeric_columns(self) -> None:
        with tempfile.TemporaryDirectory(prefix="matrix_mmap_streaming_") as temporary:
            root = Path(temporary)
            matrix = root / "matrix.mtx"
            roster = root / "barcodes.tsv"
            matrix.write_text(
                "%%MatrixMarket matrix coordinate integer general\n"
                "% sparse large-roster fixture\n"
                "2 200000 2\n"
                "1 1 1\n"
                "2 200000 3\n",
                encoding="utf-8",
            )
            with roster.open("w", encoding="utf-8") as handle:
                for index in range(200_000):
                    handle.write(f"CB{index:08d}-1\n")
            gc.collect()
            tracemalloc.start()
            rows = list(
                rna_evidence_common.iter_nonzero_matrix_barcodes(
                    matrix, roster
                )
            )
            _current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            self.assertEqual(
                rows,
                [("CB00000000-1", 1, 1), ("CB00199999-1", 3, 1)],
            )
            self.assertLess(
                peak, 8 * 1024 * 1024,
                f"disk-backed matrix streaming peak was {peak} bytes",
            )

    def test_historical_star_diagnostics_are_promoted_without_bam_scan(self) -> None:
        with tempfile.TemporaryDirectory(prefix="star_diagnostic_fixture_") as temporary:
            root = Path(temporary)
            library_dir = root / "published" / "FixtureRNA"
            library_dir.mkdir(parents=True)
            summary = library_dir / "Summary.csv"
            summary.write_text("Estimated Number of Cells,2\n", encoding="utf-8")
            task = root / "work" / "ab" / "cdef"
            task.mkdir(parents=True)
            shutil.copy2(summary, task / "Summary.csv")
            (task / "STAR_Log.out").write_text("full log\n", encoding="utf-8")
            (task / "STAR_Log.final.out").write_text(
                "final log\n", encoding="utf-8"
            )
            (task / "STAR_SJ.out.tab").write_text(
                "chr1\t1\t2\n", encoding="utf-8"
            )
            (library_dir / "STAR_Log.out").write_text(
                "stale arbitrary log\n", encoding="utf-8"
            )
            (library_dir / "STAR_Log.final.out").write_text(
                "stale arbitrary final log\n", encoding="utf-8"
            )
            (library_dir / "STAR_SJ.out.tab.gz").write_bytes(
                b"not a gzip diagnostic\n"
            )
            manifest = root / "manifest.tsv"
            runner.write_tsv(
                manifest,
                ["library", "library_dir", "summary"],
                [{
                    "library": "FixtureRNA",
                    "library_dir": library_dir,
                    "summary": summary,
                }],
            )
            audit = root / "promotion.json"
            marker = root / "promotion.ok"
            promote_rna_star_diagnostics.run(argparse.Namespace(
                manifest=str(manifest), work_root=str(root / "work"),
                audit=str(audit), marker=str(marker),
            ))
            self.assertEqual(
                (library_dir / "STAR_Log.out").read_text(encoding="utf-8"),
                "full log\n",
            )
            self.assertEqual(
                (library_dir / "STAR_Log.final.out").read_text(encoding="utf-8"),
                "final log\n",
            )
            with gzip.open(
                library_dir / "STAR_SJ.out.tab.gz", "rt", encoding="utf-8"
            ) as handle:
                self.assertEqual(handle.read(), "chr1\t1\t2\n")
            payload = json.loads(audit.read_text(encoding="utf-8"))
            self.assertFalse(payload["bam_opened"])
            self.assertFalse(payload["retired_awk_scanner_used"])
            self.assertFalse(payload["work_directory_deleted"])
            self.assertTrue(marker.is_file())
            promote_rna_star_diagnostics.validate_prior_promotion(
                manifest.resolve(), (root / "work").resolve(),
                audit.resolve(), marker.resolve(),
            )
            outputs_before = {
                name: (library_dir / name).read_bytes()
                for name in promote_rna_star_diagnostics.OUTPUT_NAMES
            }
            shutil.rmtree(root / "work")
            promote_rna_star_diagnostics.run(argparse.Namespace(
                manifest=str(manifest), work_root=str(root / "work"),
                audit=str(audit), marker=str(marker),
            ))
            self.assertEqual(
                outputs_before,
                {
                    name: (library_dir / name).read_bytes()
                    for name in promote_rna_star_diagnostics.OUTPUT_NAMES
                },
            )
            (library_dir / "STAR_Log.out").write_text(
                "tampered\n", encoding="utf-8"
            )
            with self.assertRaises(
                promote_rna_star_diagnostics.PromotionError
            ):
                promote_rna_star_diagnostics.validate_prior_promotion(
                    manifest.resolve(), (root / "work").resolve(),
                    audit.resolve(), marker.resolve(),
                )

    def test_arbitrary_star_diagnostics_without_audit_are_not_accepted(self) -> None:
        with tempfile.TemporaryDirectory(prefix="star_diagnostic_stale_") as temporary:
            root = Path(temporary)
            library_dir = root / "published" / "FixtureRNA"
            library_dir.mkdir(parents=True)
            summary = library_dir / "Summary.csv"
            summary.write_text("Estimated Number of Cells,2\n", encoding="utf-8")
            for name in promote_rna_star_diagnostics.OUTPUT_NAMES:
                (library_dir / name).write_text("arbitrary\n", encoding="utf-8")
            work = root / "work"
            work.mkdir()
            manifest = root / "manifest.tsv"
            runner.write_tsv(
                manifest,
                ["library", "library_dir", "summary"],
                [{
                    "library": "FixtureRNA", "library_dir": library_dir,
                    "summary": summary,
                }],
            )
            marker = root / "promotion.ok"
            with self.assertRaises(
                promote_rna_star_diagnostics.PromotionError
            ):
                promote_rna_star_diagnostics.run(argparse.Namespace(
                    manifest=str(manifest), work_root=str(work),
                    audit=str(root / "promotion.json"), marker=str(marker),
                ))
            self.assertFalse(marker.exists())


class CompiledFixtureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.samtools = shutil.which("samtools")
        configured = os.environ.get("RNA_BAM_EVIDENCE_BIN")
        cls.build_temporary = None
        if configured:
            cls.profiler = Path(configured).resolve()
        elif (ROOT / "rna_bam_evidence").is_file():
            cls.profiler = ROOT / "rna_bam_evidence"
        else:
            pkg_config = shutil.which("pkg-config")
            has_htslib = bool(
                pkg_config
                and subprocess.run(
                    [pkg_config, "--exists", "htslib"], check=False
                ).returncode == 0
            )
            if has_htslib:
                cls.build_temporary = tempfile.TemporaryDirectory(
                    prefix="real_htslib_build_"
                )
                cls.profiler = (
                    Path(cls.build_temporary.name) / "rna_bam_evidence"
                )
                flags = shlex.split(subprocess.check_output(
                    [pkg_config, "--cflags", "--libs", "htslib"],
                    text=True,
                ))
                subprocess.run(
                    [
                        "g++", "-std=c++17", "-O2", "-Wall", "-Wextra",
                        str(ROOT / "src" / "rna_bam_evidence.cpp"),
                        "-o", str(cls.profiler), *flags, "-lz", "-lpthread",
                    ],
                    check=True,
                )
            else:
                cls.profiler = ROOT / "rna_bam_evidence"
        if cls.samtools is None or not cls.profiler.is_file():
            raise unittest.SkipTest(
                "requires samtools plus real htslib or a built profiler "
                "(set RNA_BAM_EVIDENCE_BIN)"
            )

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.build_temporary is not None:
            cls.build_temporary.cleanup()

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rna_bam_evidence_fixture_")
        self.library = Path(self.temporary.name) / "FixtureRNA"
        (self.library / "raw").mkdir(parents=True)
        (self.library / "filtered").mkdir(parents=True)
        for destination in (
            self.library / "raw" / "barcodes.tsv.gz",
            self.library / "filtered" / "barcodes.tsv.gz",
        ):
            with (HERE / "barcodes.tsv").open("rb") as source, gzip.open(destination, "wb") as target:
                shutil.copyfileobj(source, target)
        for destination in (
            self.library / "raw" / "matrix.mtx.gz",
            self.library / "filtered" / "matrix.mtx.gz",
        ):
            with (HERE / "matrix.mtx").open("rb") as source, gzip.open(destination, "wb") as target:
                shutil.copyfileobj(source, target)
        with (HERE / "features.tsv").open("rb") as source, gzip.open(
            self.library / "raw" / "features.tsv.gz", "wb"
        ) as target:
            shutil.copyfileobj(source, target)
        shutil.copy2(HERE / "Summary.csv", self.library / "Summary.csv")
        shutil.copy2(HERE / "CellReads.stats", self.library / "CellReads.stats")
        self.bam = self.library / "gex.bam"
        subprocess.run(
            [self.samtools, "sort", "-o", str(self.bam), str(HERE / "fixture.sam")],
            check=True,
        )
        subprocess.run([self.samtools, "index", str(self.bam)], check=True)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def profiler_command(
        self, output: Path, threads: int, *, include_class_manifest: bool = True
    ) -> list[str]:
        output.mkdir()
        command = [
            str(self.profiler),
            "--bam", str(self.bam),
            "--output-dir", str(output),
            "--library", "FixtureRNA",
            "--raw-barcodes", str(self.library / "raw" / "barcodes.tsv.gz"),
            "--raw-matrix", str(self.library / "raw" / "matrix.mtx.gz"),
            "--filtered-barcodes", str(self.library / "filtered" / "barcodes.tsv.gz"),
            "--filtered-matrix", str(self.library / "filtered" / "matrix.mtx.gz"),
            "--features", str(self.library / "raw" / "features.tsv.gz"),
            "--old-raw-barcodes", str(HERE / "barcodes.tsv"),
            "--old-filtered-barcodes", str(HERE / "historical_filtered_barcodes.tsv"),
            "--rg-metadata", str(HERE / "rg_metadata.tsv"),
            "--source-order", str(HERE / "source_order.tsv"),
            "--starsolo-feature", runner.STARSOLO_FEATURE,
            "--starsolo-umi-filtering", runner.STARSOLO_UMI_FILTERING,
            "--starsolo-umi-dedup", runner.STARSOLO_UMI_DEDUP,
            "--starsolo-multimappers", runner.STARSOLO_MULTIMAPPERS,
            "--threads", str(threads),
            "--hash-bins", "100",
            "--hash-seed", str(runner.DEFAULT_HASH_SEED),
            "--expected-molecules", "9",
            # Leave enough virtual-address headroom for real htslib BGZF
            # worker stacks and shared libraries. The production memory-bound
            # behavior is asserted from the profiler audit, not by forcing the
            # tiny fixture into an unrealistically small RLIMIT_AS.
            "--max-memory-bytes", str(1024 * 1024 * 1024),
        ]
        if include_class_manifest:
            command.extend(["--class-manifest", str(HERE / "contig_classes.tsv")])
        return command

    def test_bam_preflight_uses_the_manifest_declared_index(self) -> None:
        default_index = Path(str(self.bam) + ".bai")
        declared_index = self.library / "declared-custom-index.bai"
        default_index.replace(declared_index)

        result = runner.bam_preflight(
            {"bam": self.bam, "bam_index": declared_index}
        )
        self.assertEqual(result["index_probe_region"], "chr1:1-1")
        self.assertGreaterEqual(result["index_probe_records"], 0)

    def test_bam_preflight_rejects_a_bad_declared_index(self) -> None:
        bad_index = self.library / "bad-declared-index.bai"
        bad_index.write_bytes(b"not a BAM index\n")

        with self.assertRaises(runner.EvidenceError):
            runner.bam_preflight({"bam": self.bam, "bam_index": bad_index})

    def test_exact_counters_masks_conflicts_and_threads(self) -> None:
        one = Path(self.temporary.name) / "one"
        two = Path(self.temporary.name) / "two"
        subprocess.run(self.profiler_command(one, 1), check=True)
        subprocess.run(self.profiler_command(two, 2), check=True)
        products = [*runner.COMPRESSED_PRODUCTS.keys(), *runner.PLAIN_PRODUCTS.keys()]
        for name in products:
            self.assertEqual((one / name).read_bytes(), (two / name).read_bytes(), name)

        barcode_rows = {row["CB"]: row for row in read_tsv(one / "barcode_read_metrics.tsv")}
        self.assertEqual(int(barcode_rows["AAAAAAAAAAAAAAAA-1"]["candidate_countedU_reads"]), 11)
        self.assertEqual(int(barcode_rows["AAAAAAAAAAAAAAAA-1"]["candidate_matrix_molecules"]), 4)
        self.assertEqual(int(barcode_rows["AAAAAAAAAAAAAAAA-1"]["nh_gt1_unique_gene_countedU_reads"]), 0)
        self.assertEqual(int(barcode_rows["CCCCCCCCCCCCCCCC-1"]["candidate_countedU_reads"]), 5)
        self.assertEqual(int(barcode_rows["CCCCCCCCCCCCCCCC-1"]["candidate_matrix_molecules"]), 5)
        self.assertEqual(int(barcode_rows["CCCCCCCCCCCCCCCC-1"]["nh_gt1_unique_gene_countedU_reads"]), 1)
        self.assertEqual(int(barcode_rows["AAAAAAAAAAAAAAAA-1"]["mitochondrial_reads"]), 1)
        self.assertEqual(int(barcode_rows["CCCCCCCCCCCCCCCC-1"]["rrna_reads"]), 1)
        rg_rows = read_tsv(one / "rg_summary.tsv")
        rg_by_id = {row["RG"]: row for row in rg_rows}
        self.assertEqual(rg_by_id["RG2"]["source_index"], "1")
        self.assertEqual(rg_by_id["RG3"]["source_index"], "1")
        self.assertEqual(int(rg_by_id["RG1"]["candidate_countedU_reads"]), 10)
        self.assertEqual(int(rg_by_id["RG2"]["candidate_countedU_reads"]), 5)
        self.assertEqual(int(rg_by_id["RG3"]["candidate_countedU_reads"]), 1)
        self.assertEqual(int(rg_by_id["RG1"]["candidate_matrix_molecules"]), 4)
        self.assertEqual(int(rg_by_id["RG2"]["candidate_matrix_molecules"]), 5)
        self.assertEqual(int(rg_by_id["RG3"]["candidate_matrix_molecules"]), 1)
        self.assertEqual(sum(int(row["all_records"]) for row in rg_rows), 25)
        self.assertEqual(sum(int(row["primary_mapped_reads"]) for row in rg_rows), 22)
        self.assertEqual(
            sum(
                int(row["nh_gt1_unique_gene_countedU_reads"])
                for row in rg_rows
            ),
            1,
        )
        self.assertEqual(
            sum(int(row["candidate_countedU_reads"]) for row in rg_rows),
            16,
        )
        self.assertEqual(sum(int(row["secondary_records"]) for row in rg_rows), 2)
        self.assertEqual(sum(int(row["supplementary_records"]) for row in rg_rows), 1)
        self.assertEqual(sum(int(row["qcfail_records"]) for row in rg_rows), 1)
        self.assertEqual(sum(int(row["bam_duplicate_flag_reads"]) for row in rg_rows), 1)
        self.assertEqual(sum(int(row["missing_CB"]) for row in rg_rows), 1)

        molecule_rows = read_tsv(one / "molecule_source_hash_bins.tsv")
        expected_bin = (
            min(
                stable_hash(name)
                for name in ("q01", "q02", "q04", "q11", "q12", "q18")
            )
            * 100
        ) >> 64
        cross_source = [
            row for row in molecule_rows
            if row["CB"] == "AAAAAAAAAAAAAAAA-1"
            and row["source_mask_hex"] == "0x0000000000000003"
        ]
        self.assertEqual(len(cross_source), 1)
        self.assertEqual(int(cross_source[0]["nested_min_hash_bin"]), expected_bin)
        self.assertEqual(int(cross_source[0]["candidate_matrix_molecules"]), 1)

        conflict_rows = [
            row for row in read_tsv(one / "raw_to_corrected_barcode_counts.tsv")
            if row["CR"] == "TTTTTTTTTTTTTTTT"
        ]
        self.assertEqual({row["CB"] for row in conflict_rows}, {
            "AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"
        })
        self.assertTrue(all(row["within_rg_conflict"] == "0" for row in conflict_rows))
        self.assertTrue(all(row["cross_rg_conflict"] == "1" for row in conflict_rows))
        self.assertTrue(all(row["global_conflict"] == "1" for row in conflict_rows))
        within_rows = [
            row for row in read_tsv(one / "raw_to_corrected_barcode_counts.tsv")
            if row["CR"] == "GGGGGGGGGGGGGGGG"
        ]
        self.assertEqual({row["RG"] for row in within_rows}, {"RG1", "RG2"})
        self.assertTrue(all(
            row["within_rg_conflict"] == ("1" if row["RG"] == "RG1" else "0")
            for row in within_rows
        ))
        self.assertTrue(all(row["cross_rg_conflict"] == "0" for row in within_rows))
        self.assertTrue(all(row["global_conflict"] == "1" for row in within_rows))
        summary = runner.parse_profiler_summary(one / "profiler_summary.tsv")
        self.assertEqual(summary["starsolo_feature"], runner.STARSOLO_FEATURE)
        self.assertEqual(
            summary["starsolo_umi_filtering"], runner.STARSOLO_UMI_FILTERING
        )
        self.assertEqual(summary["starsolo_umi_dedup"], runner.STARSOLO_UMI_DEDUP)
        self.assertEqual(
            summary["starsolo_multimappers"], runner.STARSOLO_MULTIMAPPERS
        )
        self.assertEqual(summary["candidate_countedU_reads"], "16")
        self.assertEqual(summary["candidate_matrix_molecules"], "9")
        self.assertEqual(summary["nh_gt1_unique_gene_countedU_reads"], "1")
        self.assertEqual(summary["nh1_logical_representative_records"], "19")
        self.assertEqual(summary["nh_gt1_logical_read_states"], "2")
        self.assertEqual(summary["nh_gt1_alignment_records_observed"], "4")
        self.assertEqual(summary["nh_gt1_singleton_feature_logical_reads"], "1")
        self.assertEqual(summary["nh_gt1_countedU_logical_reads"], "1")
        self.assertEqual(
            summary["ordinary_countedU_read_definition"],
            runner.ORDINARY_COUNTEDU_READ_DEFINITION,
        )
        self.assertEqual(
            summary["ordinary_matrix_molecule_definition"],
            runner.ORDINARY_MOLECULE_DEFINITION,
        )
        self.assertEqual(
            summary["starsolo_EM_evidence_availability"],
            runner.STARSOLO_EM_EVIDENCE_AVAILABILITY,
        )
        self.assertEqual(
            summary["ordinary_raw_matrix_exact_coordinate_reconciliation"],
            "PASS",
        )
        self.assertEqual(
            summary["ordinary_filtered_matrix_exact_coordinate_reconciliation"],
            "PASS",
        )
        self.assertEqual(summary["ordinary_raw_matrix_molecules"], "9")
        self.assertEqual(summary["ordinary_filtered_matrix_molecules"], "9")
        self.assertEqual(summary["logical_record_loops"], "1")
        self.assertEqual(
            summary["record_contribution_method"],
            "single_tag_decode_and_single_cigar_walk_reused_across_accumulators",
        )
        self.assertEqual(summary["biological_classification_status"], "complete")
        self.assertGreater(int(summary["process_MaxRSS_kib"]), 0)

    def test_missing_classification_is_explicitly_unavailable(self) -> None:
        output = Path(self.temporary.name) / "without_class_manifest"
        subprocess.run(
            self.profiler_command(
                output, 1, include_class_manifest=False
            ),
            check=True,
        )
        summary = runner.parse_profiler_summary(output / "profiler_summary.tsv")
        self.assertEqual(summary["biological_classification_status"], "unavailable")
        for row in read_tsv(output / "barcode_read_metrics.tsv"):
            self.assertEqual(row["biological_classification_status"], "unavailable")
            self.assertEqual(row["mitochondrial_reads"], "")
            self.assertEqual(row["rrna_reads"], "")

    def test_em_matrix_cannot_be_promoted_as_ordinary(self) -> None:
        output = Path(self.temporary.name) / "em_as_ordinary_must_fail"
        command = self.profiler_command(output, 1)
        matrix_index = command.index("--raw-matrix") + 1
        command[matrix_index] = str(HERE / "UniqueAndMult-EM.mtx")
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ordinary raw matrix reconciliation failed", result.stderr)

    def test_mismatched_starsolo_semantics_are_rejected_before_scan(self) -> None:
        output = Path(self.temporary.name) / "wrong_starsolo_semantics"
        command = self.profiler_command(output, 1)
        option_index = command.index("--starsolo-umi-dedup") + 1
        command[option_index] = "Exact"
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("semantically bound to STARsolo", result.stderr)

    def test_runner_reconciliation_publication_and_gather(self) -> None:
        manifest = Path(self.temporary.name) / "manifest.tsv"
        fields = [
            "library", "library_dir", "bam", "bam_index", "summary", "raw_barcodes",
            "filtered_barcodes", "raw_features", "raw_matrix", "filtered_matrix",
            "rg_metadata", "source_order", "scientific_config",
            "source_order_provenance", "old_raw_barcodes",
            "starsolo_feature", "starsolo_umi_filtering", "starsolo_umi_dedup",
            "starsolo_multimappers", "biological_classification_intent",
            "old_filtered_barcodes", "class_manifest", "native_cell_reads",
            "output_dir",
        ]
        values = {
            "library": "FixtureRNA", "library_dir": self.library, "bam": self.bam,
            "bam_index": str(self.bam) + ".bai", "summary": self.library / "Summary.csv",
            "raw_barcodes": self.library / "raw" / "barcodes.tsv.gz",
            "filtered_barcodes": self.library / "filtered" / "barcodes.tsv.gz",
            "raw_features": self.library / "raw" / "features.tsv.gz",
            "raw_matrix": self.library / "raw" / "matrix.mtx.gz",
            "filtered_matrix": self.library / "filtered" / "matrix.mtx.gz",
            "rg_metadata": HERE / "rg_metadata.tsv", "source_order": HERE / "source_order.tsv",
            "source_order_provenance": "synthetic_fixture",
            "starsolo_feature": runner.STARSOLO_FEATURE,
            "starsolo_umi_filtering": runner.STARSOLO_UMI_FILTERING,
            "starsolo_umi_dedup": runner.STARSOLO_UMI_DEDUP,
            "starsolo_multimappers": runner.STARSOLO_MULTIMAPPERS,
            "biological_classification_intent": "manifest",
            "old_raw_barcodes": HERE / "barcodes.tsv",
            "old_filtered_barcodes": HERE / "historical_filtered_barcodes.tsv",
            "class_manifest": HERE / "contig_classes.tsv",
            "native_cell_reads": self.library / "CellReads.stats",
            "output_dir": self.library / "bam_evidence",
        }
        values["scientific_config"] = Path(self.temporary.name) / "scientific_config.json"
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerow(values)
        write_fixture_scientific_config(
            Path(values["scientific_config"]), fields, values, self.profiler
        )
        run_command = [
            sys.executable, str(MAPPING / "run_rna_bam_evidence.py"), "run",
            "--manifest", str(manifest), "--row-index", "0",
            "--profiler", str(self.profiler), "--threads", "2",
            "--max-memory-gb", "1.0",
        ]
        subprocess.run(run_command, check=True)
        evidence = self.library / "bam_evidence"
        self.assertTrue((evidence / "BAM_EVIDENCE_COMPLETE.ok").is_file())
        evidence_audit = json.loads(
            (evidence / "audit.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            evidence_audit["slurm_execution"],
            runner.slurm_execution_identity(),
        )
        original_audit = evidence_audit
        trial_environment = dict(os.environ)
        for name in (
            "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_ID"
        ):
            trial_environment.pop(name, None)
        trial_environment["SLURM_JOB_ID"] = "302"
        marker = evidence / "BAM_EVIDENCE_COMPLETE.ok"
        marker_bytes = marker.read_bytes()
        marker.write_bytes(b"")
        empty_marker = subprocess.run(
            run_command, capture_output=True, text=True, check=False,
            env=trial_environment,
        )
        self.assertNotEqual(empty_marker.returncode, 0)
        self.assertIn("missing or empty evidence completion marker", empty_marker.stderr)
        marker.write_bytes(marker_bytes)
        # A normal resume validates and preserves the audit produced by the
        # original allocation.
        subprocess.run(run_command, check=True, env=trial_environment)
        self.assertEqual(
            json.loads((evidence / "audit.json").read_text(encoding="utf-8")),
            original_audit,
        )
        self.assertEqual(
            sum(int(row["countedU"]) for row in read_tsv(evidence / "CellReads.countedU.from_bam.tsv.gz")),
            16,
        )
        gather = Path(self.temporary.name) / "gather"
        subprocess.run(
            [
                sys.executable, str(MAPPING / "run_rna_bam_evidence.py"), "gather",
                "--manifest", str(manifest), "--output-dir", str(gather),
            ],
            check=True,
        )
        self.assertTrue((gather / "BAM_EVIDENCE_GATHER_COMPLETE.ok").is_file())

        summary_path = self.library / "Summary.csv"
        summary_path.write_text(
            summary_path.read_text(encoding="utf-8").replace(
                "Unique Reads in Cells Mapped to GeneFull_Ex50pAS,16",
                "Unique Reads in Cells Mapped to GeneFull_Ex50pAS,17",
            ),
            encoding="utf-8",
        )
        failed_replacement = subprocess.run(
            [*run_command, "--replace-stale"],
            capture_output=True, text=True, check=False,
        )
        self.assertNotEqual(failed_replacement.returncode, 0)
        self.assertTrue((evidence / "BAM_EVIDENCE_FAILED.json").is_file())
        self.assertFalse((evidence / "BAM_EVIDENCE_COMPLETE.ok").exists())
        self.assertFalse(
            (evidence / "CellReads.countedU.from_bam.tsv.gz").exists()
        )


if __name__ == "__main__":
    unittest.main()

