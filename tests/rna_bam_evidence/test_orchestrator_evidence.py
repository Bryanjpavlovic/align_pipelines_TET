#!/usr/bin/env python3
"""Focused tests for direct all-library RNA BAM evidence harvesting."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
MAPPING = ROOT / "scripts" / "mapping"
sys.path.insert(0, str(MAPPING))

import orchestrate_10x_mapping_qc as orchestrator


class DirectEvidenceGenerationTest(unittest.TestCase):
    def resources(self) -> orchestrator.Resources:
        return orchestrator.Resources(
            rna_driver=None,
            atac_driver=None,
            rna_workflow=None,
            atac_workflow=None,
            barcode_aggregator=None,
            trim_plotter=None,
            mapping_plotter=None,
            atac_collector=None,
            atac_plotter=None,
            bam_evidence_runner=(MAPPING / "run_rna_bam_evidence.py").resolve(),
            bam_evidence_profiler=Path(sys.executable).resolve(),
            star_diagnostic_promoter=(
                MAPPING / "promote_rna_star_diagnostics.py"
            ).resolve(),
        )

    def args(self, **overrides: object) -> argparse.Namespace:
        values: dict[str, object] = {
            "rna3_bam_evidence_source_order": ["BP1", "BP2"],
            "rna3_bam_evidence_baseline_root": None,
            "rna3_bam_evidence_class_manifest": None,
            "rna3_bam_evidence_memory_gb": 48,
            "rna3_bam_evidence_hash_bins": 100,
            "rna3_bam_evidence_cpus": 4,
            "rna3_bam_evidence_max_concurrent": 3,
            "rna3_bam_evidence_reset_failed_run": False,
        }
        values.update(overrides)
        return argparse.Namespace(**values)

    def create_libraries(self, root: Path, numbers: list[int]) -> list[str]:
        libraries: list[str] = []
        for number in numbers:
            library = f"Tet_2025_Multiome-RNA_{number}"
            libraries.append(library)
            directory = root / "mapping_output" / library
            directory.mkdir(parents=True)
            (directory / "gex.bam").write_bytes(b"bam")
        return libraries

    def generate(
        self,
        run_dir: Path,
        libraries: list[str],
        args: argparse.Namespace | None = None,
    ) -> list[orchestrator.JobSpec]:
        return orchestrator.generate_rna_bam_evidence_jobs(
            libraries,
            run_dir / "rna3",
            run_dir,
            args or self.args(),
            self.resources(),
            [],
            ["BP1", "BP2"],
        )

    def test_one_invocation_generates_the_complete_bounded_dag(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            libraries = self.create_libraries(run_dir / "rna3", [29, 2, 20, 40])
            jobs = self.generate(run_dir, libraries)

            self.assertEqual(
                [job.label for job in jobs],
                [
                    "rna3_star_diagnostics",
                    "rna3_bam_evidence",
                    "rna3_bam_evidence_gather",
                ],
            )
            self.assertEqual(
                jobs[-1].dependencies,
                ["rna3_star_diagnostics", "rna3_bam_evidence"],
            )
            array = (
                run_dir / "control" / "slurm" / "rna3_bam_evidence.sbatch"
            ).read_text(encoding="utf-8")
            self.assertIn("#SBATCH --array=0-3%3", array)
            self.assertNotIn("force-trial-rescan", array)
            self.assertFalse(
                (run_dir / "control" / "slurm" /
                 "rna3_bam_evidence_pilot.sbatch").exists()
            )
            self.assertFalse(
                (run_dir / "control" / "slurm" /
                 "rna3_bam_evidence_two_reader_benchmark.sbatch").exists()
            )
            execution = json.loads(
                (run_dir / "control" /
                 "rna3_bam_evidence_execution.json").read_text(encoding="utf-8")
            )
            self.assertEqual(execution["mode"], "direct_full_run")
            self.assertEqual(execution["selected_max_concurrent"], 3)
            self.assertEqual(
                execution["selected_libraries"],
                [
                    "Tet_2025_Multiome-RNA_2",
                    "Tet_2025_Multiome-RNA_20",
                    "Tet_2025_Multiome-RNA_29",
                    "Tet_2025_Multiome-RNA_40",
                ],
            )

    def test_reader_jobs_use_only_required_modules(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            libraries = self.create_libraries(run_dir / "rna3", [1, 2])
            self.generate(run_dir, libraries)
            script = (
                run_dir / "control" / "slurm" / "rna3_bam_evidence.sbatch"
            ).read_text(encoding="utf-8")
            self.assertIn("module purge", script)
            self.assertIn(
                "module load miniforge/3 htslib/1.20 samtools/1.20", script
            )
            self.assertNotIn("nextflow/latest", script)
            self.assertNotIn("align_pipelines/bjp", script)

    def test_value_affecting_configuration_remains_frozen(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            libraries = self.create_libraries(run_dir / "rna3", [1, 2])
            self.generate(run_dir, libraries)

            allowed = self.args(
                rna3_bam_evidence_cpus=8,
                rna3_bam_evidence_memory_gb=64,
                rna3_bam_evidence_max_concurrent=1,
            )
            self.generate(run_dir, libraries, allowed)
            array = (
                run_dir / "control" / "slurm" / "rna3_bam_evidence.sbatch"
            ).read_text(encoding="utf-8")
            self.assertIn("#SBATCH --array=0-1%1", array)
            self.assertIn("#SBATCH --cpus-per-task=8", array)
            self.assertIn("#SBATCH --mem=64G", array)

            rejected = (
                self.args(rna3_bam_evidence_hash_bins=101),
                self.args(rna3_bam_evidence_source_order=["BP2", "BP1"]),
            )
            for changed in rejected:
                with self.subTest(changed=vars(changed)):
                    with self.assertRaisesRegex(
                        orchestrator.OrchestratorError,
                        "scientific configuration differs",
                    ):
                        self.generate(run_dir, libraries, changed)

    def test_classification_change_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            libraries = self.create_libraries(run_dir / "rna3", [1])
            classes = run_dir / "classes.tsv"
            classes.write_text(
                "contig\tmitochondrial\nchr1\t0\n", encoding="utf-8"
            )
            args = self.args(rna3_bam_evidence_class_manifest=str(classes))
            self.generate(run_dir, libraries, args)
            classes.write_text(
                "contig\tmitochondrial\nchr1\t1\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(
                orchestrator.OrchestratorError,
                "scientific configuration differs",
            ):
                self.generate(run_dir, libraries, args)

    def test_failed_run_reset_archives_partial_outputs_and_regenerates_all(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            libraries = self.create_libraries(run_dir / "rna3", [1, 2])
            self.generate(run_dir, libraries)
            for library in libraries:
                evidence = (
                    run_dir / "rna3" / "mapping_output" / library / "bam_evidence"
                )
                evidence.mkdir(parents=True)
                (evidence / "BAM_EVIDENCE_COMPLETE.ok").write_text(
                    "obsolete\n", encoding="utf-8"
                )
            gathered = run_dir / "qc" / "rna3" / "bam_evidence"
            gathered.mkdir(parents=True)
            (gathered / "old.tsv").write_text("old\n", encoding="utf-8")
            orchestrator.atomic_json(
                run_dir / "control" / "job_plan.json",
                {
                    "jobs": [
                        {
                            "label": "rna3_star_diagnostics",
                            "job_id": "100",
                            "dependencies": [],
                        },
                        {
                            "label": "rna3_bam_evidence",
                            "job_id": "101",
                            "dependencies": [],
                        },
                        {
                            "label": "rna3_bam_evidence_gather",
                            "job_id": "102",
                            "dependencies": [
                                "rna3_star_diagnostics", "rna3_bam_evidence"
                            ],
                        },
                        {
                            "label": "rna3_trim_barcode_join",
                            "job_id": "103",
                            "dependencies": ["rna3_bam_evidence_gather"],
                        },
                        {
                            "label": "rna3_map",
                            "job_id": "99",
                            "dependencies": [],
                        },
                    ]
                },
            )

            with mock.patch.object(
                orchestrator,
                "slurm_accounting_states",
                return_value={"FAILED"},
            ):
                jobs = self.generate(
                    run_dir,
                    libraries,
                    self.args(rna3_bam_evidence_reset_failed_run=True),
                )

            self.assertEqual(
                [job.label for job in jobs],
                [
                    "rna3_star_diagnostics",
                    "rna3_bam_evidence",
                    "rna3_bam_evidence_gather",
                ],
            )
            migrations = list(
                (run_dir / "control" /
                 "rna3_bam_evidence_failed_run_migrations").iterdir()
            )
            self.assertEqual(len(migrations), 1)
            self.assertTrue(
                (migrations[0] / "evidence" / libraries[0] /
                 "BAM_EVIDENCE_COMPLETE.ok").is_file()
            )
            self.assertTrue(
                (migrations[0] / "gathered_bam_evidence" / "old.tsv").is_file()
            )
            self.assertFalse(
                (run_dir / "rna3" / "mapping_output" / libraries[0] /
                 "bam_evidence").exists()
            )
            retained = json.loads(
                (run_dir / "control" / "job_plan.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                [job["label"] for job in retained["jobs"]], ["rna3_map"]
            )

    def test_failed_run_reset_refuses_active_reader(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            libraries = self.create_libraries(run_dir / "rna3", [1])
            self.generate(run_dir, libraries)
            orchestrator.atomic_json(
                run_dir / "control" / "job_plan.json",
                {
                    "jobs": [
                        {
                            "label": "rna3_bam_evidence",
                            "job_id": "101",
                            "dependencies": [],
                        }
                    ]
                },
            )
            with mock.patch.object(
                orchestrator,
                "slurm_accounting_states",
                return_value={"RUNNING"},
            ):
                with self.assertRaisesRegex(
                    orchestrator.OrchestratorError,
                    "while prior RNA evidence job",
                ):
                    self.generate(
                        run_dir,
                        libraries,
                        self.args(rna3_bam_evidence_reset_failed_run=True),
                    )

    def test_mapping_configuration_immutability_is_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "existing"
            payload = {
                "release": "old",
                "inputs": {"rna3": ["/data/BP1"]},
                "references": {"rna": "/ref/a"},
                "resources": {},
                "libraries": {"rna3": ["lib1"]},
                "input_fastqs": {"rna3": []},
                "rna5_formats": {},
                "stages": ["map"],
                "options": {"rna_threads": 8},
                "scheduling": None,
                "reporting": {"rna3_bam_evidence_from_bam": False},
            }
            orchestrator.prepare_run_directory(run_dir, False, payload, None)
            enabled = json.loads(json.dumps(payload))
            enabled["reporting"] = {"rna3_bam_evidence_from_bam": True}
            orchestrator.prepare_run_directory(run_dir, True, enabled, None)

            changed = json.loads(json.dumps(enabled))
            changed["references"]["rna"] = "/ref/b"
            with self.assertRaisesRegex(
                orchestrator.OrchestratorError,
                "does not match the existing run",
            ):
                orchestrator.prepare_run_directory(run_dir, True, changed, None)

    def test_phase_and_continuation_options_are_removed(self) -> None:
        help_text = orchestrator.build_parser().format_help()
        self.assertNotIn("--rna3-bam-evidence-phase", help_text)
        self.assertNotIn("--rna3-bam-evidence-continue-audit", help_text)
        self.assertIn("--rna3-bam-evidence-reset-failed-run", help_text)


if __name__ == "__main__":
    unittest.main()
