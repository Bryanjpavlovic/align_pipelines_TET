#!/usr/bin/env python3
"""Safely promote small STAR diagnostics from a retained Nextflow work cache.

This independent step never opens a BAM and never deletes the work directory.
It matches a work task to each published library using the exact Summary.csv
SHA-256, then atomically publishes STAR_Log.out, STAR_Log.final.out, and
STAR_SJ.out.tab.gz.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


class PromotionError(RuntimeError):
    pass


RELEASE = "2026-09-05-star-diagnostics-v2-audited-reuse"
OUTPUT_NAMES = (
    "STAR_Log.out", "STAR_Log.final.out", "STAR_SJ.out.tab.gz"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def nonempty(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise PromotionError(f"missing or empty {label}: {path}")
    return path


def work_sources(task: Path, library: str) -> dict[str, Path] | None:
    choices = {
        "STAR_Log.out": (task / f"{library}Log.out", task / "STAR_Log.out"),
        "STAR_Log.final.out": (
            task / f"{library}Log.final.out", task / "STAR_Log.final.out"
        ),
        "STAR_SJ.out.tab.gz": (
            task / f"{library}SJ.out.tab",
            task / "STAR_SJ.out.tab",
            task / "STAR_SJ.out.tab.gz",
        ),
    }
    result: dict[str, Path] = {}
    for destination, candidates in choices.items():
        source = next(
            (
                path for path in candidates
                if path.is_file() and path.stat().st_size > 0
            ),
            None,
        )
        if source is None:
            return None
        result[destination] = source
    return result


def publish(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.tmp.{os.getpid()}"
    )
    try:
        if destination.suffix == ".gz" and source.suffix != ".gz":
            with source.open("rb") as incoming, temporary.open("wb") as raw:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw, mtime=0
                ) as outgoing:
                    shutil.copyfileobj(incoming, outgoing, length=1024 * 1024)
        else:
            shutil.copyfile(source, temporary)
        nonempty(temporary, "temporary STAR diagnostic")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def validate_destination(path: Path) -> dict[str, object]:
    nonempty(path, "published STAR diagnostic")
    if path.suffix == ".gz":
        uncompressed = 0
        with gzip.open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                uncompressed += len(block)
        if not uncompressed:
            raise PromotionError(f"gzip diagnostic has no content: {path}")
    else:
        uncompressed = path.stat().st_size
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "uncompressed_bytes": uncompressed,
        "sha256": digest(path),
    }


def manifest_rows(manifest: Path) -> list[dict[str, str]]:
    nonempty(manifest, "RNA evidence manifest")
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"library", "library_dir", "summary"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise PromotionError(
                f"manifest lacks columns: {', '.join(sorted(missing))}"
            )
        rows = list(reader)
    libraries = [row["library"] for row in rows]
    if (
        not rows
        or any(not library for library in libraries)
        or len(libraries) != len(set(libraries))
    ):
        raise PromotionError("manifest contains no libraries or duplicate/blank names")
    return rows


def validate_prior_promotion(
    manifest: Path,
    work_root: Path,
    audit_path: Path,
    marker: Path,
) -> dict[str, object]:
    """Revalidate a prior promotion against current inputs and every output."""
    try:
        payload = json.loads(nonempty(
            audit_path, "STAR diagnostic promotion audit"
        ).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PromotionError(f"could not read promotion audit: {audit_path}") from exc
    if not isinstance(payload, dict):
        raise PromotionError("STAR diagnostic promotion audit is not an object")
    try:
        marker_payload = json.loads(nonempty(
            marker, "STAR diagnostic promotion marker"
        ).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PromotionError(
            f"could not read promotion marker: {marker}"
        ) from exc
    rows = manifest_rows(manifest)
    expected_libraries = {row["library"] for row in rows}
    records = payload.get("libraries")
    if (
        payload.get("release") != RELEASE
        or payload.get("status") != "PASS"
        or payload.get("manifest") != str(manifest.resolve())
        or payload.get("manifest_sha256") != digest(manifest)
        or payload.get("work_root") != str(work_root.resolve())
        or payload.get("bam_opened") is not False
        or payload.get("retired_awk_scanner_used") is not False
        or payload.get("work_directory_deleted") is not False
        or not isinstance(records, dict)
        or set(records) != expected_libraries
    ):
        raise PromotionError(
            "STAR diagnostic promotion audit is stale or belongs to another run"
        )
    expected_marker = {
        "release": RELEASE,
        "status": "PASS",
        "audit": str(audit_path.resolve()),
        "audit_sha256": digest(audit_path),
        "manifest": str(manifest.resolve()),
        "manifest_sha256": digest(manifest),
    }
    if marker_payload != expected_marker:
        raise PromotionError(
            "STAR diagnostic promotion marker is stale or altered"
        )
    for row in rows:
        library = row["library"]
        record = records[library]
        summary = nonempty(Path(row["summary"]), f"{library} Summary.csv")
        if not isinstance(record, dict):
            raise PromotionError(
                f"STAR diagnostic audit has an invalid record for {library}"
            )
        outputs = record.get("outputs")
        if (
            record.get("summary_sha256") != digest(summary)
            or not isinstance(outputs, dict)
            or set(outputs) != set(OUTPUT_NAMES)
        ):
            raise PromotionError(
                f"STAR diagnostic audit is stale or incomplete for {library}"
            )
        library_dir = Path(row["library_dir"])
        for name in OUTPUT_NAMES:
            actual = validate_destination(library_dir / name)
            if outputs[name] != actual:
                raise PromotionError(
                    f"published STAR diagnostic changed for {library}: {name}"
                )
    return payload


def run(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).expanduser().resolve()
    work_root = Path(args.work_root).expanduser().resolve()
    audit_path = Path(args.audit).expanduser().resolve(strict=False)
    rows = manifest_rows(manifest)
    marker = (
        Path(args.marker).expanduser().resolve(strict=False)
        if args.marker else audit_path.with_name("STAR_DIAGNOSTIC_PROMOTION_COMPLETE.ok")
    )
    if audit_path.is_file() and marker.is_file():
        try:
            validate_prior_promotion(manifest, work_root, audit_path, marker)
        except PromotionError:
            marker.unlink(missing_ok=True)
        else:
            print(f"Revalidated STAR diagnostics for {len(rows)} libraries")
            return
    else:
        marker.unlink(missing_ok=True)
    if not work_root.is_dir():
        raise PromotionError(f"Nextflow work root does not exist: {work_root}")

    needs: dict[str, tuple[Path, str]] = {}
    records: dict[str, dict[str, object]] = {}
    for row in rows:
        library = row["library"]
        library_dir = Path(row["library_dir"])
        summary = nonempty(Path(row["summary"]), f"{library} Summary.csv")
        needs[library] = (library_dir, digest(summary))

    by_summary: dict[str, list[Path]] = {}
    if needs:
        wanted = {value[1] for value in needs.values()}
        for summary in work_root.glob("*/*/Summary.csv"):
            if not summary.is_file() or summary.stat().st_size == 0:
                continue
            summary_digest = digest(summary)
            if summary_digest in wanted:
                by_summary.setdefault(summary_digest, []).append(
                    summary.parent
                )

    for library, (library_dir, summary_digest) in sorted(needs.items()):
        candidates = [
            task for task in by_summary.get(summary_digest, [])
            if work_sources(task, library) is not None
        ]
        if not candidates:
            raise PromotionError(
                f"no retained successful STAR task matches {library} "
                f"Summary.csv and all three diagnostics"
            )
        if len(candidates) != 1:
            raise PromotionError(
                f"ambiguous retained STAR tasks for {library}: "
                + ", ".join(str(path) for path in sorted(candidates))
            )
        task = candidates[0].resolve()
        sources = work_sources(task, library)
        assert sources is not None
        for name, source in sources.items():
            publish(source, library_dir / name)
        records[library] = {
            "status": "verified_and_published_from_unique_work_task",
            "summary_sha256": summary_digest,
            "work_task": str(task),
            "candidate_task_count": len(candidates),
            "source_files": {
                name: {
                    "path": str(source),
                    "bytes": source.stat().st_size,
                    "sha256": digest(source),
                }
                for name, source in sources.items()
            },
            "outputs": {
                name: validate_destination(library_dir / name)
                for name in sources
            },
        }

    payload = {
        "release": RELEASE,
        "created_utc": utc_now(),
        "status": "PASS",
        "manifest": str(manifest),
        "manifest_sha256": digest(manifest),
        "work_root": str(work_root),
        "bam_opened": False,
        "retired_awk_scanner_used": False,
        "work_directory_deleted": False,
        "libraries": records,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = audit_path.with_name(f".{audit_path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, audit_path)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker_temporary = marker.with_name(f".{marker.name}.tmp.{os.getpid()}")
    try:
        marker_temporary.write_text(
            json.dumps(
                {
                    "release": RELEASE,
                    "status": "PASS",
                    "audit": str(audit_path.resolve()),
                    "audit_sha256": digest(audit_path),
                    "manifest": str(manifest.resolve()),
                    "manifest_sha256": digest(manifest),
                },
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        os.replace(marker_temporary, marker)
    finally:
        marker_temporary.unlink(missing_ok=True)
    print(f"Validated STAR diagnostics for {len(rows)} libraries")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--manifest", required=True)
    result.add_argument("--work-root", required=True)
    result.add_argument("--audit", required=True)
    result.add_argument("--marker")
    result.add_argument(
        "--validate-only", action="store_true",
        help="revalidate an existing audit, marker, manifest, and all outputs",
    )
    return result


if __name__ == "__main__":
    try:
        arguments = parser().parse_args()
        if arguments.validate_only:
            manifest = Path(arguments.manifest).expanduser().resolve()
            work_root = Path(arguments.work_root).expanduser().resolve()
            audit = Path(arguments.audit).expanduser().resolve()
            marker = (
                Path(arguments.marker).expanduser().resolve()
                if arguments.marker else audit.with_name(
                    "STAR_DIAGNOSTIC_PROMOTION_COMPLETE.ok"
                )
            )
            payload = validate_prior_promotion(
                manifest, work_root, audit, marker
            )
            print(
                "Validated existing STAR diagnostics for "
                f"{len(payload['libraries'])} libraries"
            )
        else:
            run(arguments)
    except (PromotionError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
