#!/usr/bin/env python3
"""Small, dependency-free readers shared by RNA evidence analyses."""

from __future__ import annotations

import csv
import gzip
import mmap
import os
import struct
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, TextIO


class AnalysisError(RuntimeError):
    pass


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield handle


def barcodes(path: Path) -> list[str]:
    if not path.is_file():
        raise AnalysisError(f"missing barcode roster: {path}")
    with open_text(path) as handle:
        result = [line.rstrip("\r\n").split("\t", 1)[0] for line in handle]
    result = [value for value in result if value]
    if not result or len(result) != len(set(result)):
        raise AnalysisError(f"empty or duplicate barcode roster: {path}")
    return result


def iter_nonzero_matrix_barcodes(
    matrix: Path, roster_path: Path
) -> Iterator[tuple[str, int, int]]:
    """Yield nonzero columns using a disk-backed fixed-width numeric vector.

    The raw barcode strings are streamed after the MatrixMarket pass. Memory is
    independent of the raw roster cardinality at Python-object level. The
    12-byte-per-column payload is an mmap-backed temporary file, rather than a
    Python list, string set, or nested dictionary. Zero-count whitelist entries
    are never yielded.
    """
    with open_text(matrix) as handle:
        banner = handle.readline().lower()
        if not banner.startswith("%%matrixmarket matrix coordinate"):
            raise AnalysisError(f"unsupported MatrixMarket file: {matrix}")
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        try:
            feature_count, column_count, declared = map(int, line.split())
        except ValueError as exc:
            raise AnalysisError(f"malformed MatrixMarket dimensions: {matrix}") from exc
        if column_count < 1:
            raise AnalysisError(f"matrix has no barcode columns: {matrix}")
        record = struct.Struct("<QI")
        scratch = os.environ.get("SLURM_TMPDIR")
        temporary_dir = scratch if scratch and Path(scratch).is_dir() else None
        try:
            backing_context = tempfile.TemporaryFile(
                prefix="rna_matrix_columns_", dir=temporary_dir
            )
        except OSError as exc:
            raise AnalysisError(
                f"could not create matrix column backing file: {exc}"
            ) from exc
        with backing_context as backing:
            try:
                backing.truncate(record.size * column_count)
                column_values = mmap.mmap(
                    backing.fileno(), record.size * column_count,
                    access=mmap.ACCESS_WRITE,
                )
            except (OSError, ValueError) as exc:
                raise AnalysisError(
                    f"could not allocate disk-backed matrix columns: {exc}"
                ) from exc
            with column_values:
                if hasattr(column_values, "madvise") and hasattr(mmap, "MADV_RANDOM"):
                    try:
                        column_values.madvise(mmap.MADV_RANDOM)
                    except OSError:
                        pass
                observed = 0
                previous_column = 0
                previous_feature = 0
                for line_number, line in enumerate(handle, start=3):
                    if not line.strip() or line.startswith("%"):
                        continue
                    try:
                        feature, column, count = map(int, line.split())
                    except ValueError as exc:
                        raise AnalysisError(
                            f"malformed matrix row {line_number}: {matrix}"
                        ) from exc
                    if (
                        not 1 <= feature <= feature_count
                        or not 1 <= column <= column_count
                        or count <= 0
                    ):
                        raise AnalysisError(
                            f"invalid matrix row {line_number}: {matrix}"
                        )
                    if (
                        column < previous_column
                        or (
                            column == previous_column
                            and feature <= previous_feature
                        )
                    ):
                        raise AnalysisError(
                            "matrix coordinates must be strictly ordered by "
                            f"barcode then feature: {matrix}:{line_number}"
                        )
                    previous_column = column
                    previous_feature = feature
                    offset = (column - 1) * record.size
                    umi, genes = record.unpack_from(column_values, offset)
                    if umi > (1 << 64) - 1 - count or genes == (1 << 32) - 1:
                        raise AnalysisError(
                            f"matrix column counter overflow at {matrix}:{line_number}"
                        )
                    record.pack_into(
                        column_values, offset, umi + count, genes + 1
                    )
                    observed += 1
                if observed != declared:
                    raise AnalysisError(f"MatrixMarket nnz mismatch: {matrix}")

                if (
                    hasattr(column_values, "madvise")
                    and hasattr(mmap, "MADV_SEQUENTIAL")
                ):
                    try:
                        column_values.madvise(mmap.MADV_SEQUENTIAL)
                    except OSError:
                        pass
                roster_rows = 0
                with open_text(roster_path) as roster_handle:
                    for roster_rows, roster_line in enumerate(
                        roster_handle, start=1
                    ):
                        barcode = roster_line.rstrip("\r\n").split("\t", 1)[0]
                        if not barcode:
                            raise AnalysisError(
                                f"blank barcode in {roster_path}:{roster_rows}"
                            )
                        if roster_rows > column_count:
                            raise AnalysisError(
                                f"matrix/barcode dimension mismatch: {matrix}"
                            )
                        offset = (roster_rows - 1) * record.size
                        umi, genes = record.unpack_from(column_values, offset)
                        if umi:
                            yield barcode, int(umi), int(genes)
                if roster_rows != column_count:
                    raise AnalysisError(
                        f"matrix/barcode dimension mismatch: {matrix}"
                    )


def find_product(directory: Path, relative: str) -> Path:
    direct = directory / relative
    if direct.is_file():
        return direct
    if direct.suffix == ".gz":
        plain = direct.with_suffix("")
        if plain.is_file():
            return plain
    else:
        compressed = direct.with_name(direct.name + ".gz")
        if compressed.is_file():
            return compressed
    raise AnalysisError(f"missing product below {directory}: {relative}")


def library_dirs(root: Path) -> dict[str, Path]:
    if not root.is_dir():
        raise AnalysisError(f"mapping root does not exist: {root}")
    result = {
        path.name: path
        for path in root.iterdir()
        if path.is_dir() and (path / "Summary.csv").is_file()
    }
    if not result:
        raise AnalysisError(f"no STARsolo library directories below {root}")
    return result


def read_tsv(path: Path, required: Iterable[str] = ()) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise AnalysisError(f"{path} lacks columns: {', '.join(sorted(missing))}")
        return list(reader)


def iter_tsv(path: Path, required: Iterable[str] = ()) -> Iterator[dict[str, str]]:
    """Stream a TSV while keeping the owning text handle open."""
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise AnalysisError(f"{path} lacks columns: {', '.join(sorted(missing))}")
        for row in reader:
            yield row


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]], *, compressed: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    opener = gzip.open if compressed else open
    with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)
