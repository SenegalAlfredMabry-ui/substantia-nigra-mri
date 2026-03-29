#!/usr/bin/env python3
"""Attach behavioral covariates to Irredeemable v9 tables and write *_with_covariates copies."""

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Sequence


DEFAULT_MANIFEST = "behavioral/uber_data_20250728.csv"
DEFAULT_ROOT = "MRI_data/analysis/tests/v9"
DEFAULT_PATTERN = "SN_nmdata_irredeemable_v9*.txt"
OUTPUT_SUFFIX = "_with_covariates"


def normalize_subject(raw: str) -> str:
    if raw is None:
        return ""
    text = str(raw).strip()
    if not text:
        return ""
    if text.upper().startswith("HUMAN_EBR-"):
        text = text.split("-", 1)[1]
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def load_manifest(path: Path, requested_cols: Optional[Sequence[str]]) -> tuple[Dict[str, Dict[str, str]], List[str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"Manifest {path} has no header")
        header = reader.fieldnames
        if not set(col.lower() for col in header) & {"sub", "subject"}:
            raise SystemExit(f"Manifest {path} is missing a subject/sub column")
        if requested_cols:
            columns = list(requested_cols)
        else:
            columns = [col for col in header if col.lower() not in ("sub", "subject")]
        roster: Dict[str, Dict[str, str]] = {}
        for row in reader:
            key = normalize_subject(row.get("sub") or row.get("subject"))
            if not key:
                continue
            roster[key] = {col: row.get(col, "") for col in columns}
    return roster, columns


def determine_tables(explicit: Sequence[str], root: Path, pattern: str) -> List[Path]:
    tables: List[Path] = []
    if explicit:
        for item in explicit:
            candidate = Path(item)
            if candidate.is_file():
                tables.append(candidate)
            else:
                print(f"[WARN] Skipping missing table {candidate}")
    else:
        if not root.is_dir():
            raise SystemExit(f"Root directory {root} does not exist")
        tables.extend(sorted(root.rglob(pattern)))
    return tables


def output_path_for(table: Path, suffix: str) -> Path:
    name = table.stem
    if name.endswith(suffix):
        return table
    return table.with_name(f"{name}{suffix}{table.suffix}")


def merge_covariates(table: Path, roster: Dict[str, Dict[str, str]], columns: Sequence[str], force: bool) -> Optional[Path]:
    output = output_path_for(table, OUTPUT_SUFFIX)
    if output == table and not force:
        print(f"[INFO] {table} already contains covariates suffix; use --force to overwrite")
        return None

    with table.open() as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if reader.fieldnames is None:
            print(f"[WARN] {table} has no header; skipping")
            return None
        header = reader.fieldnames

    missing = [col for col in columns if col not in header]
    merged_header = header + missing

    for row in rows:
        key = normalize_subject(row.get("sub") or row.get("subject"))
        extras = roster.get(key, {})
        for col in columns:
            row[col] = extras.get(col, row.get(col, ""))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=merged_header)
        writer.writeheader()
        writer.writerows(rows)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "tables",
        nargs="*",
        help="Specific Irredeemable tables to process (default: auto-discover under --root)",
    )
    parser.add_argument(
        "--manifest",
        default=DEFAULT_MANIFEST,
        help="Behavioral manifest CSV containing covariates (default: %(default)s)",
    )
    parser.add_argument(
        "--columns",
        nargs="*",
        help="Subset of manifest columns to merge (default: all columns except subject/sub)",
    )
    parser.add_argument(
        "--root",
        default=DEFAULT_ROOT,
        help="Root directory to search when tables are not explicitly provided (default: %(default)s)",
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help="Glob pattern for discovery when --tables is empty (default: %(default)s)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite files that already look like *_with_covariates",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest)
    if not manifest_path.is_file():
        raise SystemExit(f"Manifest not found at {manifest_path}")
    roster, columns = load_manifest(manifest_path, args.columns)
    tables = determine_tables(args.tables, Path(args.root), args.pattern)
    if not tables:
        print("[WARN] No Irredeemable tables matched the request")
        return
    updated = 0
    for table in tables:
        result = merge_covariates(table, roster, columns, force=args.force)
        if result:
            updated += 1
            print(f"[INFO] Wrote {result}")
    print(f"Merged covariates into {updated} table(s)")


if __name__ == "__main__":
    main()
