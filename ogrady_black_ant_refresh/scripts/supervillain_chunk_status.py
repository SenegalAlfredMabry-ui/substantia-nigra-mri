#!/usr/bin/env python3
"""Summarize supervillain warp progress per chunk."""

import argparse
import csv
import datetime as dt
from pathlib import Path as _Path
from typing import Any, Dict, List, Optional

DONE_STATUSES = {"success", "skipped_existing"}
CHUNKS: Dict[str, List[int]] = {
    "Chunk A": [117, 118, 120, 128, 130, 132, 135, 137, 140, 142, 146, 149, 151, 152, 154, 156,
                163, 165, 166, 170, 174, 201, 203, 204, 206, 207, 210, 212, 215, 217, 218, 219,
                221, 222, 223, 225, 226, 227, 228, 231],
    "Chunk B": [233, 234, 236, 237, 239, 240, 241, 242, 244, 246, 502, 503, 504, 505, 506, 508,
                510, 511, 512, 513, 515, 518, 519, 520, 522, 523, 524, 525, 526, 527, 528, 529,
                533, 534, 536, 537, 538, 539, 540, 541],
    "Chunk C": [542, 543, 544, 545, 548, 549, 551, 552, 553, 554, 555, 556, 558, 559, 560, 561,
                562, 563, 565, 566, 570, 571, 572, 574, 576, 578, 582, 584, 585, 601, 602, 604,
                605, 606, 608, 609, 610, 611, 612, 613],
    "Chunk D": [614, 530, 546, 564, 569, 575, 577, 580, 581, 596, 615, 616, 617, 619, 620, 621,
                628, 631, 632, 633, 634, 635, 637],
    "Chunk E": [217, 218, 219, 221, 222, 223, 225, 226, 227, 228, 231],
}

def parse_timestamp(value: Optional[str]) -> Optional[dt.datetime]:
    if not value:
        return None
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"
    try:
        return dt.datetime.fromisoformat(value)
    except ValueError:
        return None

def load_manifest(path: _Path) -> Dict[int, Dict[str, Any]]:
    data: Dict[int, Dict[str, Any]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        for row in reader:
            try:
                subject = int(row["subject"])
            except (KeyError, ValueError):
                continue
            data[subject] = {
                "status": row.get("status", ""),
                "timestamp": parse_timestamp(row.get("timestamp")),
            }
    return data

def summarize_chunk(name: str, ids: List[int], manifest: Dict[int, Dict[str, Any]], now: dt.datetime) -> Dict[str, Any]:
    done: List[tuple[int, Dict[str, Any]]] = []
    pending: List[int] = []
    for subject in ids:
        record = manifest.get(subject)
        if record and record.get("status") in DONE_STATUSES:
            done.append((subject, record))
        else:
            pending.append(subject)
    done.sort(key=lambda item: item[1]["timestamp"] or dt.datetime.min)
    intervals = []
    for (_, first), (_, second) in zip(done, done[1:]):
        t0 = first["timestamp"]
        t1 = second["timestamp"]
        if t0 and t1:
            intervals.append((t1 - t0).total_seconds())
    avg_interval = sum(intervals) / len(intervals) if intervals else None
    last_finish = done[-1][1]["timestamp"] if done else None
    eta = None
    if not pending and last_finish:
        eta = last_finish
    elif pending and avg_interval and last_finish:
        eta = last_finish + dt.timedelta(seconds=avg_interval * len(pending))
    elif pending and avg_interval and not last_finish:
        eta = now + dt.timedelta(seconds=avg_interval * len(pending))
    return {
        "name": name,
        "total": len(ids),
        "done": len(done),
        "pending": pending,
        "avg_interval_hours": avg_interval / 3600 if avg_interval else None,
        "last_finish": last_finish,
        "eta": eta,
    }

def format_time(value: Optional[dt.datetime]) -> str:
    if not value:
        return "n/a"
    return value.replace(tzinfo=dt.timezone.utc).isoformat().replace("+00:00", "Z")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="MRI_data/WARPS_supervillain/warp_manifest.tsv",
        help="Path to the warp manifest TSV",
    )
    parser.add_argument(
        "--chunks",
        nargs="*",
        default=sorted(CHUNKS.keys()),
        help="Subset of chunk names to report (default: all)",
    )
    parser.add_argument(
        "--now",
        help="Override the current UTC time (ISO-8601)",
    )
    args = parser.parse_args()

    manifest_path = _Path(args.manifest)
    if not manifest_path.exists():
        parser.error(f"Manifest not found: {manifest_path}")

    manifest = load_manifest(manifest_path)
    if args.now:
        override = parse_timestamp(args.now)
        if not override:
            parser.error("Could not parse --now timestamp")
        now = override
    else:
        now = dt.datetime.now(dt.timezone.utc)

    selection: List[str] = []
    for name in args.chunks:
        if name not in CHUNKS:
            parser.error(f"Unknown chunk '{name}'. Valid options: {', '.join(sorted(CHUNKS))}")
        selection.append(name)

    for name in selection:
        stats = summarize_chunk(name, CHUNKS[name], manifest, now)
        remaining = len(stats["pending"])
        print(f"{name}: {stats['done']}/{stats['total']} complete, {remaining} remaining")
        if remaining:
            ids = " ".join(str(x) for x in stats["pending"])
            print(f"  Remaining IDs: {ids}")
        avg = stats["avg_interval_hours"]
        if avg:
            print(f"  Avg interval per subject: {avg:.2f} h")
        print(f"  Last finish: {format_time(stats['last_finish'])}")
        print(f"  ETA: {format_time(stats['eta'])}")
        print()

if __name__ == "__main__":
    main()
