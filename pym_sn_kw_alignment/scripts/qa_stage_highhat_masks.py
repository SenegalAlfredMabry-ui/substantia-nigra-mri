#!/usr/bin/env python3
"""Copy Highhat v10/v9 mask files into the QA manual alignment directory."""

from __future__ import annotations

import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
HIGHHAT_LIST = REPO_ROOT / "tmp" / "highhat_with_kw_subjects.txt"
V10_BASE = REPO_ROOT / "MRI_data" / "TSE" / "probabilistic_masks_v10" / "base"
V9_BASE = REPO_ROOT / "MRI_data" / "TSE" / "probabilistic_masks_v9" / "base"
QA_ROOT = REPO_ROOT / "QA" / "manual_alignment_review"
KW_SOURCE_DIR = Path("${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/TSE/SN")


def load_subjects(list_path: Path) -> list[str]:
    if not list_path.exists():
        raise FileNotFoundError(f"Subject roster missing: {list_path}")
    return [
        line.strip()
        for line in list_path.read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def copy_v10_masks(sub: str) -> int:
    """Copy every v10 mask for a subject into QA, returns copied file count."""
    qa_dir = (QA_ROOT / sub)
    qa_dir.mkdir(parents=True, exist_ok=True)
    dest_dir = qa_dir / "black_ant_masks_v10"
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    for src in sorted(V10_BASE.glob(f"*_{sub}.nii*")):
        shutil.copy2(src, dest_dir / src.name)
        copied += 1
    if copied == 0:
        dest_dir.rmdir()
        return copied

    # Also drop a convenience copy of the SN mask at the QA root.
    sn_src = dest_dir / f"sn_groupmask_in_{sub}.nii.gz"
    if not sn_src.exists():
        sn_src = dest_dir / f"sn_groupmask_in_{sub}.nii"
    if sn_src.exists():
        suffixes = sn_src.suffixes
        ext = "".join(suffixes) if suffixes else sn_src.suffix
        sn_dest = qa_dir / f"sn_groupmask_in_{sub}_v10_hybrid{ext}"
        shutil.copy2(sn_src, sn_dest)
    return copied


def copy_v9_sn(sub: str) -> bool:
    """Copy the legacy v9 SN mask (if available) into QA."""
    qa_dir = (QA_ROOT / sub)
    qa_dir.mkdir(parents=True, exist_ok=True)
    for ext in (".nii.gz", ".nii"):
        src = V9_BASE / f"sn_groupmask_in_{sub}{ext}"
        if src.exists():
            dest = qa_dir / f"sn_groupmask_in_{sub}_v9{ext}"
            shutil.copy2(src, dest)
            return True
    return False


def copy_kw_tracings(sub: str) -> bool:
    """Copy the manual KW tracings from the MOTIP tree."""
    if not KW_SOURCE_DIR.exists():
        return False
    qa_dir = QA_ROOT / sub
    qa_dir.mkdir(parents=True, exist_ok=True)
    dest_dir = qa_dir / "kw_tracings"
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    candidates = [
        f"SN_ROI_KW_{sub}.nii",
        f"SN_ROI_KW_{sub}.nii.gz",
        f"SN_ROI_KW_{sub}.mgz",
        f"SN_ROI_KW_{sub}_ZP.nii",
        f"SN_ROI_KW_{sub}_ZP.nii.gz",
        f"SN_ROI_CNR_KW_{sub}.nii",
        f"SN_ROI_CNR_KW_{sub}.nii.gz",
    ]
    matches = []
    for name in candidates:
        src = KW_SOURCE_DIR / name
        if src.exists():
            matches.append(src)
    if not matches:
        dest_dir.rmdir()
        return False
    for src in sorted(matches):
        shutil.copy2(src, dest_dir / src.name)
    return True


def main() -> None:
    subjects = load_subjects(HIGHHAT_LIST)
    if not subjects:
        raise RuntimeError("Highhat subject list was empty")

    missing_v10: list[str] = []
    missing_v9: list[str] = []
    missing_kw: list[str] = []
    for sub in subjects:
        copied = copy_v10_masks(sub)
        if copied == 0:
            missing_v10.append(sub)
        if not copy_v9_sn(sub):
            missing_v9.append(sub)
        if not copy_kw_tracings(sub):
            missing_kw.append(sub)
    log_lines = []
    if missing_v10:
        log_lines.append(
            "[WARN] No v10 masks copied for: " + ", ".join(missing_v10)
        )
    if missing_v9:
        log_lines.append(
            "[WARN] v9 SN mask missing for: " + ", ".join(missing_v9)
        )
    if missing_kw:
        log_lines.append(
            "[WARN] KW tracings missing for: " + ", ".join(missing_kw)
        )
    if log_lines:
        qa_log = QA_ROOT / "highhat_mask_stage.log"
        qa_log.write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
