#!/usr/bin/env python3
import os
import shutil
import argparse
import json
from datetime import datetime

import numpy as np
import nibabel as nib


def backup_file(path):
    """Rename original file to file.orig (fails if .orig already exists)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    backup_path = path + ".orig"
    # The global pre-check should prevent this, but keep it here as safety:
    if os.path.exists(backup_path):
        raise FileExistsError(
            f"Refusing to overwrite existing backup file: {backup_path}"
        )
    shutil.move(path, backup_path)
    return backup_path


def load_bvals(path):
    with open(path, 'r') as f:
        txt = f.read().strip()
    parts = txt.split()
    return np.array([float(p) for p in parts])


def save_bvals(path, bvals):
    with open(path, 'w') as f:
        f.write(" ".join(
            str(int(b)) if float(b).is_integer() else str(b)
            for b in bvals
        ))


def load_bvecs(path):
    # Expect 3 rows, N columns
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr[np.newaxis, :]
    return arr


def save_bvecs(path, bvecs):
    np.savetxt(path, bvecs, fmt="%.10f")


def parse_drop_arg(drop_str):
    """
    Parse a string like '0,2,5-8' into a sorted list of unique integers.
    """
    to_drop = set()
    for part in drop_str.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-")
            start, end = int(start), int(end)
            to_drop.update(range(start, end + 1))
        else:
            to_drop.add(int(part))
    return sorted(to_drop)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Delete DWI volumes and corresponding bvals/bvecs, "
            "write a JSON log of changes, and refuse to run if .orig backups exist."
        )
    )
    parser.add_argument("nii", help="Input NIfTI file (.nii or .nii.gz)")
    parser.add_argument("bval", help="Input bval file")
    parser.add_argument("bvec", help="Input bvec file")
    parser.add_argument(
        "--drop",
        required=True,
        help=(
            "Comma-separated list of zero-based volume indices to delete, "
            "e.g. '0,5,10' or '3-7,12'."
        )
    )
    args = parser.parse_args()

    nii_path = args.nii
    bval_path = args.bval
    bvec_path = args.bvec

    # -------- FAILSAFE: refuse to run if any .orig already exists --------
    conflict_files = []
    for p in (nii_path, bval_path, bvec_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"File not found: {p}")
        backup_path = p + ".orig"
        if os.path.exists(backup_path):
            conflict_files.append(backup_path)

    if conflict_files:
        conflicts = "\n  ".join(conflict_files)
        raise SystemExit(
            "Failsafe triggered: .orig backup file(s) already exist.\n"
            "The script will NOT run to avoid overwriting previous backups.\n"
            "Conflicting backup files:\n"
            f"  {conflicts}\n\n"
            "If you really want to re-run, move or delete these .orig files first."
        )

    # ---------------- Parse indices to drop ----------------
    to_drop = parse_drop_arg(args.drop)
    print(f"Volumes to drop (0-based): {to_drop}")

    # ---------------- Load NIfTI ----------------
    img = nib.load(nii_path)
    data = img.get_fdata()
    affine = img.affine
    header = img.header

    if data.ndim != 4:
        raise ValueError(f"Expected 4D NIfTI, got shape {data.shape}")

    orig_shape = data.shape
    n_vols = orig_shape[3]
    print(f"Found {n_vols} volumes in NIfTI.")

    # Validate indices
    for idx in to_drop:
        if idx < 0 or idx >= n_vols:
            raise IndexError(f"Volume index {idx} out of range (0..{n_vols-1})")

    # ---------------- Load bvals/bvecs ----------------
    bvals = load_bvals(bval_path)
    bvecs = load_bvecs(bvec_path)

    if bvals.shape[0] != n_vols:
        raise ValueError(
            f"bvals length ({bvals.shape[0]}) != NIfTI volumes ({n_vols})"
        )
    if bvecs.shape[-1] != n_vols:
        raise ValueError(
            f"bvecs columns ({bvecs.shape[-1]}) != NIfTI volumes ({n_vols})"
        )

    # ---------------- Compute indices to keep ----------------
    keep_indices = [i for i in range(n_vols) if i not in to_drop]
    print(f"Keeping {len(keep_indices)} volumes.")

    # ---------------- Subset data ----------------
    new_data = data[..., keep_indices]
    new_bvals = bvals[keep_indices]
    new_bvecs = bvecs[..., keep_indices]
    new_shape = new_data.shape

    # ---------------- Backup originals (.orig) ----------------
    nii_backup = backup_file(nii_path)
    bval_backup = backup_file(bval_path)
    bvec_backup = backup_file(bvec_path)

    print(f"Backed up NIfTI to: {nii_backup}")
    print(f"Backed up bvals to: {bval_backup}")
    print(f"Backed up bvecs to: {bvec_backup}")

    # ---------------- Save new NIfTI with original name ----------------
    new_img = nib.Nifti1Image(new_data, affine, header)
    nib.save(new_img, nii_path)
    print(f"Saved updated NIfTI: {nii_path}")

    # ---------------- Save new bvals/bvecs with original names ----------------
    save_bvals(bval_path, new_bvals)
    print(f"Saved updated bvals: {bval_path}")

    save_bvecs(bvec_path, new_bvecs)
    print(f"Saved updated bvecs: {bvec_path}")

    # ---------------- Create JSON log ----------------
    base = os.path.splitext(os.path.basename(nii_path))[0]
    # handle .nii.gz: splitext only removes .gz
    if base.endswith(".nii"):
        base = base[:-4]
    json_name = f"{base}_desc-cleaned.json"
    json_path = os.path.join(os.path.dirname(nii_path), json_name)

    log = {
        "timestamp": datetime.now().isoformat(),
        "operation": "delete_dwi_volumes",
        "input_files_original": {
            "nii": nii_backup,
            "bval": bval_backup,
            "bvec": bvec_backup,
        },
        "output_files_new": {
            "nii": nii_path,
            "bval": bval_path,
            "bvec": bvec_path,
            "log_json": json_path,
        },
        "volumes_removed_indices_0_based": to_drop,
        "volumes_kept_indices_0_based": keep_indices,
        "original_n_volumes": int(n_vols),
        "new_n_volumes": int(len(keep_indices)),
        "original_shape": list(orig_shape),
        "new_shape": list(new_shape),
    }

    with open(json_path, "w") as f:
        json.dump(log, f, indent=4)

    print(f"Wrote JSON log: {json_path}")
    print("Done.")


if __name__ == "__main__":
    main()
