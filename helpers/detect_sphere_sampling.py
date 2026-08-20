#!/usr/bin/env python3
"""
detect_sphere_sampling.py

Detects whether diffusion MRI gradient directions (from a bvecs file)
were sampled on the whole sphere or half sphere.

Usage:
    python detect_sphere_sampling.py <bvecs_file> [--threshold 0.95] [--verbose]

Output:
    Prints "WHOLE_SPHERE" or "HALF_SPHERE" and exits with code 0 or 1 respectively.

For FSL eddy:
    WHOLE_SPHERE  ->  --slm=none  (default but can also use --slm=quadratic = recommended but slower)
    HALF_SPHERE   ->  --slm=linear
"""

import numpy as np
import argparse
import sys


def load_bvecs(bvecs_file):
    """Load bvecs file. Expected format: 3 x N text file."""
    bvecs = np.loadtxt(bvecs_file)
    # Ensure shape is (3, N)
    if bvecs.ndim == 1:
        bvecs = bvecs.reshape(-1, 1)
    if bvecs.shape[0] != 3 and bvecs.shape[1] == 3:
        bvecs = bvecs.T
    if bvecs.shape[0] != 3:
        raise ValueError(
            f"Unexpected bvecs shape: {bvecs.shape}. Expected 3 x N or N x 3."
        )
    return bvecs


def normalize_vectors(vectors):
    """Normalize each column vector to unit length."""
    norms = np.linalg.norm(vectors, axis=0)
    # Avoid division by zero for b=0 volumes
    norms[norms == 0] = 1.0
    return vectors / norms


def find_antipodal_pairs(bvecs, threshold=0.95):
    """
    For each gradient direction, check if an antipodal partner exists.

    Two vectors are antipodal if their dot product is close to -1.
    A threshold of 0.95 means the angle between them is > ~161.8°.

    Returns:
        n_with_partner: number of directions that have an antipodal partner
        n_total: total number of non-zero directions
        partner_fraction: fraction of directions with a partner
    """
    n = bvecs.shape[1]
    # Compute all pairwise dot products
    dot_products = bvecs.T @ bvecs  # N x N matrix

    # For each direction, check if any OTHER direction is antipodal
    n_with_partner = 0
    n_total = 0

    for i in range(n):
        norm_i = np.linalg.norm(bvecs[:, i])
        if norm_i < 1e-6:
            continue  # skip b=0 volumes
        n_total += 1

        # Check all other directions for antipodal match
        for j in range(n):
            if i == j:
                continue
            norm_j = np.linalg.norm(bvecs[:, j])
            if norm_j < 1e-6:
                continue
            if dot_products[i, j] <= -threshold:
                n_with_partner += 1
                break  # found a partner, move to next direction

    partner_fraction = n_with_partner / n_total if n_total > 0 else 0.0
    return n_with_partner, n_total, partner_fraction


def detect_sphere_sampling(bvecs_file, threshold=0.95, verbose=False):
    """
    Detect whether sampling was on the whole or half sphere.

    Parameters:
        bvecs_file: path to bvecs file
        threshold: dot product threshold for antipodal matching (default 0.95)
        verbose: print detailed information

    Returns:
        True if whole sphere, False if half sphere
    """
    bvecs = load_bvecs(bvecs_file)
    bvecs_norm = normalize_vectors(bvecs)

    n_with_partner, n_total, partner_fraction = find_antipodal_pairs(
        bvecs_norm, threshold
    )

    # If >50% of directions have an antipodal partner, it's whole sphere
    is_whole_sphere = partner_fraction > 0.5

    if verbose:
        print(f"Total non-zero directions: {n_total}")
        print(f"Directions with antipodal partner: {n_with_partner}")
        print(f"Fraction with partner: {partner_fraction:.2%}")
        print(f"Antipodal dot-product threshold: {threshold}")
        if is_whole_sphere:
            print("=> WHOLE_SPHERE sampling detected")
            print("=> FSL eddy recommendation: --slm=none (default but can also use --slm=quadratic = recommended but slower)")
        else:
            print("=> HALF_SPHERE sampling detected")
            print("=> FSL eddy recommendation: --slm=linear")

    return is_whole_sphere


def main():
    parser = argparse.ArgumentParser(
        description="Detect whole-sphere vs half-sphere sampling from bvecs."
    )
    parser.add_argument("bvecs", help="Path to bvecs file")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.95,
        help="Dot product threshold for antipodal matching (default: 0.95, "
        "meaning angle > ~161.8°). Lower = more tolerant.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Print detailed information"
    )
    args = parser.parse_args()

    is_whole_sphere = detect_sphere_sampling(
        args.bvecs, args.threshold, args.verbose
    )

    if is_whole_sphere:
        print("WHOLE_SPHERE")
        sys.exit(0)
    else:
        print("HALF_SPHERE")
        sys.exit(1)


if __name__ == "__main__":
    main()
