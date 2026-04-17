#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import tifffile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert STARmap goodPoints.mat and optional labels.npz into Baysor CSV/TIFF inputs."
    )
    parser.add_argument("good_points_mat", type=Path, help="Path to goodPoints.mat")
    parser.add_argument("cell_barcode_names_csv", type=Path, help="Path to cell_barcode_names.csv")
    parser.add_argument(
        "labels_npz",
        type=Path,
        nargs="?",
        help="Optional labels.npz file. If provided, writes segmentation.tiff.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for molecules.csv and segmentation.tiff (default: current directory).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    mat = scipy.io.loadmat(args.good_points_mat)
    good_points = np.asarray(mat["goodPoints"])
    good_bases = np.asarray(mat["goodBases"]).ravel()

    gene_table = pd.read_csv(args.cell_barcode_names_csv, header=None)
    gene_per_base = dict(zip(gene_table.iloc[:, 1].astype(str), gene_table.iloc[:, 2].astype(str)))

    molecules = pd.DataFrame(good_points, columns=["x", "y", "z"])
    molecules["gene"] = [gene_per_base[str(base)] for base in good_bases]
    molecules.to_csv(args.output_dir / "molecules.csv", index=False)

    if args.labels_npz is not None:
        labels = np.load(args.labels_npz)["labels"]
        tifffile.imwrite(args.output_dir / "segmentation.tiff", labels.astype(np.uint16))


if __name__ == "__main__":
    main()
