#!/usr/bin/env python3
"""
Plot stacked amino-acid frequency by protein position from an AA-frequency matrix.

Input matrix format:
  aa <tab> pos1 <tab> pos2 <tab> pos3 ...

Values are percentages (0–100).

This version fixes title/legend overlap by allocating dedicated figure space
for both elements.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY*")

PASTEL_AA_COLORS = {
    "A": "#d98ac7",
    "C": "#b7d4ea",
    "D": "#d9d9d9",
    "E": "#b7e0b1",
    "F": "#d0c1e8",
    "G": "#e8c47a",
    "H": "#c7c7e2",
    "I": "#d7d7d7",
    "K": "#b7dede",
    "L": "#cfcfcf",
    "M": "#9fd3e6",
    "N": "#e7d7d7",
    "P": "#4f8fc9",
    "Q": "#dedede",
    "R": "#e3c9c9",
    "S": "#b3b8d8",
    "T": "#9c7bc2",
    "V": "#c8c88a",
    "W": "#a36aac",
    "Y": "#9ed4df",
    "*": "#f2a6c9",
}


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def read_matrix(path: str):
    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 2 or header[0] != "aa":
            die(f"Unexpected matrix header: {header}")

        try:
            positions = [int(x) for x in header[1:]]
        except ValueError:
            die("Position columns must be integers")

        matrix = {}
        aa_rows = []

        for line_num, line in enumerate(f, start=2):
            cols = line.rstrip("\n").split("\t")
            if len(cols) != len(header):
                die(
                    f"Line {line_num} has {len(cols)} columns, expected {len(header)}"
                )

            aa = cols[0]
            try:
                vals = np.array([float(x) for x in cols[1:]], dtype=float)
            except ValueError:
                die(f"Non-numeric matrix value at line {line_num}")

            matrix[aa] = vals
            aa_rows.append(aa)

    if not aa_rows:
        die("No amino-acid rows found in matrix")

    return positions, aa_rows, matrix


def infer_ref_aa(matrix, aa_order, npos):
    ref = []
    for i in range(npos):
        best_aa = None
        best_val = -1.0
        for aa in aa_order:
            val = matrix[aa][i]
            if val > best_val:
                best_val = val
                best_aa = aa
        ref.append(best_aa)
    return ref


def main():
    ap = argparse.ArgumentParser(description="Plot stacked AA frequency by position.")
    ap.add_argument("--in-matrix", required=True, help="Input AA frequency matrix TSV")
    ap.add_argument("--out-png", required=True, help="Output PNG")
    ap.add_argument("--title", default="AA frequencies by position", help="Figure title")
    ap.add_argument("--width", type=float, default=18, help="Figure width in inches")
    ap.add_argument("--height", type=float, default=4.2, help="Figure height in inches")
    ap.add_argument(
        "--min-pct",
        type=float,
        default=0.05,
        help="Hide segments smaller than this percentage"
    )
    ap.add_argument(
        "--legend-ncol",
        type=int,
        default=8,
        help="Number of legend columns"
    )
    ap.add_argument(
        "--compress-reference",
        action="store_true",
        help="Hide dominant/reference AA fraction and only show non-reference fractions"
    )
    ap.add_argument(
        "--show-variable-rug",
        action="store_true",
        help="Draw small rug marks below variable positions"
    )
    ap.add_argument(
        "--bar-width",
        type=float,
        default=0.9,
        help="Bar width"
    )
    args = ap.parse_args()

    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)

    positions, aa_rows, matrix = read_matrix(args.in_matrix)
    aa_order = [aa for aa in DEFAULT_AA_ORDER if aa in matrix] + [
        aa for aa in aa_rows if aa not in DEFAULT_AA_ORDER
    ]

    npos = len(positions)
    x = np.array(positions, dtype=float)
    ref_aa = infer_ref_aa(matrix, aa_order, npos)

    # Build plotting values
    plot_vals = {}
    for aa in aa_order:
        vals = matrix[aa].copy()
        vals[vals < args.min_pct] = 0.0

        if args.compress_reference:
            for i in range(npos):
                if aa == ref_aa[i]:
                    vals[i] = 0.0

        plot_vals[aa] = vals

    # Use a 3-row layout: title / legend / plot
    fig = plt.figure(figsize=(args.width, args.height), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=3,
        ncols=1,
        height_ratios=[0.42, 0.55, 3.0],
        hspace=0.15
    )

    ax_title = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])
    ax = fig.add_subplot(gs[2])

    ax_title.axis("off")
    ax_legend.axis("off")

    # Title in its own row
    ax_title.text(
        0.5, 0.5, args.title,
        ha="center", va="center",
        fontsize=12
    )

    # Main plot
    bottom = np.zeros(npos, dtype=float)

    for aa in aa_order:
        vals = plot_vals[aa]
        if np.all(vals == 0):
            continue

        ax.bar(
            x,
            vals,
            bottom=bottom,
            width=args.bar_width,
            color=PASTEL_AA_COLORS.get(aa, "#cccccc"),
            label=aa,
            linewidth=0,
            align="center",
        )
        bottom += vals

    if args.compress_reference:
        ymax = max(5.0, float(np.nanmax(bottom)) * 1.08 if len(bottom) else 5.0)
        ax.set_ylim(0, ymax)
        ax.set_ylabel("Non-reference AA (%)")
    else:
        ax.set_ylim(0, 100)
        ax.set_ylabel("Percentage (%)")

    ax.set_xlabel("Amino acid position")
    ax.set_xlim(min(positions) - 1, max(positions) + 1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if args.show_variable_rug:
        y0, y1 = ax.get_ylim()
        h = (y1 - y0) * 0.018
        for i, pos in enumerate(positions):
            nonref = sum(matrix[aa][i] for aa in aa_order if aa != ref_aa[i])
            if nonref > 0:
                ax.vlines(pos, y0, y0 + h, color="black", linewidth=0.4)

    # Legend in its own row
    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    uniq_handles = []
    uniq_labels = []
    for h, l in zip(handles, labels):
        if l not in seen:
            uniq_handles.append(h)
            uniq_labels.append(l)
            seen.add(l)

    ax_legend.legend(
        uniq_handles,
        uniq_labels,
        title="AA",
        loc="center",
        ncol=args.legend_ncol,
        frameon=False,
        fontsize=8,
        title_fontsize=8,
        handlelength=1.2,
        columnspacing=0.9,
        handletextpad=0.4,
        borderaxespad=0.0,
    )

    fig.subplots_adjust(
        top=0.96,
        bottom=0.12,
        left=0.06,
        right=0.995
    )

    plt.savefig(args.out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print("== PLOT AA FREQUENCY STACKED ==")
    print(f"[OK] positions={len(positions)}")
    print(f"[OK] compress_reference={args.compress_reference}")
    print(f"[OK] wrote: {args.out_png}")


if __name__ == "__main__":
    main()
