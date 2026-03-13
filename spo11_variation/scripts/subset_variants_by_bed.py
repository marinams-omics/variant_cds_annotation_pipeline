#!/usr/bin/env python3
"""
Subset variants TSV by transcript CDS BED intervals.

Input BED:
  chrom  start0  end1  name  score  strand

Input variants TSV:
  expected at minimum:
    chrom  pos  ...
  (header optional; script autodetects and preserves header if present)

Output:
  variants TSV containing only rows whose positions fall within any BED interval

Coordinate logic:
  BED is 0-based, half-open
  variants pos is assumed 1-based
  inclusion rule:
      start0 < pos <= end1
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def ensure_parent(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_bed(path: str):
    """
    Read BED intervals into:
      intervals[chrom] = [(start0, end1), ...]
    """
    intervals = defaultdict(list)

    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 3:
                die(f"BED line {line_num} has fewer than 3 columns: {line}")

            chrom = cols[0]
            try:
                start0 = int(cols[1])
                end1 = int(cols[2])
            except ValueError:
                die(f"Invalid BED coordinates in line {line_num}: {line}")

            if start0 >= end1:
                die(f"Invalid BED interval in line {line_num}: {line}")

            intervals[chrom].append((start0, end1))

    if not intervals:
        die(f"No valid BED intervals found in: {path}")

    # sort intervals per chrom
    for chrom in intervals:
        intervals[chrom].sort()

    return intervals


def looks_like_header(line: str):
    """
    Heuristic: if col2 is not an integer, treat as header.
    """
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 2:
        return False
    try:
        int(cols[1])
        return False
    except ValueError:
        return True


def variant_in_intervals(chrom: str, pos1: int, intervals):
    """
    Check if a 1-based variant position falls inside any BED interval for that chrom.
    Inclusion rule:
      start0 < pos1 <= end1
    """
    if chrom not in intervals:
        return False

    for start0, end1 in intervals[chrom]:
        if start0 < pos1 <= end1:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description="Subset variants TSV by BED intervals.")
    ap.add_argument("--bed", required=True, help="Input BED file")
    ap.add_argument("--variants", required=True, help="Input variants TSV")
    ap.add_argument("--out", required=True, help="Output TSV with variants in BED intervals")
    args = ap.parse_args()

    print("== SUBSET VARIANTS BY BED ==")

    intervals = read_bed(args.bed)
    ensure_parent(args.out)

    total = 0
    kept = 0
    header_written = False

    with open(args.variants, "r") as fin, open(args.out, "w") as fout:
        for line_num, line in enumerate(fin, start=1):
            if not line.strip():
                continue

            # preserve header if present on first non-empty line
            if line_num == 1 and looks_like_header(line):
                fout.write(line)
                header_written = True
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                die(f"Variants line {line_num} has fewer than 2 columns: {line.strip()}")

            chrom = cols[0]
            try:
                pos1 = int(cols[1])
            except ValueError:
                die(f"Invalid variant position in line {line_num}: {line.strip()}")

            total += 1
            if variant_in_intervals(chrom, pos1, intervals):
                fout.write(line)
                kept += 1

    print(f"[OK] total_variants={total}")
    print(f"[OK] kept_variants={kept}")
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
