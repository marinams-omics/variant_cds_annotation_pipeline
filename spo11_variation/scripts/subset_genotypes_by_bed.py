#!/usr/bin/env python3
"""
Subset a genotype-preserving VCF by BED intervals.

Inputs:
  - BED file:
      chrom  start0  end1  ...
  - VCF/VCF.gz

Output:
  - filtered VCF containing only records inside BED intervals

Coordinate logic:
  BED is 0-based, half-open
  VCF positions are 1-based
  Inclusion rule:
      start0 < pos1 <= end1
"""

import argparse
import gzip
import sys
from pathlib import Path
from collections import defaultdict


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


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

    for chrom in intervals:
        intervals[chrom].sort()

    return intervals


def in_intervals(chrom: str, pos1: int, intervals):
    """
    Check if VCF 1-based position falls in any BED interval for the same chromosome.
    Rule:
      start0 < pos1 <= end1
    """
    if chrom not in intervals:
        return False

    for start0, end1 in intervals[chrom]:
        if start0 < pos1 <= end1:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description="Subset a genotype VCF by BED intervals.")
    ap.add_argument("--bed", required=True, help="Input BED file")
    ap.add_argument("--vcf", required=True, help="Input VCF/VCF.gz")
    ap.add_argument("--out", required=True, help="Output filtered VCF")
    args = ap.parse_args()

    if not Path(args.bed).exists():
        die(f"BED not found: {args.bed}")
    if not Path(args.vcf).exists():
        die(f"VCF not found: {args.vcf}")

    print("== SUBSET GENOTYPES BY BED ==")

    intervals = read_bed(args.bed)
    ensure_parent(args.out)

    header_lines = 0
    total_records = 0
    kept_records = 0

    with open_maybe_gzip(args.vcf) as fin, open(args.out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                header_lines += 1
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom = cols[0]
            try:
                pos1 = int(cols[1])
            except ValueError:
                continue

            total_records += 1

            if in_intervals(chrom, pos1, intervals):
                fout.write(line)
                kept_records += 1

    print(f"[OK] header_lines_written={header_lines}")
    print(f"[OK] total_vcf_records={total_records}")
    print(f"[OK] kept_records={kept_records}")
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
