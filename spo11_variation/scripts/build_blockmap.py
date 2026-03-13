#!/usr/bin/env python3
"""
Build CDS blockmap from a transcript CDS BED file.

Input BED format:
  chrom  start0  end1  name  score  strand

Output TSV columns:
  block  chrom  start0  end1  strand  len  cds_start  cds_end

Notes:
- BED coordinates are assumed 0-based, half-open
- CDS coordinates in output are 1-based, inclusive
- The script preserves BED row order exactly
"""

import argparse
import sys
from pathlib import Path


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def ensure_dir(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_bed(path: str):
    rows = []
    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 6:
                die(f"BED line {line_num} has fewer than 6 columns: {line}")

            chrom, start0, end1, name, score, strand = cols[:6]

            try:
                start0 = int(start0)
                end1 = int(end1)
            except ValueError:
                die(f"Invalid integer coordinates in BED line {line_num}: {line}")

            if start0 >= end1:
                die(f"Invalid BED interval (start0 >= end1) in line {line_num}: {line}")

            if strand not in {"+", "-"}:
                die(f"Invalid strand in BED line {line_num}: {line}")

            rows.append((chrom, start0, end1, name, strand))

    if not rows:
        die(f"No valid BED rows found in: {path}")

    return rows


def build_blockmap(rows):
    """
    Build blockmap preserving input row order.
    """
    out = []
    cds_cursor = 1  # 1-based CDS coordinate

    for chrom, start0, end1, name, strand in rows:
        block_len = end1 - start0
        cds_start = cds_cursor
        cds_end = cds_cursor + block_len - 1

        out.append({
            "block": name,
            "chrom": chrom,
            "start0": start0,
            "end1": end1,
            "strand": strand,
            "len": block_len,
            "cds_start": cds_start,
            "cds_end": cds_end,
        })

        cds_cursor = cds_end + 1

    return out


def write_blockmap(blockmap, out_path: str):
    ensure_dir(out_path)
    with open(out_path, "w") as out:
        out.write("block\tchrom\tstart0\tend1\tstrand\tlen\tcds_start\tcds_end\n")
        for row in blockmap:
            out.write(
                f"{row['block']}\t{row['chrom']}\t{row['start0']}\t{row['end1']}\t"
                f"{row['strand']}\t{row['len']}\t{row['cds_start']}\t{row['cds_end']}\n"
            )


def main():
    ap = argparse.ArgumentParser(description="Build CDS blockmap from BED.")
    ap.add_argument("--bed", required=True, help="Input CDS BED file")
    ap.add_argument("--out", required=True, help="Output blockmap TSV")
    args = ap.parse_args()

    print("== BUILD BLOCKMAP ==")
    rows = read_bed(args.bed)

    # Basic consistency checks
    chroms = {r[0] for r in rows}
    strands = {r[4] for r in rows}

    if len(chroms) != 1:
        die(f"BED contains multiple chromosomes: {sorted(chroms)}")
    if len(strands) != 1:
        die(f"BED contains multiple strands: {sorted(strands)}")

    blockmap = build_blockmap(rows)
    write_blockmap(blockmap, args.out)

    total_bp = sum(r["len"] for r in blockmap)
    print(f"[OK] blocks={len(blockmap)} chrom={list(chroms)[0]} strand={list(strands)[0]} total_CDS_bp={total_bp}")
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
