#!/usr/bin/env python3
"""
Extract variants from a VCF that fall within the genomic region defined by a gene GFF.

Inputs:
  - gene GFF/GFF3 defining the gene region
  - VCF/VCF.gz

Output:
  - TSV with columns:
      chrom  pos  id  ref  alt

Coordinate logic:
  GFF coordinates are 1-based inclusive
  VCF positions are 1-based
  Inclusion rule:
      start <= pos <= end

Chromosome naming is normalized:
  chr5 / Chr5 / 5 -> 5
"""

import argparse
import gzip
import sys
from pathlib import Path


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def ensure_parent(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def normalize_chromosome(chromosome: str) -> str:
    if chromosome.startswith("Chr"):
        stripped = chromosome[3:]
        return stripped if stripped else chromosome
    if chromosome.startswith("chr"):
        stripped = chromosome[3:]
        return stripped if stripped else chromosome
    return chromosome


def parse_gene_region(gff_path: str):
    """
    Parse gene GFF and return:
      chrom, start1, end1

    If multiple non-comment lines exist, they must be on the same chromosome.
    The full span is used.
    """
    chrom = None
    start1 = None
    end1 = None

    with open(gff_path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 5:
                die(f"GFF line {line_num} has fewer than 5 columns: {line}")

            this_chrom = cols[0]
            try:
                this_start = int(cols[3])
                this_end = int(cols[4])
            except ValueError:
                die(f"Invalid GFF coordinates in line {line_num}: {line}")

            if this_start > this_end:
                die(f"Invalid GFF interval in line {line_num}: {line}")

            if chrom is None:
                chrom = this_chrom
                start1 = this_start
                end1 = this_end
            else:
                if this_chrom != chrom:
                    die(f"Gene GFF contains multiple chromosomes: {chrom} and {this_chrom}")
                start1 = min(start1, this_start)
                end1 = max(end1, this_end)

    if chrom is None:
        die(f"No valid feature lines found in gene GFF: {gff_path}")

    return chrom, start1, end1


def main():
    ap = argparse.ArgumentParser(description="Extract region variants directly from VCF.")
    ap.add_argument("--gene-gff", required=True, help="Gene GFF/GFF3 defining the region")
    ap.add_argument("--vcf", required=True, help="Input VCF/VCF.gz")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    if not Path(args.gene_gff).exists():
        die(f"Gene GFF not found: {args.gene_gff}")
    if not Path(args.vcf).exists():
        die(f"VCF not found: {args.vcf}")

    print("== EXTRACT REGION VARIANTS FROM VCF ==")

    chrom, start1, end1 = parse_gene_region(args.gene_gff)
    chrom_norm = normalize_chromosome(chrom)

    print(f"[INFO] gene_region={chrom}:{start1}-{end1}")
    print(f"[INFO] normalized_gene_chrom={chrom_norm}")

    ensure_parent(args.out)

    total_records = 0
    kept_records = 0

    with open_maybe_gzip(args.vcf) as fin, open(args.out, "w") as fout:
        fout.write("chrom\tpos\tid\tref\talt\n")

        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue

            v_chrom = cols[0]
            try:
                pos1 = int(cols[1])
            except ValueError:
                continue

            vid = cols[2]
            ref = cols[3]
            alt = cols[4]

            total_records += 1

            if normalize_chromosome(v_chrom) == chrom_norm and start1 <= pos1 <= end1:
                fout.write(f"{normalize_chromosome(v_chrom)}\t{pos1}\t{vid}\t{ref}\t{alt}\n")
                kept_records += 1

    print(f"[OK] total_vcf_records_scanned={total_records}")
    print(f"[OK] kept_in_region={kept_records}")
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
