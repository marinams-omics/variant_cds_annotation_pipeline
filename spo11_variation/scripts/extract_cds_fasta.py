#!/usr/bin/env python3
"""
Extract and concatenate CDS sequence from genome FASTA using a transcript CDS BED.

Inputs:
  - genome FASTA (.fa or .fa.gz)
  - transcript CDS BED (FASTA-compatible chromosome names)

Outputs:
  - blocks TSV:
      block_name   sequence
  - CDS FASTA:
      >transcript_id_cds
      ACTG...

Assumptions:
  - BED has 6 columns:
      chrom  start0  end1  name  score  strand
  - BED is transcript-specific
  - Coordinates are 0-based, half-open
"""

import argparse
import gzip
import sys
from pathlib import Path


# -------------------------
# helpers
# -------------------------

def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def ensure_parent(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


# -------------------------
# FASTA loading
# -------------------------

def read_genome_fasta(path: str):
    """
    Load genome FASTA into memory as:
      chrom -> sequence string
    """
    genome = {}
    chrom = None
    chunks = []

    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    genome[chrom] = "".join(chunks)
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if chrom is not None:
        genome[chrom] = "".join(chunks)

    if not genome:
        die(f"No sequences found in FASTA: {path}")

    return genome


# -------------------------
# BED parsing
# -------------------------

def read_bed(path: str):
    """
    Returns list of tuples:
      (chrom, start0, end1, name, strand)
    """
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
                die(f"Invalid BED coordinates in line {line_num}: {line}")

            if start0 >= end1:
                die(f"Invalid interval in BED line {line_num}: {line}")

            if strand not in {"+", "-"}:
                die(f"Invalid strand in BED line {line_num}: {line}")

            rows.append((chrom, start0, end1, name, strand))

    if not rows:
        die(f"No valid BED rows found in: {path}")

    return rows


def sort_in_transcript_order(rows):
    """
    Sort blocks in transcript order:
      + strand -> ascending genomic start
      - strand -> descending genomic start
    """
    strands = {r[4] for r in rows}
    if len(strands) != 1:
        die(f"BED contains multiple strands: {sorted(strands)}")
    strand = next(iter(strands))

    chroms = {r[0] for r in rows}
    if len(chroms) != 1:
        die(f"BED contains multiple chromosomes: {sorted(chroms)}")

    if strand == "+":
        return sorted(rows, key=lambda x: (x[1], x[2]))
    else:
        return sorted(rows, key=lambda x: (x[1], x[2]), reverse=True)


# -------------------------
# sequence extraction
# -------------------------

def extract_blocks(rows, genome):
    """
    Extract sequences for each BED row.
    Returns list of tuples:
      (block_name, sequence)
    """
    out = []

    for chrom, start0, end1, name, strand in rows:
        if chrom not in genome:
            die(f"Chromosome {chrom} from BED not found in FASTA")

        chrom_seq = genome[chrom]
        if end1 > len(chrom_seq):
            die(f"BED interval {chrom}:{start0}-{end1} exceeds FASTA chromosome length")

        seq = chrom_seq[start0:end1]
        if strand == "-":
            seq = revcomp(seq)

        out.append((name, seq))

    return out


def transcript_id_from_bed(rows):
    """
    Derive transcript ID from BED block names like:
      GRMZM2G129913_T03|CDS1
    """
    name = rows[0][3]
    if "|CDS" in name:
        return name.split("|CDS")[0]
    return name


# -------------------------
# output
# -------------------------

def write_blocks_tsv(blocks, out_path: str):
    ensure_parent(out_path)
    with open(out_path, "w") as out:
        out.write("block\tsequence\n")
        for name, seq in blocks:
            out.write(f"{name}\t{seq}\n")


def write_cds_fasta(tx_id: str, cds_seq: str, out_path: str):
    ensure_parent(out_path)
    with open(out_path, "w") as out:
        out.write(f">{tx_id}_cds\n")
        out.write(wrap_fasta(cds_seq) + "\n")


# -------------------------
# main
# -------------------------

def main():
    ap = argparse.ArgumentParser(description="Extract and concatenate CDS sequence from genome FASTA using BED.")
    ap.add_argument("--genome", required=True, help="Genome FASTA (.fa or .fa.gz)")
    ap.add_argument("--bed", required=True, help="Transcript CDS BED (FASTA-compatible chrom names)")
    ap.add_argument("--out-blocks", required=True, help="Output blocks TSV")
    ap.add_argument("--out-cds", required=True, help="Output concatenated CDS FASTA")
    args = ap.parse_args()

    print("== EXTRACT CDS FASTA ==")

    genome = read_genome_fasta(args.genome)
    rows = read_bed(args.bed)
    rows = sort_in_transcript_order(rows)

    tx_id = transcript_id_from_bed(rows)
    blocks = extract_blocks(rows, genome)
    cds_seq = "".join(seq for _, seq in blocks)

    write_blocks_tsv(blocks, args.out_blocks)
    write_cds_fasta(tx_id, cds_seq, args.out_cds)

    print(f"[OK] transcript={tx_id}")
    print(f"[OK] blocks={len(blocks)}")
    print(f"[OK] cds_bp={len(cds_seq)}")
    print(f"[OK] wrote blocks: {args.out_blocks}")
    print(f"[OK] wrote cds:    {args.out_cds}")


if __name__ == "__main__":
    main()
