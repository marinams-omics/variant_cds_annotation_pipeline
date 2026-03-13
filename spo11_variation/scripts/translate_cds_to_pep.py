#!/usr/bin/env python3
"""
Translate CDS FASTA to peptide FASTA.

Inputs:
  - CDS FASTA (single sequence)

Outputs:
  - peptide FASTA with terminal stop if present
  - peptide FASTA without terminal stop

Checks:
  - CDS length must be multiple of 3
  - reports number of internal stops
"""

import argparse
import sys
from pathlib import Path


CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def ensure_parent(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_fasta_one(path: str):
    header = None
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:].split()[0]
                continue
            seq_parts.append(line)

    if header is None:
        die(f"No FASTA header found in: {path}")

    seq = "".join(seq_parts).upper()
    if not seq:
        die(f"No sequence found in: {path}")

    return header, seq


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def translate_cds(seq: str) -> str:
    if len(seq) % 3 != 0:
        die(f"CDS length is not divisible by 3: {len(seq)}")

    pep = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "X")
        pep.append(aa)
    return "".join(pep)


def write_fasta(header: str, seq: str, path: str):
    ensure_parent(path)
    with open(path, "w") as out:
        out.write(f">{header}\n")
        out.write(wrap_fasta(seq) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Translate CDS FASTA into peptide FASTA.")
    ap.add_argument("--cds", required=True, help="Input CDS FASTA")
    ap.add_argument("--out-pep", required=True, help="Output peptide FASTA (with stop if present)")
    ap.add_argument("--out-pep-nostop", required=True, help="Output peptide FASTA without terminal stop")
    args = ap.parse_args()

    print("== TRANSLATE CDS TO PEP ==")

    header, cds_seq = read_fasta_one(args.cds)
    pep_seq = translate_cds(cds_seq)

    internal_stops = pep_seq[:-1].count("*") if pep_seq.endswith("*") else pep_seq.count("*")
    ends_with_stop = pep_seq.endswith("*")

    pep_nostop = pep_seq[:-1] if ends_with_stop else pep_seq

    pep_header = header.replace("_cds", "_pep")

    write_fasta(pep_header, pep_seq, args.out_pep)
    write_fasta(pep_header, pep_nostop, args.out_pep_nostop)

    print(f"[OK] cds_bp={len(cds_seq)}")
    print(f"[OK] pep_aa={len(pep_seq)}")
    print(f"[OK] pep_nostop_aa={len(pep_nostop)}")
    print(f"[OK] ends_with_stop={ends_with_stop}")
    print(f"[OK] internal_stops={internal_stops}")
    print(f"[OK] wrote pep:        {args.out_pep}")
    print(f"[OK] wrote pep_nostop: {args.out_pep_nostop}")


if __name__ == "__main__":
    main()
