#!/usr/bin/env python3
"""
Summarize coding-variant results at the transcript level.

Inputs:
  - one or more CDS FASTA files
  - one or more peptide FASTA files
  - one or more per-transcript CDS variant TSV files
  - one or more per-transcript variant_effects TSV files
  - one optional site summary TSV

Output:
  - one transcript summary TSV with columns such as:
      transcript
      cds_bp
      pep_aa
      pep_nostop_aa
      n_variants
      n_sites
      n_synonymous
      n_missense
      n_nonsense
      n_non_SNP
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


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

    seq = "".join(seq_parts)
    return header, seq


def transcript_from_geneid_filename(path: str):
    """
    Extract transcript from filenames like:
      GRMZM2G129913_T03.cds.fa
      GRMZM2G129913_T03.pep.fa
      GRMZM2G129913_T03.pep.nostop.fa
    Returns:
      T03
    """
    name = Path(path).name
    parts = name.split(".")[0]  # GRMZM2G129913_T03
    if "_" not in parts:
        return parts
    return parts.split("_")[-1]


def transcript_from_effects_filename(path: str):
    """
    Extract transcript from filenames like:
      SPO11.T03.CDS.variant_effects.tsv
      SPO11.T03.CDS.variants.tsv
    Returns:
      T03
    """
    name = Path(path).name
    parts = name.split(".")
    if len(parts) >= 2:
        return parts[1]
    return name


def count_data_rows(path: str):
    with open(path, "r") as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)


def parse_effect_counts(path: str):
    counts = defaultdict(int)
    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 13 or header[12] != "effect":
            die(f"Unexpected header in effects file: {path}")

        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 13:
                continue
            effect = cols[12]
            counts[effect] += 1
    return counts


def parse_site_summary_counts(path: str):
    """
    Count number of summarized sites per transcript from:
      summarize_cds_sites.py output
    """
    counts = defaultdict(int)

    if path is None:
        return counts

    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 1 or header[0] != "transcript":
            die(f"Unexpected header in site summary file: {path}")

        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            tx = cols[0]
            counts[tx] += 1

    return counts


def main():
    ap = argparse.ArgumentParser(description="Summarize results at the transcript level.")
    ap.add_argument("--cds", nargs="+", required=True, help="CDS FASTA files")
    ap.add_argument("--pep", nargs="+", required=True, help="Peptide FASTA files (with stop)")
    ap.add_argument("--pep-nostop", nargs="+", required=True, help="Peptide FASTA files (without stop)")
    ap.add_argument("--variants", nargs="+", required=True, help="Per-transcript CDS variants TSV files")
    ap.add_argument("--effects", nargs="+", required=True, help="Per-transcript variant_effects TSV files")
    ap.add_argument("--site-summary", default=None, help="Optional CDS site summary TSV")
    ap.add_argument("--out", required=True, help="Output transcript summary TSV")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    data = defaultdict(dict)

    # CDS lengths
    for path in args.cds:
        tx = transcript_from_geneid_filename(path)
        _, seq = read_fasta_one(path)
        data[tx]["cds_bp"] = len(seq)

    # peptide lengths
    for path in args.pep:
        tx = transcript_from_geneid_filename(path)
        _, seq = read_fasta_one(path)
        data[tx]["pep_aa"] = len(seq)

    for path in args.pep_nostop:
        tx = transcript_from_geneid_filename(path)
        _, seq = read_fasta_one(path)
        data[tx]["pep_nostop_aa"] = len(seq)

    # variant counts
    for path in args.variants:
        tx = transcript_from_effects_filename(path)
        data[tx]["n_variants"] = count_data_rows(path)

    # effect counts
    for path in args.effects:
        tx = transcript_from_effects_filename(path)
        counts = parse_effect_counts(path)
        data[tx]["n_synonymous"] = counts.get("synonymous", 0)
        data[tx]["n_missense"] = counts.get("missense", 0)
        data[tx]["n_nonsense"] = counts.get("nonsense", 0)
        data[tx]["n_non_SNP"] = counts.get("non_SNP", 0)
        data[tx]["n_outside_CDS"] = counts.get("outside_CDS", 0)

    # summarized site counts
    site_counts = parse_site_summary_counts(args.site_summary)
    for tx, n in site_counts.items():
        data[tx]["n_sites"] = n

    transcripts = sorted(data.keys())

    with open(args.out, "w") as out:
        out.write(
            "\t".join([
                "transcript",
                "cds_bp",
                "pep_aa",
                "pep_nostop_aa",
                "n_variants",
                "n_sites",
                "n_synonymous",
                "n_missense",
                "n_nonsense",
                "n_non_SNP",
                "n_outside_CDS",
            ]) + "\n"
        )

        for tx in transcripts:
            row = data[tx]
            out.write(
                "\t".join([
                    tx,
                    str(row.get("cds_bp", "NA")),
                    str(row.get("pep_aa", "NA")),
                    str(row.get("pep_nostop_aa", "NA")),
                    str(row.get("n_variants", 0)),
                    str(row.get("n_sites", 0)),
                    str(row.get("n_synonymous", 0)),
                    str(row.get("n_missense", 0)),
                    str(row.get("n_nonsense", 0)),
                    str(row.get("n_non_SNP", 0)),
                    str(row.get("n_outside_CDS", 0)),
                ]) + "\n"
            )

    print("== SUMMARIZE TRANSCRIPTS ==")
    print(f"[OK] transcripts={len(transcripts)}")
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
