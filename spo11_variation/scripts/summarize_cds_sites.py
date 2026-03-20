#!/usr/bin/env python3
"""
Summarize CDS variant effects by transcript site.

Input:
  One or more *.CDS.variant_effects.tsv files with columns:
    chrom  pos  id  ref  alt  cds_pos  codon_idx  pos_in_codon
    ref_codon  alt_codon  ref_aa  alt_aa  effect

Output:
  One TSV with one row per transcript + CDS site summarizing:
    transcript
    cds_pos
    codon_idx
    ref_codon
    ref_aa
    n_variants
    effects
    alt_aas
    alt_codons
    variant_ids
    genomic_positions
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def parse_effects_file(path: str):
    """
    Parse one variant_effects.tsv file.
    Returns rows as dicts.
    """
    rows = []

    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        expected = [
            "chrom", "pos", "id", "ref", "alt",
            "cds_pos", "codon_idx", "pos_in_codon",
            "ref_codon", "alt_codon", "ref_aa", "alt_aa", "effect"
        ]
        if header[:13] != expected:
            die(f"Unexpected header in {path}: {header}")

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 13:
                die(f"Invalid line in {path}:{line_num}")

            row = dict(zip(expected, cols))
            rows.append(row)

    return rows


def transcript_from_filename(path: str):
    """
    Extract transcript label from filename like:
      SPO11.T03.CDS.variant_effects.tsv -> T03
    """
    name = Path(path).name
    parts = name.split(".")
    if len(parts) < 4:
        return name
    # usually GENE.T03.CDS.variant_effects.tsv
    return parts[1]


def sort_unique(values, numeric=False):
    vals = sorted(set(values), key=(lambda x: int(x) if numeric else x))
    return ",".join(str(v) for v in vals)


def main():
    ap = argparse.ArgumentParser(description="Summarize CDS variant effects by site.")
    ap.add_argument(
        "--effects",
        nargs="+",
        required=True,
        help="One or more *.CDS.variant_effects.tsv files"
    )
    ap.add_argument("--out", required=True, help="Output summary TSV")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    grouped = defaultdict(lambda: {
        "ref_codon": None,
        "ref_aa": None,
        "effects": [],
        "alt_aas": [],
        "alt_codons": [],
        "variant_ids": [],
        "genomic_positions": [],
    })

    for fx in args.effects:
        transcript = transcript_from_filename(fx)
        rows = parse_effects_file(fx)

        for row in rows:
            cds_pos = row["cds_pos"]
            codon_idx = row["codon_idx"]
            ref_codon = row["ref_codon"]
            ref_aa = row["ref_aa"]

            # Skip rows that are not mapped cleanly to a codon/site
            if cds_pos == "NA" or codon_idx == "NA":
                continue

            key = (transcript, int(cds_pos), int(codon_idx), ref_codon, ref_aa)

            grouped[key]["ref_codon"] = ref_codon
            grouped[key]["ref_aa"] = ref_aa
            grouped[key]["effects"].append(row["effect"])
            grouped[key]["alt_aas"].append(row["alt_aa"])
            grouped[key]["alt_codons"].append(row["alt_codon"])
            grouped[key]["variant_ids"].append(row["id"])
            grouped[key]["genomic_positions"].append(row["pos"])

    with open(args.out, "w") as out:
        out.write(
            "\t".join([
                "transcript",
                "cds_pos",
                "codon_idx",
                "ref_codon",
                "ref_aa",
                "n_variants",
                "effects",
                "alt_aas",
                "alt_codons",
                "variant_ids",
                "genomic_positions",
            ]) + "\n"
        )

        for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1])):
            transcript, cds_pos, codon_idx, ref_codon, ref_aa = key
            g = grouped[key]

            out.write(
                "\t".join([
                    transcript,
                    str(cds_pos),
                    str(codon_idx),
                    ref_codon,
                    ref_aa,
                    str(len(g["variant_ids"])),
                    sort_unique(g["effects"]),
                    sort_unique(g["alt_aas"]),
                    sort_unique(g["alt_codons"]),
                    sort_unique(g["variant_ids"]),
                    sort_unique(g["genomic_positions"], numeric=True),
                ]) + "\n"
            )

    print("== SUMMARIZE CDS SITES ==")
    print(f"[OK] input_files={len(args.effects)}")
    print(f"[OK] summarized_sites={len(grouped)}")
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
