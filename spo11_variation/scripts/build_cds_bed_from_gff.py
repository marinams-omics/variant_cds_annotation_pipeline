#!/usr/bin/env python3
"""
Build CDS BED files from project transcript definitions and full annotation.

Inputs:
  - project_gff: small GFF with target transcript entries
  - annotation_gff: full genome annotation containing CDS features

Outputs:
  - one CDS BED per transcript for FASTA-compatible chromosome naming
  - one CDS BED per transcript for VCF-compatible chromosome naming

BED format:
  chrom  start0  end1  name  score  strand

Notes:
  - BED coordinates are 0-based, half-open
  - Transcript IDs are derived from project_gff
  - CDS features are recovered from annotation_gff
"""

import argparse
import gzip
import sys
from pathlib import Path
from collections import defaultdict


# -------------------------
# basic helpers
# -------------------------

def die(msg: str):
    sys.exit(f"[ERROR] {msg}")


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)


# -------------------------
# chromosome naming
# -------------------------

def to_fasta_chrom(chrom: str) -> str:
    """
    Convert chrom name to FASTA-compatible style: Chr5
    Examples:
      chr5 -> Chr5
      5    -> Chr5
      Chr5 -> Chr5
    Leaves scaffold-like names unchanged unless they are chr/Chr-prefixed.
    """
    if chrom.startswith("Chr"):
        return chrom
    if chrom.startswith("chr"):
        return "Chr" + chrom[3:]
    if chrom.isdigit():
        return "Chr" + chrom
    return chrom


def to_vcf_chrom(chrom: str) -> str:
    """
    Convert chrom name to VCF-compatible style: 5
    Examples:
      chr5 -> 5
      Chr5 -> 5
      5    -> 5
    Leaves scaffold-like names unchanged unless they are chr/Chr-prefixed.
    """
    if chrom.startswith("Chr"):
        return chrom[3:]
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom


# -------------------------
# parsing GFFs
# -------------------------

def parse_attributes(attr_string: str):
    d = {}
    for field in attr_string.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
    return d


def parse_project_transcripts(project_gff: str):
    """
    Returns sorted target transcript IDs from the project GFF.
    """
    transcripts = set()

    with open(project_gff, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = cols
            if feature != "transcript":
                continue

            attr_dict = parse_attributes(attrs)
            tx_id = attr_dict.get("ID")
            if tx_id:
                transcripts.add(tx_id)

    if not transcripts:
        die(f"No transcript IDs found in project GFF: {project_gff}")

    return sorted(transcripts)


def parse_cds_from_annotation(annotation_gff: str, target_transcripts):
    """
    Extract CDS features for target transcripts from full annotation.

    Returns:
      cds_by_tx: dict
        tx_id -> list of tuples:
          (chrom, start1, end1, strand, phase)
    """
    target_set = set(target_transcripts)
    cds_by_tx = defaultdict(list)

    with open_maybe_gzip(annotation_gff) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = cols
            if feature != "CDS":
                continue

            attr_dict = parse_attributes(attrs)
            parent = attr_dict.get("Parent")
            if not parent:
                continue

            if parent not in target_set:
                continue

            try:
                start1 = int(start)
                end1 = int(end)
            except ValueError:
                die(f"Invalid CDS coordinates in annotation GFF: {line}")

            cds_by_tx[parent].append((chrom, start1, end1, strand, phase))

    return cds_by_tx


# -------------------------
# BED writing
# -------------------------

def sort_cds_features(features):
    """
    Sort CDS features in genomic order (ascending start).
    BED writing itself does not need transcript-order reversal for minus strand.
    For sequence concatenation later, strand-aware tools/use of bedtools -s matters.
    """
    return sorted(features, key=lambda x: (x[1], x[2]))


def write_bed(features, tx_id: str, out_path: str, chrom_mode: str):
    """
    Write one BED file for one transcript.

    chrom_mode:
      - "fa"  -> FASTA-compatible chrom naming
      - "vcf" -> VCF-compatible chrom naming
    """
    with open(out_path, "w") as out:
        for i, (chrom, start1, end1, strand, phase) in enumerate(features, start=1):
            if chrom_mode == "fa":
                bed_chrom = to_fasta_chrom(chrom)
            elif chrom_mode == "vcf":
                bed_chrom = to_vcf_chrom(chrom)
            else:
                die(f"Unknown chrom_mode: {chrom_mode}")

            start0 = start1 - 1
            name = f"{tx_id}|CDS{i}"
            score = 0
            out.write(f"{bed_chrom}\t{start0}\t{end1}\t{name}\t{score}\t{strand}\n")


# -------------------------
# main
# -------------------------

def main():
    ap = argparse.ArgumentParser(description="Build CDS BED files per transcript from GFF annotation.")
    ap.add_argument("--project-gff", required=True, help="Project transcript GFF (target transcripts only)")
    ap.add_argument("--annotation-gff", required=True, help="Full annotation GFF/GFF3 (may be .gz)")
    ap.add_argument("--out-bed-fa-dir", required=True, help="Output directory for FASTA-compatible BEDs")
    ap.add_argument("--out-bed-vcf-dir", required=True, help="Output directory for VCF-compatible BEDs")
    args = ap.parse_args()

    ensure_dir(args.out_bed_fa_dir)
    ensure_dir(args.out_bed_vcf_dir)

    print("== BUILD CDS BED FROM GFF ==")

    target_transcripts = parse_project_transcripts(args.project_gff)
    print(f"[INFO] Found {len(target_transcripts)} target transcript(s) in project GFF")

    cds_by_tx = parse_cds_from_annotation(args.annotation_gff, target_transcripts)

    missing = [tx for tx in target_transcripts if tx not in cds_by_tx]
    if missing:
        die("The following transcripts had no CDS in full annotation: " + ", ".join(missing))

    for tx in target_transcripts:
        features = sort_cds_features(cds_by_tx[tx])

        # Basic consistency checks
        strands = {x[3] for x in features}
        chroms = {x[0] for x in features}

        if len(strands) != 1:
            die(f"Transcript {tx} has CDS features on multiple strands: {sorted(strands)}")
        if len(chroms) != 1:
            die(f"Transcript {tx} has CDS features on multiple chromosomes: {sorted(chroms)}")

        out_fa = Path(args.out_bed_fa_dir) / f"{tx}.cds.bed"
        out_vcf = Path(args.out_bed_vcf_dir) / f"{tx}.cds.bed"

        write_bed(features, tx, str(out_fa), chrom_mode="fa")
        write_bed(features, tx, str(out_vcf), chrom_mode="vcf")

        total_bp = sum(end1 - start1 + 1 for (_, start1, end1, _, _) in features)
        strand = list(strands)[0]
        chrom = list(chroms)[0]

        print(
            f"[OK] {tx}: {len(features)} CDS features | "
            f"chrom={chrom} | strand={strand} | total_CDS_bp={total_bp}"
        )
        print(f"     -> {out_fa}")
        print(f"     -> {out_vcf}")

    print("✔ BED generation completed")


if __name__ == "__main__":
    main()
