#!/usr/bin/env python3
"""
CDS Variant Pipeline Input Validator
====================================

Validates the key inputs for the CDS variant annotation pipeline.

This script distinguishes between:
1. project_gff: a small project-specific GFF containing the target transcripts
2. annotation_gff: the full genome annotation GFF containing CDS features

Checks performed:
- Required files exist
- Chromosome naming styles are detectable
- Project transcript IDs can be extracted
- Each project transcript exists in the full annotation
- Each project transcript has at least one CDS in the full annotation
- Chromosomes used by project transcripts are compatible with FASTA and VCF
"""

import argparse
import gzip
import sys
from pathlib import Path
from collections import defaultdict
from typing import Set, Dict, Tuple


# -------------------------
# Error handling and utilities
# -------------------------

def die(error_message: str) -> None:
    sys.exit(f"[ERROR] {error_message}")


def open_maybe_gzip(file_path: str):
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path)


def check_exists(file_path: str, file_label: str) -> None:
    if not Path(file_path).exists():
        die(f"File not found - {file_label}: {file_path}")


# -------------------------
# Chromosome naming utilities
# -------------------------

def normalize_chromosome(chromosome: str) -> str:
    """
    Normalize chromosome names:
      chr5 -> 5
      Chr5 -> 5
      5    -> 5

    Leaves scaffold names unchanged.
    """
    if chromosome.startswith("Chr"):
        stripped = chromosome[3:]
        return stripped if stripped else chromosome
    if chromosome.startswith("chr"):
        stripped = chromosome[3:]
        return stripped if stripped else chromosome
    return chromosome


def detect_chromosome_naming_style(chromosome_names: Set[str]) -> str:
    """
    Detect predominant naming style, preferring chromosome-like names over scaffolds.
    """
    if not chromosome_names:
        return "empty"

    primary_names = [name for name in chromosome_names if not name.startswith("scaffold_")]
    names_to_check = primary_names if primary_names else list(chromosome_names)

    chrN_count = sum(1 for name in names_to_check if name.startswith("chr"))
    ChrN_count = sum(1 for name in names_to_check if name.startswith("Chr"))
    N_count = sum(1 for name in names_to_check if name.isdigit())

    if ChrN_count >= chrN_count and ChrN_count >= N_count:
        return "ChrN"
    if chrN_count >= ChrN_count and chrN_count >= N_count:
        return "chrN"
    if N_count >= chrN_count and N_count >= ChrN_count:
        return "N"
    return "mixed"


# -------------------------
# FASTA / VCF parsing
# -------------------------

def extract_fasta_chromosomes(fasta_path: str) -> Set[str]:
    chromosomes = set()
    try:
        with open_maybe_gzip(fasta_path) as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    chromosome = line[1:].split()[0]
                    chromosomes.add(chromosome)
    except Exception as e:
        die(f"Error reading FASTA file {fasta_path}: {e}")
    return chromosomes


def extract_vcf_chromosomes(vcf_path: str, max_lines: int = 200000) -> Set[str]:
    chromosomes = set()
    try:
        with open_maybe_gzip(vcf_path) as vcf_file:
            data_lines_processed = 0
            for line in vcf_file:
                if line.startswith("#"):
                    continue
                chromosome = line.split("\t")[0]
                chromosomes.add(chromosome)
                data_lines_processed += 1
                if data_lines_processed >= max_lines:
                    break
    except Exception as e:
        die(f"Error reading VCF file {vcf_path}: {e}")
    return chromosomes


# -------------------------
# GFF parsing
# -------------------------

def parse_project_transcripts(project_gff_path: str) -> Tuple[Set[str], Dict[str, str]]:
    """
    Parse the project transcript GFF and extract transcript IDs + chromosome map.

    Returns:
      - set of transcript IDs
      - dict transcript_id -> chromosome
    """
    transcripts = set()
    chromosome_map = {}

    try:
        with open(project_gff_path) as gff_file:
            for line in gff_file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                columns = line.split("\t")
                if len(columns) < 9:
                    continue

                chromosome, _, feature_type, _, _, _, _, _, attributes = columns

                if feature_type != "transcript":
                    continue
                if "ID=" not in attributes:
                    continue

                attr_dict = {}
                for field in attributes.split(";"):
                    if "=" in field:
                        key, value = field.split("=", 1)
                        attr_dict[key] = value

                transcript_id = attr_dict.get("ID")
                if transcript_id:
                    transcripts.add(transcript_id)
                    chromosome_map[transcript_id] = chromosome

    except Exception as e:
        die(f"Error parsing project GFF file {project_gff_path}: {e}")

    return transcripts, chromosome_map


def parse_full_annotation_for_transcripts(annotation_gff_path: str, target_transcripts: Set[str]) -> Tuple[Set[str], Dict[str, int], Dict[str, str]]:
    """
    Parse the full genome annotation and recover:
      - which target transcripts are found
      - CDS counts per transcript
      - chromosome map per transcript

    IMPORTANT:
    Assumes CDS Parent matches transcript ID for the transcript models of interest.
    """
    found_transcripts = set()
    cds_counts = defaultdict(int)
    chromosome_map = {}

    try:
        with open_maybe_gzip(annotation_gff_path) as gff_file:
            for line in gff_file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                columns = line.split("\t")
                if len(columns) < 9:
                    continue

                chromosome, _, feature_type, _, _, _, _, _, attributes = columns

                attr_dict = {}
                for field in attributes.split(";"):
                    if "=" in field:
                        key, value = field.split("=", 1)
                        attr_dict[key] = value

                if feature_type == "transcript":
                    transcript_id = attr_dict.get("ID")
                    if transcript_id in target_transcripts:
                        found_transcripts.add(transcript_id)
                        chromosome_map[transcript_id] = chromosome

                elif feature_type == "CDS":
                    parent_id = attr_dict.get("Parent")
                    if parent_id in target_transcripts:
                        cds_counts[parent_id] += 1

    except Exception as e:
        die(f"Error parsing annotation GFF file {annotation_gff_path}: {e}")

    return found_transcripts, cds_counts, chromosome_map


# -------------------------
# Main validation logic
# -------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate input files for CDS variant pipeline"
    )
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file (may be .gz)")
    parser.add_argument("--project-gff", required=True, help="Project transcript GFF (target transcripts only)")
    parser.add_argument("--annotation-gff", required=True, help="Full genome annotation GFF/GFF3")
    parser.add_argument("--vcf", required=True, help="VCF file with variants (may be .gz)")

    args = parser.parse_args()

    print("=" * 60)
    print("INPUT VALIDATION")
    print("=" * 60)

    # Step 1: file existence
    print("\n[1] Checking file existence...")
    check_exists(args.fasta, "FASTA")
    check_exists(args.project_gff, "project GFF")
    check_exists(args.annotation_gff, "annotation GFF")
    check_exists(args.vcf, "VCF")
    print("[✓] All required files exist")

    # Step 2: chromosome naming
    print("\n[2] Analyzing chromosome naming...")
    fasta_chromosomes = extract_fasta_chromosomes(args.fasta)
    vcf_chromosomes = extract_vcf_chromosomes(args.vcf)

    print(f"    FASTA naming style: {detect_chromosome_naming_style(fasta_chromosomes)}")
    print(f"    VCF naming style:   {detect_chromosome_naming_style(vcf_chromosomes)}")

    normalized_fasta = {normalize_chromosome(x) for x in fasta_chromosomes}
    normalized_vcf = {normalize_chromosome(x) for x in vcf_chromosomes}

    # Step 3: parse project transcripts
    print("\n[3] Reading project transcript definitions...")
    project_transcripts, project_chr_map = parse_project_transcripts(args.project_gff)

    if not project_transcripts:
        die("No transcript IDs found in project GFF")

    print(f"[✓] Found {len(project_transcripts)} target transcript(s) in project GFF")
    for tx in sorted(project_transcripts):
        print(f"    {tx}  chromosome={project_chr_map.get(tx, '?')}")

    # Step 4: validate against full annotation
    print("\n[4] Checking transcripts against full annotation...")
    found_transcripts, cds_counts, annotation_chr_map = parse_full_annotation_for_transcripts(
        args.annotation_gff,
        project_transcripts
    )

    missing_transcripts = sorted(project_transcripts - found_transcripts)
    if missing_transcripts:
        die(
            "The following project transcripts were not found in the full annotation: "
            + ", ".join(missing_transcripts)
        )

    print("[✓] All project transcripts were found in the full annotation")

    # Step 5: CDS validation
    print("\n[5] Checking CDS availability...")
    for tx in sorted(project_transcripts):
        cds_count = cds_counts.get(tx, 0)
        chromosome = annotation_chr_map.get(tx, project_chr_map.get(tx, "?"))
        print(f"    {tx:<30} chromosome={chromosome:<8} CDS features={cds_count}")
        if cds_count == 0:
            die(f"Transcript {tx} has no CDS features in the full annotation")

    print("[✓] All project transcripts have CDS features in the full annotation")

    # Step 6: chromosome compatibility
    print("\n[6] Checking chromosome compatibility...")
    transcript_normalized_chromosomes = {
        normalize_chromosome(annotation_chr_map.get(tx, project_chr_map[tx]))
        for tx in project_transcripts
    }

    missing_in_fasta = transcript_normalized_chromosomes - normalized_fasta
    if missing_in_fasta:
        die(f"Transcript chromosomes not found in FASTA: {', '.join(sorted(missing_in_fasta))}")
    else:
        print("    [✓] Transcript chromosomes are present in FASTA")

    missing_in_vcf = transcript_normalized_chromosomes - normalized_vcf
    if missing_in_vcf:
        print("    [!] WARNING: Some transcript chromosomes are not found in the VCF:")
        print(f"        {', '.join(sorted(missing_in_vcf))}")
        print("        This may indicate missing variants for these chromosomes.")
    else:
        print("    [✓] Transcript chromosomes are present in VCF")

    print("\n" + "=" * 60)
    print("✓ VALIDATION COMPLETED")
    print("=" * 60)


if __name__ == "__main__":
    main()
