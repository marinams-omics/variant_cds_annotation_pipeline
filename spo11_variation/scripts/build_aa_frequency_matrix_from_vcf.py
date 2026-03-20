#!/usr/bin/env python3
"""
Build amino-acid frequency tables per protein position from:

  - CDS FASTA
  - annotated variant effects TSV
  - genotype VCF for the same transcript

Outputs:
  1) long TSV:
       aa_pos, ref_codon, ref_aa, aa, count, pct, n_called, n_uncalled, variant_ids
  2) matrix TSV:
       aa <tab> pos1 <tab> pos2 <tab> ...

Important logic:
- Reconstruction is done per codon, per sample.
- Multiple SNPs in the same codon are applied jointly to the sample codon.
- Samples with ambiguous genotypes in any variant affecting the codon are excluded
  from the denominator for that codon.

Genotype handling:
- 0/0 or 0|0 -> reference
- 1/1 or 1|1 -> alternate
- 0/1, 1/0, 1/2, 0/2, etc. -> ambiguous (excluded)
- missing -> excluded

This conservative approach is appropriate for mostly inbred HapMap lines.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY*")


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

    seq = "".join(seq_parts).upper()
    if not seq:
        die(f"No sequence found in: {path}")

    return header, seq


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")


def parse_effects(path: str):
    """
    Parse variant_effects.tsv.

    Returns:
      variant_info[var_id] = {
        codon_idx: int,
        pos_in_codon: int,
        ref_codon: str,
        alt_codon: str,
        ref_aa: str,
        alt_aa: str,
        effect: str,
        alt_base: str,
      }

      codon_to_vars[codon_idx] = [var_id, ...]
    """
    variant_info = {}
    codon_to_vars = defaultdict(list)

    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        expected = [
            "chrom", "pos", "id", "ref", "alt",
            "cds_pos", "codon_idx", "pos_in_codon",
            "ref_codon", "alt_codon", "ref_aa", "alt_aa", "effect"
        ]
        if header[:13] != expected:
            die(f"Unexpected header in effects file: {path}")

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 13:
                die(f"Invalid line in effects file {path}:{line_num}")

            var_id = cols[2]
            codon_idx = cols[6]
            pos_in_codon = cols[7]
            ref_codon = cols[8]
            alt_codon = cols[9]
            ref_aa = cols[10]
            alt_aa = cols[11]
            effect = cols[12]

            if codon_idx == "NA" or pos_in_codon == "NA":
                continue

            try:
                codon_idx = int(codon_idx)
                pos_in_codon = int(pos_in_codon)
            except ValueError:
                continue

            if len(ref_codon) != 3 or len(alt_codon) != 3:
                continue

            alt_base = alt_codon[pos_in_codon - 1]

            variant_info[var_id] = {
                "codon_idx": codon_idx,
                "pos_in_codon": pos_in_codon,
                "ref_codon": ref_codon,
                "alt_codon": alt_codon,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "effect": effect,
                "alt_base": alt_base,
            }
            codon_to_vars[codon_idx].append(var_id)

    return variant_info, codon_to_vars


def classify_gt(gt_field: str) -> str:
    """
    Conservative genotype classification.

    Returns:
      ref       -> 0/0 or 0|0
      alt       -> 1/1 or 1|1 (or all same nonzero allele)
      missing   -> ./. or .|. or any missing allele
      ambiguous -> heterozygous or mixed nonzero alleles

    This is conservative and suitable for mostly inbred lines.
    """
    if gt_field is None or gt_field == "":
        return "missing"

    gt = gt_field.split(":")[0]

    if gt in {".", "./.", ".|."}:
        return "missing"

    sep = "/" if "/" in gt else "|" if "|" in gt else None
    if sep is None:
        # haploid-like unexpected format
        if gt == ".":
            return "missing"
        if gt == "0":
            return "ref"
        return "alt"

    alleles = gt.split(sep)

    if any(a == "." for a in alleles):
        return "missing"

    if all(a == "0" for a in alleles):
        return "ref"

    # all same nonzero allele -> alt
    if len(set(alleles)) == 1 and next(iter(set(alleles))) != "0":
        return "alt"

    # anything mixed is ambiguous
    return "ambiguous"


def parse_vcf_genotypes(vcf_path: str):
    """
    Parse VCF and return:
      samples: [sample1, sample2, ...]
      gt_by_var[var_id][sample] = one of {ref, alt, missing, ambiguous}
    """
    samples = None
    gt_by_var = defaultdict(dict)

    with open(vcf_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) < 10:
                    die(f"VCF has no sample columns: {vcf_path}")
                samples = cols[9:]
                continue

            cols = line.split("\t")
            if len(cols) < 10:
                continue

            var_id = cols[2]
            sample_fields = cols[9:]

            if samples is None:
                die(f"VCF header not found before records: {vcf_path}")

            for sample, gt_field in zip(samples, sample_fields):
                gt_by_var[var_id][sample] = classify_gt(gt_field)

    if samples is None:
        die(f"No #CHROM header found in VCF: {vcf_path}")

    return samples, gt_by_var


def main():
    ap = argparse.ArgumentParser(description="Build AA frequency matrix from CDS + effects + genotype VCF.")
    ap.add_argument("--cds", required=True, help="Transcript CDS FASTA")
    ap.add_argument("--effects", required=True, help="Transcript variant_effects.tsv")
    ap.add_argument("--vcf", required=True, help="Transcript CDS genotype VCF")
    ap.add_argument("--out-long", required=True, help="Output long TSV")
    ap.add_argument("--out-matrix", required=True, help="Output matrix TSV")
    args = ap.parse_args()

    out_long = Path(args.out_long)
    out_matrix = Path(args.out_matrix)
    out_long.parent.mkdir(parents=True, exist_ok=True)
    out_matrix.parent.mkdir(parents=True, exist_ok=True)

    _, cds_seq = read_fasta_one(args.cds)

    if len(cds_seq) % 3 != 0:
        die(f"CDS length is not divisible by 3: {len(cds_seq)}")

    n_codons = len(cds_seq) // 3
    ref_codons = {}
    ref_aas = {}

    for codon_idx in range(1, n_codons + 1):
        start = (codon_idx - 1) * 3
        codon = cds_seq[start:start + 3]
        aa = translate_codon(codon)
        ref_codons[codon_idx] = codon
        ref_aas[codon_idx] = aa

    variant_info, codon_to_vars = parse_effects(args.effects)
    samples, gt_by_var = parse_vcf_genotypes(args.vcf)

    n_samples = len(samples)

    # aa_counts[pos][aa] = count
    aa_counts = defaultdict(lambda: defaultdict(int))
    n_called_by_pos = defaultdict(int)
    n_uncalled_by_pos = defaultdict(int)
    variant_ids_by_pos = defaultdict(list)

    for codon_idx in range(1, n_codons + 1):
        ref_codon = ref_codons[codon_idx]
        vars_in_codon = codon_to_vars.get(codon_idx, [])

        if not vars_in_codon:
            # all samples have reference codon/aa
            aa_counts[codon_idx][ref_aas[codon_idx]] = n_samples
            n_called_by_pos[codon_idx] = n_samples
            n_uncalled_by_pos[codon_idx] = 0
            continue

        for var_id in vars_in_codon:
            variant_ids_by_pos[codon_idx].append(var_id)

        for sample in samples:
            sample_codon = list(ref_codon)
            sample_uncalled = False

            for var_id in vars_in_codon:
                gt_state = gt_by_var.get(var_id, {}).get(sample, "missing")

                if gt_state in {"missing", "ambiguous"}:
                    sample_uncalled = True
                    break

                if gt_state == "alt":
                    pos_in_codon = variant_info[var_id]["pos_in_codon"]
                    alt_base = variant_info[var_id]["alt_base"]
                    sample_codon[pos_in_codon - 1] = alt_base

            if sample_uncalled:
                n_uncalled_by_pos[codon_idx] += 1
                continue

            aa = translate_codon("".join(sample_codon))
            aa_counts[codon_idx][aa] += 1
            n_called_by_pos[codon_idx] += 1

    # Write long table
    with open(args.out_long, "w") as out:
        out.write(
            "\t".join([
                "aa_pos",
                "ref_codon",
                "ref_aa",
                "aa",
                "count",
                "pct",
                "n_called",
                "n_uncalled",
                "variant_ids",
            ]) + "\n"
        )

        for aa_pos in range(1, n_codons + 1):
            n_called = n_called_by_pos.get(aa_pos, 0)
            n_uncalled = n_uncalled_by_pos.get(aa_pos, 0)
            variant_ids = ",".join(sorted(set(variant_ids_by_pos.get(aa_pos, [])))) or "NA"

            # write all AA rows, including zeros, for completeness
            for aa in AA_ORDER:
                count = aa_counts[aa_pos].get(aa, 0)
                pct = (count / n_called * 100.0) if n_called > 0 else 0.0
                out.write(
                    "\t".join([
                        str(aa_pos),
                        ref_codons[aa_pos],
                        ref_aas[aa_pos],
                        aa,
                        str(count),
                        f"{pct:.4f}",
                        str(n_called),
                        str(n_uncalled),
                        variant_ids,
                    ]) + "\n"
                )

    # Write matrix (AA rows x positions columns)
    with open(args.out_matrix, "w") as out:
        header = ["aa"] + [str(i) for i in range(1, n_codons + 1)]
        out.write("\t".join(header) + "\n")

        for aa in AA_ORDER:
            row = [aa]
            for aa_pos in range(1, n_codons + 1):
                n_called = n_called_by_pos.get(aa_pos, 0)
                count = aa_counts[aa_pos].get(aa, 0)
                pct = (count / n_called * 100.0) if n_called > 0 else 0.0
                row.append(f"{pct:.4f}")
            out.write("\t".join(row) + "\n")

    variable_positions = 0
    for aa_pos in range(1, n_codons + 1):
        ref_aa = ref_aas[aa_pos]
        n_called = n_called_by_pos.get(aa_pos, 0)
        ref_count = aa_counts[aa_pos].get(ref_aa, 0)
        ref_pct = (ref_count / n_called * 100.0) if n_called > 0 else 0.0
        if ref_pct < 100.0:
            variable_positions += 1

    print("== BUILD AA FREQUENCY MATRIX FROM VCF ==")
    print(f"[OK] codon_positions={n_codons}")
    print(f"[OK] n_samples={n_samples}")
    print(f"[OK] variants_in_effects={len(variant_info)}")
    print(f"[OK] variable_positions={variable_positions}")
    print(f"[OK] wrote long:   {args.out_long}")
    print(f"[OK] wrote matrix: {args.out_matrix}")


if __name__ == "__main__":
    main()
