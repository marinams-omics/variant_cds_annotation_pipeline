#!/usr/bin/env python3
"""
Annotate CDS variants using:
  - a transcript CDS FASTA
  - a CDS blockmap
  - a variants TSV (chrom pos id ref alt ...)

Outputs TSV columns:
  chrom  pos  id  ref  alt  cds_pos  codon_idx  pos_in_codon
  ref_codon  alt_codon  ref_aa  alt_aa  effect

Notes:
- Handles transcript strand correctly
- For minus-strand transcripts:
    * genomic REF/ALT are complemented before comparison to CDS
    * genomic position is mapped to transcript CDS position
- Non simple SNPs (e.g. indels, multiallelic) are labeled as non_SNP
"""

import argparse
import gzip
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


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def complement_base(base: str) -> str:
    comp = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    return comp.get(base.upper(), "N")


def read_fasta_one(path: str):
    header = None
    seq_parts = []

    with open_maybe_gzip(path) as f:
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


def read_blockmap(path: str):
    """
    Reads blockmap TSV with columns:
      block chrom start0 end1 strand len cds_start cds_end

    Returns:
      rows: list of dicts
      total_len: int
      strand: '+' or '-'
    """
    rows = []
    header = None

    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")

            if line_num == 1:
                header = cols
                expected = ["block", "chrom", "start0", "end1", "strand", "len", "cds_start", "cds_end"]
                if cols[:8] != expected:
                    die(f"Unexpected blockmap header in {path}: {cols}")
                continue

            if len(cols) < 8:
                die(f"Invalid blockmap line {line_num}: {line}")

            block, chrom, start0, end1, strand, blen, cds_start, cds_end = cols[:8]

            row = {
                "block": block,
                "chrom": chrom,
                "start0": int(start0),
                "end1": int(end1),
                "strand": strand,
                "len": int(blen),
                "cds_start": int(cds_start),
                "cds_end": int(cds_end),
            }
            rows.append(row)

    if not rows:
        die(f"No blockmap rows found in: {path}")

    strands = {r["strand"] for r in rows}
    if len(strands) != 1:
        die(f"Blockmap contains multiple strands: {sorted(strands)}")

    strand = next(iter(strands))
    total_len = max(r["cds_end"] for r in rows)

    return rows, total_len, strand


def genomic_to_cds_pos(chrom: str, pos1: int, block_rows, total_len: int):
    """
    Map genomic 1-based position to transcript CDS 1-based position.

    block_rows are assumed to be in genomic-order blockmap coordinates.
    For plus strand:
        genomic_concat_pos = cds_start + offset
        cds_pos = genomic_concat_pos
    For minus strand:
        genomic_concat_pos = cds_start + offset
        cds_pos = total_len - genomic_concat_pos + 1
    """
    for row in block_rows:
        if row["chrom"] != chrom:
            continue

        start0 = row["start0"]
        end1 = row["end1"]

        # BED inclusion: start0 < pos1 <= end1
        if start0 < pos1 <= end1:
            offset = pos1 - (start0 + 1)  # 0-based within block
            genomic_concat_pos = row["cds_start"] + offset

            if row["strand"] == "+":
                cds_pos = genomic_concat_pos
            else:
                cds_pos = total_len - genomic_concat_pos + 1

            return cds_pos, row["strand"]

    return None, None


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")


def annotate_variant(chrom, pos1, vid, ref, alt, cds_seq, block_rows, total_len,
                     strict_ref=False, warn_counter=None, warn_limit=10):
    """
    Annotate a single variant row.
    Returns dict with output fields.
    """
    cds_pos, strand = genomic_to_cds_pos(chrom, pos1, block_rows, total_len)

    # Not in CDS: should not happen if input was subset correctly, but keep robust
    if cds_pos is None:
        return {
            "chrom": chrom,
            "pos": pos1,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "cds_pos": "NA",
            "codon_idx": "NA",
            "pos_in_codon": "NA",
            "ref_codon": "NA",
            "alt_codon": "NA",
            "ref_aa": "NA",
            "alt_aa": "NA",
            "effect": "outside_CDS"
        }

    # Only annotate simple biallelic SNPs here
    is_simple_snp = (len(ref) == 1 and len(alt) == 1 and "," not in alt)

    if not is_simple_snp:
        return {
            "chrom": chrom,
            "pos": pos1,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "cds_pos": cds_pos,
            "codon_idx": "NA",
            "pos_in_codon": "NA",
            "ref_codon": "NA",
            "alt_codon": "NA",
            "ref_aa": "NA",
            "alt_aa": "NA",
            "effect": "non_SNP"
        }

    ref_base = ref.upper()
    alt_base = alt.upper()

    # Convert to transcript/CDS orientation if minus strand
    if strand == "-":
        ref_base = complement_base(ref_base)
        alt_base = complement_base(alt_base)

    if cds_pos < 1 or cds_pos > len(cds_seq):
        die(f"Mapped CDS position out of range: cds_pos={cds_pos}, len(cds_seq)={len(cds_seq)}")

    cds_base = cds_seq[cds_pos - 1]

    if cds_base != ref_base:
        msg = (
            f"REF mismatch at {chrom}:{pos1} {vid} | "
            f"CDS base={cds_base}, expected REF={ref_base}, strand={strand}, cds_pos={cds_pos}"
        )
        if strict_ref:
            die(msg)
        else:
            if warn_counter is not None and warn_counter[0] < warn_limit:
                print(f"[WARNING] {msg}", file=sys.stderr)
                warn_counter[0] += 1

    codon_idx = (cds_pos - 1) // 3 + 1
    pos_in_codon = (cds_pos - 1) % 3 + 1

    codon_start = (codon_idx - 1) * 3
    codon_end = codon_start + 3

    if codon_end > len(cds_seq):
        return {
            "chrom": chrom,
            "pos": pos1,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "cds_pos": cds_pos,
            "codon_idx": codon_idx,
            "pos_in_codon": pos_in_codon,
            "ref_codon": "NA",
            "alt_codon": "NA",
            "ref_aa": "NA",
            "alt_aa": "NA",
            "effect": "non_SNP"
        }

    ref_codon = cds_seq[codon_start:codon_end]
    alt_codon_list = list(ref_codon)
    alt_codon_list[pos_in_codon - 1] = alt_base
    alt_codon = "".join(alt_codon_list)

    ref_aa = translate_codon(ref_codon)
    alt_aa = translate_codon(alt_codon)

    if ref_aa == alt_aa:
        effect = "synonymous"
    elif alt_aa == "*":
        effect = "nonsense"
    else:
        effect = "missense"

    return {
        "chrom": chrom,
        "pos": pos1,
        "id": vid,
        "ref": ref,
        "alt": alt,
        "cds_pos": cds_pos,
        "codon_idx": codon_idx,
        "pos_in_codon": pos_in_codon,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "effect": effect
    }


def main():
    ap = argparse.ArgumentParser(description="Annotate CDS variants from blockmap + CDS FASTA + variants TSV.")
    ap.add_argument("--blockmap", required=True, help="Input blockmap TSV")
    ap.add_argument("--cds", required=True, help="Input CDS FASTA")
    ap.add_argument("--variants", required=True, help="Input variants TSV (chrom pos id ref alt ...)")
    ap.add_argument("--strict-ref", action="store_true", help="Abort on REF/CDS mismatches")
    ap.add_argument("--warn-mismatch", type=int, default=10, help="Maximum number of mismatch warnings to print")
    args = ap.parse_args()

    _, cds_seq = read_fasta_one(args.cds)
    block_rows, total_len, strand = read_blockmap(args.blockmap)

    print("\t".join([
        "chrom", "pos", "id", "ref", "alt",
        "cds_pos", "codon_idx", "pos_in_codon",
        "ref_codon", "alt_codon", "ref_aa", "alt_aa", "effect"
    ]))

    warn_counter = [0]

    with open(args.variants, "r") as f:
        first = True
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")

            # Skip header if present
            if first:
                first = False
                try:
                    int(cols[1])
                except Exception:
                    continue

            if len(cols) < 5:
                die(f"Variants line {line_num} has fewer than 5 columns: {line}")

            chrom = cols[0]
            try:
                pos1 = int(cols[1])
            except ValueError:
                die(f"Invalid position in variants line {line_num}: {line}")

            vid = cols[2]
            ref = cols[3]
            alt = cols[4]

            ann = annotate_variant(
                chrom=chrom,
                pos1=pos1,
                vid=vid,
                ref=ref,
                alt=alt,
                cds_seq=cds_seq,
                block_rows=block_rows,
                total_len=total_len,
                strict_ref=args.strict_ref,
                warn_counter=warn_counter,
                warn_limit=args.warn_mismatch
            )

            print("\t".join(str(ann[k]) for k in [
                "chrom", "pos", "id", "ref", "alt",
                "cds_pos", "codon_idx", "pos_in_codon",
                "ref_codon", "alt_codon", "ref_aa", "alt_aa", "effect"
            ]))


if __name__ == "__main__":
    main()
