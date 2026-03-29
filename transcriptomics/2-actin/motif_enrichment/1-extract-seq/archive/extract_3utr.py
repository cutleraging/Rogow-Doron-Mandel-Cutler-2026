#!/usr/bin/env python3
"""
extract_3utr.py – fetch 3'UTR sequences by gene_id

Usage
-----
python extract_3utr.py <annotations.gtf> <reference.fa> <gene_ids.txt> > gene_3utr.fa

• annotations.gtf  – GTF file (must contain feature type “three_prime_utr” or “three_prime_UTR”)
• reference.fa     – genome FASTA (indexed with samtools faidx or `pysam.faidx`)
• gene_ids.txt     – one gene_id per line (exactly as in the GTF)

Outputs FASTA records (one per gene).  Requires `pysam`.
"""

import argparse, sys, re
from collections import defaultdict
import pysam

_RC = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("gtf")
    p.add_argument("fasta")
    p.add_argument("genes")
    return p.parse_args()

def load_genes(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}

def gene_id(attr_str):
    m = re.search(r'gene_id "([^"]+)"', attr_str)
    return m.group(1) if m else None

def main():
    args = parse_args()
    targets = load_genes(args.genes)
    utrs = defaultdict(list)          # gene_id → [(chrom,start,end,strand), …]

    # 1. collect 3'UTR intervals
    with open(args.gtf) as gtf:
        for line in gtf:
            if line.startswith("#"): 
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = line.rstrip().split("\t")
            if feature.lower() not in {"three_prime_utr", "three_prime_utr", "3utr", "3_prime_utr"}:
                continue
            gid = gene_id(attrs)
            if gid in targets:
                utrs[gid].append((chrom, int(start), int(end), strand))

    if not utrs:
        sys.exit("No matching 3'UTRs found.",)

    # 2. fetch sequences
    fasta = pysam.FastaFile(args.fasta)
    for gid, segs in utrs.items():
        # sort by genomic coordinate; maintain strand order later
        segs.sort(key=lambda x: x[1])
        seq_parts = []
        for chrom, s, e, strand in segs:
            # GTF is 1‑based inclusive; pysam is 0‑based half‑open
            seq_parts.append(fasta.fetch(chrom, s-1, e))
        seq = "".join(seq_parts)
        if segs[0][3] == "-":                      # reverse‑complement if needed
            seq = revcomp(seq)
        print(f">{gid}\n{seq}")

if __name__ == "__main__":
    main()