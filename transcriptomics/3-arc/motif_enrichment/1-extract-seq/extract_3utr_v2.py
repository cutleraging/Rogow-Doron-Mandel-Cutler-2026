#!/usr/bin/env python3
import argparse, sys, re
from collections import defaultdict
import pysam

_RC = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s): return s.translate(_RC)[::-1]

def parse_args():
    p = argparse.ArgumentParser(description="Extract 3'UTR sequences per transcript (or longest per gene).")
    p.add_argument("gtf"); p.add_argument("fasta"); p.add_argument("ids")  # ids can be gene_id or transcript_id
    p.add_argument("--id-type", choices=["gene_id","transcript_id"], default="gene_id")
    p.add_argument("--mode", choices=["per-transcript","longest-per-gene"], default="per-transcript")
    p.add_argument("--wrap", type=int, default=60)
    return p.parse_args()

def get_attr(s, key):
    m = re.search(rf'{key}\s+"([^"]+)"', s)
    return m.group(1) if m else None

def wrap(s, w): 
    return "\n".join(s[i:i+w] for i in range(0,len(s),w)) if w>0 else s

def main():
    a = parse_args()
    targets = {line.strip() for line in open(a.ids) if line.strip()}
    utrs = defaultdict(list)  # transcript_id -> [(chrom,start,end,strand)]
    t2g  = {}                 # transcript_id -> gene_id

    # Collect UTR intervals
    with open(a.gtf) as gtf:
        for line in gtf:
            if not line or line[0] == "#": 
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = line.rstrip().split("\t")
            f = feature.lower()
            # Accept common encodings of 3'UTR
            if f not in {"three_prime_utr","3_prime_utr"}: 
                continue
            tid = get_attr(attrs, "transcript_id")
            gid = get_attr(attrs, "gene_id")
            if not tid or not gid:
                continue
            keep_id = gid if a.id_type == "gene_id" else tid
            if keep_id in targets:
                utrs[tid].append((chrom, int(start), int(end), strand))
                t2g[tid] = gid

    if not utrs:
        sys.exit("No matching 3'UTRs found.")

    fa = pysam.FastaFile(a.fasta)

    def fetch_transcript_utr(tid, segs):
        # sanity: same chrom/strand
        chroms = {c for c,_,_,_ in segs}; strands = {s for _,_,_,s in segs}
        if len(chroms) != 1 or len(strands) != 1:
            return None, None  # skip malformed
        chrom = next(iter(chroms)); strand = next(iter(strands))
        segs = sorted(segs, key=lambda x: x[1])
        # merge overlaps
        merged = []
        for c,s,e,_ in segs:
            if not merged or s > merged[-1][2]:
                merged.append([c,s,e])
            else:
                merged[-1][2] = max(merged[-1][2], e)
        # fetch & orient
        seq = "".join(fa.fetch(chrom, s-1, e) for c,s,e in merged)
        if strand == "-":
            seq = revcomp(seq)
        return seq, strand

    records = []  # (name, seq)

    if a.mode == "per-transcript":
        for tid, segs in utrs.items():
            seq, strand = fetch_transcript_utr(tid, segs)
            if seq:
                name = f"{tid}|gene:{t2g[tid]}|strand:{strand}"
                records.append((name, seq))
    else:  # longest-per-gene among its transcripts present
        by_gene = defaultdict(list)
        for tid, segs in utrs.items():
            seq, strand = fetch_transcript_utr(tid, segs)
            if seq:
                by_gene[t2g[tid]].append((tid, strand, seq))
        for gid, L in by_gene.items():
            tid, strand, seq = max(L, key=lambda x: len(x[2]))
            records.append((f"{gid}|transcript:{tid}|strand:{strand}", seq))

    for name, seq in records:
        print(f">{name}\n{wrap(seq, a.wrap)}")

if __name__ == "__main__":
    main()