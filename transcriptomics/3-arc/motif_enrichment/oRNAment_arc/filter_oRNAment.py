#!/usr/bin/env python3
"""
filter_oRNAment.py

Stream-filter a huge CSV by:
- ensembl_gene_id (exact match; multiple allowed)
- region (exact match; multiple allowed)
- score >= MIN
- unpaired_probability >= MIN

Always writes a header row. Works with or without an input header. Supports
.gz/.bz2/.xz. Indices default to the canonical layout if no header.

Canonical header (for --no-header):
  0 ensembl_gene_id
  1 ensembl_transcript_id
  2 gene_biotype
  3 transcript_biotype
  4 transcript_position
  5 RBP
  6 score
  7 unpaired_probability
  8 chromosome
  9 region
 10 exon_start
 11 exon_end
"""

import argparse
import csv
import gzip
import bz2
import lzma
import sys
from typing import Iterable, Set, TextIO, Optional

CANONICAL_HEADER = [
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "gene_biotype",
    "transcript_biotype",
    "transcript_position",
    "RBP",
    "score",
    "unpaired_probability",
    "chromosome",
    "region",
    "exon_start",
    "exon_end",
]

def open_maybe_compressed(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode="rt", newline="")
    if path.endswith(".bz2"):
        return bz2.open(path, mode="rt", newline="")
    if path.endswith(".xz") or path.endswith(".lzma"):
        return lzma.open(path, mode="rt", newline="")
    return open(path, mode="r", newline="")

def load_set(cli_vals: Iterable[str], file_path: Optional[str]) -> Set[str]:
    s = set(cli_vals or [])
    if file_path:
        with open(file_path, "r") as fh:
            for line in fh:
                v = line.strip()
                if v:
                    s.add(v)
    return s

def detect_idx(header: list[str], name: str, default_idx: int) -> int:
    try:
        return header.index(name)
    except ValueError:
        if 0 <= default_idx < len(header):
            return default_idx
        raise SystemExit(f"Column '{name}' not found and default index {default_idx} is invalid for this file.")

def safe_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None

def stream_filter(
    infile: TextIO,
    out: TextIO,
    has_header: bool,
    gene_ids: Set[str],
    regions: Set[str],
    gene_col_name: Optional[str],
    region_col_name: Optional[str],
    gene_col_index: Optional[int],
    region_col_index: Optional[int],
    min_score: Optional[float],
    min_unpaired: Optional[float],
    score_col_name: Optional[str],
    unpaired_col_name: Optional[str],
    score_col_index: Optional[int],
    unpaired_col_index: Optional[int],
) -> tuple[int, int]:
    reader = csv.reader(infile)
    writer = csv.writer(out, lineterminator="\n")

    # Header + column indices
    if has_header:
        try:
            header = next(reader)
        except StopIteration:
            writer.writerow(CANONICAL_HEADER)
            return (0, 0)
        writer.writerow(header)
        g_idx = gene_col_index if gene_col_index is not None else detect_idx(header, gene_col_name or "ensembl_gene_id", 0)
        r_idx = region_col_index if region_col_index is not None else detect_idx(header, region_col_name or "region", 9)
        s_idx = score_col_index if score_col_index is not None else detect_idx(header, score_col_name or "score", 6)
        u_idx = unpaired_col_index if unpaired_col_index is not None else detect_idx(header, unpaired_col_name or "unpaired_probability", 7)
    else:
        writer.writerow(CANONICAL_HEADER)
        g_idx = 0 if gene_col_index is None else gene_col_index
        r_idx = 9 if region_col_index is None else region_col_index
        s_idx = 6 if score_col_index is None else score_col_index
        u_idx = 7 if unpaired_col_index is None else unpaired_col_index

    matched = 0
    skipped_non_numeric = 0

    use_gene = len(gene_ids) > 0
    use_region = len(regions) > 0
    use_score = min_score is not None
    use_unpaired = min_unpaired is not None

    for row in reader:
        if not row:
            continue
        max_need = max(g_idx, r_idx, s_idx, u_idx)
        if max_need >= len(row):
            continue

        if use_gene and row[g_idx] not in gene_ids:
            continue
        if use_region and row[r_idx] not in regions:
            continue

        if use_score or use_unpaired:
            sval = safe_float(row[s_idx])
            uval = safe_float(row[u_idx])
            if sval is None or uval is None:
                skipped_non_numeric += 1
                continue
            if use_score and sval < min_score:
                continue
            if use_unpaired and uval < min_unpaired:
                continue

        writer.writerow(row)
        matched += 1

    return (matched, skipped_non_numeric)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stream-filter a large CSV by gene, region, and numeric thresholds.")
    p.add_argument("csv", help="Input CSV (.csv[.gz|.bz2|.xz])")
    p.add_argument("-o", "--out", help="Output CSV path (default: stdout)")

    # Header state
    p.add_argument("--no-header", action="store_true", help="Input has no header row")

    # Gene filters
    p.add_argument("--gene-id", action="append", default=[], help="ensembl_gene_id to match (repeatable)")
    p.add_argument("--ids-file", help="File with one ensembl_gene_id per line")
    p.add_argument("--gene-col-name", help="Gene column name (default: ensembl_gene_id)")
    p.add_argument("--gene-col-index", type=int, help="Gene column index (0-based; overrides name)")

    # Region filters
    p.add_argument("--region", action="append", default=[], help="region to match (repeatable)")
    p.add_argument("--regions-file", help="File with one region per line")
    p.add_argument("--region-col-name", help="Region column name (default: region)")
    p.add_argument("--region-col-index", type=int, help="Region column index (0-based; overrides name)")

    # Thresholds
    p.add_argument("--min-score", type=float, help="Keep rows with score >= this value")
    p.add_argument("--min-unpaired", type=float, help="Keep rows with unpaired_probability >= this value")

    # Optional column overrides for numeric fields
    p.add_argument("--score-col-name", help="Score column name (default: score)")
    p.add_argument("--unpaired-col-name", help="Unpaired probability column name (default: unpaired_probability)")
    p.add_argument("--score-col-index", type=int, help="Score column index (0-based; overrides name)")
    p.add_argument("--unpaired-col-index", type=int, help="Unpaired probability column index (0-based; overrides name)")

    return p.parse_args()

def main():
    args = parse_args()

    gene_ids = load_set(args.gene_id, args.ids_file)
    regions  = load_set(args.region, args.regions_file)

    if args.out:
        outfh = open(args.out, "w", newline="")
        close_out = True
    else:
        outfh = sys.stdout
        close_out = False

    try:
        with open_maybe_compressed(args.csv) as infh:
            matched, skipped_bad = stream_filter(
                infile=infh,
                out=outfh,
                has_header=not args.no_header,
                gene_ids=gene_ids,
                regions=regions,
                gene_col_name=args.gene_col_name,
                region_col_name=args.region_col_name,
                gene_col_index=args.gene_col_index,
                region_col_index=args.region_col_index,
                min_score=args.min_score,
                min_unpaired=args.min_unpaired,
                score_col_name=args.score_col_name,
                unpaired_col_name=args.unpaired_col_name,
                score_col_index=args.score_col_index,
                unpaired_col_index=args.unpaired_col_index,
            )
    finally:
        if close_out:
            outfh.close()

    # Summary to stderr
    msg = f"Matched rows: {matched}"
    if skipped_bad:
        msg += f" | Skipped non-numeric rows: {skipped_bad}"
    print(msg, file=sys.stderr)

if __name__ == "__main__":
    main()