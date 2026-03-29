#!/usr/bin/env python3
"""
pwm_to_meme_clean_noblank.py

Concatenate RBP PWM files into a MEME-format file.
Removes:
  - 'PO A C G U' header lines
  - row indices (1,2,3,...) at line start
  - the blank line between MOTIF and letter-probability matrix
"""

import argparse, csv, os, re, sys
from collections import OrderedDict

RE_ROW_INDEX = re.compile(r"^\s*\d+\s+(.+)$")
RE_PO_HEADER = re.compile(r"^\s*PO\b", re.IGNORECASE)

def parse_args():
    p = argparse.ArgumentParser(description="Combine PWMs into MEME format (no blank line).")
    p.add_argument("--filtered", required=True)
    p.add_argument("--mapping", required=True)
    p.add_argument("--pwm-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--no-header", action="store_true")
    p.add_argument("--rbp-col-name", default="RBP")
    p.add_argument("--rbp-col-index", type=int)
    p.add_argument("--fail-on-missing", action="store_true")
    p.add_argument("--nsites", type=str, default="0")
    p.add_argument("--evalue", type=str, default="0")
    return p.parse_args()

def load_mapping(path):
    m = {}
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        for row in rdr:
            if not row:
                continue
            if len(row) < 2:
                parts = [p.strip() for p in row[0].split(",", 1)]
                if len(parts) != 2:
                    continue
                k, v = parts
            else:
                k, v = row[0].strip(), row[1].strip()
            v = v.strip().strip('"').strip("'")
            try:
                m[int(k)] = v
            except ValueError:
                continue
    return m

def iter_rbp_ids(filtered, colname, colindex, has_header=True):
    seen = OrderedDict()
    with open(filtered, newline="") as fh:
        rdr = csv.reader(fh)
        if has_header:
            header = next(rdr, [])
            idx = colindex if colindex is not None else header.index(colname)
        else:
            idx = 5 if colindex is None else colindex
        for row in rdr:
            if not row or idx >= len(row): 
                continue
            val = row[idx].strip().strip('"').strip("'")
            if not val:
                continue
            try:
                rid = int(val)
            except ValueError:
                continue
            if rid not in seen:
                seen[rid] = None
    return seen.keys()

def sanitize_name(name): 
    return "_".join(name.replace('"','').replace("'","").split())

def clean_pwm_lines(path):
    cleaned = []
    with open(path) as fh:
        for raw in fh:
            s = raw.rstrip("\n")
            if not s.strip() or RE_PO_HEADER.match(s):
                continue
            m = RE_ROW_INDEX.match(s)
            if m:
                s = m.group(1)
            if len(s.strip().split()) >= 4:
                cleaned.append(s)
    return cleaned, len(cleaned)

def main():
    a = parse_args()
    mapping = load_mapping(a.mapping)
    rbp_ids = list(iter_rbp_ids(a.filtered, a.rbp_col_name, a.rbp_col_index, not a.no_header))
    if not rbp_ids:
        print("No RBP IDs found.", file=sys.stderr)
        if a.fail_on_missing: sys.exit(1)

    with open(a.out, "w") as out:
        out.write("MEME version 5.5.8\n\n")
        out.write("ALPHABET= ACGU\n\n")
        out.write("Background letter frequencies (from uniform background):\n")
        out.write("A 0.25000 C 0.25000 G 0.25000 U 0.25000 \n\n")

        for rid in rbp_ids:
            name = mapping.get(rid)
            if name is None:
                print(f"[warn] Missing mapping for {rid}", file=sys.stderr)
                if a.fail_on_missing: sys.exit(2)
                continue
            pwm = os.path.join(a.pwm_dir, f"{rid:03d}.PWM")
            if not os.path.isfile(pwm):
                print(f"[warn] Missing PWM: {pwm}", file=sys.stderr)
                if a.fail_on_missing: sys.exit(3)
                continue

            lines, w = clean_pwm_lines(pwm)
            if w == 0:
                print(f"[warn] Empty PWM: {pwm}", file=sys.stderr)
                if a.fail_on_missing: sys.exit(4)
                continue

            motif = sanitize_name(name)
            out.write(f"MOTIF {motif}\n")
            out.write(f"letter-probability matrix: alength= 4 w= {w} nsites= {a.nsites} E= {a.evalue}\n")
            for ln in lines:
                out.write(f"{ln}\n")
            out.write("\n")

if __name__ == "__main__":
    main()