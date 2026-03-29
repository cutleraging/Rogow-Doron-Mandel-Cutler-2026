#!/usr/bin/env python3
"""
extract_by_gene.py

Extracts sequences from a FASTA file whose header contains a GN= field
matching any gene name in a supplied list.

Usage
-----
python extract_by_gene.py genes.txt proteins.fasta output.fasta

  genes.txt        one gene name per line (Ighv5-4, Igkv10-96, …)
  proteins.fasta   FASTA database to filter
  output.fasta     resulting FASTA with matching entries
"""
#!/usr/bin/env python3
"""
extract_by_gene.py

Extract sequences from a FASTA file whose header contains a GN= field
matching any gene name in a supplied list.

Gene list format
----------------
One line per entry.  If a line contains multiple names separated by
semicolons (e.g. Calm1;Calm2;Calm3), only the **first** name is used.

Usage
-----
python extract_by_gene.py genes.txt proteins.fasta output.fasta
"""
import sys
import re
from pathlib import Path


def load_gene_set(gene_file):
    """Return a set of first‐listed gene names from each non‑empty line."""
    genes = set()
    with open(gene_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            first = line.split(";", 1)[0].strip()
            if first:
                genes.add(first)
    return genes


def fasta_records(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, seq_lines = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_lines)
                header = line.rstrip("\n")
                seq_lines = []
            else:
                seq_lines.append(line.rstrip("\n"))
        if header:
            yield header, "".join(seq_lines)


def main(gene_list, fasta_in, fasta_out):
    genes = load_gene_set(gene_list)
    gn_pattern = re.compile(r"\bGN=([A-Za-z0-9._-]+)")
    written = 0

    with open(fasta_out, "w") as out_fh:
        for header, seq in fasta_records(fasta_in):
            m = gn_pattern.search(header)
            if m and m.group(1) in genes:
                out_fh.write(header + "\n")
                for i in range(0, len(seq), 60):
                    out_fh.write(seq[i : i + 60] + "\n")
                written += 1

    print(f"Written {written} records to {fasta_out}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(main.__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])