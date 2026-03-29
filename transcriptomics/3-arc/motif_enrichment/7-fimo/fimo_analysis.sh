#!/usr/bin/env bash
set -euo pipefail

#############################################################
# FIMO motif scanning (MEME Suite via Docker)
# - Inputs:
#     - MEME-format motif files (e.g., streme.meme, dreme.txt, meme.txt)
#     - FASTA sequences to scan (e.g., target UTRs, background UTRs)
#     - (optional) background model file from fasta-get-markov
# - Outputs:
#     - fimo_<motif>__on__<seq>/ (TSV+HTML per run)
#############################################################

# --------- USER INPUTS (edit these) ---------

# Motif files (absolute paths are fine; quote strings, do not use backslashes for spaces)
MOTIF_FILES=(
  "motifs_deduplicated.meme"
  "cisbp_rna_db.meme"
  "oRNAment_arc_rbp.meme"
)

# Sequences to scan (you can include both target and background FASTAs)
SEQ_FASTAS=(
  #"arc_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"
  "arc_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"
)

# Optional: background letter-frequency model (from fasta-get-markov)
# Leave empty to use the background embedded in motif files.
BFILE=""   # or "" to skip

# Where to write outputs (host path)
OUT_ROOT="$(pwd)"

# Docker image
DOCKER_IMAGE="memesuite/memesuite:latest"

# FIMO settings
FIMO_QTHRESH="0.05"   # threshold on q-values (FDR) when --qv-thresh is set
FIMO_NORC=1           # 1 = scan only the given strand (typical for 3'UTRs)
FIMO_MAX_STORED=500000  # raise if sequences are large; controls memory/IO

# -------------------------------------------

mkdir -p "$OUT_ROOT"

# Convenience flags
NORC_FLAG=()
[ "$FIMO_NORC" -eq 1 ] && NORC_FLAG=(--norc)

BG_ARGS=()
if [ -n "${BFILE}" ]; then
  [ -f "$BFILE" ] || { echo "ERROR: BFILE not found: $BFILE" >&2; exit 1; }
fi

# Iterate motif files × sequence FASTAs
for motif in "${MOTIF_FILES[@]}"; do
  if [ ! -f "$motif" ]; then
    echo "WARNING: Motif file not found, skipping: $motif" >&2
    continue
  fi
  motif_dir="$(dirname "$motif")"
  motif_base="$(basename "$motif")"
  motif_stem="${motif_base%.*}"

  for seqfa in "${SEQ_FASTAS[@]}"; do
    if [ ! -f "$seqfa" ]; then
      echo "WARNING: FASTA not found, skipping: $seqfa" >&2
      continue
    fi
    seq_dir="$(dirname "$seqfa")"
    seq_base="$(basename "$seqfa")"
    seq_stem="${seq_base%.*}"

    out_dir_host="${OUT_ROOT}/fimo_${motif_stem}__on__${seq_stem}"
    mkdir -p "$out_dir_host"

    echo ">> FIMO: ${motif_base} on ${seq_base}"
    echo "   out : ${out_dir_host}"

    # Build mounts: query motifs (/motifs), sequences (/seqs), output (/out), and optional bgfile (/bg)
    if [ -n "${BFILE}" ]; then
      bg_dir="$(dirname "$BFILE")"
      bg_base="$(basename "$BFILE")"
      docker run --rm \
        -u "$(id -u)":"$(id -g)" \
        -v "${motif_dir}:/motifs" \
        -v "${seq_dir}:/seqs" \
        -v "${bg_dir}:/bg" \
        -v "${out_dir_host}:/out" \
        "$DOCKER_IMAGE" \
        fimo \
          --verbosity 1 \
          --oc /out \
          --max-stored-scores "${FIMO_MAX_STORED}" \
          "${NORC_FLAG[@]}" \
          --bgfile "/bg/${bg_base}" \
          "/motifs/${motif_base}" \
          "/seqs/${seq_base}"
    else
      docker run --rm \
        -u "$(id -u)":"$(id -g)" \
        -v "${motif_dir}:/motifs" \
        -v "${seq_dir}:/seqs" \
        -v "${out_dir_host}:/out" \
        "$DOCKER_IMAGE" \
        fimo \
          --verbosity 1 \
          --oc /out \
          --max-stored-scores "${FIMO_MAX_STORED}" \
          "${NORC_FLAG[@]}" \
          "/motifs/${motif_base}" \
          "/seqs/${seq_base}"
    fi

    echo "   wrote: ${out_dir_host}/fimo.tsv  (plus fimo.gff, fimo.html, etc.)"
  done
done

echo ">> All FIMO scans finished. Results under: ${OUT_ROOT}"