#!/usr/bin/env bash
set -euo pipefail

#############################################################
# Run TomTom (MEME Suite via Docker) vs a ready MEME database
# - Inputs:
#     - One or more MEME-format query files (STREME/MEME/DREME outputs)
#     - A single MEME-format database file (cisbp_rna_db.meme)
# - Output:
#     - tomtom_<query_basename>/ directories (HTML + TSV)
#############################################################

# --- USER INPUTS (edit these) ---
# Put absolute paths here; don't add backslashes for spaces — just quote the string.
QUERY_FILES=(
  "motifs_deduplicated.meme"
  # add more query files if you like
)

DB_MEME="oRNAment_arc_rbp.meme"

# Where to write TomTom outputs (host path). Defaults to current dir.
OUT_ROOT="$(pwd)"

# Docker image
DOCKER_IMAGE="memesuite/memesuite:latest"

# TomTom settings
TOMTOM_DIST="ed"     # alternatives: sandelin, pearson, kullback
TOMTOM_THRESH="0.05" # q-value threshold
TOMTOM_NORC=1        # 1 = disallow reverse-complement matches (typical for oriented 3'UTRs)

# --- Sanity checks ---
[ -f "$DB_MEME" ] || { echo "ERROR: DB_MEME not found: $DB_MEME" >&2; exit 1; }
mkdir -p "$OUT_ROOT"

# Convenience flag
NORC_FLAG=()
[ "$TOMTOM_NORC" -eq 1 ] && NORC_FLAG=(--norc)

# --- Run TomTom per query file ---
for q in "${QUERY_FILES[@]}"; do
  if [ ! -f "$q" ]; then
    echo "WARNING: Query file not found, skipping: $q" >&2
    continue
  fi

  q_dir="$(dirname "$q")"
  q_base="$(basename "$q")"
  q_stem="${q_base%.*}"

  db_dir="$(dirname "$DB_MEME")"
  db_base="$(basename "$DB_MEME")"

  out_dir_host="${OUT_ROOT}/tomtom_${q_stem}"
  mkdir -p "$out_dir_host"

  echo ">> Running TomTom:"
  echo "   query: $q"
  echo "   db   : $DB_MEME"
  echo "   out  : $out_dir_host"

  # Mount the query parent as /q, the DB parent as /db, and the OUT_ROOT as /out.
  docker run --rm \
    -u "$(id -u)":"$(id -g)" \
    -v "${q_dir}:/q" \
    -v "${db_dir}:/db" \
    -v "${out_dir_host}:/out" \
    "$DOCKER_IMAGE" \
    tomtom \
      -no-ssc \
      -dist "$TOMTOM_DIST" \
      -thresh "$TOMTOM_THRESH" \
      "${NORC_FLAG[@]}" \
      -oc "/out" \
      "/q/${q_base}" \
      "/db/${db_base}"

  echo "   Open: ${out_dir_host}/tomtom.html"
done

echo ">> All TomTom runs completed under: $OUT_ROOT"