#!/usr/bin/env bash
###############################################
# DREME de novo k-mer motif discovery (MEME Suite)
# - Input: target and background FASTA files
# - Output: HTML + text + XML + logos per run
# - Runtime: Docker container (memesuite)
###############################################

# --- USER INPUTS (edit these paths) ---
TARGET_FASTA="actin_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"
BACKGROUND_FASTA="actin_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"

# Optional: a short label for this dataset (used in output folder names)
DATASET_TAG="actin_3pUTR"
# If you want DREME to search ONLY the given strand (not reverse complement),
# set this to 1 (recommended for oriented 3' UTRs). Otherwise leave 0.
DREME_NORC=1

# --- DOCKER IMAGE & MOUNT SETTINGS ---
DOCKER_IMAGE="memesuite/memesuite:latest"
WORKDIR="$(pwd)"

# --- QUICK SANITY CHECKS ---
[ -f "$TARGET_FASTA" ]     || { echo "ERROR: TARGET_FASTA not found: $TARGET_FASTA" >&2; exit 1; }
[ -f "$BACKGROUND_FASTA" ] || { echo "ERROR: BACKGROUND_FASTA not found: $BACKGROUND_FASTA" >&2; exit 1; }

# --- CREATE A RESULTS ROOT DIRECTORY ---
RESULTS_DIR="${WORKDIR}/dreme_results_${DATASET_TAG}"
mkdir -p "$RESULTS_DIR"

# --- RECORD THE DREME VERSION (useful for methods) ---
docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" dreme -version | tee "${RESULTS_DIR}/dreme_version.txt"

# Helper: NORC flag (restrict search to given strand for complementable alphabets)
NORC_FLAG=()
[ "$DREME_NORC" -eq 1 ] && NORC_FLAG=(-norc)

# ------------------------------------------------------------
# RUN 1: default short cores (k=3–7), relaxed E-value
# Good first sweep for ultra-short RBP words and U-rich cores
# ------------------------------------------------------------
OUT1="${RESULTS_DIR}/dreme_k3-7_e0.05_seed101"
mkdir -p "$OUT1"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  dreme \
  -dna \
  -p "/data/${TARGET_FASTA##*/}" \
  -n "/data/${BACKGROUND_FASTA##*/}" \
  -mink 3 \
  -maxk 7 \
  -e 0.05 \
  -m 20 \
  -s 101 \
  "${NORC_FLAG[@]}" \
  -png \
  -oc "/data/${OUT1##*/}" \
  -rna

# ------------------------------------------------------------
# RUN 2: slightly wider cores (k=4–8), balanced threshold
# Captures flanked words and a bit more specificity
# ------------------------------------------------------------
OUT2="${RESULTS_DIR}/dreme_k4-8_e0.02_seed202"
mkdir -p "$OUT2"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  dreme \
  -dna \
  -p "/data/${TARGET_FASTA##*/}" \
  -n "/data/${BACKGROUND_FASTA##*/}" \
  -mink 4 \
  -maxk 8 \
  -e 0.02 \
  -m 15 \
  -s 202 \
  "${NORC_FLAG[@]}" \
  -png \
  -oc "/data/${OUT2##*/}" \
  -rna

# ------------------------------------------------------------
# RUN 3: stringent AU-bias guardrail (k=5–7), strict E-value
# Helps reduce low-complexity mirages in AU-rich UTRs
# ------------------------------------------------------------
OUT3="${RESULTS_DIR}/dreme_k5-7_e0.005_seed303"
mkdir -p "$OUT3"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  dreme \
  -dna \
  -p "/data/${TARGET_FASTA##*/}" \
  -n "/data/${BACKGROUND_FASTA##*/}" \
  -mink 5 \
  -maxk 7 \
  -e 0.005 \
  -m 10 \
  -s 303 \
  "${NORC_FLAG[@]}" \
  -png \
  -oc "/data/${OUT3##*/}" \
  -rna

# --- WHAT YOU GET IN EACH OUT DIRECTORY ---
# - dreme.html : interactive report (open in a browser)
# - dreme.txt  : text summary
# - dreme.xml  : machine-readable results
# - logos/     : motif sequence logos (PNG)

echo ""
echo "Done. Results written under: $RESULTS_DIR"
echo "Open the HTML reports, e.g.:"
echo "  ${OUT1}/dreme.html"
echo "  ${OUT2}/dreme.html"
echo "  ${OUT3}/dreme.html"