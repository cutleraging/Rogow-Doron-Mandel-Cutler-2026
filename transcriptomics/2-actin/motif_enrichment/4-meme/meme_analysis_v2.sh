#!/usr/bin/env bash
set -euo pipefail

#############################################################
# MEME de novo motif discovery (Docker)
#############################################################

# --- USER INPUTS ---
TARGET_FASTA="actin_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"
BACKGROUND_FASTA="actin_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"
DOCKER_IMAGE="memesuite/memesuite:latest"

# --- PATHS (HOST) ---
WORKDIR="$(pwd)"
RESULTS_DIR="${WORKDIR}"   # write directly into current working dir
BFILE_DIR="${RESULTS_DIR}/background_models"

# Ensure required host dirs exist so the container can write into them
mkdir -p "$BFILE_DIR"

# --- CONTAINER MOUNTPOINT & MIRROR PATHS ---
DATA_IN="/data"                                   # inside-container root for WORKDIR
BFILE_DIR_IN="${DATA_IN}/background_models"       # mirrors ${BFILE_DIR} inside container

# --- RECORD VERSIONS ---
docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" meme -version \
  | tee "${RESULTS_DIR}/meme_version.txt"
docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" fasta-get-markov -h \
  >/dev/null 2>&1 && echo "fasta-get-markov available" >> "${RESULTS_DIR}/meme_version.txt"

# --- BUILD BACKGROUND MARKOV MODELS ---
# Write directly to host files via stdout redirection (simplest & robust)
docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" \
  fasta-get-markov -rna -m 1 "${DATA_IN}/${BACKGROUND_FASTA}" > "${BFILE_DIR}/bfile_order1.meme"

docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" \
  fasta-get-markov -rna -m 2 "${DATA_IN}/${BACKGROUND_FASTA}" > "${BFILE_DIR}/bfile_order2.meme"

# Choose which bfile to use
BFILE_CHOICE_IN="${BFILE_DIR_IN}/bfile_order2.meme"

# --- RUN 1: ZOOPS 6–12 ---
OUT1_HOST="${RESULTS_DIR}/meme_zoops_w6-12_seed11"
OUT1_IN="${DATA_IN}/meme_zoops_w6-12_seed11"
mkdir -p "$OUT1_HOST"

docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" \
  meme "${DATA_IN}/${TARGET_FASTA}" \
  -rna \
  -mod zoops \
  -nmotifs 10 \
  -minw 6 -maxw 12 \
  -bfile "${BFILE_CHOICE_IN}" \
  -evt 0.05 \
  -minsites 10 \
  -maxsites 300 \
  -seed 11 \
  -oc "${OUT1_IN}"

# --- RUN 2: ANR 6–12 ---
OUT2_HOST="${RESULTS_DIR}/meme_anr_w6-12_seed22"
OUT2_IN="${DATA_IN}/meme_anr_w6-12_seed22"
mkdir -p "$OUT2_HOST"

docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" \
  meme "${DATA_IN}/${TARGET_FASTA}" \
  -rna \
  -mod anr \
  -nmotifs 8 \
  -minw 6 -maxw 12 \
  -bfile "${BFILE_CHOICE_IN}" \
  -evt 0.02 \
  -seed 22 \
  -oc "${OUT2_IN}"

# --- RUN 3: ZOOPS 5–8 (stricter) ---
OUT3_HOST="${RESULTS_DIR}/meme_zoops_w5-8_evt0.01_seed33"
OUT3_IN="${DATA_IN}/meme_zoops_w5-8_evt0.01_seed33"
mkdir -p "$OUT3_HOST"

docker run --rm -u "$(id -u)":"$(id -g)" -v "${WORKDIR}:${DATA_IN}" "$DOCKER_IMAGE" \
  meme "${DATA_IN}/${TARGET_FASTA}" \
  -rna \
  -mod zoops \
  -nmotifs 6 \
  -minw 5 -maxw 8 \
  -bfile "${BFILE_CHOICE_IN}" \
  -evt 0.01 \
  -minsites 8 \
  -seed 33 \
  -oc "${OUT3_IN}"

echo
echo "Done. Open the reports:"
echo "  ${OUT1_HOST}/meme.html"
echo "  ${OUT2_HOST}/meme.html"
echo "  ${OUT3_HOST}/meme.html"