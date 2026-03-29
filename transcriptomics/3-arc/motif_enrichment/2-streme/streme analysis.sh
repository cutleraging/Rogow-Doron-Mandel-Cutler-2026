###############################################
# STREME de novo motif discovery (MEME Suite)
# - Input: target and background FASTA files
# - Output: HTML + MEME-format motifs per run
# - Runtime: Docker container (memesuite)
###############################################

# --- USER INPUTS (edit these paths) ---
TARGET_FASTA="arc_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"
BACKGROUND_FASTA="arc_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"

# Optional: a short label for this dataset (used in output folder names)
DATASET_TAG="arc_3pUTR"

# --- DOCKER IMAGE & MOUNT SETTINGS ---
# Use the official MEME Suite image tag you pulled (e.g., latest or a specific version)
DOCKER_IMAGE="memesuite/memesuite:latest"

# Use a working directory that contains your FASTAs (will be mounted to /data in the container)
WORKDIR="$(pwd)"

# --- QUICK SANITY CHECKS ---
[ -f "$TARGET_FASTA" ]     || { echo "ERROR: TARGET_FASTA not found: $TARGET_FASTA" >&2; exit 1; }
[ -f "$BACKGROUND_FASTA" ] || { echo "ERROR: BACKGROUND_FASTA not found: $BACKGROUND_FASTA" >&2; exit 1; }

# --- CREATE A RESULTS ROOT DIRECTORY ---
RESULTS_DIR="${WORKDIR}/streme_results_${DATASET_TAG}"
mkdir -p "$RESULTS_DIR"

# --- RECORD THE MEME SUITE VERSION (useful for methods) ---
docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" streme --version | tee "${RESULTS_DIR}/streme_version.txt"

# ------------------------------------------------------------
# RUN 1: short motifs (4–8 nt), typical for many RBP footprints
# ------------------------------------------------------------
OUT1="${RESULTS_DIR}/streme_w4-8_seed42"
mkdir -p "$OUT1"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  streme \
  --dna \
  --p "/data/${TARGET_FASTA##*/}" \
  --n "/data/${BACKGROUND_FASTA##*/}" \
  --minw 4 \
  --maxw 8 \
  --nmotifs 10 \
  --thresh 0.05 \
  --seed 42 \
  --oc "/data/${OUT1##*/}" \
  --rna
# ------------------------------------------------------------
# RUN 2: slightly wider motifs (6–10 nt) to catch bipartite/longer patterns
# ------------------------------------------------------------
OUT2="${RESULTS_DIR}/streme_w6-10_seed43"
mkdir -p "$OUT2"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  streme \
  --dna \
  --p "/data/${TARGET_FASTA##*/}" \
  --n "/data/${BACKGROUND_FASTA##*/}" \
  --minw 6 \
  --maxw 10 \
  --nmotifs 10 \
  --thresh 0.05 \
  --seed 43 \
  --oc "/data/${OUT2##*/}" \
  --rna

# ------------------------------------------------------------
# RUN 3 (optional): AU-rich guardrails
# Use a stricter threshold and fewer motifs to avoid low-complexity mirages.
# ------------------------------------------------------------
OUT3="${RESULTS_DIR}/streme_w5-7_stringent_seed44"
mkdir -p "$OUT3"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  streme \
  --dna \
  --p "/data/${TARGET_FASTA##*/}" \
  --n "/data/${BACKGROUND_FASTA##*/}" \
  --minw 5 \
  --maxw 7 \
  --nmotifs 5 \
  --thresh 0.01 \
  --seed 44 \
  --oc "/data/${OUT3##*/}" \
  --rna

# --- WHAT YOU GET IN EACH OUT DIRECTORY ---
# - streme.html  : interactive report (open in browser)
# - streme.txt   : text summary (readable)
# - streme.meme  : motifs in MEME format (feed into TomTom/FIMO/AME later)
# - logos/       : motif sequence logos

echo ""
echo "Done. Results written under: $RESULTS_DIR"
echo "Open the HTML reports, e.g.:"
echo "  ${OUT1}/streme.html"
echo "  ${OUT2}/streme.html"
echo "  ${OUT3}/streme.html"