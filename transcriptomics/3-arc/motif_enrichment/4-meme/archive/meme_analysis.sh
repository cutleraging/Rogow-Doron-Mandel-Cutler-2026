#!/usr/bin/env bash
#############################################################
# MEME de novo motif discovery (MEME Suite via Docker)
# - Inputs : target and background FASTA files
# - Outputs: HTML + text (MEME format) + XML per run
# - Notes  : Builds a background Markov model from background
#############################################################

# --- USER INPUTS (edit these paths) ---
TARGET_FASTA="arc_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"
BACKGROUND_FASTA="arc_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"

# Short label to keep results tidy
DATASET_TAG="arc_3pUTR"

# Respect UTR orientation (search only the given strand)?
# 1 = search only given strand (recommended for 3' UTRs)
# 0 = allow reverse complement
MEME_NOREVCOMP=1

# --- DOCKER IMAGE & MOUNT SETTINGS ---
DOCKER_IMAGE="memesuite/memesuite:latest"
WORKDIR="$(pwd)"

# --- QUICK SANITY CHECKS ---
[ -f "$TARGET_FASTA" ]     || { echo "ERROR: TARGET_FASTA not found: $TARGET_FASTA" >&2; exit 1; }
[ -f "$BACKGROUND_FASTA" ] || { echo "ERROR: BACKGROUND_FASTA not found: $BACKGROUND_FASTA" >&2; exit 1; }

# --- RESULTS ROOT DIRECTORY ---
RESULTS_DIR="${WORKDIR}"
mkdir -p "$RESULTS_DIR"

# --- RECORD VERSIONS (useful for Methods) ---
docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" meme -version        | tee "${RESULTS_DIR}/meme_version.txt"
docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" fasta-get-markov -h  >/dev/null 2>&1 && echo "fasta-get-markov available" >> "${RESULTS_DIR}/meme_version.txt"

# --- BUILD BACKGROUND MARKOV MODEL FROM BACKGROUND FASTA ---
# Use order-1 and order-2 models (UTRs are AU-rich; higher order helps curb composition artifacts)
BFILE_DIR="${RESULTS_DIR}/background_models"
mkdir -p "$BFILE_DIR"

# Order-1 model
docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  fasta-get-markov -dna -m 1 "/data/${BACKGROUND_FASTA##*/}" "/data/${BFILE_DIR##*/}/bfile_order1.meme"

# Order-2 model
docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  fasta-get-markov -dna -m 2 "/data/${BACKGROUND_FASTA##*/}" "/data/${BFILE_DIR##*/}/bfile_order2.meme"

# Choose which background file to use for the runs below (order-2 is stricter)
BFILE_CHOICE="${BFILE_DIR}/bfile_order2.meme"

# Helper: strand handling
REVCOMP_FLAG=()
[ "$MEME_NOREVCOMP" -eq 1 ] && REVCOMP_FLAG=(-norevcomp) || REVCOMP_FLAG=(-revcomp)

#############################################################
# RUN 1: ZOOPS model, widths 6–12 (common RBP-like footprints)
# ZOOPS = Zero or One site Per Sequence (good first assumption)
#############################################################
OUT1="${RESULTS_DIR}/meme_zoops_w6-12_seed11"
mkdir -p "$OUT1"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  meme "/data/${TARGET_FASTA##*/}" \
  -dna \
  -mod zoops \
  -nmotifs 10 \
  -minw 6 -maxw 12 \
  -bfile "/data/${BFILE_CHOICE##*/}" \
  -evt 0.05 \                # stop when next motif E-value > 0.05
  -minsites 10 \             # require some recurrence across sequences
  -maxsites 300 \            # avoid runaway repeats
  -seed 11 \
  "${REVCOMP_FLAG[@]}" \
  -oc "/data/${OUT1##*/}" \
  -rna

#############################################################
# RUN 2: ANR model, widths 6–12 (for repeated sites per UTR)
# ANR  = Any Number of Repetitions (captures multi-site RBPs)
#############################################################
OUT2="${RESULTS_DIR}/meme_anr_w6-12_seed22"
mkdir -p "$OUT2"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  meme "/data/${TARGET_FASTA##*/}" \
  -dna \
  -mod anr \
  -nmotifs 8 \
  -minw 6 -maxw 12 \
  -bfile "/data/${BFILE_CHOICE##*/}" \
  -evt 0.02 \
  -seed 22 \
  "${REVCOMP_FLAG[@]}" \
  -oc "/data/${OUT2##*/}" \
  -rna

#############################################################
# RUN 3: ZOOPS, narrower 5–8, stricter stopping (AU-bias guardrail)
# Helps isolate compact cores while resisting low-complexity noise
#############################################################
OUT3="${RESULTS_DIR}/meme_zoops_w5-8_evt0.01_seed33"
mkdir -p "$OUT3"

docker run --rm -v "${WORKDIR}:/data" "$DOCKER_IMAGE" \
  meme "/data/${TARGET_FASTA##*/}" \
  -dna \
  -mod zoops \
  -nmotifs 6 \
  -minw 5 -maxw 8 \
  -bfile "/data/${BFILE_CHOICE##*/}" \
  -evt 0.01 \
  -minsites 8 \
  -seed 33 \
  "${REVCOMP_FLAG[@]}" \
  -oc "/data/${OUT3##*/}" \
  -rna

# --- WHAT YOU GET IN EACH OUT DIRECTORY ---
# - meme.html : interactive run report (open in a browser)
# - meme.txt  : motifs in MEME format (feed this to TomTom/FIMO/AME)
# - meme.xml  : machine-readable results
# - logos/    : motif sequence logos
# - seqs.fa   : (in some versions) filtered/processed sequences

echo ""
echo "Done. Results written under: $RESULTS_DIR"
echo "Open the HTML reports, e.g.:"
echo "  ${OUT1}/meme.html"
echo "  ${OUT2}/meme.html"
echo "  ${OUT3}/meme.html"