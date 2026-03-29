library(memes)
library(universalmotif)
library(Biostrings)
library(dplyr)
library(purrr)
library(stringr)

#' End-to-end 3'UTR motif analysis with memes/MEME Suite
#'
#' @param fg_fa Path to FASTA of target 3'UTRs (foreground).
#' @param bg_fa Path to FASTA of matched background 3'UTRs.
#' @param known_db List of universalmotif objects (e.g., CISBP mouse) for known-motif tests.
#' @param do_denovo Logical; run STREME for de-novo discovery (default TRUE).
#' @param streme_args Named list of STREME args (minw, maxw, nmotifs, seed, alphaset).
#' @param ame_method AME test method ("fisher", "ranksum", "spearman"). Default "fisher".
#' @param ame_qcut Adjusted p-value cutoff to consider motifs "enriched" for scanning. Default 0.05.
#' @param fimo_thresh Per-site p-value threshold for FIMO scanning (default 1e-4).
#' @param use_markov_bg Logical; build a Markov background from bg_fa for FIMO (default TRUE).
#' @param length_normalize Logical; add site density (#sites / kb of UTR) in summaries (default TRUE).
#'
#' @return A list with components:
#' \describe{
#'   \item{streme}{STREME result (or NULL if do_denovo=FALSE).}
#'   \item{ame}{AME enrichment table for known_db.}
#'   \item{tomtom}{TomTom matches of de-novo motifs to known_db (or NULL).}
#'   \item{fimo_hits_known}{FIMO GRanges of sites from enriched known motifs.}
#'   \item{fimo_hits_denovo}{FIMO GRanges of sites from de-novo motifs (if run).}
#'   \item{per_gene_known}{Per-UTR summary for known motifs.}
#'   \item{per_gene_denovo}{Per-UTR summary for de-novo motifs (if run).}
#' }
#'
#' @details
#' - Keep alphabets consistent. If your motifs are RNA (U rows), pass STREME `alphaset="rna"`
#'   and keep FASTAs as DNA or convert T->U consistently for AME/FIMO.
#' - FIMO p-values are more realistic with a Markov background trained on bg_fa.
motif_pipeline_utr <- function(
    fg_fa,
    bg_fa,
    known_db,
    do_denovo   = TRUE,
    streme_args = list(minw=5, maxw=12, nmotifs=10, seed=1, alphaset="rna"),
    ame_method  = "fisher",
    ame_qcut    = 0.05,
    fimo_thresh = 1e-4,
    use_markov_bg   = TRUE,
    length_normalize = TRUE
) {
  # --- deps ---
  req <- c("memes", "universalmotif", "Biostrings", "dplyr", "GenomicRanges")
  miss <- req[!vapply(req, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse=", "))
  
  # --- IO ---
  fg <- Biostrings::readDNAStringSet(fg_fa)
  bg <- Biostrings::readDNAStringSet(bg_fa)
  
  # If your motifs are RNA (U row), AME/FIMO can still use DNA sequences.
  # For fully RNA semantics/logos, uncomment these lines:
  fg_rna <- Biostrings::chartr("T","U", fg)
  bg_rna <- Biostrings::chartr("T","U", bg)
  # (then use fg_rna/bg_rna in AME/FIMO calls)

  out <- list(streme=NULL, ame=NULL, tomtom=NULL,
              fimo_hits_known=NULL, fimo_hits_denovo=NULL,
              per_gene_known=NULL, per_gene_denovo=NULL)
  
  # --- 1) De novo discovery (STREME) ---
  if (isTRUE(do_denovo)) {
    out$streme <- .run_streme(fg, bg, streme_args)
  }
  
  # --- 2) Enrichment of known motifs (AME) ---
  out$ame <- memes::runAme(
    sequences = fg_rna, control = bg_rna,
    database  = known_db,
    method    = ame_method,
    outdir = file.path(getwd(), "meme_runs", "ame")
  )
  ame_sig <- dplyr::arrange(out$ame, .data$adj_pvalue)
  ame_sig <- dplyr::filter(ame_sig, .data$adj_pvalue <= ame_qcut)
  
  # Prepare a name -> motif lookup for the known_db list
  db_tbl <- tibble::tibble(
    name  = vapply(known_db, function(m) m@name, character(1)),
    motif = known_db
  )
  
  # --- 3) Optional: annotate de-novo motifs to known DB (TomTom) ---
  if (!is.null(out$streme) && length(out$streme$motif)) {
    out$tomtom <- memes::runTomTom(query = out$streme$motif, 
                                   database = known_db,
                                   outdir = file.path(getwd(), "meme_runs", "tomtom"))
  }
  
  # --- Background model for FIMO ---
  bgfile <- NULL
  if (isTRUE(use_markov_bg)) {
    bgfile <- tempfile(fileext = ".bfile")
    system2("fasta-get-markov", args = c("-m","0", bg_fa, bgfile))
    if (!file.exists(bgfile)) bgfile <- NULL
  }
  
  # --- 4a) FIMO scan for enriched known motifs ---
  if (nrow(ame_sig) > 0) {
    # join to recover the actual universalmotif objects
    top_known <- dplyr::inner_join(
      ame_sig,
      db_tbl,
      by = c("motif_name" = "name")
    )$motif
    
    out$fimo_hits_known <- memes::runFimo(sequences = fg_rna, 
                                motifs = top_known, 
                                thresh = fimo_thresh,
                                bgfile = bg_rna,
                                outdir = file.path(getwd(), "meme_runs", "fimo"))
    
    # per-UTR summary
    ktbl <- as.data.frame(out$fimo_hits_known)
    out$per_gene_known <- ktbl |>
      dplyr::group_by(.data$sequence_name, .data$motif_altname) |>
      dplyr::summarise(
        n_sites = dplyr::n(),
        min_p   = min(.data$pvalue),
        .groups = "drop"
      )
    
    # optional: normalize by UTR length (sites per kb)
    if (isTRUE(length_normalize)) {
      utr_len <- tibble::tibble(
        sequence_name = names(fg),
        utr_len = width(fg)
      )
      out$per_gene_known <- dplyr::left_join(out$per_gene_known, utr_len, by="sequence_name") |>
        dplyr::mutate(sites_per_kb = 1000 * (.data$n_sites / .data$utr_len))
    }
  }
  
  # --- 4b) FIMO scan for de-novo motifs (optional) ---
  if (!is.null(out$streme) && length(out$streme$motif)) {
    out$fimo_hits_denovo <- memes::runFimo(
      sequences = fg_rna,
      motifs    = out$streme$motif,
      thresh    = fimo_thresh,
      bgfile    = bg_rna
    )
    dtbl <- as.data.frame(out$fimo_hits_denovo)
    out$per_gene_denovo <- dtbl |>
      dplyr::group_by(.data$sequence_name, .data$motif_altname) |>
      dplyr::summarise(
        n_sites = dplyr::n(),
        min_p   = min(.data$pvalue),
        .groups = "drop"
      )
    if (isTRUE(length_normalize)) {
      utr_len <- tibble::tibble(sequence_name = names(fg), utr_len = width(fg))
      out$per_gene_denovo <- dplyr::left_join(out$per_gene_denovo, utr_len, by="sequence_name") |>
        dplyr::mutate(sites_per_kb = 1000 * (.data$n_sites / .data$utr_len))
    }
  }
  
  out
}

#' Load CISBP(-RNA) PWM files robustly
#' (skips zero-byte/empty, prefers read_cisbp, names from filename, sets strand)
#'
#' @param pwm_dir Directory containing CISBP(-RNA) PWM files (*.txt, *.pwm, *.meme)
#' @param pattern Regex for files to read. Defaults to common PWM/MEME extensions.
#' @param alphabet Either "dna" or "rna". If "rna", converts T→U and sets alphabet to RNA.
#' @param motif_strand Motif strand annotation to set on all motifs: "+", "-", or "+-".
#'   Use "+" if your database is forward-only. Default is "+-".
#' @param verbosity "message", "warn", or "quiet" — how to report skipped files.
#'
#' @return A list with:
#'   \item{motifs}{List of \code{universalmotif} objects (may be empty)}
#'   \item{log}{Data frame of skipped files and reasons}
#'
#' @details
#' - Uses `universalmotif::read_cisbp()` when possible, else falls back to other readers.
#' - Skips zero-byte and empty PWMs.
#' - Motif names are set to the file basename (without extension).
#' - Sets the `strand` slot on each motif to `motif_strand`.
#' - Never stops on error; returns empty motif list + log if nothing valid parsed.
load_cisbp_pwm_db <- function(pwm_dir,
                              pattern      = "\\.(txt|pwm|meme|tf|transfac)$",
                              alphabet     = c("dna","rna"),
                              motif_strand = c("+-", "+", "-"),
                              verbosity    = c("message","warn","quiet")) {
  
  alphabet     <- match.arg(alphabet)
  motif_strand <- match.arg(motif_strand)
  verbosity    <- match.arg(verbosity)
  
  say <- function(level = c("message","warn"), ...) {
    level <- match.arg(level)
    if (verbosity == "quiet") return(invisible(NULL))
    txt <- paste0(...)
    if (verbosity == "warn" && level == "warn") warning(txt, call. = FALSE)
    else if (verbosity == "message") message(txt)
  }
  
  files <- list.files(pwm_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) {
    say("warn", "No motif files found in ", pwm_dir)
    return(list(motifs = list(), log = data.frame(file=character(), reason=character())))
  }
  
  fi <- file.info(files)
  zero_byte <- rownames(fi)[fi$size == 0]
  non_empty <- setdiff(files, zero_byte)
  
  skip_log <- list()
  if (length(zero_byte)) {
    skip_log[[length(skip_log)+1]] <- data.frame(
      file   = basename(zero_byte),
      reason = "zero_byte_file",
      stringsAsFactors = FALSE
    )
    say("message", "Skipping ", length(zero_byte), " zero-byte files.")
  }
  
  motifs_out <- list()
  
  read_one <- function(f) {
    readers <- list(
      function(f) universalmotif::read_cisbp(f),
      function(f) universalmotif::read_matrix(f, positions = "rows"),
      function(f) universalmotif::read_transfac(f),
      function(f) universalmotif::read_jaspar(f),
      function(f) universalmotif::read_meme(f)
    )
    mot <- NULL; last_err <- NULL
    for (rd in readers) {
      mot <- tryCatch(rd(f), error = function(e) { last_err <<- e; NULL })
      if (!is.null(mot)) break
    }
    if (is.null(mot)) {
      skip_log[[length(skip_log)+1]] <<- data.frame(
        file = basename(f),
        reason = paste0("unreadable: ", conditionMessage(last_err)),
        stringsAsFactors = FALSE
      )
      return(list())
    }
    if (inherits(mot, "universalmotif")) mot <- list(mot)
    
    keep <- list()
    for (i in seq_along(mot)) {
      mm <- tryCatch(mot[[i]]@motif, error = function(e) NULL)
      if (is.null(mm) || nrow(mm)==0 || ncol(mm)==0 || any(!is.finite(mm))) {
        skip_log[[length(skip_log)+1]] <<- data.frame(
          file = basename(f), reason = "empty_motif_matrix", stringsAsFactors = FALSE
        )
        next
      }
      keep[[length(keep)+1]] <- mot[[i]]
    }
    if (!length(keep)) return(list())
    
    # set name = file basename without extension (append index if multi-motif file)
    base <- sub("\\.[^.]+$", "", basename(f))
    if (length(keep) == 1L) {
      keep[[1]]@name <- base
    } else {
      for (i in seq_along(keep)) keep[[i]]@name <- paste0(base, "_", i)
    }
    
    keep
  }
  
  for (f in non_empty) {
    got <- read_one(f)
    if (length(got)) motifs_out <- c(motifs_out, got)
  }
  
  # set strand and alphabet uniformly
  if (length(motifs_out)) {
    motifs_out <- lapply(motifs_out, function(m) {
      # strand
      m@strand <- motif_strand
      # alphabet conversion (T -> U) if RNA requested
      if (alphabet == "rna") {
        rn <- rownames(m@motif)
        if (!is.null(rn) && "T" %in% rn) {
          rn[rn == "T"] <- "U"
          rownames(m@motif) <- rn
        }
        m@alphabet <- "RNA"
      } else {
        m@alphabet <- "DNA"
      }
      m
    })
  }
  
  log_df <- if (length(skip_log)) do.call(rbind, skip_log) else
    data.frame(file = character(), reason = character(), stringsAsFactors = FALSE)
  
  if (!length(motifs_out)) {
    say("warn", "No valid motifs parsed from ", pwm_dir)
  } else {
    widths <- vapply(motifs_out, function(m) ncol(m@motif), integer(1))
    say("message", "Loaded ", length(motifs_out), " motifs. Width range: ",
        min(widths), "-", max(widths), " nt. Strand set to '", motif_strand, "'.")
  }
  
  list(motifs = motifs_out, log = log_df)
}

# robust STREME caller with safe outdir + arg-name autodetect
.run_streme <- function(fg, bg, streme_args, outdir_root = "meme_runs") {
  fs <- names(formals(memes::runStreme))
  
  # choose an explicit outdir
  od <- file.path(getwd(), outdir_root, "streme_fg_vs_bg")
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  streme_args$outdir <- od
  
  # map arg names across memes versions
  if ("alph" %in% fs && "alphaset" %in% names(streme_args)) {
    streme_args$alph <- streme_args$alphaset
    streme_args$alphaset <- NULL
  }
  
  # default to *dna* unless you converted sequences to U
  if (!("alph" %in% names(streme_args)) && !("alphaset" %in% names(streme_args))) {
    streme_args$alph <- "dna"
  }
  
  # make sure errors are printed
  if (!("silent" %in% names(streme_args))) streme_args$silent <- FALSE
  
  # also pass meme_path if you set options(memes.meme_path)
  if (is.null(streme_args$meme_path)) {
    mp <- getOption("memes.meme_path", default = NULL)
    if (!is.null(mp)) streme_args$meme_path <- mp
  }
  
  # try/catch so failure doesn’t stop pipeline
  out <- try({
    if ("input" %in% fs) {
      do.call(memes::runStreme, c(list(input = fg, control = bg), streme_args))
    } else if ("foreground" %in% fs) {
      do.call(memes::runStreme, c(list(foreground = fg, background = bg), streme_args))
    } else if ("sequences" %in% fs) {
      do.call(memes::runStreme, c(list(sequences = fg, control = bg), streme_args))
    } else {
      stop("Unknown runStreme() signature: ", paste(fs, collapse = ", "))
    }
  }, silent = FALSE)
  
  if (inherits(out, "try-error")) {
    message("\nSTREME failed. Check log files under: ", od,
            "\n(look for files like 'streme.log' or console output above).",
            "\nContinuing pipeline without de-novo motifs.")
    return(NULL)
  }
  out
}

# ---------- workflow ----------
setwd("~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment")

# Paths to your data
fg_fa <- "~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment/1-extract-seq/actin_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"     # your 3' UTR targets (DNA alphabet is fine)
bg_fa <- "~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment/1-extract-seq/actin_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"  # matched background UTRs


# 1) Load your CISBP-RNA PWMs (mouse). Keep "dna" alphabet to match DNA FASTAs.
cisbp_dir <- "~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment/CISBP-RNA_Mus_Musculus/pwms_all_motifs"
cisbp_res <- load_cisbp_pwm_db(cisbp_dir, 
                               alphabet = "rna",
                               motif_strand = "+")
cisbp_motifs <- cisbp_res$motifs


# 2) (Optional but recommended) bundle into a single MEME DB for speed/CLI compatibility
write_meme(cisbp_motifs, 
           file = "cisbp_mus_muculus_motifs.meme", 
           version = 5, 
           strand = "+")

# pipeline

# Use a short, safe working dir for MEME outputs (no spaces/special chars)
base_out <- file.path(getwd(), "meme_runs"); dir.create(base_out, showWarnings = FALSE)

# 3) Keep alphabet consistent: DNA sequences => use alph="dna".
#    (If you truly want RNA, convert FASTAs to U and use alph="rna".)

res <- motif_pipeline_utr(fg_fa = fg_fa,
                          bg_fa = bg_fa,
                          known_db = cisbp_motifs,
                          do_denovo   = TRUE,
                          streme_args = list(minw=5, maxw=12, nmotifs=10, seed=1, alphaset="rna"),
                          ame_method  = "fisher",
                          ame_qcut    = 0.05,
                          fimo_thresh = 1e-4,
                          use_markov_bg   = TRUE,
                          length_normalize = TRUE)


# 2) Use a short, safe working dir for MEME outputs (no spaces/special chars)
base_out <- file.path(getwd(), "meme_runs"); dir.create(base_out, showWarnings = FALSE)


# 3) Enrichment of known motifs in targets vs background (AME)
#    You can pass either a universalmotif list (cisbp_motifs) OR the MEME file.
#    Using the list keeps everything in-R; the MEME file uses the CLI under the hood.
ame_res <- runAme(
  sequences = readDNAStringSet(fg_fa),
  control   = readDNAStringSet(bg_fa),
  database  = cisbp_motifs,        # or database = meme_db
  method    = "fisher"              # "fisher", "ranksum", "spearman"
)
dplyr::arrange(ame_res, adj_pvalue) %>% dplyr::slice_head(n = 25)

# 4) Annotate your de novo motifs (from STREME/DREME/MEME) to CISBP with TomTom
#    Suppose you already ran STREME and have `streme_res$motif`:
# streme_res <- runStreme(readDNAStringSet(fg_fa), readDNAStringSet(bg_fa), minw=5, maxw=12, nmotifs=10)
# tom <- runTomTom(query = streme_res$motif, database = cisbp_motifs)
# head(tom$best_match)

# 5) Scan your UTRs for instances of top CISBP motifs with FIMO
#    Build a Markov background from your background FASTA for realistic p-values.
bg_model <- "bg_markov.bfile"
if (!file.exists(bg_model)) {
  system2("fasta-get-markov", c("-m","0", bg_fa, bg_model))
}

# choose a subset of enriched motifs from AME (e.g., adj_pvalue < 0.05)
top_motifs <- ame_res %>%
  filter(adj_pvalue < 0.05) %>%
  inner_join(
    tibble(name = map_chr(cisbp_motifs, ~ .x@name), motif = cisbp_motifs),
    by = c("motif_name" = "name")
  ) %>%
  pull(motif)

# fallback if nothing passes (keep a handful of strongest)
if (length(top_motifs) == 0) top_motifs <- cisbp_motifs[seq_len(min(50, length(cisbp_motifs)))]

fimo_hits <- runFimo(
  sequences = readDNAStringSet(fg_fa),
  motifs    = top_motifs,
  thresh    = 1e-4,
  bgfile    = bg_model
)

# Summarize hits per UTR and motif
fimo_tbl <- as.data.frame(fimo_hits) %>%
  group_by(sequence_name, motif_altname) %>%
  summarise(n_sites = n(), min_p = min(pvalue), .groups = "drop") %>%
  arrange(min_p)

# ---------- optional: RNA alphabet end-to-end ----------
# If you want everything in RNA-space (logos with U’s, etc.):
#   1) load motifs with alphabet="rna"
#   2) convert sequences to RNA (T->U)
#   3) use alphaset="rna" in discovery steps
#
# fg_rna <- chartr("T","U", readDNAStringSet(fg_fa))
# bg_rna <- chartr("T","U", readDNAStringSet(bg_fa))
# ame_res <- runAme(sequences = fg_rna, control = bg_rna, database = cisbp_motifs, method = "fisher")
# fimo_hits <- runFimo(sequences = fg_rna, motifs = top_motifs, thresh = 1e-4)