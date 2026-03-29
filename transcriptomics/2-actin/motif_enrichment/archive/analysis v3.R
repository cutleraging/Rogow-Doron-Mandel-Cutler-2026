# ---- functions ----

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

# ---- 0) Packages ----------------------------------------------------------------
# (Install from Bioc if needed: BiocManager::install(c("memes","universalmotif","Biostrings","dplyr","ggplot2")))
suppressPackageStartupMessages({
  library(memes)
  library(universalmotif)
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ---- 1) Inputs  ---------------------------------------------
target_fasta <- "~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment/1-extract-seq/actin_vs_no_infect_scramble_sig_enrich_3pSeqs.fasta"     # your 3' UTR targets (DNA alphabet is fine)
background_fasta <- "~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment/1-extract-seq/actin_vs_no_infect_scramble_nonSig_negEnrich_3pSeqs.fasta"  # matched background UTRs

cisbp_dir <- "~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-actin/motif_enrichment/CISBP-RNA_Mus_Musculus/pwms_all_motifs"
cisbp_res <- load_cisbp_pwm_db(cisbp_dir, 
                               alphabet = "rna",
                               motif_strand = "+")
rna_motifs_um <- cisbp_res$motifs

# ---- 2) Basic checks ------------------------------------------------------------
if (!memes::meme_is_installed()) stop("MEME Suite not found. Install MEME and/or set MEME_PATH. See check_meme_install().")
if (!exists("rna_motifs_um")) stop("Provide your CISBP-RNA motifs as a universalmotif list in object `rna_motifs_um`.")

# ---- 3) Read sequences ----------------------------------------------------------
names(targets) <- gsub("[^A-Za-z0-9._-]", "_", names(targets))
names(bg)      <- gsub("[^A-Za-z0-9._-]", "_", names(bg))

targets <- Biostrings::DNAStringSet(toupper(as.character(targets)))
bg      <- Biostrings::DNAStringSet(toupper(as.character(bg)))

# drop any empty or absurdly short sequences (shorter than 4 nt)
targets <- targets[width(targets) >= 4]
bg      <- bg[width(bg) >= 4]

# remove sequences with non-ACGT (N, R, etc.) to avoid AME choking; or mask them to N
bad_targets <- Biostrings::alphabetFrequency(targets, baseOnly = TRUE)[, c("A","C","G","T")] |> rowSums() < width(targets)
bad_bg      <- Biostrings::alphabetFrequency(bg,      baseOnly = TRUE)[, c("A","C","G","T")] |> rowSums() < width(bg)
targets <- targets[!bad_targets]
bg      <- bg[!bad_bg]

length(targets); length(bg)
stopifnot(length(targets) > 1L, length(bg) > 1L)

# ensure sequence names exist (AME/FIMO like named entries)
if (any(is.na(names(targets)) | names(targets) == "")) names(targets) <- paste0("target_", seq_along(targets))
if (any(is.na(names(bg))      | names(bg)      == "")) names(bg)      <- paste0("bg_",     seq_along(bg))

# ---- 4) Harmonize motif alphabet (RNA->DNA; U becomes T) -----------------------
# CISBP-RNA motifs are often RNA; convert to DNA alphabet for DNAStringSet inputs
motifs_dna <- universalmotif::switch_alph(rna_motifs_um)  # toggles RNA<->DNA; safe if already DNA
# optional: drop duplicated PWMs to speed things up
motifs_dna <- memes::remove_duplicate_motifs(motifs_dna)

# ---- 5) Motif enrichment with AME (targets vs background) ----------------------
# AME modes: "discriminative" here (targets vs control sequences). Default method = "fisher".
# You can tweak k-mer shuffle, scoring, thresholds via ... args (see ?runAme and vignette).
ame_res <- memes::runAme(
  input    = targets,
  control  = bg,               # you provided a real background, great
  database = motifs_dna,
  method   = "fisher",
  scoring  = "avg",
  evalue_report_threshold = 10,
  silent   = TRUE
)
# save full AME table
readr::write_csv(ame_res, "AME_targets_vs_background.csv")

# quick glance: top enriched (FDR)
ame_top <- ame_res %>%
  mutate(rank = dense_rank(adj.pvalue)) %>%
  arrange(adj.pvalue) %>%
  select(rank, motif_alt_id, name = motif_id, length, positive_hits, negative_hits, pos_frac, neg_frac, pvalue, adj.pvalue, eval)
ame_top %>% head(20)

# ---- 6) Visualize AME results as a heatmap -------------------------------------
# (Nice overview; for many motifs, consider capping scale_max)
p_heat <- ame_res %>%
  memes::plot_ame_heatmap(                      # by default, value = -log10(adj.pvalue)
    id        = motif_alt_id,
    scale_max = 7.5
  ) + ggtitle("Motif enrichment in 3'UTRs (AME: targets vs background)")
ggplot2::ggsave("AME_enrichment_heatmap.pdf", p_heat, width = 6, height = 7)

# ---- 7) Pick motifs to scan at site-level (FIMO) --------------------------------
# use FDR threshold or top-N; adjust threshold as needed
motifs_for_fimo <- ame_res %>%
  filter(adj.pvalue <= 0.05) %>%
  arrange(adj.pvalue) %>%
  pull(motif)                        # this column holds universalmotif objects
# if empty under strict FDR, relax using top-N:
if (length(motifs_for_fimo) == 0) {
  motifs_for_fimo <- ame_res %>% arrange(pvalue) %>% slice_head(n = 20) %>% pull(motif)
}

# ---- 8) Scan targets with FIMO to get binding site coordinates -----------------
# FIMO returns a GRanges with motif_id, motif_alt_id, score, pvalue, qvalue, matched_sequence, etc.
fimo_hits_targets <- memes::runFimo(
  sequences             = targets,
  motifs                = motifs_for_fimo,
  parse_genomic_coord   = FALSE,   # FASTA headers are not genomic coords here
  skip_matched_sequence = FALSE,   # set TRUE if you want smaller output; you can add sequence later
  thresh                = 1e-4,    # p-value threshold for reporting sites
  silent                = TRUE
)

# optional: also scan background to compare site densities
fimo_hits_bg <- memes::runFimo(
  sequences             = bg,
  motifs                = motifs_for_fimo,
  parse_genomic_coord   = FALSE,
  skip_matched_sequence = TRUE,
  thresh                = 1e-4,
  silent                = TRUE
)

# ---- 9) Summarize site-level enrichment (density and per-seq prevalence) -------
site_summary <- bind_rows(
  as.data.frame(mcols(fimo_hits_targets)) %>% mutate(where = "target"),
  as.data.frame(mcols(fimo_hits_bg))      %>% mutate(where = "background")
) %>%
  group_by(where, motif_id, motif_alt_id) %>%
  summarise(
    n_sites   = n(),
    median_p  = median(pvalue, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = where,
    values_from = c(n_sites, median_p),
    values_fill = 0
  ) %>%
  mutate(site_enrichment = (n_sites_target + 1) / (n_sites_background + 1)) %>%
  arrange(desc(site_enrichment))
readr::write_csv(site_summary, "FIMO_site_enrichment_summary.csv")

# ---- 10) Export per-motif, per-sequence hit counts (useful for downstream stats) -
hits_per_seq <- as.data.frame(mcols(fimo_hits_targets)) %>%
  count(motif_alt_id, sequence_name) %>%
  tidyr::pivot_wider(names_from = sequence_name, values_from = n, values_fill = 0)
readr::write_csv(hits_per_seq, "FIMO_hits_per_sequence.csv")

# ---- 11) Simple barplot of top enriched motifs by AME effect size --------------
# Here we plot pos_frac - neg_frac (prevalence difference) for quick ranking.
plot_df <- ame_res %>%
  mutate(effect = pos_frac - neg_frac) %>%
  arrange(desc(effect)) %>%
  slice_head(n = 20) %>%
  mutate(motif_label = ifelse(motif_alt_id == "", motif_id, motif_alt_id))
p_bar <- ggplot(plot_df, aes(x = reorder(motif_label, effect), y = effect)) +
  geom_col() +
  coord_flip() +
  labs(x = "Motif", y = "Prevalence difference (targets - background)",
       title = "Top enriched motifs (AME)") +
  theme_bw(base_size = 11)
ggsave("AME_top_motifs_barplot.pdf", p_bar, width = 6, height = 5)

# ---- 12) Save raw FIMO sites (targets) as BED-like TSV for genome-agnostic use --
fimo_tbl <- as.data.frame(fimo_hits_targets)
readr::write_tsv(fimo_tbl, "FIMO_sites_targets.tsv.gz")

# ---- 13) De novo discovery with STREME -----------------------------------------
# STREME tends to be fast and conservative; great first pass for de novo motifs.
streme_motifs <- memes::runStreme(
  sequences   = targets,         # 3'UTR targets
  control     = bg,               # matched background
  minw        = 5,                # min motif width
  maxw        = 12,               # max motif width (RBPs often 5–12 nt)
  nmotifs     = 15,               # number of motifs to find (tune as needed)
  thresh      = 0.05,             # E-value threshold for motif reporting
  pseudosites = 100,              # stability for small sets
  silent      = TRUE
)

# Save STREME motifs in MEME format (and as RDS) for portability
universalmotif::write_meme(streme_motifs, file = "streme_denovo.meme")
saveRDS(streme_motifs, file = "streme_denovo.rds")

# Quick visual QC: motif logos
pdf("streme_denovo_logos.pdf", width = 6, height = 8)
universalmotif::view_motifs(streme_motifs, use.names = "motif_alt_id")
dev.off()

# ---- 14) De novo discovery with MEME (complementary to STREME) -----------------
# MEME can find additional/longer motifs; 'zoops' = zero/one occurrence per sequence.
meme_motifs <- memes::runMeme(
  sequences = targets,
  control   = bg,
  mod       = "zoops",
  minw      = 5,
  maxw      = 12,
  nmotifs   = 15,
  nmotifs_search = 15,     # ensure search limit aligns with request
  evt       = 1.0,         # E-value threshold per motif
  bfile     = NULL,        # let MEME estimate background from control
  silent    = TRUE
)

universalmotif::write_meme(meme_motifs, file = "meme_denovo.meme")
saveRDS(meme_motifs, file = "meme_denovo.rds")

pdf("meme_denovo_logos.pdf", width = 6, height = 8)
universalmotif::view_motifs(meme_motifs, use.names = "motif_alt_id")
dev.off()

# ---- 15) Compare de novo motifs to CISBP-RNA with TomTom -----------------------
# This maps de novo motifs to nearest known RBPs (helps interpret families/redundancy).
# You can run separately for STREME and MEME outputs:

# STREME vs known
tomtom_streme <- memes::runTomTom(
  query    = streme_motifs,   # universalmotif list
  database = motifs_dna,      # your CISBP-RNA motifs (DNA alphabet)
  thresh   = 0.5,             # q-value threshold for reporting matches
  dist     = "ed",            # Euclidean distance on PWMs (common, robust)
  evalue   = FALSE,
  silent   = TRUE
)
readr::write_csv(tomtom_streme, "tomtom_streme_vs_cisbp.csv")

# MEME vs known
tomtom_meme <- memes::runTomTom(
  query    = meme_motifs,
  database = motifs_dna,
  thresh   = 0.5,
  dist     = "ed",
  evalue   = FALSE,
  silent   = TRUE
)
readr::write_csv(tomtom_meme, "tomtom_meme_vs_cisbp.csv")

# Convenience table: best match per de novo motif
best_tomtom_streme <- tomtom_streme |>
  dplyr::group_by(query_id) |>
  dplyr::slice_min(order_by = qvalue, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
best_tomtom_meme <- tomtom_meme |>
  dplyr::group_by(query_id) |>
  dplyr::slice_min(order_by = qvalue, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
readr::write_csv(best_tomtom_streme, "tomtom_streme_best.csv")
readr::write_csv(best_tomtom_meme,   "tomtom_meme_best.csv")

# ---- 16) Combine known + de novo for unified downstream testing ----------------
# Merge motif sets (drop exact duplicates; keep names informative)
denovo_all <- c(streme_motifs, meme_motifs)
denovo_all <- memes::remove_duplicate_motifs(denovo_all)   # exact dup removal

# Generate readable labels for de novo motifs (keep original IDs too)
denovo_all <- lapply(denovo_all, function(m) {
  if (isTRUE(is.na(m@name)) || m@name == "") m@name <- m@motif
  m
})
class(denovo_all) <- "universalmotif" # ensure list keeps class attribute

# ---- 17) Enrichment of de novo motifs with AME ---------------------------------
ame_denovo <- memes::runAme(
  input    = targets,
  control  = bg,
  database = denovo_all,
  method   = "fisher",
  scoring  = "avg",
  evalue_report_threshold = 10,
  silent   = TRUE
)
readr::write_csv(ame_denovo, "AME_denovo_targets_vs_background.csv")

# Plot heatmap of de novo enrichment
p_heat_denovo <- memes::plot_ame_heatmap(
  ame_denovo,
  id        = motif_alt_id,
  scale_max = 7.5
) + ggplot2::ggtitle("De novo motif enrichment (AME)")
ggplot2::ggsave("AME_denovo_enrichment_heatmap.pdf", p_heat_denovo, width = 6, height = 7)

# ---- 18) Site-level scanning of de novo motifs with FIMO -----------------------
# Targets
fimo_denovo_targets <- memes::runFimo(
  sequences             = targets,
  motifs                = denovo_all,
  parse_genomic_coord   = FALSE,
  skip_matched_sequence = FALSE,
  thresh                = 1e-4,
  silent                = TRUE
)
readr::write_tsv(as.data.frame(fimo_denovo_targets), "FIMO_denovo_sites_targets.tsv.gz")

# Background (for density comparison)
fimo_denovo_bg <- memes::runFimo(
  sequences             = bg,
  motifs                = denovo_all,
  parse_genomic_coord   = FALSE,
  skip_matched_sequence = TRUE,
  thresh                = 1e-4,
  silent                = TRUE
)

# Summarize site-level enrichment for de novo motifs
denovo_site_summary <- dplyr::bind_rows(
  as.data.frame(mcols(fimo_denovo_targets)) |> dplyr::mutate(where = "target"),
  as.data.frame(mcols(fimo_denovo_bg))      |> dplyr::mutate(where = "background")
) |>
  dplyr::group_by(where, motif_id, motif_alt_id) |>
  dplyr::summarise(n_sites = dplyr::n(), median_p = median(pvalue, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = where, values_from = c(n_sites, median_p), values_fill = 0) |>
  dplyr::mutate(site_enrichment = (n_sites_target + 1) / (n_sites_background + 1)) |>
  dplyr::arrange(dplyr::desc(site_enrichment))
readr::write_csv(denovo_site_summary, "FIMO_denovo_site_enrichment_summary.csv")

# ---- 19) Optional: cluster de novo motifs to reduce redundancy -----------------
# This uses PWM similarity to cluster and keep one representative per cluster.
# (Adjust 'score' threshold to be more/less aggressive.)
cmp <- universalmotif::compare_motifs(denovo_all, method = "PCC", min.overlap = 5)
hc  <- hclust(as.dist(1 - cmp), method = "average")
cl  <- cutree(hc, h = 0.25)                       # ~PCC >= 0.75 retained within clusters
keep_idx <- tapply(seq_along(denovo_all), cl, function(ix) ix[1]) |> unlist() |> sort()
denovo_representatives <- denovo_all[keep_idx]
universalmotif::write_meme(denovo_representatives, "denovo_cluster_representatives.meme")

# Re-run AME quickly on representatives (handy if you had many motifs)
ame_denovo_rep <- memes::runAme(
  input    = targets,
  control  = bg,
  database = denovo_representatives,
  method   = "fisher",
  scoring  = "avg",
  evalue_report_threshold = 10,
  silent   = TRUE
)

readr::write_csv(ame_denovo_rep, "AME_denovo_representatives.csv")