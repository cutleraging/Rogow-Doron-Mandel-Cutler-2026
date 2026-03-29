theme_Publication <- function(base_size=16, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family = "")
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#' Load RSEM processed RNA-seq data into a SummarizedExperiment
#'
#' This function reads metadata, TPM, and counts files produced by RSEM,
#' cleans and matches sample names, and returns a SummarizedExperiment
#' object containing both counts and TPM assays with sample metadata.
#'
#' @param metadata_file Metadata CSV filename (default: "sample_table.csv").
#'                      Must contain a column named "sample".
#' @param tpm_file RSEM merged gene TPM file (default: "rsem.merged.gene_tpm.tsv").
#' @param counts_file RSEM merged gene counts file (default: "rsem.merged.gene_counts.tsv").
#'
#' @return A SummarizedExperiment object with:
#'         - assays: counts (rounded integers) and TPM (numeric)
#'         - colData: sample metadata
#'         - rownames: gene IDs with version numbers stripped
#'
#' @examples
#' se <- load_rsem_data("/path/to/project/")
#'
load_rsem_data <- function(metadata_file = "sample_table.csv",
                           tpm_file = "rsem.merged.gene_tpm.tsv",
                           counts_file = "rsem.merged.gene_counts.tsv") {
  
  metadata <- read.csv(metadata_file, header = TRUE)
  rownames(metadata) <- metadata$sample
  
  # Helper to clean column names
  clean_names <- function(x) {
    x <- sub("^X", "", x)        # drop leading "X"
    x <- gsub("\\.", "-", x)     # replace "." with "-"
    return(x)
  }
  
  # Load TPM
  tpm <- read.table(tpm_file, header = TRUE, row.names = 1)
  tpm$transcript_id.s. <- NULL
  colnames(tpm) <- clean_names(colnames(tpm))
  keep <- intersect(colnames(tpm), metadata$sample)
  tpm_filtered <- tpm[, metadata$sample[metadata$sample %in% keep]]
  
  # Load counts
  counts <- read.table(counts_file, header = TRUE, row.names = 1)
  counts$transcript_id.s. <- NULL
  colnames(counts) <- clean_names(colnames(counts))
  keep <- intersect(colnames(counts), metadata$sample)
  counts_filtered <- counts[, metadata$sample[metadata$sample %in% keep]]
  
  # Build SummarizedExperiment
  se <- SummarizedExperiment(
    assays = SimpleList(
      counts = as.matrix(round(counts_filtered)),
      tpm = as.matrix(tpm_filtered)
    ),
    colData = metadata
  )
  
  # Remove version from gene names
  rownames(se) <- sub("\\.\\d+$", "", rownames(se))

  return(se)
}

#' Plot number of expressed genes per sample across multiple TPM thresholds
#'
#' This function calculates, for each sample in a SummarizedExperiment,
#' the number of genes with TPM values greater than specified thresholds
#' (default: 1, 2, 4, 8). It generates a stacked bar plot and saves it as a PDF.
#'
#' @param se A SummarizedExperiment object containing an assay named "tpm".
#' @param thresholds Numeric vector of TPM cutoffs to evaluate (default: c(1, 2, 4, 8)).
#' @param output_file Path to the PDF file where the plot will be saved
#'        (default: "genes_expressed_tpm_thresholds.pdf").
#' @param width Width of the output PDF in inches (default: 6).
#' @param height Height of the output PDF in inches (default: 5).
#'
#' @return Invisibly returns the data.frame used for plotting
#'         with columns: Sample, Threshold, Frequency.
#'
#' @examples
#' plot_genes_expressed_thresholds(se)
#' plot_genes_expressed_thresholds(se, thresholds = c(1, 5, 10))
#'
plot_genes_expressed <- function(se,
                                 thresholds = c(1, 2, 4, 8),
                                 output_file = "genes_expressed.pdf",
                                 width = 6,
                                 height = 5) {
  # Calculate number of genes above each threshold for every sample
  df_list <- lapply(thresholds, function(th) {
    data.frame(
      Sample = colnames(se),
      Threshold = paste0("TPM > ", th),
      Frequency = colSums(assay(se, "tpm") > th)
    )
  })
  
  # Combine into one data.frame
  df <- do.call(rbind, df_list)
  
  # Save stacked bar plot to PDF
  pdf(output_file, width = width, height = height)
  p <- ggplot(df, aes(x = Sample, y = Frequency, fill = Threshold)) +
    geom_col(position = "stack") +
    theme_Publication() +
    scale_fill_Publication() +
    coord_flip() +
    theme(legend.position = "right",
          legend.direction = "vertical") +
    labs(subtitle = "Number of genes detected",
         y = "# genes")
  print(p)
  dev.off()
  
  invisible(df)  # return data.frame invisibly
}

#' Plot count distributions per sample (SE or DESeq2) with median-alignment boxplot
#'
#' Creates a two-page PDF:
#'   Page 1 — ridgeline densities of log2(counts + pseudocount) per sample
#'   Page 2 — boxplots per sample with a reference line at the global median
#' Works with a SummarizedExperiment (any assay) or a DESeq2::DESeqDataSet (normalized counts).
#'
#' @param se A SummarizedExperiment (optional if \code{dds} is provided).
#' @param dds A DESeq2::DESeqDataSet (optional if \code{se} is provided).
#' @param count_assay Assay name in \code{se} to use (default: "counts").
#' @param normalized When using \code{dds}, passed to \code{counts(dds, normalized=...)} (default: TRUE).
#' @param pseudocount Value added before log transform (default: 1).
#' @param output_file PDF path (default: "count_distribution.pdf").
#' @param width,height PDF size in inches (defaults: 6, 8).
#'
#' @return Invisibly, the long-format data.frame with columns \code{Sample}, \code{Value}.
#'
#' @examples
#' # From SummarizedExperiment
#' plot_count_distribution(se)
#' # From DESeq2 normalized counts
#' plot_count_distribution(dds = dds)
plot_count_distribution <- function(se = NULL,
                                    dds = NULL,
                                    count_assay = "counts",
                                    normalized = TRUE,
                                    pseudocount = 1,
                                    output_file = "count_distribution.pdf",
                                    width = 6,
                                    height = 8) {
  if (!requireNamespace("ggridges", quietly = TRUE)) stop("Install 'ggridges'.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Install 'tidyr'.")
  if (is.null(se) && is.null(dds)) stop("Provide either 'se' or 'dds'.")
  
  # --- extract matrix ---
  if (!is.null(dds)) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install 'DESeq2'.")
    mat <- DESeq2::counts(dds, normalized = normalized)
    smp <- colnames(dds)
  } else {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Install 'SummarizedExperiment'.")
    if (!count_assay %in% SummarizedExperiment::assayNames(se)) {
      stop(sprintf("Assay '%s' not found in 'se'.", count_assay))
    }
    mat <- SummarizedExperiment::assay(se, count_assay)
    smp <- colnames(se)
  }
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  # --- transform & gather ---
  df <- as.data.frame(log2(mat + pseudocount), check.names = FALSE)
  colnames(df) <- smp
  df <- tidyr::pivot_longer(df, cols = tidyselect::everything(),
                            names_to = "Sample", values_to = "Value")
  df$Sample <- factor(df$Sample, levels = smp)
  
  # Global median of per-sample medians (for alignment check)
  sample_meds <- tapply(df$Value, df$Sample, median, na.rm = TRUE)
  global_med  <- median(sample_meds, na.rm = TRUE)
  
  grDevices::pdf(output_file, width = width, height = height)
  
  # --- Page 1: ridgeline densities ---
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = Value, y = Sample)) +
    ggridges::geom_density_ridges(scale = 2) +
    theme_Publication() +
    ggplot2::labs(
      x = "log2(counts + pseudocount)",
      y = "Sample",
      subtitle = "Count Distributions (ridgeline)"
    ) +
    ggplot2::theme(legend.position = "none")
  print(p1)
  
  # --- Page 2: boxplots with global-median reference line ---
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value)) +
    ggplot2::geom_boxplot(outlier.size = 0.4, width = 0.6) +
    ggplot2::geom_hline(yintercept = global_med, linetype = 2) +
    theme_Publication() +
    coord_flip() +
    ggplot2::labs(
      x = "Sample",
      y = "log2(counts + pseudocount)",
      subtitle = "Count Distributions (boxplots)"
    ) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  print(p2)
  
  grDevices::dev.off()
  invisible(df)
}


#' Filter SummarizedExperiment to active genes using zFPKM
#'
#' Computes zFPKM from TPM, writes a zFPKM distribution PDF,
#' and filters genes to those called active (median zFPKM > threshold
#' within each group). Returns a new SummarizedExperiment containing
#' only active genes.
#'
#' @param se A SummarizedExperiment with an assay named "tpm".
#' @param group_cols Character vector of column names in colData(se) to group by
#'   for median zFPKM (default: c("diet","sex")).
#' @param active_threshold Numeric zFPKM cutoff for activity (default: -3).
#' @param output_file PDF path for the zFPKM distribution plot
#'   (default: "zFPKM_distribution.pdf").
#' @param width,height PDF dimensions in inches (defaults: 10, 10).
#' @param facet_titles Logical; pass to zFPKM::zFPKMPlot `FacetTitles` (default: TRUE).
#'
#' @return A filtered SummarizedExperiment containing only active genes.
#'
#' @examples
#' se_active <- filter_active_genes(se)
#'
filter_active_genes <- function(
    se,
    group_cols = c("condition"),
    active_threshold = -3,
    output_file = "zFPKM_distribution.pdf",
    width = 10,
    height = 10,
    facet_titles = TRUE
) {
  # Dependencies check
  if (!requireNamespace("zFPKM", quietly = TRUE)) {
    stop("Package 'zFPKM' is required. Install via Bioconductor: BiocManager::install('zFPKM').")
  }
  for (pkg in c("dplyr","tidyr","tibble")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Please install it.", pkg))
    }
  }
  
  # Ensure the TPM assay exists
  if (!"tpm" %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay 'tpm' not found in 'se'.")
  }
  
  # compute zFPKM, store assay
  SummarizedExperiment::assay(se, "zfpkm") <- zFPKM::zFPKM(se, assayName = "tpm")
  
  # >>> FIX: coerce to base matrix before replacing Inf <<<
  zmat <- SummarizedExperiment::assay(se, "zfpkm")
  if (!is.matrix(zmat)) zmat <- as.matrix(zmat)     # handles DelayedMatrix/sparse Matrix
  storage.mode(zmat) <- "double"
  zmat[!is.finite(zmat)] <- NA_real_
  SummarizedExperiment::assay(se, "zfpkm") <- zmat
  
  # Plot zFPKM distribution to PDF
  grDevices::pdf(output_file, width = width, height = height)
  print(zFPKM::zFPKMPlot(se, assayName = "tpm", FacetTitles = facet_titles, PlotXfloor = -10))
  grDevices::dev.off()
  
  # Prepare zFPKM long table
  z_long <-
    zmat |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "zfpkm")
  
  # Bring in sample-level grouping metadata (avoid duplicating 'sample')
  sample_info <- SummarizedExperiment::colData(se) |>
    as.data.frame(check.names = FALSE)
  
  # If no 'sample' column, create it from rownames; otherwise keep existing
  if (!"sample" %in% colnames(sample_info)) {
    sample_info$sample <- rownames(sample_info)
  } else {
    # ensure character and consistent with column names of se
    sample_info$sample <- as.character(sample_info$sample)
  }
  
  # Sanity checks
  if (anyDuplicated(sample_info$sample)) {
    stop("Duplicate entries found in colData(se)$sample; sample IDs must be unique.")
  }
  if (!all(colnames(se) %in% sample_info$sample)) {
    stop("Some assay column names are missing from colData(se)$sample.")
  }
  
  # Validate grouping columns
  missing_cols <- setdiff(group_cols, colnames(sample_info))
  if (length(missing_cols) > 0) {
    stop(sprintf("Grouping columns not in colData(se): %s",
                 paste(missing_cols, collapse = ", ")))
  }
  
  # Compute median zFPKM per gene within groups
  z_grouped <-
    z_long |>
    dplyr::left_join(sample_info[, c("sample", group_cols), drop = FALSE], by = "sample") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("gene", group_cols)))) |>
    dplyr::summarise(median_zfpkm = stats::median(zfpkm, na.rm = TRUE), .groups = "drop")
  
  # Call genes active if median zFPKM exceeds threshold
  active_genes <-
    z_grouped |>
    dplyr::filter(median_zfpkm > active_threshold) |>
    dplyr::distinct(gene) |>
    dplyr::pull(gene)
  
  # Filter SE to active genes
  se_filtered <- se[rownames(se) %in% active_genes, ]
  
  return(se_filtered)
}

#' Plot within-group sample correlations using chart.Correlation
#'
#' Computes log2(TPM + 1) from a specified assay, splits samples by a grouping
#' column in `colData(se)`, and for each group with ≥ `min_reps` samples,
#' draws a correlation matrix plot (pairwise scatter, histograms, and r values)
#' using **PerformanceAnalytics::chart.Correlation**. One multi-page PDF is written,
#' with one page per group. Zero-variance genes within a group are removed to
#' avoid undefined correlations.
#'
#' @param se A \code{SummarizedExperiment} containing an assay with TPM-like values.
#' @param group_col Character scalar; name of the column in \code{colData(se)} used
#'   to group samples (e.g., "condition", "type"). Each level is plotted on its own page.
#' @param assay_name Character; assay to use (default: "tpm").
#' @param output_file Character; path to the output PDF (default: "sample_correlations_by_group.pdf").
#' @param method Correlation method passed to \code{stats::cor} ("pearson", "spearman", "kendall").
#'   Used internally by \code{chart.Correlation} (default: "pearson").
#' @param min_reps Integer; minimum number of samples required in a group to plot (default: 2).
#' @param pseudocount Numeric; value added before log2-transform (default: 1).
#' @param plot_hist Logical; whether to include histograms on the diagonal (default: TRUE).
#'
#' @return Invisibly returns a named list whose elements are the correlation matrices
#'   (\code{matrix}) for each plotted group.
#'
#' @examples
#' # Plot correlations per condition using TPM assay
#' cors <- plot_replicate_correlations(
#'   se, group_col = "condition", assay_name = "tpm",
#'   output_file = "cor_by_condition.pdf"
#' )
#'
plot_replicate_correlations <- function(
    se,
    group_col,
    assay_name = "tpm",
    output_file = "sample_correlations_by_group.pdf",
    method = c("pearson", "spearman", "kendall"),
    min_reps = 2,
    pseudocount = 1,
    plot_hist = TRUE
) {
  # --- Dependency checks ---
  if (!requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
    stop("Package 'PerformanceAnalytics' is required. Install with install.packages('PerformanceAnalytics').")
  }
  
  # --- Input validation ---
  method <- match.arg(method)
  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("Assay '%s' not found in se.", assay_name))
  }
  cd <- SummarizedExperiment::colData(se)
  if (!group_col %in% colnames(cd)) {
    stop(sprintf("Grouping column '%s' not found in colData(se).", group_col))
  }
  
  # --- Extract and transform data: log2(TPM + pseudocount) ---
  mat <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat <- log2(mat + pseudocount)  # genes x samples
  
  # --- Split samples by group ---
  groups <- as.character(cd[[group_col]])
  names(groups) <- colnames(se)
  lvl <- unique(groups)
  
  # Prepare output and collector
  grDevices::pdf(output_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  cor_list <- list()
  
  # --- Iterate over groups, plot one page per group ---
  for (g in lvl) {
    smps <- names(groups)[groups == g]
    if (length(smps) < min_reps) next  # skip undersized groups
    
    # Subset matrix to samples in this group
    X <- mat[, smps, drop = FALSE]
    
    # Drop genes with zero variance across these samples to avoid NA correlations
    # (var = 0 or all NA after transform)
    v <- apply(X, 1, stats::var, na.rm = TRUE)
    keep_genes <- is.finite(v) & v > 0
    X <- X[keep_genes, , drop = FALSE]
    
    # If still too few genes, skip
    if (nrow(X) < 2) next
    
    # Transpose to samples-as-columns for chart.Correlation (already true),
    # but ensure it's a plain data.frame
    df <- as.data.frame(X, check.names = FALSE)
    
    # Compute correlation matrix for return value
    cor_mat <- stats::cor(df, use = "pairwise.complete.obs", method = method)
    cor_list[[g]] <- cor_mat
    
    # Plot (base graphics)
    graphics::par(oma = c(0, 0, 3, 0))  # outer margin for title
    PerformanceAnalytics::chart.Correlation(
      df,
      method = method,
      histogram = plot_hist,
      pch = 16
    )
  }
  
  # Warn if nothing was plotted
  if (length(cor_list) == 0) {
    warning("No groups had at least ", min_reps, " samples to plot. No pages written.")
  }
  
  invisible(cor_list)
}

#' UpSet plot of genes detected across samples
#'
#' Builds a binary detection matrix from a SummarizedExperiment assay
#' (genes × samples), where detection = (value > threshold), and
#' draws an UpSet plot showing overlaps across samples.
#'
#' @param se A SummarizedExperiment.
#' @param assay_name Assay to use (default: "tpm").
#' @param threshold Numeric detection cutoff; values strictly greater than this
#'   are considered detected (default: 0).
#' @param output_file PDF path for the plot
#'   (default: "overlap_expressed_genes_top_10.pdf").
#' @param nintersects Number of top intersections to display (default: 10).
#' @param order_by How to order intersections; passed to UpSetR
#'   (default: "freq").
#' @param text_scale Numeric text scaling for UpSetR (default: 1.25).
#' @param show_numbers Logical; show intersection sizes on bars (default: FALSE).
#' @param width,height PDF size in inches (defaults: 10, 7).
#'
#' @return Invisibly returns the binary detection data.frame (genes × samples).
#'
#' @examples
#' plot_genes_expressed_overlap(se.filter, assay_name = "tpm", threshold = 0)
#'
plot_genes_expressed_overlap <- function(
    se,
    assay_name = "tpm",
    threshold = 0,
    output_file = "overlap_expressed_genes_top_10.pdf",
    nintersects = 10,
    order_by = "freq",
    text_scale = 1.25,
    show_numbers = FALSE,
    width = 10,
    height = 7
) {
  # Dependency check
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' is required. Install with install.packages('UpSetR').")
  }
  
  # Validate assay
  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("Assay '%s' not found in 'se'.", assay_name))
  }
  
  # Extract assay matrix and compute detection (genes × samples)
  mat <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(mat)) mat <- as.matrix(mat)     # handle DelayedMatrix/sparse
  storage.mode(mat) <- "double"
  binary_mat <- (mat > threshold)
  # Coerce logical to integer 0/1 to match your example
  binary_mat <- matrix(as.integer(binary_mat), nrow = nrow(binary_mat), dimnames = dimnames(binary_mat))
  
  # Build data.frame for UpSetR (rows = genes, columns = samples)
  binary_df <- as.data.frame(binary_mat, check.names = FALSE)
  
  # Use nicer sample labels if available in colData(se)$sample
  sample_labels <- if ("sample" %in% colnames(SummarizedExperiment::colData(se))) {
    as.character(SummarizedExperiment::colData(se)$sample)
  } else {
    colnames(se)
  }
  # Keep ordering consistent with columns
  if (length(sample_labels) == ncol(binary_df)) colnames(binary_df) <- sample_labels
  
  # Write PDF
  grDevices::pdf(output_file, width = width, height = height)
  UpSetR::upset( 
    binary_df,
    sets = colnames(binary_df),
    text.scale = text_scale,
    nintersects = nintersects,
    order.by = order_by,
    show.numbers = show_numbers
  )
  grDevices::dev.off()
  
  invisible(binary_df)
}

#' Connect to the current Ensembl release for a given species
#'
#' Queries Ensembl archives, selects the **current** release automatically,
#' finds the matching \code{*_gene_ensembl} dataset for the requested species,
#' prints a message with the chosen release/date/URL, and returns a live
#' \code{biomaRt} connection.
#'
#' @param species Character. Species name or key (e.g., "hsapiens", "Homo sapiens",
#'   "mmusculus", "Mus musculus"). Case and spacing are ignored.
#' @param biomart Character. BiomaRt to use (default: "genes").
#' @param quiet Logical. If FALSE, prints the chosen release/dataset (default: FALSE).
#'
#' @return A \code{biomaRt} Mart object connected to the current Ensembl release
#'   for the requested species.
#'
#' @examples
#' mart_hs <- get_latest_ensembl_mart("hsapiens")
#' mart_mm <- get_latest_ensembl_mart("Mus musculus")
#'
get_latest_ensembl_mart <- function(species, biomart = "genes", quiet = FALSE) {
  # --- deps ---
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Install with BiocManager::install('biomaRt').")
  }
  
  # --- helper to normalize strings for matching ---
  .norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  
  # --- get archive table and choose current release (fallback: highest numeric) ---
  arch <- biomaRt::listEnsemblArchives()
  if (!"current_release" %in% names(arch)) {
    stop("Could not retrieve Ensembl archive metadata (no 'current_release' column).")
  }
  i_cur <- which(arch$current_release == "*")
  if (length(i_cur) == 1) {
    sel <- arch[i_cur, , drop = FALSE]
  } else {
    # fallback: pick max numeric 'version' (ignore GRCh37 row)
    ver_num <- suppressWarnings(as.numeric(arch$version))
    sel <- arch[which.max(ver_num), , drop = FALSE]
  }
  ver <- suppressWarnings(as.integer(sel$version))
  rel_date <- sel$date
  rel_url  <- sel$url
  
  # --- open mart for that release ---
  mart_rel <- biomaRt::useEnsembl(biomart = biomart, version = ver)
  
  # --- find dataset for species ---
  ds <- biomaRt::listDatasets(mart_rel)
  keys <- sub("_gene_ensembl$", "", ds$dataset)
  # try exact key match (e.g., "hsapiens", "mmusculus")
  hit <- which(.norm(keys) == .norm(species))
  if (length(hit) == 0) {
    # try description contains species name (e.g., "Homo sapiens")
    hit <- grep(.norm(species), .norm(ds$description), fixed = TRUE)
  }
  if (length(hit) == 0) {
    stop(sprintf(
      "Could not find a dataset for species '%s' in Ensembl release %s.\nExamples: %s",
      species, ver %||% "?", paste(head(ds$dataset, 5), collapse = ", ")
    ))
  }
  # if multiple, take the first deterministically
  hit <- hit[1]
  dataset <- ds$dataset[hit]
  species_desc <- ds$description[hit]
  
  # --- final mart with dataset selected ---
  mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, version = ver)
  
  # --- message ---
  if (!quiet) {
    msg <- sprintf(
      "Using Ensembl release %s (%s) at %s | dataset: %s [%s]",
      if (!is.na(ver)) ver else sel$version, rel_date, rel_url, dataset, species_desc
    )
    message(msg)
  }
  
  return(mart)
}

# small infix helper (avoid importing rlang)
`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b

#' Get the latest EnsDb for a species from AnnotationHub (by date)
#'
#' Queries AnnotationHub for EnsDb resources for the requested species,
#' selects the most recently added entry (by \code{rdatadateadded}),
#' prints a message describing the chosen resource, and returns the EnsDb.
#'
#' @param species Character. Species name or key
#'   (e.g., "hsapiens", "Homo sapiens", "mmusculus", "Mus musculus").
#'   Case/spacing ignored.
#' @param prefer_genome Optional character to prefer a specific genome build
#'   when multiple entries share the latest date (e.g., "GRCh38", "GRCh37",
#'   "GRCm39"). If not found at the latest date, the function falls back to
#'   the strictly latest entry by date.
#' @param quiet Logical; if FALSE, prints the chosen date, Ensembl version,
#'   genome, and AnnotationHub ID (default: FALSE).
#'
#' @return An \code{EnsDb} object.
#'
#' @examples
#' edb_hs <- get_latest_ensdb("hsapiens")
#' edb_mm <- get_latest_ensdb("Mus musculus", prefer_genome = "GRCm39")
#'
get_latest_ensdb <- function(species, prefer_genome = NULL, quiet = FALSE) {
  # deps
  if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("Install 'AnnotationHub' (BiocManager::install('AnnotationHub')).")
  }
  if (!requireNamespace("ensembldb", quietly = TRUE)) {
    stop("Install 'ensembldb' (BiocManager::install('ensembldb')).")
  }
  
  # normalize species for matching
  .norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  species_norm <- .norm(species)
  
  ah <- AnnotationHub::AnnotationHub()
  
  # query and keep EnsDb records only
  qry <- AnnotationHub::query(ah, c("EnsDb"))
  md  <- S4Vectors::mcols(qry)
  
  # filter by species using both 'species' and 'title' fields
  keep <- (.norm(md$species) == species_norm) |
    grepl(species_norm, .norm(md$title), fixed = TRUE)
  qry <- qry[keep]
  if (length(qry) == 0L) {
    stop(sprintf("No EnsDb entries found in AnnotationHub for species '%s'.", species))
  }
  
  md <- S4Vectors::mcols(qry)
  # ensure we have dates
  dates <- as.Date(md$rdatadateadded)
  if (all(is.na(dates))) {
    stop("No 'rdatadateadded' available; cannot choose the latest by date.")
  }
  
  # pick the max date; if prefer_genome provided and present at that date, use it
  max_date <- max(dates, na.rm = TRUE)
  at_latest <- which(dates == max_date)
  
  idx <- at_latest[1]  # default choice among ties
  if (!is.null(prefer_genome)) {
    # try to find preferred genome among latest-date ties
    g <- md$genome[at_latest]
    hit <- which(.norm(g) == .norm(prefer_genome))
    if (length(hit)) idx <- at_latest[hit[1]]
  }
  
  # fallback: among ties with same date, choose highest Ensembl version if available
  if (length(at_latest) > 1) {
    ev <- suppressWarnings(as.numeric(md$ensembl_version[at_latest]))
    if (any(!is.na(ev))) {
      idx <- at_latest[which.max(ev)]
    }
  }
  
  # retrieve EnsDb
  edb <- qry[[idx]]
  
  # message
  if (!quiet) {
    msg <- sprintf(
      "AnnotationHub EnsDb selected: AH%s | species: %s | date: %s | genome: %s | Ensembl: %s",
      names(qry)[idx],
      md$species[idx],
      md$rdatadateadded[idx],
      md$genome[idx] %||% "unknown",
      tryCatch(ensembldb::ensemblVersion(edb), error = function(e) md$ensembl_version[idx] %||% "unknown")
    )
    message(msg)
  }
  
  edb
}

#' Annotate SE with Ensembl (BioMart) metadata + EnsDb ranges
#'
#' Allows different species identifiers for BioMart and EnsDb.
#' Prints how many IDs mapped in BioMart and how many ranges mapped.
#'
#' @param se SummarizedExperiment with Ensembl Gene ID rownames (versions ok).
#' @param ensemble_mart_species Character for BioMart (e.g., "hsapiens", "Homo sapiens").
#' @param ensembldb_species Character for EnsDb / AnnotationHub (e.g., "Homo sapiens").
#' @param biomart_attrs Character vector of BioMart attributes to retrieve.
#' @param prefer_genome Optional genome build preference for EnsDb (e.g., "GRCh38").
#' @param drop_unmapped_ranges Logical; drop genes without ranges (default TRUE).
#' @param quiet Logical; suppress messages.
#'
#' @return SummarizedExperiment with rowData (BioMart) and rowRanges (EnsDb).
annotate_se_ensembl <- function(
    se,
    ensemble_mart_species = "hsapiens",
    ensembldb_species     = "Homo sapiens",
    biomart_attrs = c("ensembl_gene_id","version","external_gene_name","entrezgene_id",
                      "uniprotswissprot","gene_biotype","chromosome_name",
                      "start_position","end_position","strand","description"),
    prefer_genome = NULL,
    drop_unmapped_ranges = TRUE,
    quiet = FALSE
) {
  for (pkg in c("SummarizedExperiment","S4Vectors","biomaRt","AnnotationFilter",
                "GenomicRanges","ensembldb","AnnotationHub","dplyr","IRanges")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Install '%s'.", pkg))
  }
  
  # helper: safe coerce to length-1 character
  as_scalar_chr <- function(x, nm) {
    if (length(x) != 1L) stop(sprintf("'%s' must be length 1.", nm))
    as.character(x)
  }
  ensemble_mart_species <- as_scalar_chr(ensemble_mart_species, "ensemble_mart_species")
  ensembldb_species     <- as_scalar_chr(ensembldb_species, "ensembldb_species")
  
  ids_full <- rownames(se); if (is.null(ids_full)) stop("rownames(se) must be Ensembl gene IDs.")
  ids_core <- sub("\\.\\d+$", "", ids_full)
  n_ids <- length(ids_core)
  
  # --- BioMart metadata ---
  mart <- get_latest_ensembl_mart(species = ensemble_mart_species, biomart = "genes", quiet = quiet)
  tbl <- biomaRt::getBM(
    attributes = biomart_attrs,
    filters    = "ensembl_gene_id",
    values     = ids_core,
    mart       = mart
  )
  
  m_bm <- match(ids_core, tbl$ensembl_gene_id)
  if (!quiet) message("BioMart mapped ", sum(!is.na(m_bm)), " / ", n_ids, " Ensembl IDs.")
  S4Vectors::metadata(se)$biomart_species <- ensemble_mart_species
  S4Vectors::metadata(se)$biomart_attrs <- biomart_attrs
  SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(tbl[m_bm, ])
  
  # --- EnsDb ranges ---
  edb <- get_latest_ensdb(ensembldb_species, prefer_genome = prefer_genome, quiet = quiet)
  gr <- GenomicRanges::GRanges(
    ensembldb::genes(edb, filter = AnnotationFilter::GeneIdFilter(ids_core))
  )
  
  ids_gr <- sub("\\.\\d+$", "", gr$gene_id)
  keep <- !duplicated(ids_gr)
  gr2 <- gr[keep]; ids_gr2 <- ids_gr[keep]
  
  o <- match(ids_core, ids_gr2)
  if (!quiet) message("GRanges mapped ", sum(!is.na(o)), " / ", n_ids, " Ensembl IDs.")
  
  if (drop_unmapped_ranges) {
    if (sum(!is.na(o)) == 0L) stop("No genomic ranges could be mapped.")
    keep_rows <- !is.na(o)
    if (any(!keep_rows) && !quiet) {
      miss_ids <- unique(ids_core[!keep_rows])
      message("Dropping ", sum(!keep_rows), " rows without ranges. Example IDs: ",
              paste(utils::head(miss_ids, 5), collapse = ", "), " ...")
    }
    se <- se[keep_rows, ]; o <- o[keep_rows]
    SummarizedExperiment::rowRanges(se) <- gr2[o]
  } else {
    rr <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(factor(NA, levels = GenomicRanges::seqlevels(gr2))),
      ranges   = IRanges::IRanges(start = NA_integer_, end = NA_integer_),
      strand   = S4Vectors::Rle("*")
    )
    rr <- rr[rep(1, length(ids_core))]
    matched <- which(!is.na(o)); rr[matched] <- gr2[o[matched]]
    SummarizedExperiment::rowRanges(se) <- rr
  }
  
  se
}

#' Export assay (SE) or normalized counts (DESeq2) to CSV, optionally merged with rowData
#'
#' Exports either:
#'  - an assay from a SummarizedExperiment (\code{se}, e.g., "tpm" or "counts"), or
#'  - normalized counts from a DESeq2 \code{dds} (\code{counts(dds, normalized=TRUE)}),
#' and (optionally) merges these matrices with \code{rowData} from a provided SE.
#'
#' If both \code{dds} and \code{se} are provided, the exported matrix comes from \code{dds}
#' and \code{rowData(se)} is merged by gene ID (rownames). The function preserves the
#' original gene order and collapses any list-like columns in \code{rowData}.
#'
#' @param se SummarizedExperiment providing either the assay to export (when \code{dds=NULL}),
#'   or the \code{rowData} to merge (when \code{dds} is provided).
#' @param dds DESeq2::DESeqDataSet; when provided, exports \code{counts(dds, normalized=...)}.
#' @param assay_name Assay name in \code{se} to export if \code{dds=NULL} (default: "tpm").
#' @param normalized Logical; passed to \code{counts(dds, normalized = normalized)} (default: TRUE).
#' @param file Output CSV path (default: "export.csv").
#' @param collapse_sep Separator for collapsing list-like columns in \code{rowData} (default: ";").
#'
#' @return (Invisibly) the written file path.
#'
#' @examples
#' # Export TPM from SE:
#' # export_se_assay_or_dds_csv(se, assay_name = "tpm", file = "tpm_filtered.csv")
#'
#' # Export DESeq2 normalized counts merged with rowData from SE:
#' # export_se_assay_or_dds_csv(se = se.filter, dds = dds,
#' #                            file = "counts.filtered.normalized.csv")
export_se_assay_csv <- function(se = NULL,
                               dds = NULL,
                               assay_name = "tpm",
                               normalized = TRUE,
                               file = "export.csv",
                               collapse_sep = ";") {
  # --- input checks ---
  if (is.null(dds) && is.null(se)) {
    stop("Provide either 'se' (to export an assay) or 'dds' (to export normalized counts).")
  }
  
  # --- build the expression matrix to export ---
  if (!is.null(dds)) {
    # Export normalized counts from DESeq2
    if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install 'DESeq2' to use 'dds'.")
    mat <- DESeq2::counts(dds, normalized = normalized)
    # Gene IDs from rownames(dds); samples from colnames(dds)
    gene_ids <- rownames(mat)
    sample_names <- colnames(mat)
  } else {
    # Export an assay from SE
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Install 'SummarizedExperiment' to use 'se'.")
    }
    if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
      stop(sprintf("Assay '%s' not found in 'se'.", assay_name))
    }
    mat <- SummarizedExperiment::assay(se, assay_name)
    gene_ids <- rownames(se)
    sample_names <- colnames(se)
  }
  
  # Coerce to plain matrix/data.frame for safe I/O
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  assay_df <- as.data.frame(mat, check.names = FALSE)
  colnames(assay_df) <- sample_names
  assay_df$gene_id <- gene_ids
  
  # --- prepare rowData (optional merge) ---
  rd_df <- NULL
  if (!is.null(se)) {
    rd <- SummarizedExperiment::rowData(se)
    rd_df <- as.data.frame(rd, stringsAsFactors = FALSE, check.names = FALSE)
    rd_df$gene_id <- rownames(se)
    
    # Collapse list-like columns (e.g., CharacterList) to atomic strings
    is_list_col <- vapply(rd_df, function(x) is.list(x) || inherits(x, "List"), logical(1))
    if (any(is_list_col)) {
      rd_df[is_list_col] <- lapply(rd_df[is_list_col], function(col) {
        vapply(col, function(el) paste(as.character(el), collapse = collapse_sep), character(1))
      })
    }
  }
  
  # --- merge assay/counts with rowData if available ---
  if (!is.null(rd_df)) {
    out <- merge(rd_df, assay_df, by = "gene_id", sort = FALSE, all.y = TRUE)
    # Reorder to match assay order
    ord <- match(assay_df$gene_id, out$gene_id)
    out <- out[ord, , drop = FALSE]
  } else {
    # No rowData provided: just include gene_id + counts
    out <- assay_df[, c("gene_id", sample_names), drop = FALSE]
  }
  
  # --- write CSV ---
  utils::write.csv(out, file = file, row.names = FALSE)
  invisible(file)
}

#' VST → PCA with scree, eigencor, and biplot (PCAtools)
#'
#' Runs DESeq2 VST (configurable), computes PCA with PCAtools, and writes:
#'   1) Scree plot (up to 10 PCs, or fewer if not available)
#'   2) Eigencor plot against selected metadata
#'   3) PCA biplot with user-chosen color/shape aesthetics
#'
#' @param dds DESeq2::DESeqDataSet.
#' @param blind Logical; passed to DESeq2::vst (default: TRUE).
#' @param removeVar Proportion of lowest-variance genes to remove before PCA (default: 0.1).
#' @param metavars Character vector of metadata columns to include in eigencorplot.
#'                 Use NULL to include all columns of colData(dds_vst). (Default: c("condition","type","group"))
#' @param colby Metadata column name for color in biplot (default: "type").
#' @param shape Metadata column name for point shape in biplot (default: "group").
#' @param scree_file,eigencor_file,biplot_file Output PDF paths.
#' @param scree_size,eigencor_size,biplot_size Numeric vectors c(width,height) in inches.
#'
#' @return Invisibly returns a list with \code{dds_vst}, \code{p} (PCAtools object),
#'         and \code{npc_plotted} (integer).
#'
#' @examples
#' res <- pca_with_reports(dds,
#'   blind = TRUE, removeVar = 0.1,
#'   metavars = c("condition","type","group"),
#'   colby = "type", shape = "group"
#' )
pca_with_reports <- function(
    dds,
    blind = TRUE,
    removeVar = 0.1,
    metavars = c("condition","type","group"),
    colby = "type",
    shape = "group",
    scree_file    = "screeplot.pdf",
    eigencor_file = "eigencor.pdf",
    biplot_file   = "PCA.pdf",
    scree_size    = c(5, 5),
    eigencor_size = c(5, 5),
    biplot_size   = c(6, 5)
) {
  # ---- deps ----
  for (pkg in c("DESeq2","PCAtools","SummarizedExperiment","ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Install '%s'.", pkg))
  }
  # palette
  vir <- if (requireNamespace("viridisLite", quietly = TRUE)) {
    viridisLite::viridis(100)
  } else if (requireNamespace("viridis", quietly = TRUE)) {
    viridis::viridis(100)
  } else {
    grDevices::colorRampPalette(c("#440154","#31688e","#35b779","#fde725"))(100)
  }
  
  # ---- VST ----
  dds_vst <- DESeq2::vst(dds, blind = blind)
  
  # ---- PCA (PCAtools) ----
  mat <- SummarizedExperiment::assay(dds_vst)
  meta <- as.data.frame(SummarizedExperiment::colData(dds_vst))
  
  p <- PCAtools::pca(
    mat,
    metadata  = meta,
    removeVar = removeVar
  )
  
  # How many PCs exist? cap at 10 for plots
  npc_total <- ncol(p$rotated)
  npc <- max(1, min(10, npc_total))
  
  # choose metavars
  metavars_use <- if (is.null(metavars)) colnames(meta) else {
    keep <- metavars[metavars %in% colnames(meta)]
    if (!length(keep)) stop("None of the specified 'metavars' are in colData.")
    keep
  }
  
  # ---- Scree plot ----
  grDevices::pdf(scree_file, width = scree_size[1], height = scree_size[2])
  sp <- PCAtools::screeplot(p)
  print(sp)
  grDevices::dev.off()
  
  # ---- Eigencor plot ----
  grDevices::pdf(eigencor_file, width = eigencor_size[1], height = eigencor_size[2])
  ep <- PCAtools::eigencorplot(
    p,
    components = PCAtools::getComponents(p, seq_len(npc)),
    metavars   = metavars_use,
    main       = "Principal component correlations",
    cexMain    = 1.5,
    col        = vir,
    colCorval  = "firebrick",
    fontCorval = 2,
    cexCorval  = 0.5,
    rotLabX    = 45,
    posColKey  = "top"
  )
  print(ep)
  grDevices::dev.off()
  
  # ---- Biplot ----
  grDevices::pdf(biplot_file, width = biplot_size[1], height = biplot_size[2])
  gp <- PCAtools::biplot(
    p,
    colby  = colby,
    shape  = shape,
    hline  = 0, vline = 0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    max.overlaps    = Inf,
    lab = rownames(p$metadata),
    labSize = 1,
    title = "Principle Component Analysis"
  )
  # Add your publication theme if available
  if (exists("theme_Publication")) {
    gp <- gp + theme_Publication() +
      ggplot2::theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.margin = ggplot2::margin()
      )
  }
  print(gp)
  grDevices::dev.off()
  
  invisible(list(dds_vst = dds_vst, p = p, npc_plotted = npc))
}

#' Write DESeq2 results merged with annotations (robust to list-like columns)
#'
#' Runs \code{DESeq2::results()}, merges the result with selected columns from
#' \code{rowData(dds)}, flattens any list-like columns so CSV export won't fail,
#' and writes a CSV. Also prints a compact summary of up/down counts.
#'
#' @param dds DESeq2::DESeqDataSet.
#' @param contrast Passed to \code{DESeq2::results()} (e.g., \code{c("condition","scramble_dicer","dcas_dicer")}).
#' @param alpha FDR cutoff used for the \code{padj} column (default 0.05).
#' @param annot_cols Character vector of \code{rowData(dds)} columns to merge.
#'        Missing names are ignored; duplicates are de-duplicated while preserving order.
#' @param file Output CSV path (default "de_results.csv").
#' @param collapse_sep Separator for collapsing list-like columns (default ";").
#' @param ... Additional arguments forwarded to \code{DESeq2::results()}.
#'
#' @return Invisibly returns the merged data.frame that was written.
#'
#' @examples
#' out <- write_de_results_csv(
#'   dds,
#'   contrast = c("condition","scramble_dicer","dcas_dicer"),
#'   annot_cols = c("gene_name","entrezid","gene_biotype","description"),
#'   file = "dicer_scramble_vs_dcas9.csv"
#' )
write_de_results_csv <- function(dds,
                                 contrast,
                                 alpha = 0.05,
                                 annot_cols = c("gene_name","entrezid","gene_biotype","description"),
                                 file = "de_results.csv",
                                 collapse_sep = ";",
                                 ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install 'DESeq2'.")
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Install 'SummarizedExperiment'.")
  }
  
  # 1) DE results
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha, ...)
  res_df <- as.data.frame(res, stringsAsFactors = FALSE)
  res_df$ensembl_gene_id <- rownames(res_df)
  
  # 2) Pick annotation columns safely from rowData and flatten list-like fields
  rd <- SummarizedExperiment::rowData(dds)
  rd_df <- as.data.frame(rd, stringsAsFactors = FALSE, check.names = FALSE)
  rd_df$ensembl_gene_id <- rownames(rd_df)
  
  # Use each requested column at most once, and only if present
  annot_cols <- unique(annot_cols)
  annot_cols <- annot_cols[annot_cols %in% colnames(rd_df)]
  ann <- rd_df[, c("ensembl_gene_id", annot_cols), drop = FALSE]
  
  # Flatten: CharacterList/IntegerList, list-cols, Rle, factors → character
  is_list_col <- vapply(ann, function(x) is.list(x) || inherits(x, "List"), logical(1))
  if (any(is_list_col)) {
    ann[is_list_col] <- lapply(ann[is_list_col], function(col) {
      vapply(col, function(el) paste(as.character(el), collapse = collapse_sep), character(1))
    })
  }
  ann[] <- lapply(ann, function(x) {
    if (inherits(x, "Rle")) x <- as.vector(x)
    if (is.factor(x)) x <- as.character(x)
    x
  })
  
  # 3) Merge and order by gene ID to keep stable rows
  out <- merge(ann, res_df, by = "ensembl_gene_id", sort = FALSE, all.y = TRUE)
  ord <- match(res_df$ensembl_gene_id, out$ensembl_gene_id)
  out <- out[ord, , drop = FALSE]
  
  # 4) Brief summary
  sig <- !is.na(out$padj) & out$padj < alpha
  up  <- sig & out$log2FoldChange > 0
  dn  <- sig & out$log2FoldChange < 0
  message(
    "DE summary (alpha=", alpha, "): ",
    sum(up), " up, ",
    sum(dn), " down; ",
    sum(sig), " significant out of ", nrow(out), "."
  )
  
  # 5) Write CSV (all columns now atomic)
  utils::write.csv(out, file = file, row.names = FALSE)
  invisible(out)
}

#' Differential expression report with capped labels on MA & volcano
#'
#' Runs DESeq2 `results()` for a given contrast, merges annotations, writes a CSV,
#' and saves a multi-page PDF with: (1) p-value histogram, (2) MA plot, (3) volcano,
#' (4) heatmap of significant genes. Labels on MA/volcano follow:
#'   - Only label genes that pass `padj < alpha` AND `|log2FC| >= fc_cutoff`
#'   - Label at most `max_labels_per_direction` for Up and for Down (by smallest padj)
#'   - Use `label_col` from rowData for label text (fallback to Ensembl ID)
#'
#' @param dds DESeq2::DESeqDataSet (already DESeq()-fit).
#' @param contrast length-3 vector for `results()`, e.g. c("condition","scramble_dicer","dcas_dicer").
#' @param alpha FDR threshold (default 0.05).
#' @param annot_cols rowData columns to include in the CSV merge.
#' @param results_csv CSV output filename.
#' @param plots_pdf PDF output filename (multi-page).
#' @param fc_cutoff absolute log2FC cutoff for significance labeling (default 2).
#' @param max_labels_per_direction cap per direction (Up/Down) for labels (default 10).
#' @param label_col rowData column to use for labels (default "gene_name").
#' @param top_annotation_cols optional sample-level columns for heatmap top annotation (default tries c("group","type")).
#' @param palette colors for MA/volcano points c(NS, Down, Up).
#' @param theme_publication apply `theme_Publication()` if available (default TRUE).
#' @return Invisibly: list(res_df=table, sig_ids=vector, contrast=contrast)
de_report <- function(
    dds,
    contrast,
    alpha = 0.05,
    annot_cols = c("gene_name","entrezid","gene_biotype","description"),
    results_csv = "de_results.csv",
    plots_pdf   = "de_plots.pdf",
    fc_cutoff   = 2,
    max_labels_per_direction = 10,
    label_col   = "gene_name",
    top_annotation_cols = c("group","type"),
    palette     = c("grey70","steelblue","firebrick"),
    theme_publication = TRUE
) {
  for (pkg in c("DESeq2","SummarizedExperiment","ggplot2","EnhancedVolcano",
                "ComplexHeatmap","grid","circlize","dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Install '%s'.", pkg))
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install 'ggrepel' for point labels.")
  
  # ---- results + annotations ----
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  res_df <- as.data.frame(res, stringsAsFactors = FALSE)
  res_df$ensembl_gene_id <- rownames(res_df)
  
  rd <- SummarizedExperiment::rowData(dds)
  rd_df <- as.data.frame(rd, stringsAsFactors = FALSE, check.names = FALSE)
  rd_df$ensembl_gene_id <- rownames(rd_df)
  cols <- unique(annot_cols[annot_cols %in% colnames(rd_df)])
  ann  <- rd_df[, c("ensembl_gene_id", cols), drop = FALSE]
  
  # flatten list-like rowData cols
  is_list_col <- vapply(ann, function(x) is.list(x) || inherits(x, "List"), logical(1))
  if (any(is_list_col)) {
    ann[is_list_col] <- lapply(ann[is_list_col], function(col) {
      vapply(col, function(el) paste(as.character(el), collapse = ";"), character(1))
    })
  }
  ann[] <- lapply(ann, function(x) { if (inherits(x, "Rle")) as.vector(x) else if (is.factor(x)) as.character(x) else x })
  
  merged <- merge(ann, res_df, by = "ensembl_gene_id", sort = FALSE, all.y = TRUE)
  merged <- merged[match(res_df$ensembl_gene_id, merged$ensembl_gene_id), , drop = FALSE]
  
  utils::write.csv(merged, file = results_csv, row.names = FALSE)
  
  # ---- labeling logic (shared) ----
  # choose label vector
  labvec <- if (label_col %in% names(merged)) merged[[label_col]] else merged$ensembl_gene_id
  # replace NA/blank with Ensembl ID
  labvec[is.na(labvec) | !nzchar(labvec)] <- merged$ensembl_gene_id[is.na(labvec) | !nzchar(labvec)]
  
  is_sig <- !is.na(merged$padj) & merged$padj < alpha
  is_up  <- merged$log2FoldChange >=  abs(fc_cutoff)
  is_dn  <- merged$log2FoldChange <= -abs(fc_cutoff)
  
  cand_up_idx <- which(is_sig & is_up)
  cand_dn_idx <- which(is_sig & is_dn)
  
  # order candidates by smallest padj
  up_ord <- cand_up_idx[order(merged$padj[cand_up_idx])]
  dn_ord <- cand_dn_idx[order(merged$padj[cand_dn_idx])]
  
  sel_up <- head(up_ord, max_labels_per_direction)
  sel_dn <- head(dn_ord, max_labels_per_direction)
  sel_idx <- c(sel_up, sel_dn)
  sel_labels <- labvec[sel_idx]
  
  # class for MA color
  merged <- dplyr::mutate(
    merged,
    sig_class = dplyr::case_when(
      is_sig & merged$log2FoldChange >  0 ~ "Upregulated",
      is_sig & merged$log2FoldChange <  0 ~ "Downregulated",
      TRUE ~ "NS"
    )
  )
  merged$sig_class <- factor(merged$sig_class, levels = c("NS","Downregulated","Upregulated"))
  
  # ---- heatmap prep ----
  cond_var <- contrast[1]; lev_A <- contrast[2]; lev_B <- contrast[3]
  main_title <- sprintf("%s: %s vs %s", cond_var, lev_A, lev_B)
  
  mat <- as.matrix(DESeq2::counts(dds, normalized = TRUE))
  mat_order <- mat[match(merged$ensembl_gene_id, rownames(mat)), , drop = FALSE]
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  keep_samples <- cd[[cond_var]] %in% c(lev_A, lev_B)
  mat_sub <- mat_order[, keep_samples, drop = FALSE]
  cd_sub  <- cd[keep_samples, , drop = FALSE]
  
  sig_rows <- which(is_sig)
  mat_sig <- mat_sub[sig_rows, , drop = FALSE]
  mat_sig <- mat_sig[!is.na(rownames(mat_sig)), , drop = FALSE]
  mat_sig_scaled <- if (nrow(mat_sig) > 0) t(scale(t(mat_sig))) else mat_sig
  
  ann_df <- cd_sub[, intersect(top_annotation_cols, colnames(cd_sub)), drop = FALSE]
  top_ann <- if (ncol(ann_df)) ComplexHeatmap::HeatmapAnnotation(df = ann_df) else NULL
  
  # ---- multi-page PDF ----
  grDevices::pdf(plots_pdf, width = 6, height = 6)
  
  # 1) p-value histogram
  graphics::hist(merged$pvalue,
                 main = paste(main_title),
                 xlab = "P-value", ylab = "Frequency",
                 col = "grey80", border = "grey40")
  
  # 2) MA plot with labels (ggrepel)
  gp_ma <- ggplot2::ggplot(merged, ggplot2::aes(x = log2(baseMean + 1), y = log2FoldChange, color = sig_class)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("NS" = palette[1], "Downregulated" = palette[2], "Upregulated" = palette[3])) +
    ggplot2::labs(title = main_title,
                  subtitle = sprintf("Labels: padj < %.2g & |Log2FC| ≥ %s",
                                     alpha, fc_cutoff),
                  x = "Average normalized log2(counts) + 1",
                  y = "log2 fold-change",
                  color = NULL)
  if (length(sel_idx)) {
    gp_ma <- gp_ma +
      ggrepel::geom_text_repel(
        data = merged[sel_idx, , drop = FALSE],
        ggplot2::aes(label = labvec[sel_idx]),
        size = 3, max.overlaps = Inf, min.segment.length = 0,
        box.padding = 0.25, point.padding = 0.25, segment.size = 0.3
      )
  }
  if (theme_publication && exists("theme_Publication")) gp_ma <- gp_ma + theme_Publication()
  print(gp_ma)
  
  # 3) Volcano (EnhancedVolcano) with capped labels
  keyvals <- rep(palette[1], nrow(merged)); names(keyvals) <- merged$ensembl_gene_id
  keyvals[which(is_sig & is_dn)] <- palette[2]
  keyvals[which(is_sig & is_up)] <- palette[3]
  
  # EnhancedVolcano labels by matching text in `lab` to `selectLab`
  lab_for_ev <- labvec
  select_lab <- unique(na.omit(lab_for_ev[sel_idx]))
  
  gp_vol <- EnhancedVolcano::EnhancedVolcano(
    merged,
    lab       = lab_for_ev,
    selectLab = select_lab,
    title     = main_title,
    subtitle  = sprintf("Labels: padj < %.2g & |Log2FC| ≥ %s", alpha, fc_cutoff),
    caption   = NULL,
    x         = "log2FoldChange",
    y         = "pvalue",
    pCutoffCol= "padj",
    pCutoff   = alpha,
    FCcutoff  = fc_cutoff,
    pointSize = 2.2,
    colAlpha  = 0.6,
    legendPosition = "none",
    labSize   = 3,
    labCol    = "black",
    labFace   = "bold",
    drawConnectors = TRUE,
    colConnectors  = "black",
    maxoverlapsConnectors = 25,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    colCustom = keyvals
  )
  print(gp_vol)
  
  # 4) Heatmap
  if (nrow(mat_sig_scaled) > 1 && ncol(mat_sig_scaled) > 1) {
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027"))
    ht <- ComplexHeatmap::Heatmap(
      mat_sig_scaled,
      name = "Z-score",
      column_title = paste0(main_title, "\nSignificant genes: padj < ", alpha),
      top_annotation = top_ann,
      border = TRUE,
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      row_dend_width = grid::unit(20, "mm")
    )
    draw(ht)
  } else {
    grid::grid.newpage()
    grid::grid.text("No significant genes to plot (padj < alpha).", gp = grid::gpar(cex = 1.2))
  }
  
  grDevices::dev.off()
  
  invisible(list(res_df = merged, sig_ids = merged$ensembl_gene_id[sig_rows <- which(is_sig)], contrast = contrast))
}

#' Test differences in gene biotype distributions (from annotated DESeq2 results)
#'
#' Given an annotated DE results table (e.g., output from your `de_report()`),
#' this function:
#'  1) classifies genes as Up/Down/NS using `padj < alpha` and `|log2FC| >= fc_cutoff`
#'  2) builds a contingency table: gene_biotype × {Down, NS, Up}
#'  3) runs a global chi-squared test of independence (falling back to Fisher's exact with simulation if needed)
#'  4) runs per-biotype enrichment tests:
#'       - Up vs (NS+Down)  (Fisher's exact)
#'       - Down vs (NS+Up)  (Fisher's exact)
#'     and BH-adjusts p-values
#'  5) writes a multi-page PDF of useful plots and a CSV of results
#'
#' Plots (multi-page PDF):
#'  • Stacked bar (counts) and stacked bar (proportions) by biotype × class
#'  • Lollipop of log2(OR) for Up-enrichment per biotype (with FDR)
#'  • Lollipop of log2(OR) for Down-enrichment per biotype (with FDR)
#'  • Heatmap of standardized residuals from the chi-squared test
#'
#' @param res_df Data frame with DE results and annotations (must include `padj`, `log2FoldChange`, and a biotype column).
#' @param biotype_col Column name of gene biotype in `res_df` (default: "gene_biotype").
#' @param padj_col Column name of adjusted p-values (default: "padj").
#' @param lfc_col Column name of log2 fold changes (default: "log2FoldChange").
#' @param alpha FDR threshold to call significance (default: 0.05).
#' @param fc_cutoff Absolute log2FC cutoff to call significance (default: 0).
#' @param drop_biotypes Optional character vector of biotype labels to drop (e.g., c(NA,"","unknown")).
#' @param min_biotype_count Drop biotypes with total count < this before testing/plotting (default: 10).
#' @param pdf_file Output PDF path for plots (default: "biotype_distribution_plots.pdf").
#' @param csv_file Output CSV path for per-biotype test results (default: "biotype_tests.csv").
#'
#' @return Invisibly, a list with:
#'   - counts: contingency table (biotype × class)
#'   - chisq: chisq.test (or fisher.test) object for the global test
#'   - up_enrichment: data frame of Up vs others Fisher tests per biotype
#'   - down_enrichment: data frame of Down vs others Fisher tests per biotype
#'
#' @examples
#' # res_df is the merged table from your de_report() (must include gene_biotype, padj, log2FoldChange)
#' out <- test_biotype_distribution(res_df, pdf_file = "biotype_plots.pdf", csv_file = "biotype_tests.csv")
test_biotype_distribution <- function(
    res_df,
    biotype_col = "gene_biotype",
    padj_col = "padj",
    lfc_col  = "log2FoldChange",
    alpha    = 0.05,
    fc_cutoff = 0,
    enrichment_fdr = 0.05,
    drop_biotypes = c(NA, "", "unknown", "Unknown"),
    min_biotype_count = 10,
    pdf_file = "biotype_distribution_plots.pdf",
    csv_file = "biotype_tests.csv",
    digits   = 3
) {
  for (pkg in c("dplyr","tidyr","ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Install '%s'.", pkg))
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Install 'ComplexHeatmap'.")
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Install 'circlize'.")
  
  need <- c(biotype_col, padj_col, lfc_col)
  if (!all(need %in% names(res_df))) stop("res_df must contain: ", paste(need, collapse = ", "))
  
  df <- res_df
  df[[biotype_col]] <- as.character(df[[biotype_col]])
  
  # drop unwanted biotypes
  if (length(drop_biotypes)) {
    drop_set <- setdiff(drop_biotypes, NA)
    df <- df[!(is.na(df[[biotype_col]]) & any(is.na(drop_biotypes))) &
               !(df[[biotype_col]] %in% drop_set), , drop = FALSE]
  }
  
  # call classes
  is_sig <- !is.na(df[[padj_col]]) & df[[padj_col]] < alpha &
    is.finite(df[[lfc_col]]) & abs(df[[lfc_col]]) >= abs(fc_cutoff)
  cls <- ifelse(is_sig & df[[lfc_col]] > 0, "Up",
                ifelse(is_sig & df[[lfc_col]] < 0, "Down", "NS"))
  df$.class <- factor(cls, levels = c("Down","NS","Up"))
  
  # keep biotypes with adequate counts
  keep_bio <- names(which(table(df[[biotype_col]]) >= min_biotype_count))
  df <- df[df[[biotype_col]] %in% keep_bio, , drop = FALSE]
  
  # contingency (handle empty case)
  if (nrow(df) == 0 || length(unique(df[[biotype_col]])) == 0) {
    utils::write.csv(
      data.frame(
        biotype = character(), n_down = integer(), n_ns = integer(), n_up = integer(),
        n_total = integer(), prop_down = numeric(), prop_ns = numeric(), prop_up = numeric(),
        up_log2OR = numeric(), up_ci_low = numeric(), up_ci_high = numeric(), up_p = numeric(), up_fdr = numeric(),
        down_log2OR = numeric(), down_ci_low = numeric(), down_ci_high = numeric(), down_p = numeric(), down_fdr = numeric(),
        global_test = character(), global_stat = numeric(), global_df = integer(), global_p = numeric()
      ),
      file = csv_file, row.names = FALSE
    )
    grDevices::pdf(pdf_file, width = 9, height = 5); on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage(); grid::grid.text("No biotypes pass 'min_biotype_count' filter.", gp = grid::gpar(cex = 1.2))
    return(invisible(NULL))
  }
  
  tab <- table(df[[biotype_col]], df$.class)
  
  # global test
  gl_try <- try(stats::chisq.test(tab), silent = TRUE)
  use_chisq <- !(inherits(gl_try, "try-error")) &&
    all(gl_try$expected >= 1) &&
    mean(gl_try$expected < 5) <= 0.2
  if (use_chisq) {
    gl <- gl_try
    gl_method <- gl$method
    gl_stat   <- unname(gl$statistic)
    gl_df     <- unname(gl$parameter)
    gl_p      <- unname(gl$p.value)
    stdres    <- gl$stdres
  } else {
    gl <- stats::fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
    gl_method <- "Fisher's exact (simulation)"
    gl_stat   <- NA
    gl_df     <- NA_integer_
    gl_p      <- unname(gl$p.value)
    stdres    <- NULL
  }
  
  # per-biotype tests
  bio_levels <- rownames(tab)
  fisher_block <- function(b, focus = c("Up","Down")) {
    focus <- match.arg(focus)
    if (focus == "Up") {
      a <- tab[b, "Up"];    b1 <- sum(tab[b, c("NS","Down")], na.rm = TRUE)
      c <- sum(tab[, "Up"], na.rm = TRUE) - a
      d <- sum(tab[, c("NS","Down")], na.rm = TRUE) - b1
    } else {
      a <- tab[b, "Down"];  b1 <- sum(tab[b, c("NS","Up")],   na.rm = TRUE)
      c <- sum(tab[, "Down"], na.rm = TRUE) - a
      d <- sum(tab[, c("NS","Up")],   na.rm = TRUE) - b1
    }
    ft <- stats::fisher.test(matrix(c(a,b1,c,d), nrow = 2))
    or <- unname(ifelse(is.null(ft$estimate), NA_real_, ft$estimate))
    ci <- if (!is.null(ft$conf.int)) unname(ft$conf.int) else c(NA_real_, NA_real_)
    data.frame(biotype = b, odds_ratio = or, ci_low = ci[1], ci_high = ci[2], p = ft$p.value,
               stringsAsFactors = FALSE)
  }
  up_tests   <- do.call(rbind, lapply(bio_levels, fisher_block, focus = "Up"))
  down_tests <- do.call(rbind, lapply(bio_levels, fisher_block, focus = "Down"))
  up_tests$log2OR   <- log2(up_tests$odds_ratio)
  down_tests$log2OR <- log2(down_tests$odds_ratio)
  up_tests$fdr   <- p.adjust(up_tests$p,   method = "BH")
  down_tests$fdr <- p.adjust(down_tests$p, method = "BH")
  
  # tidy table
  counts_df <- as.data.frame.matrix(tab)
  counts_df$biotype <- rownames(counts_df)
  counts_df <- dplyr::relocate(counts_df, biotype)
  props_df <- dplyr::mutate(counts_df,
                            n_total = Down + NS + Up,
                            prop_down = Down / n_total,
                            prop_ns   = NS   / n_total,
                            prop_up   = Up   / n_total)
  
  out_table <- props_df |>
    dplyr::rename(n_down = Down, n_ns = NS, n_up = Up) |>
    dplyr::left_join(dplyr::select(up_tests, biotype,
                                   up_log2OR = log2OR, up_ci_low = ci_low, up_ci_high = ci_high,
                                   up_p = p, up_fdr = fdr),
                     by = "biotype") |>
    dplyr::left_join(dplyr::select(down_tests, biotype,
                                   down_log2OR = log2OR, down_ci_low = ci_low, down_ci_high = ci_high,
                                   down_p = p, down_fdr = fdr),
                     by = "biotype") |>
    dplyr::mutate(global_test = gl_method,
                  global_stat = gl_stat, global_df = gl_df, global_p = gl_p) |>
    dplyr::arrange(dplyr::desc(n_total))
  
  num_cols <- c("prop_down","prop_ns","prop_up",
                "up_log2OR","up_ci_low","up_ci_high","up_p","up_fdr",
                "down_log2OR","down_ci_low","down_ci_high","down_p","down_fdr",
                "global_stat","global_p")
  out_table[num_cols] <- lapply(out_table[num_cols], function(x) if (is.numeric(x)) round(x, digits) else x)
  
  out_table <- out_table[, c(
    "biotype",
    "n_down","n_ns","n_up","n_total",
    "prop_down","prop_ns","prop_up",
    "up_log2OR","up_ci_low","up_ci_high","up_p","up_fdr",
    "down_log2OR","down_ci_low","down_ci_high","down_p","down_fdr",
    "global_test","global_stat","global_df","global_p"
  )]
  
  utils::write.csv(out_table, file = csv_file, row.names = FALSE)
  
  # ---------- PLOTS ----------
  grDevices::pdf(pdf_file, width = 9, height = 5); on.exit(grDevices::dev.off(), add = TRUE)
  
  long <- as.data.frame(tab) |>
    dplyr::rename(Biotype = Var1, Class = Var2, Count = Freq) |>
    dplyr::group_by(Biotype) |>
    dplyr::mutate(Prop = Count / sum(Count)) |>
    dplyr::ungroup()
  
  p_counts <- ggplot2::ggplot(long, ggplot2::aes(x = Biotype, y = Count, fill = Class)) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Gene biotype by DEG class (counts)",
                  subtitle = sprintf("Significant DEG: padj < %.2g & |log2FC| ≥ %s", alpha, fc_cutoff),
                  x = "Biotype", y = "Count") +
    ggplot2::scale_fill_manual(values = c("Down" = "#377eb8", "NS" = "grey70", "Up" = "#e41a1c")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::coord_flip()
  if (exists("theme_Publication")) p_counts <- p_counts + theme_Publication()
  print(p_counts)
  
  p_props <- ggplot2::ggplot(long, ggplot2::aes(x = Biotype, y = Prop, fill = Class)) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Gene biotype by DEG class (proportion)",
                  subtitle = sprintf("Significant DEG: padj < %.2g & |log2FC| ≥ %s", alpha, fc_cutoff),
                  x = "Biotype", y = "Proportion") +
    ggplot2::scale_fill_manual(values = c("Down" = "#377eb8", "NS" = "grey70", "Up" = "#e41a1c")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::coord_flip()
  if (exists("theme_Publication")) p_props <- p_props + theme_Publication()
  print(p_props)
  
  up_plot_df <- up_tests |>
    dplyr::mutate(sig = fdr < enrichment_fdr) |>
    dplyr::arrange(log2OR) |>
    dplyr::mutate(Biotype = factor(biotype, levels = biotype))
  
  p_up <- ggplot2::ggplot(up_plot_df, ggplot2::aes(x = log2OR, y = Biotype, color = sig)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = log2OR, y = Biotype, yend = Biotype)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "grey60")) +
    ggplot2::labs(title = "Upregulated gene biotype enrichment",
                  subtitle = sprintf("Fisher tests BH-FDR < %.2g", enrichment_fdr),
                  x = "log2(Odds Ratio)", y = NULL,
                  color = "Significant")
  if (exists("theme_Publication")) p_up <- p_up + theme_Publication()
  print(p_up)
  
  down_plot_df <- down_tests |>
    dplyr::mutate(sig = fdr < enrichment_fdr) |>
    dplyr::arrange(log2OR) |>
    dplyr::mutate(Biotype = factor(biotype, levels = biotype))
  
  p_down <- ggplot2::ggplot(down_plot_df, ggplot2::aes(x = log2OR, y = Biotype, color = sig)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = log2OR, y = Biotype, yend = Biotype)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("TRUE" = "#377eb8", "FALSE" = "grey60")) +
    ggplot2::labs(title = "Downregulated gene biotype enrichment",
                  subtitle = sprintf("Fisher tests BH-FDR < %.2g", enrichment_fdr),
                  x = "log2(Odds Ratio)", y = NULL,
                  color = "Significant")
  if (exists("theme_Publication")) p_down <- p_down + theme_Publication()
  print(p_down)
  
  if (!is.null(stdres)) {
    col_fun <- circlize::colorRamp2(c(-4, 0, 4), c("#4575b4", "white", "#d73027"))
    ht <- ComplexHeatmap::Heatmap(stdres,
                                  name = "StdResid", col = col_fun,
                                  column_title = "Chi-squared standardized residuals (biotype × class)",
                                  cluster_rows = TRUE, cluster_columns = FALSE)
    ComplexHeatmap::draw(ht)
  }
  
  invisible(list(
    counts = tab,
    global = list(method = gl_method, stat = gl_stat, df = gl_df, p = gl_p),
    up_enrichment = up_tests,
    down_enrichment = down_tests,
    table = out_table
  ))
}

#' Run GSEA (GO BP/MF/CC) from a DE result table and export CSV/PDF
#'
#' @param res_df Data frame with DE results; must contain columns for LFC and p-value.
#' @param lfc_col Column name for log2 fold-change. Default "log2FoldChange".
#' @param p_col Column name for p-value (unadjusted). Default "pvalue".
#' @param id_col Column with gene IDs matching \code{keyType} (e.g., ENTREZ). Default "entrezgene_id".
#' @param OrgDb Annotation package (e.g., \code{org.Mm.eg.db} or \code{org.Hs.eg.db}).
#' @param keyType Key type for \code{id_col} (e.g., "ENTREZID", "ENSEMBL"). Default "ENTREZID".
#' @param rank_fun Function to compute ranking score; takes numeric LFC and p, returns numeric vector.
#'   Default: \code{function(lfc, p) sign(lfc) * -log10(p)}.
#' @param file_prefix Filename prefix for exports. Default "contrast".
#' @param plot_title Title prefix for plots. Default "GSEA GO".
#' @param pvalueCutoff GSEA p-value cutoff. Default 0.1.
#' @param minGSSize,maxGSSize Gene set size bounds. Defaults 10 and 500.
#' @param simplify_go Logical; run \code{clusterProfiler::simplify()}. Default TRUE.
#' @param n_show Number of categories in dotplots. Default 10.
#' @param width,height PDF size in inches. Defaults 6, 5.
#' @details Ranking priority:
#'   1) \code{stat}; 2) \code{log2FoldChange/lfcSE}; 3) \code{sign(LFC)*-log10(pvalue)} as last resort.
#' @return List with elements \code{$bp, $mf, $cc} (each a \code{gseaResult}) and file paths.
#' @examples
#' # gsea_go_triplet(dicer_scramble_vs_dcas9, OrgDb=org.Mm.eg.db,
#' #                 id_col="entrezgene_id", file_prefix="dicer_scramble_vs_dcas9",
#' #                 plot_title="Male High-Fat vs Male Normal")
gsea_go_report <- function(
    res_df,
    id_col  = "entrezgene_id",
    OrgDb,
    keyType = "ENTREZID",
    file_prefix = "contrast",
    plot_title = "GSEA GO",
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    simplify_go = TRUE,
    n_show = 10,
    width = 6,
    height = 5
) {
  # --- deps ---
  for (pkg in c("clusterProfiler","enrichplot","ggplot2","viridis")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Install '%s'.", pkg))
  }
  if (missing(OrgDb)) stop("Provide OrgDb (e.g., org.Mm.eg.db or org.Hs.eg.db).")
  
  # --- choose ranking (prefer 'stat') ---
  has_stat  <- "stat" %in% colnames(res_df) && any(is.finite(res_df$stat))
  has_lfcse <- all(c("log2FoldChange","lfcSE") %in% colnames(res_df)) &&
    any(is.finite(res_df$log2FoldChange/res_df$lfcSE))
  has_p     <- all(c("log2FoldChange","pvalue") %in% colnames(res_df)) &&
    any(is.finite(res_df$log2FoldChange) & is.finite(res_df$pvalue))
  
  if (has_stat) {
    rank_vec <- as.numeric(res_df$stat)                            # signed test statistic
  } else if (has_lfcse) {
    rank_vec <- as.numeric(res_df$log2FoldChange / res_df$lfcSE)   # Wald stat recreation
  } else if (has_p) {
    p <- as.numeric(res_df$pvalue); p[p <= 0 | !is.finite(p)] <- NA
    rank_vec <- sign(res_df$log2FoldChange) * -log10(p)            # fallback (heuristic)
  } else {
    stop("No usable columns for ranking (need 'stat' or LFC/lfcSE or LFC+pvalue).")
  }
  
  ids <- as.character(res_df[[id_col]])
  keep <- is.finite(rank_vec) & !is.na(ids)
  ids <- ids[keep]; rank_vec <- rank_vec[keep]
  
  # Deduplicate IDs by highest absolute rank
  o <- order(abs(rank_vec), decreasing = TRUE)
  ids <- ids[o]; rank_vec <- rank_vec[o]
  keep1 <- !duplicated(ids)
  ids <- ids[keep1]; rank_vec <- rank_vec[keep1]
  
  names(rank_vec) <- ids
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  
  # --- small helper to run + simplify + export one ontology ---
  run_onto <- function(ont) {
    gse <- clusterProfiler::gseGO(
      geneList     = rank_vec,
      OrgDb        = OrgDb,
      ont          = ont,
      keyType      = keyType,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      verbose      = FALSE
    )
    if (simplify_go && inherits(gse, "gseaResult") && nrow(gse@result)) {
      # try simplify; if it fails (e.g., no similarity), keep original
      gse <- tryCatch(clusterProfiler::simplify(gse), error = function(e) gse)
    }
    
    # --- make gene IDs readable (ENTREZ/ENSEMBL -> SYMBOL) ---
    gse_readable <- tryCatch(
      clusterProfiler::setReadable(gse, OrgDb = OrgDb, keyType = keyType),
      error = function(e) gse
    )
    
    # Keep both: original IDs and SYMBOLs
    df_sym <- gse_readable@result
    df_id  <- gse@result
    if ("core_enrichment" %in% colnames(df_sym) && "core_enrichment" %in% colnames(df_id)) {
      df_sym$core_enrichment_SYMBOL <- df_sym$core_enrichment
      df_sym$core_enrichment_ID     <- df_id$core_enrichment
    }
    
    # export CSV
    csv_path <- file.path(sprintf("%s_go_%s.csv", file_prefix, tolower(ont)))
    utils::write.csv(df_sym, csv_path, row.names = FALSE)
    
    # export PDF dotplot
    pdf_path <- file.path(sprintf("%s_go_%s.pdf", file_prefix, tolower(ont)))
    grDevices::pdf(pdf_path, width = width, height = height)
    p <- enrichplot::dotplot(gse_readable, showCategory = n_show, x = "NES", font.size = 10, label_format = 50) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#444444", linewidth = 0.6) +
      viridis::scale_color_viridis(option = "viridis", direction = 1) +
      ggplot2::labs(title = sprintf("%s %s", plot_title, ont), fill = "qvalue")
    print(p)
    grDevices::dev.off()
    
    list(result = gse, csv = csv_path, pdf = pdf_path)
  }
  
  out_bp <- run_onto("BP")
  out_mf <- run_onto("MF")
  out_cc <- run_onto("CC")
  
  invisible(list(
    bp = out_bp$result, mf = out_mf$result, cc = out_cc$result,
    files = list(bp_csv = out_bp$csv, bp_pdf = out_bp$pdf,
                 mf_csv = out_mf$csv, mf_pdf = out_mf$pdf,
                 cc_csv = out_cc$csv, cc_pdf = out_cc$pdf),
    ranked_list = score
  ))
}

#' Plot a 2-set Euler diagram from gene vectors and save as PDF
#'
#' @param genes1 character. Genes in set 1.
#' @param genes2 character. Genes in set 2.
#' @param filename character. Output PDF path (e.g., "actb_arc_euler.pdf").
#' @param label1 character. Name for set 1.
#' @param label2 character. Name for set 2.
#' @param col1,col2 character. Hex fill colors.
#' @param alpha numeric. Fill transparency.
#' @param pct_mode "set" or "union". See details below.
#' @param width,height numeric. PDF size in inches.
#' @return Invisibly returns a list of counts and the fitted object.
plot_euler_2sets <- function(
    genes1, genes2, filename,
    label1 = "Set 1", label2 = "Set 2",
    col1 = "#63E6BE", col2 = "#FF41C3",
    alpha = 1,
    pct_mode = c("set","union"),
    width = 5, height = 5, margin_inch = 0
){
  stopifnot(is.character(genes1), is.character(genes2),
            is.character(filename), length(filename) == 1)
  pct_mode <- match.arg(pct_mode)
  
  A <- unique(genes1); B <- unique(genes2)
  Aonly <- length(setdiff(A, B))
  Bonly <- length(setdiff(B, A))
  AB    <- length(intersect(A, B))
  nA <- length(A); nB <- length(B); nUnion <- length(union(A, B))
  
  safe_div <- function(num, den) if (den == 0) 0 else num/den
  pct_str  <- function(p) sprintf("(%.0f%%)", 100 * p)
  
  if (pct_mode == "set") {
    lab_Aonly <- sprintf("%d\n%s", Aonly, pct_str(safe_div(Aonly, nA)))
    lab_Bonly <- sprintf("%d\n%s", Bonly, pct_str(safe_div(Bonly, nB)))
    lab_AB    <- sprintf("%d\n%s", AB,    pct_str(safe_div(AB, nUnion)))
  } else {
    lab_Aonly <- sprintf("%d\n%s", Aonly, pct_str(safe_div(Aonly, nUnion)))
    lab_Bonly <- sprintf("%d\n%s", Bonly, pct_str(safe_div(Bonly, nUnion)))
    lab_AB    <- sprintf("%d\n%s", AB,    pct_str(safe_div(AB,    nUnion)))
  }
  
  library(eulerr)
  fit <- euler(c(
    setNames(Aonly, label1),
    setNames(Bonly, label2),
    setNames(AB,    paste(label1, label2, sep = "&"))
  ))
  
  qlabs <- c(lab_Aonly, lab_Bonly, lab_AB)
  
  p <- plot(
    fit,
    fills  = list(fill = c(col1, col2), alpha = alpha),
    edges  = list(colour = "grey20", lwd = 1),
    labels = list(col = c(col1, col2), cex = 1.2, font = 2),
    quantities = list(labels = qlabs, cex = 1.0, font = 2, col = "black")
  )
  
  # add padding
  # p <- p + ggplot2::theme(plot.margin = ggplot2::unit(rep(margin_inch, 4), "in"))
  
  # Save to PDF
  ggplot2::ggsave(filename, plot = p, width = width, height = height,
                  device = grDevices::cairo_pdf)
  
  invisible(list(
    n_A_only = Aonly, n_B_only = Bonly, n_overlap = AB,
    n_A = nA, n_B = nB, n_union = nUnion, fit = fit, plot = p,
    file = normalizePath(filename, mustWork = FALSE)
  ))
}
