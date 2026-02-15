#!/usr/bin/env Rscript

# da_concordance.R - Multi-method differential accessibility concordance analysis
#
# Purpose: Compare DESeq2, edgeR, and limma-voom DA results to identify robust peaks
# Rationale: >=2/3 method concordance serves as a robustness/stability filter, not
#            statistical proof of differential accessibility
#
# Usage:
#   Rscript da_concordance.R --counts FILE --metadata FILE --deseq2 FILE --outdir DIR \
#     [--alpha 0.05] [--lfc 1] [--covariate COL]

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# Message function with timestamp
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n", sep = "")
}

# Check package availability
check_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed. Install with:\n  BiocManager::install('%s')", pkg, pkg))
  }
}

# Parse command-line arguments
option_list <- list(
  make_option(c("--counts"), type = "character", default = NULL,
              help = "featureCounts output file (TSV)", metavar = "FILE"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Sample metadata TSV (sample, condition, replicate)", metavar = "FILE"),
  make_option(c("--deseq2"), type = "character", default = NULL,
              help = "DESeq2 results CSV file", metavar = "FILE"),
  make_option(c("--outdir"), type = "character", default = NULL,
              help = "Output directory", metavar = "DIR"),
  make_option(c("--alpha"), type = "double", default = 0.05,
              help = "FDR significance threshold [default: %default]", metavar = "NUM"),
  make_option(c("--lfc"), type = "double", default = 1.0,
              help = "Log2 fold-change threshold [default: %default]", metavar = "NUM"),
  make_option(c("--covariate"), type = "character", default = NULL,
              help = "Optional covariate column from metadata", metavar = "COL")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$counts)) stop("--counts is required")
if (is.null(opt$metadata)) stop("--metadata is required")
if (is.null(opt$deseq2)) stop("--deseq2 is required")
if (is.null(opt$outdir)) stop("--outdir is required")

if (!file.exists(opt$counts)) stop(sprintf("Counts file not found: %s", opt$counts))
if (!file.exists(opt$metadata)) stop(sprintf("Metadata file not found: %s", opt$metadata))
if (!file.exists(opt$deseq2)) stop(sprintf("DESeq2 file not found: %s", opt$deseq2))

# Check required packages
check_package("edgeR")
check_package("limma")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

msg("Starting multi-method DA concordance analysis")
msg("Parameters:")
msg("  Counts:     %s", opt$counts)
msg("  Metadata:   %s", opt$metadata)
msg("  DESeq2:     %s", opt$deseq2)
msg("  Output dir: %s", opt$outdir)
msg("  Alpha:      %.3f", opt$alpha)
msg("  LFC:        %.2f", opt$lfc)
if (!is.null(opt$covariate)) {
  msg("  Covariate:  %s", opt$covariate)
}

# Load metadata
msg("Loading metadata...")
metadata <- read_tsv(opt$metadata, col_types = cols(), show_col_types = FALSE)

if (!all(c("sample", "condition") %in% colnames(metadata))) {
  stop("Metadata must contain 'sample' and 'condition' columns")
}

# Clean sample names in metadata
metadata <- metadata %>%
  mutate(sample_clean = gsub(".*/", "", sample),
         sample_clean = gsub("\\.final\\.bam$", "", sample_clean))

# Detect condition from sample names if needed
if (any(is.na(metadata$condition))) {
  metadata <- metadata %>%
    mutate(condition = case_when(
      grepl("ATF5WT", sample_clean, ignore.case = TRUE) ~ "WT",
      grepl("ATF5NULL", sample_clean, ignore.case = TRUE) ~ "KO",
      TRUE ~ condition
    ))
}

# Ensure condition is a factor with correct levels
metadata <- metadata %>%
  mutate(condition = factor(condition, levels = c("WT", "KO")))

msg("  Samples: %d", nrow(metadata))
msg("  Conditions: %s", paste(table(metadata$condition), collapse = " WT, ", sep = " KO"))

# Load counts
msg("Loading count matrix...")
counts_raw <- read_tsv(opt$counts, comment = "#", col_types = cols(), show_col_types = FALSE)

# Extract annotation columns (first 6) and count columns
annotation_cols <- counts_raw[, 1:6]
count_cols <- counts_raw[, -(1:6)]

# Clean column names to match metadata
colnames(count_cols) <- gsub(".*/", "", colnames(count_cols))
colnames(count_cols) <- gsub("\\.final\\.bam$", "", colnames(count_cols))

msg("  Peaks: %d", nrow(count_cols))
msg("  Samples in counts: %d", ncol(count_cols))

# Match samples between metadata and counts
common_samples <- intersect(metadata$sample_clean, colnames(count_cols))
if (length(common_samples) == 0) {
  stop("No matching samples between metadata and count matrix")
}

msg("  Matched samples: %d", length(common_samples))

# Subset and reorder
metadata <- metadata %>%
  filter(sample_clean %in% common_samples) %>%
  arrange(match(sample_clean, common_samples))

count_cols <- count_cols[, common_samples]

# Verify order
if (!all(colnames(count_cols) == metadata$sample_clean)) {
  stop("Sample order mismatch between counts and metadata")
}

# Apply low-count filter
msg("Applying low-count filter (rowSums >= 10)...")
keep <- rowSums(count_cols) >= 10
count_matrix <- as.matrix(count_cols[keep, ])
annotation_filtered <- annotation_cols[keep, ]

msg("  Peaks retained: %d (%.1f%%)", sum(keep), 100 * mean(keep))

# Prepare design matrix
msg("Preparing design matrix...")
if (!is.null(opt$covariate)) {
  if (!opt$covariate %in% colnames(metadata)) {
    stop(sprintf("Covariate '%s' not found in metadata", opt$covariate))
  }
  design_formula <- as.formula(sprintf("~ condition + %s", opt$covariate))
  msg("  Design: ~condition + %s", opt$covariate)
} else {
  design_formula <- ~ condition
  msg("  Design: ~condition")
}

design <- model.matrix(design_formula, data = metadata)

# ============================================================================
# edgeR analysis
# ============================================================================
msg("Running edgeR analysis...")

dge <- DGEList(counts = count_matrix, group = metadata$condition)
dge <- calcNormFactors(dge)

msg("  Estimating dispersions...")
dge <- estimateDisp(dge, design)

msg("  Fitting GLM...")
fit_edger <- glmQLFit(dge, design)

msg("  Performing QL F-test...")
qlf <- glmQLFTest(fit_edger, coef = "conditionKO")

results_edger <- as.data.frame(topTags(qlf, n = Inf)) %>%
  rownames_to_column("peak_id") %>%
  as_tibble() %>%
  rename(
    baseMean = logCPM,
    log2FoldChange = logFC,
    pvalue = PValue,
    padj = FDR
  ) %>%
  select(peak_id, baseMean, log2FoldChange, pvalue, padj)

msg("  edgeR complete: %d peaks tested", nrow(results_edger))
msg("  Significant (padj < %.3f): %d", opt$alpha, sum(results_edger$padj < opt$alpha, na.rm = TRUE))

# Write edgeR results
edger_out <- file.path(opt$outdir, "DA_edgeR.csv")
write_csv(results_edger, edger_out)
msg("  Saved: %s", edger_out)

# ============================================================================
# limma-voom analysis
# ============================================================================
msg("Running limma-voom analysis...")

v <- voom(dge, design, plot = FALSE)

msg("  Fitting linear model...")
fit_limma <- lmFit(v, design)
fit_limma <- eBayes(fit_limma)

msg("  Extracting results...")
results_limma <- topTable(fit_limma, coef = "conditionKO", number = Inf, sort.by = "none") %>%
  rownames_to_column("peak_id") %>%
  as_tibble() %>%
  rename(
    baseMean = AveExpr,
    log2FoldChange = logFC,
    pvalue = P.Value,
    padj = adj.P.Val
  ) %>%
  select(peak_id, baseMean, log2FoldChange, pvalue, padj)

msg("  limma-voom complete: %d peaks tested", nrow(results_limma))
msg("  Significant (padj < %.3f): %d", opt$alpha, sum(results_limma$padj < opt$alpha, na.rm = TRUE))

# Write limma results
limma_out <- file.path(opt$outdir, "DA_limma.csv")
write_csv(results_limma, limma_out)
msg("  Saved: %s", limma_out)

# ============================================================================
# Load DESeq2 results
# ============================================================================
msg("Loading DESeq2 results...")

results_deseq2 <- read_csv(opt$deseq2, col_types = cols(), show_col_types = FALSE)

# Extract unshrunken LFC if available, otherwise use shrunken
if ("log2FoldChange_unshrunken" %in% colnames(results_deseq2)) {
  msg("  Using unshrunken log2FoldChange for comparison")
  results_deseq2 <- results_deseq2 %>%
    rename(log2FoldChange = log2FoldChange_unshrunken)
} else if (!"log2FoldChange" %in% colnames(results_deseq2)) {
  stop("DESeq2 results must contain 'log2FoldChange' or 'log2FoldChange_unshrunken'")
}

# Ensure peak_id column exists
if (!"peak_id" %in% colnames(results_deseq2)) {
  if ("Geneid" %in% colnames(results_deseq2)) {
    results_deseq2 <- results_deseq2 %>% rename(peak_id = Geneid)
  } else {
    stop("DESeq2 results must contain 'peak_id' or 'Geneid' column")
  }
}

results_deseq2 <- results_deseq2 %>%
  select(peak_id, baseMean, log2FoldChange, pvalue, padj)

msg("  DESeq2 loaded: %d peaks", nrow(results_deseq2))
msg("  Significant (padj < %.3f): %d", opt$alpha, sum(results_deseq2$padj < opt$alpha, na.rm = TRUE))

# ============================================================================
# Merge results and calculate concordance
# ============================================================================
msg("Merging results from all three methods...")

merged <- results_deseq2 %>%
  rename_with(~ paste0("deseq2_", .), -peak_id) %>%
  full_join(
    results_edger %>% rename_with(~ paste0("edger_", .), -peak_id),
    by = "peak_id"
  ) %>%
  full_join(
    results_limma %>% rename_with(~ paste0("limma_", .), -peak_id),
    by = "peak_id"
  )

msg("  Total peaks in merged set: %d", nrow(merged))

# Calculate significance flags
msg("Calculating concordance metrics...")

merged <- merged %>%
  mutate(
    deseq2_sig = !is.na(deseq2_padj) & deseq2_padj < opt$alpha & abs(deseq2_log2FoldChange) >= opt$lfc,
    edger_sig = !is.na(edger_padj) & edger_padj < opt$alpha & abs(edger_log2FoldChange) >= opt$lfc,
    limma_sig = !is.na(limma_padj) & limma_padj < opt$alpha & abs(limma_log2FoldChange) >= opt$lfc,
    n_methods_sig = deseq2_sig + edger_sig + limma_sig,
    high_confidence = n_methods_sig >= 2
  )

# Write merged results with header comment
msg("Writing merged results...")
comparison_out <- file.path(opt$outdir, "DA_comparison.csv")

# Add framing comment
comment_line <- sprintf("# Multi-method DA comparison. n_methods_sig >= 2 indicates >=2/3 concordance, used as a robustness/stability filter (not statistical proof). Parameters: alpha=%.3f, lfc=%.2f", opt$alpha, opt$lfc)

con <- file(comparison_out, "w")
writeLines(comment_line, con)
close(con)

write_csv(merged, comparison_out, append = TRUE)
msg("  Saved: %s", comparison_out)

# Write high-confidence subset
high_conf <- merged %>%
  filter(high_confidence) %>%
  arrange(desc(n_methods_sig), deseq2_padj)

highconf_out <- file.path(opt$outdir, "DA_high_confidence.csv")
con <- file(highconf_out, "w")
writeLines(comment_line, con)
close(con)
write_csv(high_conf, highconf_out, append = TRUE)
msg("  Saved: %s (%d peaks)", highconf_out, nrow(high_conf))

# ============================================================================
# Generate summary statistics
# ============================================================================
msg("Generating concordance summary...")

summary_lines <- c(
  "Multi-method DA Concordance Summary",
  "====================================",
  sprintf("Analysis date: %s", Sys.time()),
  sprintf("Parameters: alpha=%.3f, lfc=%.2f", opt$alpha, opt$lfc),
  "",
  "Total peaks tested:",
  sprintf("  DESeq2:     %d", sum(!is.na(merged$deseq2_padj))),
  sprintf("  edgeR:      %d", sum(!is.na(merged$edger_padj))),
  sprintf("  limma-voom: %d", sum(!is.na(merged$limma_padj))),
  "",
  "Significant peaks per method:",
  sprintf("  DESeq2:     %d (%.1f%%)", sum(merged$deseq2_sig, na.rm = TRUE), 100 * mean(merged$deseq2_sig, na.rm = TRUE)),
  sprintf("  edgeR:      %d (%.1f%%)", sum(merged$edger_sig, na.rm = TRUE), 100 * mean(merged$edger_sig, na.rm = TRUE)),
  sprintf("  limma-voom: %d (%.1f%%)", sum(merged$limma_sig, na.rm = TRUE), 100 * mean(merged$limma_sig, na.rm = TRUE)),
  "",
  "Concordance breakdown:",
  sprintf("  3/3 methods: %d peaks", sum(merged$n_methods_sig == 3, na.rm = TRUE)),
  sprintf("  2/3 methods: %d peaks", sum(merged$n_methods_sig == 2, na.rm = TRUE)),
  sprintf("  1/3 methods: %d peaks", sum(merged$n_methods_sig == 1, na.rm = TRUE)),
  sprintf("  0/3 methods: %d peaks", sum(merged$n_methods_sig == 0, na.rm = TRUE)),
  "",
  sprintf("High-confidence peaks (>=2/3): %d (%.1f%% of union)",
          sum(merged$high_confidence, na.rm = TRUE),
          100 * mean(merged$high_confidence, na.rm = TRUE)),
  "",
  "Pairwise overlap (of significant peaks):",
  sprintf("  DESeq2 & edgeR:      %d", sum(merged$deseq2_sig & merged$edger_sig, na.rm = TRUE)),
  sprintf("  DESeq2 & limma-voom: %d", sum(merged$deseq2_sig & merged$limma_sig, na.rm = TRUE)),
  sprintf("  edgeR & limma-voom:  %d", sum(merged$edger_sig & merged$limma_sig, na.rm = TRUE)),
  "",
  "Pairwise Spearman correlation (unshrunken log2FC):",
  {
    complete_de <- !is.na(merged$deseq2_log2FoldChange) & !is.na(merged$edger_log2FoldChange)
    complete_dl <- !is.na(merged$deseq2_log2FoldChange) & !is.na(merged$limma_log2FoldChange)
    complete_el <- !is.na(merged$edger_log2FoldChange) & !is.na(merged$limma_log2FoldChange)
    rho_de <- if (sum(complete_de) > 2) cor(merged$deseq2_log2FoldChange[complete_de], merged$edger_log2FoldChange[complete_de], method = "spearman") else NA
    rho_dl <- if (sum(complete_dl) > 2) cor(merged$deseq2_log2FoldChange[complete_dl], merged$limma_log2FoldChange[complete_dl], method = "spearman") else NA
    rho_el <- if (sum(complete_el) > 2) cor(merged$edger_log2FoldChange[complete_el], merged$limma_log2FoldChange[complete_el], method = "spearman") else NA
    c(
      sprintf("  DESeq2 vs edgeR:      rho = %.4f (n = %d)", rho_de, sum(complete_de)),
      sprintf("  DESeq2 vs limma-voom: rho = %.4f (n = %d)", rho_dl, sum(complete_dl)),
      sprintf("  edgeR vs limma-voom:  rho = %.4f (n = %d)", rho_el, sum(complete_el))
    )
  },
  "",
  "Interpretation:",
  "  >= 2/3 concordance serves as a robustness/stability filter.",
  "  It does NOT constitute statistical proof of differential accessibility.",
  "  Use high-confidence peaks for downstream prioritization."
)

summary_out <- file.path(opt$outdir, "concordance_summary.txt")
writeLines(summary_lines, summary_out)
msg("  Saved: %s", summary_out)

# Print summary to console
cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n\n")

# ============================================================================
# Generate LFC scatter plots
# ============================================================================
msg("Generating LFC scatter plots...")

# Prepare data for scatter (significant peaks only for highlighting)
scatter_data <- merged %>%
  filter(!is.na(deseq2_log2FoldChange) & !is.na(edger_log2FoldChange) & !is.na(limma_log2FoldChange))

plot_scatter <- function(x_col, y_col, x_lab, y_lab, data) {
  ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = high_confidence), alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
    scale_color_manual(
      values = c("FALSE" = "gray70", "TRUE" = "firebrick"),
      labels = c("FALSE" = "< 2/3", "TRUE" = ">= 2/3"),
      name = "Concordance"
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      title = sprintf("%s vs %s", y_lab, x_lab)
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    coord_fixed()
}

p1 <- plot_scatter("deseq2_log2FoldChange", "edger_log2FoldChange",
                   "DESeq2 log2FC", "edgeR log2FC", scatter_data)
p2 <- plot_scatter("deseq2_log2FoldChange", "limma_log2FoldChange",
                   "DESeq2 log2FC", "limma-voom log2FC", scatter_data)
p3 <- plot_scatter("edger_log2FoldChange", "limma_log2FoldChange",
                   "edgeR log2FC", "limma-voom log2FC", scatter_data)

# Check for patchwork package
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

scatter_out <- file.path(opt$outdir, "method_lfc_scatter.pdf")

if (has_patchwork) {
  library(patchwork)
  combined <- (p1 | p2 | p3) +
    plot_annotation(
      title = "Multi-method log2 Fold-Change Comparison",
      subtitle = sprintf("High-confidence: >= 2/3 methods (alpha=%.3f, lfc=%.2f)", opt$alpha, opt$lfc)
    )
  ggsave(scatter_out, combined, width = 15, height = 5)
  msg("  Saved: %s (3-panel)", scatter_out)
} else {
  pdf(scatter_out, width = 15, height = 5)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  msg("  Saved: %s (3 separate pages)", scatter_out)
  msg("  Note: Install 'patchwork' for combined layout")
}

# ============================================================================
# Generate Venn diagram
# ============================================================================
msg("Generating Venn diagram...")

venn_data <- list(
  DESeq2 = merged %>% filter(deseq2_sig) %>% pull(peak_id),
  edgeR = merged %>% filter(edger_sig) %>% pull(peak_id),
  limma = merged %>% filter(limma_sig) %>% pull(peak_id)
)

has_venn <- requireNamespace("ggVennDiagram", quietly = TRUE)

venn_out <- file.path(opt$outdir, "method_venn.pdf")

if (has_venn) {
  library(ggVennDiagram)

  p_venn <- ggVennDiagram(venn_data, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = "Significant Peak Overlap Across Methods",
      subtitle = sprintf("alpha=%.3f, lfc=%.2f", opt$alpha, opt$lfc)
    ) +
    theme(legend.position = "none")

  ggsave(venn_out, p_venn, width = 7, height = 6)
  msg("  Saved: %s", venn_out)
} else {
  msg("  Skipping Venn diagram (install 'ggVennDiagram' for this feature)")
}

# ============================================================================
# Complete
# ============================================================================
msg("Concordance analysis complete!")
msg("Output directory: %s", opt$outdir)
msg("")
msg("Key files:")
msg("  - DA_comparison.csv       : Full merged results")
msg("  - DA_high_confidence.csv  : High-confidence peaks (>= 2/3 methods)")
msg("  - concordance_summary.txt : Summary statistics")
msg("  - method_lfc_scatter.pdf  : LFC comparison scatter plots")
if (has_venn) {
  msg("  - method_venn.pdf         : Venn diagram of significant peaks")
}
