# DESeq2 differential accessibility analysis for ATAC-seq
# Andres Chousal Cantu
# Icahn School of Medicine
#
# Usage:
#   Rscript deseq_atac.R <featureCounts.txt> <output_dir> [samples.tsv]
#   Rscript deseq_atac.R <deseq_results.csv> <output_dir> --figures-only [gene_map.tsv]
#
# Arguments:
#   1. featureCounts output file (from union peaks SAF) OR DESeq2 results CSV
#   2. Output directory for results
#   3. (Optional) Sample metadata TSV with columns: sample, condition, replicate
#      OR gene mapping TSV when --figures-only is used
#   4. (Optional) --figures-only flag to skip DESeq2 analysis and generate figures only

# =============================================================================
# DEPENDENCY CHECKS
# =============================================================================
check_package <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      stop(sprintf("[ERROR] Package '%s' required. Install with:\n  BiocManager::install('%s')",
                   pkg, pkg), call. = FALSE)
    } else {
      stop(sprintf("[ERROR] Package '%s' required. Install with:\n  install.packages('%s')",
                   pkg, pkg), call. = FALSE)
    }
  }
}

check_package("readr")
check_package("dplyr")
check_package("ggplot2")
check_package("stringr")
check_package("tibble")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript deseq_atac.R <featureCounts.txt> <output_dir> [samples.tsv]\n       Rscript deseq_atac.R <deseq_results.csv> <output_dir> --figures-only [gene_map.tsv]",
       call. = FALSE)
}

input_file <- args[1]
out_dir <- args[2]
figures_only <- "--figures-only" %in% args
gene_map_file <- if (figures_only && length(args) >= 4 && args[4] != "--figures-only") args[4] else NULL
meta_file <- if (!figures_only && length(args) >= 3 && args[3] != "--figures-only") args[3] else NULL

if (!file.exists(input_file)) {
  stop(sprintf("[ERROR] Input file not found: %s", input_file), call. = FALSE)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# FIGURES-ONLY MODE: Skip DESeq2, read results and generate figures
# =============================================================================
if (figures_only) {
  message("[DESeq2-Figures] Running in figures-only mode")
  message("[DESeq2-Figures] Reading DESeq2 results from: ", input_file)
  
  # Read DESeq2 results CSV
  da <- read_csv(input_file, show_col_types = FALSE)
  
  # Validate required columns
  required_cols <- c("peak_id", "log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(da))
  if (length(missing_cols) > 0) {
    stop(sprintf("[ERROR] DESeq2 results missing columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  
  message("[DESeq2-Figures] Columns in results: ", paste(colnames(da), collapse = ", "))
  
  # Try to load gene mapping for labels (peak -> gene name)
  gene_map <- NULL
  if (!is.null(gene_map_file) && file.exists(gene_map_file)) {
    message("[DESeq2-Figures] Loading gene mapping from: ", gene_map_file)
    gene_map <- tryCatch(
      read_tsv(gene_map_file, show_col_types = FALSE,
               col_names = c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist")),
      error = function(e) { message("[WARN] Could not read gene map: ", e$message); NULL }
    )
  }
  
  # Process results for visualization
  da <- da %>%
    mutate(
      sig = case_when(
        is.na(padj) | padj > 0.05 | abs(log2FoldChange) < 1 ~ "NS",
        log2FoldChange >= 1 ~ "Up",
        log2FoldChange <= -1 ~ "Down",
        TRUE ~ "NS"
      ),
      sig = factor(sig, levels = c("NS", "Up", "Down")),
      neg_log10_padj = -log10(pmax(padj, 1e-300))
    )
  
  # Add gene labels if mapping available
  if (!is.null(gene_map) && "gene_name" %in% colnames(gene_map)) {
    da <- da %>%
      left_join(gene_map %>% select(peak_id, gene_name) %>% distinct(),
                by = "peak_id")
    da$label <- da$gene_name
  }
  
  # Fallback: use peak ID if no gene name
  if (is.null(gene_map) || !"label" %in% colnames(da) || all(is.na(da$label))) {
    da$label <- sub("-[0-9]+$", "", da$peak_id)
  }
  
  # Select significant peaks to label
  n_sig_total <- sum(da$sig %in% c("Up", "Down"), na.rm = TRUE)
  label_cap <- 20
  da <- da %>% mutate(row_id = row_number())
  
  sig_for_labels <- da %>%
    filter(sig %in% c("Up", "Down"), !is.na(label), nzchar(label)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    group_by(sig, label) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  if (nrow(sig_for_labels) <= label_cap) {
    label_rows <- sig_for_labels$row_id
  } else {
    label_rows <- sig_for_labels %>%
      group_by(sig) %>%
      arrange(padj, desc(abs(log2FoldChange)), .by_group = TRUE) %>%
      slice_head(n = ceiling(label_cap / 2)) %>%
      ungroup() %>%
      arrange(padj, desc(abs(log2FoldChange))) %>%
      slice_head(n = label_cap) %>%
      pull(row_id)
  }
  
  da <- da %>%
    mutate(
      show_label = ifelse(row_id %in% label_rows, label, NA_character_),
      show_dot = !is.na(show_label)
    ) %>%
    select(-row_id)
  
  # Volcano plot
  message("[DESeq2-Figures] Generating volcano plot")
  p <- ggplot(da, aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
    geom_point(aes(alpha = sig), size = 1.05) +
    geom_point(
      data = da %>% filter(sig %in% c("Up", "Down")),
      size = 1.5, alpha = 0.85, show.legend = FALSE
    ) +
    geom_point(
      data = da %>% filter(show_dot),
      aes(fill = sig),
      shape = 21, size = 2.2, stroke = 0.35, color = "black",
      alpha = 0.95, show.legend = FALSE
    ) +
    scale_color_manual(values = c("NS" = "grey70", "Up" = "#D55E00", "Down" = "#0072B2")) +
    scale_alpha_manual(values = c("NS" = 0.14, "Up" = 0.70, "Down" = 0.70), guide = "none") +
    scale_fill_manual(values = c("NS" = "grey70", "Up" = "#D55E00", "Down" = "#0072B2"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    labs(title = "ATAC-seq DA (ATF5 KO vs WT)",
         subtitle = sprintf("Up: %d | Down: %d (padj<0.05, |LFC|>1); labels: %d",
                            sum(da$sig == "Up"), sum(da$sig == "Down"), sum(da$show_dot)),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Status") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")
  
  # Add labels with ggrepel if available
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = da %>% filter(!is.na(show_label)),
      aes(label = show_label),
      size = 2.8, fontface = "italic", max.overlaps = Inf, segment.color = "grey55",
      min.segment.length = 0.1, box.padding = 0.3, point.padding = 0.2,
      show.legend = FALSE
    )
  }
  
  ggsave(file.path(out_dir, "volcano_ATAC_DA.png"), p, width = 7, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "volcano_ATAC_DA.pdf"), p, width = 7, height = 7)
  message("[DESeq2-Figures] Volcano plot saved")
  
  # MA plot (if baseMean available)
  if ("baseMean" %in% colnames(da) && !all(is.na(da$baseMean))) {
    message("[DESeq2-Figures] Generating MA plot")
    p_ma <- ggplot(da, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
      geom_point(aes(alpha = sig), size = 1.05) +
      geom_point(
        data = da %>% filter(sig %in% c("Up", "Down")),
        size = 1.5, alpha = 0.85, show.legend = FALSE
      ) +
      geom_point(
        data = da %>% filter(show_dot),
        aes(fill = sig),
        shape = 21, size = 2.2, stroke = 0.35, color = "black",
        alpha = 0.95, show.legend = FALSE
      ) +
      scale_color_manual(values = c("NS" = "grey70", "Up" = "#D55E00", "Down" = "#0072B2")) +
      scale_alpha_manual(values = c("NS" = 0.14, "Up" = 0.70, "Down" = 0.70), guide = "none") +
      scale_fill_manual(values = c("NS" = "grey70", "Up" = "#D55E00", "Down" = "#0072B2"), guide = "none") +
      geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"), color = "grey50") +
      labs(title = "MA Plot (ATAC-seq DA)",
           subtitle = sprintf("Up: %d | Down: %d; labels: %d",
                              sum(da$sig == "Up"), sum(da$sig == "Down"), sum(da$show_dot)),
           x = "log10(mean counts + 1)", y = "log2 Fold Change", color = "Status") +
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom")
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_ma <- p_ma + ggrepel::geom_text_repel(
        data = da %>% filter(!is.na(show_label)),
        aes(label = show_label),
        size = 2.8, fontface = "italic", max.overlaps = Inf, segment.color = "grey55",
        min.segment.length = 0.1, box.padding = 0.3, point.padding = 0.2,
        show.legend = FALSE
      )
    }
    
    ggsave(file.path(out_dir, "MA_plot_ATAC_DA.png"), p_ma, width = 7, height = 6, dpi = 300)
    ggsave(file.path(out_dir, "MA_plot_ATAC_DA.pdf"), p_ma, width = 7, height = 6)
    message("[DESeq2-Figures] MA plot saved")
  }
  
  message("[DESeq2-Figures] Figures generation complete")
  quit(status = 0)
}

# =============================================================================
# FULL MODE: Read featureCounts and run DESeq2 analysis
# =============================================================================
check_package("DESeq2", bioc = TRUE)
check_package("apeglm", bioc = TRUE)

suppressPackageStartupMessages({
  library(DESeq2)
})

message("[DESeq2] Reading featureCounts file: ", input_file)
fc <- read.delim(input_file, comment.char = "#", check.names = FALSE)

# Diagnostic: check column count
message("[DESeq2] featureCounts columns: ", paste(colnames(fc), collapse = ", "))

peak_ids <- fc$Geneid
counts <- as.matrix(fc[, -(1:6)])
rownames(counts) <- peak_ids

# Verify counts matrix is numeric
if (!is.numeric(counts)) {
  warning(sprintf("[WARN] Counts matrix has mode '%s', attempting coercion to numeric", mode(counts)))
  counts <- apply(counts, 2, as.numeric)
  if (!is.numeric(counts)) {
    stop(sprintf("[ERROR] Could not coerce counts to numeric. Mode: %s. Check featureCounts output format.",
                 mode(counts)), call. = FALSE)
  }
}

# Verify we have samples
if (ncol(counts) == 0) {
  stop(sprintf("[ERROR] No sample columns found in featureCounts output. Total columns: %d (expected >= 7 for 1+ samples)",
               ncol(fc)), call. = FALSE)
}

# Clean sample names (remove path, keep basename without .final.bam)
samples <- colnames(counts)
samples_clean <- gsub(".*/", "", samples)
samples_clean <- gsub("\\.final\\.bam$", "", samples_clean)
colnames(counts) <- samples_clean
samples <- samples_clean

message("[DESeq2] Samples in count matrix: ", paste(samples, collapse = ", "))

# =============================================================================
# SAMPLE METADATA
# =============================================================================
# Try to read metadata file
use_metadata <- FALSE
if (!is.null(meta_file) && file.exists(meta_file)) {
  message("[DESeq2] Reading metadata file: ", meta_file)
  meta <- read_tsv(meta_file, show_col_types = FALSE)

  # Validate metadata is not empty
  if (nrow(meta) == 0) {
    warning("[WARN] Metadata file is empty, using pattern-based detection")
  } else {
    # Validate metadata has required columns
    required_cols <- c("sample", "condition")
    missing_cols <- setdiff(required_cols, colnames(meta))
    if (length(missing_cols) > 0) {
      stop(sprintf("[ERROR] Metadata missing required columns: %s",
                   paste(missing_cols, collapse = ", ")), call. = FALSE)
    }

    message(sprintf("[DESeq2] Metadata loaded: %d samples", nrow(meta)))
    use_metadata <- TRUE

    # Match samples in count matrix to metadata
    matched <- samples %in% meta$sample
    if (!all(matched)) {
      missing <- samples[!matched]
      warning(sprintf("[WARN] Samples not in metadata (using pattern fallback): %s",
                      paste(missing, collapse = ", ")))
    }
  }
} else {
  message("[DESeq2] No metadata file provided, using pattern-based detection")
}

# Build coldata from metadata if available, otherwise use patterns
coldata <- data.frame(sample = samples)

if (use_metadata) {
  # Match each sample to metadata
  coldata$condition <- sapply(samples, function(s) {
    m <- meta$condition[meta$sample == s]
    if (length(m) == 1) return(m)
    # Fallback to pattern
    if (grepl("ATF5WT", s, ignore.case = TRUE)) return("WT")
    return("KO")
  })
  coldata$replicate <- sapply(samples, function(s) {
    if ("replicate" %in% colnames(meta)) {
      m <- meta$replicate[meta$sample == s]
      if (length(m) == 1) return(m)
    }
    str_extract(s, "n[0-9]+")
  })
} else {
  # Pattern-based detection (fallback)
  coldata$condition <- sapply(samples, function(s) {
    if (grepl("ATF5WT", s, ignore.case = TRUE)) return("WT")
    if (grepl("ATF5NULL|ATF5KO|KO", s, ignore.case = TRUE)) return("KO")
    return("KO")
  })
  coldata$replicate <- str_extract(samples, "n[0-9]+")
}

# Convert to factors with explicit levels (WT vs KO)
coldata$condition <- factor(coldata$condition, levels = c("WT", "KO"))
coldata$replicate <- factor(coldata$replicate)
rownames(coldata) <- samples

message("[DESeq2] Sample metadata:")
print(coldata)

# Validate all samples have condition assigned
if (any(is.na(coldata$condition))) {
  stop("[ERROR] Some samples have NA condition - check metadata/patterns", call. = FALSE)
}

# =============================================================================
# DESeq2 ANALYSIS
# =============================================================================
message("[DESeq2] Creating DESeqDataSet")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~condition
)

# Filter low-count peaks
message("[DESeq2] Filtering low-count peaks (rowSum >= 10)")
n_before <- nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]
n_after <- nrow(dds)
message(sprintf("[DESeq2] Peaks: %d -> %d (removed %d)", n_before, n_after, n_before - n_after))

# Run DESeq2
message("[DESeq2] Running DESeq2")
dds <- DESeq(dds)

# LFC shrinkage with apeglm
message("[DESeq2] Shrinking log2 fold changes (apeglm)")
res <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm")

# =============================================================================
# SAVE RESULTS
# =============================================================================
res_tbl <- as.data.frame(res) %>%
  rownames_to_column("peak_id") %>%
  arrange(padj, desc(abs(log2FoldChange)))

out_file <- file.path(out_dir, "DA_results_DESeq2.csv")
write_csv(res_tbl, out_file)
message("[DESeq2] Full results saved to: ", out_file)

# Subsets: significant differential peaks
up <- res_tbl %>% filter(!is.na(padj), padj <= 0.05, log2FoldChange >= 1)
down <- res_tbl %>% filter(!is.na(padj), padj <= 0.05, log2FoldChange <= -1)

write_csv(up, file.path(out_dir, "DA_KO_up.csv"))
write_csv(down, file.path(out_dir, "DA_KO_down.csv"))

message(sprintf("[DESeq2] Differential peaks (padj<=0.05, |LFC|>=1): Up=%d, Down=%d",
                nrow(up), nrow(down)))

# =============================================================================
# PCA PLOT
# =============================================================================
message("[DESeq2] Generating PCA plot")
vsd <- vst(dds, blind = TRUE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$sample <- rownames(pcaData)

p <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  labs(title = "PCA of ATAC-seq samples",
       subtitle = "VST-transformed counts")

pca_file <- file.path(out_dir, "PCA_samples.pdf")
ggsave(
  filename = pca_file,
  plot = p,
  width = 6,
  height = 5
)
message("[DESeq2] PCA plot saved to: ", pca_file)

# Save PCA coordinates
pca_coords_file <- file.path(out_dir, "PCA_coordinates.csv")
pca_coords <- pcaData %>%
  select(sample, condition, PC1, PC2) %>%
  mutate(PC1_var = percentVar[1], PC2_var = percentVar[2])
write_csv(pca_coords, pca_coords_file)
message("[DESeq2] PCA coordinates saved to: ", pca_coords_file)

# =============================================================================
# SUMMARY
# =============================================================================
message("\n[DESeq2] Analysis complete")
message("  Results:    ", out_file)
message("  Up peaks:   ", file.path(out_dir, "DA_KO_up.csv"))
message("  Down peaks: ", file.path(out_dir, "DA_KO_down.csv"))
message("  PCA plot:   ", pca_file)
message("  PCA coords: ", pca_coords_file)
