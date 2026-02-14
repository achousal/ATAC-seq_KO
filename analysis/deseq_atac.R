# DESeq2 differential accessibility analysis for ATAC-seq
# Andres Chousal Cantu
# Icahn School of Medicine
#
# Usage:
#   Rscript deseq_atac.R <featureCounts.txt> <output_dir> [samples.tsv]
#
# Arguments:
#   1. featureCounts output file (from union peaks SAF)
#   2. Output directory for results
#   3. (Optional) Sample metadata TSV with columns: sample, condition, replicate

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

check_package("DESeq2", bioc = TRUE)
check_package("apeglm", bioc = TRUE)  # Required for lfcShrink type="apeglm"
check_package("readr")
check_package("dplyr")
check_package("ggplot2")
check_package("stringr")
check_package("tibble")

suppressPackageStartupMessages({
  library(DESeq2)
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
  stop("Usage: Rscript deseq_atac.R <featureCounts.txt> <output_dir> [samples.tsv]",
       call. = FALSE)
}

fc_file <- args[1]
out_dir <- args[2]
meta_file <- if (length(args) >= 3) args[3] else NULL

if (!file.exists(fc_file)) {
  stop(sprintf("[ERROR] featureCounts file not found: %s", fc_file), call. = FALSE)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# READ COUNT MATRIX
# =============================================================================
message("[DESeq2] Reading featureCounts file: ", fc_file)
fc <- read.delim(fc_file, comment.char = "#", check.names = FALSE)

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
    return("KO")  # default to KO if ambiguous
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
