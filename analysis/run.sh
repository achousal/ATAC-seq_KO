#!/usr/bin/env bash
# =============================================================================
# run.sh - ATAC-seq Pipeline (ATF5 KO Study)
# =============================================================================
# Single entry point for all pipeline stages.
#
# Usage:
#   bash run.sh [options]
#
# Stages (run in order unless skipped):
#   1. atac     - QC, alignment, filtering, peaks, DESeq2
#   2. integrate - Peak-to-gene mapping (requires --rna-deg and --gtf)
#   3. figures  - Volcano, MA, QC, heatmaps, integration panels
#
# Options:
#   --skip-atac         Skip ATAC-seq processing
#   --skip-integrate    Skip peak-to-gene integration
#   --skip-figures      Skip figure generation
#   --force             Re-run all steps even if outputs exist
#   --rna-deg FILE      Path to RNA-seq DEG CSV (for integration)
#   --gtf FILE          Path to GTF annotation (for integration)
#   --dry-run           Show what would run without executing
#   --validate          Run output validation only
#   -h, --help          Show this help
#
# Examples:
#   bash run.sh                          # Run full pipeline
#   bash run.sh --skip-atac              # Skip upstream, run figures only
#   bash run.sh --rna-deg degs.csv --gtf anno.gtf  # Enable integration
#   bash run.sh --validate               # Validate existing outputs
#
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
log()  { printf "[%(%F %T)T] %s\n" -1 "$*"; }
skip() { printf "[%(%F %T)T] [SKIP] %s\n" -1 "$*"; }
warn() { printf "[%(%F %T)T] [WARN] %s\n" -1 "$*" >&2; }
die()  { printf "[%(%F %T)T] [ERROR] %s\n" -1 "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

checkpoint() {
    local outfile="$1" desc="$2"
    if [[ -s "$outfile" ]]; then
        skip "$desc exists: $outfile"
        return 0
    fi
    return 1
}

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh" || die "Failed to load config.sh"

# Derived paths
DESEQ_DIR="${ATAC_BASE}/analysis/subread/deseq"
MACS_DIR="${ATAC_BASE}/analysis/macs"
DEEPTOOLS_DIR="${ATAC_BASE}/analysis/deeptools"
SAMTOOLS_DIR="${ATAC_BASE}/analysis/samtools"
INTEGRATION_DIR="${ATAC_BASE}/analysis/gene_integration"
RESULTS_DIR="${ATAC_BASE}/results"
FIGURES_DIR="${RESULTS_DIR}/figures"
TABLES_DIR="${RESULTS_DIR}/tables"
LOG_DIR="${ATAC_BASE}/analysis/logs"

# Defaults
SKIP_ATAC=0
SKIP_INTEGRATE=0
SKIP_FIGURES=0
FORCE=0
DRY_RUN=0
VALIDATE_ONLY=0
RNA_DEG_FILE=""
GTF_FILE=""

# -----------------------------------------------------------------------------
# Parse arguments
# -----------------------------------------------------------------------------
usage() { sed -n '/^# Usage:/,/^# ====/p' "$0" | grep -v '^# ===' | sed 's/^# //'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-atac)      SKIP_ATAC=1; shift ;;
        --skip-integrate) SKIP_INTEGRATE=1; shift ;;
        --skip-figures)   SKIP_FIGURES=1; shift ;;
        --force)          FORCE=1; shift ;;
        --rna-deg)        RNA_DEG_FILE="$2"; shift 2 ;;
        --gtf)            GTF_FILE="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=1; shift ;;
        --validate)       VALIDATE_ONLY=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)                die "Unknown argument: $1 (use --help)" ;;
    esac
done

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
mkdir -p "$LOG_DIR" "$RESULTS_DIR" "$FIGURES_DIR" "$TABLES_DIR" \
         "${FIGURES_DIR}/qc" "${FIGURES_DIR}/da" "${FIGURES_DIR}/integration" \
         "${FIGURES_DIR}/heatmaps"

log "=========================================="
log "ATAC-seq Pipeline"
log "=========================================="
log "ATAC_BASE:   $ATAC_BASE"
log "SKIP_ATAC:   $SKIP_ATAC"
log "SKIP_INTEGRATE: $SKIP_INTEGRATE"
log "SKIP_FIGURES: $SKIP_FIGURES"
log "FORCE:       $FORCE"
log "=========================================="

if [[ $DRY_RUN -eq 1 ]]; then
    log "[DRY-RUN] Would execute pipeline with above settings"
    exit 0
fi

# -----------------------------------------------------------------------------
# Validation-only mode
# -----------------------------------------------------------------------------
if [[ $VALIDATE_ONLY -eq 1 ]]; then
    log "Running validation only"
    bash "${SCRIPT_DIR}/validate_outputs.sh" "${ATAC_BASE}/analysis"
    exit $?
fi

# =============================================================================
# STAGE 1: ATAC-seq Processing
# =============================================================================
log "STAGE 1: ATAC-seq Processing"

if [[ $SKIP_ATAC -eq 1 ]]; then
    skip "ATAC-seq processing (--skip-atac)"
else
    DA_RESULTS="${DESEQ_DIR}/DA_results_DESeq2.csv"
    
    if [[ $FORCE -eq 0 ]] && checkpoint "$DA_RESULTS" "DESeq2 results"; then
        log "ATAC-seq outputs exist, skipping"
    else
        log "Running ATAC-seq pipeline"
        bash "${SCRIPT_DIR}/atacseq_analysis.sh" 2>&1 | tee "${LOG_DIR}/atac.$(date +%Y%m%d_%H%M%S).log"
    fi
fi

# Validate ATAC outputs
DA_RESULTS="${DESEQ_DIR}/DA_results_DESeq2.csv"
UNION_PEAKS="${MACS_DIR}/union_peaks.bed"
[[ -s "$DA_RESULTS" ]] || die "DA results not found: $DA_RESULTS"
[[ -s "$UNION_PEAKS" ]] || die "Union peaks not found: $UNION_PEAKS"
log "ATAC-seq outputs validated"

# =============================================================================
# STAGE 2: Peak-to-Gene Integration
# =============================================================================
log "STAGE 2: Peak-to-Gene Integration"

if [[ $SKIP_INTEGRATE -eq 1 ]]; then
    skip "Integration (--skip-integrate)"
elif [[ -z "$RNA_DEG_FILE" || -z "$GTF_FILE" ]]; then
    warn "Integration requires --rna-deg and --gtf; skipping"
    SKIP_INTEGRATE=1
else
    [[ -f "$RNA_DEG_FILE" ]] || die "RNA DEG file not found: $RNA_DEG_FILE"
    [[ -f "$GTF_FILE" ]] || die "GTF file not found: $GTF_FILE"
    
    mkdir -p "$INTEGRATION_DIR"
    MERGED_CSV="${INTEGRATION_DIR}/genelevel_RNA_Prom_Enh_merged.csv"
    
    # Step 2a: Peak-to-gene mapping
    if [[ $FORCE -eq 0 ]] && checkpoint "${INTEGRATION_DIR}/promoter_DA.csv" "Promoter DA"; then
        log "Gene integration prep already done"
    else
        log "Running gene integration prep"
        bash "${SCRIPT_DIR}/gene_integration_prep.sh" \
            --macs "$MACS_DIR" \
            --deseq "$DA_RESULTS" \
            --gtf "$GTF_FILE" \
            --outdir "$INTEGRATION_DIR" \
            2>&1 | tee "${LOG_DIR}/integration_prep.$(date +%Y%m%d_%H%M%S).log"
    fi
    
    # Step 2b: RNA+ATAC merge and figures
    if [[ $FORCE -eq 0 ]] && checkpoint "$MERGED_CSV" "Gene-level merge"; then
        log "RNA+ATAC merge already done"
    else
        if have Rscript && [[ -f "${SCRIPT_DIR}/rna_atac_figures.R" ]]; then
            have module && module load R/4.2.0 2>/dev/null || true
            
            log "Running RNA+ATAC integration"
            Rscript "${SCRIPT_DIR}/rna_atac_figures.R" \
                --rna_deg_file "$RNA_DEG_FILE" \
                --datadir "$INTEGRATION_DIR" \
                --figdir "${FIGURES_DIR}/integration" \
                2>&1 | tee "${LOG_DIR}/rna_atac.$(date +%Y%m%d_%H%M%S).log"

            # Run exploratory figure panel
            if [[ -f "${SCRIPT_DIR}/rna_atac_exploratory.R" ]]; then
                log "Running exploratory RNA+ATAC figures"
                Rscript "${SCRIPT_DIR}/rna_atac_exploratory.R" \
                    --rna_deg_file "$RNA_DEG_FILE" \
                    --datadir "$INTEGRATION_DIR" \
                    --figdir "${FIGURES_DIR}/integration" \
                    2>&1 | tee "${LOG_DIR}/rna_atac_exploratory.$(date +%Y%m%d_%H%M%S).log"
            fi

            # Copy tables
            for f in "$INTEGRATION_DIR"/*.csv; do
                [[ -f "$f" ]] && cp "$f" "$TABLES_DIR/"
            done
        else
            warn "rna_atac_figures.R not found or Rscript unavailable"
        fi
    fi
fi

# =============================================================================
# STAGE 3: Figures
# =============================================================================
log "STAGE 3: Figure Generation"

if [[ $SKIP_FIGURES -eq 1 ]]; then
    skip "Figure generation (--skip-figures)"
else
    # Load R if on HPC
    have module && module load R/4.2.0 2>/dev/null || true
    
    # -------------------------------------------------------------------------
    # 3a: ATAC DA visualizations
    # -------------------------------------------------------------------------
    VOLCANO_PNG="${FIGURES_DIR}/da/volcano_ATAC_DA.png"
    
    if [[ $FORCE -eq 1 ]] || ! checkpoint "$VOLCANO_PNG" "Volcano plot"; then
        log "Generating DA figures"
        GENE_MAP="${INTEGRATION_DIR}/promoter_closest.tsv"
        Rscript "${SCRIPT_DIR}/deseq_atac.R" "$DA_RESULTS" "${FIGURES_DIR}/da" --figures-only "$GENE_MAP" \
            2>&1 || {
            # Fallback: generate figures inline if deseq_atac.R fails
            Rscript --vanilla - "$DA_RESULTS" "${FIGURES_DIR}/da" "$GENE_MAP" <<'EOF'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    message("ggrepel not installed - labels will overlap. Install with: install.packages('ggrepel')")
  }
})
args <- commandArgs(trailingOnly = TRUE)
da <- read_csv(args[1], show_col_types = FALSE)
outdir <- args[2]; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
gene_map_file <- if (length(args) >= 3) args[3] else NULL

# Standardize column names
if (!"log2FoldChange" %in% names(da) && "log2FC" %in% names(da)) da <- rename(da, log2FoldChange = log2FC)
if (!"padj" %in% names(da) && "FDR" %in% names(da)) da <- rename(da, padj = FDR)

# Try to load gene mapping for labels (peak -> gene name)
# File format: TSV with columns peak_id, gene_id, gene_name, strand, dist_bp, abs_dist
gene_map <- NULL
if (!is.null(gene_map_file) && file.exists(gene_map_file)) {
  message("Loading gene mapping from: ", gene_map_file)
  gene_map <- tryCatch(
    read_tsv(gene_map_file, show_col_types = FALSE,
             col_names = c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist")),
    error = function(e) { message("Could not read gene map: ", e$message); NULL }
  )
}

# Compute significance and labeling columns
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

# Add gene labels if mapping available, otherwise use peak IDs
if (!is.null(gene_map) && "gene_name" %in% names(gene_map)) {
  # Join on peak_id
  peak_col <- intersect(c("peak_id", "Geneid", "peak"), names(gene_map))[1]
  da_col <- intersect(c("peak_id", "Geneid", "peak"), names(da))[1]
  if (!is.na(peak_col) && !is.na(da_col)) {
    da <- da %>%
      left_join(gene_map %>% select(all_of(c(peak_col, "gene_name"))) %>% distinct(),
                by = setNames(peak_col, da_col))
    da$label <- da$gene_name
  }
}

# Fallback: use truncated peak ID if no gene name
if (!"label" %in% names(da) || all(is.na(da$label))) {
  da_col <- intersect(c("peak_id", "Geneid", "peak"), names(da))[1]
  if (!is.na(da_col)) {
    # Shorten peak IDs: chr1:1000-2000 -> chr1:1000
    da$label <- sub("-[0-9]+$", "", da[[da_col]])
  } else {
    da$label <- NA_character_
  }
}

# Select significant peaks to label:
# label all if few; otherwise keep strongest significant hits (adaptive cap).
n_sig_total <- sum(da$sig %in% c("Up", "Down"), na.rm = TRUE)
label_cap <- min(max(10, ceiling(0.25 * n_sig_total)), 20)
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

# Build volcano plot
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
ggsave(file.path(outdir, "volcano_ATAC_DA.png"), p, width = 7, height = 7, dpi = 300)
ggsave(file.path(outdir, "volcano_ATAC_DA.pdf"), p, width = 7, height = 7)
message("Volcano plot saved")

# MA plot
if ("baseMean" %in% names(da) && !all(is.na(da$baseMean))) {
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
  ggsave(file.path(outdir, "MA_plot_ATAC_DA.png"), p_ma, width = 7, height = 6, dpi = 300)
  ggsave(file.path(outdir, "MA_plot_ATAC_DA.pdf"), p_ma, width = 7, height = 6)
  message("MA plot saved")
}
message("Figures saved to: ", outdir)
EOF
        }
    fi
    
    # -------------------------------------------------------------------------
    # 3b: Copy PCA plot
    # -------------------------------------------------------------------------
    PCA_SRC="${DESEQ_DIR}/PCA_samples.pdf"
    PCA_DST="${FIGURES_DIR}/qc/PCA_samples.pdf"
    [[ -f "$PCA_SRC" && ! -f "$PCA_DST" ]] && cp "$PCA_SRC" "$PCA_DST"
    
    # -------------------------------------------------------------------------
    # 3c: deepTools heatmaps (if available)
    # -------------------------------------------------------------------------
    if have computeMatrix && have plotHeatmap; then
        HEATMAP="${FIGURES_DIR}/heatmaps/heatmap_all_peaks.pdf"
        BW_FILES=("${DEEPTOOLS_DIR}"/*.CPM.bw)
        
        if [[ ${#BW_FILES[@]} -gt 0 && -f "${BW_FILES[0]}" ]]; then
            if [[ $FORCE -eq 1 ]] || ! checkpoint "$HEATMAP" "Heatmap"; then
                log "Generating deepTools heatmap"
                MATRIX="${DEEPTOOLS_DIR}/matrix_all.gz"
                LABELS=$(printf "%s " "${BW_FILES[@]}" | sed 's|.*/||g; s|\.CPM\.bw||g')
                
                computeMatrix reference-point -R "$UNION_PEAKS" -S "${BW_FILES[@]}" \
                    -o "$MATRIX" -b 3000 -a 3000 --referencePoint center \
                    -p "${ATAC_THREADS:-8}" --skipZeros --samplesLabel $LABELS 2>/dev/null || true
                
                [[ -s "$MATRIX" ]] && plotHeatmap -m "$MATRIX" -out "$HEATMAP" \
                    --dpi 300 --colorMap Reds --heatmapHeight 12 2>/dev/null || true
            fi
        fi
    fi
    
    # -------------------------------------------------------------------------
    # 3d: Copy tables to results
    # -------------------------------------------------------------------------
    log "Copying tables to results"
    [[ -f "$DA_RESULTS" ]] && cp "$DA_RESULTS" "$TABLES_DIR/"
    [[ -f "${DESEQ_DIR}/DA_KO_up.csv" ]] && cp "${DESEQ_DIR}/DA_KO_up.csv" "$TABLES_DIR/"
    [[ -f "${DESEQ_DIR}/DA_KO_down.csv" ]] && cp "${DESEQ_DIR}/DA_KO_down.csv" "$TABLES_DIR/"
    [[ -f "${DESEQ_DIR}/PCA_coordinates.csv" ]] && cp "${DESEQ_DIR}/PCA_coordinates.csv" "$TABLES_DIR/"
fi

# =============================================================================
# Summary
# =============================================================================
log "=========================================="
log "Pipeline Complete"
log "=========================================="
log "Results: $RESULTS_DIR"
log ""
log "Key outputs:"
log "  DESeq2:     $DA_RESULTS"
log "  Peaks:      $UNION_PEAKS"
log "  MultiQC:    ${ATAC_BASE}/analysis/multiqc/multiqc_report.html"
[[ -f "${INTEGRATION_DIR}/genelevel_RNA_Prom_Enh_merged.csv" ]] && \
    log "  Integration: ${INTEGRATION_DIR}/genelevel_RNA_Prom_Enh_merged.csv"
log ""
log "Figures: $FIGURES_DIR"
log "Tables:  $TABLES_DIR"
log "=========================================="
