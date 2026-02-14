#!/usr/bin/env Rscript
# RNA vs ATAC Exploratory Figure Panel
# Generates comprehensive exploratory visualizations for RNA/Promoter/Enhancer integration
#
# Usage:
#   Rscript rna_atac_exploratory.R \
#     --rna_deg_file /path/to/degs.csv \
#     --datadir /path/to/gene_integration \
#     --figdir /path/to/output/figures

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
  library(patchwork)
  library(pheatmap)
  library(UpSetR)
  library(ggVennDiagram)
})

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n", sep = "")

# ============================================================================
# Data Loading (reused from rna_atac_figures.R)
# ============================================================================

clean_ens <- function(x) {
  x <- as.character(x)
  gsub("\\.\\d+$", "", x)
}

collapse_peak_to_gene <- function(da_df, map_df, alpha = 0.05) {
  da2 <- da_df %>%
    mutate(
      peak_id = as.character(peak_id),
      padj = as.numeric(padj),
      log2FoldChange = as.numeric(log2FoldChange)
    )

  map2 <- map_df %>%
    transmute(
      peak_id = as.character(peak_id),
      gene_id_clean = clean_ens(gene_id),
      gene_name = as.character(gene_name)
    )

  joined <- inner_join(map2, da2, by = "peak_id")
  if (nrow(joined) == 0) {
    return(tibble(
      gene_id_clean = character(),
      n_peaks = integer(),
      n_sig = integer(),
      best_peak_id = character(),
      best_log2FC = double(),
      best_padj = double(),
      mean_log2FC = double(),
      median_log2FC = double()
    ))
  }

  best_rows <- joined %>%
    group_by(gene_id_clean) %>%
    arrange(is.na(padj), padj, desc(abs(log2FoldChange))) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene_id_clean,
      best_peak_id = peak_id,
      best_log2FC = log2FoldChange,
      best_padj = padj
    )

  summ <- joined %>%
    group_by(gene_id_clean) %>%
    summarise(
      n_peaks = n(),
      n_sig = sum(!is.na(padj) & padj < alpha),
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      median_log2FC = median(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    )

  left_join(summ, best_rows, by = "gene_id_clean")
}

load_and_merge_data <- function(opt) {
  # Load mapping files
  prom_map <- read_tsv(
    file.path(opt$datadir, "promoter_closest.tsv"),
    col_names = c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist"),
    show_col_types = FALSE
  )
  enh_map <- read_tsv(
    file.path(opt$datadir, "enhancer_closest.tsv"),
    col_names = c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist"),
    show_col_types = FALSE
  )

  # Load DA results
  prom_da <- read_csv(file.path(opt$datadir, "promoter_DA.csv"), show_col_types = FALSE)
  enh_da <- read_csv(file.path(opt$datadir, "enhancer_DA.csv"), show_col_types = FALSE)

  # Collapse to gene level
  prom_gene <- collapse_peak_to_gene(prom_da, prom_map, alpha = opt$alpha) %>%
    rename_with(~paste0("prom_", .x), -gene_id_clean)
  enh_gene <- collapse_peak_to_gene(enh_da, enh_map, alpha = opt$alpha) %>%
    rename_with(~paste0("enh_", .x), -gene_id_clean)

  # Load RNA DEGs
  rna <- read_csv(opt$rna_deg_file, show_col_types = FALSE)

  gene_id_col <- intersect(names(rna), c("gene_id", "GeneID", "ensembl_gene_id", "id"))
  gene_name_col <- intersect(names(rna), c("gene_name", "symbol", "GeneName", "external_gene_name", "gene"))
  lfc_col <- intersect(names(rna), c("log2FoldChange", "log2FC", "lfc"))
  padj_col <- intersect(names(rna), c("padj", "adj_pval", "FDR", "qval", "q_value"))

  if (length(gene_id_col) == 0) stop("RNA DEGs missing gene_id-like column.")
  if (length(lfc_col) == 0) stop("RNA DEGs missing log2FoldChange-like column.")
  if (length(padj_col) == 0) stop("RNA DEGs missing padj/FDR-like column.")

  gene_id_col <- gene_id_col[1]
  gene_name_col <- if (length(gene_name_col) > 0) gene_name_col[1] else NULL
  lfc_col <- lfc_col[1]
  padj_col <- padj_col[1]

  rna2 <- rna %>%
    transmute(
      gene_id_clean = clean_ens(.data[[gene_id_col]]),
      gene_label = if (!is.null(gene_name_col)) as.character(.data[[gene_name_col]]) else clean_ens(.data[[gene_id_col]]),
      rna_log2FC = as.numeric(.data[[lfc_col]]),
      rna_padj = as.numeric(.data[[padj_col]])
    ) %>%
    distinct(gene_id_clean, .keep_all = TRUE)

  # Merge
  gene_merged <- rna2 %>%
    left_join(prom_gene, by = "gene_id_clean") %>%
    left_join(enh_gene, by = "gene_id_clean") %>%
    mutate(
      rna_sig = !is.na(rna_padj) & rna_padj < opt$alpha,
      prom_sig = !is.na(prom_best_padj) & prom_best_padj < opt$alpha,
      enh_sig = !is.na(enh_best_padj) & enh_best_padj < opt$alpha,
      any_sig = rna_sig | prom_sig | enh_sig,
      sig_class = case_when(
        rna_sig & (prom_sig | enh_sig) ~ "RNA + ATAC",
        rna_sig ~ "RNA only",
        (prom_sig | enh_sig) ~ "ATAC only",
        TRUE ~ "not sig"
      ),
      # Direction categories
      rna_dir = case_when(
        !rna_sig ~ "NS",
        rna_log2FC > 0 ~ "Up",
        rna_log2FC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      prom_dir = case_when(
        !prom_sig ~ "NS",
        prom_best_log2FC > 0 ~ "Up",
        prom_best_log2FC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      enh_dir = case_when(
        !enh_sig ~ "NS",
        enh_best_log2FC > 0 ~ "Up",
        enh_best_log2FC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )

  gene_merged
}

# ============================================================================
# Figure 1: Hexbin Density Plots
# ============================================================================

make_hexbin_panel <- function(df, x, y, xlab, ylab, title) {
  xvals <- df[[x]]
  yvals <- df[[y]]
  ok <- is.finite(xvals) & is.finite(yvals)
  rho <- suppressWarnings(cor(xvals[ok], yvals[ok], method = "spearman"))
  n_ok <- sum(ok)

  df_filt <- df[ok, ]

  p <- ggplot(df_filt, aes(x = .data[[x]], y = .data[[y]])) +
    geom_hex(bins = 40) +
    scale_fill_viridis_c(option = "plasma", trans = "log10", name = "Count") +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    geom_smooth(method = "lm", color = "white", linewidth = 0.8, se = FALSE, linetype = "solid") +
    labs(
      title = title,
      subtitle = sprintf("Spearman rho = %.3f (n = %d)", rho, n_ok),
      x = xlab, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(color = "grey40", size = 9),
      legend.position = "right"
    ) +
    coord_cartesian(expand = TRUE)

  p
}

fig1_hexbin <- function(df, figdir) {
  msg("Figure 1: Hexbin density plots")

  df_rna_prom <- df %>% filter(!is.na(rna_log2FC), !is.na(prom_best_log2FC))
  df_rna_enh <- df %>% filter(!is.na(rna_log2FC), !is.na(enh_best_log2FC))
  df_prom_enh <- df %>% filter(!is.na(prom_best_log2FC), !is.na(enh_best_log2FC))

  p1 <- make_hexbin_panel(df_rna_prom, "rna_log2FC", "prom_best_log2FC",
                          "RNA log2FC", "Promoter log2FC", "RNA vs Promoter")
  p2 <- make_hexbin_panel(df_rna_enh, "rna_log2FC", "enh_best_log2FC",
                          "RNA log2FC", "Enhancer log2FC", "RNA vs Enhancer")
  p3 <- make_hexbin_panel(df_prom_enh, "prom_best_log2FC", "enh_best_log2FC",
                          "Promoter log2FC", "Enhancer log2FC", "Promoter vs Enhancer")

  p_combined <- (p1 | p2 | p3) +
    plot_annotation(title = "Hexbin Density: RNA vs ATAC Fold Changes",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))

  ggsave(file.path(figdir, "Fig1_hexbin_density.png"), p_combined,
         width = 14, height = 5, dpi = 300)

  list(p1 = p1, p2 = p2, p3 = p3)
}

# ============================================================================
# Figure 2: Quadrant Concordance Plots
# ============================================================================

make_quadrant_panel <- function(df, x, y, xlab, ylab, title) {
  xvals <- df[[x]]
  yvals <- df[[y]]
  ok <- is.finite(xvals) & is.finite(yvals)
  df_filt <- df[ok, ]

  # Quadrant counts
  q_uu <- sum(df_filt[[x]] > 0 & df_filt[[y]] > 0)
  q_ud <- sum(df_filt[[x]] > 0 & df_filt[[y]] < 0)
  q_du <- sum(df_filt[[x]] < 0 & df_filt[[y]] > 0)
  q_dd <- sum(df_filt[[x]] < 0 & df_filt[[y]] < 0)

  # Concordant vs discordant
  concordant <- q_uu + q_dd
  discordant <- q_ud + q_du
  total <- concordant + discordant

  # Fisher test for enrichment
  mat <- matrix(c(q_uu, q_ud, q_du, q_dd), nrow = 2)
  fisher_p <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA)

  # Axis limits (symmetric around 0)
  max_abs_x <- max(abs(df_filt[[x]]), na.rm = TRUE)
  max_abs_y <- max(abs(df_filt[[y]]), na.rm = TRUE)

  # Colors by sig_class
  pal <- c(
    "not sig" = "grey70",
    "ATAC only" = "#D55E00",
    "RNA only" = "#CC79A7",
    "RNA + ATAC" = "#0072B2"
  )

  p <- ggplot(df_filt, aes(x = .data[[x]], y = .data[[y]], color = sig_class)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = pal, drop = FALSE) +

    # Quadrant labels
    annotate("text", x = max_abs_x * 0.7, y = max_abs_y * 0.85,
             label = sprintf("Up-Up\n(n=%d)", q_uu), size = 3, fontface = "bold", color = "grey30") +
    annotate("text", x = max_abs_x * 0.7, y = -max_abs_y * 0.85,
             label = sprintf("Up-Down\n(n=%d)", q_ud), size = 3, fontface = "bold", color = "grey30") +
    annotate("text", x = -max_abs_x * 0.7, y = max_abs_y * 0.85,
             label = sprintf("Down-Up\n(n=%d)", q_du), size = 3, fontface = "bold", color = "grey30") +
    annotate("text", x = -max_abs_x * 0.7, y = -max_abs_y * 0.85,
             label = sprintf("Down-Down\n(n=%d)", q_dd), size = 3, fontface = "bold", color = "grey30") +

    labs(
      title = title,
      subtitle = sprintf("Concordant: %d (%.1f%%) | Fisher p = %.2e",
                         concordant, 100 * concordant / total, fisher_p),
      x = xlab, y = ylab, color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(color = "grey40", size = 9),
      legend.position = "bottom"
    ) +
    coord_cartesian(xlim = c(-max_abs_x, max_abs_x), ylim = c(-max_abs_y, max_abs_y))

  p
}

fig2_quadrant <- function(df, figdir) {
  msg("Figure 2: Quadrant concordance plots")

  # Filter to genes with at least one significant result
  df_sig <- df %>% filter(any_sig)

  df_rna_prom <- df_sig %>% filter(!is.na(rna_log2FC), !is.na(prom_best_log2FC))
  df_rna_enh <- df_sig %>% filter(!is.na(rna_log2FC), !is.na(enh_best_log2FC))

  p1 <- make_quadrant_panel(df_rna_prom, "rna_log2FC", "prom_best_log2FC",
                            "RNA log2FC", "Promoter log2FC",
                            "RNA vs Promoter (sig genes)")
  p2 <- make_quadrant_panel(df_rna_enh, "rna_log2FC", "enh_best_log2FC",
                            "RNA log2FC", "Enhancer log2FC",
                            "RNA vs Enhancer (sig genes)")

  p_combined <- (p1 | p2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  p_combined <- p_combined +
    plot_annotation(title = "Quadrant Analysis: Direction Concordance",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))

  ggsave(file.path(figdir, "Fig2_quadrant_concordance.png"), p_combined,
         width = 12, height = 6, dpi = 300)

  list(p1 = p1, p2 = p2)
}

# ============================================================================
# Figure 3: UpSet + Venn
# ============================================================================

fig3_overlap <- function(df, figdir) {
  msg("Figure 3: UpSet plot and Venn diagram")

  # Prepare sets
  rna_genes <- df$gene_id_clean[df$rna_sig]
  prom_genes <- df$gene_id_clean[df$prom_sig]
  enh_genes <- df$gene_id_clean[df$enh_sig]

  gene_list <- list(
    RNA = rna_genes,
    Promoter = prom_genes,
    Enhancer = enh_genes
  )

  # UpSet plot
  png(file.path(figdir, "Fig3a_upset.png"), width = 8, height = 6, units = "in", res = 300)
  print(upset(
    fromList(gene_list),
    order.by = "freq",
    sets.bar.color = c("#CC79A7", "#D55E00", "#0072B2"),
    main.bar.color = "grey30",
    text.scale = 1.3,
    point.size = 3
  ))
  dev.off()
  msg("  Saved Fig3a_upset.png")

  # Venn diagram
  p_venn <- ggVennDiagram(
    gene_list,
    label = "count",
    label_alpha = 0,
    set_color = c("#CC79A7", "#D55E00", "#0072B2")
  ) +
    scale_fill_gradient(low = "white", high = "#0072B2", name = "Count") +
    labs(title = "Overlap of Significant Genes") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(file.path(figdir, "Fig3b_venn.png"), p_venn, width = 7, height = 6, dpi = 300)
  msg("  Saved Fig3b_venn.png")

  invisible(gene_list)
}

# ============================================================================
# Figure 4: Focused Scatter (sig genes only, better labeling)
# ============================================================================

fig4_focused_scatter <- function(df, figdir, lfc_threshold = 0.5, max_labels = 35) {
  msg("Figure 4: Focused scatter (sig genes only)")

  select_focus_labels <- function(panel_df, padj_a, padj_b, sig_a, sig_b, lfc_a, lfc_b,
                                  lfc_threshold = 0.5, max_labels = 35) {
    ranked <- panel_df %>%
      mutate(
        sig_a_val = dplyr::coalesce(as.logical(.data[[sig_a]]), FALSE),
        sig_b_val = dplyr::coalesce(as.logical(.data[[sig_b]]), FALSE),
        axis_sig = sig_a_val | sig_b_val,
        both_sig = sig_a_val & sig_b_val,
        best_padj = pmin(
          ifelse(is.na(.data[[padj_a]]), 1, .data[[padj_a]]),
          ifelse(is.na(.data[[padj_b]]), 1, .data[[padj_b]])
        ),
        max_abs_lfc = pmax(abs(.data[[lfc_a]]), abs(.data[[lfc_b]]), na.rm = TRUE),
        max_abs_lfc = ifelse(is.finite(max_abs_lfc), max_abs_lfc, 0)
      ) %>%
      filter(!is.na(gene_label), nzchar(gene_label), axis_sig)

    candidates <- ranked %>% filter(max_abs_lfc >= lfc_threshold)
    if (nrow(candidates) < min(10, max_labels)) {
      candidates <- ranked
    }

    candidates %>%
      arrange(desc(both_sig), best_padj, desc(max_abs_lfc)) %>%
      distinct(gene_label, .keep_all = TRUE) %>%
      slice_head(n = max_labels)
  }

  df_sig <- df %>%
    filter(any_sig) %>%
    mutate(
      # Concordance category for coloring
      concordance = case_when(
        !rna_sig | (!prom_sig & !enh_sig) ~ "One modality",
        (rna_log2FC > 0 & (prom_best_log2FC > 0 | enh_best_log2FC > 0)) |
          (rna_log2FC < 0 & (prom_best_log2FC < 0 | enh_best_log2FC < 0)) ~ "Concordant",
        TRUE ~ "Discordant"
      )
    )

  pal <- c("Concordant" = "#0072B2", "Discordant" = "#D55E00", "One modality" = "grey60")

  # RNA vs Promoter
  df_rna_prom <- df_sig %>% filter(!is.na(rna_log2FC), !is.na(prom_best_log2FC))
  lab_rna_prom <- select_focus_labels(
    df_rna_prom,
    padj_a = "rna_padj", padj_b = "prom_best_padj",
    sig_a = "rna_sig", sig_b = "prom_sig",
    lfc_a = "rna_log2FC", lfc_b = "prom_best_log2FC",
    lfc_threshold = lfc_threshold, max_labels = max_labels
  )

  p1 <- ggplot(df_rna_prom, aes(x = rna_log2FC, y = prom_best_log2FC, color = concordance)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_point(alpha = 0.45, size = 1.5) +
    geom_point(
      data = lab_rna_prom,
      aes(fill = concordance),
      shape = 21, size = 2.4, stroke = 0.35, color = "black",
      alpha = 0.95, show.legend = FALSE
    ) +
    geom_text_repel(
      data = lab_rna_prom,
      aes(label = gene_label),
      size = 3.0, fontface = "italic", max.overlaps = Inf,
      box.padding = 0.3, point.padding = 0.2, segment.color = "grey50"
    ) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      title = "RNA vs Promoter (significant genes)",
      x = "RNA log2FC", y = "Promoter log2FC", color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

  # RNA vs Enhancer
  df_rna_enh <- df_sig %>% filter(!is.na(rna_log2FC), !is.na(enh_best_log2FC))
  lab_rna_enh <- select_focus_labels(
    df_rna_enh,
    padj_a = "rna_padj", padj_b = "enh_best_padj",
    sig_a = "rna_sig", sig_b = "enh_sig",
    lfc_a = "rna_log2FC", lfc_b = "enh_best_log2FC",
    lfc_threshold = lfc_threshold, max_labels = max_labels
  )

  p2 <- ggplot(df_rna_enh, aes(x = rna_log2FC, y = enh_best_log2FC, color = concordance)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_point(alpha = 0.45, size = 1.5) +
    geom_point(
      data = lab_rna_enh,
      aes(fill = concordance),
      shape = 21, size = 2.4, stroke = 0.35, color = "black",
      alpha = 0.95, show.legend = FALSE
    ) +
    geom_text_repel(
      data = lab_rna_enh,
      aes(label = gene_label),
      size = 3.0, fontface = "italic", max.overlaps = Inf,
      box.padding = 0.3, point.padding = 0.2, segment.color = "grey50"
    ) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      title = "RNA vs Enhancer (significant genes)",
      x = "RNA log2FC", y = "Enhancer log2FC", color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

  p_combined <- (p1 | p2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(file.path(figdir, "Fig4_focused_scatter.png"), p_combined,
         width = 12, height = 6, dpi = 300)

  list(p1 = p1, p2 = p2)
}

# ============================================================================
# Figure 5: Stratified Boxplots
# ============================================================================

fig5_stratified_boxplots <- function(df, figdir) {
  msg("Figure 5: Stratified boxplots")

  # Use RNA direction as stratification
  df_strat <- df %>%
    filter(!is.na(rna_dir)) %>%
    mutate(rna_dir = factor(rna_dir, levels = c("Down", "NS", "Up")))

  # Promoter boxplot
  df_prom <- df_strat %>% filter(!is.na(prom_best_log2FC))

  # Kruskal-Wallis test
  kw_prom <- kruskal.test(prom_best_log2FC ~ rna_dir, data = df_prom)

  p1 <- ggplot(df_prom, aes(x = rna_dir, y = prom_best_log2FC, fill = rna_dir)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Down" = "#3182bd", "NS" = "grey70", "Up" = "#e6550d")) +
    labs(
      title = "Promoter ATAC by RNA Direction",
      subtitle = sprintf("Kruskal-Wallis p = %.2e", kw_prom$p.value),
      x = "RNA Direction (KO vs WT)", y = "Promoter log2FC"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")

  # Enhancer boxplot
  df_enh <- df_strat %>% filter(!is.na(enh_best_log2FC))

  kw_enh <- kruskal.test(enh_best_log2FC ~ rna_dir, data = df_enh)

  p2 <- ggplot(df_enh, aes(x = rna_dir, y = enh_best_log2FC, fill = rna_dir)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Down" = "#3182bd", "NS" = "grey70", "Up" = "#e6550d")) +
    labs(
      title = "Enhancer ATAC by RNA Direction",
      subtitle = sprintf("Kruskal-Wallis p = %.2e", kw_enh$p.value),
      x = "RNA Direction (KO vs WT)", y = "Enhancer log2FC"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")

  p_combined <- (p1 | p2) +
    plot_annotation(title = "ATAC Changes Stratified by RNA Direction",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))

  ggsave(file.path(figdir, "Fig5_stratified_boxplots.png"), p_combined,
         width = 10, height = 5, dpi = 300)

  list(p1 = p1, p2 = p2, kw_prom = kw_prom, kw_enh = kw_enh)
}

# ============================================================================
# Figure 6: Top DEG Heatmap
# ============================================================================

fig6_heatmap <- function(df, figdir, topN = 50) {
  msg("Figure 6: Top DEG heatmap")

  # Top N by abs(RNA log2FC)
  top_genes <- df %>%
    filter(!is.na(rna_log2FC)) %>%
    arrange(desc(abs(rna_log2FC))) %>%
    slice_head(n = topN) %>%
    filter(!is.na(prom_best_log2FC) | !is.na(enh_best_log2FC))

  if (nrow(top_genes) < 5) {
    msg("  Not enough genes with both RNA and ATAC data for heatmap")
    return(NULL)
  }

  # Build matrix
  mat <- top_genes %>%
    select(gene_label, rna_log2FC, prom_best_log2FC, enh_best_log2FC) %>%
    column_to_rownames("gene_label") %>%
    as.matrix()

  colnames(mat) <- c("RNA", "Promoter", "Enhancer")

  # Annotation for significance
  annot_row <- top_genes %>%
    transmute(
      gene_label,
      RNA_sig = ifelse(rna_sig, "sig", "ns"),
      Prom_sig = ifelse(prom_sig, "sig", "ns"),
      Enh_sig = ifelse(enh_sig, "sig", "ns")
    ) %>%
    column_to_rownames("gene_label")

  annot_colors <- list(
    RNA_sig = c("sig" = "#0072B2", "ns" = "grey90"),
    Prom_sig = c("sig" = "#D55E00", "ns" = "grey90"),
    Enh_sig = c("sig" = "#009E73", "ns" = "grey90")
  )

  # Color scale centered at 0
  max_abs <- max(abs(mat), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("#3182bd", "white", "#e6550d"))(100)

  png(file.path(figdir, "Fig6_top_DEG_heatmap.png"), width = 8, height = 12, units = "in", res = 300)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,  # Keep ordered by RNA effect size
    cluster_cols = FALSE,
    breaks = breaks,
    color = colors,
    na_col = "grey95",
    annotation_row = annot_row,
    annotation_colors = annot_colors,
    main = sprintf("Top %d DEGs: RNA vs ATAC log2FC", nrow(mat)),
    fontsize = 9,
    fontsize_row = 7,
    cellwidth = 30,
    cellheight = 10
  )
  dev.off()

  msg("  Saved Fig6_top_DEG_heatmap.png")

  invisible(mat)
}

# ============================================================================
# Figure 7: Promoter vs Enhancer Comparison (Barbell)
# ============================================================================

fig7_barbell <- function(df, figdir, topN = 40) {
  msg("Figure 7: Promoter vs Enhancer barbell")

  # Filter genes with both promoter and enhancer data, and RNA sig
  df_both <- df %>%
    filter(rna_sig, !is.na(prom_best_log2FC), !is.na(enh_best_log2FC)) %>%
    arrange(desc(abs(rna_log2FC))) %>%
    slice_head(n = topN) %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_label)),
      rna_direction = ifelse(rna_log2FC > 0, "RNA Up", "RNA Down")
    )

  if (nrow(df_both) < 5) {
    msg("  Not enough genes for barbell plot")
    return(NULL)
  }

  # Prepare for geom_segment
  df_long <- df_both %>%
    select(gene_label, rna_direction, prom_best_log2FC, enh_best_log2FC) %>%
    pivot_longer(cols = c(prom_best_log2FC, enh_best_log2FC),
                 names_to = "region", values_to = "log2FC") %>%
    mutate(region = ifelse(region == "prom_best_log2FC", "Promoter", "Enhancer"))

  p <- ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    # Connecting segments
    geom_segment(
      data = df_both,
      aes(x = prom_best_log2FC, xend = enh_best_log2FC, y = gene_label, yend = gene_label),
      color = "grey60", linewidth = 0.6
    ) +
    # Points for promoter and enhancer
    geom_point(
      data = df_long,
      aes(x = log2FC, y = gene_label, color = region),
      size = 2.5
    ) +
    scale_color_manual(values = c("Promoter" = "#D55E00", "Enhancer" = "#009E73")) +
    facet_wrap(~rna_direction, scales = "free_y", ncol = 2) +
    labs(
      title = "Promoter vs Enhancer ATAC Changes for Top RNA DEGs",
      subtitle = "Line connects promoter and enhancer log2FC for each gene",
      x = "ATAC log2FC", y = NULL, color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )

  ggsave(file.path(figdir, "Fig7_prom_vs_enh_barbell.png"), p,
         width = 12, height = 10, dpi = 300)

  p
}

# ============================================================================
# Main
# ============================================================================

option_list <- list(
  make_option("--rna_deg_file", type = "character",
              help = "RNA DEGs CSV (must include gene_id + log2FoldChange + padj/FDR)"),
  make_option("--datadir", type = "character",
              help = "Directory containing promoter/enhancer outputs from prep script"),
  make_option("--figdir", type = "character", default = NULL,
              help = "Directory to save figures (default: <datadir>/figures_exploratory)"),
  make_option("--alpha", type = "double", default = 0.05,
              help = "Significance threshold [default: 0.05]"),
  make_option("--topN", type = "integer", default = 50,
              help = "Number of top DEGs for heatmap/barbell [default: 50]")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$rna_deg_file), file.exists(opt$rna_deg_file))
stopifnot(!is.null(opt$datadir), dir.exists(opt$datadir))

figdir <- if (is.null(opt$figdir) || opt$figdir == "") {
  file.path(opt$datadir, "figures_exploratory")
} else {
  opt$figdir
}
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

msg("=== RNA vs ATAC Exploratory Figures ===")
msg("RNA DEG file: %s", opt$rna_deg_file)
msg("Data directory: %s", opt$datadir)
msg("Figure directory: %s", figdir)
msg("Alpha: %.3f", opt$alpha)

# Load and merge data
df <- load_and_merge_data(opt)
msg("Loaded %d genes", nrow(df))
msg("  RNA sig: %d", sum(df$rna_sig, na.rm = TRUE))
msg("  Promoter sig: %d", sum(df$prom_sig, na.rm = TRUE))
msg("  Enhancer sig: %d", sum(df$enh_sig, na.rm = TRUE))

# Generate figures
fig1_hexbin(df, figdir)
fig2_quadrant(df, figdir)
fig3_overlap(df, figdir)
fig4_focused_scatter(df, figdir)
fig5_stratified_boxplots(df, figdir)
fig6_heatmap(df, figdir, topN = opt$topN)
fig7_barbell(df, figdir, topN = opt$topN)

msg("=== DONE ===")
msg("All figures saved to: %s", figdir)
