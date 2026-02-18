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
  library(tibble)
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

stouffer_combine <- function(pvals, lfcs, weights = NULL) {
  p <- as.numeric(pvals)
  fc <- as.numeric(lfcs)
  keep <- is.finite(p) & p > 0 & p < 1 & is.finite(fc)
  p <- p[keep]; fc <- fc[keep]
  if (length(p) == 0) return(list(pvalue = NA_real_, direction = NA_character_))
  z <- sign(fc) * qnorm(1 - p / 2)
  if (!is.null(weights)) {
    w <- as.numeric(weights)[keep]
    z_combined <- sum(w * z) / sqrt(sum(w^2))
  } else {
    z_combined <- sum(z) / sqrt(length(z))
  }
  p_combined <- 2 * pnorm(-abs(z_combined))
  dir <- if (z_combined > 0) "up" else if (z_combined < 0) "down" else "ambiguous"
  list(pvalue = p_combined, direction = dir)
}

collapse_peak_to_gene <- function(da_df, map_df, alpha = 0.05) {
  # Validate raw pvalue column
  if (!"pvalue" %in% colnames(da_df)) {
    pval_col <- intersect(colnames(da_df), c("pvalue", "pval", "p_value"))
    if (length(pval_col) == 0) {
      stop("DA results missing raw 'pvalue' column. Re-run deseq_atac.R (updated version) to generate it.",
           call. = FALSE)
    }
    da_df$pvalue <- da_df[[pval_col[1]]]
  }

  da2 <- da_df %>%
    mutate(
      peak_id = as.character(peak_id),
      padj = as.numeric(padj),
      pvalue = as.numeric(pvalue),
      log2FoldChange = as.numeric(log2FoldChange)
    )

  # Detect 7-column (weighted) vs 6-column (legacy) map files
  has_weight <- ncol(map_df) >= 7
  if (has_weight) {
    map2 <- map_df %>%
      transmute(
        peak_id = as.character(.[[1]]),
        gene_id_clean = clean_ens(.[[2]]),
        gene_name = as.character(.[[3]]),
        weight = as.numeric(.[[7]])
      )
  } else {
    map2 <- map_df %>%
      transmute(
        peak_id = as.character(.[[1]]),
        gene_id_clean = clean_ens(.[[2]]),
        gene_name = as.character(.[[3]])
      )
  }

  joined <- inner_join(map2, da2, by = "peak_id")
  if (nrow(joined) == 0) {
    return(tibble(
      gene_id_clean = character(), n_peaks = integer(), n_sig = integer(),
      best_peak_id = character(), best_log2FC = double(), best_padj = double(),
      combined_pvalue = double(), combined_padj = double(),
      combined_direction = character(), weighted_log2FC = double(),
      mean_log2FC = double(), median_log2FC = double()
    ))
  }

  # Stouffer combination per gene
  gene_stats <- joined %>%
    group_by(gene_id_clean) %>%
    summarise(
      n_peaks = n(),
      n_sig = sum(!is.na(padj) & padj < alpha),
      stouffer = {
        w <- if (has_weight && "weight" %in% names(cur_data())) weight else NULL
        list(stouffer_combine(pvalue, log2FoldChange, w))
      },
      weighted_log2FC = if (has_weight && "weight" %in% names(cur_data())) {
        wt <- weight / sum(weight)
        sum(wt * log2FoldChange, na.rm = TRUE)
      } else {
        mean(log2FoldChange, na.rm = TRUE)
      },
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      median_log2FC = median(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      combined_pvalue = sapply(stouffer, function(x) x$pvalue),
      combined_direction = sapply(stouffer, function(x) x$direction)
    ) %>%
    select(-stouffer)

  # BH correction AFTER combination
  gene_stats$combined_padj <- p.adjust(gene_stats$combined_pvalue, method = "BH")

  # Best peak per gene (legacy compat)
  best_rows <- joined %>%
    group_by(gene_id_clean) %>%
    arrange(is.na(padj), padj, desc(abs(log2FoldChange))) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene_id_clean, best_peak_id = peak_id, best_log2FC = log2FoldChange, best_padj = padj)

  left_join(gene_stats, best_rows, by = "gene_id_clean")
}

load_and_merge_data <- function(opt) {
  # Load mapping files (auto-detect 7-col weighted vs 6-col legacy)
  prom_map_raw <- read_tsv(
    file.path(opt$datadir, "promoter_closest.tsv"),
    col_names = FALSE, show_col_types = FALSE
  )
  n_prom_cols <- ncol(prom_map_raw)
  if (n_prom_cols >= 7) {
    colnames(prom_map_raw) <- c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist", "weight")
  } else {
    colnames(prom_map_raw) <- c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist")[seq_len(n_prom_cols)]
  }
  prom_map <- prom_map_raw

  enh_map_raw <- read_tsv(
    file.path(opt$datadir, "enhancer_closest.tsv"),
    col_names = FALSE, show_col_types = FALSE
  )
  n_enh_cols <- ncol(enh_map_raw)
  if (n_enh_cols >= 7) {
    colnames(enh_map_raw) <- c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist", "weight")
  } else {
    colnames(enh_map_raw) <- c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist")[seq_len(n_enh_cols)]
  }
  enh_map <- enh_map_raw

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

  pval_raw_col <- intersect(names(rna), c("pvalue", "pval", "p_value"))

  rna2 <- rna %>%
    transmute(
      gene_id_clean = clean_ens(.data[[gene_id_col]]),
      gene_label = if (!is.null(gene_name_col)) as.character(.data[[gene_name_col]]) else clean_ens(.data[[gene_id_col]]),
      rna_log2FC = as.numeric(.data[[lfc_col]]),
      rna_padj = as.numeric(.data[[padj_col]]),
      rna_pvalue = if (length(pval_raw_col) > 0) {
        as.numeric(.data[[pval_raw_col[1]]])
      } else {
        warning("RNA DEGs missing raw pvalue column; using padj as fallback")
        as.numeric(.data[[padj_col]])
      }
    ) %>%
    distinct(gene_id_clean, .keep_all = TRUE)

  # Merge
  gene_merged <- rna2 %>%
    left_join(prom_gene, by = "gene_id_clean") %>%
    left_join(enh_gene, by = "gene_id_clean") %>%
    mutate(
      rna_sig = !is.na(rna_padj) & rna_padj < opt$alpha,
      prom_sig = !is.na(prom_combined_padj) & prom_combined_padj < opt$alpha,
      enh_sig = !is.na(enh_combined_padj) & enh_combined_padj < opt$alpha,
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
    padj_a = "rna_padj", padj_b = "prom_combined_padj",
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
      box.padding = 0.3, point.padding = 0.2, segment.color = "grey50",
      show.legend = FALSE
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
    padj_a = "rna_padj", padj_b = "enh_combined_padj",
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
      box.padding = 0.3, point.padding = 0.2, segment.color = "grey50",
      show.legend = FALSE
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

fig6_heatmap <- function(df, figdir, counts_file = NULL, metadata_file = NULL,
                         datadir = NULL, topN = 50) {
  msg("Figure 6: Sample-level ATAC accessibility heatmap")

  if (is.null(counts_file) || !file.exists(counts_file)) {
    msg("  No counts file provided or not found; skipping Fig6")
    return(NULL)
  }
  if (is.null(metadata_file) || !file.exists(metadata_file)) {
    msg("  No metadata file provided or not found; skipping Fig6")
    return(NULL)
  }

  # Load sample metadata
  meta <- readr::read_tsv(metadata_file, show_col_types = FALSE)
  meta <- meta %>%
    mutate(sample_clean = gsub("[^A-Za-z0-9]", "", sample))

  # Load count matrix (featureCounts format: 6 annotation cols + sample cols)
  counts_raw <- readr::read_tsv(counts_file, comment = "#", show_col_types = FALSE)
  peak_ids <- counts_raw[[1]]
  count_cols <- counts_raw[, -(1:6)]
  colnames(count_cols) <- gsub(".*/", "", colnames(count_cols))
  colnames(count_cols) <- gsub("\\.final\\.bam$", "", colnames(count_cols))

  # Match samples
  common <- intersect(meta$sample_clean, colnames(count_cols))
  if (length(common) < 2) {
    msg("  Fewer than 2 matching samples; skipping Fig6")
    return(NULL)
  }
  meta <- meta %>% filter(sample_clean %in% common) %>% arrange(condition, sample_clean)
  count_mat <- as.matrix(count_cols[, meta$sample_clean])
  rownames(count_mat) <- peak_ids

  # Select top DA peaks with gene labels
  top_genes <- df %>%
    filter(any_sig, !is.na(gene_label), nzchar(gene_label)) %>%
    mutate(
      n_sig = as.integer(rna_sig) + as.integer(prom_sig) + as.integer(enh_sig),
      best_padj = pmin(
        ifelse(is.na(rna_padj), 1, rna_padj),
        ifelse(is.na(prom_combined_padj), 1, prom_combined_padj),
        ifelse(is.na(enh_combined_padj), 1, enh_combined_padj)
      )
    ) %>%
    arrange(desc(n_sig), best_padj) %>%
    distinct(gene_label, .keep_all = TRUE) %>%
    slice_head(n = topN)

  if (nrow(top_genes) < 5) {
    msg("  Not enough significant genes for heatmap")
    return(NULL)
  }

  # Map genes to their best promoter peak (closest, most significant)
  prom_map_file <- NULL
  candidates <- c(
    if (!is.null(datadir)) file.path(datadir, "promoter_multigene.tsv"),
    file.path(dirname(counts_file), "..", "gene_integration", "promoter_multigene.tsv")
  )
  for (cand in candidates) {
    if (file.exists(cand)) { prom_map_file <- cand; break }
  }

  if (!is.null(prom_map_file) && file.exists(prom_map_file)) {
    prom_map <- readr::read_tsv(
      prom_map_file,
      col_names = c("peak_id", "gene_id", "gene_name", "strand", "dist_bp", "abs_dist", "weight"),
      show_col_types = FALSE
    )
    # Best peak per gene (closest to TSS)
    peak_gene <- prom_map %>%
      filter(gene_name %in% top_genes$gene_label) %>%
      arrange(abs_dist) %>%
      distinct(gene_name, .keep_all = TRUE)
  } else {
    # Fall back: try to match peak IDs from df if available
    msg("  Promoter map not found; attempting peak_id match from DA results")
    peak_gene <- top_genes %>%
      filter(!is.na(prom_best_peak_id)) %>%
      transmute(peak_id = prom_best_peak_id, gene_name = gene_label)
  }

  # Subset count matrix to selected peaks
  valid_peaks <- intersect(peak_gene$peak_id, rownames(count_mat))
  if (length(valid_peaks) < 5) {
    msg("  Too few peaks matched (%d); skipping Fig6", length(valid_peaks))
    return(NULL)
  }

  peak_gene_sub <- peak_gene %>% filter(peak_id %in% valid_peaks)
  mat_sub <- count_mat[peak_gene_sub$peak_id, , drop = FALSE]

  # Normalize: log2(CPM + 1)
  lib_sizes <- colSums(count_mat)
  mat_cpm <- t(t(mat_sub) / lib_sizes * 1e6)
  mat_log <- log2(mat_cpm + 1)

  # Z-score per peak (row)
  mat_z <- t(scale(t(mat_log)))
  rownames(mat_z) <- peak_gene_sub$gene_name

  # Remove rows with zero variance (all-NA z-scores)
  keep_rows <- apply(mat_z, 1, function(x) !all(is.na(x)))
  mat_z <- mat_z[keep_rows, , drop = FALSE]

  if (nrow(mat_z) < 5) {
    msg("  Too few valid rows after z-scoring; skipping Fig6")
    return(NULL)
  }

  # Column annotation: condition
  col_annot <- data.frame(
    Condition = meta$condition,
    row.names = meta$sample_clean
  )
  annot_colors <- list(Condition = c("WT" = "#0072B2", "KO" = "#D55E00"))

  # Color scale: blue-white-red centered at 0
  max_z <- min(max(abs(mat_z), na.rm = TRUE), 3)
  breaks <- seq(-max_z, max_z, length.out = 101)
  colors <- colorRampPalette(c("#3182bd", "white", "#e6550d"))(100)

  png(file.path(figdir, "Fig6_sample_heatmap.png"),
      width = 8, height = max(8, nrow(mat_z) * 0.2 + 2),
      units = "in", res = 300)
  pheatmap::pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    breaks = breaks,
    color = colors,
    na_col = "grey95",
    annotation_col = col_annot,
    annotation_colors = annot_colors,
    main = sprintf("ATAC Accessibility: Top %d DA Genes (z-scored log2CPM)", nrow(mat_z)),
    fontsize = 9,
    fontsize_row = 7,
    show_colnames = TRUE
  )
  dev.off()

  msg("  Saved Fig6_sample_heatmap.png (%d genes x %d samples)", nrow(mat_z), ncol(mat_z))
  invisible(mat_z)
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
# Figure 8: Sensitivity Analysis Comparison
# ============================================================================

fig8_sensitivity <- function(datadir, figdir) {
  msg("Figure 8: Sensitivity analysis comparison")

  summary_file <- file.path(datadir, "sensitivity", "summary.tsv")
  if (!file.exists(summary_file)) {
    msg("  No sensitivity summary found; skipping Fig8")
    return(NULL)
  }

  sens <- read_tsv(summary_file, show_col_types = FALSE)
  if (nrow(sens) == 0) {
    msg("  Empty sensitivity summary; skipping Fig8")
    return(NULL)
  }

  sens_long <- sens %>%
    pivot_longer(cols = starts_with("n_"), names_to = "metric", values_to = "count") %>%
    mutate(
      metric = gsub("n_", "", metric),
      metric = gsub("_", " ", metric),
      window = factor(window, levels = sort(unique(window)))
    )

  p1 <- ggplot(sens_long, aes(x = window, y = count, fill = metric)) +
    geom_col(position = "dodge") +
    labs(
      title = "Gene/Peak Counts by Enhancer Window Size",
      x = "Enhancer window (bp)", y = "Count", fill = "Metric"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

  # Jaccard overlap heatmap across window sizes
  sens_dirs <- list.dirs(file.path(datadir, "sensitivity"), recursive = FALSE)
  if (length(sens_dirs) >= 2) {
    window_genes <- list()
    for (d in sens_dirs) {
      w <- basename(d)
      enh_file <- file.path(d, "enhancer_multigene.tsv")
      if (file.exists(enh_file)) {
        genes <- unique(read_tsv(enh_file, col_names = FALSE, show_col_types = FALSE)[[3]])
        window_genes[[w]] <- genes
      }
    }

    if (length(window_genes) >= 2) {
      windows <- names(window_genes)
      jaccard_mat <- matrix(0, nrow = length(windows), ncol = length(windows),
                            dimnames = list(windows, windows))
      for (i in seq_along(windows)) {
        for (j in seq_along(windows)) {
          a <- window_genes[[windows[i]]]
          b <- window_genes[[windows[j]]]
          jaccard_mat[i, j] <- length(intersect(a, b)) / length(union(a, b))
        }
      }

      png(file.path(figdir, "Fig8b_sensitivity_jaccard.png"),
          width = 7, height = 6, units = "in", res = 300)
      pheatmap::pheatmap(
        jaccard_mat,
        main = "Jaccard Overlap: Enhancer Genes by Window Size",
        display_numbers = TRUE,
        number_format = "%.2f",
        color = colorRampPalette(c("white", "#0072B2"))(50)
      )
      dev.off()
      msg("  Saved Fig8b_sensitivity_jaccard.png")
    }
  }

  ggsave(file.path(figdir, "Fig8_sensitivity_comparison.png"), p1,
         width = 8, height = 5, dpi = 300)
  msg("  Saved Fig8_sensitivity_comparison.png")

  p1
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
  make_option("--counts", type = "character", default = NULL,
              help = "featureCounts output file for sample-level heatmap (optional)"),
  make_option("--metadata", type = "character", default = NULL,
              help = "Sample metadata TSV (sample, condition, replicate) for heatmap (optional)"),
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
fig6_heatmap(df, figdir, counts_file = opt$counts, metadata_file = opt$metadata,
             datadir = opt$datadir, topN = opt$topN)
fig7_barbell(df, figdir, topN = opt$topN)
fig8_sensitivity(opt$datadir, figdir)

msg("=== DONE ===")
msg("All figures saved to: %s", figdir)
