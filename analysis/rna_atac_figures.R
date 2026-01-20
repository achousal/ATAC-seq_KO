#!/usr/bin/env Rscript

# example run:
# Rscript rna_atac_figures.R \
#   --rna_deg_file "/sc/arion/projects/Chipuk_Laboratory/henaoj02/20251016_Bulk_RNAseq_Fibroblast_ATF5_WT_vs_KO/BiNGS-Bulk-RNA-seq-pipeline/data_rna/processed/analysis/differential_expression/condition_ATF5KO_vs_WT_degs.csv" \
#   --datadir /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq/ATF5_KO/analysis/gene_integration \
#   --alpha 0.05 \
#   --topN 30

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
})

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n", sep = "")

clean_ens <- function(x) {
  x <- as.character(x)
  gsub("\\.\\d+$", "", x)
}

# Collapse peak-level DA to gene-level using the "best peak" rule:
# best = min padj (NA last); tie-break = max abs(log2FC)
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

qrange <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 10) return(NULL)
  as.numeric(quantile(v, probs = c(0.01, 0.99), na.rm = TRUE))
}

make_panel <- function(df, x, y, xlab, ylab, title,
                       label_genes = NULL, xlim = NULL, ylim = NULL) {

  stopifnot(all(c(x, y, "sig_class") %in% colnames(df)))

  # enforce consistent legend order (adjust if your strings differ)
  levs <- c("not sig", "ATAC only", "RNA only", "RNA + ATAC")
  df <- df %>%
    dplyr::mutate(sig_class = factor(sig_class, levels = levs))

  # Colorblind-friendly-ish palette (Okabe-Ito variants)
  pal <- c(
    "not sig"    = "grey80",
    "ATAC only"  = "#D55E00",
    "RNA only"   = "#CC79A7",
    "RNA + ATAC" = "#0072B2"
  )

  # Make non-sig quiet without losing them
  alpha_map <- c(
    "not sig"    = 0.12,
    "ATAC only"  = 0.55,
    "RNA only"   = 0.60,
    "RNA + ATAC" = 0.75
  )

  # Spearman correlation (robust to outliers)
  xvals <- df[[x]]; yvals <- df[[y]]
  ok <- is.finite(xvals) & is.finite(yvals)
  rho <- suppressWarnings(stats::cor(xvals[ok], yvals[ok], method = "spearman"))
  n_ok <- sum(ok)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +

    # points (alpha varies by class; legend only for color)
    ggplot2::geom_point(
      ggplot2::aes(color = sig_class, alpha = sig_class),
      size = 1.35
    ) +

    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = alpha_map, guide = "none", drop = FALSE) +

    ggplot2::labs(
      title = title,
      subtitle = sprintf("Spearman ρ = %.2f; n = %d", rho, n_ok),
      x = xlab, y = ylab, color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      legend.position = "bottom"
    )

  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE)
  }

  # Optional: add y=x dashed line for promoter vs enhancer panels (helps interpretation)
  if (grepl("Promoter", xlab) && grepl("Enhancer", ylab)) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                                 linewidth = 0.35, color = "grey75")
  }

  if (!is.null(label_genes) && length(label_genes) > 0 && "gene_label" %in% names(df)) {
    lab_df <- df %>% dplyr::filter(!is.na(gene_label), gene_label %in% label_genes)

    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = gene_label),
      size = 3.1,
      max.overlaps = 30,          # <- key: don’t let labels swamp the panel
      box.padding = 0.30,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.color = "grey45",
      segment.alpha = 0.6,
      seed = 1
    )
  }

  p
}

# ---------------- args ----------------
option_list <- list(
  make_option("--rna_deg_file", type="character",
              help="RNA DEGs CSV (must include gene_id + log2FoldChange + padj/FDR)"),
  make_option("--datadir", type="character",
              help="Directory containing promoter/enhancer outputs from prep script; also where CSV tables will be written"),
  make_option("--figdir", type="character", default=NULL,
              help="Directory to save figures (default: <datadir>/figures)"),
  make_option("--alpha", type="double", default=0.05),
  make_option("--topN", type="integer", default=30),

  # produced by prep .sh (defaults assume datadir)
  make_option("--promoter_closest_tsv", type="character", default=NULL),
  make_option("--enhancer_closest_tsv", type="character", default=NULL),
  make_option("--promoter_DA_csv", type="character", default=NULL),
  make_option("--enhancer_DA_csv", type="character", default=NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$rna_deg_file), file.exists(opt$rna_deg_file))
stopifnot(!is.null(opt$datadir), dir.exists(opt$datadir))

prom_map_file <- ifelse(is.null(opt$promoter_closest_tsv),
                        file.path(opt$datadir, "promoter_closest.tsv"),
                        opt$promoter_closest_tsv)
enh_map_file  <- ifelse(is.null(opt$enhancer_closest_tsv),
                        file.path(opt$datadir, "enhancer_closest.tsv"),
                        opt$enhancer_closest_tsv)
prom_da_file  <- ifelse(is.null(opt$promoter_DA_csv),
                        file.path(opt$datadir, "promoter_DA.csv"),
                        opt$promoter_DA_csv)
enh_da_file   <- ifelse(is.null(opt$enhancer_DA_csv),
                        file.path(opt$datadir, "enhancer_DA.csv"),
                        opt$enhancer_DA_csv)

stopifnot(file.exists(prom_map_file), file.exists(enh_map_file))
stopifnot(file.exists(prom_da_file), file.exists(enh_da_file))

msg("datadir: %s", opt$datadir)

figdir <- if (is.null(opt$figdir) || opt$figdir == "") {
  file.path(opt$datadir, "figures")
} else {
  opt$figdir
}
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
msg("Figdir: %s", figdir)

# ---------------- load prep outputs ----------------
prom_map <- read_tsv(
  prom_map_file,
  col_names = c("peak_id","gene_id","gene_name","strand","dist_bp","abs_dist"),
  show_col_types = FALSE
)
enh_map <- read_tsv(
  enh_map_file,
  col_names = c("peak_id","gene_id","gene_name","strand","dist_bp","abs_dist"),
  show_col_types = FALSE
)

prom_da <- read_csv(prom_da_file, show_col_types = FALSE)
enh_da  <- read_csv(enh_da_file, show_col_types = FALSE)

# ---------------- gene-level collapse ----------------
prom_gene <- collapse_peak_to_gene(prom_da, prom_map, alpha = opt$alpha) %>%
  rename_with(~paste0("prom_", .x), -gene_id_clean)

enh_gene <- collapse_peak_to_gene(enh_da, enh_map, alpha = opt$alpha) %>%
  rename_with(~paste0("enh_", .x), -gene_id_clean)

# ---------------- load RNA DEGs (robust columns) ----------------
rna <- read_csv(opt$rna_deg_file, show_col_types = FALSE)

gene_id_col <- intersect(names(rna), c("gene_id","GeneID","ensembl_gene_id","id"))
gene_name_col <- intersect(names(rna), c("gene_name","symbol","GeneName","external_gene_name","gene"))
lfc_col  <- intersect(names(rna), c("log2FoldChange","log2FC","lfc"))
padj_col <- intersect(names(rna), c("padj","adj_pval","FDR","qval","q_value"))

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

# ---------------- merge gene-level ----------------
gene_merged <- rna2 %>%
  left_join(prom_gene, by = "gene_id_clean") %>%
  left_join(enh_gene,  by = "gene_id_clean") %>%
  mutate(
    rna_sig  = !is.na(rna_padj) & rna_padj < opt$alpha,
    prom_sig = !is.na(prom_best_padj) & prom_best_padj < opt$alpha,
    enh_sig  = !is.na(enh_best_padj) & enh_best_padj < opt$alpha,
    sig_class = case_when(
      rna_sig & (prom_sig | enh_sig) ~ "RNA + ATAC",
      rna_sig ~ "RNA only",
      (prom_sig | enh_sig) ~ "ATAC only",
      TRUE ~ "not sig"
    )
  )

write_csv(gene_merged, file.path(opt$datadir, "genelevel_RNA_Prom_Enh_merged.csv"))

top_rna <- gene_merged %>%
  arrange(is.na(rna_padj), rna_padj, desc(abs(rna_log2FC))) %>%
  slice(1:opt$topN)

write_csv(top_rna, file.path(opt$datadir, sprintf("top%d_RNA_genes.csv", opt$topN)))

# ---------------- summary stats ----------------
summ <- tibble(
  n_genes_rna = nrow(rna2),
  n_genes_merged = nrow(gene_merged),
  n_rna_sig = sum(gene_merged$rna_sig, na.rm = TRUE),
  n_prom_sig = sum(gene_merged$prom_sig, na.rm = TRUE),
  n_enh_sig = sum(gene_merged$enh_sig, na.rm = TRUE),
  n_rna_and_prom = sum(gene_merged$rna_sig & gene_merged$prom_sig, na.rm = TRUE),
  n_rna_and_enh  = sum(gene_merged$rna_sig & gene_merged$enh_sig, na.rm = TRUE),
  n_prom_and_enh = sum(gene_merged$prom_sig & gene_merged$enh_sig, na.rm = TRUE),
  n_all_three    = sum(gene_merged$rna_sig & gene_merged$prom_sig & gene_merged$enh_sig, na.rm = TRUE),
  cor_rna_prom = suppressWarnings(cor(gene_merged$rna_log2FC, gene_merged$prom_best_log2FC, use = "complete.obs")),
  cor_rna_enh  = suppressWarnings(cor(gene_merged$rna_log2FC, gene_merged$enh_best_log2FC,  use = "complete.obs")),
  cor_prom_enh = suppressWarnings(cor(gene_merged$prom_best_log2FC, gene_merged$enh_best_log2FC, use = "complete.obs"))
)
write_csv(summ, file.path(opt$datadir, "summary_stats.csv"))

# ---------------- figures ----------------
label_genes <- unique(top_rna$gene_label)

df_rna_prom <- gene_merged %>% filter(!is.na(rna_log2FC), !is.na(prom_best_log2FC))
df_rna_enh  <- gene_merged %>% filter(!is.na(rna_log2FC), !is.na(enh_best_log2FC))
df_prom_enh <- gene_merged %>% filter(!is.na(prom_best_log2FC), !is.na(enh_best_log2FC))

p1 <- make_panel(df_rna_prom, "rna_log2FC","prom_best_log2FC",
                 "RNA log2FC (ATF5KO vs WT)", "Promoter best-peak log2FC",
                 "RNA vs Promoter", label_genes)

p2 <- make_panel(df_rna_enh, "rna_log2FC","enh_best_log2FC",
                 "RNA log2FC (ATF5KO vs WT)", "Enhancer best-peak log2FC",
                 "RNA vs Enhancer", label_genes)

p3 <- make_panel(df_prom_enh, "prom_best_log2FC","enh_best_log2FC",
                 "Promoter best-peak log2FC", "Enhancer best-peak log2FC",
                 "Promoter vs Enhancer", label_genes)

have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (have_patchwork) {
  p_full <- (p1 + p2 + p3) +
    patchwork::plot_layout(ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  ggsave(file.path(figdir, "Figure1_threepanel_fullrange_RNA_Prom_Enh.png"),
         plot = p_full, width = 15, height = 5, dpi = 300)
} else {
  ggsave(file.path(figdir, "Panel1_RNA_vs_Prom_full.png"), p1, width = 5, height = 5, dpi = 300)
  ggsave(file.path(figdir, "Panel2_RNA_vs_Enh_full.png"),  p2, width = 5, height = 5, dpi = 300)
  ggsave(file.path(figdir, "Panel3_Prom_vs_Enh_full.png"), p3, width = 5, height = 5, dpi = 300)
}

x1 <- qrange(df_rna_prom$rna_log2FC); y1 <- qrange(df_rna_prom$prom_best_log2FC)
x2 <- qrange(df_rna_enh$rna_log2FC);  y2 <- qrange(df_rna_enh$enh_best_log2FC)
x3 <- qrange(df_prom_enh$prom_best_log2FC); y3 <- qrange(df_prom_enh$enh_best_log2FC)

p1z <- make_panel(df_rna_prom, "rna_log2FC","prom_best_log2FC",
                  "RNA log2FC (ATF5KO vs WT)", "Promoter best-peak log2FC",
                  "RNA vs Promoter (DE-scale)", label_genes, xlim = x1, ylim = y1)

p2z <- make_panel(df_rna_enh, "rna_log2FC","enh_best_log2FC",
                  "RNA log2FC (ATF5KO vs WT)", "Enhancer best-peak log2FC",
                  "RNA vs Enhancer (DE-scale)", label_genes, xlim = x2, ylim = y2)

p3z <- make_panel(df_prom_enh, "prom_best_log2FC","enh_best_log2FC",
                  "Promoter best-peak log2FC", "Enhancer best-peak log2FC",
                  "Promoter vs Enhancer (DE-scale)", label_genes, xlim = x3, ylim = y3)

if (have_patchwork) {
  p_full <- (p1z + p2z + p3z) +
    patchwork::plot_layout(ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  ggsave(file.path(figdir, "Figure2_threepanel_DEscale_RNA_Prom_Enh.png"),
         plot = p_full, width = 15, height = 5, dpi = 300)
} else {
  ggsave(file.path(figdir, "Panel1_RNA_vs_Prom_DEscale.png"), p1z, width = 5, height = 5, dpi = 300)
  ggsave(file.path(figdir, "Panel2_RNA_vs_Enh_DEscale.png"),  p2z, width = 5, height = 5, dpi = 300)
  ggsave(file.path(figdir, "Panel3_Prom_vs_Enh_DEscale.png"), p3z, width = 5, height = 5, dpi = 300)
}

msg("DONE. Outputs in: %s", opt$datadir)
msg("Figures in: %s", figdir)
msg("Key outputs:")
msg("  genelevel_RNA_Prom_Enh_merged.csv")
msg("  top%d_RNA_genes.csv", opt$topN)
msg("  summary_stats.csv")
msg("  Figure1_threepanel_fullrange_RNA_Prom_Enh.png")
msg("  Figure2_threepanel_DEscale_RNA_Prom_Enh.png")