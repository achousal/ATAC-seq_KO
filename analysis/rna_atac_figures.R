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

# Collapse peak-level DA to gene-level using signed Stouffer's Z
collapse_peak_to_gene <- function(da_df, map_df, alpha = 0.05) {
  # Validate raw pvalue column
  if (!"pvalue" %in% colnames(da_df)) {
    # Try fallback columns
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

modality_rank <- function(padj, log2fc) {
  p <- as.numeric(padj)
  fc <- as.numeric(log2fc)

  p[!is.finite(p)] <- 1
  fc[!is.finite(fc)] <- 0
  p <- pmin(pmax(p, 1e-300), 1)

  ord <- order(p, -abs(fc))
  ranks <- integer(length(p))
  ranks[ord] <- seq_along(ord)
  ranks
}

harmonic_mean_p <- function(pvals) {
  p <- as.numeric(pvals)
  p <- p[is.finite(p)]
  if (length(p) == 0) return(NA_real_)
  p <- pmin(pmax(p, 1e-300), 1)
  length(p) / sum(1 / p)
}

select_panel_labels <- function(df, sig_cols, padj_cols, lfc_cols, top_n = 30) {
  sig_cols <- intersect(sig_cols, names(df))
  padj_cols <- intersect(padj_cols, names(df))
  lfc_cols <- intersect(lfc_cols, names(df))
  if (!("gene_label" %in% names(df)) || length(sig_cols) == 0 || length(padj_cols) == 0 || length(lfc_cols) == 0) {
    return(character(0))
  }

  n <- nrow(df)
  any_sig <- rep(FALSE, n)
  both_sig <- rep(TRUE, n)
  for (col in sig_cols) {
    v <- as.logical(df[[col]])
    v[is.na(v)] <- FALSE
    any_sig <- any_sig | v
    both_sig <- both_sig & v
  }

  padj_mat <- sapply(padj_cols, function(col) {
    v <- as.numeric(df[[col]])
    v[is.na(v)] <- 1
    v
  })
  if (is.null(dim(padj_mat))) padj_mat <- matrix(padj_mat, ncol = 1)
  best_padj <- apply(padj_mat, 1, min)

  lfc_mat <- sapply(lfc_cols, function(col) abs(as.numeric(df[[col]])))
  if (is.null(dim(lfc_mat))) lfc_mat <- matrix(lfc_mat, ncol = 1)
  lfc_mat[!is.finite(lfc_mat)] <- NA_real_
  max_abs_lfc <- apply(lfc_mat, 1, function(x) {
    m <- suppressWarnings(max(x, na.rm = TRUE))
    if (is.finite(m)) m else 0
  })

  rank_df <- tibble(
    gene_label = as.character(df$gene_label),
    any_sig = any_sig,
    both_sig = both_sig,
    best_padj = best_padj,
    max_abs_lfc = max_abs_lfc
  ) %>%
    filter(!is.na(gene_label), nzchar(gene_label)) %>%
    arrange(desc(any_sig), desc(both_sig), best_padj, desc(max_abs_lfc)) %>%
    distinct(gene_label, .keep_all = TRUE)

  rank_df %>% slice_head(n = top_n) %>% pull(gene_label)
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
    ggplot2::scale_fill_manual(values = pal, guide = "none", drop = FALSE) +
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
    lab_df <- df %>% dplyr::filter(
      !is.na(gene_label),
      gene_label %in% label_genes,
      sig_class != "not sig"
    )

    p <- p + ggplot2::geom_point(
      data = lab_df,
      ggplot2::aes(fill = sig_class),
      shape = 21,
      size = 2.2,
      stroke = 0.3,
      color = "black",
      alpha = 0.95,
      show.legend = FALSE
    )

    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = gene_label),
      size = 3.1,
      fontface = "italic",
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
# Auto-detect 7-col (weighted) vs 6-col (legacy) map files
load_map <- function(path) {
  tmp <- read_tsv(path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(tmp) >= 7) {
    colnames(tmp) <- c("peak_id","gene_id","gene_name","strand","dist_bp","abs_dist","weight")[seq_len(ncol(tmp))]
  } else {
    colnames(tmp) <- c("peak_id","gene_id","gene_name","strand","dist_bp","abs_dist")[seq_len(ncol(tmp))]
  }
  tmp
}
prom_map <- load_map(prom_map_file)
enh_map  <- load_map(enh_map_file)

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
pval_raw_col <- intersect(names(rna), c("pvalue", "pval", "p_value"))

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
    rna_padj = as.numeric(.data[[padj_col]]),
    rna_pvalue = if (length(pval_raw_col) > 0) {
      as.numeric(.data[[pval_raw_col[1]]])
    } else {
      warning("RNA DEGs missing raw pvalue column; using padj as fallback for HMP (conservative)")
      as.numeric(.data[[padj_col]])
    }
  ) %>%
  distinct(gene_id_clean, .keep_all = TRUE)

# ---------------- merge gene-level ----------------
gene_merged <- rna2 %>%
  left_join(prom_gene, by = "gene_id_clean") %>%
  left_join(enh_gene,  by = "gene_id_clean") %>%
  mutate(
    rna_sig  = !is.na(rna_padj) & rna_padj < opt$alpha,
    prom_sig = !is.na(prom_combined_padj) & prom_combined_padj < opt$alpha,
    enh_sig  = !is.na(enh_combined_padj) & enh_combined_padj < opt$alpha,
    sig_class = case_when(
      rna_sig & (prom_sig | enh_sig) ~ "RNA + ATAC",
      rna_sig ~ "RNA only",
      (prom_sig | enh_sig) ~ "ATAC only",
      TRUE ~ "not sig"
    )
  )

# Modality-specific ranks for robust rank aggregation
gene_merged <- gene_merged %>%
  mutate(
    rna_rank = modality_rank(rna_padj, rna_log2FC),
    prom_rank = modality_rank(prom_combined_padj, prom_weighted_log2FC),
    enh_rank = modality_rank(enh_combined_padj, enh_weighted_log2FC)
  )

write_csv(gene_merged, file.path(opt$datadir, "genelevel_RNA_Prom_Enh_merged.csv"))

top_rna <- gene_merged %>%
  arrange(is.na(rna_padj), rna_padj, desc(abs(rna_log2FC))) %>%
  slice(1:opt$topN)

write_csv(top_rna, file.path(opt$datadir, sprintf("top%d_RNA_genes.csv", opt$topN)))

n_genes_total <- nrow(gene_merged)

integrated_ranked <- gene_merged %>%
  filter(rna_sig & (prom_sig | enh_sig)) %>%
  rowwise() %>%
  mutate(
    n_sig_modalities = 1L + as.integer(prom_sig) + as.integer(enh_sig),
    rank_product = exp(mean(log(c(
      rna_rank,
      if (isTRUE(prom_sig)) prom_rank else NULL,
      if (isTRUE(enh_sig)) enh_rank else NULL
    )))),
    rank_product_norm = rank_product / n_genes_total,
    integrated_rank_score = -log10(pmax(rank_product_norm, 1e-300)),
    hmp_raw = harmonic_mean_p(c(
      rna_pvalue,
      if (isTRUE(prom_sig)) prom_combined_pvalue else NULL,
      if (isTRUE(enh_sig)) enh_combined_pvalue else NULL
    )),
    hmp_neglog10 = -log10(pmax(hmp_raw, 1e-300)),
    atac_best_layer = case_when(
      isTRUE(prom_sig) & isTRUE(enh_sig) ~ ifelse(
        ifelse(is.na(prom_combined_padj), 1, prom_combined_padj) <= ifelse(is.na(enh_combined_padj), 1, enh_combined_padj),
        "promoter",
        "enhancer"
      ),
      isTRUE(prom_sig) ~ "promoter",
      isTRUE(enh_sig) ~ "enhancer",
      TRUE ~ NA_character_
    ),
    atac_best_log2FC = case_when(
      atac_best_layer == "promoter" ~ prom_best_log2FC,
      atac_best_layer == "enhancer" ~ enh_best_log2FC,
      TRUE ~ NA_real_
    ),
    direction_relation = case_when(
      is.na(rna_log2FC) | is.na(atac_best_log2FC) ~ NA_character_,
      rna_log2FC == 0 | atac_best_log2FC == 0 ~ "ambiguous",
      sign(rna_log2FC) == sign(atac_best_log2FC) ~ "concordant",
      TRUE ~ "discordant"
    )
  ) %>%
  ungroup() %>%
  mutate(hmp_evidence_score = p.adjust(hmp_raw, method = "BH")) %>%
  arrange(desc(n_sig_modalities), rank_product, hmp_evidence_score, desc(abs(rna_log2FC))) %>%
  mutate(integrated_rank = row_number())

# Write integrated ranking with header comment
integrated_csv_path <- file.path(opt$datadir, "integrated_genes_ranked.csv")
cat("# Integrated gene ranking (RNA+ATAC)\n", file = integrated_csv_path)
cat("# hmp_evidence_score: BH-adjusted harmonic mean p-value (evidence ranking, not formal significance)\n", file = integrated_csv_path, append = TRUE)
cat("# Formal significance calls: rna_sig, prom_sig, enh_sig (based on modality-specific padj)\n", file = integrated_csv_path, append = TRUE)
cat("# Combined p-values: prom_combined_padj, enh_combined_padj (Stouffer's Z on raw p-values, then BH)\n", file = integrated_csv_path, append = TRUE)
write_csv(integrated_ranked, integrated_csv_path, append = TRUE)

top_integrated <- integrated_ranked %>%
  slice_head(n = opt$topN)

write_csv(top_integrated, file.path(opt$datadir, "top_integrated_genes.csv"))
write_csv(top_integrated, file.path(opt$datadir, sprintf("top%d_integrated_genes.csv", opt$topN)))

# ---------------- summary stats ----------------
summ <- tibble(
  n_genes_rna = nrow(rna2),
  n_genes_merged = nrow(gene_merged),
  n_integrated_ranked = nrow(integrated_ranked),
  n_integrated_concordant = sum(integrated_ranked$direction_relation == "concordant", na.rm = TRUE),
  n_integrated_discordant = sum(integrated_ranked$direction_relation == "discordant", na.rm = TRUE),
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
df_rna_prom <- gene_merged %>% filter(!is.na(rna_log2FC), !is.na(prom_best_log2FC))
df_rna_enh  <- gene_merged %>% filter(!is.na(rna_log2FC), !is.na(enh_best_log2FC))
df_prom_enh <- gene_merged %>% filter(!is.na(prom_best_log2FC), !is.na(enh_best_log2FC))

label_genes_rna_prom <- select_panel_labels(
  df_rna_prom,
  sig_cols = c("rna_sig", "prom_sig"),
  padj_cols = c("rna_padj", "prom_combined_padj"),
  lfc_cols = c("rna_log2FC", "prom_best_log2FC"),
  top_n = opt$topN
)

label_genes_rna_enh <- select_panel_labels(
  df_rna_enh,
  sig_cols = c("rna_sig", "enh_sig"),
  padj_cols = c("rna_padj", "enh_combined_padj"),
  lfc_cols = c("rna_log2FC", "enh_best_log2FC"),
  top_n = opt$topN
)

label_genes_prom_enh <- select_panel_labels(
  df_prom_enh,
  sig_cols = c("prom_sig", "enh_sig"),
  padj_cols = c("prom_combined_padj", "enh_combined_padj"),
  lfc_cols = c("prom_best_log2FC", "enh_best_log2FC"),
  top_n = opt$topN
)

p1 <- make_panel(df_rna_prom, "rna_log2FC","prom_best_log2FC",
                 "RNA log2FC (ATF5KO vs WT)", "Promoter best-peak log2FC",
                 "RNA vs Promoter", label_genes_rna_prom)

p2 <- make_panel(df_rna_enh, "rna_log2FC","enh_best_log2FC",
                 "RNA log2FC (ATF5KO vs WT)", "Enhancer best-peak log2FC",
                 "RNA vs Enhancer", label_genes_rna_enh)

p3 <- make_panel(df_prom_enh, "prom_best_log2FC","enh_best_log2FC",
                 "Promoter best-peak log2FC", "Enhancer best-peak log2FC",
                 "Promoter vs Enhancer", label_genes_prom_enh)

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

msg("DONE. Outputs in: %s", opt$datadir)
msg("Figures in: %s", figdir)
msg("Key outputs:")
msg("  genelevel_RNA_Prom_Enh_merged.csv")
msg("  top%d_RNA_genes.csv", opt$topN)
msg("  integrated_genes_ranked.csv (RNA+ATAC rank product)")
msg("  top_integrated_genes.csv")
msg("  top%d_integrated_genes.csv", opt$topN)
msg("  summary_stats.csv")
msg("  Figure1_threepanel_fullrange_RNA_Prom_Enh.png")
