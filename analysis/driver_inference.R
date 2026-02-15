#!/usr/bin/env Rscript
# driver_inference.R
# Minimal driver inference: integrates ATAC DA, motif enrichment, and RNA DE
# This is hypothesis-generating evidence integration, not causal proof of driver status.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Logging helper
msg <- function(...) {
  message(sprintf("[%s] %s", Sys.time(), paste0(...)))
}

# Package check helper
check_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
}

# Strip Ensembl version suffixes
clean_ens <- function(x) {
  gsub("\\.\\d+$", "", x)
}

# Parse HOMER knownResults.txt
parse_homer <- function(homer_dir, direction, pval_thresh) {
  known_results <- file.path(homer_dir, "knownResults.txt")

  if (!file.exists(known_results)) {
    msg(sprintf("Warning: HOMER results not found at %s", known_results))
    return(tibble(
      TF = character(0),
      Direction = character(0),
      Motif_Pval = numeric(0),
      Motif_Qval = numeric(0),
      Full_Motif_Name = character(0)
    ))
  }

  msg(sprintf("Parsing HOMER results from %s", known_results))

  # HOMER knownResults.txt is tab-delimited with header
  # Col 1: Motif Name, Col 3: P-value, Col 5: q-value
  homer <- read_tsv(known_results, show_col_types = FALSE)
  colnames(homer) <- c("MotifName", "Consensus", "Pval", "LogPval", "Qval",
                       "NumTargetSeq", "PctTargetSeq", "NumBgSeq", "PctBgSeq")

  # Extract TF name (text before first parenthesis)
  homer <- homer %>%
    mutate(
      TF = str_extract(MotifName, "^[^(]+"),
      TF = str_trim(TF),
      Direction = direction,
      Motif_Pval = as.numeric(Pval),
      Motif_Qval = as.numeric(Qval),
      Full_Motif_Name = MotifName
    ) %>%
    filter(Motif_Pval <= pval_thresh) %>%
    select(TF, Direction, Motif_Pval, Motif_Qval, Full_Motif_Name)

  msg(sprintf("  Found %d enriched motifs (p <= %.3f)", nrow(homer), pval_thresh))

  return(homer)
}

# Parse RNA-seq DEGs with flexible column detection
parse_rna_deg <- function(rna_file, alpha) {
  msg(sprintf("Parsing RNA-seq DEGs from %s", rna_file))

  rna <- read_csv(rna_file, show_col_types = FALSE)

  # Detect gene ID column
  gene_id_col <- intersect(c("gene_id", "GeneID", "ensembl_gene_id"), colnames(rna))
  if (length(gene_id_col) == 0) {
    stop("Could not find gene ID column (expected: gene_id, GeneID, or ensembl_gene_id)")
  }
  gene_id_col <- gene_id_col[1]

  # Detect gene name column
  gene_name_col <- intersect(c("gene_name", "symbol", "gene_symbol", "GeneName"), colnames(rna))
  if (length(gene_name_col) == 0) {
    msg("Warning: No gene name column found; using gene ID only")
    gene_name_col <- NULL
  } else {
    gene_name_col <- gene_name_col[1]
  }

  # Detect log2FC column
  lfc_col <- intersect(c("log2FoldChange", "log2FC", "logFC"), colnames(rna))
  if (length(lfc_col) == 0) {
    stop("Could not find log2FC column (expected: log2FoldChange, log2FC, or logFC)")
  }
  lfc_col <- lfc_col[1]

  # Detect adjusted p-value column
  padj_col <- intersect(c("padj", "FDR", "adj.P.Val", "qvalue"), colnames(rna))
  if (length(padj_col) == 0) {
    stop("Could not find adjusted p-value column (expected: padj, FDR, adj.P.Val, or qvalue)")
  }
  padj_col <- padj_col[1]

  # Build output
  rna_out <- rna %>%
    mutate(
      gene_id_clean = clean_ens(.data[[gene_id_col]]),
      log2FC = as.numeric(.data[[lfc_col]]),
      padj = as.numeric(.data[[padj_col]]),
      RNA_DE = !is.na(padj) & padj <= alpha
    )

  if (!is.null(gene_name_col)) {
    rna_out <- rna_out %>%
      mutate(TF = .data[[gene_name_col]]) %>%
      select(gene_id_clean, TF, log2FC, padj, RNA_DE)
  } else {
    rna_out <- rna_out %>%
      mutate(TF = gene_id_clean) %>%
      select(gene_id_clean, TF, log2FC, padj, RNA_DE)
  }

  msg(sprintf("  Found %d DE genes (padj <= %.3f)", sum(rna_out$RNA_DE, na.rm = TRUE), alpha))

  return(rna_out)
}

# Integrate evidence
integrate_evidence <- function(motif_up, motif_down, rna_deg, ranked_file) {
  msg("Integrating evidence across ATAC motifs and RNA DE")

  # Parse ranked peaks for ATAC DA evidence
  ranked <- read_csv(ranked_file, show_col_types = FALSE)

  # Detect nearest_gene column
  gene_col <- intersect(c("nearest_gene", "gene_name", "symbol"), colnames(ranked))
  if (length(gene_col) == 0) {
    msg("Warning: No nearest_gene column in ranked file; skipping ATAC DA evidence")
    atac_genes <- character(0)
  } else {
    gene_col <- gene_col[1]
    atac_genes <- unique(ranked[[gene_col]])
  }

  # Combine motif evidence from both directions
  motif_all <- bind_rows(motif_up, motif_down) %>%
    group_by(TF) %>%
    summarize(
      Motif_Directions = paste(unique(Direction), collapse = ";"),
      Best_Motif_Pval = min(Motif_Pval, na.rm = TRUE),
      Best_Motif_Qval = min(Motif_Qval, na.rm = TRUE),
      Full_Motif_Names = paste(unique(Full_Motif_Name), collapse = " | "),
      .groups = "drop"
    ) %>%
    mutate(Motif_Evidence = TRUE)

  # Merge with RNA DE
  integrated <- rna_deg %>%
    full_join(motif_all, by = "TF") %>%
    replace_na(list(
      Motif_Evidence = FALSE,
      RNA_DE = FALSE
    )) %>%
    mutate(
      ATAC_DA = TF %in% atac_genes,
      Evidence_Count = as.integer(ATAC_DA) + as.integer(Motif_Evidence) + as.integer(RNA_DE),
      Tier = case_when(
        Evidence_Count == 3 ~ "Tier 1: ATAC DA + Motif + RNA DE",
        Evidence_Count == 2 ~ "Tier 2: Two evidence lines",
        Evidence_Count == 1 ~ "Tier 3: Single evidence (exploratory)",
        TRUE ~ "No evidence"
      )
    ) %>%
    filter(Evidence_Count > 0) %>%
    arrange(desc(Evidence_Count), Best_Motif_Pval, padj)

  msg(sprintf("  Tier 1 candidates: %d", sum(integrated$Tier == "Tier 1: ATAC DA + Motif + RNA DE")))
  msg(sprintf("  Tier 2 candidates: %d", sum(integrated$Tier == "Tier 2: Two evidence lines")))
  msg(sprintf("  Tier 3 candidates: %d", sum(integrated$Tier == "Tier 3: Single evidence (exploratory)")))

  return(integrated)
}

# Write summary
write_summary <- function(integrated, outfile, motif_pval, alpha) {
  msg(sprintf("Writing evidence summary to %s", outfile))

  sink(outfile)
  cat("Driver Inference Evidence Summary\n")
  cat("=====================================\n\n")
  cat("This is hypothesis-generating evidence integration, not causal proof of driver status.\n\n")
  cat(sprintf("Parameters:\n"))
  cat(sprintf("  - Motif p-value threshold: %.3f\n", motif_pval))
  cat(sprintf("  - RNA DE significance threshold: %.3f\n", alpha))
  cat("\n")

  cat("Evidence Tiers:\n")
  cat("  - Tier 1: ATAC DA + Motif Enrichment + RNA DE (strongest candidates)\n")
  cat("  - Tier 2: Two of three evidence lines\n")
  cat("  - Tier 3: Single evidence line (exploratory)\n\n")

  cat("Results:\n")
  tier_counts <- integrated %>%
    count(Tier) %>%
    arrange(desc(n))

  for (i in seq_len(nrow(tier_counts))) {
    cat(sprintf("  %s: %d candidates\n", tier_counts$Tier[i], tier_counts$n[i]))
  }

  cat("\n")
  cat("Top 10 Tier 1 Candidates (if any):\n")
  tier1 <- integrated %>%
    filter(Tier == "Tier 1: ATAC DA + Motif + RNA DE") %>%
    head(10)

  if (nrow(tier1) > 0) {
    for (i in seq_len(nrow(tier1))) {
      cat(sprintf("  %d. %s (RNA log2FC: %.2f, Motif p: %.2e, RNA padj: %.2e)\n",
                  i, tier1$TF[i], tier1$log2FC[i], tier1$Best_Motif_Pval[i], tier1$padj[i]))
    }
  } else {
    cat("  (No Tier 1 candidates found)\n")
  }

  sink()
}

# Main function
main <- function() {
  option_list <- list(
    make_option(c("--ranked"), type = "character", default = NULL,
                help = "Integrated ranked genes CSV", metavar = "FILE"),
    make_option(c("--homer_up"), type = "character", default = NULL,
                help = "HOMER results directory for KO-up peaks", metavar = "DIR"),
    make_option(c("--homer_down"), type = "character", default = NULL,
                help = "HOMER results directory for KO-down peaks", metavar = "DIR"),
    make_option(c("--rna_deg"), type = "character", default = NULL,
                help = "RNA-seq DEG CSV", metavar = "FILE"),
    make_option(c("--outdir"), type = "character", default = NULL,
                help = "Output directory", metavar = "DIR"),
    make_option(c("--motif_pval"), type = "double", default = 0.05,
                help = "Motif p-value threshold [default: %default]", metavar = "NUM"),
    make_option(c("--alpha"), type = "double", default = 0.05,
                help = "RNA DE significance threshold [default: %default]", metavar = "NUM")
  )

  parser <- OptionParser(
    usage = "%prog --ranked FILE --homer_up DIR --homer_down DIR --rna_deg FILE --outdir DIR",
    option_list = option_list,
    description = "\nMinimal driver inference: integrates ATAC DA, motif enrichment, and RNA DE.\nThis is hypothesis-generating evidence integration, not causal proof of driver status."
  )

  opt <- parse_args(parser)

  # Validate required arguments
  required_args <- c("ranked", "homer_up", "homer_down", "rna_deg", "outdir")
  missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

  if (length(missing_args) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")))
  }

  # Validate input files
  if (!file.exists(opt$ranked)) {
    stop(sprintf("Ranked genes file not found: %s", opt$ranked))
  }
  if (!file.exists(opt$rna_deg)) {
    stop(sprintf("RNA DEG file not found: %s", opt$rna_deg))
  }
  if (!dir.exists(opt$homer_up)) {
    msg(sprintf("Warning: HOMER up directory not found: %s", opt$homer_up))
  }
  if (!dir.exists(opt$homer_down)) {
    msg(sprintf("Warning: HOMER down directory not found: %s", opt$homer_down))
  }

  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
    msg(sprintf("Created output directory: %s", opt$outdir))
  }

  msg("Starting driver inference")
  msg(sprintf("  Ranked genes: %s", opt$ranked))
  msg(sprintf("  HOMER up: %s", opt$homer_up))
  msg(sprintf("  HOMER down: %s", opt$homer_down))
  msg(sprintf("  RNA DEG: %s", opt$rna_deg))
  msg(sprintf("  Output dir: %s", opt$outdir))
  msg(sprintf("  Motif p-value threshold: %.3f", opt$motif_pval))
  msg(sprintf("  RNA DE alpha: %.3f", opt$alpha))

  # Parse inputs
  motif_up <- parse_homer(opt$homer_up, "KO_up", opt$motif_pval)
  motif_down <- parse_homer(opt$homer_down, "KO_down", opt$motif_pval)
  rna_deg <- parse_rna_deg(opt$rna_deg, opt$alpha)

  # Integrate evidence
  integrated <- integrate_evidence(motif_up, motif_down, rna_deg, opt$ranked)

  # Write outputs
  out_csv <- file.path(opt$outdir, "driver_candidates_tiered.csv")
  out_summary <- file.path(opt$outdir, "driver_evidence_summary.txt")

  write_csv(integrated, out_csv)
  msg(sprintf("Wrote driver candidates to %s", out_csv))

  write_summary(integrated, out_summary, opt$motif_pval, opt$alpha)

  msg("Driver inference complete")
}

# Run main
main()
