# CLAUDE.md - ATAC-seq Pipeline (ATF5 KO Study)

## Overview

ATAC-seq analysis pipeline for MEF samples: ATF5 wild-type vs ATF5 knockout.
Runs on Minerva HPC (LSF scheduler) or locally.

**Pipeline Architecture:**
- Single entry point ([analysis/run.sh](analysis/run.sh)) orchestrates three stages
- [analysis/run.lsf](analysis/run.lsf) wraps run.sh for HPC submission
- Modular stage execution with checkpoint logic (skips existing outputs)
- Optional RNA+ATAC integration when provided with RNA-seq DEGs and GTF

## Quick Start

```bash
# HPC: Submit full pipeline (ATAC + integration + figures)
bsub < analysis/run.lsf

# Local: Run full pipeline
bash analysis/run.sh

# Local: Run figures only (requires existing ATAC outputs)
bash analysis/run.sh --skip-atac

# Enable RNA+ATAC integration (edit run.lsf or pass via CLI)
bash analysis/run.sh --rna-deg /path/to/degs.csv --gtf /path/to/anno.gtf

# Validate outputs without re-running
bash analysis/run.sh --validate

# Preview what would run
bash analysis/run.sh --dry-run
```

## Project Structure

```
ATAC-seq_KO/
  analysis/
    run.sh                    # Main pipeline orchestrator (entry point)
    run.lsf                   # LSF submission wrapper for HPC
    config.sh                 # Centralized configuration (paths, params)
    samples.tsv               # Sample metadata (sample, condition, replicate)
    atacseq_analysis.sh       # Stage 1: ATAC processing (QC → peaks → DESeq2)
    deseq_atac.R              # DESeq2 differential accessibility analysis
    idr_analysis.sh           # IDR analysis for peak reproducibility
    compute_qc_metrics.sh     # Extended QC metrics computation
    da_concordance.R          # Multi-method DA concordance (DESeq2, edgeR, limma)
    gene_integration_prep.sh  # Stage 2a: Peak-to-gene mapping (promoter/enhancer)
    rna_atac_figures.R        # Stage 2b: RNA+ATAC integration and correlation plots
    rna_atac_exploratory.R    # Stage 2c: Exploratory figure panel (hexbin, quadrant, heatmap, etc.)
    driver_inference.R        # Driver inference with tiered TF candidates
    validate_outputs.sh       # Output validation and integrity checks
    logs/                     # Pipeline logs (created on first run)
    BINGS/                    # (optional) Additional analysis scripts
  raw/                        # Reference files (user-provided)
    bowtie2_index/            # Bowtie2 genome index (mm10)
    mm10-blacklist.v2.bed     # ENCODE blacklist regions
    mm10.chrom.sizes          # Chromosome sizes
    gencode.vM25.annotation.gtf  # Gene annotation (for integration)
  results/                    # Final organized outputs
    figures/
      qc/                     # PCA, QC plots
      da/                     # Volcano, MA plots
      integration/            # RNA+ATAC correlation plots
      heatmaps/               # deepTools heatmaps
    tables/                   # CSV/TSV results (DA, merged gene-level)
  data/                       # Input data (optional, for local runs)
  .claude/                    # Claude Code configuration
  .serena/                    # Serena agent state
```

## Pipeline Stages

Orchestrated by [analysis/run.sh](analysis/run.sh) with checkpoint logic to skip completed stages:

| Stage | Script | Outputs | Can Skip? |
|-------|--------|---------|-----------|
| 1. ATAC Processing | [atacseq_analysis.sh](analysis/atacseq_analysis.sh) | BAMs, BigWigs, peaks, DESeq2 | Yes (`--skip-atac`) |
| 1.5. Multi-method DA Concordance | [da_concordance.R](analysis/da_concordance.R) | edgeR + limma + DESeq2 comparison | Yes (`--skip-concordance`) |
| 2. RNA+ATAC Integration | [gene_integration_prep.sh](analysis/gene_integration_prep.sh) + [rna_atac_figures.R](analysis/rna_atac_figures.R) + [rna_atac_exploratory.R](analysis/rna_atac_exploratory.R) | Peak-gene maps, merged tables, correlation + exploratory plots | Yes (`--skip-integrate`) |
| 2.5. Driver Inference | [driver_inference.R](analysis/driver_inference.R) | Tiered TF driver candidates | Yes (`--skip-drivers`) |
| 3. Figure Generation | inline in [run.sh](analysis/run.sh) | Volcano, MA, PCA, heatmaps | Yes (`--skip-figures`) |

**Checkpoint behavior:** Existing outputs are automatically skipped unless `--force` is used.

## Configuration

**Main config:** [analysis/config.sh](analysis/config.sh)
- `ATAC_BASE` - Project root (`/sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO`)
- `ATAC_RAW_FASTQ` - FASTQ location
- `ATAC_REFERENCE_DIR` - Reference files (blacklist, bowtie2 index, chrom.sizes)
- `ATAC_THREADS` - Thread count (auto-detected from LSF or defaults to 8)
- `MACS_BIN` - MACS3 path (`/hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3`)
- `ATAC_CONDA_ENV` - Conda environment (`csb`)
- `ATAC_IDR_ENABLED` - Enable IDR peak reproducibility analysis (default: false)
- `ATAC_IDR_THRESHOLD` - IDR threshold for consensus peaks (default: 0.05)
- `ATAC_IDR_MODE` - IDR filtering mode: strict or permissive (default: strict)
- `ATAC_QC_FAIL_POLICY` - QC failure policy: warn or halt (default: warn)
- `ATAC_QC_FRIP_MIN` - Minimum FRiP threshold (default: 0.15)
- `ATAC_QC_TSS_ENRICH_MIN` - Minimum TSS enrichment threshold (default: 5.0)
- `ATAC_DISTANCE_DECAY` - Peak-to-gene distance decay function: none, linear, exponential (default: exponential)
- `ATAC_MAX_DISTAL_GENES` - Max genes linked per distal peak (default: 3)
- `ATAC_CONCORDANCE_ENABLED` - Enable multi-method DA concordance (default: false)
- `ATAC_COVARIATE` - Optional design covariate column (default: none)
- `ATAC_DRIVER_INFERENCE` - Enable driver inference (default: false)

**HPC submission:** [analysis/run.lsf](analysis/run.lsf)
- Edit `SKIP_ATAC`, `SKIP_INTEGRATE`, `SKIP_FIGURES`, `FORCE` flags
- Set `RNA_DEG_FILE` and `GTF_FILE` to enable integration
- LSF resources: 8 cores, 24h, 4GB/core

## Key Outputs

### ATAC-seq Results

| Output | Location |
|--------|----------|
| Final BAMs | `analysis/samtools/*.final.bam` |
| BigWigs (CPM) | `analysis/deeptools/*.CPM.bw` |
| Union peaks | `analysis/macs/union_peaks.bed` |
| DA results (all) | `analysis/subread/deseq/DA_results_DESeq2.csv` |
| DA up (KO) | `analysis/subread/deseq/DA_KO_up.csv` |
| DA down (KO) | `analysis/subread/deseq/DA_KO_down.csv` |
| MultiQC | `analysis/multiqc/multiqc_report.html` |

### IDR (in analysis/macs/idr/)

| Output | Location |
|--------|----------|
| IDR summary | `analysis/macs/idr/idr_summary.tsv` |
| IDR consensus peaks | `analysis/macs/idr/idr_consensus_peaks.bed` |
| IDR method documentation | `analysis/macs/idr/idr_method.txt` |

### Extended QC

| Output | Location |
|--------|----------|
| Extended QC metrics | `analysis/subread/qc_metrics_extended.tsv` |

### Concordance (in analysis/subread/concordance/)

| Output | Location |
|--------|----------|
| Method comparison | `analysis/subread/concordance/DA_comparison.csv` |
| High-confidence peaks | `analysis/subread/concordance/DA_high_confidence.csv` |
| Concordance summary | `analysis/subread/concordance/concordance_summary.txt` |

### Figures (in results/)

| Figure | Location |
|--------|----------|
| Volcano plot | `results/figures/da/volcano_ATAC_DA.png` (labeled with top 15 up/down genes) |
| MA plot | `results/figures/da/MA_plot_ATAC_DA.png` (labeled with top 15 up/down genes) |
| PCA plot | `results/figures/qc/PCA_samples.pdf` |
| Heatmap | `results/figures/heatmaps/heatmap_all_peaks.pdf` |

**Labeling behavior:** If gene integration has been run, volcano/MA plots label peaks with their nearest gene name. Otherwise, peaks are labeled by their genomic coordinates (e.g., `chr1:1000`).

### Integration (optional)

| Output | Location |
|--------|----------|
| Gene-level merge | `analysis/gene_integration/genelevel_RNA_Prom_Enh_merged.csv` |
| Integrated ranking (all) | `analysis/gene_integration/integrated_genes_ranked.csv` |
| Top integrated genes | `analysis/gene_integration/top_integrated_genes.csv` |
| Correlation plots | `results/figures/integration/Figure{1,2}_threepanel_*.png` |
| Top DEGs | `analysis/gene_integration/top30_RNA_genes.csv` |

**Integrated ranking behavior:** genes significant in RNA and at least one ATAC layer are ranked by multi-modal rank product (RNA/promoter/enhancer) with harmonic-mean p-value as secondary evidence.

### Driver Inference

| Output | Location |
|--------|----------|
| Driver candidates | `analysis/gene_integration/driver_candidates_tiered.csv` |
| Evidence summary | `analysis/gene_integration/driver_evidence_summary.txt` |

**Exploratory figures** (from [rna_atac_exploratory.R](analysis/rna_atac_exploratory.R)):

| Figure | Description |
|--------|-------------|
| `Fig1_hexbin_density.png` | Hexbin density plots (avoids overplotting) |
| `Fig2_quadrant_concordance.png` | Concordance analysis with quadrant counts |
| `Fig3a_upset.png` | UpSet plot of sig gene overlaps |
| `Fig3b_venn.png` | Venn diagram of sig gene overlaps |
| `Fig4_focused_scatter.png` | Scatter of sig genes only, labeled by effect size |
| `Fig5_stratified_boxplots.png` | ATAC changes by RNA direction |
| `Fig6_top_DEG_heatmap.png` | Heatmap of top 50 DEGs (RNA, Promoter, Enhancer) |
| `Fig7_prom_vs_enh_barbell.png` | Promoter vs Enhancer comparison per gene |

## Samples

| Sample | Condition | Replicate |
|--------|-----------|-----------|
| MEF-ATF5WT-n2 | WT | 2 |
| MEF-ATF5WT-n3 | WT | 3 |
| MEF-ATF5WT-n4 | WT | 4 |
| MEF-ATF5NULL-n2 | KO | 2 |
| MEF-ATF5NULL-n3 | KO | 3 |
| MEF-ATF5NULL-n4 | KO | 4 |

See [analysis/samples.tsv](analysis/samples.tsv) for metadata.

## Environment

```bash
# HPC environment
module load anaconda
conda activate csb
module load R/4.2.0  # for DESeq2 and integration

# csb environment contains:
# - fastqc, trim_galore, bowtie2, samtools, bedtools
# - deeptools, subread (featureCounts), multiqc, picard
```

**R packages required:**
- Core: DESeq2, apeglm, readr, dplyr, ggplot2, ggrepel, optparse, patchwork
- Concordance: edgeR, limma
- Exploratory figures: pheatmap, UpSetR, ggVennDiagram, tidyr

## Reproducibility Features

The pipeline ensures reproducibility through:

1. **Version tracking:** Software versions recorded in `analysis/logs/versions.txt`
2. **Checkpoint logic:** Outputs are timestamped and never overwritten unless `--force` is used
3. **Configuration as code:** All parameters centralized in [analysis/config.sh](analysis/config.sh)
4. **Sample metadata:** Explicit sample-to-condition mapping in [analysis/samples.tsv](analysis/samples.tsv)
5. **Deterministic processing:** Fixed parameters (MAPQ thresholds, genome size, peak calling settings)
6. **Run provenance:** LSF job IDs and timestamps in all log files
7. **Modular stages:** Can re-run individual stages without full pipeline re-execution

## Pipeline Options

```
--skip-atac          Skip ATAC processing (use existing BAMs/peaks)
--skip-concordance   Skip multi-method DA concordance
--skip-integrate     Skip RNA+ATAC integration
--skip-drivers       Skip driver inference
--skip-figures       Skip figure generation
--force              Re-run all steps (ignore checkpoints)
--rna-deg FILE       RNA-seq DEG file for integration
--gtf FILE           GTF annotation for integration
--covariate COL      Include covariate in design formula
--validate           Run validation only (no execution)
--dry-run            Preview what would run without executing
```

## Current Project State

**Status:** Pipeline ready for execution (not yet run on full dataset)

**Git History:**
- `2f3e666` - dir fix
- `240e8a1` - feat: initial ATAC-seq pipeline setup

**What's Configured:**
- ✅ Project structure established
- ✅ All scripts written and validated
- ✅ Configuration file set for Minerva HPC
- ✅ RNA integration paths configured in run.lsf
- ✅ Sample metadata in [analysis/samples.tsv](analysis/samples.tsv)
- ✅ LSF submission wrapper ready

**Next Steps:**
1. Ensure reference files exist in `raw/` (bowtie2 index, blacklist, chrom.sizes, GTF)
2. Create logs directory: `mkdir -p analysis/logs`
3. Submit pipeline: `bsub < analysis/run.lsf`
4. Monitor progress: `bjobs -w` and `tail -f analysis/logs/*.pipeline.log`

**Integration Status:**
- RNA-seq DEG file path configured in run.lsf
- GTF annotation path configured in run.lsf
- Integration will run automatically when pipeline is submitted

## Troubleshooting

### Check Pipeline Status

```bash
# Monitor LSF job
bjobs -w

# View pipeline log (most informative)
tail -f analysis/logs/*.pipeline.log

# View LSF stdout/stderr
tail -f analysis/logs/*.out
tail -f analysis/logs/*.err

# Search for errors
grep -i error analysis/logs/*.log
```

### Common Issues

**Issue: `logs/` directory doesn't exist**
```bash
mkdir -p analysis/logs
```

**Issue: Reference files not found**
```bash
# Verify all required files exist
ls raw/bowtie2_index/index.*.bt2
ls raw/mm10-blacklist.v2.bed
ls raw/mm10.chrom.sizes
ls raw/gencode.vM25.annotation.gtf  # only needed for integration
```

**Issue: Conda environment `csb` not found**
```bash
# Create environment with required tools
conda create -n csb -c bioconda -c conda-forge \
  fastqc trim-galore bowtie2 samtools bedtools \
  deeptools subread multiqc picard
```

**Issue: MACS3 not found**
```bash
# Update path in config.sh
ls -l /hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3
# If different, edit MACS_BIN in analysis/config.sh
```

**Issue: R packages missing**
```bash
module load R/4.2.0
R
> if (!require("BiocManager")) install.packages("BiocManager")
> BiocManager::install(c("DESeq2", "apeglm", "edgeR", "limma"))
> install.packages(c("readr", "dplyr", "ggplot2", "stringr", "tibble", "ggrepel", "optparse", "patchwork"))

# For exploratory figures (rna_atac_exploratory.R)
> install.packages(c("pheatmap", "UpSetR", "ggVennDiagram", "tidyr"))
```

### Re-running Stages

The pipeline uses checkpoint logic—existing outputs are skipped unless `--force` is used.

**Re-run figures only:**
```bash
bash analysis/run.sh --skip-atac --skip-integrate
```

**Force re-run DESeq2 and downstream:**
```bash
rm analysis/subread/deseq/DA_results_DESeq2.csv
bash analysis/run.sh --skip-atac --skip-integrate
```

**Full re-run from scratch:**
```bash
bash analysis/run.sh --force
```

## Related Documentation

- Full pipeline details: [README.md](README.md)
- Configuration reference: [analysis/config.sh](analysis/config.sh)
- Sample metadata format: [analysis/samples.tsv](analysis/samples.tsv)
