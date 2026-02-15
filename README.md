# ATAC-seq Pipeline: ATF5 KO Study

Complete workflow for analyzing ATAC-seq data from MEF (mouse embryonic fibroblast) samples comparing ATF5 wild-type vs ATF5 knockout conditions on Minerva HPC.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Workflow Order of Operations](#workflow-order-of-operations)
- [Detailed Step-by-Step Guide](#detailed-step-by-step-guide)
- [Troubleshooting and Re-runs](#troubleshooting-and-re-runs)
- [Output Files](#output-files)
- [Validation](#validation)

## Prerequisites

### Required Software (on Minerva)

- Anaconda/Miniconda with `csb` environment
- Tools in `csb` environment: fastqc, trim_galore, bowtie2, samtools, bedtools, deepTools, featureCounts, multiqc, picard
- R (4.2.0) with packages: DESeq2, apeglm, ggplot2, dplyr, readr, stringr, tibble
- MACS3 (hardcoded path: `/hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3`)
- HOMER (module: homer/4.10)

### Required Reference Files

Located in `/sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/raw/`:
- `mm10-blacklist.v2.bed`: ENCODE blacklist regions
- `bowtie2_index/index`: Bowtie2 genome index (mm10)
- `mm10.chrom.sizes`: Chromosome sizes file

### Sample Data

Raw FASTQ files (paired-end):
- MEF-ATF5WT-n2, n3, n4 (3 biological replicates, WT condition)
- MEF-ATF5NULL-n2, n3, n4 (3 biological replicates, KO condition)

## Quick Start

```bash
# 1. Navigate to project directory
cd /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO

# 2. Create log directory (required before submission)
mkdir -p analysis/logs

# 3. Submit pipeline to LSF scheduler
bsub < analysis/run_atac.lsf

# 4. Monitor job status
bjobs -w

# 5. View live pipeline log
tail -f analysis/logs/*.pipeline.log
```

## Workflow Order of Operations

### High-Level Overview

1. **Setup and Configuration** (manual)
2. **Job Submission** (manual: `bsub`)
3. **QC and Preprocessing** (automated)
4. **Alignment and Filtering** (automated)
5. **Peak Calling** (automated)
6. **Quantification** (automated)
7. **Differential Analysis** (automated)
8. **Validation** (manual)
9. **RNA-ATAC Integration** (optional, manual)
10. **Integration Visualization** (optional, manual)

### Detailed Execution Flow

```
run_atac.lsf (LSF wrapper)
    |
    +-- Load environment (LSF job configuration)
    +-- Set resource limits (8 cores, 24h, 4GB/core)
    +-- Create logs directory
    |
    v
atacseq_analysis.sh (main pipeline)
    |
    +-- Load config.sh (paths, parameters)
    +-- Activate conda environment (csb)
    +-- Validate reference files
    |
    +-- STEP 1: QC + Preprocessing + Alignment + Filtering + Dedup
    |   |
    |   +-- For each sample:
    |       +-- FastQC (raw reads)
    |       +-- TrimGalore (Nextera adapters)
    |       +-- Bowtie2 (align with --very-sensitive -X 2000)
    |       +-- Samtools (sort, filter: MAPQ>=20, proper pairs, canonical chr only)
    |       +-- Picard MarkDuplicates (REMOVE_DUPLICATES=true)
    |       +-- Samtools (final sort and index)
    |       +-- deepTools bamCoverage (CPM-normalized BigWigs)
    |       +-- Record QC metrics (alignment rate, chrM fraction, dup rate, etc.)
    |
    +-- STEP 1b: MultiQC (aggregate all QC reports)
    |
    +-- STEP 2: Per-replicate Peak Calling
    |   |
    |   +-- For each sample:
    |       +-- MACS3 callpeak (BAMPE mode, --nomodel)
    |       +-- Bedtools intersect (remove blacklisted regions)
    |
    +-- STEP 3: IDR Analysis (if idr available) or Union Peak Set
    |   |
    |   +-- For each condition (WT, KO):
    |   |   +-- Pool replicate BAMs
    |   |   +-- Call peaks on pooled sample
    |   |   +-- Generate two pseudoreplicates (random 50/50 split)
    |   |   +-- Call peaks on each pseudoreplicate
    |   |   +-- Call peaks on true replicates
    |   |   +-- Run IDR on true replicates
    |   |   +-- Run IDR on pooled pseudoreplicates
    |   |   +-- Determine Nt and Np (peaks passing IDR threshold)
    |   |   +-- Select top max(Nt, Np) peaks from pooled peak set
    |   +-- Merge WT-optimal and KO-optimal peaks to form consensus
    |   +-- Fallback to union peak calling if IDR unavailable
    |   +-- Convert to SAF format (for featureCounts)
    |
    +-- STEP 3.5: Extended QC Metrics
    |   |
    |   +-- For each sample:
    |       +-- Compute FRiP (fraction of reads in peaks)
    |       +-- Compute TSS enrichment score
    |       +-- Compute fragment size distribution (mean, median)
    |       +-- Compute nucleosome-free ratio (NFR)
    |       +-- Compute blacklist fraction
    |       +-- Compute library complexity (PBC1, PBC2)
    |   +-- Write qc_metrics_extended.tsv
    |
    +-- STEP 3.6: HOMER Annotation and Motif Analysis
    |   |
    |   +-- Load HOMER module
    |   +-- Center peaks (1bp centers for consistent -size window)
    |   +-- annotatePeaks.pl (annotate union peaks, size=200)
    |   +-- findMotifsGenome.pl (motif enrichment, size=200, len=8,10,12)
    |   +-- makeTagDirectory (create tag directories for each sample)
    |   +-- Unload HOMER module
    |
    +-- STEP 4: Quantification
    |   |
    |   +-- featureCounts (count reads in union peaks, paired-end mode)
    |
    +-- STEP 5: Differential Accessibility Analysis
        |
        +-- Load R/4.2.0 module
        +-- Run deseq_atac.R
            |
            +-- Load count matrix
            +-- Parse sample metadata (from samples.tsv or filenames)
            +-- Create DESeq2 dataset (condition: NULL vs WT)
            +-- Run differential analysis
            +-- Apply lfcShrink (apeglm)
            +-- Generate outputs:
                +-- DA_results_DESeq2.csv (all peaks with stats)
                +-- DA_KO_up.csv (peaks with higher accessibility in KO, padj<=0.05, LFC>=1)
                +-- DA_KO_down.csv (peaks with lower accessibility in KO, padj<=0.05, LFC<=-1)
                +-- PCA_samples.pdf (sample PCA plot, VST-transformed)
                +-- PCA_coordinates.csv (PC1/PC2 coordinates + variance explained)
        +-- Unload R module
    |
    +-- STEP 5.5: Multi-Method Concordance Analysis
        |
        +-- Load R/4.2.0 module
        +-- Run da_concordance.R
            |
            +-- Load count matrix and metadata
            +-- Run DESeq2 (unshrunken LFC)
            +-- Run edgeR (TMM normalization)
            +-- Run limma-voom
            +-- Compare results across methods
            +-- Flag high-confidence peaks (significant in >=2 of 3 methods)
            +-- Generate outputs:
                +-- DA_comparison.csv (side-by-side method comparison)
                +-- DA_high_confidence.csv (peaks passing concordance filter)
        +-- Unload R module

### Optional Post-Pipeline Steps (Manual)

After the main pipeline completes, you can optionally integrate ATAC-seq results with RNA-seq data:

```
OPTIONAL: RNA-ATAC Integration (manual)
    |
    +-- gene_integration_prep.sh
        |
        +-- Extract TSS coordinates from GTF
        +-- Map peaks to nearest TSS (bedtools closest)
        +-- Split peaks into promoter (≤2kb) vs enhancer (2kb-50kb)
        +-- Subset DA results by promoter/enhancer
        +-- Outputs:
            +-- analysis/gene_integration/tss.bed
            +-- analysis/gene_integration/union_peaks_closestTSS.tsv
            +-- analysis/gene_integration/promoter_closest.tsv
            +-- analysis/gene_integration/enhancer_closest.tsv
            +-- analysis/gene_integration/promoter_DA.csv
            +-- analysis/gene_integration/enhancer_DA.csv
    |
    +-- rna_atac_figures.R
        |
        +-- Load RNA DEGs (from separate RNA-seq analysis)
        +-- Load promoter/enhancer DA results
        +-- Collapse peak-level DA to gene-level (best peak per gene)
        +-- Merge RNA + promoter + enhancer at gene level
        +-- Classify genes: RNA+ATAC, RNA-only, ATAC-only, not sig
        +-- Generate correlation plots:
            +-- RNA vs Promoter log2FC
            +-- RNA vs Enhancer log2FC
            +-- Promoter vs Enhancer log2FC
        +-- Outputs:
            +-- analysis/gene_integration/genelevel_RNA_Prom_Enh_merged.csv
            +-- analysis/gene_integration/top30_RNA_genes.csv
            +-- analysis/gene_integration/integrated_genes_ranked.csv
            +-- analysis/gene_integration/top_integrated_genes.csv
            +-- analysis/gene_integration/summary_stats.csv
            +-- analysis/gene_integration/figures/Figure1_threepanel_fullrange_RNA_Prom_Enh.png
            +-- analysis/gene_integration/figures/Figure2_threepanel_DEscale_RNA_Prom_Enh.png
    |
    +-- Driver Inference (optional, requires RNA-seq and motif data)
        |
        +-- Load R/4.2.0 module
        +-- Run driver_inference.R
            |
            +-- Load RNA DEGs
            +-- Load ATAC DA results (promoter + enhancer)
            +-- Load HOMER motif enrichment results
            +-- Identify TFs with DA at own locus
            +-- Identify TFs with motif enrichment in DA peaks
            +-- Identify TFs with concordant RNA expression changes
            +-- Integrate evidence and tier candidates:
                +-- Tier 1: all three evidence lines
                +-- Tier 2: two of three evidence lines
                +-- Tier 3: single evidence line (exploratory)
            +-- Generate outputs:
                +-- driver_candidates_tiered.csv
                +-- driver_evidence_summary.txt
```

## Detailed Step-by-Step Guide

### 1. Initial Setup (One-time)

Verify configuration file paths:

```bash
# Edit if needed
nano analysis/config.sh

# Required variables:
# - ATAC_BASE: Project root
# - ATAC_RAW_FASTQ: Raw FASTQ location
# - ATAC_REFERENCE_DIR: Reference files directory
# - MACS_BIN: MACS3 executable path
```

Verify sample metadata:

```bash
# Check samples.tsv exists and is correct
cat analysis/samples.tsv

# Expected format (tab-separated):
# sample	condition	replicate
# MEF-ATF5WT-n2	WT	n2
# MEF-ATF5WT-n3	WT	n3
# ...
```

### 2. Pre-submission Checklist

```bash
# Ensure log directory exists
mkdir -p /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/analysis/logs

# Verify reference files exist
ls /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/raw/mm10-blacklist.v2.bed
ls /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/raw/bowtie2_index/index.*.bt2
ls /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/raw/mm10.chrom.sizes

# Verify raw FASTQ files accessible
ls /sc/arion/projects/Chipuk_Laboratory/henaoj02/20251107_ATAC_seq_Fibroblasts_WT_vs_ATF5KO/30-1227913735/00_fastq/*_R1_001.fastq.gz
```

### 3. Submit Pipeline

```bash
cd /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO
bsub < analysis/run_atac.lsf
```

Expected output:
```
Job <JOBID> is submitted to queue <premium>.
```

### 4. Monitor Execution

```bash
# Check job status
bjobs -w

# View LSF stdout/stderr
tail -f analysis/logs/<JOBID>.out
tail -f analysis/logs/<JOBID>.err

# View pipeline log (more detailed)
tail -f analysis/logs/<JOBID>.pipeline.log

# Check resource usage
bjobs -l <JOBID>
```

### 5. Verify Completion

Pipeline completion indicators:

```bash
# Check for final DESeq2 output
ls -lh analysis/subread/deseq/DA_results_DESeq2.csv

# Check for MultiQC report
ls -lh analysis/multiqc/multiqc_report.html

# Check pipeline log for "[DONE]" message
tail analysis/logs/*.pipeline.log
```

### 6. Validate Outputs

Run validation script:

```bash
bash analysis/validate_outputs.sh /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/analysis
```

Expected output: All checks should PASS.

### 7. Review Visualizations

**PCA plot** (generated automatically by pipeline):
```bash
# View PCA plot
open analysis/subread/deseq/PCA_samples.pdf

# Or copy to local machine for viewing
scp user@minerva:/path/to/analysis/subread/deseq/PCA_samples.pdf .
```

**PCA plot features:**
- VST-transformed counts (variance-stabilizing transformation)
- PC1 vs PC2 with percent variance explained
- Color-coded by condition (WT vs NULL/KO)
- Sample names labeled on each point
- Should show clear separation between conditions if biological signal is strong

**MultiQC report** (interactive HTML):
```bash
# View MultiQC report (comprehensive QC dashboard)
open analysis/multiqc/multiqc_report.html

# Includes:
# - Read quality scores (FastQC)
# - Adapter content
# - Alignment rates (Bowtie2)
# - Duplication rates (Picard)
# - Library complexity estimates
# - Feature assignment stats (featureCounts)
```

**HOMER motif enrichment** (interactive HTML):
```bash
# View motif analysis
open analysis/homer/motifs_union.size200/homerResults.html

# Shows:
# - Enriched transcription factor motifs
# - Motif logos
# - P-values and enrichment statistics
# - Known TF matches
```

### 8. RNA-ATAC Integration (Optional)

**Prerequisites:**
- RNA-seq DEG results (CSV with columns: gene_id, log2FoldChange, padj)
- GTF annotation file (e.g., gencode.vM25.annotation.gtf)

**Step 7a: Prepare peak-to-gene mapping**

```bash
# Run gene integration prep script
cd /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO

# Option 1: Use default paths (edit script first to set paths)
bash analysis/gene_integration_prep.sh

# Option 2: Provide paths via CLI
bash analysis/gene_integration_prep.sh \
  --base /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO \
  --gtf /path/to/gencode.vM25.annotation.gtf \
  --outdir analysis/gene_integration \
  --prom_win 2000 \
  --enh_max 50000
```

This creates:
- Promoter peak mappings (peaks within ±2kb of TSS)
- Enhancer peak mappings (peaks 2kb-50kb from TSS)
- Promoter/enhancer-specific DA results

**Step 7b: Generate RNA-ATAC integration figures**

```bash
# Load R module
module load R/4.2.0

# Run integration visualization script
Rscript analysis/rna_atac_figures.R \
  --rna_deg_file /path/to/RNA_DEGs.csv \
  --datadir analysis/gene_integration \
  --alpha 0.05 \
  --topN 30

# Outputs will be in:
# - analysis/gene_integration/genelevel_RNA_Prom_Enh_merged.csv
# - analysis/gene_integration/integrated_genes_ranked.csv
# - analysis/gene_integration/top_integrated_genes.csv
# - analysis/gene_integration/figures/*.png
```

**Expected outputs:**
- Gene-level merged table (RNA + promoter + enhancer DA)
- Top N DEGs table
- Ranked integrated genes table (RNA + ATAC significance gate, rank-product ordering)
- Summary statistics (overlap counts, correlations)
- Three-panel scatter plots:
  - RNA log2FC vs Promoter log2FC
  - RNA log2FC vs Enhancer log2FC
  - Promoter log2FC vs Enhancer log2FC
- Both full-range and DE-scale versions

**Integration figure features:**
- Color-coded by significance class (RNA+ATAC, RNA-only, ATAC-only, not sig)
- Spearman correlation coefficients
- Top N genes labeled
- Colorblind-friendly palette

## Troubleshooting and Re-runs

### Check for Errors

```bash
# Search for ERROR messages in pipeline log
grep -i error analysis/logs/*.pipeline.log

# Check LSF stderr
cat analysis/logs/*.err

# Check DESeq2 stderr log
cat analysis/logs/deseq2.stderr.log
```

### Partial Re-runs

Pipeline uses checkpoint logic: already-generated files are skipped (marked with `[SKIP]`). To force regeneration of specific steps:

**Re-run counting and differential analysis only:**

```bash
# Preview what will be deleted
bash analysis/cleanup_for_rerun.sh --dry-run

# Delete counting/DESeq2 outputs
bash analysis/cleanup_for_rerun.sh

# Re-submit
bsub < analysis/run_atac.lsf
```

**Re-run entire pipeline from scratch:**

```bash
# WARNING: This deletes ALL analysis outputs
rm -rf analysis/fastqc analysis/trimmgalore analysis/bowtie2 \
       analysis/samtools analysis/picard analysis/deeptools \
       analysis/macs analysis/subread analysis/homer analysis/multiqc

# Re-submit
mkdir -p analysis/logs
bsub < analysis/run_atac.lsf
```

### Common Issues

**Issue: MACS3 not found**
```bash
# Verify path
ls -l /hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3

# If different, edit config.sh:
export MACS_BIN="/path/to/macs3"
```

**Issue: R packages missing (main pipeline)**
```bash
# Load R and install packages
module load R/4.2.0
R
> if (!require("BiocManager")) install.packages("BiocManager")
> BiocManager::install(c("DESeq2", "apeglm"))
> install.packages(c("readr", "dplyr", "ggplot2", "stringr", "tibble"))
```

**Issue: R packages missing (RNA-ATAC integration)**
```bash
module load R/4.2.0
R
> install.packages(c("readr", "dplyr", "stringr", "ggplot2", "ggrepel", "optparse"))
> # Optional but recommended for better figure layouts:
> install.packages("patchwork")
```

**Issue: Conda environment not found**
```bash
# List available environments
conda env list

# If csb doesn't exist, create it
conda create -n csb -c bioconda -c conda-forge \
  fastqc trim-galore bowtie2 samtools bedtools \
  deeptools subread multiqc picard
```

## Output Files

### Key Outputs

| Output | Path | Description |
|--------|------|-------------|
| Final BAMs | `analysis/samtools/*.final.bam` | Filtered, deduplicated BAMs (MAPQ>=20, proper pairs, no chrM) |
| BigWigs | `analysis/deeptools/*.CPM.bw` | CPM-normalized coverage tracks (for IGV) |
| IDR summary | `analysis/macs/idr/idr_summary.tsv` | IDR analysis summary (Nt, Np, optimal peak counts) |
| IDR consensus peaks | `analysis/macs/idr/idr_consensus_peaks.bed` | Optimal peak set from IDR analysis (WT + KO merged) |
| IDR method flag | `analysis/macs/idr/idr_method.txt` | Indicates whether IDR or union fallback was used |
| Union peaks | `analysis/macs/union_peaks.bed` | Merged peak set (blacklist-filtered, used if no IDR) |
| Count matrix | `analysis/subread/union_peaks_featureCounts.txt` | Read counts per peak per sample |
| DA results | `analysis/subread/deseq/DA_results_DESeq2.csv` | All peaks with log2FC, p-values, padj (apeglm shrinkage) |
| DA up (KO) | `analysis/subread/deseq/DA_KO_up.csv` | Sig. peaks: padj<=0.05, log2FC>=1 |
| DA down (KO) | `analysis/subread/deseq/DA_KO_down.csv` | Sig. peaks: padj<=0.05, log2FC<=-1 |
| DA comparison | `analysis/subread/concordance/DA_comparison.csv` | Side-by-side comparison of DESeq2, edgeR, limma-voom |
| DA high-confidence | `analysis/subread/concordance/DA_high_confidence.csv` | Peaks significant in >=2 of 3 methods |
| PCA plot | `analysis/subread/deseq/PCA_samples.pdf` | Sample PCA (VST-transformed, labeled) |
| PCA coordinates | `analysis/subread/deseq/PCA_coordinates.csv` | PC1/PC2 values + variance explained |
| HOMER annotation | `analysis/homer/union_peaks.annotatePeaks.size200.txt` | Peak genomic context |
| HOMER motifs | `analysis/homer/motifs_union.size200/homerResults.html` | Enriched motifs |
| QC table (basic) | `analysis/subread/qc_metrics.tsv` | Alignment and duplication stats |
| QC table (extended) | `analysis/subread/qc_metrics_extended.tsv` | FRiP, TSS enrichment, fragment size, NFR, blacklist fraction, PBC1/PBC2 |
| MultiQC | `analysis/multiqc/multiqc_report.html` | Aggregated QC report |
| Versions | `analysis/logs/versions.txt` | Software versions used |

### Visualizations Generated by Main Pipeline

The main ATAC-seq pipeline automatically generates the following visualizations:

| Visualization | Path | Description |
|---------------|------|-------------|
| PCA plot | `analysis/subread/deseq/PCA_samples.pdf` | Sample clustering by PC1 vs PC2 (VST-transformed counts), colored by condition, labeled with sample names. Shows percent variance explained by each PC. |
| MultiQC report | `analysis/multiqc/multiqc_report.html` | Interactive HTML report aggregating QC metrics from all pipeline steps (FastQC, Bowtie2, Picard, featureCounts, etc.). Includes read quality, alignment rates, duplication rates, GC content, and more. |
| HOMER motifs | `analysis/homer/motifs_union.size200/homerResults.html` | Interactive HTML report of enriched motifs in accessible regions, including motif logos, p-values, and target percentages. |

**BigWig files for genome browser visualization:**
- `analysis/deeptools/*.CPM.bw` - CPM-normalized coverage tracks for each sample (load into IGV/UCSC Genome Browser)

### Optional RNA-ATAC Integration Outputs

| Output | Path | Description |
|--------|------|-------------|
| TSS coordinates | `analysis/gene_integration/tss.bed` | Gene TSS positions (1bp) |
| Peak-to-TSS map | `analysis/gene_integration/union_peaks_closestTSS.tsv` | All peaks mapped to nearest TSS |
| Promoter map | `analysis/gene_integration/promoter_closest.tsv` | Peaks within ±2kb of TSS |
| Enhancer map | `analysis/gene_integration/enhancer_closest.tsv` | Peaks 2kb-50kb from TSS |
| Promoter DA | `analysis/gene_integration/promoter_DA.csv` | DA results for promoter peaks |
| Enhancer DA | `analysis/gene_integration/enhancer_DA.csv` | DA results for enhancer peaks |
| Gene-level merged | `analysis/gene_integration/genelevel_RNA_Prom_Enh_merged.csv` | RNA + ATAC merged at gene level |
| Top DEGs | `analysis/gene_integration/top30_RNA_genes.csv` | Top N differentially expressed genes |
| Integrated ranking (all) | `analysis/gene_integration/integrated_genes_ranked.csv` | Genes significant in RNA and >=1 ATAC layer, ranked by multi-modal rank product |
| Top integrated genes | `analysis/gene_integration/top_integrated_genes.csv` | Top N integrated genes from combined RNA+ATAC ranking |
| Summary stats | `analysis/gene_integration/summary_stats.csv` | Overlap counts and correlations |
| Integration figures | `analysis/gene_integration/figures/*.png` | RNA-ATAC correlation plots |
| Driver candidates (tiered) | `analysis/gene_integration/driver_candidates_tiered.csv` | TF candidates with evidence tier (1/2/3) and supporting evidence |
| Driver evidence summary | `analysis/gene_integration/driver_evidence_summary.txt` | Human-readable summary of driver inference results |

### Directory Structure

```
analysis/
├── fastqc/              # Raw read QC reports
├── trimmgalore/         # Trimmed FASTQs + QC
├── bowtie2/             # Alignment logs
├── samtools/            # BAM files (aligned, filtered, final)
├── picard/              # Deduplication metrics
├── deeptools/           # BigWig coverage tracks
├── macs/                # Peak files (per-replicate and union)
│   ├── rep_peaks/       # Per-replicate peaks
│   ├── rep_peaks_bl/    # Blacklist-filtered peaks
│   ├── idr/             # IDR analysis outputs
│   │   ├── idr_summary.tsv
│   │   ├── idr_consensus_peaks.bed
│   │   └── idr_method.txt
│   ├── union_peaks.bed
│   └── union_peaks.saf
├── homer/               # Annotation and motifs
│   ├── union_peaks.annotatePeaks.size200.txt
│   ├── motifs_union.size200/
│   └── tagdirs/         # Tag directories per sample
├── subread/             # Quantification
│   ├── union_peaks_featureCounts.txt
│   ├── qc_metrics.tsv
│   ├── qc_metrics_extended.tsv
│   ├── deseq/           # Differential analysis
│   │   ├── DA_results_DESeq2.csv
│   │   ├── DA_KO_up.csv
│   │   ├── DA_KO_down.csv
│   │   ├── PCA_samples.pdf
│   │   └── PCA_coordinates.csv
│   └── concordance/     # Multi-method comparison
│       ├── DA_comparison.csv
│       └── DA_high_confidence.csv
├── gene_integration/    # Optional: RNA-ATAC integration
│   ├── tss.bed
│   ├── peaks.center.bed
│   ├── union_peaks_closestTSS.tsv
│   ├── promoter_closest.tsv
│   ├── enhancer_closest.tsv
│   ├── promoter_peak_ids.txt
│   ├── enhancer_peak_ids.txt
│   ├── promoter_DA.csv
│   ├── enhancer_DA.csv
│   ├── genelevel_RNA_Prom_Enh_merged.csv
│   ├── top30_RNA_genes.csv
│   ├── integrated_genes_ranked.csv
│   ├── top_integrated_genes.csv
│   ├── summary_stats.csv
│   ├── driver_candidates_tiered.csv
│   ├── driver_evidence_summary.txt
│   └── figures/
│       ├── Figure1_threepanel_fullrange_RNA_Prom_Enh.png
│       └── Figure2_threepanel_DEscale_RNA_Prom_Enh.png
├── multiqc/             # Aggregated QC report
└── logs/                # Pipeline and tool logs
    ├── <JOBID>.out
    ├── <JOBID>.err
    ├── <JOBID>.pipeline.log
    ├── versions.txt
    └── deseq2.stderr.log
```

## Validation

### Quick Validation

```bash
# Check BAM integrity
for bam in analysis/samtools/*.final.bam; do
  samtools quickcheck "$bam" && echo "OK: $(basename "$bam")"
done

# Count union peaks
wc -l analysis/macs/union_peaks.bed

# Count DA peaks
wc -l analysis/subread/deseq/DA_KO_*.csv

# Preview QC metrics
column -t analysis/subread/qc_metrics.tsv | head
```

### Comprehensive Validation

```bash
bash analysis/validate_outputs.sh /path/to/analysis
```

Checks:
- Final BAM files exist and are indexed
- BigWig files exist
- Peak files exist (IDR consensus or union peaks)
- featureCounts output exists
- DESeq2 outputs exist (DA results, PCA, up/down peaks)
- Concordance outputs exist (DA comparison, high-confidence peaks)
- HOMER outputs exist (annotation, motifs)
- MultiQC report exists
- QC tables exist (basic and extended metrics)
- IDR outputs exist (if IDR was used)

## Analysis Notes

### Biological Context

- **Study**: ATF5 knockout in mouse embryonic fibroblasts (MEFs)
- **Genome**: mm10 (GRCm38)
- **Library prep**: Nextera-based ATAC-seq (Tn5 transposase)
- **Replicates**: 3 biological replicates per condition (WT and KO)
- **Batch structure**: Single batch (no batch correction needed)

### Methods

#### IDR (Irreproducible Discovery Rate)

Peak reproducibility was assessed using the ENCODE IDR framework. For each condition, replicate BAMs were pooled and peaks called on the pooled sample. Two pseudoreplicates were generated by random 50/50 read splitting. The optimal peak set was determined using the ENCODE max(Nt, Np) rule: Nt = peaks passing IDR threshold (0.05) in true replicate comparison; Np = peaks passing threshold in pooled pseudoreplicate comparison. The top max(Nt, Np) peaks from the pooled peak set constitute the optimal set per condition. Consensus peaks were formed by merging WT-optimal and KO-optimal peak sets. If IDR software is unavailable, the pipeline falls back to union peak calling.

#### QC Metrics

Extended quality control metrics are computed per sample: FRiP (Fraction of Reads in Peaks), TSS enrichment score, mean/median fragment size, nucleosome-free ratio (NFR; fragments less than 150bp divided by fragments 150-300bp), blacklist fraction, and library complexity (PBC1, PBC2). ENCODE minimum thresholds: FRiP greater than or equal to 0.1, TSS enrichment greater than or equal to 5, NFR ratio greater than or equal to 1.0.

#### Multi-Method Concordance

Differential accessibility is independently assessed by DESeq2, edgeR, and limma-voom using unshrunken log2 fold changes for harmonized comparison. Peaks significant (padj less than 0.05, absolute LFC greater than or equal to 1) in at least 2 of 3 methods are flagged as high-confidence. This concordance filter serves as a robustness/stability check, not a formal statistical proof of differential accessibility.

#### Peak-to-Gene Assignment

Peaks are assigned to genes using multi-gene window assignment within a configurable window (default: 50kb from TSS). Assignments are distance-weighted using an exponential decay function (weight = 1/(1 + dist/decay), default decay = 10kb). Promoter peaks (within 2kb of TSS) are prioritized. Distal (enhancer) assignments are capped at 5 genes per peak by weight rank.

#### Integration Statistics

Peak-to-gene integration uses signed Stouffer Z-score combination on raw p-values with post-hoc BH correction. Cross-modality evidence is ranked by harmonic mean p-value (HMP) of raw per-modality p-values. HMP serves as an evidence ranking metric, not a formal joint FDR.

#### Driver Inference

Candidate driver transcription factors are identified by integrating three evidence lines: (1) differential chromatin accessibility at the TF's own locus, (2) enrichment of the TF's binding motif in directionally-changed peaks (HOMER findMotifsGenome.pl), and (3) concordant RNA expression changes. Candidates are tiered: Tier 1 (all three evidence lines), Tier 2 (two of three), Tier 3 (single evidence, exploratory). This is hypothesis-generating evidence integration, not causal proof of driver status.

#### Interpretation Boundaries

Integrated rankings identify candidates for experimental validation, not confirmed regulatory relationships. Statistical associations between chromatin accessibility and gene expression do not establish causality.

### Pipeline Parameters

**Main ATAC-seq Pipeline:**
- **Alignment**: Bowtie2 `--very-sensitive`, max fragment size 2000bp
- **Filtering**: MAPQ>=20, proper pairs only, canonical chromosomes only (chr1-19,X,Y), no chrM
- **Deduplication**: Picard MarkDuplicates with `REMOVE_DUPLICATES=true`
- **Normalization**: CPM (counts per million) for BigWigs
- **Peak calling**: MACS3 BAMPE mode, no model, genome size 1.87e9
- **Blacklist**: ENCODE mm10 blacklist v2
- **Motif search**: HOMER, window size 200bp, motif lengths 8,10,12

**Optional RNA-ATAC Integration:**
- **Promoter definition**: Peaks within ±2kb of TSS
- **Enhancer definition**: Peaks 2kb-50kb from TSS
- **Gene-level collapse**: Best peak per gene (min padj, tie-break by max abs(log2FC))
- **Integrated gene ranking**: RNA-significant + >=1 ATAC-significant genes ranked by geometric mean of modality ranks (RNA/promoter/enhancer); harmonic mean p-value retained as secondary evidence
- **Significance threshold**: padj < 0.05 (default, adjustable)
- **Top genes for labeling**: Up to top 30, prioritizing significant genes across RNA/promoter/enhancer by best padj and effect size (default, adjustable)

### Reproducibility

All run metadata recorded in `analysis/logs/versions.txt`:
- Software versions
- Configuration file path
- Conda environment
- Tool paths

## References

- ENCODE ATAC-seq pipeline: https://www.encodeproject.org/atac-seq/
- mm10 blacklist v2: https://github.com/Boyle-Lab/Blacklist
- MACS3: https://github.com/macs3-project/MACS
- HOMER: http://homer.ucsd.edu/homer/

## Contact

Andres Chousal Cantu
Icahn School of Medicine at Mount Sinai
