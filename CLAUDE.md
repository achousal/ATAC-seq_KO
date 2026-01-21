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
    gene_integration_prep.sh  # Stage 2a: Peak-to-gene mapping (promoter/enhancer)
    rna_atac_figures.R        # Stage 2b: RNA+ATAC integration and correlation plots
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
| 2. RNA+ATAC Integration | [gene_integration_prep.sh](analysis/gene_integration_prep.sh) + [rna_atac_figures.R](analysis/rna_atac_figures.R) | Peak-gene maps, merged tables, correlation plots | Yes (`--skip-integrate`) |
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

### Figures (in results/)

| Figure | Location |
|--------|----------|
| Volcano plot | `results/figures/da/volcano_ATAC_DA.png` |
| MA plot | `results/figures/da/MA_plot_ATAC_DA.png` |
| PCA plot | `results/figures/qc/PCA_samples.pdf` |
| Heatmap | `results/figures/heatmaps/heatmap_all_peaks.pdf` |

### Integration (optional)

| Output | Location |
|--------|----------|
| Gene-level merge | `analysis/gene_integration/genelevel_RNA_Prom_Enh_merged.csv` |
| Correlation plots | `results/figures/integration/Figure{1,2}_threepanel_*.png` |
| Top DEGs | `analysis/gene_integration/top30_RNA_genes.csv` |

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
--skip-atac         Skip ATAC processing (use existing BAMs/peaks)
--skip-integrate    Skip RNA+ATAC integration
--skip-figures      Skip figure generation
--force             Re-run all steps (ignore checkpoints)
--rna-deg FILE      RNA-seq DEG file for integration
--gtf FILE          GTF annotation for integration
--validate          Run validation only (no execution)
--dry-run           Preview what would run without executing
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
> BiocManager::install(c("DESeq2", "apeglm"))
> install.packages(c("readr", "dplyr", "ggplot2", "stringr", "tibble", "ggrepel", "optparse"))
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
