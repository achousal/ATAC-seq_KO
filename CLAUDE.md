# CLAUDE.md - ATAC-seq Pipeline (ATF5 KO Study)

## Project Overview

ATAC-seq analysis pipeline for MEF (mouse embryonic fibroblast) samples comparing ATF5 wild-type vs ATF5 knockout conditions. Pipeline runs on Minerva HPC (LSF scheduler).

## Quick Reference

```bash
# Submit pipeline
bsub < code/run_atac.lsf

# Check job status
bjobs -l $JOBID

# View logs
tail -f /sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq/ATF5_KO/analysis/logs/*.pipeline.log
```

## Directory Structure

```
ATAC-seq/
  code/
    config.sh              # Pipeline configuration (paths, params)
    samples.tsv            # Sample metadata
    run_atac.lsf           # LSF submission wrapper
    atacseq_analysis.sh    # Main pipeline script
    deseq_atac.R           # Differential accessibility analysis
    validate_outputs.sh    # Output validation script
    cleanup_for_rerun.sh   # Cleanup script for re-runs
    rna_atac_figures.R     # Integration figures
    gene_integration_prep.sh
  raw/                     # Reference files (on Minerva)
    mm10-blacklist.v2.bed
    bowtie2_index/
    mm10.chrom.sizes
  analysis/                # All outputs (on Minerva)
    fastqc/
    trimmgalore/
    bowtie2/
    samtools/              # *.final.bam (main BAMs)
    picard/
    deeptools/             # *.CPM.bw (BigWigs)
    macs/                  # Peaks
    homer/                 # Annotation + motifs
    subread/               # featureCounts + DESeq2
    multiqc/
    logs/
```

## Samples

| Sample | Condition | Replicate |
|--------|-----------|-----------|
| MEF-ATF5WT-n2 | WT | n2 |
| MEF-ATF5WT-n3 | WT | n3 |
| MEF-ATF5WT-n4 | WT | n4 |
| MEF-ATF5NULL-n2 | NULL (KO) | n2 |
| MEF-ATF5NULL-n3 | NULL (KO) | n3 |
| MEF-ATF5NULL-n4 | NULL (KO) | n4 |

**Library prep**: Nextera-based (standard ATAC-seq)
**Batch structure**: Single batch (no batch correction needed)
**Genome**: mm10 (GRCm38)

## Pipeline Steps

1. **QC + Trim**: FastQC -> TrimGalore (Nextera adapters)
2. **Align**: Bowtie2 --very-sensitive -X 2000
3. **Filter**: MAPQ>=20, proper pairs, canonical chr only (no chrM)
4. **Dedup**: Picard MarkDuplicates (REMOVE_DUPLICATES=true)
5. **BigWigs**: deepTools bamCoverage --normalizeUsing CPM
6. **Peaks**: MACS3 per-replicate -> blacklist filter -> union
7. **Counts**: featureCounts on union peaks
8. **DA**: DESeq2 (NULL vs WT comparison)
9. **Annotation**: HOMER annotatePeaks + findMotifsGenome

## Key Outputs

| Output | Path |
|--------|------|
| Final BAMs | `analysis/samtools/*.final.bam` |
| BigWigs | `analysis/deeptools/*.CPM.bw` |
| Union peaks | `analysis/macs/union_peaks.bed` |
| DA results | `analysis/subread/deseq/DA_results_DESeq2.csv` |
| DA up (KO) | `analysis/subread/deseq/DA_KO_up.csv` |
| DA down (KO) | `analysis/subread/deseq/DA_KO_down.csv` |
| PCA coordinates | `analysis/subread/deseq/PCA_coordinates.csv` |
| Motifs | `analysis/homer/motifs_union.size200/` |
| QC table | `analysis/subread/qc_metrics.tsv` |
| MultiQC | `analysis/multiqc/multiqc_report.html` |

## Environment

```bash
# Minerva
module load anaconda
conda activate csb

# Required tools (in csb env)
fastqc, trim_galore, bowtie2, samtools, bedtools,
deepTools, featureCounts, multiqc, picard

# MACS3 (hardcoded path - Minerva Rocky9)
/hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3

# R (loaded separately)
module load R/4.2.0
# Packages: DESeq2, apeglm, ggplot2, dplyr, readr, stringr, tibble
```

## Configuration

Edit `code/config.sh` to customize:
- `ATAC_BASE`: Project root directory
- `ATAC_RAW_FASTQ`: Location of raw FASTQ files
- `ATAC_REFERENCE_DIR`: Reference files (blacklist, index, chrom sizes)
- `MACS_BIN`: Path to MACS3 executable
- `ATAC_THREADS`: Number of threads (auto-detected from LSF if available)

## Re-running After Changes

If you need to regenerate counting/DESeq2 results:

```bash
# Preview what will be deleted
bash code/cleanup_for_rerun.sh --dry-run

# Delete affected files
bash code/cleanup_for_rerun.sh

# Re-submit pipeline
mkdir -p /sc/arion/.../analysis/logs
bsub < code/run_atac.lsf
```

## Validation

```bash
bash code/validate_outputs.sh /path/to/analysis
```

## Testing

No automated tests currently. To validate outputs:

```bash
# Check BAM integrity
for bam in analysis/samtools/*.final.bam; do
  samtools quickcheck "$bam" && echo "OK: $bam"
done

# Check peak counts
wc -l analysis/macs/union_peaks.bed

# Check DESeq2 ran
head analysis/subread/deseq/DA_results_DESeq2.csv
```

## References

- ENCODE ATAC-seq pipeline: https://www.encodeproject.org/atac-seq/
- mm10 blacklist v2: https://github.com/Boyle-Lab/Blacklist
- MACS3: https://github.com/macs3-project/MACS
