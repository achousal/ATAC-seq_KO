#!/usr/bin/env bash
# =============================================================================
# config.sh - ATAC-seq Pipeline Configuration
# =============================================================================
# Edit paths and parameters for your environment.
# This file is sourced by run.sh and atacseq_analysis.sh.
# =============================================================================

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

# Project root (where analysis/, raw/, results/ live)
export ATAC_BASE="/sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO"

# Raw FASTQ location
export ATAC_RAW_FASTQ="/sc/arion/projects/Chipuk_Laboratory/henaoj02/20251107_ATAC_seq_Fibroblasts_WT_vs_ATF5KO/30-1227913735/00_fastq"

# Reference files
export ATAC_REFERENCE_DIR="${ATAC_BASE}/raw"
export ATAC_BLACKLIST="${ATAC_REFERENCE_DIR}/mm10-blacklist.v2.bed"
export ATAC_BOWTIE_INDEX="${ATAC_REFERENCE_DIR}/bowtie2_index/index"
export ATAC_CHROM_SIZES="${ATAC_REFERENCE_DIR}/mm10.chrom.sizes"

# Sample metadata (TSV: sample, condition, replicate)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export ATAC_SAMPLE_META="${SCRIPT_DIR}/samples.tsv"

# -----------------------------------------------------------------------------
# TOOLS
# -----------------------------------------------------------------------------

# MACS3 path (Minerva Rocky9)
export MACS_BIN="/hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3"

# Conda environment name
export ATAC_CONDA_ENV="csb"

# -----------------------------------------------------------------------------
# PARAMETERS
# -----------------------------------------------------------------------------

# Threads (auto-detect from LSF or default to 8)
export ATAC_THREADS="${LSB_MAX_NUM_PROCESSORS:-8}"

# Alignment quality filter
export ATAC_MAPQ=20

# MACS3 genome size (mm10 effective)
export ATAC_GENOME_SIZE="1.87e9"

# BigWig bin size
export ATAC_BW_BINSIZE=10

# HOMER settings
export ATAC_HOMER_SIZE=200
export ATAC_HOMER_MOTIF_LEN="8,10,12"

# Keep intermediate files (0 = delete, 1 = keep)
export ATAC_KEEP_INTERMEDIATES=0

# -----------------------------------------------------------------------------
# CONDITION PATTERNS (fallback if samples.tsv missing)
# -----------------------------------------------------------------------------
export ATAC_WT_PATTERN="ATF5WT"
export ATAC_KO_PATTERN="ATF5NULL"
