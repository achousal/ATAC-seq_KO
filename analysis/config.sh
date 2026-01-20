#!/usr/bin/env bash
# ATAC-seq pipeline configuration
# Source this file from atacseq_analysis.sh

# =============================================================================
# PATHS - Edit for your environment
# =============================================================================

# Project base directory (where analysis/ and raw/ live)
export ATAC_BASE="/sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO"

# Raw FASTQ location
export ATAC_RAW_FASTQ="/sc/arion/projects/Chipuk_Laboratory/henaoj02/20251107_ATAC_seq_Fibroblasts_WT_vs_ATF5KO/30-1227913735/00_fastq"

# Reference files directory
export ATAC_REFERENCE_DIR="${ATAC_BASE}/raw"

# Reference files (mm10)
export ATAC_BLACKLIST="${ATAC_REFERENCE_DIR}/mm10-blacklist.v2.bed"
export ATAC_BOWTIE_INDEX="${ATAC_REFERENCE_DIR}/bowtie2_index/index"
export ATAC_CHROM_SIZES="${ATAC_REFERENCE_DIR}/mm10.chrom.sizes"

# =============================================================================
# TOOLS
# =============================================================================

# MACS3 binary (Minerva Rocky9 module path)
export MACS_BIN="/hpc/packages/minerva-rocky9/macs/3.0.2/MACS3/bin/macs3"

# =============================================================================
# PARAMETERS
# =============================================================================

# Threads: use LSF allocation if available, otherwise default
export ATAC_THREADS="${LSB_MAX_NUM_PROCESSORS:-8}"

# Alignment filtering
export ATAC_MAPQ=20

# Keep intermediate files (0=delete, 1=keep)
export ATAC_KEEP_INTERMEDIATES=0

# MACS3 genome size (mm10 effective)
export ATAC_GENOME_SIZE="1.87e9"

# BigWig bin size
export ATAC_BW_BINSIZE=10

# HOMER motif search parameters
export ATAC_HOMER_SIZE=200
export ATAC_HOMER_MOTIF_LEN="8,10,12"

# =============================================================================
# SAMPLE METADATA
# =============================================================================

# Sample metadata file (TSV with columns: sample, condition, replicate)
# If not provided or missing, pipeline falls back to pattern-based detection
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export ATAC_SAMPLE_META="${SCRIPT_DIR}/samples.tsv"

# Fallback patterns for condition detection (if no metadata file)
export ATAC_WT_PATTERN="ATF5WT"
export ATAC_KO_PATTERN="ATF5NULL"

# =============================================================================
# CONDA ENVIRONMENT
# =============================================================================

export ATAC_CONDA_ENV="csb"
