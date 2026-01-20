#!/usr/bin/env bash
# Validate ATAC-seq pipeline outputs
# Usage: validate_outputs.sh <analysis_dir>

set -euo pipefail

ANALYSIS_DIR="${1:?Usage: validate_outputs.sh <analysis_dir>}"

echo "=== ATAC-seq Output Validation ==="
echo "Analysis directory: $ANALYSIS_DIR"
echo ""

errors=0
warnings=0

check_pass() { echo "[PASS] $1"; }
check_fail() { echo "[FAIL] $1"; ((errors++)); }
check_warn() { echo "[WARN] $1"; ((warnings++)); }

# -----------------------------------------------------------------------------
# 1. Directory structure
# -----------------------------------------------------------------------------
echo "--- Checking directory structure ---"

for dir in samtools picard deeptools macs subread homer multiqc logs; do
    if [[ -d "$ANALYSIS_DIR/$dir" ]]; then
        check_pass "Directory exists: $dir/"
    else
        check_fail "Directory missing: $dir/"
    fi
done

# -----------------------------------------------------------------------------
# 2. BAM integrity
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking BAM integrity ---"

bam_count=0
for bam in "$ANALYSIS_DIR"/samtools/*.final.bam; do
    [[ -f "$bam" ]] || continue
    ((bam_count++))
    if samtools quickcheck "$bam" 2>/dev/null; then
        check_pass "BAM valid: $(basename "$bam")"
    else
        check_fail "BAM corrupt: $(basename "$bam")"
    fi

    # Check index exists
    if [[ -f "${bam}.bai" ]]; then
        check_pass "Index exists: $(basename "$bam").bai"
    else
        check_fail "Index missing: $(basename "$bam").bai"
    fi
done

if [[ $bam_count -eq 0 ]]; then
    check_fail "No .final.bam files found"
fi

# -----------------------------------------------------------------------------
# 3. BigWig files
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking BigWig files ---"

bw_count=0
for bw in "$ANALYSIS_DIR"/deeptools/*.CPM.bw; do
    [[ -f "$bw" ]] || continue
    ((bw_count++))
    size=$(stat -f%z "$bw" 2>/dev/null || stat -c%s "$bw" 2>/dev/null || echo 0)
    if [[ $size -gt 1000 ]]; then
        check_pass "BigWig valid: $(basename "$bw") (${size} bytes)"
    else
        check_fail "BigWig too small: $(basename "$bw") (${size} bytes)"
    fi
done

if [[ $bw_count -eq 0 ]]; then
    check_fail "No .CPM.bw files found"
fi

# -----------------------------------------------------------------------------
# 4. Peak files
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking peak files ---"

union_bed="$ANALYSIS_DIR/macs/union_peaks.bed"
union_saf="$ANALYSIS_DIR/macs/union_peaks.saf"

if [[ -s "$union_bed" ]]; then
    peak_count=$(wc -l < "$union_bed")
    check_pass "Union peaks: $peak_count peaks"
else
    check_fail "Union peaks missing or empty"
fi

if [[ -s "$union_saf" ]]; then
    # Verify 1-based coordinates (column 3 should be > 0)
    zero_starts=$(awk -F'\t' '$3 <= 0 {count++} END {print count+0}' "$union_saf")
    if [[ $zero_starts -eq 0 ]]; then
        check_pass "SAF coordinates are 1-based"
    else
        check_fail "SAF has $zero_starts rows with 0-based or negative starts"
    fi
else
    check_fail "SAF file missing or empty"
fi

# Check per-replicate peaks
rep_peak_count=$(find "$ANALYSIS_DIR/macs/rep_peaks_bl" -name "*.bl.narrowPeak" 2>/dev/null | wc -l)
if [[ $rep_peak_count -gt 0 ]]; then
    check_pass "Per-replicate BL-filtered peaks: $rep_peak_count files"
else
    check_warn "No per-replicate peaks found"
fi

# -----------------------------------------------------------------------------
# 5. Count matrix
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking count matrix ---"

fc_out="$ANALYSIS_DIR/subread/union_peaks_featureCounts.txt"
if [[ -s "$fc_out" ]]; then
    samples_in_fc=$(head -2 "$fc_out" | tail -1 | awk -F'\t' '{print NF-6}')
    check_pass "featureCounts matrix: $samples_in_fc samples"
else
    check_fail "featureCounts output missing"
fi

# -----------------------------------------------------------------------------
# 6. DESeq2 results
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking DESeq2 results ---"

deseq_dir="$ANALYSIS_DIR/subread/deseq"
da_results="$deseq_dir/DA_results_DESeq2.csv"

if [[ -s "$da_results" ]]; then
    da_rows=$(wc -l < "$da_results")
    check_pass "DESeq2 results: $((da_rows - 1)) peaks tested"
else
    check_fail "DESeq2 results missing"
fi

for subset in DA_KO_up.csv DA_KO_down.csv; do
    if [[ -s "$deseq_dir/$subset" ]]; then
        count=$(($(wc -l < "$deseq_dir/$subset") - 1))
        check_pass "$subset: $count differential peaks"
    else
        check_warn "$subset missing or empty"
    fi
done

if [[ -s "$deseq_dir/PCA_samples.pdf" ]]; then
    check_pass "PCA plot exists"
else
    check_warn "PCA plot missing"
fi

# -----------------------------------------------------------------------------
# 7. HOMER outputs
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking HOMER outputs ---"

homer_anno="$ANALYSIS_DIR/homer/union_peaks.annotatePeaks.size200.txt"
homer_motifs="$ANALYSIS_DIR/homer/motifs_union.size200/homerResults.html"

if [[ -s "$homer_anno" ]]; then
    check_pass "HOMER annotation exists"
else
    check_warn "HOMER annotation missing"
fi

if [[ -s "$homer_motifs" ]]; then
    check_pass "HOMER motif results exist"
else
    check_warn "HOMER motif results missing"
fi

# -----------------------------------------------------------------------------
# 8. QC files
# -----------------------------------------------------------------------------
echo ""
echo "--- Checking QC files ---"

qc_tsv="$ANALYSIS_DIR/subread/qc_metrics.tsv"
if [[ -s "$qc_tsv" ]]; then
    qc_samples=$(($(wc -l < "$qc_tsv") - 1))
    check_pass "QC metrics: $qc_samples samples"
else
    check_warn "QC metrics missing"
fi

multiqc="$ANALYSIS_DIR/multiqc/multiqc_report.html"
if [[ -s "$multiqc" ]]; then
    check_pass "MultiQC report exists"
else
    check_warn "MultiQC report missing"
fi

versions="$ANALYSIS_DIR/logs/versions.txt"
if [[ -s "$versions" ]]; then
    check_pass "Versions log exists"
else
    check_warn "Versions log missing"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=== Summary ==="
echo "Errors:   $errors"
echo "Warnings: $warnings"

if [[ $errors -gt 0 ]]; then
    echo ""
    echo "VALIDATION FAILED - $errors critical issues found"
    exit 1
else
    echo ""
    echo "VALIDATION PASSED"
    exit 0
fi
