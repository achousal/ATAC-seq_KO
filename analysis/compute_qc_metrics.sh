#!/usr/bin/env bash
# =============================================================================
# compute_qc_metrics.sh - Extended QC metrics for ATAC-seq
# =============================================================================
# Computes per-sample ATAC-seq QC metrics:
#   - FRiP (Fraction of Reads in Peaks)
#   - TSS enrichment
#   - Fragment size distribution (mean/median)
#   - NFR ratio (nucleosome-free vs mononucleosome)
#   - Blacklist fraction
#   - PBC1/PBC2 (library complexity)
#
# Outputs: qc_metrics_extended.tsv with pass/fail flags per metric
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

# =============================================================================
# HELPERS
# =============================================================================
log(){ printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*"; }
skip(){ log "[SKIP] $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"; }

usage(){
  cat <<EOF
Usage: $(basename "$0") [options]

Compute extended QC metrics for ATAC-seq samples.

Options:
  --samtools-dir DIR    Directory containing final BAMs (*.final.bam)
  --peaks-bed FILE      Union peaks BED file
  --blacklist FILE      Blacklist regions BED
  --tss-bed FILE        TSS BED file (chr start end gene_id gene_name strand)
  --bigwig-dir DIR      Directory containing BigWig files (optional, for deepTools)
  --gtf FILE            GTF annotation (auto-generates TSS BED if --tss-bed empty)
  --threads INT         Thread count (default: from config or 8)
  --outdir DIR          Output directory (default: analysis/qc)
  --policy POLICY       QC failure policy: warn|exclude (default: warn)
  -h, --help            Show this help

Environment:
  Sources config.sh from same directory if available.
  Uses ATAC_QC_FRIP_MIN, ATAC_QC_TSS_ENRICH_MIN, ATAC_QC_NFR_RATIO_MIN,
  ATAC_QC_FAIL_POLICY from config.

Output:
  qc_metrics_extended.tsv - One row per sample with metrics and pass/fail flags
  qc_excluded_samples.txt - List of samples failing QC (if --policy exclude)
EOF
}

# =============================================================================
# LOAD CONFIGURATION
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.sh"

if [[ -f "$CONFIG_FILE" ]]; then
    log "Loading config: $CONFIG_FILE"
    source "$CONFIG_FILE"
fi

# =============================================================================
# DEFAULT PARAMETERS (can be overridden by config or CLI)
# =============================================================================
SAMTOOLS_DIR="${SCRIPT_DIR}/samtools"
PEAKS_BED="${SCRIPT_DIR}/macs/union_peaks.bed"
BLACKLIST="${ATAC_BLACKLIST:-}"
TSS_BED="${ATAC_TSS_BED:-}"
BIGWIG_DIR="${SCRIPT_DIR}/deeptools"
GTF="${ATAC_REFERENCE_DIR:-${SCRIPT_DIR}/../raw}/gencode.vM25.annotation.gtf"
THREADS="${ATAC_THREADS:-8}"
OUTDIR="${SCRIPT_DIR}/qc"
POLICY="${ATAC_QC_FAIL_POLICY:-warn}"

# QC thresholds (ENCODE minimums)
FRIP_MIN="${ATAC_QC_FRIP_MIN:-0.1}"
TSS_ENRICH_MIN="${ATAC_QC_TSS_ENRICH_MIN:-5}"
NFR_RATIO_MIN="${ATAC_QC_NFR_RATIO_MIN:-1.0}"

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    --samtools-dir) SAMTOOLS_DIR="$2"; shift 2;;
    --peaks-bed) PEAKS_BED="$2"; shift 2;;
    --blacklist) BLACKLIST="$2"; shift 2;;
    --tss-bed) TSS_BED="$2"; shift 2;;
    --bigwig-dir) BIGWIG_DIR="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --policy) POLICY="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

# =============================================================================
# VALIDATE INPUTS
# =============================================================================
[[ "$POLICY" =~ ^(warn|exclude)$ ]] || die "Invalid --policy: $POLICY (must be warn|exclude)"
[[ -d "$SAMTOOLS_DIR" ]] || die "Missing samtools directory: $SAMTOOLS_DIR"
[[ -s "$PEAKS_BED" ]] || die "Missing peaks BED: $PEAKS_BED"
[[ -s "$BLACKLIST" ]] || die "Missing blacklist: $BLACKLIST"

need samtools
need bedtools
need awk

mkdir -p "$OUTDIR"

# =============================================================================
# AUTO-GENERATE TSS BED IF NEEDED
# =============================================================================
if [[ -z "$TSS_BED" || ! -s "$TSS_BED" ]]; then
    if [[ -s "$GTF" ]]; then
        log "Auto-generating TSS BED from GTF: $GTF"
        TSS_BED="${OUTDIR}/tss_autogen.bed"

        awk -F'\t' 'BEGIN{OFS="\t"}
          $0 !~ /^#/ && $3=="gene" {
            chr=$1; start=$4; end=$5; strand=$7; attr=$9;

            gene_id="NA"; gene_name="NA";
            if (match(attr, /gene_id "([^"]+)"/, a)) gene_id=a[1];
            if (match(attr, /gene_name "([^"]+)"/, b)) gene_name=b[1];
            if (gene_name=="NA") gene_name=gene_id;

            tss = (strand=="+") ? start : end;
            tss0 = tss - 1;

            print chr, tss0, tss0+1, gene_id, gene_name, strand
          }' "$GTF" | LC_ALL=C sort -k1,1 -k2,2n > "$TSS_BED"

        log "Generated TSS BED: $TSS_BED ($(wc -l < "$TSS_BED") genes)"
    else
        die "No TSS BED provided and GTF not found: $GTF"
    fi
fi

[[ -s "$TSS_BED" ]] || die "TSS BED is empty or missing: $TSS_BED"

# =============================================================================
# INITIALIZE OUTPUT
# =============================================================================
QC_TSV="${OUTDIR}/qc_metrics_extended.tsv"
QC_EXCLUDED="${OUTDIR}/qc_excluded_samples.txt"

# Header
printf "sample\tcondition\tfrip\ttss_enrichment\tfragment_mean\tfragment_median\tnfr_ratio\tblacklist_fraction\tpbc1\tpbc2\tfrip_pass\ttss_pass\tnfr_pass\toverall_pass\n" > "$QC_TSV"

# Clear exclusion list
> "$QC_EXCLUDED"

log "Output: $QC_TSV"
log "Thresholds: FRiP>=${FRIP_MIN}, TSS_enrich>=${TSS_ENRICH_MIN}, NFR_ratio>=${NFR_RATIO_MIN}"

# =============================================================================
# HELPER: Get sample condition
# =============================================================================
get_sample_condition() {
    local samplename="$1"
    local condition="NULL"

    # Try metadata file first
    if [[ -s "${ATAC_SAMPLE_META:-}" ]]; then
        condition=$(awk -F'\t' -v s="$samplename" 'NR>1 && $1==s {print $2; exit}' "$ATAC_SAMPLE_META")
        if [[ -n "$condition" ]]; then
            echo "$condition"
            return
        fi
    fi

    # Fallback to pattern matching
    if [[ "$samplename" =~ ATF5WT ]]; then
        condition="WT"
    elif [[ "$samplename" =~ ATF5NULL ]]; then
        condition="KO"
    fi

    echo "$condition"
}

# =============================================================================
# COMPUTE METRICS PER SAMPLE
# =============================================================================
log "Processing samples from: $SAMTOOLS_DIR"

for final_bam in "${SAMTOOLS_DIR}"/*.final.bam; do
    [[ -f "$final_bam" ]] || continue

    samplename=$(basename "$final_bam" .final.bam)
    condition=$(get_sample_condition "$samplename")

    log "[*] ${samplename} (${condition})"

    # Locate filtered BAM (pre-dedup, for PBC/blacklist metrics)
    filt_bam="${SAMTOOLS_DIR}/${samplename}.mq${ATAC_MAPQ:-20}.canon.bam"
    [[ -f "$filt_bam" ]] || {
        log "[WARN] Filtered BAM not found: $filt_bam (skipping PBC/blacklist)"
        filt_bam=""
    }

    # -------------------------------------------------------------------------
    # 1. FRiP (Fraction of Reads in Peaks)
    # -------------------------------------------------------------------------
    log "  Computing FRiP..."
    reads_in_peaks=$(bedtools intersect -u -abam "$final_bam" -b "$PEAKS_BED" | samtools view -c -)
    total_reads=$(samtools view -c "$final_bam")

    if [[ $total_reads -gt 0 ]]; then
        frip=$(awk -v p="$reads_in_peaks" -v t="$total_reads" 'BEGIN{printf "%.4f", p/t}')
    else
        frip="0.0000"
    fi

    # -------------------------------------------------------------------------
    # 2. Fragment size distribution (mean/median)
    # -------------------------------------------------------------------------
    log "  Computing fragment size distribution..."
    frag_stats_tmp="${OUTDIR}/${samplename}.frag_stats.tmp"

    samtools stats -@ "$THREADS" "$final_bam" | grep "^IS" > "$frag_stats_tmp" || true

    if [[ -s "$frag_stats_tmp" ]]; then
        # Parse insert size distribution: IS column format is: IS <size> <count>
        # Compute mean and median from distribution
        frag_mean=$(awk '{
            sum += $2 * $3;
            count += $3;
        } END {
            if (count > 0) printf "%.1f", sum/count; else print "NA"
        }' "$frag_stats_tmp")

        frag_median=$(awk '{
            sizes[NR] = $2;
            counts[NR] = $3;
            total += $3;
        } END {
            if (total == 0) { print "NA"; exit }
            half = total / 2;
            cum = 0;
            for (i = 1; i <= NR; i++) {
                cum += counts[i];
                if (cum >= half) {
                    print sizes[i];
                    exit;
                }
            }
        }' "$frag_stats_tmp")
    else
        frag_mean="NA"
        frag_median="NA"
    fi

    rm -f "$frag_stats_tmp"

    # -------------------------------------------------------------------------
    # 3. NFR ratio (nucleosome-free vs mononucleosome)
    # -------------------------------------------------------------------------
    log "  Computing NFR ratio..."

    if [[ -s "$frag_stats_tmp.is" ]]; then
        rm -f "$frag_stats_tmp.is"
    fi

    samtools stats -@ "$THREADS" "$final_bam" | grep "^IS" > "$frag_stats_tmp.is" || true

    if [[ -s "$frag_stats_tmp.is" ]]; then
        nfr_count=$(awk '$2 < 150 {sum += $3} END {print sum+0}' "$frag_stats_tmp.is")
        mono_count=$(awk '$2 >= 150 && $2 <= 300 {sum += $3} END {print sum+0}' "$frag_stats_tmp.is")

        if [[ $mono_count -gt 0 ]]; then
            nfr_ratio=$(awk -v n="$nfr_count" -v m="$mono_count" 'BEGIN{printf "%.4f", n/m}')
        else
            nfr_ratio="0.0000"
        fi
    else
        nfr_ratio="NA"
    fi

    rm -f "$frag_stats_tmp.is"

    # -------------------------------------------------------------------------
    # 4. Blacklist fraction (on filtered pre-dedup BAM)
    # -------------------------------------------------------------------------
    log "  Computing blacklist fraction..."

    if [[ -n "$filt_bam" && -f "$filt_bam" ]]; then
        bl_reads=$(bedtools intersect -u -abam "$filt_bam" -b "$BLACKLIST" | samtools view -c -)
        filt_total=$(samtools view -c "$filt_bam")

        if [[ $filt_total -gt 0 ]]; then
            bl_frac=$(awk -v b="$bl_reads" -v t="$filt_total" 'BEGIN{printf "%.6f", b/t}')
        else
            bl_frac="0.000000"
        fi
    else
        bl_frac="NA"
    fi

    # -------------------------------------------------------------------------
    # 5. PBC1 (library complexity: distinct locations / total reads)
    # -------------------------------------------------------------------------
    log "  Computing PBC1..."

    if [[ -n "$filt_bam" && -f "$filt_bam" ]]; then
        pbc_tmp="${OUTDIR}/${samplename}.pbc.tmp"

        # Extract genomic positions (chr:start-end), count unique
        bedtools bamtobed -i "$filt_bam" 2>/dev/null \
            | cut -f1-3 \
            | LC_ALL=C sort \
            | LC_ALL=C uniq -c \
            > "$pbc_tmp" || true

        if [[ -s "$pbc_tmp" ]]; then
            n_distinct=$(wc -l < "$pbc_tmp")
            n_total=$(awk '{sum += $1} END {print sum+0}' "$pbc_tmp")

            if [[ $n_total -gt 0 ]]; then
                pbc1=$(awk -v d="$n_distinct" -v t="$n_total" 'BEGIN{printf "%.4f", d/t}')
            else
                pbc1="0.0000"
            fi
        else
            pbc1="NA"
        fi

        rm -f "$pbc_tmp"
    else
        pbc1="NA"
    fi

    # -------------------------------------------------------------------------
    # 6. PBC2 (N_1_read_locations / N_2_read_locations)
    # -------------------------------------------------------------------------
    log "  Computing PBC2..."

    if [[ -n "$filt_bam" && -f "$filt_bam" ]]; then
        pbc_tmp="${OUTDIR}/${samplename}.pbc.tmp"

        bedtools bamtobed -i "$filt_bam" 2>/dev/null \
            | cut -f1-3 \
            | LC_ALL=C sort \
            | LC_ALL=C uniq -c \
            > "$pbc_tmp" || true

        if [[ -s "$pbc_tmp" ]]; then
            n_one=$(awk '$1 == 1 {count++} END {print count+0}' "$pbc_tmp")
            n_two=$(awk '$1 == 2 {count++} END {print count+0}' "$pbc_tmp")

            if [[ $n_two -gt 0 ]]; then
                pbc2=$(awk -v o="$n_one" -v t="$n_two" 'BEGIN{printf "%.4f", o/t}')
            else
                pbc2="NA"
            fi
        else
            pbc2="NA"
        fi

        rm -f "$pbc_tmp"
    else
        pbc2="NA"
    fi

    # -------------------------------------------------------------------------
    # 7. TSS enrichment (simplified: reads near TSS vs flanks)
    # -------------------------------------------------------------------------
    log "  Computing TSS enrichment..."

    tss_tmp="${OUTDIR}/${samplename}.tss.tmp"

    # TSS regions: +/- 100bp around TSS
    bedtools slop -i "$TSS_BED" -g "${ATAC_CHROM_SIZES}" -b 100 2>/dev/null \
        | bedtools intersect -u -abam "$final_bam" -b - \
        | samtools view -c - > "$tss_tmp.center" || echo "0" > "$tss_tmp.center"

    # Flanking regions: 1000bp upstream and downstream (exclude center +/- 100bp)
    # Upstream: TSS - 1000 to TSS - 100
    # Downstream: TSS + 100 to TSS + 1000
    awk 'BEGIN{OFS="\t"}{
        chr=$1; tss_start=$2; tss_end=$3; id=$4; name=$5; strand=$6;

        # TSS is at tss_start (0-based); upstream/downstream are strand-aware
        if (strand == "+") {
            up_start = tss_start - 1000;
            up_end = tss_start - 100;
            down_start = tss_start + 100;
            down_end = tss_start + 1000;
        } else {
            up_start = tss_start + 100;
            up_end = tss_start + 1000;
            down_start = tss_start - 1000;
            down_end = tss_start - 100;
        }

        if (up_start < 0) up_start = 0;
        if (up_end < 0) up_end = 0;
        if (down_start < 0) down_start = 0;
        if (down_end < 0) down_end = 0;

        if (up_end > up_start) print chr, up_start, up_end, id, name, strand;
        if (down_end > down_start) print chr, down_start, down_end, id, name, strand;
    }' "$TSS_BED" | LC_ALL=C sort -k1,1 -k2,2n \
        | bedtools intersect -u -abam "$final_bam" -b - \
        | samtools view -c - > "$tss_tmp.flanks" || echo "0" > "$tss_tmp.flanks"

    center_reads=$(cat "$tss_tmp.center")
    flank_reads=$(cat "$tss_tmp.flanks")

    if [[ $flank_reads -gt 0 ]]; then
        # TSS enrichment = (center reads / center bp) / (flank reads / flank bp)
        # Center: 200bp per TSS, Flanks: 1800bp per TSS
        n_tss=$(wc -l < "$TSS_BED")
        center_bp=$((n_tss * 200))
        flank_bp=$((n_tss * 1800))

        tss_enrichment=$(awk -v c="$center_reads" -v cb="$center_bp" \
                             -v f="$flank_reads" -v fb="$flank_bp" \
                             'BEGIN{
                                 if (fb > 0 && cb > 0) {
                                     center_density = c / cb;
                                     flank_density = f / fb;
                                     if (flank_density > 0)
                                         printf "%.4f", center_density / flank_density;
                                     else
                                         print "0.0000";
                                 } else print "0.0000";
                             }')
    else
        tss_enrichment="0.0000"
    fi

    rm -f "$tss_tmp.center" "$tss_tmp.flanks"

    # -------------------------------------------------------------------------
    # 8. Pass/fail flags
    # -------------------------------------------------------------------------
    frip_pass="FALSE"
    tss_pass="FALSE"
    nfr_pass="FALSE"

    if [[ "$frip" != "NA" ]]; then
        if awk -v f="$frip" -v m="$FRIP_MIN" 'BEGIN{exit (f >= m ? 0 : 1)}'; then
            frip_pass="TRUE"
        fi
    fi

    if [[ "$tss_enrichment" != "NA" ]]; then
        if awk -v t="$tss_enrichment" -v m="$TSS_ENRICH_MIN" 'BEGIN{exit (t >= m ? 0 : 1)}'; then
            tss_pass="TRUE"
        fi
    fi

    if [[ "$nfr_ratio" != "NA" ]]; then
        if awk -v n="$nfr_ratio" -v m="$NFR_RATIO_MIN" 'BEGIN{exit (n >= m ? 0 : 1)}'; then
            nfr_pass="TRUE"
        fi
    fi

    overall_pass="TRUE"
    if [[ "$frip_pass" == "FALSE" || "$tss_pass" == "FALSE" || "$nfr_pass" == "FALSE" ]]; then
        overall_pass="FALSE"
    fi

    # -------------------------------------------------------------------------
    # 9. Write row to output
    # -------------------------------------------------------------------------
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$samplename" "$condition" "$frip" "$tss_enrichment" \
        "$frag_mean" "$frag_median" "$nfr_ratio" "$bl_frac" \
        "$pbc1" "$pbc2" \
        "$frip_pass" "$tss_pass" "$nfr_pass" "$overall_pass" >> "$QC_TSV"

    # -------------------------------------------------------------------------
    # 10. Apply policy
    # -------------------------------------------------------------------------
    if [[ "$overall_pass" == "FALSE" ]]; then
        if [[ "$POLICY" == "warn" ]]; then
            log "[WARN] Sample failed QC: $samplename (FRiP=${frip}, TSS=${tss_enrichment}, NFR=${nfr_ratio})"
        elif [[ "$POLICY" == "exclude" ]]; then
            log "[EXCLUDE] Sample failed QC: $samplename (FRiP=${frip}, TSS=${tss_enrichment}, NFR=${nfr_ratio})"
            echo "$samplename" >> "$QC_EXCLUDED"
        fi
    else
        log "  PASS (FRiP=${frip}, TSS=${tss_enrichment}, NFR=${nfr_ratio})"
    fi

done

# =============================================================================
# SUMMARY
# =============================================================================
log "[DONE] Extended QC metrics:"
log "  Output: $QC_TSV"

n_samples=$(tail -n +2 "$QC_TSV" | wc -l)
n_passed=$(awk -F'\t' 'NR>1 && $NF=="TRUE" {count++} END{print count+0}' "$QC_TSV")
n_failed=$((n_samples - n_passed))

log "  Samples: $n_samples total, $n_passed passed, $n_failed failed"

if [[ "$POLICY" == "exclude" && $n_failed -gt 0 ]]; then
    log "  Excluded samples list: $QC_EXCLUDED"
fi
