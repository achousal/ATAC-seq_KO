#!/usr/bin/env bash
# =============================================================================
# idr_analysis.sh - ENCODE-faithful IDR analysis for ATAC-seq
# =============================================================================
# Implements ENCODE IDR protocol with pooled pseudoreplicates.
# Falls back to union peaks if idr command not available.
#
# Usage:
#   bash idr_analysis.sh [OPTIONS]
#
# Options:
#   --macs-dir DIR      MACS output directory (default: analysis/macs)
#   --samtools-dir DIR  Samtools BAM directory (default: analysis/samtools)
#   --threads N         Thread count (default: from config or 8)
#   --outdir DIR        IDR output directory (default: macs-dir/idr)
#
# Outputs:
#   WT_optimal_peaks.bed, KO_optimal_peaks.bed
#   idr_consensus_peaks.bed, idr_consensus_peaks.saf
#   idr_summary.tsv
#   idr_method.txt
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

log(){ printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*"; }
skip(){ log "[SKIP] $*"; }
die(){ log "[ERROR] $*"; exit 1; }
have(){ command -v "$1" >/dev/null 2>&1; }

# =============================================================================
# LOAD CONFIGURATION
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.sh"

if [[ -f "$CONFIG_FILE" ]]; then
    log "Loading config: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    die "Config not found: $CONFIG_FILE"
fi

# =============================================================================
# PARSE CLI ARGS
# =============================================================================
MACS_DIR="${SCRIPT_DIR}/macs"
SAMTOOLS_DIR="${SCRIPT_DIR}/samtools"
THREADS="${ATAC_THREADS:-8}"
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --macs-dir)
            MACS_DIR="$2"
            shift 2
            ;;
        --samtools-dir)
            SAMTOOLS_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

# Default outdir
if [[ -z "$OUTDIR" ]]; then
    OUTDIR="${MACS_DIR}/idr"
fi

mkdir -p "$OUTDIR"

# =============================================================================
# VALIDATE CONFIG
# =============================================================================
IDR_THRESHOLD="${ATAC_IDR_THRESHOLD:-0.05}"
IDR_MIN_REPS="${ATAC_IDR_MIN_REPS:-2}"
IDR_MODE="${ATAC_IDR_MODE:-encode}"
WT_PATTERN="${ATAC_WT_PATTERN:-ATF5WT}"
KO_PATTERN="${ATAC_KO_PATTERN:-ATF5NULL}"
GENOME_SIZE_MACS="${ATAC_GENOME_SIZE:-1.87e9}"

if [[ -z "${MACS_BIN:-}" ]]; then
    die "MACS_BIN not set in config"
fi

# =============================================================================
# CHECK IDR AVAILABILITY
# =============================================================================
if ! have idr; then
    log "[WARNING] idr command not found. Falling back to union peaks."
    echo "fallback" > "${OUTDIR}/idr_method.txt"
    echo "IDR command not available; using union peaks as fallback." >> "${OUTDIR}/idr_method.txt"
    exit 0
fi

# =============================================================================
# CHECKPOINT
# =============================================================================
FINAL_OUTPUT="${OUTDIR}/idr_consensus_peaks.saf"
if [[ -f "$FINAL_OUTPUT" ]]; then
    skip "IDR consensus already exists: $FINAL_OUTPUT"
    exit 0
fi

log "Starting ENCODE IDR analysis (mode: $IDR_MODE, threshold: $IDR_THRESHOLD)"

# =============================================================================
# COLLECT REPLICATES PER CONDITION
# =============================================================================
declare -a WT_BAMS KO_BAMS

while IFS= read -r bam; do
    bam_name=$(basename "$bam" .final.bam)
    if [[ "$bam_name" =~ $WT_PATTERN ]]; then
        WT_BAMS+=("$bam")
    elif [[ "$bam_name" =~ $KO_PATTERN ]]; then
        KO_BAMS+=("$bam")
    fi
done < <(find "$SAMTOOLS_DIR" -name "*.final.bam" | sort)

log "Found ${#WT_BAMS[@]} WT replicates, ${#KO_BAMS[@]} KO replicates"

if [[ ${#WT_BAMS[@]} -lt $IDR_MIN_REPS ]] && [[ ${#KO_BAMS[@]} -lt $IDR_MIN_REPS ]]; then
    log "[WARNING] Insufficient replicates for both conditions (min $IDR_MIN_REPS). Falling back to union peaks."
    echo "fallback" > "${OUTDIR}/idr_method.txt"
    echo "Insufficient replicates: WT=${#WT_BAMS[@]}, KO=${#KO_BAMS[@]} (min $IDR_MIN_REPS required)" >> "${OUTDIR}/idr_method.txt"
    exit 0
fi

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Pool BAMs
pool_bams() {
    local -n bams=$1
    local out_bam=$2
    log "Pooling ${#bams[@]} BAMs -> $out_bam"
    samtools merge -@ "$THREADS" -f "$out_bam" "${bams[@]}"
    samtools index -@ "$THREADS" "$out_bam"
}

# Create pseudoreplicates (50/50 split)
create_pseudoreps() {
    local pool_bam=$1
    local pr1=$2
    local pr2=$3
    log "Creating pseudoreplicates from $pool_bam"

    # Count total reads
    local total_reads=$(samtools view -c "$pool_bam")
    local half=$((total_reads / 2))

    # Shuffle and split
    samtools view -h "$pool_bam" | \
        awk -v seed="$RANDOM" 'BEGIN{srand(seed)} /^@/{print; next} {print rand()" "$0}' | \
        sort -k1,1n | \
        cut -d' ' -f2- | \
        awk -v half="$half" 'NR<=half' | \
        samtools view -@ "$THREADS" -bS - > "$pr1"

    samtools view -h "$pool_bam" | \
        awk -v seed="$RANDOM" 'BEGIN{srand(seed)} /^@/{print; next} {print rand()" "$0}' | \
        sort -k1,1n | \
        cut -d' ' -f2- | \
        awk -v half="$half" 'NR>half' | \
        samtools view -@ "$THREADS" -bS - > "$pr2"

    samtools index -@ "$THREADS" "$pr1"
    samtools index -@ "$THREADS" "$pr2"
}

# Call peaks with MACS3
call_peaks() {
    local bam=$1
    local name=$2
    local outdir=$3

    log "Calling peaks: $name"
    mkdir -p "$outdir"

    "$MACS_BIN" callpeak \
        -t "$bam" \
        -n "$name" \
        --outdir "$outdir" \
        -f BAMPE \
        -g "$GENOME_SIZE_MACS" \
        --nomodel \
        --keep-dup all \
        --slocal 1000 \
        2>&1 | grep -v "INFO" || true
}

# Run IDR
run_idr() {
    local peak1=$1
    local peak2=$2
    local out=$3

    log "Running IDR: $(basename $peak1) vs $(basename $peak2)"

    idr --samples "$peak1" "$peak2" \
        --input-file-type narrowPeak \
        --rank p.value \
        --output-file "$out" \
        --plot \
        --idr-threshold "$IDR_THRESHOLD" \
        --soft-idr-threshold 0.1 \
        --log-output-file "${out}.log" \
        2>&1 | grep -v "INFO" || true
}

# Count peaks passing IDR
count_idr_peaks() {
    local idr_file=$1
    awk -v thresh="$IDR_THRESHOLD" '$5 >= -log(thresh)/log(10)' "$idr_file" | wc -l
}

# Get top N peaks from pooled peak file
get_top_peaks() {
    local pooled_peaks=$1
    local n=$2
    local out=$3

    log "Extracting top $n peaks from pooled set"
    sort -k8,8nr "$pooled_peaks" | head -n "$n" > "$out"
}

# =============================================================================
# ENCODE IDR PROTOCOL (per condition)
# =============================================================================

process_condition_encode() {
    local condition=$1
    local -n rep_bams=$2
    local workdir="${OUTDIR}/${condition}"

    mkdir -p "$workdir"

    log "=== Processing condition: $condition (ENCODE mode) ==="

    # Check replicate count
    if [[ ${#rep_bams[@]} -lt $IDR_MIN_REPS ]]; then
        log "[WARNING] Insufficient replicates for $condition (${#rep_bams[@]} < $IDR_MIN_REPS). Skipping IDR."
        return 1
    fi

    # Step 1: Pool replicates
    local pooled_bam="${workdir}/pooled.bam"
    pool_bams rep_bams "$pooled_bam"

    # Step 2: Call peaks on pooled
    local pooled_peaks_dir="${workdir}/pooled_peaks"
    call_peaks "$pooled_bam" "pooled" "$pooled_peaks_dir"
    local pooled_peaks="${pooled_peaks_dir}/pooled_peaks.narrowPeak"

    # Step 3: Create pseudoreplicates
    local pr1_bam="${workdir}/pseudorep1.bam"
    local pr2_bam="${workdir}/pseudorep2.bam"
    create_pseudoreps "$pooled_bam" "$pr1_bam" "$pr2_bam"

    # Step 4: Call peaks on pseudoreplicates
    call_peaks "$pr1_bam" "pseudorep1" "${workdir}/pr1_peaks"
    call_peaks "$pr2_bam" "pseudorep2" "${workdir}/pr2_peaks"

    local pr1_peaks="${workdir}/pr1_peaks/pseudorep1_peaks.narrowPeak"
    local pr2_peaks="${workdir}/pr2_peaks/pseudorep2_peaks.narrowPeak"

    # Step 5: True replicate IDR (all pairwise combinations)
    local idr_dir="${workdir}/idr_true_reps"
    mkdir -p "$idr_dir"

    local -a idr_files
    for ((i=0; i<${#rep_bams[@]}-1; i++)); do
        for ((j=i+1; j<${#rep_bams[@]}; j++)); do
            local rep1_name=$(basename "${rep_bams[$i]}" .final.bam)
            local rep2_name=$(basename "${rep_bams[$j]}" .final.bam)

            # Need per-rep peaks
            local rep1_peaks="${MACS_DIR}/rep_peaks_bl/${rep1_name}.bl.narrowPeak"
            local rep2_peaks="${MACS_DIR}/rep_peaks_bl/${rep2_name}.bl.narrowPeak"

            if [[ ! -f "$rep1_peaks" ]] || [[ ! -f "$rep2_peaks" ]]; then
                log "[WARNING] Per-rep peaks not found for $rep1_name or $rep2_name. Skipping pair."
                continue
            fi

            local idr_out="${idr_dir}/idr_${rep1_name}_vs_${rep2_name}.txt"
            run_idr "$rep1_peaks" "$rep2_peaks" "$idr_out"
            idr_files+=("$idr_out")
        done
    done

    if [[ ${#idr_files[@]} -eq 0 ]]; then
        log "[ERROR] No valid true replicate IDR comparisons for $condition"
        return 1
    fi

    # Max Nt across all pairwise comparisons
    local Nt=0
    for idr_file in "${idr_files[@]}"; do
        local n=$(count_idr_peaks "$idr_file")
        if [[ $n -gt $Nt ]]; then
            Nt=$n
        fi
    done

    log "True replicate IDR: Nt = $Nt"

    # Step 6: Pooled pseudoreplicate IDR
    local pr_idr="${workdir}/idr_pseudoreps.txt"
    run_idr "$pr1_peaks" "$pr2_peaks" "$pr_idr"

    local Np=$(count_idr_peaks "$pr_idr")
    log "Pooled pseudoreplicate IDR: Np = $Np"

    # Step 7: ENCODE decision rule: max(Nt, Np)
    local optimal_n=$Nt
    if [[ $Np -gt $Nt ]]; then
        optimal_n=$Np
    fi

    log "Optimal peak count: max(Nt=$Nt, Np=$Np) = $optimal_n"

    # Step 8: Rescue ratio
    local rescue_ratio=0
    if [[ $Nt -gt 0 ]]; then
        rescue_ratio=$(awk -v np="$Np" -v nt="$Nt" 'BEGIN{printf "%.2f", np/nt}')
    fi

    if (( $(awk -v rr="$rescue_ratio" 'BEGIN{print (rr > 2.0)}') )); then
        log "[WARNING] High rescue ratio ($rescue_ratio > 2.0) for $condition. Possible quality issue."
    fi

    # Step 9: Extract top optimal_n peaks from pooled
    local optimal_peaks="${OUTDIR}/${condition}_optimal_peaks.bed"
    get_top_peaks "$pooled_peaks" "$optimal_n" "$optimal_peaks"

    # Write summary
    echo -e "${condition}\t${Nt}\t${Np}\t${optimal_n}\t${rescue_ratio}" >> "${OUTDIR}/idr_summary.tsv"

    log "=== $condition complete: $optimal_n optimal peaks ==="
    return 0
}

# =============================================================================
# PAIRWISE MODE (simplified)
# =============================================================================

process_condition_pairwise() {
    local condition=$1
    local -n rep_bams=$2
    local workdir="${OUTDIR}/${condition}"

    mkdir -p "$workdir"

    log "=== Processing condition: $condition (pairwise mode) ==="

    if [[ ${#rep_bams[@]} -lt $IDR_MIN_REPS ]]; then
        log "[WARNING] Insufficient replicates for $condition (${#rep_bams[@]} < $IDR_MIN_REPS). Skipping IDR."
        return 1
    fi

    # Pairwise IDR on true replicates
    local idr_dir="${workdir}/idr_pairwise"
    mkdir -p "$idr_dir"

    declare -A peak_pass_count

    for ((i=0; i<${#rep_bams[@]}-1; i++)); do
        for ((j=i+1; j<${#rep_bams[@]}; j++)); do
            local rep1_name=$(basename "${rep_bams[$i]}" .final.bam)
            local rep2_name=$(basename "${rep_bams[$j]}" .final.bam)

            local rep1_peaks="${MACS_DIR}/rep_peaks_bl/${rep1_name}.bl.narrowPeak"
            local rep2_peaks="${MACS_DIR}/rep_peaks_bl/${rep2_name}.bl.narrowPeak"

            if [[ ! -f "$rep1_peaks" ]] || [[ ! -f "$rep2_peaks" ]]; then
                log "[WARNING] Per-rep peaks not found. Skipping."
                continue
            fi

            local idr_out="${idr_dir}/idr_${rep1_name}_vs_${rep2_name}.txt"
            run_idr "$rep1_peaks" "$rep2_peaks" "$idr_out"

            # Mark peaks passing IDR
            awk -v thresh="$IDR_THRESHOLD" '$5 >= -log(thresh)/log(10) {print $1":"$2"-"$3}' "$idr_out" | \
                while read peak_id; do
                    peak_pass_count["$peak_id"]=$((${peak_pass_count["$peak_id"]:-0} + 1))
                done
        done
    done

    # Filter peaks passing in >=2 comparisons
    local consensus="${OUTDIR}/${condition}_optimal_peaks.bed"
    > "$consensus"

    for peak_id in "${!peak_pass_count[@]}"; do
        if [[ ${peak_pass_count[$peak_id]} -ge 2 ]]; then
            IFS=':' read -r chr rest <<< "$peak_id"
            IFS='-' read -r start end <<< "$rest"
            echo -e "${chr}\t${start}\t${end}\t${peak_id}\t${peak_pass_count[$peak_id]}" >> "$consensus"
        fi
    done

    local n_consensus=$(wc -l < "$consensus")
    log "=== $condition complete: $n_consensus consensus peaks ==="

    echo -e "${condition}\tNA\tNA\t${n_consensus}\tNA" >> "${OUTDIR}/idr_summary.tsv"
    return 0
}

# =============================================================================
# PROCESS BOTH CONDITIONS
# =============================================================================

echo -e "condition\tNt\tNp\tmax_Nt_Np\trescue_ratio" > "${OUTDIR}/idr_summary.tsv"

if [[ "$IDR_MODE" == "encode" ]]; then
    process_condition_encode "WT" WT_BAMS || log "[WARNING] WT IDR failed or skipped"
    process_condition_encode "KO" KO_BAMS || log "[WARNING] KO IDR failed or skipped"
elif [[ "$IDR_MODE" == "pairwise" ]]; then
    process_condition_pairwise "WT" WT_BAMS || log "[WARNING] WT IDR failed or skipped"
    process_condition_pairwise "KO" KO_BAMS || log "[WARNING] KO IDR failed or skipped"
else
    die "Unknown IDR_MODE: $IDR_MODE (must be 'encode' or 'pairwise')"
fi

# =============================================================================
# CREATE CONSENSUS PEAK SET
# =============================================================================

log "Creating IDR consensus peak set"

WT_OPTIMAL="${OUTDIR}/WT_optimal_peaks.bed"
KO_OPTIMAL="${OUTDIR}/KO_optimal_peaks.bed"

if [[ ! -f "$WT_OPTIMAL" ]] && [[ ! -f "$KO_OPTIMAL" ]]; then
    die "No optimal peaks generated for either condition. IDR failed."
fi

# Union of WT + KO optimal peaks
CONSENSUS_BED="${OUTDIR}/idr_consensus_peaks.bed"
cat "$WT_OPTIMAL" "$KO_OPTIMAL" 2>/dev/null | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - | \
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak_"NR}' > "$CONSENSUS_BED"

log "Consensus peaks: $(wc -l < $CONSENSUS_BED)"

# =============================================================================
# CONVERT TO SAF FORMAT
# =============================================================================

CONSENSUS_SAF="${OUTDIR}/idr_consensus_peaks.saf"
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
     {print $4,$1,$2+1,$3,"+"}' "$CONSENSUS_BED" > "$CONSENSUS_SAF"

log "Created SAF: $CONSENSUS_SAF"

# =============================================================================
# DOCUMENT METHOD
# =============================================================================

METHOD_FILE="${OUTDIR}/idr_method.txt"
cat > "$METHOD_FILE" <<EOF
ENCODE IDR Analysis
===================
Mode: $IDR_MODE
IDR threshold: $IDR_THRESHOLD
Min replicates: $IDR_MIN_REPS

Decision rule (ENCODE):
- Nt = peaks passing IDR in true replicate comparisons
- Np = peaks passing IDR in pooled pseudoreplicate comparison
- Optimal peak count per condition = max(Nt, Np)
- Rescue ratio = Np/Nt (warn if > 2.0)

Consensus: union of WT + KO optimal peaks, merged and assigned IDs.
EOF

log "IDR analysis complete. Outputs in: $OUTDIR"

# =============================================================================
# CLEANUP
# =============================================================================

if [[ "${ATAC_KEEP_INTERMEDIATES:-0}" -eq 0 ]]; then
    log "Cleaning up intermediate BAMs"
    rm -f "${OUTDIR}"/*/pooled.bam* "${OUTDIR}"/*/pseudorep*.bam*
fi

log "Done."
