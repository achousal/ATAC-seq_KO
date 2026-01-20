#!/usr/bin/env bash
# Andres Chousal Cantu
# Icahn School of Medicine
#
# atacseq_analysis.sh
# ATAC-seq pipeline (mm10) with:
#   - Nextera adapter trimming
#   - Read filters: canonical chr only + MAPQ>=20 + proper pairs
#   - BigWigs via deepTools bamCoverage --normalizeUsing CPM
#   - Per-replicate MACS3 peaks -> blacklist filter -> union peaks -> featureCounts -> DESeq2
#   - HOMER: union peak annotation + motif enrichment
#   - MultiQC aggregation
#   - QC table: bowtie2 overall alignment rate, mapped reads, chrM frac (pre-filter),
#               filtered mapped, final mapped, Picard dup rate, est library size
#
# Raw FASTQs: MEF-ATF5WT-n{2,3,4}_R{1,2}_001.fastq.gz
#             MEF-ATF5NULL-n{2,3,4}_R{1,2}_001.fastq.gz

set -euo pipefail
IFS=$'\n\t'

log(){ printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*"; }
skip(){ log "[SKIP] $*"; }
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
    log "[ERROR] Config not found: $CONFIG_FILE"
    exit 1
fi

# Validate required config vars
for var in ATAC_BASE ATAC_RAW_FASTQ ATAC_BLACKLIST ATAC_BOWTIE_INDEX ATAC_CHROM_SIZES; do
    if [[ -z "${!var:-}" ]]; then
        log "[ERROR] Required config var not set: $var"
        exit 1
    fi
done

# =============================================================================
# OPTIONAL: modules + conda
# =============================================================================
if have module; then
    module load anaconda 2>/dev/null || true
fi

if have conda; then
    eval "$(conda shell.bash hook)"
    conda activate "${ATAC_CONDA_ENV:-csb}"
else
    echo "[ERROR] conda not found in PATH. Load anaconda/miniconda first." >&2
    exit 1
fi

# =============================================================================
# DERIVED PATHS FROM CONFIG
# =============================================================================
THREADS="${ATAC_THREADS:-8}"
MAPQ="${ATAC_MAPQ:-20}"
KEEP_INTERMEDIATES="${ATAC_KEEP_INTERMEDIATES:-0}"
MACS_BIN="${MACS_BIN}"
GENOME_SIZE_MACS="${ATAC_GENOME_SIZE:-1.87e9}"

raw_atac="$ATAC_RAW_FASTQ"
raw_dir="$ATAC_REFERENCE_DIR"
bl_regions="$ATAC_BLACKLIST"
bowtie_index="$ATAC_BOWTIE_INDEX"
chrom_sizes="$ATAC_CHROM_SIZES"

analysis_dir="${ATAC_BASE}/analysis"
fastqc_dir="${analysis_dir}/fastqc"
trimg_dir="${analysis_dir}/trimmgalore"
bowtie2_dir="${analysis_dir}/bowtie2"
samtools_dir="${analysis_dir}/samtools"
picard_dir="${analysis_dir}/picard"
deeptools_dir="${analysis_dir}/deeptools"
macs_dir="${analysis_dir}/macs"
subread_dir="${analysis_dir}/subread"
homer_dir="${analysis_dir}/homer"
deseq_dir="${subread_dir}/deseq"
multiqc_dir="${analysis_dir}/multiqc"
log_dir="${analysis_dir}/logs"

deseq_R="${SCRIPT_DIR}/deseq_atac.R"
sample_meta="${ATAC_SAMPLE_META:-${SCRIPT_DIR}/samples.tsv}"

mkdir -p "$fastqc_dir" "$trimg_dir" "$bowtie2_dir" "$samtools_dir" \
         "$picard_dir" "$deeptools_dir" "$macs_dir" "$subread_dir" \
         "$homer_dir" "$deseq_dir" "$multiqc_dir" "$log_dir"

# =============================================================================
# SANITY CHECKS
# =============================================================================
[[ -s "$chrom_sizes" ]] || { log "[ERROR] Missing $chrom_sizes"; exit 1; }
[[ -s "$bl_regions"  ]] || { log "[ERROR] Missing $bl_regions"; exit 1; }
[[ -d "$raw_atac"    ]] || { log "[ERROR] Missing FASTQ dir $raw_atac"; exit 1; }
[[ -x "$MACS_BIN"    ]] || { log "[ERROR] MACS3 not executable at $MACS_BIN"; exit 1; }

# Genome version check (warn if index path doesn't contain mm10)
if [[ ! "$bowtie_index" =~ mm10 ]]; then
    log "[WARN] Bowtie index path does not contain 'mm10' - verify correct genome"
fi

# Canonical chromosomes filter (mouse): chr1-19, chrX, chrY
awk_keep_canonical='BEGIN{OFS="\t"} /^@/ {print; next} $3 ~ /^chr([0-9]+|X|Y)$/ {print}'

# =============================================================================
# TOOL PRESENCE
# =============================================================================
req_tools=(fastqc trim_galore bowtie2 samtools bedtools featureCounts bamCoverage multiqc)
for t in "${req_tools[@]}"; do
    have "$t" || { log "[ERROR] Required tool not found in PATH: $t"; exit 1; }
done

# Picard: either wrapper "picard" OR $PICARD jar
if ! have picard && [[ -z "${PICARD:-}" ]]; then
    log "[ERROR] Picard not found (no 'picard' wrapper and \$PICARD not set)"
    exit 1
fi

# =============================================================================
# RECORD VERSIONS
# =============================================================================
versions_txt="${log_dir}/versions.txt"
if [[ ! -s "$versions_txt" ]]; then
    {
        echo "=== Versions ($(date)) ==="
        echo "Config: $CONFIG_FILE"
        conda --version || true
        fastqc --version 2>&1 | head -n1 || true
        trim_galore --version 2>&1 | head -n2 || true
        bowtie2 --version 2>&1 | head -n1 || true
        samtools --version 2>&1 | head -n2 || true
        bedtools --version 2>&1 | head -n1 || true
        featureCounts -v 2>&1 | head -n1 || true
        bamCoverage --version 2>&1 | head -n1 || true
        multiqc --version 2>&1 | head -n1 || true
        if have picard; then
            picard MarkDuplicates --version 2>&1 | head -n1 || true
        else
            java -jar "$PICARD" MarkDuplicates --version 2>&1 | head -n1 || true
        fi
        echo "MACS3: $("$MACS_BIN" --version 2>&1 | head -n1)"
    } > "$versions_txt"
fi

# =============================================================================
# QC TABLE
# Note: picard_percent_dup is computed AFTER MAPQ/proper-pair/canonical filtering
# =============================================================================
qc_tsv="${subread_dir}/qc_metrics.tsv"
if [[ ! -s "$qc_tsv" ]]; then
    printf "sample\tcondition\treplicate\tbowtie2_overall_alignment_rate\taligned_mapped\taligned_chrM_mapped\taligned_chrM_frac\tfiltered_mapped\tfinal_mapped\tpicard_percent_dup_postfilter\tpicard_est_library_size\n" > "$qc_tsv"
fi

# =============================================================================
# HELPER: Get sample metadata
# =============================================================================
get_sample_condition() {
    local samplename="$1"
    local condition="NULL"

    # Try metadata file first
    if [[ -s "$sample_meta" ]]; then
        condition=$(awk -F'\t' -v s="$samplename" 'NR>1 && $1==s {print $2; exit}' "$sample_meta")
        if [[ -n "$condition" ]]; then
            echo "$condition"
            return
        fi
    fi

    # Fallback to pattern matching
    if [[ "$samplename" =~ ${ATAC_WT_PATTERN:-ATF5WT} ]]; then
        condition="WT"
    fi
    echo "$condition"
}

get_sample_replicate() {
    local samplename="$1"
    local replicate="NA"

    # Try metadata file first
    if [[ -s "$sample_meta" ]]; then
        replicate=$(awk -F'\t' -v s="$samplename" 'NR>1 && $1==s {print $3; exit}' "$sample_meta")
        if [[ -n "$replicate" ]]; then
            echo "$replicate"
            return
        fi
    fi

    # Fallback to pattern extraction
    replicate=$(echo "$samplename" | grep -oE 'n[0-9]+' || echo "NA")
    echo "$replicate"
}

# =============================================================================
# STEP 1: PREPROCESS -> ALIGN -> FILTER -> DEDUP -> BW
# =============================================================================
log "STEP 1: Preprocess/Align/Filter/Dedup/BigWig"
found_any=0

for r1 in "${raw_atac}"/*_R1_001.fastq.gz; do
    [[ -e "$r1" ]] || continue
    found_any=1
    r2="${r1/_R1_/_R2_}"
    [[ -f "$r2" ]] || { log "Missing mate for $(basename "$r1")"; continue; }

    samplename=$(basename "$r1" _R1_001.fastq.gz)

    tr1="${trimg_dir}/${samplename}_R1_001_val_1.fq.gz"
    tr2="${trimg_dir}/${samplename}_R2_001_val_2.fq.gz"
    sam="${bowtie2_dir}/${samplename}.sam"

    aligned_bam="${samtools_dir}/${samplename}.aligned.sorted.bam"
    filt_bam="${samtools_dir}/${samplename}.mq${MAPQ}.canon.bam"
    final_bam="${samtools_dir}/${samplename}.final.bam"

    log "[*] ${samplename}"

    # FastQC (raw)
    fq1_html="${fastqc_dir}/$(basename "$r1" .fastq.gz)_fastqc.html"
    fq2_html="${fastqc_dir}/$(basename "$r2" .fastq.gz)_fastqc.html"
    if [[ ! -s "$fq1_html" || ! -s "$fq2_html" ]]; then
        fastqc --threads "$THREADS" -o "$fastqc_dir" "$r1" "$r2"
    else
        skip "FastQC exists: ${samplename}"
    fi

    # Trim Galore (with Nextera adapters for ATAC-seq)
    if [[ ! -s "$tr1" || ! -s "$tr2" ]]; then
        trim_galore --paired --cores "$THREADS" --gzip --fastqc \
            --nextera --output_dir "$trimg_dir" "$r1" "$r2"
    else
        skip "TrimGalore exists: ${samplename}"
    fi

    # Bowtie2 (add read groups)
    if [[ ! -s "$sam" ]]; then
        bowtie2 --very-sensitive -p "$THREADS" -X 2000 -x "$bowtie_index" \
            --rg-id "$samplename" \
            --rg "SM:$samplename" \
            --rg "LB:lib1" \
            --rg "PL:ILLUMINA" \
            --rg "PU:$samplename" \
            -1 "$tr1" -2 "$tr2" -S "$sam" 2> "${bowtie2_dir}/${samplename}.log"
    else
        skip "SAM exists: ${samplename}"
    fi

    # Aligned sorted BAM (pre-filter)
    if [[ ! -s "$aligned_bam" || ! -s "${aligned_bam}.bai" ]]; then
        samtools view -@ "$THREADS" -bS "$sam" | samtools sort -@ "$THREADS" -o "$aligned_bam" -
        samtools index -@ "$THREADS" "$aligned_bam"
        if [[ "$KEEP_INTERMEDIATES" -eq 0 ]]; then
            # Verify BAM before deleting SAM
            if samtools quickcheck "$aligned_bam"; then
                rm -f "$sam"
            else
                log "[ERROR] BAM failed quickcheck, keeping SAM: $aligned_bam"
            fi
        fi
    else
        skip "Aligned BAM exists: ${samplename}"
    fi

    # Read filters: MAPQ>=MAPQ, proper pairs, no secondary/supplementary, canonical chr only
    if [[ ! -s "$filt_bam" || ! -s "${filt_bam}.bai" ]]; then
        samtools view -@ "$THREADS" -h -q "$MAPQ" -f 2 -F 256 -F 2048 "$aligned_bam" \
            | awk "$awk_keep_canonical" \
            | samtools sort -@ "$THREADS" -o "$filt_bam" -
        samtools index -@ "$THREADS" "$filt_bam"
    else
        skip "Filtered BAM exists: ${samplename}"
    fi

    # Picard MarkDuplicates
    dedup_bam="${picard_dir}/${samplename}.dedup.bam"
    pic_metrics="${picard_dir}/${samplename}.metrics.txt"
    if [[ ! -s "$dedup_bam" || ! -s "$pic_metrics" ]]; then
        if have picard; then
            picard MarkDuplicates \
                I="$filt_bam" O="$dedup_bam" \
                M="$pic_metrics" \
                REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
        else
            java -jar "$PICARD" MarkDuplicates \
                I="$filt_bam" O="$dedup_bam" \
                M="$pic_metrics" \
                REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
        fi
    else
        skip "Dedup BAM exists: ${samplename}"
    fi

    # Final sort + index
    if [[ ! -s "$final_bam" || ! -s "${final_bam}.bai" ]]; then
        samtools sort -@ "$THREADS" -o "$final_bam" "$dedup_bam"
        samtools index -@ "$THREADS" "$final_bam"
        if [[ "$KEEP_INTERMEDIATES" -eq 0 ]]; then
            # Verify final BAM before cleanup
            if samtools quickcheck "$final_bam"; then
                rm -f "$dedup_bam"
            else
                log "[ERROR] Final BAM failed quickcheck: $final_bam"
            fi
        fi
    else
        skip "Final BAM exists: ${samplename}"
    fi

    # BigWig via deepTools (CPM, paired-end fragment coverage)
    bw="${deeptools_dir}/${samplename}.CPM.bw"
    if [[ ! -s "$bw" ]]; then
        bamCoverage --bam "$final_bam" \
            --outFileName "$bw" \
            --outFileFormat bigwig \
            --binSize "${ATAC_BW_BINSIZE:-10}" \
            --normalizeUsing CPM \
            --extendReads \
            --centerReads \
            --numberOfProcessors "$THREADS"
    else
        skip "BigWig exists: ${samplename}"
    fi

    # QC row
    if ! awk -v s="$samplename" 'BEGIN{FS=OFS="\t"} NR==1{next} $1==s{found=1} END{exit(found?0:1)}' "$qc_tsv"; then

        condition=$(get_sample_condition "$samplename")
        replicate=$(get_sample_replicate "$samplename")

        aln_rate=$(grep -m1 "overall alignment rate" "${bowtie2_dir}/${samplename}.log" | awk '{print $1}' || echo "NA")

        aligned_mapped=$(samtools idxstats "$aligned_bam" | awk '{sum+=$3} END{print sum+0}')
        aligned_chrM=$(samtools idxstats "$aligned_bam" | awk '$1=="chrM"{print $3+0}')
        aligned_chrM_frac=$(awk -v m="$aligned_mapped" -v c="$aligned_chrM" 'BEGIN{ if(m>0) printf("%.6f", c/m); else printf("NA") }')

        filtered_mapped=$(samtools idxstats "$filt_bam" | awk '{sum+=$3} END{print sum+0}')
        final_mapped=$(samtools idxstats "$final_bam" | awk '{sum+=$3} END{print sum+0}')

        # Robust Picard metrics parsing using header
        pic_pct_dup=$(awk -F'\t' '
            /^LIBRARY/ { for(i=1;i<=NF;i++) if($i=="PERCENT_DUPLICATION") col=i; next }
            col && NF>=col && !/^#/ && !/^$/ { print $col; exit }
        ' "$pic_metrics" 2>/dev/null || echo "NA")

        pic_lib_size=$(awk -F'\t' '
            /^LIBRARY/ { for(i=1;i<=NF;i++) if($i=="ESTIMATED_LIBRARY_SIZE") col=i; next }
            col && NF>=col && !/^#/ && !/^$/ { print $col; exit }
        ' "$pic_metrics" 2>/dev/null || echo "NA")

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$samplename" "$condition" "$replicate" \
            "$aln_rate" "$aligned_mapped" "$aligned_chrM" "$aligned_chrM_frac" \
            "$filtered_mapped" "$final_mapped" \
            "$pic_pct_dup" "$pic_lib_size" >> "$qc_tsv"

    else
        skip "QC already recorded: ${samplename}"
    fi

done

[[ $found_any -eq 1 ]] || { log "No *_R1_001.fastq.gz found in ${raw_atac}"; exit 1; }

# MultiQC
log "STEP 1b: MultiQC"
if [[ ! -s "${multiqc_dir}/multiqc_report.html" ]]; then
    multiqc -f -o "$multiqc_dir" "$analysis_dir" 2>&1 | tee "${log_dir}/multiqc.log" || true
else
    skip "MultiQC exists: ${multiqc_dir}/multiqc_report.html"
fi

# =============================================================================
# STEP 2: PER-REPLICATE PEAKS (MACS3) + BLACKLIST FILTER
# =============================================================================
log "STEP 2: Per-replicate MACS3 peaks + blacklist filter"
log "[INFO] Using MACS: $MACS_BIN"
"$MACS_BIN" --version

mkdir -p "${macs_dir}/rep_peaks" "${macs_dir}/rep_peaks_bl"

for fbam in "${samtools_dir}"/*.final.bam; do
    [[ -f "$fbam" ]] || continue
    tag=$(basename "$fbam" .final.bam)

    out_np="${macs_dir}/rep_peaks/${tag}_peaks.narrowPeak"
    if [[ ! -s "$out_np" ]]; then
        "$MACS_BIN" callpeak -t "$fbam" -n "$tag" --outdir "${macs_dir}/rep_peaks" \
            -f BAMPE -g "$GENOME_SIZE_MACS" --nomodel --keep-dup all --slocal 1000 \
            2> "${macs_dir}/rep_peaks/${tag}.log"
    else
        skip "MACS per-rep exists: $tag"
    fi

    bl_np="${macs_dir}/rep_peaks_bl/${tag}.bl.narrowPeak"
    if [[ ! -s "$bl_np" ]]; then
        bedtools intersect -a "$out_np" -b "$bl_regions" -v > "$bl_np"
    else
        skip "BL-filtered peaks exist: $tag"
    fi
done

# =============================================================================
# STEP 3: UNION PEAKS -> SAF
# =============================================================================
log "STEP 3: Union of BL-filtered replicate peaks"
union_bed="${macs_dir}/union_peaks.bed"
union_saf="${macs_dir}/union_peaks.saf"

if [[ ! -s "$union_bed" ]]; then
    cat "${macs_dir}/rep_peaks_bl/"*.bl.narrowPeak \
        | cut -f1-3 \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"peak_"NR}' > "$union_bed"
else
    skip "Union BED exists"
fi

# SAF format: GeneID, Chr, Start (1-based), End, Strand
# BED is 0-based half-open, SAF is 1-based inclusive
# CRITICAL FIX: Convert start from 0-based to 1-based
if [[ ! -s "$union_saf" ]]; then
    awk 'BEGIN{OFS="\t"}{print $4,$1,$2+1,$3,"."}' "$union_bed" > "$union_saf"
else
    skip "Union SAF exists"
fi

# =============================================================================
# STEP 3.5: HOMER: PEAK ANNOTATION + MOTIF ENRICHMENT
# =============================================================================
log "STEP 3.5: HOMER annotation + motifs (size=${ATAC_HOMER_SIZE:-200})"

# Load HOMER via module (Minerva)
if have ml; then
    ml homer/4.10 2>/dev/null || ml homer 2>/dev/null || true
elif have module; then
    module load homer/4.10 2>/dev/null || module load homer 2>/dev/null || true
fi

# Project-local preparsed dir for speed/reproducibility
preparsedDir="${homer_dir}/preparsed_dir_mm10"
mkdir -p "$homer_dir" "$preparsedDir"

# Center peaks for consistent -size window
union_center_bed="${macs_dir}/union_peaks.center1bp.bed"
if [[ ! -s "$union_center_bed" ]]; then
    awk 'BEGIN{OFS="\t"}{
        mid=int(($2+$3)/2);
        s=mid; e=mid+1;
        if(s<0) s=0;
        print $1,s,e,$4
    }' "$union_bed" > "$union_center_bed"
else
    skip "Union center BED exists"
fi

homer_size="${ATAC_HOMER_SIZE:-200}"
homer_anno="${homer_dir}/union_peaks.annotatePeaks.size${homer_size}.txt"
homer_motif_dir="${homer_dir}/motifs_union.size${homer_size}"
mkdir -p "$homer_motif_dir"

# Peak annotation
if [[ ! -s "$homer_anno" ]]; then
    annotatePeaks.pl "$union_center_bed" mm10 -size "$homer_size" > "$homer_anno"
else
    skip "HOMER annotation exists: $homer_anno"
fi

# Motif enrichment
if [[ ! -s "${homer_motif_dir}/homerResults.html" ]]; then
    findMotifsGenome.pl "$union_center_bed" mm10 "$homer_motif_dir" \
        -preparsedDir "$preparsedDir" \
        -size "$homer_size" \
        -len "${ATAC_HOMER_MOTIF_LEN:-8,10,12}" \
        -p "$THREADS" \
        -mask
else
    skip "HOMER motifs exist: ${homer_motif_dir}"
fi

# Tag directories
tag_root="${homer_dir}/tagdirs"
mkdir -p "$tag_root"
for fbam in "${samtools_dir}"/*.final.bam; do
    [[ -f "$fbam" ]] || continue
    tag=$(basename "$fbam" .final.bam)
    td="${tag_root}/${tag}"
    if [[ ! -d "$td" ]]; then
        makeTagDirectory "$td" "$fbam" -sspe -keepOne -tbp 1
    else
        skip "TagDir exists: $td"
    fi
done

if command -v module >/dev/null 2>&1; then
    module unload homer/4.10 2>/dev/null || module unload homer 2>/dev/null || true
fi

# =============================================================================
# STEP 4: featureCounts (UNION PEAKS)
# =============================================================================
log "STEP 4: featureCounts (union peaks)"
fc_out="${subread_dir}/union_peaks_featureCounts.txt"
if [[ ! -s "$fc_out" ]]; then
    featureCounts -p -B -C -T "$THREADS" \
        -a "$union_saf" -F SAF \
        -o "$fc_out" \
        "${samtools_dir}/"*.final.bam
else
    skip "featureCounts exists: $fc_out"
fi

# =============================================================================
# STEP 5: DESeq2
# =============================================================================
log "STEP 5: DESeq2"
if [[ ! -s "${deseq_dir}/DA_results_DESeq2.csv" ]]; then
    if have ml; then
        ml R/4.2.0 2>/dev/null || true
    elif have module; then
        module load R/4.2.0 2>/dev/null || true
    fi

    have Rscript || { log "[ERROR] Rscript not found in PATH"; exit 1; }
    [[ -s "$deseq_R" ]] || { log "[ERROR] Missing DESeq2 script: $deseq_R"; exit 1; }

    # Pass sample metadata file as third argument if it exists
    if [[ -s "$sample_meta" ]]; then
        Rscript "$deseq_R" "$fc_out" "$deseq_dir" "$sample_meta" 2> "${log_dir}/deseq2.stderr.log"
    else
        Rscript "$deseq_R" "$fc_out" "$deseq_dir" 2> "${log_dir}/deseq2.stderr.log"
    fi

    if have ml; then
        ml unload R/4.2.0 2>/dev/null || true
    elif have module; then
        module unload R/4.2.0 2>/dev/null || true
    fi
else
    skip "DESeq2 outputs exist: ${deseq_dir}"
fi

log "[DONE] Key outputs:"
log "  Versions:       ${versions_txt}"
log "  MultiQC:        ${multiqc_dir}/multiqc_report.html"
log "  QC table:       ${qc_tsv}"
log "  Final BAMs:     ${samtools_dir}/*.final.bam"
log "  BigWigs (CPM):  ${deeptools_dir}/*.CPM.bw"
log "  Per-rep peaks:  ${macs_dir}/rep_peaks/*_peaks.narrowPeak"
log "  Union peaks:    ${macs_dir}/union_peaks.bed"
log "  HOMER anno:     ${homer_dir}/union_peaks.annotatePeaks.size${homer_size}.txt"
log "  HOMER motifs:   ${homer_dir}/motifs_union.size${homer_size}/homerResults.html"
log "  Tag dirs:       ${homer_dir}/tagdirs/*"
log "  Count matrix:   ${subread_dir}/union_peaks_featureCounts.txt"
log "  DESeq2:         ${deseq_dir}/DA_results_DESeq2.csv + PCA_samples.pdf + PCA_coordinates.csv"
