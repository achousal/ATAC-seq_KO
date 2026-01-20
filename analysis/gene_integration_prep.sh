#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ---------------- helpers ----------------
log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"; }

# ---------------- config (edit if needed) ----------------
BASE="/sc/arion/projects/Chipuk_Laboratory/chousa01/ATAC-seq_KO/analysis/atacseq_analysis.sh"
MACS="${BASE}/analysis/macs"
DESEQ="${BASE}/analysis/subread/deseq/DA_results_DESeq2.csv"
GTF="${BASE}/raw/gencode.vM25.annotation.gtf"

# output folder
OUT="${BASE}/analysis/gene_integration"

# windows
PROM_WIN=2000      # promoter: within +/- 2kb of TSS
ENH_MAX=50000      # enhancer: 2kb–50kb from TSS

# If GTF uses "1" but peaks use "chr1", set to 1 (or run with --add_chr_prefix)
ADD_CHR_PREFIX=0

# ---------------- optional CLI overrides ----------------
usage(){
  cat <<EOF
Usage:
  $(basename "$0") [options]

Options (override defaults in script):
  --base PATH
  --macs PATH
  --deseq CSV
  --gtf  GTF
  --outdir DIR
  --prom_win INT
  --enh_max INT
  --add_chr_prefix
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base) BASE="$2"; MACS="${BASE}/analysis/macs"; DESEQ="${BASE}/analysis/subread/deseq/DA_results_DESeq2.csv"; GTF="${BASE}/raw/gencode.vM25.annotation.gtf"; shift 2;;
    --macs) MACS="$2"; shift 2;;
    --deseq) DESEQ="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --prom_win) PROM_WIN="$2"; shift 2;;
    --enh_max) ENH_MAX="$2"; shift 2;;
    --add_chr_prefix) ADD_CHR_PREFIX=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

# ---------------- environment/tools ----------------
# Load tools if modules exist
if command -v module >/dev/null 2>&1; then
  module load bedtools >/dev/null 2>&1 || true
fi

need awk
need sort
need bedtools

# ---------------- inputs ----------------
UNION_PEAKS_BED="${MACS}/union_peaks.bed"

[[ -f "$GTF" ]] || die "Missing GTF: $GTF"
[[ -f "$UNION_PEAKS_BED" ]] || die "Missing union peaks BED: $UNION_PEAKS_BED"
[[ -f "$DESEQ" ]] || die "Missing DESeq2 peak DA CSV: $DESEQ"

mkdir -p "$OUT"
log "OUT=$OUT"
log "UNION_PEAKS_BED=$UNION_PEAKS_BED"
log "DESEQ=$DESEQ"
log "GTF=$GTF"
log "PROM_WIN=$PROM_WIN ENH_MAX=$ENH_MAX ADD_CHR_PREFIX=$ADD_CHR_PREFIX"

# -----------------------------
# 1) GTF -> TSS BED (1bp)
# columns: chr start0 end1 gene_id gene_name strand
# -----------------------------
log "[1/5] Build TSS BED -> $OUT/tss.bed"
awk -F'\t' -v addchr="$ADD_CHR_PREFIX" 'BEGIN{OFS="\t"}
  $0 !~ /^#/ && $3=="gene" {
    chr=$1; start=$4; end=$5; strand=$7; attr=$9;

    if (addchr==1 && chr !~ /^chr/) chr="chr"chr;

    gene_id="NA"; gene_name="NA";
    if (match(attr, /gene_id "([^"]+)"/, a)) gene_id=a[1];
    if (match(attr, /gene_name "([^"]+)"/, b)) gene_name=b[1];
    if (gene_name=="NA") gene_name=gene_id;

    tss = (strand=="+") ? start : end;  # 1-based
    tss0 = tss - 1;                     # BED start 0-based

    print chr, tss0, tss0+1, gene_id, gene_name, strand
  }' "$GTF" \
| LC_ALL=C sort -k1,1 -k2,2n \
> "$OUT/tss.bed"

# -----------------------------
# 2) Union peaks -> centers (1bp)
# peaks.center.bed: chr center center+1 peak_id
# -----------------------------
log "[2/5] Build peak centers -> $OUT/peaks.center.bed"
awk 'BEGIN{OFS="\t"}
  {
    chr=$1; s=$2; e=$3; id=$4;
    c = int((s+e)/2);
    print chr, c, c+1, id
  }' "$UNION_PEAKS_BED" \
| LC_ALL=C sort -k1,1 -k2,2n \
> "$OUT/peaks.center.bed"

# -----------------------------
# 3) Closest mapping: center -> TSS
# union_peaks_closestTSS.tsv columns:
#   chr start end peak_id gene_id gene_name strand dist_bp abs_dist
# dist_bp = peak_center_start - tss_start0 (signed)
# -----------------------------
log "[3/5] bedtools closest -> $OUT/union_peaks_closestTSS.tsv"
bedtools closest -a "$OUT/peaks.center.bed" -b "$OUT/tss.bed" -t first \
| awk 'BEGIN{OFS="\t"}
  {
    # a: 1-4  chr aS aE peak_id
    # b: 5-10 chr bS bE gene_id gene_name strand
    chr=$1; aS=$2; aE=$3; peak=$4;
    gene_id=$8; gene_name=$9; strand=$10;
    dist = aS - $6;             # center - tss0
    absd = (dist<0)? -dist:dist;
    print chr, aS, aE, peak, gene_id, gene_name, strand, dist, absd
  }' \
> "$OUT/union_peaks_closestTSS.tsv"

# -----------------------------
# 4) Promoter/enhancer split
# promoter_closest.tsv columns:
#   peak_id gene_id gene_name strand dist_bp abs_dist
# -----------------------------
log "[4/5] Split promoter/enhancer"
awk -v W="$PROM_WIN" 'BEGIN{OFS="\t"}
  { if ($9 <= W) print $4,$5,$6,$7,$8,$9 }' \
  "$OUT/union_peaks_closestTSS.tsv" > "$OUT/promoter_closest.tsv"

awk -v W="$PROM_WIN" -v M="$ENH_MAX" 'BEGIN{OFS="\t"}
  { if ($9 > W && $9 <= M) print $4,$5,$6,$7,$8,$9 }' \
  "$OUT/union_peaks_closestTSS.tsv" > "$OUT/enhancer_closest.tsv"

cut -f1 "$OUT/promoter_closest.tsv" | LC_ALL=C sort -u > "$OUT/promoter_peak_ids.txt"
cut -f1 "$OUT/enhancer_closest.tsv" | LC_ALL=C sort -u > "$OUT/enhancer_peak_ids.txt"

log "Promoter peaks: $(wc -l < "$OUT/promoter_peak_ids.txt")"
log "Enhancer peaks: $(wc -l < "$OUT/enhancer_peak_ids.txt")"

# -----------------------------
# 5) Subset ATAC DA CSV into promoter/enhancer DA
# Robustly finds the 'peak_id' column in the CSV header.
# -----------------------------
log "[5/5] Subset DESeq2 DA -> promoter_DA.csv / enhancer_DA.csv"

subset_csv_by_ids(){
  local csv="$1"
  local ids="$2"
  local out="$3"

  awk -v FS=',' -v OFS=',' -v idfile="$ids" '
    BEGIN{
      while ((getline line < idfile) > 0) { ids[line]=1 }
      close(idfile)
    }
    NR==1{
      peak_col=0
      for (i=1;i<=NF;i++){
        if ($i=="peak_id") { peak_col=i; break }
      }
      if (peak_col==0){
        print "ERROR: peak_id column not found in header" > "/dev/stderr"
        exit 2
      }
      print $0
      next
    }
    {
      pid = $peak_col
      if (pid in ids) print $0
    }
  ' "$csv" > "$out"
}

subset_csv_by_ids "$DESEQ" "$OUT/promoter_peak_ids.txt" "$OUT/promoter_DA.csv"
subset_csv_by_ids "$DESEQ" "$OUT/enhancer_peak_ids.txt"  "$OUT/enhancer_DA.csv"

log "Wrote: $OUT/tss.bed"
log "Wrote: $OUT/union_peaks_closestTSS.tsv"
log "Wrote: $OUT/promoter_closest.tsv / $OUT/enhancer_closest.tsv"
log "Wrote: $OUT/promoter_DA.csv / $OUT/enhancer_DA.csv"
log "DONE."