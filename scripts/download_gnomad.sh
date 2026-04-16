#!/bin/bash
# scripts/download_gnomad.sh — Download and merge gnomAD v4.1 per-chromosome VCFs.
# gnomAD v4.1 ships with chr-prefix (UCSC); handles conversion to nochr (Ensembl).
# Prints final merged VCF path to stdout (captured by install.sh).
# Requires: gsutil, bcftools, tabix

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="$(dirname "$SCRIPT_DIR")/install.log"
log()  { echo "[$(date '+%F %T')] $*" | tee -a "$LOG" >&2; }
err()  { echo "[$(date '+%F %T')] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

OUT_DIR="${OUT_DIR:-$HOME/references/GRCh38}"
THREADS=8
CHR_STYLE="nochr"
GNOMAD_VERSION="4.1"
GNOMAD_GCS="gs://gcp-public-data--gnomad/release/${GNOMAD_VERSION}/vcf/genomes"
CHROMS=( {1..22} X Y )

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir)     OUT_DIR="$2";   shift 2 ;;
    --threads)     THREADS="$2";   shift 2 ;;
    --chr-style)   CHR_STYLE="$2"; shift 2 ;;
    --version)
      GNOMAD_VERSION="$2"
      GNOMAD_GCS="gs://gcp-public-data--gnomad/release/${GNOMAD_VERSION}/vcf/genomes"
      shift 2 ;;
    *) err "Unknown option: $1" ;;
  esac
done

[[ "$CHR_STYLE" == "nochr" || "$CHR_STYLE" == "chr" ]] \
  || err "--chr-style must be 'nochr' or 'chr'"

command -v gsutil   &>/dev/null || err "gsutil not found. Install Google Cloud SDK."
command -v bcftools &>/dev/null || err "bcftools not found."
command -v tabix    &>/dev/null || err "tabix not found."

GNOMAD_DIR="$OUT_DIR/gnomad_v${GNOMAD_VERSION}_per_chrom"
mkdir -p "$GNOMAD_DIR"

MERGED_VCF="$OUT_DIR/gnomad.genomes.v${GNOMAD_VERSION}.sites.merged.${CHR_STYLE}.vcf.gz"

if [[ -f "$MERGED_VCF" && -f "${MERGED_VCF}.tbi" ]]; then
  log "gnomAD v${GNOMAD_VERSION} (${CHR_STYLE}) already merged: $MERGED_VCF"
  echo "$MERGED_VCF"; exit 0
fi

# Build chr->nochr rename map (gnomAD ships chr-prefix)
RENAME_MAP="$OUT_DIR/gnomad_chr_to_nochr.map"
if [[ "$CHR_STYLE" == "nochr" ]]; then
  : > "$RENAME_MAP"
  for c in {1..22} X Y; do
    printf "chr%s\t%s\n" "$c" "$c" >> "$RENAME_MAP"
  done
  printf "chrM\tMT\n"  >> "$RENAME_MAP"
  printf "chrMT\tMT\n" >> "$RENAME_MAP"
fi

VCF_LIST=()
for C in "${CHROMS[@]}"; do
  GCS_VCF="${GNOMAD_GCS}/gnomad.genomes.v${GNOMAD_VERSION}.sites.chr${C}.vcf.bgz"
  LOCAL_RAW="$GNOMAD_DIR/gnomad.v${GNOMAD_VERSION}.chr${C}.raw.vcf.bgz"

  if [[ ! -f "$LOCAL_RAW" ]]; then
    log "gnomAD chr${C}: downloading..."
    gsutil -m cp "$GCS_VCF"      "$LOCAL_RAW"       2>>"$LOG" || err "gsutil failed for chr${C}"
    gsutil -m cp "${GCS_VCF}.tbi" "${LOCAL_RAW}.tbi" 2>>"$LOG" \
      || { log "Rebuilding tbi for chr${C}..."; tabix -p vcf "$LOCAL_RAW" 2>>"$LOG"; }
  fi

  if [[ "$CHR_STYLE" == "nochr" ]]; then
    LOCAL_OUT="$GNOMAD_DIR/gnomad.v${GNOMAD_VERSION}.chr${C}.nochr.vcf.gz"
    if [[ ! -f "$LOCAL_OUT" ]]; then
      log "gnomAD chr${C}: converting chr -> nochr..."
      bcftools annotate \
        --rename-chrs "$RENAME_MAP" \
        --threads 2 -Oz -o "$LOCAL_OUT" \
        "$LOCAL_RAW" 2>>"$LOG" || err "bcftools rename failed for chr${C}"
      tabix -p vcf "$LOCAL_OUT" 2>>"$LOG"
    fi
  else
    LOCAL_OUT="$GNOMAD_DIR/gnomad.v${GNOMAD_VERSION}.chr${C}.chr.vcf.gz"
    if [[ ! -f "$LOCAL_OUT" ]]; then
      bcftools view --threads 2 -Oz -o "$LOCAL_OUT" "$LOCAL_RAW" 2>>"$LOG" \
        || err "bcftools view failed for chr${C}"
      tabix -p vcf "$LOCAL_OUT" 2>>"$LOG"
    fi
  fi
  VCF_LIST+=( "$LOCAL_OUT" )
done

log "Merging ${#VCF_LIST[@]} gnomAD VCFs -> $MERGED_VCF ..."
printf '%s\n' "${VCF_LIST[@]}" > "$GNOMAD_DIR/concat_list.txt"
bcftools concat \
  --file-list "$GNOMAD_DIR/concat_list.txt" \
  --allow-overlaps \
  --threads "$THREADS" \
  -Oz -o "$MERGED_VCF" 2>>"$LOG" || err "bcftools concat failed"
tabix -p vcf "$MERGED_VCF" 2>>"$LOG"
[[ -f "$RENAME_MAP" ]] && rm -f "$RENAME_MAP"
log "gnomAD merged: $MERGED_VCF"
echo "$MERGED_VCF"
