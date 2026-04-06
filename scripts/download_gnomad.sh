#!/bin/bash
# scripts/download_gnomad.sh
# Downloads gnomAD v4.1 genome VCFs per chromosome via gsutil,
# merges with bcftools concat, indexes, and cleans up chr files.
# If merged file already exists, skips download entirely.
#
# Usage: download_gnomad.sh --out-dir DIR [--threads INT] [--skip-if-exists]
# Output: <OUT_DIR>/gnomad.genomes.v4.1.sites.merged.vcf.gz (+ .tbi)

set -euo pipefail
log() { echo "[$(date '+%F %T')] [gnomAD] $*"; }
err() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

OUT_DIR=""
THREADS=8
CONCAT_THREADS=20
SKIP_IF_EXISTS=true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir)          OUT_DIR="$2"; shift 2 ;;
    --threads)          THREADS="$2"; CONCAT_THREADS="$2"; shift 2 ;;
    --no-skip)          SKIP_IF_EXISTS=false; shift ;;
    -h|--help)
      echo "Usage: download_gnomad.sh --out-dir DIR [--threads INT] [--no-skip]"
      exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

[[ -n "$OUT_DIR" ]] || err "--out-dir is required"
mkdir -p "$OUT_DIR"

MERGED="$OUT_DIR/gnomad.genomes.v4.1.sites.merged.vcf.gz"

# Skip entirely if merged file already present
if $SKIP_IF_EXISTS && [[ -f "$MERGED" && -f "${MERGED}.tbi" ]]; then
  log "Merged gnomAD already present — skipping download: $MERGED"
  echo "$MERGED"
  exit 0
fi

# Check gsutil
if ! command -v gsutil &>/dev/null; then
  err "gsutil not found. Install Google Cloud SDK: https://cloud.google.com/sdk/docs/install"
fi
command -v bcftools &>/dev/null || err "bcftools not found. Run install.sh first."

WORK_DIR="$OUT_DIR/gnomad_genomes_v4_1_tmp"
mkdir -p "$WORK_DIR"

log "Downloading gnomAD v4.1 per-chromosome VCFs to $WORK_DIR ..."
CHROMS=( {1..22} X Y )
for chr in "${CHROMS[@]}"; do
  BGZ="$WORK_DIR/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
  TBI="${BGZ}.tbi"
  if [[ -f "$BGZ" && -f "$TBI" ]]; then
    log "  chr${chr}: already present, skipping."
  else
    log "  chr${chr}: downloading..."
    gsutil -q cp \
      "gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz" \
      "$BGZ"
    gsutil -q cp \
      "gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz.tbi" \
      "$TBI"
  fi
done

# Build ordered VCF list
VCF_LIST="$WORK_DIR/vcf.list"; : > "$VCF_LIST"
ls "$WORK_DIR"/gnomad.genomes.v4.1.sites.chr*.vcf.bgz | sort -V >> "$VCF_LIST"
[[ -s "$VCF_LIST" ]] || err "No per-chromosome VCFs found in $WORK_DIR"

# Merge
log "Merging chromosomes -> $MERGED (concat threads=$CONCAT_THREADS)"
bcftools concat -Oz --threads "$CONCAT_THREADS" -f "$VCF_LIST" -o "$MERGED"

# Index
log "Indexing merged VCF (threads=$THREADS)"
bcftools index --tbi --threads "$THREADS" "$MERGED"

# Cleanup per-chromosome files
log "Cleaning up per-chromosome files..."
rm -f "$WORK_DIR"/gnomad.genomes.v4.1.sites.chr*.vcf.bgz \
       "$WORK_DIR"/gnomad.genomes.v4.1.sites.chr*.vcf.bgz.tbi \
       "$VCF_LIST"
rmdir "$WORK_DIR" 2>/dev/null || true

log "gnomAD ready: $MERGED"
echo "$MERGED"
