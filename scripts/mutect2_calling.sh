#!/bin/bash

#Author: Babul Pradhan, Ph.D. Student, Institute of Bioinformatics (IOB), Bangalore, India
#Date: 09/02/2025
#Purpose: Somatic variant calling from preprocessed tumor RNA-seq BAM files using GATK Mutect2 in tumor-only mode with Panel of Normals (PoN) and gnomAD population allele frequency database

# scripts/mutect2_calling.sh
# Mutect2 somatic calling: tumor-normal paired OR tumor-only fallback.
# Tumor-only mode is used automatically when --normal-bam / --normal-sm are absent,
# or when --tumor-only is explicitly passed.
#
# Usage:
#   Paired:      mutect2_calling.sh [opts] TUMOR_BAM NORMAL_BAM OUT_DIR
#   Tumor-only:  mutect2_calling.sh --tumor-only [opts] TUMOR_BAM OUT_DIR

set -euo pipefail
log() { echo "[$(date '+%F %T')] $*"; }
err() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

GATK=""; REF=""; PON=""; GNOMAD=""
TUMOR_SM=""; NORMAL_SM=""
CLEANUP=false
TUMOR_ONLY=false

while [[ "$#" -gt 0 && "$1" =~ ^- ]]; do
  case "$1" in
    --gatk)        GATK="$2";      shift 2 ;;
    --ref)         REF="$2";       shift 2 ;;
    --pon)         PON="$2";       shift 2 ;;
    --gnomad)      GNOMAD="$2";    shift 2 ;;
    --tumor-sm)    TUMOR_SM="$2";  shift 2 ;;
    --normal-sm)   NORMAL_SM="$2"; shift 2 ;;
    --tumor-only)  TUMOR_ONLY=true; shift ;;
    --cleanup|-c)  CLEANUP=true;   shift ;;
    -h|--help)
      echo "Usage: mutect2_calling.sh [--tumor-only] [opts] TUMOR_BAM [NORMAL_BAM] OUT_DIR"
      exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

# Positional args: paired = 3, tumor-only = 2
if $TUMOR_ONLY || [[ $# -eq 2 ]]; then
  [[ $# -eq 2 ]] || err "Tumor-only mode requires: TUMOR_BAM OUT_DIR"
  TUMOR_BAM="$1"; OUT_DIR="$2"
  TUMOR_ONLY=true
elif [[ $# -eq 3 ]]; then
  TUMOR_BAM="$1"; NORMAL_BAM="$2"; OUT_DIR="$3"
else
  err "Usage: mutect2_calling.sh [--tumor-only] TUMOR_BAM [NORMAL_BAM] OUT_DIR"
fi

[[ -n "$GATK" && -n "$REF" ]] || err "--gatk and --ref are required."
[[ -f "$TUMOR_BAM" ]] || err "Tumor BAM not found: $TUMOR_BAM"
! $TUMOR_ONLY && [[ -f "$NORMAL_BAM" ]] || { $TUMOR_ONLY || err "Normal BAM not found: $NORMAL_BAM"; }
mkdir -p "$OUT_DIR"

# Derive sample names from BAM filenames if not provided
[[ -n "$TUMOR_SM" ]] || TUMOR_SM=$(basename "$TUMOR_BAM" | sed -E 's/(_recal)?(Aligned\.sortedByCoord\.out)?\.bam$//')
if ! $TUMOR_ONLY; then
  [[ -n "$NORMAL_SM" ]] || NORMAL_SM=$(basename "$NORMAL_BAM" | sed -E 's/(_recal)?(Aligned\.sortedByCoord\.out)?\.bam$//')
fi

UNF="$OUT_DIR/unfiltered.vcf.gz"
FIL="$OUT_DIR/mutect2.vcf.gz"

# Idempotency check
if [[ -f "$FIL" && -f "${FIL}.tbi" ]]; then
  log "Skip Mutect2 (filtered VCF exists): $FIL"
  exit 0
fi

# Build Mutect2 command
CMD=( "$GATK" Mutect2 -R "$REF" -I "$TUMOR_BAM" -tumor "$TUMOR_SM" )

if $TUMOR_ONLY; then
  log "Mutect2 mode: TUMOR-ONLY | sample=$TUMOR_SM"
else
  CMD+=( -I "$NORMAL_BAM" -normal "$NORMAL_SM" )
  log "Mutect2 mode: PAIRED | tumor=$TUMOR_SM | normal=$NORMAL_SM"
fi

[[ -n "$PON"    && -f "$PON"    ]] && CMD+=( --panel-of-normals "$PON" )
[[ -n "$GNOMAD" && -f "$GNOMAD" ]] && CMD+=( --germline-resource "$GNOMAD" )
CMD+=( -O "$UNF" )

log "Running: ${CMD[*]}"
"${CMD[@]}"

log "FilterMutectCalls -> $FIL"
"$GATK" FilterMutectCalls -R "$REF" -V "$UNF" -O "$FIL"

$CLEANUP && { log "Cleanup: removing $UNF"; rm -f "$UNF" "${UNF}.tbi" || true; }
log "Done: $FIL"
