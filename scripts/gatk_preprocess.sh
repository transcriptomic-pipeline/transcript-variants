#!/bin/bash
# scripts/gatk_preprocess.sh — AddOrReplaceReadGroups, optional MarkDuplicates,
# SplitNCigarReads, BaseRecalibrator, ApplyBQSR

set -euo pipefail
log() { echo "[$(date '+%F %T')] [$SAMPLE] $*"; }
err() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

GATK=""; PICARD_JAR=""; REF=""; DBSNP=""; MILLS=""
SAMPLE=""; INPUT=""; OUT_DIR=""
MARK_DUP=true
RGID=""; RGLB="lib1"; RGPL="ILLUMINA"; RGPU=""; RGSM=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gatk)      GATK="$2"; shift 2 ;;
    --picard)    PICARD_JAR="$2"; shift 2 ;;
    --ref)       REF="$2"; shift 2 ;;
    --dbsnp)     DBSNP="$2"; shift 2 ;;
    --mills)     MILLS="$2"; shift 2 ;;
    --sample)    SAMPLE="$2"; shift 2 ;;
    --input)     INPUT="$2"; shift 2 ;;
    --out-dir)   OUT_DIR="$2"; shift 2 ;;
    --mark-dup)  MARK_DUP="$2"; shift 2 ;;
    --rgid)      RGID="$2"; shift 2 ;;
    --rglb)      RGLB="$2"; shift 2 ;;
    --rgpl)      RGPL="$2"; shift 2 ;;
    --rgpu)      RGPU="$2"; shift 2 ;;
    --rgsm)      RGSM="$2"; shift 2 ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

[[ -n "$GATK" && -n "$REF" && -n "$SAMPLE" && -n "$INPUT" && -n "$OUT_DIR" ]] \
  || err "Missing required args."

# Defaults for RG fields
[[ -z "$RGID" ]] && RGID="$SAMPLE"
[[ -z "$RGPU" ]] && RGPU="${SAMPLE}.unit1"
[[ -z "$RGSM" ]] && RGSM="$SAMPLE"

PICARD_CMD=""
if [[ -n "$PICARD_JAR" && -f "$PICARD_JAR" ]]; then
  PICARD_CMD="java -Xmx8G -jar $PICARD_JAR"
elif command -v picard &>/dev/null; then
  PICARD_CMD="picard"
else
  err "Picard not found. Provide --picard or install it."
fi

cd "$OUT_DIR"

# 1. AddOrReplaceReadGroups
log "AddOrReplaceReadGroups"
$PICARD_CMD AddOrReplaceReadGroups \
  -I "$INPUT" -O "${SAMPLE}_rg.bam" \
  -SO coordinate \
  -RGID "$RGID" -RGLB "$RGLB" -RGPL "$RGPL" -RGPU "$RGPU" -RGSM "$RGSM" \
  --CREATE_INDEX true

CURRENT_BAM="${OUT_DIR}/${SAMPLE}_rg.bam"

# 2. MarkDuplicates (optional)
if [[ "$MARK_DUP" == "true" || "$MARK_DUP" == "1" || "$MARK_DUP" == "yes" ]]; then
  log "MarkDuplicates"
  $PICARD_CMD MarkDuplicates \
    -I "$CURRENT_BAM" \
    -O "${SAMPLE}_dedup.bam" \
    -M "${SAMPLE}_dedup_metrics.txt" \
    --CREATE_INDEX true
  rm -f "$CURRENT_BAM" "${CURRENT_BAM%.bam}.bai"
  CURRENT_BAM="${OUT_DIR}/${SAMPLE}_dedup.bam"
else
  log "Skipping MarkDuplicates (--no-mark-dup)"
fi

# 3. SplitNCigarReads
log "SplitNCigarReads"
"$GATK" --java-options -Xmx8G SplitNCigarReads \
  -R "$REF" -I "$CURRENT_BAM" -O "${SAMPLE}_split.bam"
rm -f "$CURRENT_BAM" "${CURRENT_BAM%.bam}.bai"
CURRENT_BAM="${OUT_DIR}/${SAMPLE}_split.bam"

# 4. BaseRecalibrator
KNOWN_SITES=()
[[ -f "$DBSNP" ]] && KNOWN_SITES+=( --known-sites "$DBSNP" )
[[ -f "$MILLS" ]] && KNOWN_SITES+=( --known-sites "$MILLS" )

if (( ${#KNOWN_SITES[@]} )); then
  log "BaseRecalibrator"
  "$GATK" --java-options -Xmx8G BaseRecalibrator \
    -R "$REF" -I "$CURRENT_BAM" \
    "${KNOWN_SITES[@]}" \
    -O "${SAMPLE}_recal.table"

  # 5. ApplyBQSR
  log "ApplyBQSR"
  "$GATK" --java-options -Xmx8G ApplyBQSR \
    -R "$REF" -I "$CURRENT_BAM" \
    --bqsr-recal-file "${SAMPLE}_recal.table" \
    -O "${SAMPLE}_recal.bam"
  rm -f "$CURRENT_BAM" "${CURRENT_BAM%.bam}.bai"
else
  log "No known-sites VCFs found; skipping BQSR — renaming split BAM as recal BAM."
  mv "$CURRENT_BAM" "${OUT_DIR}/${SAMPLE}_recal.bam"
  [[ -f "${CURRENT_BAM%.bam}.bai" ]] && mv "${CURRENT_BAM%.bam}.bai" "${OUT_DIR}/${SAMPLE}_recal.bai"
fi

log "Preprocessing complete: ${OUT_DIR}/${SAMPLE}_recal.bam"
