#!/bin/bash
# scripts/build_pon.sh — per-normal Mutect2 tumor-only -> GenomicsDBImport -> CreateSomaticPanelOfNormals

set -euo pipefail
log() { echo "[$(date '+%F %T')] $*" >&2; }
err() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

GATK=""; REF=""; JOBS=1; MODE="missing"; TMPDIR=""; WORKSPACE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gatk)      GATK="$2"; shift 2 ;;
    --ref)       REF="$2"; shift 2 ;;
    --jobs)      JOBS="$2"; shift 2 ;;
    --mode)      MODE="$2"; shift 2 ;;
    --tmp-dir)   TMPDIR="$2"; shift 2 ;;
    --workspace) WORKSPACE="$2"; shift 2 ;;
    -h|--help)   echo "Usage: build_pon.sh --gatk P --ref P [--jobs N] [--mode missing|all|ask] [--workspace DIR] PON_LIST OUTPUT_VCF"; exit 0 ;;
    *) break ;;
  esac
done
[[ -n "$GATK" && -n "$REF" && $# -eq 2 ]] || err "Usage: build_pon.sh --gatk P --ref P PON_LIST OUTPUT_VCF"

PON_LIST="$1"; OUT_VCF="$2"
OUT_DIR="$(dirname "$(readlink -f "$OUT_VCF")")"
[[ -n "$TMPDIR" ]] || TMPDIR="$OUT_DIR/pon_work"
[[ -n "$WORKSPACE" ]] || WORKSPACE="$OUT_DIR/pon_db"
LOG_DIR="$OUT_DIR/logs/pon"; mkdir -p "$TMPDIR" "$WORKSPACE" "$LOG_DIR"

# Early reuse
if [[ -f "$OUT_VCF" && -f "${OUT_VCF}.tbi" ]]; then
  case "$MODE" in
    all)     log "PoN exists — rebuilding (mode=all)" ;;
    missing) log "PoN exists — reusing (mode=missing): $OUT_VCF"; exit 0 ;;
    ask)
      if [[ -t 0 ]]; then
        read -p "Reuse existing PoN ($OUT_VCF)? (Y/n) " -r ans; ans=${ans:-Y}
        [[ "$ans" =~ ^[Yy]$ ]] && { log "Reusing PoN: $OUT_VCF"; exit 0; }
      else log "Reusing PoN (non-interactive): $OUT_VCF"; exit 0; fi ;;
  esac
fi

mapfile -t BAMS < <(grep -v '^[[:space:]]*$' "$PON_LIST")
(( ${#BAMS[@]} )) || err "PON list is empty: $PON_LIST"

SAMPLES=(); VCFS=(); LOGS=()
for bam in "${BAMS[@]}"; do
  bn="$(basename "$bam")"; sm="${bn%.bam}"; sm="${sm%_recal}"
  SAMPLES+=( "$sm" ); VCFS+=( "$TMPDIR/${sm}_pon.vcf.gz" ); LOGS+=( "$LOG_DIR/${sm}.log" )
done

MISSING_VCF=()
for vcf in "${VCFS[@]}"; do [[ -f "$vcf" && -f "${vcf}.tbi" ]] || MISSING_VCF+=( "$vcf" ); done
log "Per-normal VCFs: ${#VCFS[@]} total | missing=${#MISSING_VCF[@]} (mode=$MODE)"

supports_waitn() { local a=${BASH_VERSINFO[0]:-0} b=${BASH_VERSINFO[1]:-0}; (( a>4 || (a==4 && b>=3) )); }

CMDS=()
case "$MODE" in
  all)
    for i in "${!SAMPLES[@]}"; do
      CMDS+=( "$GATK Mutect2 -R \"$REF\" --max-mnp-distance 0 -I \"${BAMS[$i]}\" --tumor \"${SAMPLES[$i]}\" -O \"${VCFS[$i]}\" >\"${LOGS[$i]}\" 2>&1" )
    done ;;
  missing|ask)
    for i in "${!SAMPLES[@]}"; do
      [[ -f "${VCFS[$i]}" && -f "${VCFS[$i]}.tbi" ]] && continue
      CMDS+=( "$GATK Mutect2 -R \"$REF\" --max-mnp-distance 0 -I \"${BAMS[$i]}\" --tumor \"${SAMPLES[$i]}\" -O \"${VCFS[$i]}\" >\"${LOGS[$i]}\" 2>&1" )
    done ;;
esac

if (( ${#CMDS[@]} )); then
  log "Running ${#CMDS[@]} per-normal Mutect2 calls (JOBS=$JOBS)..."
  _pids=(); _use_waitn=false; supports_waitn && _use_waitn=true
  for cmd in "${CMDS[@]}"; do
    bash -c "$cmd" & _pids+=( $! )
    if (( ${#_pids[@]} >= JOBS )); then
      $_use_waitn && { wait -n || true; } || { wait "${_pids[0]}" || true; _pids=( "${_pids[@]:1}" ); }
    fi
  done
  for pid in "${_pids[@]}"; do wait "$pid" || true; done
fi

MISS=()
for vcf in "${VCFS[@]}"; do [[ -f "$vcf" && -f "${vcf}.tbi" ]] || MISS+=( "$vcf" ); done
(( ${#MISS[@]} == 0 )) || err "Missing per-normal VCFs: ${MISS[*]}"

# Intervals from dict
DICT="${REF%.*}.dict"; [[ -f "$DICT" ]] || err "Reference dict not found: $DICT"
INT_FILE="$TMPDIR/pon.intervals.list"; : > "$INT_FILE"
awk '/^@SQ/{for(i=1;i<=NF;i++){if($i~/^SN:/){split($i,a,":");print a[2]}}}' "$DICT" >> "$INT_FILE"
[[ -s "$INT_FILE" ]] || err "No contigs in dict: $DICT"

# GenomicsDBImport
if [[ -d "$WORKSPACE" && -s "$WORKSPACE/MANIFEST.json" ]] && [[ "$MODE" != "all" ]]; then
  log "Reusing GenomicsDB workspace: $WORKSPACE"
else
  log "GenomicsDBImport -> $WORKSPACE"
  args=( GenomicsDBImport -R "$REF" --genomicsdb-workspace-path "$WORKSPACE" --reader-threads "$JOBS" -L "$INT_FILE" )
  for vcf in "${VCFS[@]}"; do args+=( -V "$vcf" ); done
  "$GATK" "${args[@]}" >"$LOG_DIR/genomicsdbimport.log" 2>&1
fi

# CreateSomaticPanelOfNormals
log "CreateSomaticPanelOfNormals -> $OUT_VCF"
"$GATK" CreateSomaticPanelOfNormals -R "$REF" -V "gendb://$WORKSPACE" -O "$OUT_VCF" >"$LOG_DIR/create_pon.log" 2>&1
[[ -f "$OUT_VCF" ]] || err "PoN creation failed — see $LOG_DIR/create_pon.log"
log "PoN complete: $OUT_VCF"
