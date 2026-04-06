#!/bin/bash
# run_transcript_variants.sh — Main entry point for transcript-variants somatic pipeline.
# Repo: https://github.com/transcriptomic-pipeline/transcript-variants
#
# Stages:
#   1. Preprocessing  (AddOrReplaceReadGroups, optional MarkDuplicates,
#                      SplitNCigarReads, BaseRecalibrator, ApplyBQSR)
#   2. PoN            (optional: custom path | build from matched normals | build from all normals)
#   3. Mutect2        (tumor-normal pairs, optional PoN + gnomAD germline resource)
#   4. Filter         (FilterMutectCalls -> PASS-only extraction)
#
# Tools are auto-installed via install.sh if not found (unless --no-auto-install).
# gnomAD is downloaded via scripts/download_gnomad.sh during install (default).
# If user provides --gnomad PATH, that path is used directly without any download.

set -euo pipefail
set -E
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/.pipeline_env"

log()  { echo "[$(date '+%F %T')] $*"; }
warn() { echo "[$(date '+%F %T')] WARN: $*" >&2; }
err()  { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }
trap 'err "Unexpected error at line $LINENO: $BASH_COMMAND"' ERR

# ============================================================
# DEFAULTS  (overridden by .pipeline_env then CLI)
# ============================================================
GATK_BIN=""; PICARD_JAR=""
DEFAULT_REF=""; DEFAULT_DBSNP=""; DEFAULT_MILLS=""; DEFAULT_GNOMAD=""

REF=""; DBSNP=""; MILLS=""; GNOMAD=""
OUT_DIR="$(pwd)/somatic_output"
JOBS=1
PREPROC_MODE="missing"       # missing | all | ask
MARK_DUPLICATES=true
USE_PON=false
PON_MODE=""                  # custom | build
PON_BUILD_NORMALS="all"      # all | matched
CUSTOM_PON=""
SAMPLE_MAP="$SCRIPT_DIR/config/sample_map.txt"
RG_FILE="$SCRIPT_DIR/config/samples_rg_details.txt"
TUMOR_PATTERNS="tumor,_T"
NORMAL_PATTERNS="normal,_N"
CLEANUP_PREPROC=false
CLEANUP_SOMATIC=false
AUTO_INSTALL=true
INSTALL_ONLY=false

print_usage() {
  cat <<'EOF'
Usage: run_transcript_variants.sh [OPTIONS]

  Repo: https://github.com/transcriptomic-pipeline/transcript-variants

Input:
  -m, --sample-map FILE    Sample map (default: config/sample_map.txt)
                           Columns: SAMPLE_ID  TYPE  BAM_PATH  [MATCHED_WITH]
  -rg, --read-group FILE   RG details file (default: config/samples_rg_details.txt)

Reference resources:
  -r,  --ref PATH          Reference FASTA (default: from .pipeline_env / install)
  --dbsnp PATH             dbSNP VCF (default: from .pipeline_env)
  --mills PATH             Mills indels VCF (default: from .pipeline_env)
  --gnomad PATH            gnomAD VCF — if provided, used directly (no download)
                           (default: ~/references/GRCh38/gnomad.genomes.v4.1.sites.merged.vcf.gz)

Output:
  -o, --out-dir DIR        Output directory (default: ./somatic_output)

Processing:
  --no-mark-dup            Skip MarkDuplicates (default: run)
  -j, --jobs INT           Parallel jobs (default: 1)
  --pre-proc-missing       Reuse existing outputs, run only missing (default)
  --pre-proc-all           Force rerun everything
  --pre-proc-ask           Ask before reusing (interactive)
  -cp, --cleanup-preproc   Remove preprocessing intermediates after completion
  -cs, --cleanup-somatic   Remove unfiltered.vcf.gz after filtering

PoN options:
  --use-pon                Enable PoN filtering (default: off)
  --pon-custom PATH        Use existing custom/external PoN (implies --use-pon)
  --pon-build              Build PoN from all normals (implies --use-pon)
  --pon-build-matched      Build PoN from matched/adjacent normals only

Tool overrides:
  --gatk PATH              GATK binary path
  --picard PATH            picard.jar path

Install:
  --no-auto-install        Do not auto-invoke install.sh
  --install-only           Run install.sh and exit

  -h, --help
EOF
}

# ============================================================
# Parse CLI
# ============================================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--sample-map)      SAMPLE_MAP="$2"; shift 2 ;;
    -rg|--read-group)     RG_FILE="$2"; shift 2 ;;
    -r|--ref)             REF="$2"; shift 2 ;;
    --dbsnp)              DBSNP="$2"; shift 2 ;;
    --mills)              MILLS="$2"; shift 2 ;;
    --gnomad)             GNOMAD="$2"; shift 2 ;;
    -o|--out-dir)         OUT_DIR="$2"; shift 2 ;;
    --no-mark-dup)        MARK_DUPLICATES=false; shift ;;
    -j|--jobs)            JOBS="$2"; shift 2 ;;
    --pre-proc-missing)   PREPROC_MODE="missing"; shift ;;
    --pre-proc-all)       PREPROC_MODE="all"; shift ;;
    --pre-proc-ask)       PREPROC_MODE="ask"; shift ;;
    -cp|--cleanup-preproc) CLEANUP_PREPROC=true; shift ;;
    -cs|--cleanup-somatic) CLEANUP_SOMATIC=true; shift ;;
    --use-pon)            USE_PON=true; shift ;;
    --pon-custom)         CUSTOM_PON="$2"; USE_PON=true; PON_MODE="custom"; shift 2 ;;
    --pon-build)          USE_PON=true; PON_MODE="build"; PON_BUILD_NORMALS="all"; shift ;;
    --pon-build-matched)  USE_PON=true; PON_MODE="build"; PON_BUILD_NORMALS="matched"; shift ;;
    --gatk)               GATK_BIN="$2"; shift 2 ;;
    --picard)             PICARD_JAR="$2"; shift 2 ;;
    --no-auto-install)    AUTO_INSTALL=false; shift ;;
    --install-only)       INSTALL_ONLY=true; shift ;;
    -h|--help)            print_usage; exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

# ============================================================
# Source environment from install
# ============================================================
[[ -f "$ENV_FILE" ]] && source "$ENV_FILE"

# Fill in defaults from .pipeline_env only if not set via CLI
[[ -z "$REF" ]]   && REF="${DEFAULT_REF:-}"
[[ -z "$DBSNP" ]] && DBSNP="${DEFAULT_DBSNP:-}"
[[ -z "$MILLS" ]] && MILLS="${DEFAULT_MILLS:-}"
[[ -z "$GNOMAD" ]] && GNOMAD="${DEFAULT_GNOMAD:-}"

# ============================================================
# Auto-install check
# ============================================================
needs_install() {
  # Check GATK
  local gbin="${GATK_BIN:-$(command -v gatk 2>/dev/null || true)}"
  [[ -x "$gbin" ]] || return 0
  # Check Picard
  local pjar="${PICARD_JAR:-}"
  [[ -f "$pjar" ]] || command -v picard &>/dev/null || return 0
  # Check samtools + bcftools
  command -v samtools &>/dev/null || return 0
  command -v bcftools &>/dev/null || return 0
  return 1
}

if $INSTALL_ONLY; then
  log "Running install.sh..."
  bash "$SCRIPT_DIR/install.sh"
  exit 0
fi

if $AUTO_INSTALL && needs_install; then
  log "Required tools not found — auto-invoking install.sh..."
  bash "$SCRIPT_DIR/install.sh"
  [[ -f "$ENV_FILE" ]] && source "$ENV_FILE"
  [[ -z "$REF" ]]    && REF="${DEFAULT_REF:-}"
  [[ -z "$DBSNP" ]]  && DBSNP="${DEFAULT_DBSNP:-}"
  [[ -z "$MILLS" ]]  && MILLS="${DEFAULT_MILLS:-}"
  [[ -z "$GNOMAD" ]] && GNOMAD="${DEFAULT_GNOMAD:-}"
fi

# Resolve GATK
[[ -z "$GATK_BIN" ]] && GATK_BIN="${DEFAULT_GATK_BIN:-$(command -v gatk 2>/dev/null || true)}"
[[ -x "$GATK_BIN" ]] || err "GATK not found. Run install.sh or provide --gatk PATH."

# Resolve Picard
[[ -z "$PICARD_JAR" ]] && PICARD_JAR="${DEFAULT_PICARD_JAR:-}"

# Reference validation
[[ -f "$REF" ]] || err "Reference FASTA not found: '${REF}'. Run install.sh or provide --ref PATH."
[[ -f "$DBSNP" ]] || warn "dbSNP not found at '${DBSNP}' — BQSR known-sites will be incomplete."
[[ -f "$MILLS" ]] || warn "Mills VCF not found at '${MILLS}' — BQSR known-sites will be incomplete."

# gnomAD: warn but don't abort; it is optional for Mutect2
if [[ -z "$GNOMAD" || ! -f "$GNOMAD" ]]; then
  GNOMAD_DEFAULT="$HOME/references/GRCh38/gnomad.genomes.v4.1.sites.merged.vcf.gz"
  if [[ -f "$GNOMAD_DEFAULT" ]]; then
    GNOMAD="$GNOMAD_DEFAULT"
    log "gnomAD auto-detected: $GNOMAD"
  else
    warn "gnomAD VCF not found. Mutect2 will run without germline resource (--gnomad not set)."
    warn "To install: bash install.sh (downloads gnomAD v4.1 by default)."
    warn "To skip:    run without --gnomad (already the case)."
    GNOMAD=""
  fi
fi

# ============================================================
# Setup output + logging
# ============================================================
mkdir -p "$OUT_DIR/logs"
LOG_DIR="$OUT_DIR/logs"
exec > >(tee -a "$LOG_DIR/master.stdout.log") 2>&1

log "==========================================="
log "  transcript-variants somatic pipeline"
log "  https://github.com/transcriptomic-pipeline/transcript-variants"
log "==========================================="
log "  GATK        : $GATK_BIN"
log "  Picard      : ${PICARD_JAR:-system}"
log "  Ref         : $REF"
log "  dbSNP       : ${DBSNP:-none}"
log "  Mills       : ${MILLS:-none}"
log "  gnomAD      : ${GNOMAD:-none}"
log "  MarkDup     : $MARK_DUPLICATES"
log "  Use PoN     : $USE_PON | mode=${PON_MODE:-n/a} | normals=${PON_BUILD_NORMALS}"
log "  Preproc     : $PREPROC_MODE"
log "  Jobs        : $JOBS"
log "  Out dir     : $OUT_DIR"

# ============================================================
# Load sample map
# ============================================================
declare -A SAMPLE_TYPE SAMPLE_BAM SAMPLE_MATCH
ALL_SAMPLES=()

classify_default() {
  local sm="$1" lo; lo=$(tr '[:upper:]' '[:lower:]' <<<"$sm")
  IFS=',' read -ra TP <<<"$TUMOR_PATTERNS"
  for p in "${TP[@]}"; do
    [[ "$lo" == *"$(tr '[:upper:]' '[:lower:]' <<<"$p")"* ]] && { echo tumor; return; }
  done
  IFS=',' read -ra NP <<<"$NORMAL_PATTERNS"
  for p in "${NP[@]}"; do
    [[ "$lo" == *"$(tr '[:upper:]' '[:lower:]' <<<"$p")"* ]] && { echo normal; return; }
  done
  echo unknown
}

[[ -f "$SAMPLE_MAP" ]] || err "Sample map not found: $SAMPLE_MAP — copy config/sample_map.txt and fill it in."

log "Loading sample map: $SAMPLE_MAP"
while IFS=$'\t' read -r sid stype sbam smatch; do
  [[ -z "$sid" || "$sid" =~ ^# ]] && continue
  sid="${sid// /}"
  [[ -z "$stype" || "$stype" =~ ^# ]] && stype=$(classify_default "$sid")
  stype=$(tr '[:upper:]' '[:lower:]' <<<"${stype%% *}")
  [[ -f "$sbam" ]] || { warn "BAM not found for $sid: '$sbam' — skipping."; continue; }
  SAMPLE_TYPE[$sid]="$stype"
  SAMPLE_BAM[$sid]="$sbam"
  [[ -n "${smatch:-}" && ! "${smatch}" =~ ^# ]] && SAMPLE_MATCH[$sid]="${smatch%% *}"
  ALL_SAMPLES+=( "$sid" )
done <"$SAMPLE_MAP"

(( ${#ALL_SAMPLES[@]} )) || err "No valid samples found in $SAMPLE_MAP."

TUMORS=(); NORMALS=()
for sm in "${ALL_SAMPLES[@]}"; do
  case "${SAMPLE_TYPE[$sm]}" in
    tumor)  TUMORS+=( "$sm" ) ;;
    normal) NORMALS+=( "$sm" ) ;;
    *)      warn "Could not classify '$sm' as tumor/normal — skipping." ;;
  esac
done
log "Tumors  (${#TUMORS[@]}): ${TUMORS[*]}"
log "Normals (${#NORMALS[@]}): ${NORMALS[*]}"

# ============================================================
# Load RG details
# ============================================================
declare -A RG_MAP  # SAMPLE -> "RGID|RGLB|RGPL|RGPU|RGSM"
if [[ -f "$RG_FILE" ]]; then
  log "Loading RG details: $RG_FILE"
  while IFS=$'\t' read -r sample rgid rglb rgpl rgpu rgsm; do
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue
    RG_MAP["$sample"]="${rgid}|${rglb}|${rgpl}|${rgpu}|${rgsm}"
  done <"$RG_FILE"
else
  warn "RG file not found: $RG_FILE — default RG values will be used."
fi

get_rg() {
  local sm="$1"
  [[ -n "${RG_MAP[$sm]:-}" ]] && echo "${RG_MAP[$sm]}" \
                               || echo "${sm}|lib1|ILLUMINA|${sm}.unit1|${sm}"
}

# ============================================================
# Parallel runner helpers
# ============================================================
supports_waitn() {
  local a=${BASH_VERSINFO[0]:-0} b=${BASH_VERSINFO[1]:-0}
  (( a>4 || (a==4 && b>=3) ))
}
_run_par_items() {
  local func=$1; shift
  local -a pids=(); local use_waitn=false; supports_waitn && use_waitn=true
  for item in "$@"; do
    ( $func "$item" ) & pids+=( $! )
    if (( ${#pids[@]} >= JOBS )); then
      $use_waitn && { wait -n || true; } \
                 || { wait "${pids[0]}" || true; pids=( "${pids[@]:1}" ); }
    fi
  done
  for pid in "${pids[@]}"; do wait "$pid" || true; done
}
_run_par_pairs() {
  local -a pids=(); local use_waitn=false; supports_waitn && use_waitn=true
  for pair in "$@"; do
    ( IFS=':' read -r t n <<<"$pair"; run_mutect2_pair "$t" "$n" ) & pids+=( $! )
    if (( ${#pids[@]} >= JOBS )); then
      $use_waitn && { wait -n || true; } \
                 || { wait "${pids[0]}" || true; pids=( "${pids[@]:1}" ); }
    fi
  done
  for pid in "${pids[@]}"; do wait "$pid" || true; done
}

# ============================================================
# STEP 1: Preprocessing
# ============================================================
preprocess_sample() {
  local sm="$1"
  local bam="${SAMPLE_BAM[$sm]}"
  local sdir="$OUT_DIR/$sm"; mkdir -p "$sdir"
  local out_recal="$sdir/${sm}_recal.bam"

  case "$PREPROC_MODE" in
    missing)
      [[ -f "$out_recal" ]] && { log "Skip preprocess (exists): $out_recal"; return; } ;;
    ask)
      if [[ -f "$out_recal" ]]; then
        if [[ -t 0 ]]; then
          read -p "Reuse existing preprocessed BAM for $sm? (Y/n) " -r ans; ans=${ans:-Y}
          [[ "$ans" =~ ^[Yy]$ ]] && { log "Reuse preprocess: $out_recal"; return; }
        else
          log "Reuse preprocess (non-interactive): $out_recal"; return
        fi
      fi ;;
    all) : ;;
  esac

  # Build RG flag string safely
  IFS='|' read -ra RGP <<<"$(get_rg "$sm")"
  log "Preprocessing $sm ..."
  bash "$SCRIPT_DIR/scripts/gatk_preprocess.sh" \
    --gatk      "$GATK_BIN" \
    --picard    "${PICARD_JAR:-}" \
    --ref       "$REF" \
    --dbsnp     "${DBSNP:-}" \
    --mills     "${MILLS:-}" \
    --mark-dup  "$MARK_DUPLICATES" \
    --sample    "$sm" \
    --input     "$bam" \
    --out-dir   "$sdir" \
    --rgid      "${RGP[0]}" \
    --rglb      "${RGP[1]}" \
    --rgpl      "${RGP[2]}" \
    --rgpu      "${RGP[3]}" \
    --rgsm      "${RGP[4]}" \
    >"$LOG_DIR/${sm}_preprocess.log" 2>&1 \
    || { err "Preprocessing FAILED for $sm — see $LOG_DIR/${sm}_preprocess.log"; }

  [[ -f "$out_recal" ]] || err "Preprocessing produced no recal BAM for $sm."
  log "Preprocessing OK: $out_recal"
}

log "--- STEP 1: Preprocessing (${#ALL_SAMPLES[@]} samples) ---"
_run_par_items preprocess_sample "${ALL_SAMPLES[@]}"

# ============================================================
# STEP 2: Panel of Normals (optional)
# ============================================================
PON_VCF=""
if $USE_PON; then

  if [[ "$PON_MODE" == "custom" ]]; then
    [[ -f "$CUSTOM_PON" ]] || err "Custom PoN not found: $CUSTOM_PON"
    PON_VCF="$CUSTOM_PON"
    log "Using custom PoN: $PON_VCF"

  elif [[ "$PON_MODE" == "build" ]]; then
    PON_OUT="$OUT_DIR/pon/custom_pon.vcf.gz"
    mkdir -p "$OUT_DIR/pon"

    # Resolve which normals to include
    declare -a PON_NORMALS=()
    if [[ "$PON_BUILD_NORMALS" == "matched" ]]; then
      log "PoN build: using matched normals only."
      for t in "${TUMORS[@]}"; do
        local_match="${SAMPLE_MATCH[$t]:-}"
        if [[ -z "$local_match" ]]; then
          # Auto-pair by stripping known suffixes and matching
          pid="${t%_T}"; pid="${pid%_tumor}"
          for n in "${NORMALS[@]}"; do
            nstem="${n%_N}"; nstem="${nstem%_normal}"
            [[ "$nstem" == "$pid" ]] && { local_match="$n"; break; }
          done
        fi
        if [[ -n "$local_match" ]]; then
          PON_NORMALS+=( "$local_match" )
        else
          warn "No matched normal found for $t — excluded from PoN."
        fi
      done
      # Deduplicate
      mapfile -t PON_NORMALS < <(printf '%s\n' "${PON_NORMALS[@]}" | sort -u)
    else
      PON_NORMALS=( "${NORMALS[@]}" )
      log "PoN build: using all ${#PON_NORMALS[@]} normals."
    fi

    (( ${#PON_NORMALS[@]} )) || err "No normals available to build PoN."

    PON_BAM_LIST="$OUT_DIR/pon/pon_bams.list"; : > "$PON_BAM_LIST"
    for n in "${PON_NORMALS[@]}"; do
      echo "$OUT_DIR/$n/${n}_recal.bam" >> "$PON_BAM_LIST"
    done

    log "Building PoN from ${#PON_NORMALS[@]} normals (normals=${PON_BUILD_NORMALS})..."
    bash "$SCRIPT_DIR/scripts/build_pon.sh" \
      --gatk  "$GATK_BIN" \
      --ref   "$REF" \
      --jobs  "$JOBS" \
      --mode  "$PREPROC_MODE" \
      "$PON_BAM_LIST" "$PON_OUT" >"$LOG_DIR/build_pon.log" 2>&1

    [[ -f "$PON_OUT" ]] || err "PoN build failed — see $LOG_DIR/build_pon.log"
    PON_VCF="$PON_OUT"
    log "PoN ready: $PON_VCF"
  fi
fi

# ============================================================
# STEP 3: Tumor-Normal pairing
# ============================================================
build_pairs() {
  for t in "${TUMORS[@]}"; do
    local match="${SAMPLE_MATCH[$t]:-}"
    if [[ -z "$match" ]]; then
      local pid="${t%_T}"; pid="${pid%_tumor}"
      for n in "${NORMALS[@]}"; do
        local nstem="${n%_N}"; nstem="${nstem%_normal}"
        [[ "$nstem" == "$pid" ]] && { match="$n"; break; }
      done
    fi
    if [[ -n "$match" ]]; then
      echo "${t}:${match}"
    else
      warn "No normal match for tumor $t — skipping."
    fi
  done
}

# ============================================================
# STEP 3: Mutect2 per pair
# ============================================================
run_mutect2_pair() {
  local tumor="$1" normal="$2"
  local pdir="$OUT_DIR/${tumor}_vs_${normal}"; mkdir -p "$pdir"
  local final_vcf="$pdir/mutect2.vcf.gz"
  local pass_vcf="$pdir/mutect2.pass.vcf.gz"
  local sdir="$OUT_DIR/logs/status"; mkdir -p "$sdir"

  case "$PREPROC_MODE" in
    missing)
      if [[ -f "$pass_vcf" && -f "${pass_vcf}.tbi" ]]; then
        log "Skip Mutect2 (pass VCF exists): $pass_vcf"
        : > "$sdir/${tumor}_vs_${normal}.skip"; return
      fi ;;
    ask)
      if [[ -f "$pass_vcf" && -f "${pass_vcf}.tbi" ]]; then
        if [[ -t 0 ]]; then
          read -p "Reuse Mutect2 result for ${tumor}_vs_${normal}? (Y/n) " -r ans; ans=${ans:-Y}
          [[ "$ans" =~ ^[Yy]$ ]] && { log "Reuse: $pass_vcf"; : > "$sdir/${tumor}_vs_${normal}.reuse"; return; }
        else
          log "Reuse (non-interactive): $pass_vcf"; : > "$sdir/${tumor}_vs_${normal}.reuse"; return
        fi
      fi ;;
    all) : ;;
  esac

  IFS='|' read -ra TRG <<<"$(get_rg "$tumor")"
  IFS='|' read -ra NRG <<<"$(get_rg "$normal")"
  local tsm="${TRG[4]:-$tumor}" nsm="${NRG[4]:-$normal}"

  log "Mutect2: $tumor (SM=$tsm) vs $normal (SM=$nsm)"

  local mutect_args=(
    --gatk     "$GATK_BIN"
    --ref      "$REF"
    --tumor-sm "$tsm"
    --normal-sm "$nsm"
  )
  [[ -n "$PON_VCF" && -f "$PON_VCF" ]] && mutect_args+=( --pon "$PON_VCF" )
  [[ -n "$GNOMAD"  && -f "$GNOMAD"  ]] && mutect_args+=( --gnomad "$GNOMAD" )
  $CLEANUP_SOMATIC && mutect_args+=( --cleanup )

  bash "$SCRIPT_DIR/scripts/mutect2_calling.sh" \
    "${mutect_args[@]}" \
    "$OUT_DIR/$tumor/${tumor}_recal.bam" \
    "$OUT_DIR/$normal/${normal}_recal.bam" \
    "$pdir" >"$LOG_DIR/${tumor}_vs_${normal}_mutect2.log" 2>&1

  # PASS filtering
  if [[ -f "$final_vcf" && -f "${final_vcf}.tbi" ]]; then
    bash "$SCRIPT_DIR/scripts/filter_pass.sh" "$final_vcf" "$pass_vcf" \
      >>"$LOG_DIR/${tumor}_vs_${normal}_mutect2.log" 2>&1
    : > "$sdir/${tumor}_vs_${normal}.ok"
    log "Mutect2 OK -> $pass_vcf"
  else
    : > "$sdir/${tumor}_vs_${normal}.fail"
    log "Mutect2 FAILED: ${tumor}_vs_${normal} — see $LOG_DIR/${tumor}_vs_${normal}_mutect2.log"
  fi
}

log "--- STEP 3: Somatic Variant Calling ---"
rm -rf "$OUT_DIR/logs/status"; mkdir -p "$OUT_DIR/logs/status"
mapfile -t PAIRS < <(build_pairs)
(( ${#PAIRS[@]} )) || err "No tumor-normal pairs could be formed. Check sample_map.txt."
log "Pairs to process: ${#PAIRS[@]}"
_run_par_pairs "${PAIRS[@]}"

# ============================================================
# Final summary
# ============================================================
ok=0; fail=0; skip=0; failed_list=()
sdir="$OUT_DIR/logs/status"
for pair in "${PAIRS[@]}"; do
  IFS=: read -r t n <<<"$pair"
  pass_vcf="$OUT_DIR/${t}_vs_${n}/mutect2.pass.vcf.gz"
  if [[ -f "$pass_vcf" && -f "${pass_vcf}.tbi" ]]; then
    (( ok+=1 ))
  elif [[ -f "$sdir/${t}_vs_${n}.skip" || -f "$sdir/${t}_vs_${n}.reuse" ]]; then
    (( skip+=1 ))
  else
    (( fail+=1 )); failed_list+=( "${t} vs ${n}" )
  fi
done

log "==========================================="
log "  Pipeline Summary"
log "  Total pairs : ${#PAIRS[@]}"
log "  OK (PASS)   : $ok"
log "  Skipped     : $skip"
log "  Failed      : $fail"
(( fail > 0 )) && {
  log "  Failed pairs:"
  printf "    %s\n" "${failed_list[@]}"
  log "  Rerun with --pre-proc-missing to retry failed pairs."
}
log "  Results : $OUT_DIR"
log "  Logs    : $LOG_DIR"
log "==========================================="
