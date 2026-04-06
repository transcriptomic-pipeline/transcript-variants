#!/bin/bash
# run_transcript_variants.sh — Main entry point for transcript-variants somatic pipeline.
# Repo: https://github.com/transcriptomic-pipeline/transcript-variants
#
# Mutect2 calling strategy (auto-detected from sample_map.txt):
#   - Matched normal found    -> tumor-normal paired Mutect2
#   - Matched normal missing  -> tumor-only Mutect2 (fallback, with warning)
#   - No normals at all       -> all tumors run in tumor-only mode (auto-detected)
# PoN is applied in both modes if --use-pon / --pon-custom / --pon-build is set.

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
# DEFAULTS
# ============================================================
GATK_BIN=""; PICARD_JAR=""
DEFAULT_REF=""; DEFAULT_DBSNP=""; DEFAULT_MILLS=""; DEFAULT_GNOMAD=""

REF=""; DBSNP=""; MILLS=""; GNOMAD=""
OUT_DIR="$(pwd)/somatic_output"
JOBS=1
PREPROC_MODE="missing"
MARK_DUPLICATES=true
USE_PON=false
PON_MODE=""
PON_BUILD_NORMALS="all"
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
  -r,  --ref PATH          Reference FASTA
  --dbsnp PATH             dbSNP VCF for BQSR
  --mills PATH             Mills indels VCF for BQSR
  --gnomad PATH            gnomAD VCF (used directly, no download triggered)
                           Default: ~/references/GRCh38/gnomad.genomes.v4.1.sites.merged.vcf.gz

Output:
  -o, --out-dir DIR        Output directory (default: ./somatic_output)

Processing:
  --no-mark-dup            Skip MarkDuplicates (default: run)
  -j, --jobs INT           Parallel jobs (default: 1)
  --pre-proc-missing       Reuse existing, run only missing (default)
  --pre-proc-all           Force rerun everything
  --pre-proc-ask           Ask before reusing (interactive)
  -cp, --cleanup-preproc   Remove preprocessing intermediates
  -cs, --cleanup-somatic   Remove unfiltered.vcf.gz after filtering

PoN options (default: no PoN):
  --use-pon                Enable PoN filtering
  --pon-custom PATH        Use existing custom PoN (implies --use-pon)
  --pon-build              Build PoN from all normals (implies --use-pon)
  --pon-build-matched      Build PoN from matched normals only

Calling mode (auto-detected — no flag needed):
  Paired mode:             Used when matched normal is found in sample_map
  Tumor-only mode:         Auto-fallback when no matched normal for a tumor,
                           or when no normals exist at all in sample_map

Tool overrides:
  --gatk PATH
  --picard PATH

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
    -m|--sample-map)       SAMPLE_MAP="$2";       shift 2 ;;
    -rg|--read-group)      RG_FILE="$2";          shift 2 ;;
    -r|--ref)              REF="$2";              shift 2 ;;
    --dbsnp)               DBSNP="$2";            shift 2 ;;
    --mills)               MILLS="$2";            shift 2 ;;
    --gnomad)              GNOMAD="$2";           shift 2 ;;
    -o|--out-dir)          OUT_DIR="$2";          shift 2 ;;
    --no-mark-dup)         MARK_DUPLICATES=false; shift ;;
    -j|--jobs)             JOBS="$2";             shift 2 ;;
    --pre-proc-missing)    PREPROC_MODE="missing"; shift ;;
    --pre-proc-all)        PREPROC_MODE="all";    shift ;;
    --pre-proc-ask)        PREPROC_MODE="ask";    shift ;;
    -cp|--cleanup-preproc) CLEANUP_PREPROC=true;  shift ;;
    -cs|--cleanup-somatic) CLEANUP_SOMATIC=true;  shift ;;
    --use-pon)             USE_PON=true;           shift ;;
    --pon-custom)          CUSTOM_PON="$2"; USE_PON=true; PON_MODE="custom"; shift 2 ;;
    --pon-build)           USE_PON=true; PON_MODE="build"; PON_BUILD_NORMALS="all";     shift ;;
    --pon-build-matched)   USE_PON=true; PON_MODE="build"; PON_BUILD_NORMALS="matched"; shift ;;
    --gatk)                GATK_BIN="$2";         shift 2 ;;
    --picard)              PICARD_JAR="$2";       shift 2 ;;
    --no-auto-install)     AUTO_INSTALL=false;    shift ;;
    --install-only)        INSTALL_ONLY=true;     shift ;;
    -h|--help)             print_usage; exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

# ============================================================
# Source .pipeline_env
# ============================================================
[[ -f "$ENV_FILE" ]] && source "$ENV_FILE"
[[ -z "$REF" ]]    && REF="${DEFAULT_REF:-}"
[[ -z "$DBSNP" ]]  && DBSNP="${DEFAULT_DBSNP:-}"
[[ -z "$MILLS" ]]  && MILLS="${DEFAULT_MILLS:-}"
[[ -z "$GNOMAD" ]] && GNOMAD="${DEFAULT_GNOMAD:-}"

# ============================================================
# Auto-install check
# ============================================================
needs_install() {
  local gbin="${GATK_BIN:-$(command -v gatk 2>/dev/null || true)}"
  [[ -x "$gbin" ]]          || return 0
  local pjar="${PICARD_JAR:-}"
  [[ -f "$pjar" ]]          || command -v picard &>/dev/null || return 0
  command -v samtools &>/dev/null || return 0
  command -v bcftools &>/dev/null || return 0
  return 1
}

if $INSTALL_ONLY; then
  log "Running install.sh..."; bash "$SCRIPT_DIR/install.sh"; exit 0
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

[[ -z "$GATK_BIN" ]] && GATK_BIN="${DEFAULT_GATK_BIN:-$(command -v gatk 2>/dev/null || true)}"
[[ -z "$PICARD_JAR" ]] && PICARD_JAR="${DEFAULT_PICARD_JAR:-}"
[[ -x "$GATK_BIN" ]] || err "GATK not found. Run install.sh or provide --gatk PATH."
[[ -f "$REF" ]] || err "Reference FASTA not found: '${REF}'. Run install.sh or --ref PATH."
[[ -f "$DBSNP" ]] || warn "dbSNP not found — BQSR known-sites incomplete."
[[ -f "$MILLS" ]] || warn "Mills VCF not found — BQSR known-sites incomplete."

# gnomAD: try default path, warn if still missing
if [[ -z "$GNOMAD" || ! -f "$GNOMAD" ]]; then
  GNOMAD_DEFAULT="$HOME/references/GRCh38/gnomad.genomes.v4.1.sites.merged.vcf.gz"
  if [[ -f "$GNOMAD_DEFAULT" ]]; then
    GNOMAD="$GNOMAD_DEFAULT"
    log "gnomAD auto-detected: $GNOMAD"
  else
    warn "gnomAD not found — Mutect2 will run without germline resource."
    warn "To install: bash install.sh  (downloads gnomAD v4.1 by default)"
    GNOMAD=""
  fi
fi

# ============================================================
# Output + logging
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

[[ -f "$SAMPLE_MAP" ]] || err "Sample map not found: $SAMPLE_MAP"

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

(( ${#ALL_SAMPLES[@]} )) || err "No valid samples in $SAMPLE_MAP."

TUMORS=(); NORMALS=()
for sm in "${ALL_SAMPLES[@]}"; do
  case "${SAMPLE_TYPE[$sm]}" in
    tumor)  TUMORS+=( "$sm" ) ;;
    normal) NORMALS+=( "$sm" ) ;;
    *)      warn "Cannot classify '$sm' — skipping." ;;
  esac
done

log "Tumors  (${#TUMORS[@]}): ${TUMORS[*]:-none}"
log "Normals (${#NORMALS[@]}): ${NORMALS[*]:-none}"

# ============================================================
# Load RG details
# ============================================================
declare -A RG_MAP
if [[ -f "$RG_FILE" ]]; then
  log "Loading RG details: $RG_FILE"
  while IFS=$'\t' read -r sample rgid rglb rgpl rgpu rgsm; do
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue
    RG_MAP["$sample"]="${rgid}|${rglb}|${rgpl}|${rgpu}|${rgsm}"
  done <"$RG_FILE"
else
  warn "RG file not found: $RG_FILE — defaults will be used."
fi

get_rg() {
  local sm="$1"
  [[ -n "${RG_MAP[$sm]:-}" ]] && echo "${RG_MAP[$sm]}" \
                               || echo "${sm}|lib1|ILLUMINA|${sm}.unit1|${sm}"
}

# ============================================================
# Parallel runners
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
_run_par_jobs() {
  # Runs arbitrary bash -c strings in parallel
  local -a pids=(); local use_waitn=false; supports_waitn && use_waitn=true
  for job_str in "$@"; do
    bash -c "$job_str" & pids+=( $! )
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
          read -p "Reuse preprocessed BAM for $sm? (Y/n) " -r ans; ans=${ans:-Y}
          [[ "$ans" =~ ^[Yy]$ ]] && { log "Reuse: $out_recal"; return; }
        else
          log "Reuse (non-interactive): $out_recal"; return
        fi
      fi ;;
    all) : ;;
  esac

  IFS='|' read -ra RGP <<<"$(get_rg "$sm")"
  log "Preprocessing $sm ..."
  bash "$SCRIPT_DIR/scripts/gatk_preprocess.sh" \
    --gatk     "$GATK_BIN"          \
    --picard   "${PICARD_JAR:-}"    \
    --ref      "$REF"               \
    --dbsnp    "${DBSNP:-}"         \
    --mills    "${MILLS:-}"         \
    --mark-dup "$MARK_DUPLICATES"   \
    --sample   "$sm"                \
    --input    "$bam"               \
    --out-dir  "$sdir"              \
    --rgid     "${RGP[0]}"          \
    --rglb     "${RGP[1]}"          \
    --rgpl     "${RGP[2]}"          \
    --rgpu     "${RGP[3]}"          \
    --rgsm     "${RGP[4]}"          \
    >"$LOG_DIR/${sm}_preprocess.log" 2>&1 \
    || err "Preprocessing FAILED for $sm — see $LOG_DIR/${sm}_preprocess.log"

  [[ -f "$out_recal" ]] || err "No recal BAM produced for $sm."
  log "Preprocessing OK: $out_recal"
}

log "--- STEP 1: Preprocessing (${#ALL_SAMPLES[@]} samples) ---"
(( ${#ALL_SAMPLES[@]} )) && _run_par_items preprocess_sample "${ALL_SAMPLES[@]}"

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
    if (( ${#NORMALS[@]} == 0 )); then
      warn "No normal samples found in sample_map — cannot build PoN. Proceeding without PoN."
      USE_PON=false
    else
      PON_OUT="$OUT_DIR/pon/custom_pon.vcf.gz"; mkdir -p "$OUT_DIR/pon"

      # Determine which normals to include
      declare -a PON_NORMALS=()
      if [[ "$PON_BUILD_NORMALS" == "matched" ]]; then
        log "PoN build: matched normals only."
        for t in "${TUMORS[@]}"; do
          local_match="${SAMPLE_MATCH[$t]:-}"
          if [[ -z "$local_match" ]]; then
            local pid="${t%_T}"; pid="${pid%_tumor}"
            for n in "${NORMALS[@]}"; do
              local nstem="${n%_N}"; nstem="${nstem%_normal}"
              [[ "$nstem" == "$pid" ]] && { local_match="$n"; break; }
            done
          fi
          [[ -n "$local_match" ]] && PON_NORMALS+=( "$local_match" ) \
                                  || warn "No matched normal for $t — excluded from PoN."
        done
        mapfile -t PON_NORMALS < <(printf '%s\n' "${PON_NORMALS[@]}" | sort -u)
      else
        PON_NORMALS=( "${NORMALS[@]}" )
        log "PoN build: all ${#PON_NORMALS[@]} normals."
      fi

      if (( ${#PON_NORMALS[@]} == 0 )); then
        warn "No qualifying normals for PoN build — proceeding without PoN."
        USE_PON=false
      else
        PON_BAM_LIST="$OUT_DIR/pon/pon_bams.list"; : > "$PON_BAM_LIST"
        for n in "${PON_NORMALS[@]}"; do
          echo "$OUT_DIR/$n/${n}_recal.bam" >> "$PON_BAM_LIST"
        done

        bash "$SCRIPT_DIR/scripts/build_pon.sh" \
          --gatk "$GATK_BIN" --ref "$REF" \
          --jobs "$JOBS" --mode "$PREPROC_MODE" \
          "$PON_BAM_LIST" "$PON_OUT" >"$LOG_DIR/build_pon.log" 2>&1

        [[ -f "$PON_OUT" ]] || err "PoN build failed — see $LOG_DIR/build_pon.log"
        PON_VCF="$PON_OUT"
        log "PoN ready: $PON_VCF"
      fi
    fi
  fi
fi

# ============================================================
# STEP 3: Build calling jobs (paired OR tumor-only per tumor)
# ============================================================

# Returns matched normal for a tumor, or empty string if none
find_normal_for_tumor() {
  local t="$1"
  # Explicit mapping from sample_map
  local match="${SAMPLE_MATCH[$t]:-}"
  if [[ -n "$match" ]]; then
    # Validate the normal exists in our sample list
    [[ -n "${SAMPLE_BAM[$match]:-}" ]] && { echo "$match"; return; }
    warn "MATCHED_WITH '$match' for $t not found in sample_map — trying auto-pair."
  fi
  # Auto-pair by stripping known suffixes and matching stem
  local pid="${t%_T}"; pid="${pid%_tumor}"
  for n in "${NORMALS[@]}"; do
    local nstem="${n%_N}"; nstem="${nstem%_normal}"
    [[ "$nstem" == "$pid" ]] && { echo "$n"; return; }
  done
  echo ""   # no match
}

# ---- Decide calling mode for each tumor ----
# CALL_JOBS is an array of strings: "TUMOR:NORMAL" or "TUMOR:" (tumor-only)
declare -a CALL_JOBS=()
PAIRED_COUNT=0; TUMOR_ONLY_COUNT=0

for t in "${TUMORS[@]}"; do
  normal=$(find_normal_for_tumor "$t")
  if [[ -n "$normal" ]]; then
    CALL_JOBS+=( "${t}:${normal}" )
    (( PAIRED_COUNT+=1 ))
    log "  Calling plan: $t vs $normal  [PAIRED]"
  else
    CALL_JOBS+=( "${t}:" )
    (( TUMOR_ONLY_COUNT+=1 ))
    warn "  Calling plan: $t  [TUMOR-ONLY — no matched normal found]"
  fi
done

log "Calling plan: ${PAIRED_COUNT} paired | ${TUMOR_ONLY_COUNT} tumor-only"

# ============================================================
# STEP 3: Run Mutect2 per job
# ============================================================
run_calling_job() {
  local job="$1"
  local tumor="${job%%:*}"
  local normal="${job##*:}"    # empty string if tumor-only
  local is_paired=false
  [[ -n "$normal" ]] && is_paired=true

  local pdir
  if $is_paired; then
    pdir="$OUT_DIR/${tumor}_vs_${normal}"
  else
    pdir="$OUT_DIR/${tumor}_tumor_only"
  fi
  mkdir -p "$pdir"

  local final_vcf="$pdir/mutect2.vcf.gz"
  local pass_vcf="$pdir/mutect2.pass.vcf.gz"
  local sdir="$OUT_DIR/logs/status"; mkdir -p "$sdir"
  local status_key="${job//\//_}"   # safe filename from job key

  case "$PREPROC_MODE" in
    missing)
      if [[ -f "$pass_vcf" && -f "${pass_vcf}.tbi" ]]; then
        log "Skip (pass VCF exists): $pass_vcf"
        : > "$sdir/${status_key}.skip"; return
      fi ;;
    ask)
      if [[ -f "$pass_vcf" && -f "${pass_vcf}.tbi" ]]; then
        local label="$tumor"; $is_paired && label+=" vs $normal"
        if [[ -t 0 ]]; then
          read -p "Reuse Mutect2 for $label? (Y/n) " -r ans; ans=${ans:-Y}
          [[ "$ans" =~ ^[Yy]$ ]] && { log "Reuse: $pass_vcf"; : > "$sdir/${status_key}.reuse"; return; }
        else
          log "Reuse (non-interactive): $pass_vcf"; : > "$sdir/${status_key}.reuse"; return
        fi
      fi ;;
    all) : ;;
  esac

  IFS='|' read -ra TRG <<<"$(get_rg "$tumor")"
  local tsm="${TRG[4]:-$tumor}"

  local mutect_args=(
    --gatk      "$GATK_BIN"
    --ref       "$REF"
    --tumor-sm  "$tsm"
  )
  [[ -n "$PON_VCF" && -f "$PON_VCF" ]]   && mutect_args+=( --pon "$PON_VCF" )
  [[ -n "$GNOMAD"  && -f "$GNOMAD"  ]]   && mutect_args+=( --gnomad "$GNOMAD" )
  $CLEANUP_SOMATIC                        && mutect_args+=( --cleanup )

  local log_label
  if $is_paired; then
    IFS='|' read -ra NRG <<<"$(get_rg "$normal")"
    local nsm="${NRG[4]:-$normal}"
    mutect_args+=( --normal-sm "$nsm" )
    log_label="${tumor}_vs_${normal}"
    log "Mutect2 PAIRED: $tsm vs $nsm"
    bash "$SCRIPT_DIR/scripts/mutect2_calling.sh" \
      "${mutect_args[@]}" \
      "$OUT_DIR/$tumor/${tumor}_recal.bam" \
      "$OUT_DIR/$normal/${normal}_recal.bam" \
      "$pdir" >"$LOG_DIR/${log_label}_mutect2.log" 2>&1
  else
    mutect_args+=( --tumor-only )
    log_label="${tumor}_tumor_only"
    log "Mutect2 TUMOR-ONLY: $tsm"
    bash "$SCRIPT_DIR/scripts/mutect2_calling.sh" \
      "${mutect_args[@]}" \
      "$OUT_DIR/$tumor/${tumor}_recal.bam" \
      "$pdir" >"$LOG_DIR/${log_label}_mutect2.log" 2>&1
  fi

  if [[ -f "$final_vcf" && -f "${final_vcf}.tbi" ]]; then
    bash "$SCRIPT_DIR/scripts/filter_pass.sh" "$final_vcf" "$pass_vcf" \
      >>"$LOG_DIR/${log_label}_mutect2.log" 2>&1
    : > "$sdir/${status_key}.ok"
    log "Mutect2 OK -> $pass_vcf"
  else
    : > "$sdir/${status_key}.fail"
    log "Mutect2 FAILED: $log_label — see $LOG_DIR/${log_label}_mutect2.log"
  fi
}

log "--- STEP 3: Somatic Variant Calling (${#CALL_JOBS[@]} jobs) ---"
rm -rf "$OUT_DIR/logs/status"; mkdir -p "$OUT_DIR/logs/status"

(( ${#CALL_JOBS[@]} )) || err "No calling jobs to run."

# Run all jobs in parallel (up to $JOBS at once)
declare -a pids=()
use_waitn=false; supports_waitn && use_waitn=true
for job in "${CALL_JOBS[@]}"; do
  ( run_calling_job "$job" ) & pids+=( $! )
  if (( ${#pids[@]} >= JOBS )); then
    $use_waitn && { wait -n || true; } \
               || { wait "${pids[0]}" || true; pids=( "${pids[@]:1}" ); }
  fi
done
for pid in "${pids[@]}"; do wait "$pid" || true; done

# ============================================================
# Final summary — based on actual output files
# ============================================================
ok_paired=0; ok_tonly=0; fail=0; skip=0; failed_list=()
sdir="$OUT_DIR/logs/status"

for job in "${CALL_JOBS[@]}"; do
  local tumor="${job%%:*}"
  local normal="${job##*:}"
  local is_paired=false; [[ -n "$normal" ]] && is_paired=true
  local status_key="${job//\//_}"

  local pdir
  $is_paired && pdir="$OUT_DIR/${tumor}_vs_${normal}" \
             || pdir="$OUT_DIR/${tumor}_tumor_only"
  local pass_vcf="$pdir/mutect2.pass.vcf.gz"

  if [[ -f "$pass_vcf" && -f "${pass_vcf}.tbi" ]]; then
    $is_paired && (( ok_paired+=1 )) || (( ok_tonly+=1 ))
  elif [[ -f "$sdir/${status_key}.skip" || -f "$sdir/${status_key}.reuse" ]]; then
    (( skip+=1 ))
  else
    (( fail+=1 ))
    $is_paired && failed_list+=( "$tumor vs $normal [paired]" ) \
               || failed_list+=( "$tumor [tumor-only]" )
  fi
done

log "==========================================="
log "  Pipeline Summary"
log "  Total jobs      : ${#CALL_JOBS[@]}"
log "  OK paired       : $ok_paired"
log "  OK tumor-only   : $ok_tonly"
log "  Skipped/reused  : $skip"
log "  Failed          : $fail"
(( fail > 0 )) && {
  log "  Failed jobs:"
  printf "    %s\n" "${failed_list[@]}"
  log "  Rerun with --pre-proc-missing to retry only failed jobs."
}
log "  Results : $OUT_DIR"
log "  Logs    : $LOG_DIR"
log "==========================================="
