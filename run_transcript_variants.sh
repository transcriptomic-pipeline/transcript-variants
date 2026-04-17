#!/bin/bash
# run_transcript_variants.sh — Somatic variant calling pipeline
# Key fixes:
#   1. No exported-function/xargs dependency
#   2. Native bash parallel job control
#   3. Mutect2 uses SM from BAM headers, not directory/sample IDs
#   4. Single calling-plan build reused across steps
#   5. Adds --skip-preprocessed to avoid re-running finished preprocessing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="output"
LOG_DIR=""
MAIN_LOG=""

log()  { echo "[$(date '+%F %T')] $*" | tee -a "$MAIN_LOG"; }
warn() { echo "[$(date '+%F %T')] WARN: $*" | tee -a "$MAIN_LOG" >&2; }
err()  { echo "[$(date '+%F %T')] ERROR: $*" | tee -a "$MAIN_LOG" >&2; exit 1; }

GATK_BIN=""
PICARD_JAR=""
REF_FA=""
DBSNP_VCF=""
MILLS_VCF=""
GNOMAD_VCF=""
SAMTOOLS_BIN=""

JOBS=4
MARKDUP=true
SKIP_PREPROCESSED=false

USE_PON=false
PON_MODE="n/a"
PON_NORMALS="all"
PON_VCF=""

SAMPLE_MAP="$SCRIPT_DIR/config/sample_map.txt"
RG_FILE="$SCRIPT_DIR/config/samples_rg_details.txt"

ENV_FILE="$SCRIPT_DIR/.pipeline_env"
if [[ -f "$ENV_FILE" ]]; then
    # shellcheck source=/dev/null
    source "$ENV_FILE"
    GATK_BIN="${GATK_BIN:-${GATK_BIN:-}}"
    PICARD_JAR="${PICARD_JAR:-${PICARD_JAR:-}}"
    REF_FA="${REF_FA:-${DEFAULT_REF:-}}"
    DBSNP_VCF="${DBSNP_VCF:-${DEFAULT_DBSNP:-}}"
    MILLS_VCF="${MILLS_VCF:-${DEFAULT_MILLS:-}}"
    GNOMAD_VCF="${GNOMAD_VCF:-${DEFAULT_GNOMAD:-}}"
fi

print_usage() {
    cat <<EOF
Usage: run_transcript_variants.sh [OPTIONS]

  -o, --out-dir DIR        Output directory (default: output)
  -j, --jobs N             Parallel jobs (default: 4)

  --ref FA                 Reference FASTA
  --dbsnp VCF              dbSNP VCF
  --mills VCF              Mills indels VCF
  --gnomad VCF             gnomAD VCF
  --gatk PATH              GATK binary
  --picard JAR             Picard JAR
  --samtools PATH          samtools binary

  --sample-map FILE        Sample map (default: config/sample_map.txt)
  --rg-file FILE           RG details file (default: config/samples_rg_details.txt)

  --no-markdup             Skip Picard MarkDuplicates
  --skip-preprocessed      Do not schedule samples that already have *_recal.bam + *.bai
  --pon VCF                Existing panel-of-normals VCF
  --pon-mode create|use    Only 'use' is supported in this script
  --pon-normals all|FILE   Reserved for future use

  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--out-dir)          OUT_DIR="$2"; shift 2 ;;
        -j|--jobs)             JOBS="$2"; shift 2 ;;

        --ref)                 REF_FA="$2"; shift 2 ;;
        --dbsnp)               DBSNP_VCF="$2"; shift 2 ;;
        --mills)               MILLS_VCF="$2"; shift 2 ;;
        --gnomad)              GNOMAD_VCF="$2"; shift 2 ;;
        --gatk)                GATK_BIN="$2"; shift 2 ;;
        --picard)              PICARD_JAR="$2"; shift 2 ;;
        --samtools)            SAMTOOLS_BIN="$2"; shift 2 ;;

        --sample-map)          SAMPLE_MAP="$2"; shift 2 ;;
        --rg-file)             RG_FILE="$2"; shift 2 ;;

        --no-markdup)          MARKDUP=false; shift ;;
        --skip-preprocessed)   SKIP_PREPROCESSED=true; shift ;;

        --pon)                 PON_VCF="$2"; USE_PON=true; shift 2 ;;
        --pon-mode)            PON_MODE="$2"; shift 2 ;;
        --pon-normals)         PON_NORMALS="$2"; shift 2 ;;

        -h|--help)             print_usage; exit 0 ;;
        *)                     echo "Unknown option: $1" >&2; print_usage; exit 1 ;;
    esac
done

[[ "$JOBS" =~ ^[0-9]+$ ]] || err "--jobs must be an integer"
(( JOBS >= 1 )) || err "--jobs must be >= 1"

if [[ -z "$GATK_BIN" ]]; then
    GATK_BIN="$(command -v gatk 2>/dev/null || true)"
fi
if [[ -z "$SAMTOOLS_BIN" ]]; then
    SAMTOOLS_BIN="$(command -v samtools 2>/dev/null || true)"
fi

mkdir -p "$OUT_DIR"
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$LOG_DIR"
MAIN_LOG="$OUT_DIR/pipeline.log"
: > "$MAIN_LOG"

trap 'echo "[$(date "+%F %T")] ERROR: Unexpected error at line $LINENO: $BASH_COMMAND" | tee -a "$MAIN_LOG" >&2' ERR

[[ -n "$GATK_BIN" && -x "$GATK_BIN" ]] || err "GATK binary not found/executable: ${GATK_BIN:-<unset>}"
[[ -n "$PICARD_JAR" && -f "$PICARD_JAR" ]] || err "Picard JAR not found: ${PICARD_JAR:-<unset>}"
[[ -n "$SAMTOOLS_BIN" && -x "$SAMTOOLS_BIN" ]] || err "samtools not found/executable: ${SAMTOOLS_BIN:-<unset>}"

[[ -n "$REF_FA" && -f "$REF_FA" ]] || err "Reference FASTA not found: ${REF_FA:-<unset>}"
[[ -f "${REF_FA}.fai" ]] || err "Reference FASTA index missing: ${REF_FA}.fai"

dict_path="${REF_FA%.*}.dict"
[[ -f "$dict_path" ]] || warn "Sequence dictionary not found: $dict_path"

[[ -n "$DBSNP_VCF" && -f "$DBSNP_VCF" ]] || err "dbSNP VCF not found: ${DBSNP_VCF:-<unset>}"
[[ -f "${DBSNP_VCF}.tbi" ]] || err "dbSNP index missing: ${DBSNP_VCF}.tbi"

[[ -n "$MILLS_VCF" && -f "$MILLS_VCF" ]] || err "Mills VCF not found: ${MILLS_VCF:-<unset>}"
[[ -f "${MILLS_VCF}.tbi" ]] || err "Mills index missing: ${MILLS_VCF}.tbi"

if [[ -n "${GNOMAD_VCF:-}" ]]; then
    [[ -f "$GNOMAD_VCF" ]] || err "gnomAD VCF not found: $GNOMAD_VCF"
    [[ -f "${GNOMAD_VCF}.tbi" ]] || err "gnomAD index missing: ${GNOMAD_VCF}.tbi"
fi

[[ -f "$SAMPLE_MAP" ]] || err "Sample map not found: $SAMPLE_MAP"
[[ -f "$RG_FILE" ]] || err "RG details file not found: $RG_FILE"

if $USE_PON; then
    [[ "$PON_MODE" == "use" || "$PON_MODE" == "n/a" ]] || err "Only --pon-mode use is supported"
    [[ -n "$PON_VCF" && -f "$PON_VCF" ]] || err "PoN VCF not found: ${PON_VCF:-<unset>}"
    [[ -f "${PON_VCF}.tbi" ]] || err "PoN index missing: ${PON_VCF}.tbi"
fi

log "==========================================="
log "  transcript-variants somatic pipeline"
log "==========================================="
log "  GATK              : $GATK_BIN"
log "  Picard            : $PICARD_JAR"
log "  samtools          : $SAMTOOLS_BIN"
log "  Ref               : $REF_FA"
log "  dbSNP             : $DBSNP_VCF"
log "  Mills             : $MILLS_VCF"
log "  gnomAD            : ${GNOMAD_VCF:-<none>}"
log "  MarkDup           : $MARKDUP"
log "  Skip preprocessed : $SKIP_PREPROCESSED"
log "  Use PoN           : $USE_PON | mode=$PON_MODE | normals=$PON_NORMALS"
log "  Jobs              : $JOBS"
log "  Out dir           : $OUT_DIR"

declare -A SAMPLE_BAM=()
declare -A SAMPLE_TYPE=()
declare -A SAMPLE_PAIR=()

declare -A RG_ID=()
declare -A RG_SM=()
declare -A RG_LB=()
declare -A RG_PL=()
declare -A RG_PU=()

declare -a TUMORS=()
declare -a NORMALS=()
declare -a PAIRED_JOBS=()
declare -a TUMOR_ONLY_JOBS=()

FORCE_TUMOR_ONLY=false

infer_type() {
    local id="${1,,}"
    if [[ "$id" =~ (^|[_.-])t($|[_.-]) ]] \
        || [[ "$id" =~ _tumor$ ]] \
        || [[ "$id" =~ ^tumor([_.-]|$) ]] \
        || [[ "$id" =~ _tumour$ ]] \
        || [[ "$id" =~ ^tumour([_.-]|$) ]] \
        || [[ "$id" =~ [_.-]tum($|[_.-]) ]]; then
        echo "tumor"
        return
    fi

    if [[ "$id" =~ (^|[_.-])n($|[_.-]) ]] \
        || [[ "$id" =~ _normal$ ]] \
        || [[ "$id" =~ ^normal([_.-]|$) ]] \
        || [[ "$id" =~ _norm($|[_.-]) ]]; then
        echo "normal"
        return
    fi

    echo ""
}

strip_suffix() {
    local id="$1"
    id="${id%_T}"
    id="${id%_N}"
    id="${id%_tumor}"
    id="${id%_normal}"
    id="${id%_tumour}"
    id="${id%_norm}"
    id="${id%_tum}"
    echo "$id"
}

get_bam_sm() {
    local bam="$1"
    local -a sms=()

    mapfile -t sms < <(
        "$SAMTOOLS_BIN" view -H "$bam" \
        | awk -F'\t' '
            /^@RG/ {
                for (i=1; i<=NF; i++) {
                    if ($i ~ /^SM:/) {
                        sm=$i
                        sub(/^SM:/, "", sm)
                        print sm
                    }
                }
            }
        ' | sort -u
    )

    if (( ${#sms[@]} == 0 )); then
        return 2
    elif (( ${#sms[@]} > 1 )); then
        printf '%s\n' "${sms[@]}" >&2
        return 3
    fi

    echo "${sms[0]}"
}

wait_for_slot() {
    while (( $(jobs -rp | wc -l) >= JOBS )); do
        wait -n || true
    done
}

run_parallel_items() {
    local func_name="$1"
    local fail_file="$2"
    shift 2

    : > "$fail_file"

    local item
    for item in "$@"; do
        wait_for_slot
        (
            if ! "$func_name" "$item"; then
                echo "$item" >> "$fail_file"
                exit 1
            fi
        ) &
    done

    wait || true
}

log "Loading RG details: $RG_FILE"
# Expected columns (tab-separated):
# SAMPLE  RGID  RGLB  RGPL  RGPU  RGSM
while read -r sample_id rgid rglb rgpl rgpu rgsm _rest; do
  [[ -z "${sample_id:-}" ]] && continue
  [[ "$sample_id" =~ ^# ]] && continue

  RG_ID["$sample_id"]="$rgid"   # Read Group ID
  RG_LB["$sample_id"]="$rglb"   # Library ID (LIBxxxxx)
  RG_PL["$sample_id"]="$rgpl"   # Platform (ILLUMINA, etc.)
  RG_PU["$sample_id"]="$rgpu"   # Platform unit
  RG_SM["$sample_id"]="$rgsm"   # Sample name in BAM header (1493_13_T, etc.)
done < "$RG_FILE"

log "Loading sample map: $SAMPLE_MAP"
while IFS=$'\t' read -r sample_id type bam_path matched_with || [[ -n "${sample_id:-}" ]]; do
    [[ -z "${sample_id:-}" ]] && continue
    [[ "$sample_id" =~ ^# ]] && continue

    if [[ -z "${bam_path:-}" ]]; then
        warn "No BAM path for '$sample_id' — skipping"
        continue
    fi
    if [[ ! -f "$bam_path" ]]; then
        warn "BAM not found for '$sample_id': $bam_path — skipping"
        continue
    fi

    role=""
    lower_type="${type,,}"

    if [[ "$lower_type" == "tumor" || "$lower_type" == "tumour" ]]; then
        role="tumor"
    elif [[ "$lower_type" == "normal" || "$lower_type" == "norm" ]]; then
        role="normal"
    else
        role="$(infer_type "$sample_id")"
        if [[ -z "$role" ]]; then
            warn "Cannot determine role for '$sample_id' (TYPE='$type') — skipping"
            continue
        fi
        log "  Auto-detected role for '$sample_id': $role"
    fi

    SAMPLE_BAM["$sample_id"]="$bam_path"
    SAMPLE_TYPE["$sample_id"]="$role"

    if [[ -n "${matched_with:-}" ]]; then
        SAMPLE_PAIR["$sample_id"]="$matched_with"
    fi

    if [[ "$role" == "tumor" ]]; then
        TUMORS+=("$sample_id")
    else
        NORMALS+=("$sample_id")
    fi
done < "$SAMPLE_MAP"

(( ${#TUMORS[@]} > 0 )) || err "No tumor samples loaded from $SAMPLE_MAP"

if (( ${#NORMALS[@]} == 0 )); then
    FORCE_TUMOR_ONLY=true
    warn "No normal samples found — forcing tumor-only mode for all tumors"
fi

log "Loaded  : ${#TUMORS[@]} tumor(s), ${#NORMALS[@]} normal(s)"

find_normal_match() {
    local tumor="$1"
    local candidate=""
    local base_t=""
    local base_n=""
    local n=""

    if $FORCE_TUMOR_ONLY; then
        echo ""
        return 0
    fi

    if [[ -n "${SAMPLE_PAIR[$tumor]:-}" ]]; then
        candidate="${SAMPLE_PAIR[$tumor]}"
        if [[ -n "${SAMPLE_BAM[$candidate]:-}" && "${SAMPLE_TYPE[$candidate]:-}" == "normal" ]]; then
            echo "$candidate"
            return 0
        else
            warn "Explicit MATCHED_WITH '$candidate' for '$tumor' not found or not normal"
        fi
    fi

    base_t="$(strip_suffix "$tumor")"
    for n in "${NORMALS[@]}"; do
        base_n="$(strip_suffix "$n")"
        if [[ -n "$base_t" && "$base_t" == "$base_n" ]]; then
            echo "$n"
            return 0
        fi
    done

    echo ""
}

log "--- Building calling plan ---"
for tumor in "${TUMORS[@]}"; do
    matched="$(find_normal_match "$tumor")"
    if [[ -n "$matched" ]]; then
        PAIRED_JOBS+=("${tumor}:${matched}")
        log "  $tumor -> $matched [PAIRED]"
    else
        TUMOR_ONLY_JOBS+=("$tumor")
        log "  $tumor -> [TUMOR-ONLY]"
    fi
done

log "Calling plan: ${#PAIRED_JOBS[@]} paired | ${#TUMOR_ONLY_JOBS[@]} tumor-only"
(( ${#PAIRED_JOBS[@]} + ${#TUMOR_ONLY_JOBS[@]} > 0 )) || err "Calling plan is empty"

preprocess_sample() {
    local sample="$1"
    local input_bam="${SAMPLE_BAM[$sample]:-}"
    local sdir="$OUT_DIR/$sample"
    local slog="$LOG_DIR/${sample}_preprocess.log"

    local sorted_bam="$sdir/${sample}_sorted.bam"
    local dedup_bam="$sdir/${sample}_dedup.bam"
    local bqsr_table="$sdir/${sample}_bqsr.table"
    local recal_bam="$sdir/${sample}_recal.bam"
    local dup_metrics="$sdir/${sample}_dup_metrics.txt"

    local rg_id="${RG_ID[$sample]:-}"
    local rg_sm="${RG_SM[$sample]:-}"
    local rg_lb="${RG_LB[$sample]:-}"
    local rg_pl="${RG_PL[$sample]:-}"
    local rg_pu="${RG_PU[$sample]:-}"

    mkdir -p "$sdir"
    : > "$slog"

    if [[ -f "$recal_bam" && -f "${recal_bam}.bai" ]]; then
        log "Preprocessing already done: $sample"
        return 0
    fi

    [[ -n "$input_bam" ]] || { echo "[$(date '+%F %T')] ERROR: No BAM path for $sample" | tee -a "$MAIN_LOG" >&2; return 1; }
    [[ -f "$input_bam" ]] || { echo "[$(date '+%F %T')] ERROR: Input BAM not found for $sample: $input_bam" | tee -a "$MAIN_LOG" >&2; return 1; }

    [[ -n "$rg_id" ]] || { echo "[$(date '+%F %T')] ERROR: No RG_ID for $sample in $RG_FILE" | tee -a "$MAIN_LOG" >&2; return 1; }
    [[ -n "$rg_sm" ]] || { echo "[$(date '+%F %T')] ERROR: No RG_SM for $sample in $RG_FILE" | tee -a "$MAIN_LOG" >&2; return 1; }
    [[ -n "$rg_lb" ]] || { echo "[$(date '+%F %T')] ERROR: No RG_LB for $sample in $RG_FILE" | tee -a "$MAIN_LOG" >&2; return 1; }
    [[ -n "$rg_pl" ]] || { echo "[$(date '+%F %T')] ERROR: No RG_PL for $sample in $RG_FILE" | tee -a "$MAIN_LOG" >&2; return 1; }
    [[ -n "$rg_pu" ]] || { echo "[$(date '+%F %T')] ERROR: No RG_PU for $sample in $RG_FILE" | tee -a "$MAIN_LOG" >&2; return 1; }

    log "Preprocessing $sample ..."

    if [[ ! -f "$sorted_bam" ]]; then
        java -jar "$PICARD_JAR" AddOrReplaceReadGroups \
            I="$input_bam" \
            O="$sorted_bam" \
            SORT_ORDER=coordinate \
            RGID="$rg_id" \
            RGSM="$rg_sm" \
            RGLB="$rg_lb" \
            RGPL="$rg_pl" \
            RGPU="$rg_pu" \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT \
            >> "$slog" 2>&1 || {
                echo "[$(date '+%F %T')] ERROR: AddOrReplaceReadGroups failed for $sample" | tee -a "$MAIN_LOG" >&2
                return 1
            }
    fi

    if $MARKDUP; then
        if [[ ! -f "$dedup_bam" ]]; then
            java -jar "$PICARD_JAR" MarkDuplicates \
                I="$sorted_bam" \
                O="$dedup_bam" \
                M="$dup_metrics" \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=LENIENT \
                >> "$slog" 2>&1 || {
                    echo "[$(date '+%F %T')] ERROR: MarkDuplicates failed for $sample" | tee -a "$MAIN_LOG" >&2
                    return 1
                }
        fi
    else
        dedup_bam="$sorted_bam"
        [[ -f "${dedup_bam}.bai" ]] || "$SAMTOOLS_BIN" index -@ 2 "$dedup_bam" >> "$slog" 2>&1
    fi

    if [[ ! -f "$bqsr_table" ]]; then
        "$GATK_BIN" BaseRecalibrator \
            -I "$dedup_bam" \
            -R "$REF_FA" \
            --known-sites "$DBSNP_VCF" \
            --known-sites "$MILLS_VCF" \
            -O "$bqsr_table" \
            >> "$slog" 2>&1 || {
                echo "[$(date '+%F %T')] ERROR: BaseRecalibrator failed for $sample" | tee -a "$MAIN_LOG" >&2
                return 1
            }
    fi

    if [[ ! -f "$recal_bam" ]]; then
        "$GATK_BIN" ApplyBQSR \
            -I "$dedup_bam" \
            -R "$REF_FA" \
            --bqsr-recal-file "$bqsr_table" \
            -O "$recal_bam" \
            >> "$slog" 2>&1 || {
                echo "[$(date '+%F %T')] ERROR: ApplyBQSR failed for $sample" | tee -a "$MAIN_LOG" >&2
                return 1
            }
    fi

    if [[ ! -f "${recal_bam}.bai" ]]; then
        "$SAMTOOLS_BIN" index -@ 2 "$recal_bam" >> "$slog" 2>&1 || {
            echo "[$(date '+%F %T')] ERROR: samtools index failed for $sample recal BAM" | tee -a "$MAIN_LOG" >&2
            return 1
        }
    fi

    log "Preprocessing done: $sample"
}

ALL_SAMPLES=( "${TUMORS[@]}" "${NORMALS[@]}" )
PREPROCESS_QUEUE=()

if $SKIP_PREPROCESSED; then
    for sample in "${ALL_SAMPLES[@]}"; do
        recal_bam="$OUT_DIR/$sample/${sample}_recal.bam"
        if [[ -f "$recal_bam" && -f "${recal_bam}.bai" ]]; then
            log "Skipping completed preprocessing: $sample"
        else
            PREPROCESS_QUEUE+=( "$sample" )
        fi
    done
else
    PREPROCESS_QUEUE=( "${ALL_SAMPLES[@]}" )
fi

log "--- STEP 1: Preprocessing (${#PREPROCESS_QUEUE[@]} scheduled / ${#ALL_SAMPLES[@]} total) ---"

if (( ${#PREPROCESS_QUEUE[@]} > 0 )); then
    PREPROC_FAIL_FILE="$OUT_DIR/preproc_failures.txt"
    run_parallel_items preprocess_sample "$PREPROC_FAIL_FILE" "${PREPROCESS_QUEUE[@]}"

    if [[ -s "$PREPROC_FAIL_FILE" ]]; then
        while IFS= read -r sample; do
            log "ERROR: Preprocessing FAILED for $sample — see $LOG_DIR/${sample}_preprocess.log"
        done < "$PREPROC_FAIL_FILE"
        err "$(wc -l < "$PREPROC_FAIL_FILE") sample(s) failed preprocessing."
    fi
    rm -f "$PREPROC_FAIL_FILE"
else
    log "No samples need preprocessing."
fi

run_mutect2_pair() {
    local job="$1"
    local tumor="${job%%:*}"
    local normal="${job##*:}"

    local tumor_bam="$OUT_DIR/$tumor/${tumor}_recal.bam"
    local normal_bam="$OUT_DIR/$normal/${normal}_recal.bam"

    local pdir="$OUT_DIR/${tumor}_vs_${normal}"
    local label="${tumor}_vs_${normal}"
    local mlog="$LOG_DIR/${label}_mutect2.log"
    local raw_vcf="$pdir/${label}.unfiltered.vcf.gz"

    local tumor_sm=""
    local normal_sm=""
    local -a args=()

    mkdir -p "$pdir"
    : > "$mlog"

    [[ -f "$tumor_bam" ]] || { echo "[$(date '+%F %T')] ERROR: Tumor BAM not found: $tumor_bam" | tee -a "$MAIN_LOG" >&2; return 1; }
    [[ -f "$normal_bam" ]] || { echo "[$(date '+%F %T')] ERROR: Normal BAM not found: $normal_bam" | tee -a "$MAIN_LOG" >&2; return 1; }

    if [[ -f "$raw_vcf" && -f "${raw_vcf}.tbi" ]]; then
        log "Mutect2 already done: $tumor vs $normal"
        return 0
    fi

    tumor_sm="$(get_bam_sm "$tumor_bam")" || {
        echo "[$(date '+%F %T')] ERROR: Could not determine unique tumor SM from $tumor_bam" | tee -a "$MAIN_LOG" >&2
        return 1
    }
    normal_sm="$(get_bam_sm "$normal_bam")" || {
        echo "[$(date '+%F %T')] ERROR: Could not determine unique normal SM from $normal_bam" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    [[ "$tumor_sm" != "$normal_sm" ]] || {
        echo "[$(date '+%F %T')] ERROR: Tumor and normal have identical SM ($tumor_sm)" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    log "Mutect2 PAIRED: $tumor vs $normal"

    args=(
        -R "$REF_FA"
        -I "$tumor_bam"  --tumor-sample "$tumor_sm"
        -I "$normal_bam" --normal-sample "$normal_sm"
        -O "$raw_vcf"
        --f1r2-tar-gz "$pdir/${label}_f1r2.tar.gz"
    )

    if [[ -n "${GNOMAD_VCF:-}" && -f "$GNOMAD_VCF" ]]; then
        args+=( --germline-resource "$GNOMAD_VCF" )
    fi

    if $USE_PON && [[ -n "${PON_VCF:-}" && -f "$PON_VCF" ]]; then
        args+=( --panel-of-normals "$PON_VCF" )
    fi

    "$GATK_BIN" Mutect2 "${args[@]}" >> "$mlog" 2>&1 || {
        echo "[$(date '+%F %T')] ERROR: Mutect2 FAILED for $tumor vs $normal — $mlog" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    log "Mutect2 done: $tumor vs $normal"
}

run_mutect2_tumor_only() {
    local tumor="$1"
    local tumor_bam="$OUT_DIR/$tumor/${tumor}_recal.bam"

    local pdir="$OUT_DIR/${tumor}_tumoronly"
    local mlog="$LOG_DIR/${tumor}_tumoronly_mutect2.log"
    local raw_vcf="$pdir/${tumor}.unfiltered.vcf.gz"

    local tumor_sm=""
    local -a args=()

    mkdir -p "$pdir"
    : > "$mlog"

    [[ -f "$tumor_bam" ]] || { echo "[$(date '+%F %T')] ERROR: Tumor BAM not found: $tumor_bam" | tee -a "$MAIN_LOG" >&2; return 1; }

    if [[ -f "$raw_vcf" && -f "${raw_vcf}.tbi" ]]; then
        log "Mutect2 already done: $tumor [tumor-only]"
        return 0
    fi

    tumor_sm="$(get_bam_sm "$tumor_bam")" || {
        echo "[$(date '+%F %T')] ERROR: Could not determine unique tumor SM from $tumor_bam" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    log "Mutect2 TUMOR-ONLY: $tumor"

    args=(
        -R "$REF_FA"
        -I "$tumor_bam"
        --tumor-sample "$tumor_sm"
        -O "$raw_vcf"
        --f1r2-tar-gz "$pdir/${tumor}_f1r2.tar.gz"
    )

    if [[ -n "${GNOMAD_VCF:-}" && -f "$GNOMAD_VCF" ]]; then
        args+=( --germline-resource "$GNOMAD_VCF" )
    fi

    "$GATK_BIN" Mutect2 "${args[@]}" >> "$mlog" 2>&1 || {
        echo "[$(date '+%F %T')] ERROR: Mutect2 FAILED for $tumor [tumor-only] — $mlog" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    log "Mutect2 done: $tumor [tumor-only]"
}

log "--- STEP 2: Somatic Variant Calling (${#PAIRED_JOBS[@]} paired + ${#TUMOR_ONLY_JOBS[@]} tumor-only) ---"

CALL_FAIL_FILE="$OUT_DIR/call_failures.txt"
: > "$CALL_FAIL_FILE"

if (( ${#PAIRED_JOBS[@]} > 0 )); then
    run_parallel_items run_mutect2_pair "$CALL_FAIL_FILE" "${PAIRED_JOBS[@]}"
fi

if (( ${#TUMOR_ONLY_JOBS[@]} > 0 )); then
    tmp_fail="$OUT_DIR/call_failures_tumoronly.txt"
    run_parallel_items run_mutect2_tumor_only "$tmp_fail" "${TUMOR_ONLY_JOBS[@]}"
    [[ -s "$tmp_fail" ]] && cat "$tmp_fail" >> "$CALL_FAIL_FILE"
    rm -f "$tmp_fail"
fi

if [[ -s "$CALL_FAIL_FILE" ]]; then
    while IFS= read -r item; do
        log "ERROR: Variant calling FAILED for $item"
    done < "$CALL_FAIL_FILE"
    err "$(wc -l < "$CALL_FAIL_FILE") Mutect2 job(s) failed."
fi
rm -f "$CALL_FAIL_FILE"

filter_pair() {
    local job="$1"
    local tumor="${job%%:*}"
    local normal="${job##*:}"

    local pdir="$OUT_DIR/${tumor}_vs_${normal}"
    local label="${tumor}_vs_${normal}"
    local flog="$LOG_DIR/${label}_filter.log"

    local raw_vcf="$pdir/${label}.unfiltered.vcf.gz"
    local stats_file="${raw_vcf}.stats"
    local f1r2_gz="$pdir/${label}_f1r2.tar.gz"
    local rom_gz="$pdir/${label}_read-orientation-model.tar.gz"
    local filtered_vcf="$pdir/${label}.filtered.vcf.gz"

    local pileup_tumor="$pdir/${tumor}_pileup.table"
    local pileup_normal="$pdir/${normal}_pileup.table"
    local contam_table="$pdir/${tumor}_contamination.table"
    local seg_table="$pdir/${tumor}_segments.table"

    local -a filter_args=()

    : > "$flog"

    [[ -f "$raw_vcf" ]] || { echo "[$(date '+%F %T')] WARN: Raw VCF missing for $tumor vs $normal — skipping filter" | tee -a "$MAIN_LOG"; return 0; }

    if [[ -f "$filtered_vcf" && -f "${filtered_vcf}.tbi" ]]; then
        log "Filtering already done: $tumor vs $normal"
        return 0
    fi

    if [[ -f "$f1r2_gz" && ! -f "$rom_gz" ]]; then
        "$GATK_BIN" LearnReadOrientationModel \
            -I "$f1r2_gz" \
            -O "$rom_gz" \
            >> "$flog" 2>&1 || warn "LearnReadOrientationModel failed for $tumor vs $normal"
    fi

    if [[ -n "${GNOMAD_VCF:-}" && -f "$GNOMAD_VCF" ]]; then
        if [[ ! -f "$pileup_tumor" ]]; then
            "$GATK_BIN" GetPileupSummaries \
                -I "$OUT_DIR/$tumor/${tumor}_recal.bam" \
                -V "$GNOMAD_VCF" \
                -L "$GNOMAD_VCF" \
                -O "$pileup_tumor" \
                >> "$flog" 2>&1 || true
        fi

        if [[ ! -f "$pileup_normal" ]]; then
            "$GATK_BIN" GetPileupSummaries \
                -I "$OUT_DIR/$normal/${normal}_recal.bam" \
                -V "$GNOMAD_VCF" \
                -L "$GNOMAD_VCF" \
                -O "$pileup_normal" \
                >> "$flog" 2>&1 || true
        fi

        if [[ -f "$pileup_tumor" && -f "$pileup_normal" && ! -f "$contam_table" ]]; then
            "$GATK_BIN" CalculateContamination \
                -I "$pileup_tumor" \
                -matched "$pileup_normal" \
                -O "$contam_table" \
                --tumor-segmentation "$seg_table" \
                >> "$flog" 2>&1 || true
        fi
    fi

    filter_args=(
        -R "$REF_FA"
        -V "$raw_vcf"
        --stats "$stats_file"
        -O "$filtered_vcf"
    )

    [[ -f "$rom_gz" ]] && filter_args+=( --ob-priors "$rom_gz" )
    [[ -f "$contam_table" ]] && filter_args+=( --contamination-table "$contam_table" )
    [[ -f "$seg_table" ]] && filter_args+=( --tumor-segmentation "$seg_table" )

    "$GATK_BIN" FilterMutectCalls "${filter_args[@]}" >> "$flog" 2>&1 || {
        echo "[$(date '+%F %T')] ERROR: FilterMutectCalls failed for $tumor vs $normal — $flog" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    log "Filtering done: $tumor vs $normal"
}

filter_tumor_only() {
    local tumor="$1"

    local pdir="$OUT_DIR/${tumor}_tumoronly"
    local flog="$LOG_DIR/${tumor}_tumoronly_filter.log"

    local raw_vcf="$pdir/${tumor}.unfiltered.vcf.gz"
    local stats_file="${raw_vcf}.stats"
    local f1r2_gz="$pdir/${tumor}_f1r2.tar.gz"
    local rom_gz="$pdir/${tumor}_read-orientation-model.tar.gz"
    local filtered_vcf="$pdir/${tumor}.filtered.vcf.gz"

    local pileup_tumor="$pdir/${tumor}_pileup.table"
    local contam_table="$pdir/${tumor}_contamination.table"
    local seg_table="$pdir/${tumor}_segments.table"

    local -a filter_args=()

    : > "$flog"

    [[ -f "$raw_vcf" ]] || { echo "[$(date '+%F %T')] WARN: Raw VCF missing for $tumor [tumor-only] — skipping filter" | tee -a "$MAIN_LOG"; return 0; }

    if [[ -f "$filtered_vcf" && -f "${filtered_vcf}.tbi" ]]; then
        log "Filtering already done: $tumor [tumor-only]"
        return 0
    fi

    if [[ -f "$f1r2_gz" && ! -f "$rom_gz" ]]; then
        "$GATK_BIN" LearnReadOrientationModel \
            -I "$f1r2_gz" \
            -O "$rom_gz" \
            >> "$flog" 2>&1 || warn "LearnReadOrientationModel failed for $tumor [tumor-only]"
    fi

    if [[ -n "${GNOMAD_VCF:-}" && -f "$GNOMAD_VCF" ]]; then
        if [[ ! -f "$pileup_tumor" ]]; then
            "$GATK_BIN" GetPileupSummaries \
                -I "$OUT_DIR/$tumor/${tumor}_recal.bam" \
                -V "$GNOMAD_VCF" \
                -L "$GNOMAD_VCF" \
                -O "$pileup_tumor" \
                >> "$flog" 2>&1 || true
        fi

        if [[ -f "$pileup_tumor" && ! -f "$contam_table" ]]; then
            "$GATK_BIN" CalculateContamination \
                -I "$pileup_tumor" \
                -O "$contam_table" \
                --tumor-segmentation "$seg_table" \
                >> "$flog" 2>&1 || true
        fi
    fi

    filter_args=(
        -R "$REF_FA"
        -V "$raw_vcf"
        --stats "$stats_file"
        -O "$filtered_vcf"
    )

    [[ -f "$rom_gz" ]] && filter_args+=( --ob-priors "$rom_gz" )
    [[ -f "$contam_table" ]] && filter_args+=( --contamination-table "$contam_table" )
    [[ -f "$seg_table" ]] && filter_args+=( --tumor-segmentation "$seg_table" )

    "$GATK_BIN" FilterMutectCalls "${filter_args[@]}" >> "$flog" 2>&1 || {
        echo "[$(date '+%F %T')] ERROR: FilterMutectCalls failed for $tumor [tumor-only] — $flog" | tee -a "$MAIN_LOG" >&2
        return 1
    }

    log "Filtering done: $tumor [tumor-only]"
}

log "--- STEP 3: Filtering Mutect2 calls ---"

FILTER_FAIL_FILE="$OUT_DIR/filter_failures.txt"
: > "$FILTER_FAIL_FILE"

if (( ${#PAIRED_JOBS[@]} > 0 )); then
    run_parallel_items filter_pair "$FILTER_FAIL_FILE" "${PAIRED_JOBS[@]}"
fi

if (( ${#TUMOR_ONLY_JOBS[@]} > 0 )); then
    tmp_fail="$OUT_DIR/filter_failures_tumoronly.txt"
    run_parallel_items filter_tumor_only "$tmp_fail" "${TUMOR_ONLY_JOBS[@]}"
    [[ -s "$tmp_fail" ]] && cat "$tmp_fail" >> "$FILTER_FAIL_FILE"
    rm -f "$tmp_fail"
fi

if [[ -s "$FILTER_FAIL_FILE" ]]; then
    while IFS= read -r item; do
        log "ERROR: FilterMutectCalls FAILED for $item"
    done < "$FILTER_FAIL_FILE"
    err "$(wc -l < "$FILTER_FAIL_FILE") filter job(s) failed."
fi
rm -f "$FILTER_FAIL_FILE"

log "==========================================="
log "  Pipeline complete."
log "  Output : $OUT_DIR"
log "  Log    : $MAIN_LOG"
log "==========================================="
