#!/bin/bash
# install.sh — Install tools and reference resources for transcript-variants pipeline.
# Skips any tool/resource already installed.
# gnomAD v4.1 is downloaded and merged by default unless --no-gnomad is passed.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="$SCRIPT_DIR/install.log"
log()  { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
warn() { echo "[$(date '+%F %T')] WARN: $*" | tee -a "$LOG" >&2; }
err()  { echo "[$(date '+%F %T')] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

# ---------- Defaults ----------
TOOLS_DIR="${TOOLS_DIR:-$HOME/tools}"
REF_DIR="${REF_DIR:-$HOME/references/GRCh38}"
GATK_VERSION="4.4.0.0"
GATK_DIR="$TOOLS_DIR/gatk-${GATK_VERSION}"
PICARD_DIR="$TOOLS_DIR/picard"
SAMTOOLS_DIR="$TOOLS_DIR/samtools"
BCFTOOLS_DIR="$TOOLS_DIR/bcftools"

INSTALL_REF=true
INSTALL_GNOMAD=true      # download+merge gnomAD v4.1 by default
GNOMAD_THREADS=8
NO_APT=false

print_usage() {
  cat <<EOF
Usage: install.sh [OPTIONS]
  --no-ref           Skip reference FASTA + known-sites download
  --no-gnomad        Skip gnomAD v4.1 download/merge
  --gnomad-threads N Threads for gnomAD merge/index (default: 8)
  --tools-dir DIR    Install tools here (default: ~/tools)
  --ref-dir DIR      Install references here (default: ~/references/GRCh38)
  --no-apt           Skip apt-get (assume system dependencies present)
  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-ref)          INSTALL_REF=false; shift ;;
    --no-gnomad)       INSTALL_GNOMAD=false; shift ;;
    --gnomad-threads)  GNOMAD_THREADS="$2"; shift 2 ;;
    --tools-dir)
      TOOLS_DIR="$2"
      GATK_DIR="$TOOLS_DIR/gatk-${GATK_VERSION}"
      PICARD_DIR="$TOOLS_DIR/picard"
      SAMTOOLS_DIR="$TOOLS_DIR/samtools"
      BCFTOOLS_DIR="$TOOLS_DIR/bcftools"
      shift 2 ;;
    --ref-dir)         REF_DIR="$2"; shift 2 ;;
    --no-apt)          NO_APT=true; shift ;;
    -h|--help)         print_usage; exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

mkdir -p "$TOOLS_DIR" "$REF_DIR"
log "========================================="
log "  transcript-variants :: install.sh"
log "========================================="
log "Tools dir : $TOOLS_DIR"
log "Ref dir   : $REF_DIR"
log "Install gnomAD: $INSTALL_GNOMAD"

# ---------- System dependencies ----------
if ! $NO_APT && command -v apt-get &>/dev/null; then
  log "Checking system packages..."
  PKGS=()
  for p in wget curl unzip tabix openjdk-17-jre-headless; do
    dpkg -s "$p" &>/dev/null || PKGS+=( "$p" )
  done
  if (( ${#PKGS[@]} )); then
    log "Installing apt packages: ${PKGS[*]}"
    sudo apt-get update -qq && sudo apt-get install -y "${PKGS[@]}" &>>"$LOG"
  else
    log "System packages OK."
  fi
fi

# ---------- samtools ----------
if command -v samtools &>/dev/null; then
  log "samtools: $(command -v samtools) [OK]"
elif [[ -x "$SAMTOOLS_DIR/bin/samtools" ]]; then
  log "samtools: $SAMTOOLS_DIR [OK]"
  export PATH="$SAMTOOLS_DIR/bin:$PATH"
else
  log "Installing samtools..."
  SM_VER="1.19.2"
  cd "$TOOLS_DIR"
  wget -q "https://github.com/samtools/samtools/releases/download/${SM_VER}/samtools-${SM_VER}.tar.bz2"
  tar -xjf "samtools-${SM_VER}.tar.bz2"
  cd "samtools-${SM_VER}"
  ./configure --prefix="$SAMTOOLS_DIR" --without-curses --disable-lzma &>>"$LOG"
  make -j4 &>>"$LOG" && make install &>>"$LOG"
  export PATH="$SAMTOOLS_DIR/bin:$PATH"
  log "samtools installed: $SAMTOOLS_DIR"
fi

# ---------- bcftools ----------
if command -v bcftools &>/dev/null; then
  log "bcftools: $(command -v bcftools) [OK]"
elif [[ -x "$BCFTOOLS_DIR/bin/bcftools" ]]; then
  log "bcftools: $BCFTOOLS_DIR [OK]"
  export PATH="$BCFTOOLS_DIR/bin:$PATH"
else
  log "Installing bcftools..."
  BC_VER="1.19"
  cd "$TOOLS_DIR"
  wget -q "https://github.com/samtools/bcftools/releases/download/${BC_VER}/bcftools-${BC_VER}.tar.bz2"
  tar -xjf "bcftools-${BC_VER}.tar.bz2"
  cd "bcftools-${BC_VER}"
  ./configure --prefix="$BCFTOOLS_DIR" &>>"$LOG"
  make -j4 &>>"$LOG" && make install &>>"$LOG"
  export PATH="$BCFTOOLS_DIR/bin:$PATH"
  log "bcftools installed: $BCFTOOLS_DIR"
fi

# ---------- GATK ----------
if command -v gatk &>/dev/null; then
  log "GATK: $(command -v gatk) [OK]"
elif [[ -f "$GATK_DIR/gatk" ]]; then
  log "GATK: $GATK_DIR [OK]"
else
  log "Downloading GATK ${GATK_VERSION}..."
  wget -q -P "$TOOLS_DIR" \
    "https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip"
  unzip -q "$TOOLS_DIR/gatk-${GATK_VERSION}.zip" -d "$TOOLS_DIR"
  rm -f "$TOOLS_DIR/gatk-${GATK_VERSION}.zip"
  log "GATK installed: $GATK_DIR"
fi

# ---------- Picard ----------
if command -v picard &>/dev/null; then
  log "picard: $(command -v picard) [OK]"
elif [[ -f "$PICARD_DIR/picard.jar" ]]; then
  log "picard: $PICARD_DIR/picard.jar [OK]"
else
  log "Downloading Picard..."
  mkdir -p "$PICARD_DIR"
  wget -q -O "$PICARD_DIR/picard.jar" \
    "https://github.com/broadinstitute/picard/releases/latest/download/picard.jar"
  log "Picard installed: $PICARD_DIR/picard.jar"
fi

# ---------- Reference resources ----------
GNOMAD_VCF=""
if $INSTALL_REF; then
  ENSEMBL_REL="113"
  REF_FA="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  DBSNP_VCF="$REF_DIR/dbsnp_nochr.vcf.gz"
  MILLS_VCF="$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.nochr.vcf.gz"
  GNOMAD_VCF="$REF_DIR/gnomad.genomes.v4.1.sites.merged.vcf.gz"

  # Reference FASTA
  if [[ -f "$REF_FA" ]]; then
    log "Reference FASTA: $REF_FA [OK]"
  else
    log "Downloading GRCh38 primary assembly (Ensembl r${ENSEMBL_REL})..."
    wget -q -P "$REF_DIR" \
      "https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    log "Decompressing reference..."
    gunzip "$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    log "Indexing (fai + dict)..."
    samtools faidx "$REF_FA"
    GATK_BIN_LOCAL="${GATK_DIR}/gatk"
    [[ -f "$GATK_BIN_LOCAL" ]] || GATK_BIN_LOCAL="$(command -v gatk 2>/dev/null || true)"
    [[ -x "$GATK_BIN_LOCAL" ]] && "$GATK_BIN_LOCAL" CreateSequenceDictionary -R "$REF_FA" &>>"$LOG" \
      || warn "GATK not found for dict creation; run: gatk CreateSequenceDictionary -R $REF_FA"
    log "Reference ready: $REF_FA"
  fi
  [[ -f "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

  # dbSNP
  if [[ -f "$DBSNP_VCF" ]]; then
    log "dbSNP: $DBSNP_VCF [OK]"
  else
    warn "dbSNP not found: $DBSNP_VCF"
    warn "Download from GATK resource bundle and place at: $DBSNP_VCF"
  fi

  # Mills
  if [[ -f "$MILLS_VCF" ]]; then
    log "Mills indels: $MILLS_VCF [OK]"
  else
    warn "Mills VCF not found: $MILLS_VCF"
    warn "Download from GATK resource bundle and place at: $MILLS_VCF"
  fi

  # gnomAD v4.1
  if $INSTALL_GNOMAD; then
    if [[ -f "$GNOMAD_VCF" && -f "${GNOMAD_VCF}.tbi" ]]; then
      log "gnomAD v4.1 merged: $GNOMAD_VCF [OK]"
    else
      log "Downloading and merging gnomAD v4.1 (this will take a long time)..."
      GNOMAD_VCF=$(bash "$SCRIPT_DIR/scripts/download_gnomad.sh" \
        --out-dir "$REF_DIR" \
        --threads "$GNOMAD_THREADS")
      log "gnomAD ready: $GNOMAD_VCF"
    fi
  else
    log "Skipping gnomAD download (--no-gnomad). Provide --gnomad at runtime if needed."
  fi
fi

# ---------- Write environment config ----------
ENV_FILE="$SCRIPT_DIR/.pipeline_env"
GATK_BIN_FINAL="${GATK_DIR}/gatk"
[[ -f "$GATK_BIN_FINAL" ]] || GATK_BIN_FINAL="$(command -v gatk 2>/dev/null || echo '')"
PICARD_JAR_FINAL="${PICARD_DIR}/picard.jar"
[[ -f "$PICARD_JAR_FINAL" ]] || PICARD_JAR_FINAL="$(command -v picard 2>/dev/null || echo '')"

cat >"$ENV_FILE" <<ENVEOF
# Auto-generated by install.sh — sourced by run_transcript_variants.sh
GATK_BIN="${GATK_BIN_FINAL}"
PICARD_JAR="${PICARD_JAR_FINAL}"
DEFAULT_REF="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DEFAULT_DBSNP="${REF_DIR}/dbsnp_nochr.vcf.gz"
DEFAULT_MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.nochr.vcf.gz"
DEFAULT_GNOMAD="${GNOMAD_VCF:-${REF_DIR}/gnomad.genomes.v4.1.sites.merged.vcf.gz}"
ENVEOF

log "Environment config: $ENV_FILE"
log "========================================="
log "  Installation complete."
log "========================================="
