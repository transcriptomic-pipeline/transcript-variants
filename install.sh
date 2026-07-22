#!/bin/bash

#Author: Babul Pradhan, Ph.D. Student, Institute of Bioinformatics (IOB), Bangalore, India
#Date: 09/02/2025
#Purpose: Installation script for setting up GATK, SAMtools, and Picard tools required for the transcript variant calling module of the transcriptomic analysis pipeline

# install.sh — Install tools and references for transcript-variants pipeline.
# Downloads dbSNP 146 (hg38) and Mills indels from Broad FTP bundle.
# Optionally downloads gnomAD v4.1.
# Supports --chr-style chr|nochr (default: nochr, matches Ensembl GRCh38).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="$SCRIPT_DIR/install.log"
log()  { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
warn() { echo "[$(date '+%F %T')] WARN: $*"  | tee -a "$LOG" >&2; }
err()  { echo "[$(date '+%F %T')] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

# ─── Defaults ────────────────────────────────────────────────
TOOLS_DIR="${TOOLS_DIR:-$HOME/tools}"
REF_DIR="${REF_DIR:-$HOME/references/GRCh38}"
GATK_VERSION="4.4.0.0"
GATK_DIR="$TOOLS_DIR/gatk-${GATK_VERSION}"
PICARD_DIR="$TOOLS_DIR/picard"
SAMTOOLS_DIR="$TOOLS_DIR/samtools"
BCFTOOLS_DIR="$TOOLS_DIR/bcftools"

INSTALL_REF=true
INSTALL_GNOMAD=true
GNOMAD_THREADS=8
CHR_STYLE="nochr"   # nochr (Ensembl, default) | chr (UCSC)
NO_APT=false

# Broad Institute public FTP — anonymous, no credentials needed
# Files are pre-bgzipped .vcf.gz — NO bgzip step required
BROAD_FTP="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38"

print_usage() {
  cat <<EOF
Usage: install.sh [OPTIONS]

  --chr-style STR      'nochr' (Ensembl, default) or 'chr' (UCSC/hg38)
                       nochr: strip 'chr' prefix — matches Ensembl GRCh38
                       chr:   keep 'chr' prefix  — matches UCSC hg38
  --no-ref             Skip reference FASTA + known-sites download
  --no-gnomad          Skip gnomAD v4.1 download/merge
  --gnomad-threads N   Threads for gnomAD concat/index (default: 8)
  --tools-dir DIR      Install tools here (default: ~/tools)
  --ref-dir DIR        Install references here (default: ~/references/GRCh38)
  --no-apt             Skip apt-get installs
  -h, --help

Source URLs (Broad FTP, anonymous):
  dbSNP : ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
  Mills : ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr-style)       CHR_STYLE="$2"; shift 2 ;;
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

[[ "$CHR_STYLE" == "nochr" || "$CHR_STYLE" == "chr" ]] \
  || err "--chr-style must be 'nochr' or 'chr' (got: $CHR_STYLE)"

mkdir -p "$TOOLS_DIR" "$REF_DIR"
log "========================================="
log "  transcript-variants :: install.sh"
log "========================================="
log "  Tools dir     : $TOOLS_DIR"
log "  Ref dir       : $REF_DIR"
log "  Chr style     : $CHR_STYLE"
log "  Install gnomAD: $INSTALL_GNOMAD"

# ─── Helper: safe wget with visible progress ──────────────────
# wget -c enables resume on re-run; works natively for FTP + HTTP
safe_wget() {
  local url="$1" dest="$2"
  log "  Downloading: $(basename "$dest")"
  log "  Source: $url"
  wget -c --progress=dot:giga -O "$dest" "$url" 2>&1 | tee -a "$LOG" || {
    rm -f "$dest"
    err "Download failed: $url  (see $LOG for details)"
  }
  [[ -s "$dest" ]] || { rm -f "$dest"; err "Downloaded file is empty: $dest"; }
  log "  Done: $dest ($(du -sh "$dest" | cut -f1))"
}

# ─── Helper: build chr<->nochr rename map ─────────────────────
# Broad FTP: UCSC chr-prefix (chr1..chr22, chrX, chrY, chrM)
# Ensembl:   no prefix        (1..22, X, Y, MT)
build_rename_map() {
  local direction="$1" outfile="$2"
  : > "$outfile"
  for c in {1..22} X Y; do
    if [[ "$direction" == "chr_to_nochr" ]]; then
      printf "chr%s\t%s\n" "$c" "$c" >> "$outfile"
    else
      printf "%s\tchr%s\n" "$c" "$c" >> "$outfile"
    fi
  done
  if [[ "$direction" == "chr_to_nochr" ]]; then
    printf "chrM\tMT\n"  >> "$outfile"
    printf "chrMT\tMT\n" >> "$outfile"
  else
    printf "MT\tchrM\n"  >> "$outfile"
  fi
}

# ─── Helper: convert downloaded VCF to target chr style ───────
# $1 = input .vcf.gz  (chr-style, from Broad FTP)
# $2 = output .vcf.gz (target chr style)
# $3 = rename_map file  OR  "" to skip renaming (keep chr-style)
convert_vcf_chrstyle() {
  local src="$1" dst="$2" rename_map="${3:-}"
  if [[ -n "$rename_map" && -f "$rename_map" ]]; then
    log "  Converting chr -> nochr: $(basename "$src") -> $(basename "$dst")"
    bcftools annotate \
      --rename-chrs "$rename_map" \
      --threads 4 \
      -Oz -o "$dst" \
      "$src" 2>>"$LOG" \
      || err "bcftools --rename-chrs failed on $(basename "$src") — see $LOG"
    tabix -p vcf "$dst" 2>>"$LOG" \
      || err "tabix failed on $(basename "$dst")"
    rm -f "$src" "${src}.tbi"
  else
    mv "$src"       "$dst"
    mv "${src}.tbi" "${dst}.tbi" 2>/dev/null \
      || tabix -p vcf "$dst" 2>>"$LOG"
  fi
  log "  Ready (${CHR_STYLE}): $dst [OK]"
}

# ============================================================
# System dependencies
# ============================================================
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

# ============================================================
# samtools / bcftools / GATK / Picard
# ============================================================
if command -v samtools &>/dev/null; then
  log "samtools: $(command -v samtools) [OK]"
elif [[ -x "$SAMTOOLS_DIR/bin/samtools" ]]; then
  log "samtools: $SAMTOOLS_DIR [OK]"
  export PATH="$SAMTOOLS_DIR/bin:$PATH"
else
  log "Installing samtools 1.19.2..."
  cd "$TOOLS_DIR"
  safe_wget "https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2" \
            "samtools-1.19.2.tar.bz2"
  tar -xjf "samtools-1.19.2.tar.bz2"
  cd "samtools-1.19.2"
  ./configure --prefix="$SAMTOOLS_DIR" --without-curses --disable-lzma &>>"$LOG"
  make -j4 &>>"$LOG" && make install &>>"$LOG"
  export PATH="$SAMTOOLS_DIR/bin:$PATH"
  cd "$SCRIPT_DIR"
  log "samtools installed: $SAMTOOLS_DIR"
fi

if command -v bcftools &>/dev/null; then
  log "bcftools: $(command -v bcftools) [OK]"
elif [[ -x "$BCFTOOLS_DIR/bin/bcftools" ]]; then
  log "bcftools: $BCFTOOLS_DIR [OK]"
  export PATH="$BCFTOOLS_DIR/bin:$PATH"
else
  log "Installing bcftools 1.19..."
  cd "$TOOLS_DIR"
  safe_wget "https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2" \
            "bcftools-1.19.tar.bz2"
  tar -xjf "bcftools-1.19.tar.bz2"
  cd "bcftools-1.19"
  ./configure --prefix="$BCFTOOLS_DIR" &>>"$LOG"
  make -j4 &>>"$LOG" && make install &>>"$LOG"
  export PATH="$BCFTOOLS_DIR/bin:$PATH"
  cd "$SCRIPT_DIR"
  log "bcftools installed: $BCFTOOLS_DIR"
fi

if command -v gatk &>/dev/null; then
  log "GATK: $(command -v gatk) [OK]"
elif [[ -f "$GATK_DIR/gatk" ]]; then
  log "GATK: $GATK_DIR/gatk [OK]"
else
  log "Downloading GATK ${GATK_VERSION}..."
  safe_wget \
    "https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip" \
    "$TOOLS_DIR/gatk-${GATK_VERSION}.zip"
  unzip -q "$TOOLS_DIR/gatk-${GATK_VERSION}.zip" -d "$TOOLS_DIR"
  rm -f "$TOOLS_DIR/gatk-${GATK_VERSION}.zip"
  log "GATK installed: $GATK_DIR"
fi

if command -v picard &>/dev/null; then
  log "picard: $(command -v picard) [OK]"
elif [[ -f "$PICARD_DIR/picard.jar" ]]; then
  log "picard: $PICARD_DIR/picard.jar [OK]"
else
  log "Downloading Picard..."
  mkdir -p "$PICARD_DIR"
  safe_wget "https://github.com/broadinstitute/picard/releases/latest/download/picard.jar" \
            "$PICARD_DIR/picard.jar"
  log "Picard installed: $PICARD_DIR/picard.jar"
fi

# ============================================================
# Reference resources
# ============================================================
DBSNP_VCF=""; MILLS_VCF=""; GNOMAD_VCF=""

if $INSTALL_REF; then

  ENSEMBL_REL="113"
  REF_FA="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

  # Filenames encode chr style — both can coexist in the same ref dir
  if [[ "$CHR_STYLE" == "nochr" ]]; then
    DBSNP_VCF="$REF_DIR/dbsnp_146.hg38.nochr.vcf.gz"
    MILLS_VCF="$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.nochr.vcf.gz"
  else
    DBSNP_VCF="$REF_DIR/dbsnp_146.hg38.chr.vcf.gz"
    MILLS_VCF="$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.chr.vcf.gz"
  fi
  GNOMAD_VCF="$REF_DIR/gnomad.genomes.v4.1.sites.merged.${CHR_STYLE}.vcf.gz"

  # ─── Reference FASTA ─────────────────────────────────────
  if [[ -f "$REF_FA" ]]; then
    log "Reference FASTA: $REF_FA [OK]"
  else
    log "Downloading GRCh38 primary assembly (Ensembl r${ENSEMBL_REL})..."
    safe_wget \
      "https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
      "$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    log "Decompressing reference..."
    gunzip "$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    log "Indexing (fai)..."
    samtools faidx "$REF_FA"
    GATK_BIN_TMP="${GATK_DIR}/gatk"
    [[ -f "$GATK_BIN_TMP" ]] || GATK_BIN_TMP="$(command -v gatk 2>/dev/null || true)"
    [[ -x "$GATK_BIN_TMP" ]] && {
      log "Creating sequence dictionary..."
      "$GATK_BIN_TMP" CreateSequenceDictionary -R "$REF_FA" &>>"$LOG"
    } || warn "Run manually: gatk CreateSequenceDictionary -R $REF_FA"
    log "Reference ready: $REF_FA"
  fi
  [[ -f "${REF_FA}.fai" ]] || { log "Indexing reference..."; samtools faidx "$REF_FA"; }

  # ─── dbSNP 146 ───────────────────────────────────────────
  # Source: Broad FTP, pre-bgzipped .vcf.gz, chr-prefix (UCSC)
  if [[ -f "$DBSNP_VCF" && -f "${DBSNP_VCF}.tbi" ]]; then
    log "dbSNP: $DBSNP_VCF [OK]"
  else
    DBSNP_RAW="$REF_DIR/dbsnp_146.hg38.raw.chr.vcf.gz"
    if [[ ! -f "$DBSNP_RAW" || ! -f "${DBSNP_RAW}.tbi" ]]; then
      log "Downloading dbSNP 146 hg38 from Broad FTP..."
      safe_wget "${BROAD_FTP}/dbsnp_146.hg38.vcf.gz"     "$DBSNP_RAW"
      safe_wget "${BROAD_FTP}/dbsnp_146.hg38.vcf.gz.tbi" "${DBSNP_RAW}.tbi"
    else
      log "dbSNP raw (chr-style) already present: $DBSNP_RAW"
    fi
    if [[ "$CHR_STYLE" == "nochr" ]]; then
      RENAME_MAP="$REF_DIR/chr_to_nochr.map"
      build_rename_map "chr_to_nochr" "$RENAME_MAP"
      convert_vcf_chrstyle "$DBSNP_RAW" "$DBSNP_VCF" "$RENAME_MAP"
      rm -f "$RENAME_MAP"
    else
      convert_vcf_chrstyle "$DBSNP_RAW" "$DBSNP_VCF" ""
    fi
  fi

  # ─── Mills and 1000G indels ───────────────────────────────
  # Source: Broad FTP, pre-bgzipped .vcf.gz, chr-prefix (UCSC)
  if [[ -f "$MILLS_VCF" && -f "${MILLS_VCF}.tbi" ]]; then
    log "Mills indels: $MILLS_VCF [OK]"
  else
    MILLS_RAW="$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.raw.chr.vcf.gz"
    if [[ ! -f "$MILLS_RAW" || ! -f "${MILLS_RAW}.tbi" ]]; then
      log "Downloading Mills and 1000G indels (hg38) from Broad FTP..."
      safe_wget \
        "${BROAD_FTP}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        "$MILLS_RAW"
      safe_wget \
        "${BROAD_FTP}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
        "${MILLS_RAW}.tbi"
    else
      log "Mills raw (chr-style) already present: $MILLS_RAW"
    fi
    if [[ "$CHR_STYLE" == "nochr" ]]; then
      RENAME_MAP="$REF_DIR/chr_to_nochr.map"
      [[ -f "$RENAME_MAP" ]] || build_rename_map "chr_to_nochr" "$RENAME_MAP"
      convert_vcf_chrstyle "$MILLS_RAW" "$MILLS_VCF" "$RENAME_MAP"
      rm -f "$RENAME_MAP"
    else
      convert_vcf_chrstyle "$MILLS_RAW" "$MILLS_VCF" ""
    fi
  fi

  # ─── gnomAD v4.1 ─────────────────────────────────────────
  if $INSTALL_GNOMAD; then
    if [[ -f "$GNOMAD_VCF" && -f "${GNOMAD_VCF}.tbi" ]]; then
      log "gnomAD v4.1 (${CHR_STYLE}): $GNOMAD_VCF [OK]"
    else
      command -v gsutil &>/dev/null \
        || err "gnomAD requires gsutil. Install: https://cloud.google.com/sdk/docs/install  Or use --no-gnomad."
      log "Downloading and merging gnomAD v4.1 (chr-style=${CHR_STYLE}) ..."
      GNOMAD_VCF=$(bash "$SCRIPT_DIR/scripts/download_gnomad.sh" \
        --out-dir   "$REF_DIR" \
        --threads   "$GNOMAD_THREADS" \
        --chr-style "$CHR_STYLE") \
        || err "gnomAD download/merge failed"
      log "gnomAD ready: $GNOMAD_VCF"
    fi
  else
    log "Skipping gnomAD (--no-gnomad). Provide --gnomad PATH at runtime if needed."
    GNOMAD_VCF=""
  fi

else
  log "Skipping reference installation (--no-ref)."
fi

# ============================================================
# Write .pipeline_env
# ============================================================
GATK_BIN_FINAL="${GATK_DIR}/gatk"
[[ -f "$GATK_BIN_FINAL" ]] || GATK_BIN_FINAL="$(command -v gatk 2>/dev/null || echo '')"

ENV_FILE="$SCRIPT_DIR/.pipeline_env"
cat >"$ENV_FILE" <<ENVEOF
# Auto-generated by install.sh — re-run to regenerate
GATK_BIN="${GATK_BIN_FINAL}"
PICARD_JAR="${PICARD_DIR}/picard.jar"
DEFAULT_REF="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DEFAULT_DBSNP="${DBSNP_VCF:-}"
DEFAULT_MILLS="${MILLS_VCF:-}"
DEFAULT_GNOMAD="${GNOMAD_VCF:-}"
CHR_STYLE="${CHR_STYLE}"
ENVEOF

log "Environment config written: $ENV_FILE"
log "========================================="
log "  Installation complete."
log "  Chr style : $CHR_STYLE"
log "  dbSNP     : ${DBSNP_VCF:-skipped (--no-ref)}"
log "  Mills     : ${MILLS_VCF:-skipped (--no-ref)}"
log "  gnomAD    : ${GNOMAD_VCF:-(deferred — provide --gnomad at runtime)}"
log "========================================="
