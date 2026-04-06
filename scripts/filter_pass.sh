#!/bin/bash
# scripts/filter_pass.sh — Extract PASS-only variants from a Mutect2 filtered VCF

set -euo pipefail
[[ $# -eq 2 ]] || { echo "Usage: filter_pass.sh <input.vcf.gz> <output.pass.vcf.gz>"; exit 1; }

INFILE="$1"; OUTFILE="$2"
command -v bcftools &>/dev/null || { echo "ERROR: bcftools not found"; exit 1; }

bcftools view -f PASS "$INFILE" -Oz -o "$OUTFILE"
tabix -p vcf "$OUTFILE"
echo "PASS-only VCF written to: $OUTFILE"
echo "Index written to: ${OUTFILE}.tbi"
echo "PASS variant count: $(bcftools view -H "$OUTFILE" | wc -l)"
