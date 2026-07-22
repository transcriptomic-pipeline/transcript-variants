# Transcript Variants Module — Transcriptomic Pipeline

Somatic variant calling from RNA-seq BAM files of the tumor samples using GATK Mutect2. Includes GATK best-practice preprocessing (duplicate marking, SplitNCigarReads, BQSR), Panel of Normals (PoN) construction, gnomAD population filtering, and PASS variant extraction.

**Author:** Babul Pradhan, Ph.D. Student, Institute of Bioinformatics (IOB), Bangalore, India  
**Date:** 09/02/2025  
**Part of:** Automated Transcriptomic Analysis Pipeline

---

## Overview

This module identifies somatic single nucleotide variants (SNVs) and small insertions/deletions (indels) from RNA-seq data of the tumors, following GATK best practices adapted for RNA-seq variant calling.

**Pipeline position:** `Sorted BAM → [Transcript Variants] → Filtered Somatic VCF`

---

## Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| GATK | ≥ 4.3 | Mutect2 somatic variant calling, BQSR, FilterMutectCalls |
| Picard | ≥ 2.27 | Duplicate marking (MarkDuplicates) |
| SAMtools | ≥ 1.15 | BAM manipulation and indexing |

---

## Repository Structure

transcript-variants/
├── run_transcript_variants.sh # Main orchestration script
├── install.sh # Tool installation script
├── config/ # Tool path configuration
├── scripts/
│ ├── gatk_preprocess.sh # Duplicate marking, SplitNCigarReads, BQSR
│ ├── mutect2_calling.sh # Somatic variant calling with Mutect2
│ ├── build_pon.sh # Panel of Normals construction
│ ├── download_gnomad.sh # gnomAD allele frequency VCF download
│ └── filter_pass.sh # PASS variant filtering post-FilterMutectCalls
├── README.md
└── LICENSE


---

## Pipeline Steps

Preprocessing (gatk_preprocess.sh)
└── MarkDuplicates (Picard) → SplitNCigarReads → BQSR → Recalibrated BAM

Panel of Normals (build_pon.sh)
└── Mutect2 (tumor-only on normals) → GenomicsDBImport → CreateSomaticPanelOfNormals → PoN VCF

gnomAD Download (download_gnomad.sh)
└── Download + index gnomAD v3.1.2 VCF

Variant Calling (mutect2_calling.sh)
└── Mutect2 (tumor-only + PoN + gnomAD) → GetPileupSummaries → CalculateContamination → FilterMutectCalls

PASS Filtering (filter_pass.sh)
└── bcftools/GATK filter → PASS-only somatic VCF


---

## Installation

```bash
bash install.sh
```

---

## Usage

```bash
bash run_transcript_variants.sh -i <bam_dir> -o <output_dir> -r <reference.fa> [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `-i, --input DIR` | Directory containing preprocessed BAM files |
| `-o, --output DIR` | Output directory for VCF results |
| `-r, --reference FASTA` | Reference genome FASTA (GRCh38) |

### Key Optional Arguments

| Argument | Description |
|----------|-------------|
| `--pon VCF` | Pre-built Panel of Normals VCF |
| `--gnomad VCF` | gnomAD allele frequency VCF |
| `--normal-bam-dir DIR` | Normal BAM directory for PoN construction |
| `-t, --threads NUM` | CPU threads (default: 8) |

### Example

```bash
# Full pipeline with PoN and gnomAD
bash run_transcript_variants.sh \
  -i alignment_results/ \
  -o variant_results/ \
  -r /ref/GRCh38.fa \
  --pon pon.vcf.gz \
  --gnomad gnomad.genomes.v3.1.2.sites.vcf.gz

# Build PoN from normal samples
bash scripts/build_pon.sh \
  -i normal_bams/ \
  -o pon_output/ \
  -r /ref/GRCh38.fa
```

---

## Output Structure

variant_results/
├── preprocessed/
│ └── [sample]_recal.bam # BQSR-recalibrated BAM
├── vcf/
│ ├── [sample]_mutect2.vcf.gz # Raw Mutect2 calls
│ ├── [sample]_filtered.vcf.gz # After FilterMutectCalls
│ └── [sample]_PASS.vcf.gz # Final PASS-only somatic variants
├── contamination/
│ └── [sample]_contamination.table
└── logs/


---

## Notes

- gnomAD v3.1.2 (GRCh38) is used for population allele frequency filtering
- PoN reduces germline and technical artefact calls
- Designed for **tumor-only** calling (no matched normal RNA-seq required)
