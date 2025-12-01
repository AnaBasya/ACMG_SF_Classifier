# ACMG Secondary Findings Classifier

## Overview
This repository contains a production-oriented implementation of an automated-first ACMG Secondary Findings (SF) classifier. 

### Goals
- **Maximise safe automation:** Automatically classify and report variants where evidence is robust and conforms to decision rules.
- **Minimise manual workload:** Use high-quality DB evidence (ClinVar, gnomAD) and conservative computational thresholds to push variants to automated conclusions where appropriate.
- **Preserve safety:** Gene-specific exceptions (e.g., *TTN* meta-exons, LOF mechanism validity) trigger specific logic or manual review.

## Contents
- `acmg_sf_classifier.py`
  - Main pipeline script. Performs VCF parsing, VEP annotation (optional), applies gene-specific rules and ACMG criteria, consults databases, and writes outputs.
- `scripts/`
  - Helper utilities (e.g., `make_ttn_meta_csv.py`, download scripts).
- `README.md`
  - This document.

## Requirements

### Python Environment
- Python 3.9+
- Required pip packages:
  - `pandas`
  - `pysam`
  - `cyvcf2`
  - `intervaltree` (required for TTN meta-exon support)

### System Tools
- **VEP (Variant Effect Predictor):** Required if running annotation within the pipeline.
- **tabix / bgzip:** Required for handling VCF/TSV indexes.

### Databases
Provide paths via CLI arguments or a JSON config (`--db-paths-json`).
- **ClinVar:** VCF (bgzipped & indexed). Recommended for PS1/PM5/Disease matching.
- **gnomAD:** VCF or TSV (bgzipped & indexed). Supports v2/v3/v4.
- **dbNSFP:** Used for in-silico scores (REVEL, SpliceAI, CADD, AlphaMissense).
- **ACMG Table:** CSV/TSV defining reportable genes and rules.

## The ACMG Rules File (`--acmg-table`)
The pipeline requires a gene-specific rules file. This file acts as the "control center," dictating which genes are analyzed and how specific variants are treated.

**Critical Columns:**
1.  **Gene Symbol:** HUGO symbol (e.g., *BRCA1*, *TTN*).
2.  **Reportable:** Only genes marked suitable are processed.
3.  **MOI (Inheritance):** `AD`, `AR`, `X-linked`. Used for PM3 (in-trans) and zygosity checks.
4.  **Disease/Phenotype:** Used to check against ClinVar traits.
5.  **Variants to report/Notes:** Scanned for keywords like:
    - *"Truncating variants only"*: Disables reporting of missense.
    - *"Biallelic"*: Triggers recessive (PM3) logic.
    - *"GOF"*: Disables automated PVS1 (LOF) criteria.

## Usage

### 1. Single Sample Analysis
Process a single proband VCF (and optional parents for trio analysis).

```bash
python acmg_sf_classifier.py \
  --proband sample.vcf.gz \
  --father father.vcf.gz \
  --mother mother.vcf.gz \
  --outdir results/sample_id \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --run-vep-dbnsfp \
  --vep /path/to/vep \
  --vep-cache /path/to/.vep \
  --fasta /path/to/hg38.fa \
  --db-paths-json config/db_paths.json
```

### 2. Batch Mode
Recursively searches a directory for `*_proband.vcf` files and aggregates results.

```bash
python acmg_sf_classifier.py \
  --batch-input-dir data/data_batch \
  --outdir batch_results \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --db-paths-json config/db_paths.json
```

## Command Line Arguments

- `--proband`: Proband VCF (bgz recommended).
- `--acmg-table` **[REQUIRED]**: ACMG rules CSV/TSV.
- `--outdir`: Output directory.
- `--ttn-meta`: TTN meta-exon CSV (optional, restricts PVS1 in TTN).
- `--gnomad-v2` / `--gnomad-v3` / `--gnomad-v4`: Paths to population DBs.
- `--clinvar-tabix`: Path to ClinVar VCF.
- `--dbnsfp`: Path to dbNSFP database.
- `--run-vep-dbnsfp`: Trigger VEP execution.
- `--aggressive`: Enable stricter filtering for VUS (use high-quality DB evidence to upgrade VUS -> P/LP).

## Design & Rules Summary

### PVS1 (Loss-of-function)
- Checks LOF mechanism validity per gene.
- NMD detection via exon coordinates.
- Canonical splice site (+/- 1-2bp).
- *TTN* meta-exon PSI checks (if `--ttn-meta` provided).

### PS1 / PM5
- Exact protein match (PS1) or codon match (PM5) to ClinVar Pathogenic entries (Requires >= 2 Stars or no conflicts).
- Validates Phenotype match between ACMG table and ClinVar trait.

### PM3 (In-trans)
- Uses trio data (if provided) to confirm `in-trans` status.
- In Batch/Singleton mode without parents, uses statistical heuristics to assign PM3_Supporting/Moderate.

### Pre-computation
The pipeline caches the ClinVar protein index to speed up subsequent runs.

## License & Clinical Note
**This software is a decision-support tool.** It is validated and integrated into the workflow of the **"Genome" Shared Resource Centre** at the Research Centre for Medical Genetics (Moscow, Russia).

All findings must be confirmed by a certified clinical geneticist and orthogonal methods (e.g., Sanger sequencing).

All rights reserved. **Patent Pending**
