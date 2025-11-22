# ACMG Secondary Findings Classifier

## Overview
This repository contains a production-oriented implementation of an automated-first ACMG Secondary Findings (SF) classifier.

### Goals
- **Maximise safe automation:** automatically classify and report variants where evidence is robust and conforms to decision rules.
- **Minimise manual workload for clinical geneticists:** use high-quality DB evidence and conservative computational thresholds to push variants to automated conclusions where appropriate.
- **Preserve safety:** gene-specific exceptions (e.g., TTN meta-exon, genes where LOF is not the disease mechanism) still trigger manual review.

## Contents
- `acmg_sf_pipeline.py`
  - Main pipeline script. Performs VEP annotation (if necessary), applies gene-specific rules and ACMG criteria, consults databases, and writes outputs.
- `make_ttn_meta_csv.py`
  - Utility to create a TTN meta-exon CSV used by the classifier. Normalises PSI values and produces an ordered CSV with `chr,start,end,exon_id,psi_meta,lof_allowed`.
- `README.md`
  - This document.

## Requirements

### Python Environment
- Python 3.8+
- Required pip packages:
  - `pandas`
  - `pysam`
  - `cyvcf2`
  - `intervaltree` (only required if using TTN meta-exon support)

### System Tools
- **VEP (local)** + installed plugins (SpliceAI, REVEL, CADD, etc.) if you want the script to run VEP.
- **tabix** (optional, for indexing VEP output).

### Databases
Provide paths via `--db-paths-json` or use defaults:
- ClinGen gene validity CSV
- ClinGen variant assertions CSV
- ClinVar VCF (recommended; required for ClinVar-based evidence)
- HGMD Pro VCF (optional)
- gnomAD VCF (optional)
- INTERNAL_DB CSV (optional)

## The ACMG Rules File (`ACMG_SF_v3.3.csv` / `.tsv`)
The pipeline requires a gene-specific rules file (passed via `--acmg-tsv`). This file acts as the "control center" for the classifier, dictating which genes are analyzed and how specific variants are treated.

**Format:** Tab-separated (TSV) or Comma-separated (CSV) file containing the following critical columns:
1.  **Gene Symbol:** The HUGO symbol (e.g., *BRCA1*, *TTN*).
2.  **Reportable as SF:** Values `Yes`/`No`. Only genes marked `Yes` are processed.
3.  **MOI (or Inheritance):** Mode of inheritance (e.g., `AD`, `AR`, `X-linked`). Used for PM3 and zygosity checks.
4.  **Variants to report:** Free text describing scope (e.g., "All P/LP", "Truncating only").
5.  **Reporting Guidance / Mechanism:** text fields scanned for keywords to trigger special flags:
    - **"Truncating variants only":** Disables reporting of missense variants for this gene.
    - **"Biallelic" / "2 variants":** triggers specific logic for Recessive conditions.
    - **"GOF" / "Gain of function":** Disables automated PVS1 (LOF) criteria, as LOF is not the mechanism.

## VEP Configuration Guide
Since VEP is often installed in a separate Conda environment or specific system location, correct configuration is crucial.

### 1. Setting the VEP Executable Path (`--vep`)
You do not need to add VEP to your global `PATH`. Instead, find the full path to the executable and pass it to the script.
*Example:* `--vep /home/user/miniconda3/envs/vep/bin/vep`

### 2. Locating the VEP Cache (`--vep-cache`)
The cache directory is required for annotation. It usually resides in a hidden folder in your home directory or within the Conda environment.
*Example:* `--vep-cache /home/user/.vep`

## Command Line Arguments
The script accepts the following arguments. Flags marked **[REQUIRED]** must be provided for the script to run.

### Input Files
*   `--proband` **[REQUIRED]**
    *   **Description:** Path to the Proband's VCF file. BGZIP compressed (`.vcf.gz`) is recommended. If the VCF is not annotated, the script will attempt to run VEP on it.
*   `--father` *(Optional)*
    *   **Description:** Path to the Father's VCF file. Used for phasing (PM3) and de novo (PS2/PM6) analysis.
*   `--mother` *(Optional)*
    *   **Description:** Path to the Mother's VCF file. Used for phasing (PM3) and de novo (PS2/PM6) analysis.

### Configuration & Rules
*   `--acmg-tsv` **[REQUIRED]**
    *   **Description:** Path to the ACMG Secondary Findings rules file (e.g., `ACMG_SF_v3.3.csv` or `.tsv`). This file defines reportable genes and specific handling rules (see section above).
*   `--ttn-meta` *(Optional)*
    *   **Description:** Path to the TTN meta-exon CSV.
    *   **Columns:** `chr`, `start`, `end`, `psi_meta`, `lof_allowed`.
    *   **Function:** If provided, restricts PVS1 application in *TTN* to exons with high PSI (percent spliced-in) values.
*   `--db-paths-json` *(Optional)*
    *   **Description:** Path to a JSON file mapping database keys to file paths. Use this to override default paths for ClinVar, HGMD, gnomAD, etc.

### VEP Settings
*   `--vep` *(Default: "vep")*
    *   **Description:** The command or full path to the VEP executable. If VEP is in a specific conda environment, provide the full path (e.g., `/envs/vep/bin/vep`).
*   `--vep-cache` *(Optional)*
    *   **Description:** Path to the VEP cache directory (containing species folders). Required if running VEP locally.
*   `--fasta` *(Optional)*
    *   **Description:** Path to the reference genome FASTA file. Required by VEP for HGVS calculation.
*   `--vep-extra` *(Default: [])*
    *   **Description:** Extra arguments to pass to VEP, typically for plugins.
    *   **Usage:** Tokenized list. Example: `--vep-extra "--plugin" "SpliceAI,snv,0.2" "--plugin" "REVEL"`

### Outputs & Execution
*   `--outdir` *(Default: "acmg_sf_results")*
    *   **Description:** Directory where output CSVs and logs will be saved. Created if it does not exist.
*   `--aggressive` *(Flag)*
    *   **Description:** Enables "Aggressive" DB-driven overrides.
    *   **Function:** If set, variants calculated as VUS by the algorithm will be automatically upgraded to **Pathogenic** or **Likely Pathogenic** if strong evidence exists in high-quality databases (ClinVar 2+ stars, HGMD DM with publications, etc.). This reduces manual review volume but relies heavily on database accuracy.

## Design & Rules Summary

### PVS1: Loss-of-function decision tree
Implemented using:
- Gene LOF mechanism confirmation via ClinGen validity or TSV guidance.
- NMD tags (from VEP CSQ).
- Canonical splice site handling with SpliceAI >= 0.20.
- Start-loss heuristics and last-exon handling (PM4).
- TTN meta-exon handling using provided TTN meta-exon CSV.

### PS1/PM5
- Uses ClinGen variant assertions, Internal DB, ClinVar protein-level evidence (indexed and cached), and HGMD DM (>=2 PMIDs).
- ClinVar protein-level matching requires review status **>=2★** for inclusion.

### PP3/BP4
- **Missense:** REVEL thresholds (support 0.644, moderate 0.932), MutPred2 optional.
- **Splicing:** SpliceAI single threshold 0.20 (>= 0.2 supports splice impact).
- CADD >=25 recorded as informative.

### PM3
- Batch processing across variants in the same sample/gene to determine compound-heterozygosity evidence.
- Uses parental genotypes for **in-trans** determination where available; otherwise assigns supporting/moderate PM3.

### PS2/PM6
- De novo calling uses parental genotypes and allele-balance (AB) + depth (DP) heuristics to reduce false positives.

## ClinVar protein-level index & cache
- Building an index from ClinVar for protein-level matches (p.HGVS) is costly.
- The pipeline builds this index once and stores it as a gzipped JSON cache (path configured in `DB_PATHS` or db_paths JSON under `CLINVAR_INDEX_CACHE`).
- On subsequent runs, the cache is loaded if the ClinVar VCF mtime equals cached mtime, saving significant startup time.
- If ClinVar VCF is updated, cache is automatically detected as stale and rebuilt.

## Usage Example

```bash
python acmg_sf__pipeline.py \
  --proband sample.vcf.gz \
  --acmg-tsv "ACMG_SF_v3.3.csv" \
  --outdir results \
  --vep /home/user/miniconda3/envs/vep/bin/vep \
  --vep-cache /home/user/.vep \
  --fasta /path/to/human.fasta \
  --vep-extra "--plugin" "SpliceAI,snv,0.2" "--plugin" "REVEL" \
  --ttn-meta ttn_meta.csv \
  --db-paths-json db_paths.json \
  --father father.vcf.gz --mother mother.vcf.gz \
  --aggressive
```

## License & Clinical Note
This code is a clinical tool template, it is validated and integrated into workflow of the "Genome" Shared Resource Centre at the Research Centre for Medical Genetics (Moscow, Russia).

All rights reserved + **Patent Pending**
