# ACMG Secondary Findings Classifier

## Overview
This repository implements a production-grade bioinformatics pipeline for automating **ACMG Secondary Findings (SF v3.2/v3.3)** classification. 

It is designed to intake VCF files (Single or Trio), annotate them, cross-reference against clinical and population databases, and apply rigorous ACMG/AMP criteria (PVS1, PS1, PM2, PP3, etc.) to identify Reportable Variants.

## Repository Structure

```
.
├── acmg_sf_classifier.py        # Main pipeline script
├── ACMG_SF_v3.3_full.csv        # Primary Rules Engine (Genes, MOI, Phenotypes)
├── config/                      # JSON configurations (DB paths overrides)
├── data/                        # Input VCFs and helper data
│   ├── exons.csv                # Transcript definitions for NMD analysis
│   ├── header.vcf               # Template headers
│   └── ...
├── databases/                   # Reference Databases
│   ├── clinvar/                 # ClinVar VCF (for PS1/PM5)
│   ├── gnomad/                  # Population frequencies (v2/v3/v4)
│   ├── vep/                     # VEP Cache
│   ├── dbnsfp/                  # dbNSFP database (REVEL, SpliceAI, CADD)
│   └── ...
├── scripts/                     # Helper utilities
│   ├── requirements.txt         # Python dependencies
│   ├── make_ttn_meta_csv.py     # Utility for TTN transcripts
│   └── ...
├── results/                     # Output for Single-Proband runs (Singletons/Trios)
└── batch_results/               # Aggregated output for Batch runs
    ├── individual_results/      # Detailed logs per sample from batch run
    ├── FINAL_auto_conclusions.csv
    └── FINAL_manual_review.csv
```


## Detailed Workflow Logic

### 1. Robust Gene Normalization
To resolve discrepancies between HGNC symbols used in VCFs, ClinVar, and the ACMG table, the pipeline uses a strict normalization strategy:
- Strips whitespace/quotes.
- Uppercases strings.
- Replaces non-alphanumeric characters with underscores.
- **Impact:** Ensures `BRCA1`, `brca1`, and ` BRCA-1 ` are treated identically for database lookups and rule application.

### 2. Multi-Source gnomAD Aggregation
The pipeline supports simultaneous input of **gnomAD v2, v3, and v4**.
- **Formats:** Accepts both `.vcf.bgz` (slow, exact) and TSV files (`.bgz` + `tabix`, fast).
- **Logic:** For every variant, it extracts frequencies from *all* provided versions.
- **Aggregation:** Calculates the **Global Max AF** across all datasets.
- **Purpose:** Used conservatively for:
  - **BA1/BS1 (Benign):** Requires AF > 5% / 1% in *any* major population to exclude a variant.
  - **PM2 (Rare):** Requires AF < threshold (e.g., 0.0001) in *all* databases to apply the criteria.

### 3. Execution Modes

#### A. Single Proband Mode (`--proband`)
Analyzes one individual (optionally with parents).
- **Trio Analysis:** If `--father` and `--mother` are provided, the script uses parental genotypes to check for **De Novo (PS2)** and **In-Trans (PM3)** status.
- **Outputs:** Creates a dedicated directory defined by `--outdir`. Generates:
  - `auto_conclusions.csv` (Pathogenic/Likely Pathogenic).
  - `manual_review_list.csv` (VUS needing attention, AR single hits).
  - `all_candidates.csv` (Raw list of all variants in ACMG genes passing QC).

#### B. Batch Cohort Mode (`--batch-input-dir`)
Recursively searches a directory tree to process whole cohorts.
- **Discovery:** Scans the folder structure for files ending in `*_proband.vcf`.
- **Family Inference:** Automatically detects `*_father.vcf` and `*_mother.vcf` in the same subdirectory as the proband.
- **Aggregation:** 
  1. Runs the pipeline for each sample individually (results stored in `individual_results/`).
  2. Collapses all findings into master files: `FINAL_auto_conclusions.csv` and `FINAL_manual_review.csv` in the output root.

### 4. Database & Criteria Engine

| Criteria | Logic / Database Used |
| :--- | :--- |
| **PVS1** | **Loss-of-Function.** Uses `data/exons.csv` to check for NMD triggers (Last Exon rule). Supports gene-specific overrides (e.g., TTN meta-exons via `--ttn-meta`) and excludes GOF genes defined in the ACMG table. |
| **PS1 / PM5** | **ClinVar.** Built-in indexing of ClinVar VCF (Protein/Splice index). Requires >= 2 Stars (or 1 star w/o conflicts) for P/LP entries to trigger these points. Matches Phenotype strings between ACMG table and ClinVar. |
| **PM2** | **Rareness.** Strict absence/rarity in aggregate gnomAD data. |
| **PP3** | **In-Silico.** VEP + dbNSFP. Triggers based on thresholds for REVEL, AlphaMissense, and SpliceAI. |
| **BA1 / BS1** | **Common Variant.** Automated exclusion if Max AF > 0.05 (BA1) or > 0.01 (BS1). |

## Usage Examples

### 1. Single Sample (Manual paths)
```bash
python acmg_sf_classifier.py \\
  --proband data/31233/proband.vcf.gz \\
  --father data/31233/father.vcf.gz \\
  --mother data/31233/mother.vcf.gz \\
  --outdir results/31233 \\
  --acmg-table ACMG_SF_v3.3_full.csv \\
  --run-vep-dbnsfp \\
  --db-paths-json config/db_paths.json
2. Batch Mode (Processing a Cohort)

This will recursively process all samples found in data/data_batch/ and save aggregated results to batch_results/.

code
Bash
download
content_copy
expand_less
python acmg_sf_classifier.py \\
  --batch-input-dir data/data_batch \\
  --outdir batch_results \\
  --acmg-table ACMG_SF_v3.3_full.csv \\
  --db-paths-json config/db_paths.json \\
  --exons-file data/exons.csv
Setup & Requirements

Environment:

code
Bash
download
content_copy
expand_less
pip install -r scripts/requirements.txt
# Requires: pandas, pysam, cyvcf2, intervaltree

System Tools:

vep (if using --run-vep-dbnsfp)

bcftools + tabix (required for VCF manipulation)

Outputs

auto_conclusions.csv: High-confidence Pathogenic/Likely Pathogenic variants. Ready for report drafting.

manual_review_list.csv: Complex cases (Recessive single hits, VUS with strong in-silico scores, ClinVar conflicts).

run_info.json: Metadata about the run (counts, paths, timestamps).

Developed for clinical NGS workflows at RCMG ("Genome" Centre).