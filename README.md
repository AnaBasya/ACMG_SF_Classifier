# ACMG_SF_Classifier

**Automated Prioritization Algorithm for Secondary Findings in ACMG SF v3.3 Genes**

## Overview

**ACMG_SF_Classifier** is a production-grade, fully automated pipeline for the prioritization and interpretation of secondary findings in genes from the American College of Medical Genetics and Genomics (ACMG) Secondary Findings version 3.3 (SF v3.3) list. The algorithm translates ACMG/AMP guidelines into a rigorous, quantitative framework based on the Bayesian model proposed by Tavtigian et al. (2018) and refined by the Association for Clinical Genomic Science (ACGS 2024).

The pipeline has been validated on over 7,800 clinical samples at the Research Centre For Medical Genetics (RCMG), Moscow, and is registered with the Russian Federal Service for Intellectual Property (certificate No. 2026613737). It achieves **98.78% sensitivity** and **99.73% specificity** on real-world clinical data, reducing manual expert review effort by >1,000-fold and processing a whole-genome sample in approximately **12 seconds**.

### Key Features

- **Weighted Bayesian scoring system** implementing ACMG/AMP criteria with empirically derived weights
- **Novel PP5ext/BP5ext scoring** that mathematically integrates database evidence (ClinVar, HGMD) with phenotype concordance
- **Transcript-aware NMD prediction** for precise PVS1 assignment using exon structure and the 50–55 nucleotide rule
- **Phasing-aware PM3 batching** for recessive gene analysis
- **Automated carrier filtering** for heterozygous variants in autosomal-recessive genes
- **Delins cluster detection** for complex allele identification
- **Multi-level triage logic** distributing variants into automatic reporting, manual review, or filtering
- **Batch processing support** for high-throughput clinical workflows



## Intellectual Property Notice

**ACMG_SF_Classifier** is registered intellectual property of the Research Centre For Medical Genetics (RCMG), Moscow, Russian Federation.

- **Certificate No.:** 2026613737
- **Application No.:** 2026611977 (January 30, 2026)
- **Date of Registration:** February 9, 2026

**Authors:** Anna Basova, Timofei Vizerov, Artem Borovikov, Kseniya Zabudskaya, Nikita Beskorovainiy, Ekaterina Malysheva, Oxana Ryzhkova, Sergey Kutsev

**Copyright © 2026 Research Centre For Medical Genetics. All rights reserved.**

Unauthorized commercial use, reproduction, or distribution of this software or its associated algorithms without prior written permission from the copyright holder is strictly prohibited. For licensing inquiries, please contact the Research Centre For Medical Genetics.

## Repository Structure

```
.
├── acmg_sf_classifier_POINTS_3_points.py   # RECOMMENDED: Production version (threshold = 3 points)
├── acmg_sf_classifier_POINTS.py            # Points-based version (older)
├── acmg_sf_classifier_NO_VEP.py            # Legacy: no VEP, no extended scoring
├── acmg_sf_classifier_VEP.py               # Legacy: VEP, no extended scoring
├── ACMG_SF_v3_3_full.csv                   # Primary gene rules table
├── nmd_data/                               # NMD/exon resources
├── scripts/                                # Helper utilities
```

## Version Comparison

The repository contains multiple versions reflecting the iterative development process.

**The recommended version is `acmg_sf_classifier_POINTS_3_points.py`**, which demonstrated optimal performance metrics on the validation cohort (sensitivity: 99.99%, specificity: 99.99%) and the main clinical cohort (sensitivity: 98.78%, specificity: 99.73%).

## Weighted ACMG/AMP Scoring System

The algorithm implements the Bayesian framework proposed by Tavtigian et al. (2018) and refined by ACGS (2024). Each ACMG/AMP criterion is assigned a numerical weight based on its strength level.

### Pathogenic Criteria

| Criterion | Strength | Points |
|--|-|--|
| PVS1 | Very Strong | +8 |
| PVS1 | Strong | +4 |
| PVS1 | Moderate | +2 |
| PVS1 | Supporting | +1 |
| PS1 | Strong | +4 |
| PS1 | Moderate | +2 |
| PS1 | Supporting | +1 |
| PS2 | Very Strong | +8* |
| PS2 | Moderate | +2* |
| PM1 | — | +2 |
| PM2 | Moderate | +2 |
| PM3 | Very Strong | +8 |
| PM3 | Strong | +4 |
| PM3 | Moderate | +2 |
| PM3 | Supporting | +1 |
| PM4 | Strong | +4 |
| PM4 | Moderate | +2 |
| PM4 | Supporting | +1 |
| PM5 | Moderate | +2 |
| PM5 | Supporting | +1 |
| PP1 | Strong | +4 |
| PP1 | Moderate | +2 |
| PP1 | Supporting | +1 |
| PP3 | Moderate | +2 |
| PP3 | Supporting | +1 |
| PP5ext | Very Strong | +8 |
| PP5ext | Strong | +4 |
| PP5ext | Moderate | +2 |
| PP5ext | Supporting | +1 |

\* PS2 criteria are **completely suppressed** in the automated branch to prevent false positives in monogenic opportunistic screening.

### Benign Criteria

| Criterion | Strength | Points |
|--|-|--|
| BA1 | Stand-alone | −8 |
| BS1 | Strong | −4 |
| BS2 | Strong | −4 |
| BP4 | Supporting | −1 |
| BP5ext | Supporting | −1 |
| BP5ext | Moderate | −2 |
| BP5ext | Strong | −4 |
| BP6 | Supporting | −1 |
| BP7 | Supporting | −1 |

### Classification Thresholds

| Classification | Points |
|-|--|
| Pathogenic | ≥10 |
| Likely Pathogenic | 6–9 |
| VUS | 0–5 |
| Likely Benign | −1 to −5 |
| Benign | ≤−6 |

### Safety Guards

- **Minimum pathogenic criteria:** At least 2 independent pathogenic criteria (PVS/PS/PM/PP) required for P/LP classification
- **Pathogenic downgrade:** If exactly 2 pathogenic criteria present → downgrade to Likely Pathogenic
- **VUS downgrade:** If <1 strong pathogenic criterion → downgrade to VUS
- **Carrier suppression:** Single heterozygous P/LP variants in AR genes are automatically filtered



## Extended Database Scoring System (PP5ext/BP5ext)

The PP5ext/BP5ext system provides nuanced, mathematically weighted evidence from ClinVar and HGMD. Coefficients were empirically derived by expert interpreters at RCMG and validated on clinical data.

### ClinVar Scoring

| Classification | Stars | Raw Block | Description |
|-|-|--|-|
| Pathogenic | 4★ | 4.5 | Practice guideline (e.g., ClinGen) |
| Pathogenic | 3★ | 3.0 | Expert panel (e.g., ClinGen VCEP) |
| Pathogenic | 2★ | 3.0 | Multiple submitters, no conflict |
| Likely Pathogenic | 3★ | 2.5 | Expert panel |
| Likely Pathogenic | 2★ | 2.0 | Multiple submitters, no conflict |
| Pathogenic | 1★ | 2.25 | Single submitter with criteria |
| Likely Pathogenic | 1★ | 1.5 | Single submitter with criteria |
| Pathogenic | 0★ | 1.0 | No criteria provided |

### HGMD Scoring

| Class | Raw Block | Description |
|-|--|-|
| DM (dmsupported ≥0.9) | 6.0 | Disease-causing mutation, strong support |
| DM (dmsupported 0.7–0.89) | 4.0 | Disease-causing mutation, good support |
| DM (dmsupported 0.4–0.69) | 3.0 | Disease-causing mutation, moderate support |
| DM? (dmsupported ≥0.7) | 1.5 | Likely disease-causing |
| DM? (dmsupported 0.4–0.69) | 1.0 | Likely disease-causing, limited support |
| DFP | −1.0 | Disease-associated polymorphism |
| DP | −2.0 | Disease-associated polymorphism (functional) |
| FP | −2.0 | Functional polymorphism |

### Strength Mapping

| Raw Total | PP5ext Strength | Points | BP5ext Strength | Points |
|--|--|--|--|--|
| ≥8.0 | Very Strong | +8 | — | — |
| 4.0–7.9 | Strong | +4 | — | — |
| 2.0–3.9 | Moderate | +2 | — | — |
| 1.0–1.9 | Supporting | +1 | — | — |
| −1.0 to −1.9 | — | — | Supporting | −1 |
| −2.0 to −3.9 | — | — | Moderate | −2 |
| ≤−4.0 | — | — | Strong | −4 |

**Phenotype Match Modifier:** If no phenotype match exists between the database entry and the ACMG disease, the raw total is capped at ±1.0 to prevent over-interpretation of phenotype-mismatched evidence.



## Manual Review Threshold Selection

The optimal threshold of **3 points** was empirically determined during validation. The table below illustrates the sensitivity vs. manual review burden trade-off:

| Threshold | Sensitivity | Specificity | Variants to Manual Review | F1-Score |
|--|-|-|-|-|
| 6 points | 91.67% | 99.87% | 1,847 | 68.75% |
| 5 points | 93.06% | 99.84% | 2,480 | 70.53% |
| **3 points** | **98.78%** | **99.73%** | **3,889** | **59.11%** |
| 2 points | 99.32% | 99.63% | 5,673 | 54.78% |
| 1 point | 99.53% | 99.51% | 8,475 | 47.23% |

**At the 3-point threshold:**
- **Sensitivity:** 98.78% (only 3 false negatives out of 430 secondary findings)
- **Specificity:** 99.73%
- **Manual review burden:** 3,889 variants across 7,827 samples (0.10% of all ACMG-filtered variants)
- **Processing time:** ~12 seconds per whole-genome sample



## Performance Metrics

### Validation Cohort (43 samples, 47 complex variants)

| Metric | Value |
|--|-|
| Accuracy | 99.99% |
| Sensitivity (Recall) | 99.99% |
| Specificity | 99.99% |
| Precision | 69.50% |
| F1-Score | 71.89% |

### Main Clinical Cohort (7,827 samples, 430 secondary findings)

| Metric | Value |
|--|-|
| Accuracy | 99.93% |
| Sensitivity (Recall) | 98.78% |
| Specificity | 99.73% |
| Precision | 53.31% |
| F1-Score | 59.11% |

### Secondary Finding Frequency

- **Overall frequency:** 5.27% (430 findings in 413 patients)
- **Most frequent genes:** BRCA1 (11.9%), LDLR (9.5%), BRCA2 (8.6%), TTN (7.4%), APOB (4.2%)



## Quick Start

### Prerequisites

```bash
# Python 3.11 or higher required
python3 --version

# Install dependencies
pip install -r scripts/requirements.txt
```

### Single Sample 

```bash
python3 acmg_sf_classifier_POINTS_3_points.py \
  --proband sample.vcf.gz \
  --outdir results/sample \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --gnomad-v4 /ngs/gnomad/gnomad.genomes.v4.0.vcf.gz \
  --gnomad-v3 /ngs/gnomad/v3.1.2/gnomad.genomes.v3.1.2.hg38.tsv.bgz \
  --gnomad-v2 /ngs/gnomad/gnomad.genomes.r2.1.sites.tsv.bgz \
  --hgmd /ngs/hgmd/hgmd_2022_vars_dmsupport.csv \
  --exons-file nmd_data/acmg_simple_exons.csv \
  --nmd-table nmd_data/acmg_nmd_table.tsv
```

### Batch Mode (Cohort Processing)

```bash
python3 acmg_sf_classifier_POINTS_3_points.py \
  --batch-input-dir ./data/all_converted_vcf/  \
  --outdir ./3pts_15_9_results/ \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --gnomad-v4 /ngs/gnomad/gnomad.genomes.v4.0.vcf.gz \
  --gnomad-v3 /ngs/gnomad/v3.1.2/gnomad.genomes.v3.1.2.hg38.tsv.bgz \
  --gnomad-v2 /ngs/gnomad/gnomad.genomes.r2.1.sites.tsv.bgz \
  --hgmd /ngs/hgmd/hgmd_2022_vars_dmsupport.csv \
  --exons-file nmd_data/acmg_simple_exons.csv /
  --nmd-table nmd_data/acmg_nmd_table.tsv
```



## Command-Line Arguments

| Argument | Description | Required |
|-|-|-|
| `--proband` | Proband VCF file (bgzipped) | Yes* |
| `--father` | Father VCF file (optional) | No |
| `--mother` | Mother VCF file (optional) | No |
| `--batch-input-dir` | Directory containing samples for batch processing | Yes* |
| `--outdir` | Output directory | No (default: acmg_sf_results) |
| `--acmg-table` | ACMG SF v3.3 gene rules table | Yes |
| `--db-paths-json` | JSON file with database paths | No |
| `--exons-file` | Transcript/exons CSV for NMD prediction | No |
| `--nmd-table` | Precomputed NMD table (TSV) | No |
| `--transcript-map` | Ensembl ↔ RefSeq transcript mapping | No |
| `--gene-bed` | BED file for pre-filtering (expanded by 20kb) | No |
| `--gnomad-v2` | gnomAD v2.1.1 site file | No |
| `--gnomad-v3` | gnomAD v3.1.2 site file | No |
| `--gnomad-v4` | gnomAD v4.1.0 site file | No |
| `--internal-db` | Internal database CSV | No |
| `--hgmd` | HGMD Professional CSV | No |
| `--acmg-recs` | Additional gene recommendations | No |
| `--proband-sample` | Sample name (for multi-sample VCFs) | No |
| `--no-strict-qc` | Disable strict QC filtering | No |

\* Either `--proband` or `--batch-input-dir` must be provided.



## Output Files

| File | Description |
|-|-|
| `all_candidates.csv` | Complete annotation of all variants in ACMG SF v3.3 genes |
| `auto_conclusions.csv` | Variants directed to automatic reporting (Pathogenic/Likely Pathogenic) |
| `manual_review_list.csv` | Variants requiring expert manual review |

### Output Columns

Key columns include:
- Variant identifiers (chrom, pos, ref, alt, gene)
- Functional annotation (HGVSc, HGVSp, consequence, impact)
- Population frequencies (gnomAD AF, version)
- Computational predictions (REVEL, SpliceAI, AlphaMissense, Dscore)
- Database evidence (ClinVar, HGMD, internal DB)
- ACMG criteria assigned and per-criterion explanations
- Total points and classification
- Triage flags (`in_auto_report`, `in_manual_review`, `filtered_out`)



## Citation

If you use ACMG_SF_Classifier in your research, please cite:

**Software:**

```bibtex
@software{Basova_ACMG_SF_Classifier_2026,
  author = {Anna Basova, Timofei Vizerov, Artem Borovikov, Kseniya Zabudskaya, Nikita Beskorovainiy, Ekaterina Malysheva, Oxana Ryzhkova, Sergey Kutsev},
  title = {ACMG_SF_Classifier: Automated Prioritization Algorithm for Secondary Findings in ACMG SF v3.3 Genes},
  year = {2026},
  url = {https://github.com/AnaBasya/ACMG_SF_Classifier},
  note = {Software registered with the Russian Federal Service for Intellectual Property (Certificate No. 2026613737)}
}
```



## License

```
MIT License

Copyright (c) 2026 Research Centre For Medical Genetics (RCMG)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

**Note:** While the software is open-source under the MIT License, the associated algorithm and methodology are registered intellectual property of the Research Centre For Medical Genetics (certificate No. 2026613737). Commercial use of the algorithm or its derivatives without prior written permission is prohibited.



## Contact

**Anna Basova**  
Email: anna.m.basova@gmail.com  
GitHub: [AnaBasya](https://github.com/AnaBasya)

**Research Centre For Medical Genetics (RCMG)**  
Moscow, Russian Federation  
Website: [https://med-gen.ru](https://med-gen.ru)



## Acknowledgments

This work was supported by the Research Centre For Medical Genetics. The authors thank all members of the bioinformatics and genetics teams for their contributions to data generation and expert interpretation. 


*Developed for clinical NGS workflows at RCMG ("Genome" Centre).*  
