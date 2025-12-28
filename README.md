# ACMG Secondary Findings Classifier — README (updated)

This README describes the current production-ready pipeline implemented in `acmg_sf_classifier.py` for automated classification of ACMG Secondary Findings (SF v3.x).  
The pipeline annotates VCFs, aggregates multi-source evidence (ClinVar, HGMD, internal DB, gnomAD), predicts NMD/last-exon effects, applies gene-specific rules from the ACMG SF table, and performs conservative automatic vs manual triage.


# ACMG Secondary Findings Classifier

## Overview

This repository implements a production-grade bioinformatics pipeline for automating **ACMG Secondary Findings (SF v3.3)** classification. 

The pipeline:
- Accepts single-sample or trio VCFs (supports batch/cohort mode).
- Optionally annotates with VEP + dbNSFP and SpliceAI plugin.
- Aggregates gnomAD v2/v3/v4 allele frequencies and uses the maximum observed AF for population-based criteria.
- Integrates ClinVar (including a protein-level index for PS1/PM5), HGMD, and an internal database (CSV).
- Uses transcript/exon resources (exons CSV, transcript map, NMD table) to make transcript-aware PVS1/NMD decisions.
- Implements PM3 batching to evaluate AR pairing and phasing evidence.
- Detects delins/adjacent indel clusters and applies phasing-based rules.
- Produces 3 main outputs per run: `all_candidates.csv`, `auto_conclusions.csv`, `manual_review_list.csv` plus `run_info.json`.

## Repository Structure

```
.
├── acmg_sf_classifier.py        # Main pipeline script
├── ACMG_SF_v3.3_full.csv        # Primary Rules Engine (Genes, MOI, Phenotypes, Rules)
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
└── batch_results_final/         # Aggregated output for Batch runs
    ├── individual_results/      # Detailed logs per sample from batch run
    ├── FINAL_all_candidates.csv
    ├── FINAL_auto_conclusions.csv
    └── FINAL_manual_review.csv
```

## Quick start

Install dependencies (recommended in a dedicated environment):

```bash
pip install -r scripts/requirements.txt
# required: pandas, pysam, cyvcf2 (optional), intervaltree (optional)
```

Single/proband example:

```bash
python3 acmg_sf_classifier.py \
  --proband data/31233/proband.vcf.gz \
  --father data/31233/father.vcf.gz \
  --mother data/31233/mother.vcf.gz \
  --outdir results/31233 \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --run-vep-dbnsfp \
  --db-paths-json config/db_paths.json
```

Batch mode (cohort):

```bash
python3 acmg_sf_classifier.py \
  --batch-input-dir data/data_batch \
  --outdir batch_results \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --db-paths-json config/db_paths.json \
  --exons-file data/exons.csv
```

Your team’s batch example (as provided):

```bash
python3 acmg_sf_classifier.py \
  --batch-input-dir ./data/data_batch \
  --outdir ./batch_results_final \
  --acmg-table ACMG_SF_v3.3_full.csv \
  --vep vep \
  --vep-cache /home/anna/anna/ACMG_SF_Classifier/databases/vep/.vep \
  --fasta ~/ngs-data/ngs/ref/hg38.fa \
  --gene-bed ./data/mane_out/mane_intervals.sorted.bed \
  --run-vep-dbnsfp \
  --dbnsfp /home/anna/ngs-data/ngs/dbNSFP/ver4.9/dbNSFP_hg38.gz \
  --dbnsfp-fields REVEL_score,CADD_phred,AlphaMissense_score \
  --gnomad-v4 ~/ngs-data/ngs/gnomad/v4/gnomad.genomes.v4.0.sites.vcf.gz \
  --gnomad-v3 ~/ngs-data/ngs/gnomad/v3.1.2/gnomad.genomes.v3.1.2.hg38.tsv.bgz \
  --gnomad-v2 ~/ngs-data/ngs/gnomad/gnomad.genomes.r2.1.sites.tsv.bgz \
  --internal-db ./databases/internal/export_all_annotations.csv \
  --acmg-recs ./ACMG_SF_Annotations.tsv \
  --hgmd /home/nik/share/ccu-ngs/ngs/bases/hgmd/dmsupport/hgmd_2025_vars_dmsupport.csv \
  --exons-file nmd_data/acmg_simple_exons.csv \
  --nmd-table nmd_data/acmg_nmd_table.tsv \
  --transcript-map refseq2enst_mane.tsv
```

---

## Inputs / Repository layout (key files)

- `acmg_sf_classifier.py` — main pipeline script (entry point).
- `ACMG_SF_v3.3_full.csv` — ACMG SF gene table (gene, diseases, MOI, rules_text, positional guidance).
- `databases/clinvar/` — ClinVar VCF (used for PS1/PM5 and ClinVar evidence).
- `databases/vep/` — VEP cache and plugin data (SpliceAI).
- `databases/dbnsfp/` — dbNSFP archive (for REVEL, CADD, AlphaMissense, SpliceAI fields).
- `databases/internal/` — optional internal DB CSV (vid, anntype, date etc).
- `exons.csv` / `nmd_table.tsv` — exon/transcript files used by NMD predictor.
- `results/`, `batch_results/` — outputs.

---

## Default thresholds and configuration (used by code)

Key thresholds used by the pipeline (these are configurable in the code):

- Computational evidence:
  - REVEL_SUPPORTING = 0.644
  - REVEL_MODERATE = 0.932
  - ALPHA_MISSENSE_SUPPORTING = 0.564
  - CADD_SUPPORTING = 30

- Splice predictors:
  - SPLICEAI_PVS1 = 0.70
  - SPLICEAI_MODERATE = 0.20
  - SPLICEAI_PS1_THRESHOLD = 0.10

- Population frequency:
  - BA1_AF = 0.05 (stand‑alone benign)
  - BS1_AF = 0.01 (strong benign)

- PM2 (MOI-aware):
  - PM2_AD_AF = 0.0001
  - PM2_AR_AF = 0.005
  - PM2_XLD_AF = 0.0001
  - PM2_XLR_AF = 0.003

- BA1 exceptions (example):
  - HFE: p.Cys282Tyr
  - BTD: p.Asp444His

- Strict QC filters:
  - STRICT_MIN_DP = 10 (heterozygous)
  - STRICT_MIN_AF = 0.20 (sample-level alternate allele fraction)
  - Homozygotes/hemizygotes: lower DP allowed (e.g. 5)

---

## Points map and classification thresholds

This pipeline uses a point-based (Tavtigian et al., (2018)) aggregation of criteria. The `POINTS_MAP` used in code:

- PVS1_VeryStrong = +8
- PVS1_Strong = +4
- PVS1_Moderate = +2
- PVS1_Supporting = +1

- PS1_Strong = +4
- PS1_Moderate = +2
- PS1_Supporting = +1

- PS2_VeryStrong = +8
- PS2_Moderate = +2

- PM1 = +2
- PM2_Moderate = +2

- PM3_VeryStrong = +8
- PM3_Strong = +4
- PM3_Moderate = +2
- PM3_Supporting = +1

- PM4_Strong = +4
- PM4_Moderate = +2
- PM4_Supporting = +1

- PM5_Moderate = +2
- PM5_Supporting = +1

- PP1_Strong = +4
- PP1_Moderate = +2
- PP1_Supporting = +1

- PP3_Moderate = +2
- PP3_Supporting = +1

- PP5_Supporting = +1

Benign/negative points:
- BA1 = −8
- BS1 = −4
- BS2 = −4
- BP4 = −1
- BP5 = −1
- BP6 = −1
- BP7 = −1

Classification mapping from summed points:
- Pathogenic: total_points >= 10
- Likely pathogenic: total_points in [6 .. 9]
- VUS: total_points in [0 .. 5]
- Likely benign: total_points in [−5 .. −1]
- Benign: total_points <= −6

Additional safeguard:
- MIN_PATHOGENIC_CRITERIA = 2 (minimum number of distinct pathogenic criteria PVS/PS/PM/PP required). If the point sum reaches Pathogenic/Likely pathogenic thresholds but there are fewer than `MIN_PATHOGENIC_CRITERIA` distinct pathogenic criteria, the automated class is downgraded to VUS.

---

## Detailed logic for assigning each ACMG criterion

Below is a precise description of how each major criterion is evaluated in the pipeline (mirrors the implementation).

### BA1 (Stand-alone Benign)
- If aggregated gnomAD max AF ≥ BA1_AF (0.05) and variant NOT in gene-specific BA1 exception list → assign BA1 (−8).
- BA1 short-circuits classification: variant is returned as Benign.

### BS1 (Strong Benign by frequency)
- If aggregated gnomAD max AF ≥ BS1_AF (0.01) → assign BS1 (−4).

### PVS1 (Loss of Function — decision tree)
Applied to obvious LOF consequences: stop_gained, frameshift, canonical splice (+/−1–2), start_lost, etc.

Decision flow:
1. Respect gene-specific override: if gene flagged `lof_not_reportable` in ACMG table → do not apply PVS1.
2. Internal/external NMD annotations:
   - If internal annotation (e.g., VEP NMD or `nmd` field) indicates escape (last exon) → treat as escape.
   - Else run transcript-aware predictor `variant_triggers_nmd()`:
     - Uses HGVSc cDNA offsets when present (c.123+1, c.456-2).
     - Maps cDNA → exon via `nmd_map` (RefSeq exon CDS lengths) and `tx_map` (ENST→NM).
     - Uses `exons_df` (if provided) to find last exon coordinates and penultimate exon terminal window (default 50 nt).
     - Predicts NMD = True only if variant lies upstream of penultimate terminal window for all relevant transcripts (conservative logic).
3. Assignments:
   - If predicted to trigger NMD → PVS1_VeryStrong (+8).
   - If predicted to escape NMD:
     - canonical splice in last exon → PM4_Moderate (+2)
     - stop/frameshift in last exon → PM4_Moderate (or PM4_Strong depending on truncation %; the pipeline prefers PM4_Moderate by default)
   - Start-loss: default PVS1_Moderate (+2). The pipeline searches for in-frame downstream ATGs; presence of plausible downstream start may affect level.
4. Internal NMD overrides are honored (if internal predictor marks last exon escape, PVS1 not applied).

### PM4 (Protein length change)
- In-frame indels and stop-loss:
  - small in-frame indels (<10 aa) outside critical domains → PM4_Supporting (+1)
  - moderate/important in-frame indels → PM4_Moderate (+2)
  - stop-loss predicted to significantly alter C‑terminus → PM4_Strong (+4)

### PS1 / PM5 (Protein-level evidence)
- PS1 (exact same AA change): uses ClinVar protein-level index (HGVS protein). Only high-quality ClinVar entries (stars ≥ 2) are used for strong PS1.
  - exact AA match Pathogenic (not likely) → PS1_Strong (+4).
  - likely pathogenic match → PS1_Supporting/PS1_Moderate (conservatively supporting).
- PM5 (different AA at same residue reported P/LP):
  - Pathogenic at same AA position (different AA) → PM5_Moderate (+2)
  - LP at same position → PM5_Supporting (+1)

### PM2 (Absence / rarity)
- MOI-aware frequency thresholds:
  - AD (dominant): PM2_AD_AF (0.0001)
  - AR (recessive): PM2_AR_AF (0.005)
  - XLR/XLD: use PM2_XLR_AF / PM2_XLD_AF
- If gnomAD max AF < MOI-specific threshold → PM2_Moderate (+2)

### PP3 / BP4 (Computational predictions)
- PP3 assignment:
  - SpliceAI ≥ 0.20 → PP3_Moderate (+2)
  - REVEL ≥ 0.932 → PP3_Moderate (+2)
  - CADD ≥ 30 → treated as strong supporting computational evidence (commonly used as PP3)
  - REVEL ≥ 0.644 or AlphaMissense ≥ 0.564 → PP3_Supporting (+1)
- IMPORTANT: PP3 is NOT applied if PVS1 was applied or if variant is a clear LOF/canonical splice consequence (to avoid stacking computational evidence with PVS1).

### PP5 (Reputable source)
- If high-quality database annotation exists (ClinVar ≥ 2★ P/LP, or HGMD DM with ≥ 2 pubs/dmsupported) and there is no HQ conflict → PP5_Supporting (+1).

### PM3 (Recessive / trans observations) — batching and phasing-aware
- PM3 scoring is batch-aggregated per sample/gene by `apply_pm3_batch()`:
  - Confirmed in trans with a P/LP partner (phased 0|1 vs 1|0, or parental phasing confirming trans) → +1.0 per confirmed trans (counts toward PM3 levels).
  - Unphased unconfirmed pairings contribute smaller amounts (0.5 for certain unphased P/LP) and are capped for unphased contributions (max 1.0 aggregated unphased).
  - Homozygous P/LP contributes to recessive evidence (counts as partial).
- PM3 levels (aggregated):
  - total >= 4.0 → PM3_VeryStrong (+8)
  - total >= 2.0 → PM3_Strong (+4)
  - total >= 1.0 → PM3_Moderate (+2)
  - total >= 0.5 → PM3_Supporting (+1)
- PM3 is downgraded to supporting if parental phasing evidence is absent for otherwise high raw scores (policy in apply_pm3_batch).
- Sources of trans confirmation:
  - proband phased GT (0|1 & 1|0),
  - parental phasing patterns (one variant inherited from father and other from mother),
  - cooccurrence evidence (provided TSV or DB lookup; pipeline accepts `--gnomad-cooccurrence` table).

### PP1 / PS2 (Segregation / De novo)
- PS2 (confirmed de novo) assigned when proband genotype indicates de novo and parental genotypes support absence in parents; strength depends on confirmation.
- PP1: segregation evidence used conservatively; strengths (supporting/moderate/strong) may require manual curation.

### PM1 (Mutational hotspot / functional domain)
- If variant lies in a documented mutational hotspot or critical domain per gene rules → PM1 (+2).

### HGMD / Internal DB / Gene validation
- HGMD DM with ≥ 2 publications or DM support count ≥ 2 is considered high-quality pathogenic evidence (treated similarly to high-quality ClinVar for PP5).
- Internal DB entries are included in database evidence but not automatically considered HQ unless marked/curated as such.
- Gene name validation: DB entries are validated (normalize_gene_name) to avoid spuriously applying DB evidence for mismatched gene annotations.

### Conflicts handling
- If HQ sources disagree between Pathogenic and Benign (substantive conflict) → variant is marked as HQ conflict and routed to manual review.
- Mild disagreements (Pathogenic vs Likely pathogenic) may be resolved conservatively towards Pathogenic for auto-assignment.

---

## Triage logic: automatic vs manual vs filtered

The triage aims to:
- Place only high-confidence results into `auto_conclusions.csv`
- Route uncertain or complex situations to `manual_review_list.csv`
- Filter clear carriers or low-quality calls out (remain recorded in `all_candidates.csv` and audit outputs)

Major rules:

1. Strict QC filtering
   - Low coverage or low allelic fraction will mark variant as failed QC (kept for auditing but not for reporting).
   - Heterozygous calls require DP ≥ STRICT_MIN_DP and AF ≥ STRICT_MIN_AF (defaults above). Homozygous/hemizygous thresholds are more permissive.

2. Carrier suppression (important protective rule)
   - For a sample/gene where MOI is recessive (AR/XLR) and the sample has exactly:
     - 1 heterozygous P/LP,
     - NO homozygous/hemizygous P/LP,
     - NO VUS with total_points ≥ 5,
     → then all variants in that gene for that sample are marked as filtered (treated as carrier) and are NOT escalated to manual or auto. This prevents "naked het" from being elevated.

3. AR Pairing & Rescue
   - 1 P/LP het + 1 VUS(≥5): both are routed to manual (phasing/pairing review).
   - 1 P/LP het + multiple VUS≥5: all relevant variants are routed to manual.
   - ≥2 P/LP hets:
     - If in trans (phased or parental pattern or cooccurrence) → manual review (and may become auto if other criteria satisfy).
     - If phasing unknown → manual (phasing required).

4. HQ DB-driven automatic assignment
   - If variant has unanimous HQ evidence (ClinVar ≥ 2★ Pathogenic, or HGMD DM with strong support) with no substantive conflict and MOI/zygosity consistent → assign automatically (auto), add PP5_Supporting and output to `auto_conclusions.csv`. AR heterozygotes remain blocked until pairing confirmed.

5. Algorithmic scoring automatic assignment
   - If total points >= 10 (Pathogenic) or >= 6 (Likely pathogenic) AND `MIN_PATHOGENIC_CRITERIA` satisfied → variant becomes auto (unless DB conflict or AR-het blocking applies).
   - Otherwise variant remains VUS or is routed to manual if flagged.

6. Delins / complex indels
   - Adjacent variants within 10 bp cluster. If phased evidence shows cis (same haplotype) → cluster considered delins and sent to manual for evaluation. If trans → not considered delins. If unphased and cluster meets P/LP thresholds → manual.

7. XL male handling
   - Male heterozygotes on chrX with AF <= 0.9 → manual for hemizygosity check. If AF >= 0.9 treat as hemizygous.

8. DB gene-mismatch
   - If DB evidence originates from a different gene (gene mismatch) that evidence is ignored and a note is added. If both HGMD and ClinVar mismatch disease annotation, escalate to manual (with carrier suppression guard for AR single hets).

---

## Outputs

- `all_candidates.csv` — full list of candidates passing gene-level filters and QC (with per-criterion explanations, DB summary, disease match metadata).
- `auto_conclusions.csv` — variants selected for automatic reporting (Pathogenic / Likely pathogenic) with reasons.
- `manual_review_list.csv` — variants flagged for manual curation (AR pairs, delins, HQ conflicts, ambiguous PM3 evidence, phenotype mismatches).
- `run_info.json` — run metadata (timestamps, counts, file paths).

Columns include variant identifiers, gene, HGVSc/HGVSp, REVEL/CADD/SpliceAI values, ClinVar/HGMD/internal DB summary, criteria assigned and per‑criterion textual explanations, total points and automated class, and triage flags (`in_auto_report`, `in_manual_review`, `filtered_out`).

---
*Developed for clinical NGS workflows at RCMG ("Genome" Centre).*