# ACMG Secondary Findings Classifier

Overview
--------
This repository contains a production-oriented implementation of an automated-first ACMG Secondary Findings (SF) classifier.

Goals
- Maximise safe automation: automatically classify and report variants where evidence is robust and conforms to decision rules.
- Minimise manual workload for clinical geneticists: use high-quality DB evidence and conservative computational thresholds to push variants to automated conclusions where appropriate.
- Preserve safety: gene-specific exceptions (e.g., TTN meta-exon, genes where LOF is not the disease mechanism) still trigger manual review.

Contents
--------
- `acmg_sf_clinical_pipeline.py`
  - Main pipeline script. Performs VEP annotation (if necessary), applies gene-specific rules and ACMG criteria, consults databases, and writes outputs.
- `make_ttn_meta_csv.py`
  - Utility to create a TTN meta-exon CSV used by the classifier. Normalises PSI values and produces an ordered CSV with `chr,start,end,exon_id,psi_meta,lof_allowed`.
- `README.md`
  - This document.

Requirements
------------
- Python 3.8+
- pip packages:
  - pandas
  - pysam
  - cyvcf2
  - intervaltree (only required if using TTN meta-exon support)
- System tools:
  - VEP (local) + installed plugins (SpliceAI, REVEL, CADD, etc.) if you want the script to run VEP.
  - tabix (optional, for indexing VEP output)
- Databases (provide paths via --db-paths-json or defaults):
  - ClinGen gene validity CSV
  - ClinGen variant assertions CSV
  - ClinVar VCF (recommended; required for ClinVar-based evidence)
  - HGMD Pro VCF (optional)
  - gnomAD VCF (optional)
  - INTERNAL_DB CSV (optional)
- The ACMG SF TSV (ACMG_SF_v3.3.tsv for example")

Design & Rules Summary
----------------------
- PVS1: Loss-of-function decision tree implemented, using:
  - gene LOF mechanism confirmation via ClinGen validity or TSV guidance
  - NMD tags (from VEP CSQ)
  - canonical splice site handling with SpliceAI >= 0.20
  - start-loss heuristics and last-exon handling (PM4)
  - TTN meta-exon handling using provided TTN meta-exon CSV
- PS1/PM5:
  - Uses ClinGen variant assertions, Internal DB, ClinVar protein-level evidence (indexed and cached), and HGMD DM (>=2 PMIDs)
  - ClinVar protein-level matching requires review status >=2★ for inclusion
- PP3/BP4:
  - Missense: REVEL thresholds (support 0.644, moderate 0.932), MutPred2 optional
  - Splicing: SpliceAI single threshold 0.20 (>= 0.2 supports splice impact)
  - CADD >=25 recorded as informative
- PM3:
  - Batch processing across variants in same sample/gene to determine compound-heterozygosity evidence
  - Uses parental genotypes for in-trans determination where available; otherwise assigns supporting/moderate PM3
- PS2/PM6:
  - De novo calling uses parental genotypes and allele-balance (AB) + depth (DP) heuristics to reduce false positives
- DB-driven override:
  - If algorithmic scoring yields VUS but a high-quality DB (ClinGen, ClinVar >=2★, HGMD DM>=2 pubs, Internal P/LP) contains P/LP evidence, script can automatically set final classification to DB classification to reduce manual workload. Enabled by --aggressive flag.

ClinVar protein-level index & cache
----------------------------------
- Building an index from ClinVar for protein-level matches (p.HGVS) is costly.
- The pipeline builds this index once and stores it as a gzipped JSON cache (path configured in DB_PATHS or db_paths JSON under "CLINVAR_INDEX_CACHE").
- On subsequent runs, the cache is loaded if the ClinVar VCF mtime equals cached mtime, saving significant startup time.
- If ClinVar VCF is updated, cache is automatically detected as stale and rebuilt.

Usage
-----
1. Prepare inputs:
   - Annotated VCF with VEP CSQ OR a raw VCF and local VEP available.
   - ACMG SF TSV (gene-specific guidance file).
   - Optional parental VCFs for trios, TTN meta CSV, and DB files.
2. Optionally prepare a db_paths JSON to override defaults (including CLINVAR_INDEX_CACHE path).
3. Run the classifier:
   ```
   python acmg_sf_clinical_pipeline.py \
     --proband sample.vcf.gz \
     --acmg-tsv "ACMG_SF_v3.3.tsv" \
     --outdir results \
     --vep /path/to/vep \
     --vep-cache /path/to/vep/cache \
     --fasta /path/to/human.fasta \
     --vep-extra "--plugin" "SpliceAI,snv,0.2" "--plugin" "REVEL" \
     --ttn-meta ttn_meta.csv \
     --db-paths-json db_paths.json \
     --father father.vcf.gz --mother mother.vcf.gz \
     --aggressive
   ```

Outputs
-------
- `all_candidates.csv` — comprehensive table of candidates with annotations, criteria, and scoring details.
- `auto_conclusions.csv` — variants assigned automated conclusions (Pathogenic/Likely pathogenic) and not requiring review.
- `manual_review_list.csv` — minimized list of variants that require manual inspection (gene-specific exceptions, TTN outside meta-exon, homozygous AR checks, ambiguous phasing).
- `run_info.json` — metadata about the run (inputs, parameters, counts).

TTN meta-exon CSV creation
--------------------------
Use `make_ttn_meta_csv.py` to produce the TTN meta CSV:
```
python make_ttn_meta_csv.py --input ttn_input.tsv --output ttn_meta.csv
```
Input must contain chr,start,end and a psi column (name containing "psi"). PSI may be 0-1 or 0-100.

License & Clinical Cote
-----------------------
This code is a clinical tool template, it is validated and integrated into workflow of the "Genome" Shared Resource Centre at the Research Centre for Medical Genetics (Moscow, Russia).

All rights reserved + Patent Pending
