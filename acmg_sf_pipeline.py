#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACMG Secondary Findings Classifier

This script:
- Annotates with local VEP if needed (user must have VEP + plugins configured).
- Applies gene-specific rules from ACMG SF TSV.
- Implements PVS1, PS1/PM5 (including ClinVar protein-level matching with cache), PP3/BP4,
  PM3 batch, PS2/PM6 (allele balance), TTN meta-exon handling.
- Minimizes manual_review flags; only P/LP variants are returned as secondary findings.
- Outputs:
    - all_candidates.csv  (all candidates, auditing)
    - auto_conclusions.csv (automated P/LP conclusions only)
    - manual_review_list.csv (P/LP variants that require manual review)
    - run_info.json

Author: 
Clinical Bioinformatics Team of the Research Centre for Medical Genetics (Code by Anna Basova)

"""

from __future__ import annotations
import os
import sys
import argparse
import json
import gzip
import logging
import subprocess
import time
import re
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any
from collections import defaultdict
import pandas as pd
import pysam
from cyvcf2 import VCF
try:
    from intervaltree import IntervalTree
except Exception:
    IntervalTree = None  # will error later if TTN required

# ---------------------------
# Configuration / constants
# ---------------------------
DB_PATHS_DEFAULT = {
    "CLINGEN_VALIDITY": "databases/clingen_gene_disease.csv",
    "CLINGEN_VARIANTS": "databases/clingen_variant_pathogenicity.csv",
    "INTERNAL_DB": "databases/internal_db.csv",
    "CLINVAR_VCF": "databases/clinvar.vcf.gz",
    # optional, used if you have access
    "HGMD_VCF": "databases/hgmd_pro.vcf.gz",
    "GNOMAD_VCF": "databases/gnomad.vcf.gz",
    "CLINVAR_INDEX_CACHE": "databases/clinvar_protein_index.json.gz"
}

# thresholds (tunable)
THRESH = {
    "REVEL_SUPPORT": 0.644,
    "REVEL_MODERATE": 0.932,
    "MUTPRED_SUPPORT": 0.737,
    "MUTPRED_MODERATE": 0.829,
    "SPLICEAI": 0.20,
    "BA1_AF": 0.05,
    "BS1_AF": 0.01,
    "PM2_STRICT": 1e-5
}

# Points mapping (ACGS 2024)
POINTS_MAP = {
    "PVS1_VeryStrong": 8, "PVS1_Strong": 4, "PVS1_Moderate": 2, "PVS1_Supporting": 1,
    "PS1_Strong": 4, "PS1_Moderate": 2, "PS1_Supporting": 1,
    "PS2_VeryStrong": 8, "PS2_Moderate": 2,
    "PM1": 2, "PM2_Moderate": 2, "PM2_Supporting": 1,
    "PM3_Supporting": 1, "PM3_Moderate": 2, "PM3_Strong": 4, "PM3_VeryStrong": 8,
    "PM4_Moderate": 2, "PM4_Supporting": 1, "PM5_Moderate": 2, "PM5_Supporting": 1,
    "PP1_Supporting": 1, "PP1_Moderate": 2, "PP1_Strong": 4,
    "PP3_Moderate": 2, "PP3_Supporting": 1,
    # benign
    "BA1": -8, "BS1": -4, "BS2": -4, "BP4": -1, "BP7": -1
}

CLASSIFICATION_THRESHOLDS = [
    ("Pathogenic", 10, None),
    ("Likely pathogenic", 6, 9),
    ("VUS", 0, 5),
    ("Likely benign", -5, -1),
    ("Benign", None, -6)
]

LOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
CANONICAL_SPLICE = {"splice_acceptor_variant", "splice_donor_variant"}

MIN_DP_FOR_PS2 = 15
MIN_AB_FOR_HET = 0.20

# ---------------------------
# Logging
# ---------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("acmg_sf_pipeline")

# ---------------------------
# Data model
# ---------------------------
@dataclass
class VariantRecord:
    """
    Lightweight container for per-variant data used throughout the pipeline.
    Only store fields required for scoring and output.
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    sample: str

    gene: str = ""
    transcript: str = ""
    consequence: str = ""
    hgvsc: str = ""
    hgvsp: str = ""
    exon: str = ""
    is_last_exon: bool = False

    # population / computational
    gnomad_af: Optional[float] = None
    gnomad_hom: Optional[int] = None
    revel: Optional[float] = None
    mutpred2: Optional[float] = None
    spliceai: Optional[float] = None
    cadd: Optional[float] = None
    nmd: Optional[str] = ""

    # genotype info
    proband_gt: str = "./."
    father_gt: str = "./."
    mother_gt: str = "./."
    proband_dp: Optional[int] = None
    proband_ad: Optional[Tuple[int,int]] = None
    proband_ab: Optional[float] = None

    db_hits: List[Tuple[str,str]] = field(default_factory=list)  # (classification, source)
    criteria_assigned: List[str] = field(default_factory=list)
    criteria_points: Dict[str,int] = field(default_factory=dict)
    total_points: int = 0
    automated_class: str = "VUS"
    manual_review: bool = False
    manual_reasons: List[str] = field(default_factory=list)

# ---------------------------
# Helper functions
# ---------------------------
def run_cmd(cmd: List[str], check=True) -> subprocess.CompletedProcess:
    """
    Run a system command, log stdout/stderr on failure.
    Returns CompletedProcess.
    """
    logger.debug("Running command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if check and proc.returncode != 0:
        logger.error("Command failed (%d): %s\nstdout:\n%s\nstderr:\n%s", proc.returncode, " ".join(cmd), proc.stdout, proc.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return proc

def vep_available(vep_cmd: str = "vep") -> bool:
    try:
        p = subprocess.run([vep_cmd, "--version"], capture_output=True, text=True, timeout=10)
        return p.returncode == 0
    except Exception:
        return False

def vcf_has_csq(vcf_path: str) -> bool:
    """
    Quick check for CSQ header in VCF using pysam.
    """
    try:
        vf = pysam.VariantFile(vcf_path)
        return vf.header.info.get("CSQ") is not None
    except Exception:
        return False

def parse_csq_header_from_vcf(vcf_path: str) -> List[str]:
    """
    Extract the CSQ field order from VEP-annotated VCF header.
    """
    vf = pysam.VariantFile(vcf_path)
    csq = vf.header.info.get("CSQ")
    if csq and csq.description and "Format:" in csq.description:
        header = csq.description.split("Format:")[1].strip()
        return [h.strip() for h in header.split("|")]
    for rec in vf.header.records:
        if rec.key == "INFO" and rec.get("ID") == "CSQ":
            desc = rec.get("Description") or ""
            if "Format:" in desc:
                header = desc.split("Format:")[1].strip()
                return [h.strip() for h in header.split("|")]
    return []

def parse_vep_csq_line(line: str, header: List[str]) -> Dict[str,str]:
    parts = line.split("|")
    out = {}
    for i, key in enumerate(header):
        out[key] = parts[i] if i < len(parts) else ""
    return out

def parse_protein_pos(hgvsp: str) -> Optional[int]:
    if not hgvsp:
        return None
    m = re.search(r'p\.[A-Za-z]{1,3}(\d+)', hgvsp)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return None
    return None

# ---------------------------
# Gene-specific rules loader
# ---------------------------
def load_acmg_tsv(tsv_path: str) -> Dict[str, Dict[str,Any]]:
    """
    Load ACMG SF TSV and derive common per-gene flags.

    The TSV contains many columns and gene-specific guidance.
    We parse out fields that affect automated decision rules:
      - reportable: whether gene is reportable as a secondary finding (Reportable as SF == Yes)
      - moi/inheritance
      - variants_to_report, reporting_guidance, mechanism
      - special_flags derived heuristically:
         require_biallelic, report_truncating_only, report_missense_only,
         lof_not_reportable, truncating_only
    The rules dict is keyed by gene symbol.
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    rules: Dict[str, Dict[str,Any]] = {}
    for _, row in df.iterrows():
        gene = row.get("Gene Symbol") or row.get("Gene") or ""
        if not gene:
            continue
        gene = str(gene).strip()
        reportable_field = row.get("Reportable as SF") or ""
        reportable = str(reportable_field).strip().lower() == "yes"
        moi = row.get("MOI") or row.get("Inheritance") or ""
        variants_to_report = row.get("Variants to report") or ""
        reporting_guidance = row.get("Reporting Guidance Comment") or ""
        mechanism = row.get("Механизм") or row.get("Mechanism") or ""
        guidance_text = " ".join([str(reporting_guidance or ""), str(mechanism or ""), str(variants_to_report or "")]).lower()
        special_flags = {
            "require_biallelic": bool("biallelic" in guidance_text or "2 variants" in variants_to_report.lower() or "2 het" in variants_to_report.lower()),
            "report_truncating_only": ("truncating" in guidance_text and "only" in guidance_text) or ("truncating variants only" in variants_to_report.lower()),
            "report_missense_only": ("missense" in guidance_text and ("only" in guidance_text or "только" in guidance_text)),
            # treat common GOF phrases as forbidding LOF
            "lof_not_reportable": any(x in guidance_text for x in ("gain-of-function", "gof", "go f", "go f")) or ("gof" in guidance_text),
            "truncating_only": "truncating variants only" in guidance_text or "truncating" in variants_to_report.lower()
        }
        rules[gene] = {
            "reportable": reportable,
            "moi": moi,
            "variants_to_report": variants_to_report,
            "reporting_guidance": reporting_guidance,
            "mechanism": mechanism,
            "special_flags": special_flags,
            "raw_row": row.to_dict()
        }
    logger.info("Loaded gene-specific rules for %d genes", len(rules))
    return rules

# ---------------------------
# Database manager with ClinVar protein-index caching
# ---------------------------
class DatabaseManager:
    """
    Loads ClinGen, ClinVar, HGMD (optional), gnomAD and internal DB (optional).
    Builds and caches a ClinVar protein-level index (gene -> list of protein annotations) to avoid repeated indexing.
    Provide "CLINVAR_INDEX_CACHE" path in db_paths to enable cache writing/loading.
    """
    def __init__(self, db_paths: Dict[str,str]):
        self.paths = dict(db_paths or {})
        self.clingen_validity = pd.DataFrame()
        self.clingen_variants = pd.DataFrame()
        self.internal_db = pd.DataFrame()
        self.clinvar_vcf: Optional[pysam.VariantFile] = None
        self.hgmd_vcf: Optional[pysam.VariantFile] = None
        self.gnomad_vf: Optional[pysam.VariantFile] = None
        self.clinvar_protein_index: Dict[str, List[dict]] = defaultdict(list)
        self._load_all()

    # ---------------------------
    # Cache helpers
    # ---------------------------
    def _clinvar_mtime(self, clinvar_path: str) -> Optional[float]:
        try:
            return float(os.path.getmtime(clinvar_path))
        except Exception:
            return None

    def _load_clinvar_index_cache(self, cache_path: str, clinvar_path: str) -> Optional[dict]:
        """
        Load gzipped JSON cache and validate mtime of clinvar VCF.
        Return index dict or None if stale/unusable.
        """
        if not os.path.exists(cache_path):
            return None
        try:
            with gzip.open(cache_path, "rt", encoding="utf-8") as fh:
                payload = json.load(fh)
            cached_mtime = payload.get("clinvar_mtime")
            current_mtime = self._clinvar_mtime(clinvar_path)
            if cached_mtime is None or current_mtime is None:
                logger.info("ClinVar cache present but cannot validate mtime -> rebuild")
                return None
            if float(cached_mtime) != float(current_mtime):
                logger.info("ClinVar cache stale (cached=%s current=%s) -> rebuild", cached_mtime, current_mtime)
                return None
            index = payload.get("index", {})
            logger.info("Loaded ClinVar protein-level index from cache: %s", cache_path)
            return index
        except Exception as e:
            logger.warning("Failed to load clinvar index cache %s: %s", cache_path, e)
            return None

    def _save_clinvar_index_cache(self, cache_path: str, clinvar_path: str, index: dict):
        """
        Write gzipped JSON cache with clinvar mtime and index dictionary.
        Atomic write to avoid partial files.
        """
        try:
            mtime = self._clinvar_mtime(clinvar_path)
            payload = {"clinvar_mtime": mtime, "index": index}
            tmp = cache_path + ".tmp"
            with gzip.open(tmp, "wt", encoding="utf-8") as fh:
                json.dump(payload, fh, ensure_ascii=False)
            os.replace(tmp, cache_path)
            logger.info("Wrote ClinVar protein index cache to %s (mtime=%s)", cache_path, mtime)
        except Exception as e:
            logger.warning("Failed to write clinvar index cache %s: %s", cache_path, e)
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass

    # ---------------------------
    # Load all DBs and optionally build or load clinvar index
    # ---------------------------
    def _load_all(self):
        # ClinGen validity
        p = self.paths.get("CLINGEN_VALIDITY")
        if p and os.path.exists(p):
            try:
                self.clingen_validity = pd.read_csv(p, dtype=str)
                logger.info("Loaded ClinGen validity")
            except Exception as e:
                logger.warning("Failed to load ClinGen validity: %s", e)

        # ClinGen variants
        p = self.paths.get("CLINGEN_VARIANTS")
        if p and os.path.exists(p):
            try:
                self.clingen_variants = pd.read_csv(p, dtype=str)
                logger.info("Loaded ClinGen variants")
            except Exception as e:
                logger.warning("Failed to load ClinGen variants: %s", e)

        # Internal DB
        p = self.paths.get("INTERNAL_DB")
        if p and os.path.exists(p):
            try:
                self.internal_db = pd.read_csv(p, dtype=str)
                if 'chrom' in self.internal_db.columns and 'pos' in self.internal_db.columns:
                    def make_key(r):
                        try:
                            pos = int(r['pos'])
                        except Exception:
                            pos = r.get('pos')
                        return f"{r.get('chrom')}-{pos}-{r.get('ref')}-{r.get('alt')}"
                    self.internal_db['key'] = self.internal_db.apply(make_key, axis=1)
                logger.info("Loaded Internal DB")
            except Exception as e:
                logger.warning("Failed to load Internal DB: %s", e)

        # ClinVar VCF and protein index
        p = self.paths.get("CLINVAR_VCF")
        cache_path = self.paths.get("CLINVAR_INDEX_CACHE")
        if p and os.path.exists(p):
            try:
                self.clinvar_vcf = pysam.VariantFile(p)
                logger.info("Opened ClinVar VCF: %s", p)
            except Exception as e:
                logger.warning("Failed to open ClinVar VCF: %s", e)
                self.clinvar_vcf = None
            # attempt to load index from cache
            loaded_index = None
            if self.clinvar_vcf and cache_path:
                try:
                    loaded_index = self._load_clinvar_index_cache(cache_path, p)
                except Exception:
                    loaded_index = None
            if loaded_index:
                self.clinvar_protein_index = defaultdict(list, loaded_index)
            else:
                # build index and cache if possible
                if self.clinvar_vcf:
                    try:
                        index = self._index_clinvar_protein_level(build_only=True)
                        self.clinvar_protein_index = defaultdict(list, index)
                        if cache_path:
                            try:
                                self._save_clinvar_index_cache(cache_path, p, index)
                            except Exception:
                                logger.warning("Failed to save clinvar protein index cache")
                    except Exception as e:
                        logger.warning("Failed to build ClinVar protein index: %s", e)
        else:
            logger.info("ClinVar VCF not provided or not found; ClinVar-based protein matching disabled")

        # HGMD
        p = self.paths.get("HGMD_VCF")
        if p and os.path.exists(p):
            try:
                self.hgmd_vcf = pysam.VariantFile(p)
                logger.info("Opened HGMD VCF")
            except Exception as e:
                logger.warning("Failed to open HGMD VCF: %s", e)

        # gnomAD
        p = self.paths.get("GNOMAD_VCF")
        if p and os.path.exists(p):
            try:
                self.gnomad_vf = pysam.VariantFile(p)
                logger.info("Opened gnomAD VCF")
            except Exception:
                logger.info("gnomAD provided but not readable; rely on CSQ AFs")

    def _index_clinvar_protein_level(self, build_only: bool = False) -> dict:
        """
        Build an index mapping gene -> list of protein-level ClinVar entries.

        Each entry is a dict: {"hgvsp": str, "clnsig": str, "stars": int, "chrom": str, "pos": int}

        This operation is performed once and cached to reduce startup cost for repeated runs.
        """
        if not self.clinvar_vcf:
            return {}

        index = defaultdict(list)
        for rec in self.clinvar_vcf.fetch():
            geneinfo = rec.info.get('GENEINFO') or rec.info.get('Gene') or rec.info.get('GENE')
            if not geneinfo:
                continue
            if isinstance(geneinfo, (list, tuple)):
                geneinfo = geneinfo[0]
            gene = str(geneinfo).split(':')[0] if ':' in str(geneinfo) else str(geneinfo)

            # compute review stars using CLNREVSTAT if present
            rev = rec.info.get('CLNREVSTAT') or ''
            revs = ";".join(rev) if isinstance(rev, (list,tuple)) else str(rev).lower()
            stars = 0
            if 'practice_guideline' in revs:
                stars = 4
            elif 'expert_panel' in revs:
                stars = 3
            elif 'criteria_provided' in revs and 'multiple_submitters' in revs and ('no_conflicts' in revs or 'no_conflict' in revs):
                stars = 2
            elif 'criteria_provided' in revs:
                stars = 1

            clnsig = rec.info.get('CLNSIG') or ''
            clnsig_str = ";".join(clnsig) if isinstance(clnsig, (list,tuple)) else str(clnsig)

            # scan INFO fields for p.<...> patterns to extract protein-level HGVS
            p_candidates = set()
            for key in rec.info.keys():
                try:
                    val = rec.info.get(key)
                    if not val:
                        continue
                    vals = val if isinstance(val, (list,tuple)) else [val]
                    for v in vals:
                        if not isinstance(v, str):
                            v = str(v)
                        for m in re.finditer(r'p\.[A-Za-z0-9_\-\*\(\)]+', v):
                            token = m.group(0).rstrip(';,')
                            p_candidates.add(token)
                except Exception:
                    continue

            for p_hgvs in p_candidates:
                entry = {
                    "hgvsp": p_hgvs,
                    "clnsig": clnsig_str,
                    "stars": stars,
                    "chrom": getattr(rec, "chrom", None),
                    "pos": getattr(rec, "pos", None)
                }
                index[gene].append(entry)
        # convert to plain dict (json-serializable) when returning
        return {k: v for k, v in index.items()}

    # ---------------------------
    # Public helpers used by scoring engine
    # ---------------------------
    def get_gene_mechanism(self, gene_symbol: str) -> str:
        """
        Return "Valid" if ClinGen validity table marks this gene with definitive/strong LOF mechanism,
        otherwise "Disputed" or "Unknown".
        """
        if self.clingen_validity.empty:
            return "Unknown"
        try:
            lowcols = [c.lower() for c in self.clingen_validity.columns]
            if 'gene_symbol' in lowcols:
                row = self.clingen_validity[self.clingen_validity['gene_symbol'] == gene_symbol]
                if not row.empty:
                    classification = str(row.iloc[0].get('classification','')).lower()
                    if classification in ('definitive','strong'):
                        return "Valid"
                    return "Disputed"
            row = self.clingen_validity[self.clingen_validity.iloc[:,0].astype(str).str.contains(gene_symbol, na=False)]
            if row.empty:
                return "Unknown"
            classification = str(row.iloc[0].get('classification','')).lower()
            if classification in ('definitive','strong'):
                return "Valid"
            return "Disputed"
        except Exception:
            return "Unknown"

    def get_variant_db_evidence(self, vr: VariantRecord) -> List[Tuple[str,str]]:
        """
        Return list of high-quality DB evidence tuples (classification_label, source_name)
        Sources considered:
          - ClinGen variant-level assertions (hgvs_c or hgvs_p)
          - Internal DB exact coordinate matches (classification)
          - ClinVar coordinate-level P/LP with review >=2★
          - HGMD DM with >=2 PMIDs
        """
        results: List[Tuple[str,str]] = []

        # ClinGen variant-level
        if not self.clingen_variants.empty:
            try:
                if vr.hgvsc:
                    hits = self.clingen_variants[self.clingen_variants.get('hgvs_c','') == vr.hgvsc]
                    for _, r in hits.iterrows():
                        results.append((r.get('pathogenicity') or r.get('classification') or "", "ClinGen"))
                if vr.hgvsp:
                    hits = self.clingen_variants[self.clingen_variants.get('hgvs_p','') == vr.hgvsp]
                    for _, r in hits.iterrows():
                        results.append((r.get('pathogenicity') or r.get('classification') or "", "ClinGen"))
            except Exception:
                pass

        # Internal DB
        if not self.internal_db.empty and 'key' in self.internal_db.columns:
            try:
                key = f"{vr.chrom}-{vr.pos}-{vr.ref}-{vr.alt}"
                hits = self.internal_db[self.internal_db['key'] == key]
                for _, r in hits.iterrows():
                    results.append((r.get('classification') or "", "InternalDB"))
            except Exception:
                pass

        # ClinVar coordinate-level with >=2★
        if self.clinvar_vcf:
            try:
                for rec in self.clinvar_vcf.fetch(vr.chrom, vr.pos - 1, vr.pos):
                    if getattr(rec, 'ref', None) != vr.ref:
                        continue
                    if not (rec.alts and vr.alt in rec.alts):
                        continue
                    rev = rec.info.get('CLNREVSTAT') or ''
                    revs = ";".join(rev) if isinstance(rev, (list,tuple)) else str(rev).lower()
                    stars = 0
                    if 'practice_guideline' in revs:
                        stars = 4
                    elif 'expert_panel' in revs:
                        stars = 3
                    elif 'criteria_provided' in revs and 'multiple_submitters' in revs and ('no_conflicts' in revs or 'no_conflict' in revs):
                        stars = 2
                    elif 'criteria_provided' in revs:
                        stars = 1
                    if stars >= 2:
                        clnsig = rec.info.get('CLNSIG') or ''
                        s = ";".join(clnsig) if isinstance(clnsig, (list,tuple)) else str(clnsig)
                        s = s.lower()
                        if 'pathogenic' in s:
                            results.append(("Pathogenic","ClinVar"))
                        elif 'likely_pathogenic' in s or 'likely pathogenic' in s:
                            results.append(("Likely pathogenic","ClinVar"))
            except Exception:
                pass

        # HGMD DM with >=2 PMIDs
        if self.hgmd_vcf:
            try:
                for rec in self.hgmd_vcf.fetch(vr.chrom, vr.pos - 1, vr.pos):
                    if getattr(rec,'ref', None) != vr.ref:
                        continue
                    if not (rec.alts and vr.alt in rec.alts):
                        continue
                    hgmd_class = rec.info.get('CLASS') or rec.info.get('HGMD_CLASS') or ""
                    pmids = rec.info.get('PMID') or rec.info.get('PMIDs') or []
                    pmid_count = len(set(pmids)) if pmids else 0
                    if str(hgmd_class).upper() == 'DM' and pmid_count >= 2:
                        results.append(("Pathogenic","HGMD_Pro"))
            except Exception:
                pass

        return results

    def annotate_gnomad_for_variant(self, vr: VariantRecord):
        """
        Populate gnomAD AF/hom counts from provided gnomAD VCF if available.
        """
        if not self.gnomad_vf:
            return
        try:
            for rec in self.gnomad_vf.fetch(vr.chrom, vr.pos - 1, vr.pos):
                if getattr(rec,'ref', None) != vr.ref:
                    continue
                if not (rec.alts and vr.alt in rec.alts):
                    continue
                af = 0.0
                if 'AF' in rec.info:
                    a = rec.info.get('AF')
                    try:
                        af = float(a[0] if isinstance(a,(list,tuple)) else a)
                    except Exception:
                        af = 0.0
                else:
                    try:
                        ac = rec.info.get('AC'); an = rec.info.get('AN')
                        ac0 = ac[0] if isinstance(ac,(list,tuple)) else ac
                        an0 = an[0] if isinstance(an,(list,tuple)) else an
                        if ac0 is not None and an0:
                            af = float(ac0) / float(an0)
                    except Exception:
                        af = 0.0
                vr.gnomad_af = af
                try:
                    vr.gnomad_hom = int(rec.info.get('nhomalt') or 0)
                except Exception:
                    pass
                break
        except Exception:
            pass

    def find_protein_variant_by_hgvsp(self, gene: str, hgvsp: Optional[str], pos: Optional[int] = None) -> List[Tuple[str,str]]:
        """
        Return matches from ClinGen, InternalDB, and ClinVar protein-level index.
        ClinVar matches are accepted only when review stars >= 2.
        """
        results: List[Tuple[str,str]] = []

        # ClinGen
        if not self.clingen_variants.empty:
            try:
                if hgvsp:
                    hits = self.clingen_variants[self.clingen_variants.get('hgvs_p','') == hgvsp]
                    for _, r in hits.iterrows():
                        results.append((r.get('pathogenicity') or r.get('classification') or "", "ClinGen"))
                if pos:
                    if 'hgvs_p' in self.clingen_variants.columns:
                        try:
                            pos_series = self.clingen_variants['hgvs_p'].fillna('').astype(str).str.extract(r'p\.[A-Za-z]+(\d+)', expand=False)
                            mask = pos_series.dropna().astype(float) == float(pos)
                            sel = self.clingen_variants[mask.fillna(False)]
                            for _, r in sel.iterrows():
                                results.append((r.get('pathogenicity') or r.get('classification') or "", "ClinGen"))
                        except Exception:
                            pass
            except Exception:
                pass

        # Internal DB
        if not self.internal_db.empty:
            try:
                if 'hgvs_p' in self.internal_db.columns and hgvsp:
                    hits = self.internal_db[self.internal_db['hgvs_p'] == hgvsp]
                    for _, r in hits.iterrows():
                        results.append((r.get('classification') or "", "InternalDB"))
                if 'hgvs_p' in self.internal_db.columns and pos:
                    try:
                        pos_series = self.internal_db['hgvs_p'].fillna('').astype(str).str.extract(r'p\.[A-Za-z]+(\d+)', expand=False)
                        mask = pos_series.dropna().astype(float) == float(pos)
                        sel = self.internal_db[mask.fillna(False)]
                        for _, r in sel.iterrows():
                            results.append((r.get('classification') or "", "InternalDB"))
                    except Exception:
                        pass
            except Exception:
                pass

        # ClinVar protein-level index
        try:
            entries = self.clinvar_protein_index.get(gene) or []
            if hgvsp:
                for e in entries:
                    if e.get("hgvsp") and e["hgvsp"].strip().lower() == hgvsp.strip().lower() and (e.get("stars",0) >= 2):
                        clnsig = e.get("clnsig","")
                        if "pathogen" in str(clnsig).lower():
                            results.append(("Pathogenic","ClinVar"))
                        elif "likely" in str(clnsig).lower():
                            results.append(("Likely pathogenic","ClinVar"))
            if pos:
                for e in entries:
                    p = e.get("hgvsp")
                    if not p:
                        continue
                    m = re.search(r'p\.[A-Za-z]{1,3}(\d+)', p)
                    if not m:
                        continue
                    try:
                        ppos = int(m.group(1))
                    except Exception:
                        continue
                    if ppos == pos and (e.get("stars",0) >= 2):
                        clnsig = e.get("clnsig","")
                        if "pathogen" in str(clnsig).lower():
                            results.append(("Pathogenic","ClinVar"))
                        elif "likely" in str(clnsig).lower():
                            results.append(("Likely pathogenic","ClinVar"))
        except Exception:
            pass

        return results

# ---------------------------
# ACMG scoring engine
# ---------------------------
class ACMGEngine:
    """
    Implements the primary scoring logic. Designed to be conservative but automated-first:
    - Apply strong automated criteria when supported by databases and computational evidence.
    - Minimize manual_review flags; set them only for gene-specific critical cases (e.g., TTN outside meta-exon)
      or when phenotype/zygosity/phasing confirmation is required for safe reporting (e.g., homozygous AR).
    """
    def __init__(self, dbm: DatabaseManager, gene_rules: Dict[str, Dict[str,Any]], ttn_meta: Dict[str, IntervalTree] = None, spliceai_thr: float = 0.20):
        self.dbm = dbm
        self.gene_rules = gene_rules
        self.ttn_meta = ttn_meta or {}
        self.spliceai_thr = spliceai_thr

    def _get_gene_rule(self, gene: str) -> Dict[str,Any]:
        return self.gene_rules.get(gene, {"reportable": False, "special_flags": {}, "raw_row": {}})

    def _ttn_overlap(self, chrom: str, pos: int) -> Optional[dict]:
        chromn = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
        if chromn not in self.ttn_meta:
            return None
        hits = self.ttn_meta[chromn].at(pos - 1)
        if not hits:
            return None
        return max(hits, key=lambda it: float(it.data.get("psi_meta", 0) or 0)).data

    # ---- PVS1 decision tree ----
    def _apply_pvs1(self, vr: VariantRecord):
        rule = self._get_gene_rule(vr.gene)
        special = rule.get("special_flags", {})
        gene_mech = self.dbm.get_gene_mechanism(vr.gene)
        lof_allowed = True
        # apply local guidance and ClinGen gene mechanism
        if special.get("lof_not_reportable"):
            lof_allowed = False
        if gene_mech != "Valid" and not rule.get("mechanism"):
            lof_allowed = False

        cons = (vr.consequence or "").lower()
        is_lof = any(x in cons for x in LOF_CONSEQUENCES)
        if not is_lof or not lof_allowed:
            return

        # TTN handling: only apply PVS1 automatically if variant falls into high-PSI meta-exon and guidance allows
        if vr.gene.upper() == "TTN":
            ttnover = self._ttn_overlap(vr.chrom, vr.pos)
            if not ttnover or float(ttnover.get("psi_meta", 0) or 0) < 70 or str(ttnover.get("lof_allowed","1")).lower() in ("0","false","no"):
                # TTN interpretation is gene-specific and clinically sensitive: keep manual_review for safety
                vr.manual_reasons.append("TTN truncating variant outside high-PSI meta-exon -> PVS1 not applied automatically; manual review recommended")
                vr.manual_review = True
                return

        # Last exon -> PM4 (protein-length change)
        if vr.is_last_exon:
            vr.criteria_assigned.append("PM4_Moderate")
            vr.criteria_points["PM4_Moderate"] = POINTS_MAP.get("PM4_Moderate", 2)
            return

        # NMD tag (from VEP) if present
        nmd_tag = (vr.nmd or "").lower()
        if nmd_tag:
            if "triggers_nmd" in nmd_tag or "yes" in nmd_tag:
                vr.criteria_assigned.append("PVS1_VeryStrong")
                vr.criteria_points["PVS1_VeryStrong"] = POINTS_MAP.get("PVS1_VeryStrong", 8)
                return
            elif "escapes_nmd" in nmd_tag or "no" in nmd_tag:
                vr.criteria_assigned.append("PVS1_Moderate")
                vr.criteria_points["PVS1_Moderate"] = POINTS_MAP.get("PVS1_Moderate", 2)
                return

        # canonical splice site handling using unified SpliceAI threshold
        if any(x in cons for x in CANONICAL_SPLICE):
            if vr.spliceai is not None and vr.spliceai >= self.spliceai_thr:
                vr.criteria_assigned.append("PVS1_VeryStrong")
                vr.criteria_points["PVS1_VeryStrong"] = POINTS_MAP.get("PVS1_VeryStrong", 8)
                return
            else:
                vr.criteria_assigned.append("PVS1_Strong")
                vr.criteria_points["PVS1_Strong"] = POINTS_MAP.get("PVS1_Strong", 4)
                return

        # start-loss 
        if "start_lost" in cons:
            pos_aa = parse_protein_pos(vr.hgvsp) if vr.hgvsp else None
            if pos_aa and pos_aa > 20:
                vr.criteria_assigned.append("PVS1_Moderate")
                vr.criteria_points["PVS1_Moderate"] = POINTS_MAP.get("PVS1_Moderate", 2)
            else:
                vr.criteria_assigned.append("PVS1_VeryStrong")
                vr.criteria_points["PVS1_VeryStrong"] = POINTS_MAP.get("PVS1_VeryStrong", 8)
            return

        # frameshift / stop gained -> PVS1_VeryStrong (unless last exon handled earlier)
        if any(x in cons for x in ("frameshift_variant", "stop_gained")):
            vr.criteria_assigned.append("PVS1_VeryStrong")
            vr.criteria_points["PVS1_VeryStrong"] = POINTS_MAP.get("PVS1_VeryStrong", 8)
            return

    # ---- PS1 & PM5 (protein and splice equivalence) ----
    def _apply_ps1_pm5(self, vr: VariantRecord):
        pos_aa = parse_protein_pos(vr.hgvsp) if vr.hgvsp else None
        prot_matches = self.dbm.find_protein_variant_by_hgvsp(vr.gene, vr.hgvsp, pos_aa)
        if prot_matches:
            has_path = any(cls and 'pathogen' in str(cls).lower() for cls,_ in prot_matches if cls)
            has_likely = any(cls and 'likely' in str(cls).lower() for cls,_ in prot_matches if cls)
            if has_path:
                vr.criteria_assigned.append("PS1_Strong")
                vr.criteria_points["PS1_Strong"] = POINTS_MAP.get("PS1_Strong", 4)
                vr.manual_reasons.append("PS1_Strong: exact protein match to pathogenic variant in high-quality DB")
                return
            if has_likely:
                vr.criteria_assigned.append("PS1_Moderate")
                vr.criteria_points["PS1_Moderate"] = POINTS_MAP.get("PS1_Moderate", 2)
                vr.manual_reasons.append("PS1_Moderate: exact protein match to likely pathogenic variant in high-quality DB")
                return

        # PM5 same-residue different AA
        if pos_aa and ("missense" in (vr.consequence or "").lower() or "inframe" in (vr.consequence or "").lower()):
            residue_matches = self.dbm.find_protein_variant_by_hgvsp(vr.gene, None, pos_aa)
            if residue_matches:
                has_path = any(cls and 'pathogen' in str(cls).lower() for cls,_ in residue_matches if cls)
                has_likely = any(cls and 'likely' in str(cls).lower() for cls,_ in residue_matches if cls)
                if has_path:
                    vr.criteria_assigned.append("PM5_Moderate")
                    vr.criteria_points["PM5_Moderate"] = POINTS_MAP.get("PM5_Moderate", 2)
                    vr.manual_reasons.append("PM5_Moderate: residue has pathogenic missense in high-quality DB")
                    return
                if has_likely:
                    vr.criteria_assigned.append("PM5_Supporting")
                    vr.criteria_points["PM5_Supporting"] = POINTS_MAP.get("PM5_Supporting", 1)
                    vr.manual_reasons.append("PM5_Supporting: residue has likely pathogenic missense in high-quality DB")
                    return

        # splice-equivalence PS1 (SpliceAI + DB)
        if vr.spliceai is not None and vr.spliceai >= self.spliceai_thr:
            s_matches = self.dbm.find_splice_equivalent_variants(vr.gene, vr.pos, radius=6)
            if any(cls and ('pathogen' in str(cls).lower() or 'likely' in str(cls).lower()) for cls,_ in s_matches if cls):
                vr.criteria_assigned.append("PS1_Supporting")
                vr.criteria_points["PS1_Supporting"] = POINTS_MAP.get("PS1_Supporting", 1)
                vr.manual_reasons.append("PS1_Supporting: splice-equivalent P/LP variant found in high-quality DB")
                return

    # ---- PP3 / BP4 computational metrics ----
    def _apply_pp3_bp4(self, vr: VariantRecord):
        cons = (vr.consequence or "").lower()
        is_missense = any(x in cons for x in ("missense", "inframe_variant", "inframe_insertion", "inframe_deletion"))
        if is_missense:
            if vr.revel is not None:
                if vr.revel >= THRESH["REVEL_MODERATE"]:
                    vr.criteria_assigned.append("PP3_Moderate")
                    vr.criteria_points["PP3_Moderate"] = POINTS_MAP.get("PP3_Moderate", 2)
                elif vr.revel >= THRESH["REVEL_SUPPORT"]:
                    vr.criteria_assigned.append("PP3_Supporting")
                    vr.criteria_points["PP3_Supporting"] = POINTS_MAP.get("PP3_Supporting", 1)
                elif vr.revel < 0.4:
                    vr.criteria_assigned.append("BP4")
                    vr.criteria_points["BP4"] = POINTS_MAP.get("BP4", -1)
            if vr.mutpred2 is not None:
                if vr.mutpred2 >= THRESH["MUTPRED_MODERATE"]:
                    vr.criteria_assigned.append("PP3_Moderate")
                    vr.criteria_points["PP3_Moderate"] = POINTS_MAP.get("PP3_Moderate", 2)
                elif vr.mutpred2 >= THRESH["MUTPRED_SUPPORT"]:
                    vr.criteria_assigned.append("PP3_Supporting")
                    vr.criteria_points["PP3_Supporting"] = POINTS_MAP.get("PP3_Supporting", 1)
        # SpliceAI unified threshold
        if vr.spliceai is not None:
            if vr.spliceai >= self.spliceai_thr:
                vr.criteria_assigned.append("PP3_Supporting")
                vr.criteria_points["PP3_Supporting"] = POINTS_MAP.get("PP3_Supporting", 1)
            else:
                vr.criteria_assigned.append("BP4")
                vr.criteria_points["BP4"] = POINTS_MAP.get("BP4", -1)
        # CADD used as informative note only
        if vr.cadd is not None and vr.cadd >= 25:
            vr.manual_reasons.append(f"CADD {vr.cadd} >= 25 (informative high deleteriousness)")

    # ---- PS2 / PM6 (de novo) ----
    def _apply_ps2_pm6(self, vr: VariantRecord):
        pg = vr.proband_gt
        fg = vr.father_gt
        mg = vr.mother_gt
        dp = vr.proband_dp or 0
        ab = vr.proband_ab or 0.0
        if pg in ("0/1","1/0","0|1","1|0"):
            # robust de novo: both parents hom-ref, allele balance & depth meet thresholds
            if dp >= MIN_DP_FOR_PS2 and ab >= MIN_AB_FOR_HET and fg in ("0/0","0|0") and mg in ("0/0","0|0"):
                vr.criteria_assigned.append("PS2_VeryStrong")
                vr.criteria_points["PS2_VeryStrong"] = POINTS_MAP.get("PS2_VeryStrong", 8)
                return
            # parents negative but AB/DP not robust -> PM6 assumed de novo
            if fg in ("0/0","0|0") and mg in ("0/0","0|0"):
                vr.criteria_assigned.append("PM6")
                vr.criteria_points["PM6"] = 1
                return

    # ---- PM3 batch for AR genes (compound het) ----
    def apply_pm3_batch(self, candidates: List[VariantRecord]):
        by_sample_gene = defaultdict(lambda: defaultdict(list))
        for v in candidates:
            by_sample_gene[v.sample][v.gene].append(v)
        for sample, genes in by_sample_gene.items():
            for gene, varlist in genes.items():
                if len(varlist) == 0:
                    continue
                # homozygous rare -> PM3_Strong
                if len(varlist) == 1:
                    v = varlist[0]
                    if v.proband_gt in ("1/1","1|1"):
                        v.criteria_assigned.append("PM3_Strong")
                        v.criteria_points["PM3_Strong"] = POINTS_MAP.get("PM3_Strong", 4)
                        v.manual_reasons.append("PM3_Strong: homozygous proband (AR scenario); phenotype correlation recommended")
                        v.manual_review = True
                    continue
                # pairwise checks
                for i in range(len(varlist)):
                    for j in range(i+1, len(varlist)):
                        v1 = varlist[i]
                        v2 = varlist[j]
                        def is_het(gt): return gt in ("0/1","1/0","0|1","1|0")
                        in_trans = False
                        if is_het(v1.father_gt) and is_het(v2.mother_gt):
                            in_trans = True
                        if is_het(v2.father_gt) and is_het(v1.mother_gt):
                            in_trans = True
                        if in_trans:
                            for vv in (v1, v2):
                                vv.criteria_assigned.append("PM3_Strong")
                                vv.criteria_points["PM3_Strong"] = POINTS_MAP.get("PM3_Strong", 4)
                                vv.manual_reasons.append("PM3_Strong: likely in trans by parental genotypes")
                                vv.manual_review = True
                        else:
                            db_p1 = any(cls and ('pathogen' in str(cls).lower() or 'likely' in str(cls).lower()) for cls,_ in v1.db_hits) if v1.db_hits else False
                            db_p2 = any(cls and ('pathogen' in str(cls).lower() or 'likely' in str(cls).lower()) for cls,_ in v2.db_hits) if v2.db_hits else False
                            if db_p1 and db_p2:
                                for vv in (v1, v2):
                                    vv.criteria_assigned.append("PM3_Moderate")
                                    vv.criteria_points["PM3_Moderate"] = POINTS_MAP.get("PM3_Moderate", 2)
                                    vv.manual_reasons.append("PM3_Moderate: two P/LP variants in same gene, phasing unknown")
                                    vv.manual_review = True
                            else:
                                for vv in (v1, v2):
                                    vv.criteria_assigned.append("PM3_Supporting")
                                    vv.criteria_points["PM3_Supporting"] = POINTS_MAP.get("PM3_Supporting", 1)
                                    vv.manual_reasons.append("PM3_Supporting: two rare variants in same gene, phasing unknown")
                                    vv.manual_review = True

    # ---- classification helper ----
    def _classify_by_points(self, pts: int) -> str:
        for name, low, high in CLASSIFICATION_THRESHOLDS:
            if low is None and high is not None:
                if pts <= high:
                    return name
            elif high is None and low is not None:
                if pts >= low:
                    return name
            elif low is not None and high is not None and low <= pts <= high:
                return name
        return "VUS"

    # ---- score single variant (main entry) ----
    def score_variant(self, vr: VariantRecord):
        # clear old
        vr.criteria_assigned = []
        vr.criteria_points = {}
        vr.manual_reasons = []

        # gather DB evidence early; used for PM3 decisions below
        vr.db_hits = self.dbm.get_variant_db_evidence(vr)

        # population filters (BA1/BS1/PM2)
        af = vr.gnomad_af if vr.gnomad_af is not None else 0.0
        if af >= THRESH["BA1_AF"]:
            vr.criteria_assigned.append("BA1")
            vr.criteria_points["BA1"] = POINTS_MAP.get("BA1", -8)
            vr.total_points = sum(vr.criteria_points.values())
            vr.automated_class = self._classify_by_points(vr.total_points)
            vr.manual_reasons.append(f"BA1: population AF {af} >= {THRESH['BA1_AF']}")
            return
        if af >= THRESH["BS1_AF"]:
            vr.criteria_assigned.append("BS1")
            vr.criteria_points["BS1"] = POINTS_MAP.get("BS1", -4)
            vr.manual_reasons.append(f"BS1: population AF {af} >= {THRESH['BS1_AF']}")

        if af == 0.0:
            vr.criteria_assigned.append("PM2_Moderate")
            vr.criteria_points["PM2_Moderate"] = POINTS_MAP.get("PM2_Moderate", 2)
        elif 0.0 < af < THRESH["PM2_STRICT"]:
            vr.criteria_assigned.append("PM2_Supporting")
            vr.criteria_points["PM2_Supporting"] = POINTS_MAP.get("PM2_Supporting", 1)

        # gene-specific rules may exclude variant types from automated reporting
        rule = self._get_gene_rule(vr.gene)
        special = rule.get("special_flags", {})
        cons = (vr.consequence or "").lower()
        if special.get("report_missense_only") and "missense" not in cons:
            vr.manual_reasons.append("Gene rule: only missense variants reportable -> excluded from automated SF")
            vr.manual_review = True
            return
        if special.get("report_truncating_only") and not any(x in cons for x in LOF_CONSEQUENCES):
            vr.manual_reasons.append("Gene rule: only truncating variants reportable -> excluded from automated SF")
            vr.manual_review = True
            return

        # apply main criteria
        self._apply_pvs1(vr)
        self._apply_ps1_pm5(vr)
        self._apply_pp3_bp4(vr)
        self._apply_ps2_pm6(vr)

        # sum points
        total = 0
        for _, pts in vr.criteria_points.items():
            try:
                total += int(pts)
            except Exception:
                try:
                    total += float(pts)
                except Exception:
                    pass
        vr.total_points = int(total)
        vr.automated_class = self._classify_by_points(vr.total_points)

        # (??? still thinking) automation-first DB-driven override to reduce manual workload:
        # If algorithm yields VUS but there's high-quality DB P/LP evidence, set automated_class to DB classification.
        if vr.automated_class == "VUS" and vr.db_hits:
            priority = {"ClinGen": 4, "ClinVar": 3, "HGMD_Pro": 3, "InternalDB": 2}
            selected = None
            best_score = 0
            for cls, src in vr.db_hits:
                score = priority.get(src, 0)
                if score > best_score:
                    best_score = score
                    selected = (cls, src)
            if selected:
                cls, src = selected
                cls_str = str(cls).lower()
                if "pathogen" in cls_str:
                    vr.automated_class = "Pathogenic"
                    vr.manual_reasons.append(f"Automated override to Pathogenic ({src}) to reduce manual review")
                elif "likely" in cls_str:
                    vr.automated_class = "Likely pathogenic"
                    vr.manual_reasons.append(f"Automated override to Likely Pathogenic ({src}) to reduce manual review")

        # enforce gene-specific prohibition of LOF if present
        if special.get("lof_not_reportable"):
            removed = False
            for code in list(vr.criteria_points.keys()):
                if code.startswith("PVS1"):
                    vr.criteria_points.pop(code, None)
                    removed = True
            if removed:
                vr.criteria_assigned = [c for c in vr.criteria_assigned if not c.startswith("PVS1")]
                vr.manual_reasons.append("Removed PVS1 per gene-specific rule forbidding LOF")
                # update totals & class
                total = sum(vr.criteria_points.values()) if vr.criteria_points else 0
                vr.total_points = int(total)
                vr.automated_class = self._classify_by_points(vr.total_points)

# ---------------------------
# VEP annotation helper
# ---------------------------
def run_local_vep(in_vcf: str, out_vcf: str, vep_cmd: str, vep_cache: Optional[str], fasta: Optional[str], extra_args: List[str], threads: int = 4) -> str:
    """
    Run VEP locally to annotate VCF. The user must provide vep executable path and optionally cache and fasta.
    extra_args: list of extra arguments (e.g. ["--plugin", "SpliceAI,snv,0.2", "--plugin", "REVEL"])
    """
    cmd = [vep_cmd, "--input_file", in_vcf, "--output_file", out_vcf, "--vcf", "--compress_output", "bgzip", "--offline", "--cache", "--everything", "--fork", str(threads)]
    if vep_cache:
        cmd += ["--dir_cache", vep_cache]
    if fasta:
        cmd += ["--fasta", fasta]
    if extra_args:
        cmd += extra_args
    logger.info("Running VEP. Ensure VEP cache and plugins are configured correctly.")
    run_cmd(cmd)
    gz = out_vcf if out_vcf.endswith(".gz") else out_vcf + ".gz"
    try:
        run_cmd(["tabix", "-p", "vcf", gz], check=False)
    except Exception:
        logger.debug("tabix missing or indexing failed.")
    return gz

# ---------------------------
# Build VariantRecord from cyvcf2 record
# ---------------------------
def build_variant_record_from_cyvcf2(rec, csq_header: List[str], acmg_genes: set, dad_vf: Optional[pysam.VariantFile], mom_vf: Optional[pysam.VariantFile]) -> Optional[VariantRecord]:
    """
    Parse cyvcf2 record and select MANE annotation for ACMG SF genes.
    Extract computational annotations (REVEL, SpliceAI, CADD), gnomAD AF (from CSQ), and genotype fields (AD/DP).
    Also fetch parental genotypes from provided pysam VariantFile objects if present.
    """
    if 'CSQ' not in rec.INFO:
        return None
    csq_raw = rec.INFO.get('CSQ')
    csq_lines = csq_raw if isinstance(csq_raw, (list,tuple)) else str(csq_raw).split(',')
    best_ann = None
    for line in csq_lines:
        ann = parse_vep_csq_line(line, csq_header)
        gene = ann.get('SYMBOL') or ann.get('Gene') or ""
        if not gene:
            continue
        if gene not in acmg_genes:
            continue
        # prefer MANE
        if ann.get('MANE_SELECT','') or ann.get('MANE','') or ann.get('MANE_PLUS_CLINICAL',''):
            best_ann = ann
            break
        if best_ann is None:
            best_ann = ann
    if not best_ann:
        return None

    gene = best_ann.get('SYMBOL') or best_ann.get('Gene') or ""
    transcript = best_ann.get('Feature') or best_ann.get('Feature_id') or ""
    hgvsc = best_ann.get('HGVSc') or ""
    hgvsp = best_ann.get('HGVSp') or ""
    consequence = best_ann.get('Consequence') or best_ann.get('Consequence_terms') or ""
    exon = best_ann.get('EXON') or ""
    is_last_exon = False
    if exon and '/' in exon:
        try:
            a,b = exon.split('/')
            is_last_exon = (a == b)
        except Exception:
            is_last_exon = False

    def get_num(keys):
        for k in keys:
            v = best_ann.get(k)
            if v not in (None, "", "."):
                try:
                    return float(v)
                except Exception:
                    if isinstance(v, str) and "|" in v:
                        parts = v.split("|")
                        nums = []
                        for p in parts[1:]:
                            try:
                                nums.append(float(p))
                            except Exception:
                                pass
                        if nums:
                            return max(nums)
        return None

    revel = get_num(['REVEL_score','REVEL'])
    mutpred2 = get_num(['MutPred2_score','MutPred2'])
    spliceai = get_num(['SpliceAI_pred_DS_max','SpliceAI','spliceai'])
    cadd = get_num(['CADD_PHRED','CADD'])
    gnomad_af = get_num(['gnomAD_AF','gnomad_AF','AF'])

    sample_name = rec.samples[0].name if rec.samples else ""
    gt = rec.genotypes[0][:2] if rec.genotypes else (None, None)
    prob_gt = "/".join("." if x is None else str(x) for x in gt)

    ad = None
    dp = None
    ab = None
    try:
        adf = rec.format('AD')
        if adf is not None and len(adf) > 0:
            advals = adf[0]
            if isinstance(advals, (list,tuple)) and len(advals) >= 2:
                refad = int(advals[0]); altad = int(advals[1])
                ad = (refad, altad)
                dp = refad + altad
                if dp > 0:
                    ab = altad / float(dp)
    except Exception:
        try:
            dpvals = rec.format('DP')
            if dpvals is not None and len(dpvals) > 0:
                dp = int(dpvals[0])
        except Exception:
            dp = None

    father_gt = "./."
    mother_gt = "./."
    if dad_vf:
        try:
            for dr in dad_vf.fetch(rec.CHROM, rec.POS - 1, rec.POS):
                if getattr(dr, 'ref', None) == rec.REF and (dr.alts and rec.ALT[0] in dr.alts):
                    p = next(iter(dr.samples))
                    g = dr.samples[p].get('GT')
                    father_gt = "/".join("." if x is None else str(x) for x in g)
                    break
        except Exception:
            pass
    if mom_vf:
        try:
            for mr in mom_vf.fetch(rec.CHROM, rec.POS - 1, rec.POS):
                if getattr(mr, 'ref', None) == rec.REF and (mr.alts and rec.ALT[0] in mr.alts):
                    p = next(iter(mr.samples))
                    g = mr.samples[p].get('GT')
                    mother_gt = "/".join("." if x is None else str(x) for x in g)
                    break
        except Exception:
            pass

    vr = VariantRecord(
        chrom=rec.CHROM,
        pos=int(rec.POS),
        ref=rec.REF,
        alt=str(rec.ALT[0]) if rec.ALT else "",
        sample=sample_name,
        gene=gene,
        transcript=transcript,
        consequence=consequence,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        exon=exon,
        is_last_exon=is_last_exon,
        gnomad_af=float(gnomad_af) if gnomad_af is not None else None,
        revel=float(revel) if revel is not None else None,
        mutpred2=float(mutpred2) if mutpred2 is not None else None,
        spliceai=float(spliceai) if spliceai is not None else None,
        cadd=float(cadd) if cadd is not None else None,
        nmd=best_ann.get('NMD','') or "",
        proband_gt=prob_gt,
        father_gt=father_gt,
        mother_gt=mother_gt,
        proband_dp=dp,
        proband_ad=ad,
        proband_ab=ab
    )
    return vr

# ---------------------------
# Pipeline orchestration
# ---------------------------
def process_vcf(proband_vcf: str,
                father_vcf: Optional[str],
                mother_vcf: Optional[str],
                outdir: str,
                vep_cmd: str,
                vep_cache: Optional[str],
                fasta: Optional[str],
                vep_extra: List[str],
                acmg_tsv: str,
                db_paths: Dict[str,str],
                ttn_meta_csv: Optional[str],
                aggressive: bool = True):
    """
    Full run:
      - ensure VEP annotation
      - parse CSQ and collect candidate variants restricted to reportable genes
      - per-variant scoring (PVS1/PS1/PP3/PS2)
      - PM3 batch pass
      - apply DB-driven overrides (if aggressive)
      - write outputs CSVs and run_info
    """
    os.makedirs(outdir, exist_ok=True)

    # load gene-specific rules
    gene_rules = load_acmg_tsv(acmg_tsv)
    acmg_genes = set(g for g, r in gene_rules.items() if r.get("reportable", False))
    logger.info("Reportable ACMG SF genes loaded: %d", len(acmg_genes))

    # ensure VEP annotation present
    annotated_vcf = proband_vcf
    if not vcf_has_csq(proband_vcf):
        logger.info("Input VCF missing VEP CSQ. Running local VEP (requires VEP and plugins configured).")
        if not vep_available(vep_cmd):
            raise RuntimeError("VEP not found; provide correct --vep path or pre-annotate VCF with CSQ.")
        out_v = os.path.join(outdir, os.path.basename(proband_vcf).rstrip(".vcf") + ".vep.vcf")
        annotated_vcf = run_local_vep(proband_vcf, out_v, vep_cmd, vep_cache, fasta, vep_extra)
        logger.info("VEP produced: %s", annotated_vcf)
    else:
        logger.info("Input VCF already contains CSQ from VEP.")

    # load DBs and TTN meta
    dbm = DatabaseManager(db_paths)
    ttn_trees = {}
    if ttn_meta_csv:
        if not IntervalTree:
            logger.error("intervaltree python package is required for TTN meta handling. Install it and retry.")
            raise RuntimeError("Missing intervaltree dependency")
        if os.path.exists(ttn_meta_csv):
            ttn_df = pd.read_csv(ttn_meta_csv)
            for _, r in ttn_df.iterrows():
                chrom = r.get('chr') or r.get('chrom') or r.get('chromosome')
                if not chrom:
                    continue
                if not str(chrom).startswith("chr"):
                    chrom = "chr" + str(chrom)
                if chrom not in ttn_trees:
                    ttn_trees[chrom] = IntervalTree()
                start = int(r.get('start') or 0)
                end = int(r.get('end') or 0)
                ttn_trees[chrom].addi(start, end, r.to_dict())
            logger.info("Loaded TTN meta-exon data for %d chromosomes", len(ttn_trees))
        else:
            logger.warning("TTN meta-exon CSV not found at %s", ttn_meta_csv)

    engine = ACMGEngine(dbm, gene_rules, ttn_meta=ttn_trees, spliceai_thr=THRESH["SPLICEAI"])

    csq_header = parse_csq_header_from_vcf(annotated_vcf)
    if not csq_header:
        logger.warning("CSQ header not parsed; VEP output may not contain required fields.")

    dad_vf = pysam.VariantFile(father_vcf) if father_vcf and os.path.exists(father_vcf) else None
    mom_vf = pysam.VariantFile(mother_vcf) if mother_vcf and os.path.exists(mother_vcf) else None

    vcf_reader = VCF(annotated_vcf)
    candidates: List[VariantRecord] = []
    for rec in vcf_reader:
        vr = build_variant_record_from_cyvcf2(rec, csq_header, acmg_genes, dad_vf, mom_vf)
        if not vr:
            continue
        # annotate gnomAD from provided gnomAD VCF if available
        if (vr.gnomad_af is None or vr.gnomad_af == 0.0) and dbm.gnomad_vf:
            dbm.annotate_gnomad_for_variant(vr)
        candidates.append(vr)
    logger.info("Collected %d candidate variants in ACMG SF genes", len(candidates))

    # First pass scoring per variant
    for v in candidates:
        engine.score_variant(v)

    # PM3 (batch across sample/gene)
    engine.apply_pm3_batch(candidates)

    # finalize: recompute totals and apply DB-driven overrides if aggressive
    for v in candidates:
        total = 0
        for _, pts in v.criteria_points.items():
            try:
                total += int(pts)
            except Exception:
                try:
                    total += float(pts)
                except Exception:
                    pass
        v.total_points = int(total)
        v.automated_class = engine._classify_by_points(v.total_points)
        if aggressive and v.db_hits:
            priority = {"ClinGen":4, "ClinVar":3, "HGMD_Pro":3, "InternalDB":2}
            selected = None
            best_score = 0
            for cls, src in v.db_hits:
                score = priority.get(src, 0)
                if score > best_score:
                    best_score = score
                    selected = (cls, src)
            if selected:
                cls, src = selected
                cls_str = str(cls).lower()
                if "pathogen" in cls_str:
                    v.automated_class = "Pathogenic"
                    v.manual_reasons.append(f"Automated DB-driven override to Pathogenic ({src})")
                elif "likely" in cls_str:
                    v.automated_class = "Likely pathogenic"
                    v.manual_reasons.append(f"Automated DB-driven override to Likely Pathogenic ({src})")

# ---------------------------
# outputs
# ---------------------------
rows = []
manual_rows = []
auto_rows = []

for v in candidates:
    # recompute/ensure totals & class are up to date
    total = 0
    for _, pts in v.criteria_points.items():
        try:
            total += int(pts)
        except Exception:
            try:
                total += float(pts)
            except Exception:
                pass
    v.total_points = int(total)
    v.automated_class = engine._classify_by_points(v.total_points)

    # If "aggressive" DB-driven override is desired it was already applied above;
    # here we rely on v.automated_class and v.manual_review values.

    row = {
        "sample": v.sample,
        "chrom": v.chrom,
        "pos": v.pos,
        "ref": v.ref,
        "alt": v.alt,
        "gene": v.gene,
        "transcript": v.transcript,
        "consequence": v.consequence,
        "hgvsc": v.hgvsc,
        "hgvsp": v.hgvsp,
        "exon": v.exon,
        "is_last_exon": v.is_last_exon,
        "gnomad_af": v.gnomad_af,
        "revel": v.revel,
        "mutpred2": v.mutpred2,
        "spliceai": v.spliceai,
        "cadd": v.cadd,
        "nmd": v.nmd,
        "proband_gt": v.proband_gt,
        "father_gt": v.father_gt,
        "mother_gt": v.mother_gt,
        "proband_dp": v.proband_dp,
        "proband_ad": v.proband_ad,
        "proband_ab": v.proband_ab,
        "db_hits": ";".join([f"{c}|{s}" for (c,s) in v.db_hits]) if v.db_hits else "",
        "criteria_assigned": ";".join(v.criteria_assigned),
        "criteria_points": json.dumps(v.criteria_points, ensure_ascii=False),
        "total_points": v.total_points,
        "automated_class": v.automated_class,
        "manual_review": v.manual_review,
        "manual_reasons": ";".join(v.manual_reasons)
    }
    rows.append(row)

    # Only Pathogenic / Likely pathogenic are considered secondary findings.
    # Put them into manual or automated lists depending on manual_review flag.
    if v.automated_class in ("Pathogenic", "Likely pathogenic"):
        if v.manual_review:
            manual_rows.append(row)
        else:
            auto_rows.append(row)

# write outputs
out_all = os.path.join(outdir, "all_candidates.csv")
out_manual = os.path.join(outdir, "manual_review_list.csv")
out_auto = os.path.join(outdir, "auto_conclusions.csv")

pd.DataFrame(rows).to_csv(out_all, index=False)
pd.DataFrame(manual_rows).to_csv(out_manual, index=False)
pd.DataFrame(auto_rows).to_csv(out_auto, index=False)

# run metadata 
run_info = {
    "proband_vcf": proband_vcf,
    "annotated_vcf": annotated_vcf,
    "vep_cmd": vep_cmd,
    "vep_cache": vep_cache,
    "fasta": fasta,
    "acmg_tsv": acmg_tsv,
    "db_paths": db_paths,
    "num_candidates_total": len(rows),
    "num_manual_plp": len(manual_rows),
    "num_auto_plp": len(auto_rows),
    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
}
with open(os.path.join(outdir, "run_info.json"), "w", encoding="utf-8") as fh:
    json.dump(run_info, fh, indent=2)

logger.info("Wrote outputs: %s, %s, %s", out_all, out_manual, out_auto)
logger.info("Run summary: candidates_total=%d automated_PLP=%d manual_PLP=%d", len(rows), len(auto_rows), len(manual_rows))

# ---------------------------
# CLI
# ---------------------------
def main():
    parser = argparse.ArgumentParser(description="ACMG SF classifier (VEP-based, gene-specific rules, ClinVar protein-index caching)")
    parser.add_argument("--proband", required=True, help="Proband VCF (bgz recommended)")
    parser.add_argument("--father", help="Father VCF (optional)")
    parser.add_argument("--mother", help="Mother VCF (optional)")
    parser.add_argument("--outdir", default="acmg_sf_results", help="Output directory")
    parser.add_argument("--vep", default="vep", help="VEP executable path")
    parser.add_argument("--vep-cache", help="VEP cache directory")
    parser.add_argument("--fasta", help="Reference FASTA for VEP")
    parser.add_argument("--vep-extra", nargs="*", default=[], help="Extra tokenized args to pass to VEP (plugins etc.)")
    parser.add_argument("--acmg-tsv", required=True, help="ACMG SF TSV with gene-specific rules (your file)")
    parser.add_argument("--ttn-meta", help="TTN meta-exon CSV (columns: chr,start,end,psi_meta,lof_allowed)")
    parser.add_argument("--db-paths-json", help="JSON file to override default DB paths")
    parser.add_argument("--aggressive", action="store_true", help="Favor DB-driven automated overrides to reduce manual review")
    args = parser.parse_args()

    db_paths = dict(DB_PATHS_DEFAULT)
    if args.db_paths_json:
        if os.path.exists(args.db_paths_json):
            with open(args.db_paths_json) as fh:
                db_paths.update(json.load(fh))
        else:
            logger.warning("--db-paths-json not found; using defaults")

    process_vcf(
        proband_vcf=args.proband,
        father_vcf=args.father,
        mother_vcf=args.mother,
        outdir=args.outdir,
        vep_cmd=args.vep,
        vep_cache=args.vep_cache,
        fasta=args.fasta,
        vep_extra=args.vep_extra,
        acmg_tsv=args.acmg_tsv,
        db_paths=db_paths,
        ttn_meta_csv=args.ttn_meta,
        aggressive=args.aggressive
    )

if __name__ == "__main__":
    main()
