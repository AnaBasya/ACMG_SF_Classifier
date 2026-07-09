#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACMG Secondary Findings Classifier 

Purpose:
  Automates ACMG Secondary Findings (SF v3.x) classification for single samples or cohorts.
  Uses annotation, aggregates ClinVar/HGMD/internal and gnomAD evidence, 
  applies gene-specific ACMG SF rules (including transcript-aware PVS1/NMD logic),
  evaluates PS1/PM5, PM3 batching (phasing-aware), delins detection, and triages results
  into automatic reports and manual-review lists.

Inputs:
  - Proband VCF (optionally father/mother VCFs), ACMG SF table (CSV/TSV),
    optional: dbNSFP, ClinVar VCF, HGMD CSV, internal DB CSV, gnomAD v2/v3/v4, exons/NMD/transcript maps.

Outputs:
  - all_candidates.csv, auto_conclusions.csv, manual_review_list.csv, run_info.json

Usage hint:
  Run `python3 acmg_sf_classifier.py --help` for full CLI options and recommended workflows.

Author: 
  Clinical Bioinformatics Team of the Research Centre for Medical Genetics (Code by Anna Basova)
"""

from __future__ import annotations
from urllib.parse import unquote
import os
import sys
import json
import argparse
import logging
import subprocess
import shutil
import shlex
import time
import gzip
import re
import sqlite3
import glob
import concurrent.futures
import csv
import ast 
import traceback
import difflib
from dataclasses import dataclass, field
from typing import Tuple, List, Dict, Any, Optional, Set
from collections import defaultdict
# from cooccurrence_helper import query_pair_gnomad, fetch_cooccurrence_for_pairs (cooccurrence analysis handled by external script)

# optional libs
try:
    import pandas as pd
    HAVE_PANDAS = True
except Exception:
    pd = None
    HAVE_PANDAS = False

try:
    import pysam
    HAVE_PYSAM = True
except Exception:
    pysam = None
    HAVE_PYSAM = False

try:
    from cyvcf2 import VCF
    HAVE_CYVCF2 = True
except Exception:
    VCF = None
    HAVE_CYVCF2 = False

try:
    from intervaltree import IntervalTree
    HAVE_INTERVALTREE = True
except Exception:
    IntervalTree = None
    HAVE_INTERVALTREE = False

# ---------------------------
# Logging
# ---------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("acmg_sf_classifier")

# ---------------------------
# Defaults & thresholds
# ---------------------------
DB_PATHS_DEFAULT = {
    "CLINVAR_VCF": "databases/clinvar/clinvar_with_variant.mane_filtered.vcf.gz",
    "CLINVAR_EXPANDED_CSV": "databases/clinvar/clinvar_conditions_summary.filtered.with_hgvsp.csv",
    "CLINVAR_FILTERED_CSV": "databases/clinvar/clinvar_conditions_summary.filtered.with_hgvsp.csv",
    

}

THRESH = {
    # Computational Evidence (PP3)
    "REVEL_SUPPORTING": 0.64,
    "REVEL_MODERATE": 0.93,
    "ALPHA_MISSENSE_SUPPORTING": 0.56,
    # "CADD_SUPPORTING": 30,  
    "SQUIRLS_SUPPORTING": 0.70,             

    # Splicing Predictors
    "SPLICEAI_PVS1": 0.70,
    "SPLICEAI_MODERATE": 0.20,
    "SPLICEAI_PS1_THRESHOLD": 0.10,

    # Population Frequency
    "BA1_AF": 0.05,
    "BS1_AF": 0.01,

    # PM2 Thresholds
    "PM2_AD_AF": 0.0001,
    "PM2_AR_AF": 0.005,
    "PM2_XLD_AF": 0.0001,
    "PM2_XLR_AF": 0.003,

    # BA1 Exception List
    "BA1_EXCEPTIONS": {
        "HFE": ["p.Cys282Tyr"],
        "BTD": ["p.Asp444His"]
    },

    "DSCORE_MANUAL_THRESHOLD": 0.7,      # If Dscore >= 0.7 and Impact MODERATE -> manual for VUS 1-5
    "DSCORE_PERFECT": 1.0,               # If Dscore == 1.0 -> manual for VUS 1-5 (any Impact)
    "VUS_POINTS_MIN_FOR_MANUAL": 1,      # Minimum points for VUS to trigger manual review
    "VUS_POINTS_MAX_FOR_MANUAL": 5       # Maximum points for VUS to trigger manual review
}

# -----------------------------------------------------------
# Points map and classification thresholds
# -----------------------------------------------------------
POINTS_MAP = {
    "PVS1_VeryStrong": 8,
    "PVS1_Strong": 4,
    "PVS1_Moderate": 2,
    "PVS1_Supporting": 1,

    "PS1_Strong": 4,
    "PS1_Moderate": 2,
    "PS1_Supporting": 1,

    "PS2_VeryStrong": 8,
    "PS2_Moderate": 2,

    "PM1": 2,
    "PM2_Moderate": 2,

    "PM3_VeryStrong": 8,
    "PM3_Strong": 4,
    "PM3_Moderate": 2,
    "PM3_Supporting": 1,

    "PM4_Strong": 4,
    "PM4_Moderate": 2,
    "PM4_Supporting": 1,

    "PM5_Moderate": 2,
    "PM5_Supporting": 1,

    "PP1_Strong": 4,
    "PP1_Moderate": 2,
    "PP1_Supporting": 1,

    "PP3_Moderate": 2,
    "PP3_Supporting": 1,

    "BA1": -8,
    "BS1": -4,
    "BS2": -4,
    "BP4": -1,
    "BP6": -1,
    "BP7": -1,

    # extended PP5/BP5 strength levels (PP5ext / BP5ext)
    "PP5ext_VeryStrong": 8,
    "PP5ext_Strong": 4,
    "PP5ext_Moderate": 2,
    "PP5ext_Supporting": 1,

    "BP5ext_Supporting": -1,
    "BP5ext_Moderate": -2,
    "BP5ext_Strong": -4,
}

CLASSIFICATION_THRESHOLDS = [
    ("Pathogenic", 10, None),
    ("Likely pathogenic", 6, 9),
    ("VUS", 0, 5),
    ("Likely benign", -5, -1),
    ("Benign", None, -6)
]

MIN_PATHOGENIC_CRITERIA = 2
MIN_BENIGN_CRITERIA = 2

LOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "frameshift", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
CANONICAL_SPLICE = {"splice_acceptor_variant", "splice_donor_variant"}
STRICT_QC_ENABLED = True
STRICT_MIN_DP = 10
STRICT_MIN_AF = 0.20

# ---------------------------
# Data model
# ---------------------------
@dataclass
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    sample: str

    query_variant: str = ""
    gene: str = ""
    gene_raw: str = ""
    transcript: str = ""
    consequence: str = ""
    hgvsc: str = ""
    hgvsp: str = ""
    tid: Optional[str] = ""
    exon: str = ""
    is_last_exon: bool = False
    nmd: Optional[str] = ""

    sample_sex: str = "Unknown"

    gnomad_af: Optional[float] = None
    gnomad_details: Dict[str, Any] = field(default_factory=dict)

    clinvar_sig: str = ""
    clinvar_trait: str = ""
    clinvar_review: str = ""
    clinvar_stars: str = "0"
    clinvar_phenotype_info: str = ""            
    clinvar_consensus_star: Optional[int] = None  
    clinvar_phenotype_match: str = "No"   

    internal_db_sig: str = ""
    internal_db_details: str = ""

    hgmd_id: str = ""
    hgmd_class: str = ""
    hgmd_phen: str = ""
    hgmd_rankscore: str = "Not provided"
    hgmd_support_positive: int = 0
    hgmd_support_total: int = 0

    pathogenic_criteria_count: int = 0
    benign_criteria_count: int = 0
    hgmd_publication_count: int = 0

    dscore: Optional[float] = None

    gnomad_hom: Optional[int] = None
    revel: Optional[float] = None
    alpha_missense: Optional[float] = None
    alpha_missense_category: Optional[str] = None  # 'P', 'B', 'A', etc.
    spliceai: Optional[float] = None
    #cadd: Optional[float] = None
    nmdesc_predictor: Optional[str] = ""
    squirls_score: Optional[float] = None
    squirls_category: Optional[str] = None

    proband_gt: str = "./."
    father_gt: str = "./."
    mother_gt: str = "./."
    proband_dp: Optional[int] = None
    proband_ad: Optional[Tuple[int,int]] = None
    proband_af: Optional[float] = None

    db_hits: List[Tuple[str,str]] = field(default_factory=list)
    criteria_assigned: List[str] = field(default_factory=list)
    criteria_points: Dict[str,int] = field(default_factory=dict)
    criteria_explanations: Dict[str,str] = field(default_factory=dict)
    total_points: int = 0
    automated_class: str = "VUS"
    manual_review: bool = False
    manual_reasons: List[str] = field(default_factory=list)
    _raw_ann: Dict[str,Any] = field(default_factory=dict)

    has_strong_pathogenic_evidence: bool = False
    has_strong_benign_evidence: bool = False
    has_weak_pathogenic_evidence: bool = False

    is_hq_conflict_db: bool = False
    conflict_details: str = ""

    is_hq_pathogenic_db: bool = False
    is_hq_benign_db: bool = False

    _disease_match: bool = False
    _triage_group: str = ""

    ps1_matches: List[Dict] = field(default_factory=list)
    pm5_matches: List[Dict] = field(default_factory=list)
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    aa_pos: Optional[int] = None
    phase: str = "unknown"
    is_delins_part: bool = False

    filtered_reasons: List[str] = field(default_factory=list)
    auto_report_reasons: List[str] = field(default_factory=list)

# ---------------------------
# Utilities
# ---------------------------
def run_cmd(cmd: List[str], check=True, timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    logger.debug("Running: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if check and proc.returncode != 0:
        logger.error("Command failed (%d): %s\nstdout:\n%s\nstderr:\n%s", proc.returncode, " ".join(cmd), proc.stdout, proc.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return proc

_RE_LEADING_TRAILING = re.compile(r'^[\s"\'`]+|[\s"\'`]+$')
_RE_NON_ALNUM = re.compile(r'[^A-Za-z0-9_]')

def normalize_clinvar_label(sig: str) -> str:
    if not sig:
        return "Not provided"
    
    s_lower = sig.lower().replace("_", " ")
    
    if "conflict" in s_lower:
        return "Conflicting classifications of pathogenicity"
    
    if "likely pathogenic" in s_lower:
        return "Likely pathogenic"
    
    if "pathogenic" in s_lower and "likely" not in s_lower:
        return "Pathogenic"
    
    if "likely benign" in s_lower:
        return "Likely benign"
    
    if "benign" in s_lower and "pathogen" not in s_lower:
        return "Benign"
    
    if "uncertain" in s_lower or "vus" in s_lower:
        return "Uncertain significance"
    
    return sig

def calculate_consensus_star(submissions: List[Dict]) -> int:
    """
    Calculate consensus star rating following clinvar_merge.py logic:
    - If any submission has 3-4 stars -> return max star
    - If >=2 one-star submissions agree on vector -> promote to 2 stars
    - Otherwise return max star among all matching submissions
    """
    if not submissions:
        return 0
    
    stars_list = [int(sub.get("stars", 0) or 0) for sub in submissions]
    max_star = max(stars_list) if stars_list else 0
    
    # 3-4 stars have priority
    if max_star >= 3:
        return max_star
    
    # Count one-star submissions that are NOT conflicts
    one_star_subs = []
    for sub in submissions:
        if int(sub.get("stars", 0) or 0) == 1:
            sig = str(sub.get("clnsig", "")).lower()
            # Skip conflict submissions for promotion (they don't count toward agreement)
            if "conflict" not in sig:
                one_star_subs.append(sub)
    
    if len(one_star_subs) >= 2:
        # Check vector agreement among one-star submissions
        sigs_one_star = [str(sub.get("clnsig", "")).lower() for sub in one_star_subs]
        
        has_path = any('pathogenic' in s and 'benign' not in s and 'conflict' not in s for s in sigs_one_star)
        has_benign = any('benign' in s and 'pathogen' not in s for s in sigs_one_star)
        has_vus = any('uncertain' in s or 'vus' in s for s in sigs_one_star)
        
        categories = sum([has_path, has_benign, has_vus])
        
        if categories <= 1:
            # No conflict -> promote to 2 stars
            return 2
    
    return max_star

def _uniq_preserve(items):
    """
    Remove exact duplicates and empty strings while preserving order.
    
    Args:
        items: List of items (strings)
    
    Returns:
        List with duplicates removed, order preserved
    """
    out = []
    seen = set()
    for it in items:
        if not it:
            continue
        s = str(it).strip()
        if not s or s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out

def aggregate_clinvar_signatures(vr):
    try:
        _clin_sigs_all = []

        # 1) from expanded submissions stored on variant record (vr.clinvar_submissions)
        for s in getattr(vr, "clinvar_submissions", []) or []:
            try:
                cs = ""
                if isinstance(s, dict):
                    cs = str(s.get("clnsig", "") or "").strip()
                else:
                    cs = str(s or "").strip()
                if not cs:
                    continue
                parts = re.split(r'[;,\|/]+', cs)
                for p in parts:
                    p2 = p.strip()
                    if p2:
                        _clin_sigs_all.append(p2)
            except Exception:
                continue

        # 2) from any existing vr.clinvar_sig (e.g. earlier VCF-based lookup)
        existing_sig = getattr(vr, "clinvar_sig", "") or ""
        if existing_sig and str(existing_sig).strip() and str(existing_sig).strip().lower() not in ("not provided", ".", "null", "none"):
            try:
                parts = re.split(r'[;,\|/]+', str(existing_sig))
                for p in parts:
                    p2 = p.strip()
                    if p2:
                        _clin_sigs_all.append(p2)
            except Exception:
                pass

        # Deduplicate preserving order (case-insensitive uniqueness, keep original case)
        uniq = []
        seen = set()
        for item in _clin_sigs_all:
            key = item.strip().lower()
            if not key:
                continue
            if key in seen:
                continue
            seen.add(key)
            uniq.append(item.strip())

        # Final assignment
        if uniq:
            vr.clinvar_sig = ";".join(uniq)
        else:
            if existing_sig and existing_sig.strip() and existing_sig.strip().lower() not in ("not provided", ".", "null", "none"):
                vr.clinvar_sig = existing_sig
            else:
                vr.clinvar_sig = "Not provided"
    except Exception:
        try:
            if not getattr(vr, "clinvar_sig", None):
                vr.clinvar_sig = "Not provided"
        except Exception:
            pass

def get_sample_names_from_vcf_header(vcf_path: str) -> List[str]:
    """
    Extract sample names from VCF header.
    """
    samples = []
    try:
        opener = gzip.open if (str(vcf_path).endswith((".gz", ".bgz")) or is_gzipped(vcf_path)) else open
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith("#CHROM"):
                    parts = re.split(r'\t+', ln.rstrip("\n"))
                    if len(parts) > 9:
                        samples = parts[9:]
                    break
    except Exception as e:
        logger.warning(f"Could not extract sample names from {vcf_path}: {e}")
    return samples

def normalize_gene_name(g: str) -> str:
    if not g:
        return ""
    s = str(g)
    s = _RE_LEADING_TRAILING.sub("", s)
    if (s.startswith("(") and s.endswith(")")) or (s.startswith("[") and s.endswith("]")):
        s = s[1:-1].strip()
    s = s.replace('"','').replace("'", "").replace("`","")
    s = re.sub(r'\s+', '_', s)
    s = _RE_NON_ALNUM.sub('_', s)
    s = re.sub(r'_+', '_', s)
    s = s.strip('_').upper()
    return s


def _normalize_spelling_token(tok: str) -> str:
    """Normalize token for spelling variations (UK/US, ligatures, etc.)"""
    if not tok:
        return ""
    s = str(tok).lower()
    
    # Extended replacements for medical terms
    replacements = {
        # UK/US variations
        'ae': 'e', 'oe': 'e', 'leuka': 'leuka', 'leuk': 'leuk',
        'haem': 'hem', 'haemo': 'hemo', 'tumour': 'tumor', 'tumours': 'tumors',
        'colour': 'color', 'favour': 'favor', 'metre': 'meter',
        'centre': 'center', 'fibre': 'fiber', 'anaemia': 'anemia',
        'anaesthesia': 'anesthesia', 'diarrhoea': 'diarrhea',
        'oedema': 'edema', 'oesophagus': 'esophagus',
        'foetus': 'fetus', 'manoeuvre': 'maneuver',
        # Medical term variations
        'adrenoleukodystrophy': 'adrenoleukodystrophy',
        'adrenoleucodystrophy': 'adrenoleukodystrophy',  # Normalize both to same
        'leucodystrophy': 'leukodystrophy',
        'leuco': 'leuko',
        # Common ligatures
        'æ': 'ae', 'œ': 'oe',
    }
    
    for old, new in replacements.items():
        s = s.replace(old, new)
    
    # Remove punctuation and extra spaces
    s = re.sub(r'[^a-z0-9\s]', '', s)
    s = re.sub(r'\s+', ' ', s).strip()
    
    return s

def _token_similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    na = _normalize_spelling_token(a)
    nb = _normalize_spelling_token(b)
    if not na or not nb:
        return 0.0
    # exact normalized equality
    if na == nb:
        return 1.0
    # substring containment (helpful for long disease names)
    if len(na) >= 5 and na in nb:
        return 1.0
    if len(nb) >= 5 and nb in na:
        return 1.0
    # difflib fallback
    return difflib.SequenceMatcher(None, na, nb).ratio()

def check_disease_match(acmg_disease_string: str, variant_trait: str, 
                       extra_recs_list: List[Dict] = None,
                       gene_rules: Dict[str, Any] = None) -> Tuple[bool, str, str]:
    """
    Improved disease matching:
     - preserves prior behavior (direct substring / token intersection)
     - excludes generic disease words (syndrome, disorder, disease, phenotype...) from token-matching
     - adds token-level fuzzy matching (with AE/OE normalization) with a conservative similarity threshold
       so we accept small orthographic differences (e.g., 'cholesterolemia' vs 'cholesterolaemia')
     - extra_recs_list (whitelist/exclusion) is kept as before
    Returns (matched: bool, match_level: str, reason: str)
    """
    vt_raw = str(variant_trait or "").strip()
    acmg_raw = str(acmg_disease_string or "").strip()
    # if no ACMG disease -> cannot match
    if not acmg_raw:
        return False, "no_acmg_disease", "No ACMG disease defined"

    # Cleaning & normalization helper (preserve previous normalization but keep tokens)
    def norm(s: str) -> str:
        s2 = str(s).lower()
        s2 = s2.replace("_", " ").replace("|", " ").replace("/", " ")
        # collapse punctuation to spaces
        s2 = re.sub(r'[^\w\s]', ' ', s2)
        s2 = re.sub(r'\s+', ' ', s2).strip()
        s2 = re.sub(r'^[^a-z0-9]+|[^a-z0-9]+$', '', s2)
        return s2

    acmg_norm = norm(acmg_raw)
    vt_norm = norm(vt_raw)

    # If DB trait empty -> mark as unknown DB phenotype (do not treat as explicit absence)
    if not vt_norm:
        return False, "unknown_db_phenotype", "No phenotype text in DB (unknown)"

    # Whitelist/exclusion matching (unchanged logic)
    if extra_recs_list:
        for rec in extra_recs_list:
            target = str(rec.get("disease","") or "").lower().replace("_"," ").strip()
            if not target:
                continue
            target_norm = norm(target)
            if target_norm and (target_norm in vt_norm or target_norm in acmg_norm):
                if rec.get("inclusion", True):
                    return True, "whitelist_match", f"Matched recommendation: {target}"
                else:
                    return False, "whitelist_exclude", f"Excluded by recommendation: {target}"

    # 1) Direct substring / containment (fast, exact)
    if acmg_norm and (acmg_norm in vt_norm or vt_norm in acmg_norm):
        return True, "direct_match", "Direct substring match"

    # 2) Tokenize and remove generic disease words (stopwords)
    # We exclude short tokens (<=3 chars) and generic disease descriptors that cause false positives (e.g., 'syndrome')
    GENERIC_TERMS = {
        "syndrome", "syndromes", "disease", "disorder", "disorders", 
        "phenotype", "phenotypes", "trait", "traits", "susceptibility", "tumour", "tumor", "tumours"
    }

    def tokens(s: str) -> List[str]:
        toks = [t for t in re.split(r'\s+', s) if t and len(t) > 3]
        # filter out pure generic disease words
        toks = [t for t in toks if t not in GENERIC_TERMS]
        return toks

    acmg_tokens = set(tokens(acmg_norm))
    vt_tokens = set(tokens(vt_norm))

    # exact token overlap (after removing generic words)
    common = acmg_tokens.intersection(vt_tokens)
    if common:
        return True, "token_match", f"Shared terms: {', '.join(sorted(common))}"

    # 3) fuzzy token matching:
    # allow small orthographic differences; be conservative w/ threshold to avoid matching generic tokens
    # thresholds:
    #   - high similarity threshold for short words (len<=8): >=0.90
    #   - slightly lower for longer words: >=0.80
    # Also allow direct normalized-equality using our AE/OE normalization
    for at in acmg_tokens:
        for vt in vt_tokens:
            # normalized equality check (handles ae/oe -> e)
            if _normalize_spelling_token(at) == _normalize_spelling_token(vt) and len(_normalize_spelling_token(at)) >= 4:
                return True, "fuzzy_token_norm_equal", f"Normalized-equal tokens: {at} ~ {vt}"
            sim = _token_similarity(at, vt)
            # determine threshold by token length (use longer token length)
            L = max(len(_normalize_spelling_token(at)), len(_normalize_spelling_token(vt)))
            if L <= 8:
                thr = 0.88  
            else:
                thr = 0.77   
            if sim >= thr:
                return True, "fuzzy_token_match", f"Fuzzy match tokens: {at} ~ {vt} (sim={sim:.2f})"

    # 3.5) Advanced medical term matching with stemming
    # This handles variations like "adrenoleukodystrophy" vs "adrenoleucodystrophy"
    medical_stems = {
        'adrenoleukodystrophy': ['adrenoleukodystroph', 'adrenoleucodystroph', 'ald'],
        'leukodystrophy': ['leukodystroph', 'leucodystroph'],
        'cardiomyopathy': ['cardiomyopath', 'cardiomyopathy'],
        'channelopathy': ['channelopath'],
    }
    
    # Check if any ACMG token maps to a medical stem that appears in variant trait
    for acmg_tok in acmg_tokens:
        acmg_norm = _normalize_spelling_token(acmg_tok)
        for canonical, variants in medical_stems.items():
            if acmg_norm in variants or acmg_norm == canonical:
                for vt_tok in vt_tokens:
                    vt_norm = _normalize_spelling_token(vt_tok)
                    if vt_norm in variants or vt_norm == canonical:
                        return True, "medical_stem_match", f"Medical term match: {acmg_tok} ~ {vt_tok}"

    # 4) stem/token based similarity (fallback)
    def stem_token(tok: str) -> str:
        return re.sub(r'(itis|opathy|osis|ic|al|us|a|ia|ism)$', '', tok)
    acmg_stems = {stem_token(t) for t in acmg_tokens}
    vt_stems = {stem_token(t) for t in vt_tokens}
    if acmg_stems.intersection(vt_stems):
        return True, "stem_token_match", f"Stem overlap: {', '.join(sorted(acmg_stems.intersection(vt_stems)))}"

    return False, "mismatch", "No overlap found"

def normalize_chrom(chrom: str) -> str:
    c = str(chrom)
    if c.startswith("chr"):
        return c
    if re.match(r'^\d+$', c) or c in ("X","Y","MT","M"):
        return "chr" + c
    return c

def is_gzipped(path: str) -> bool:
    # Check file magic bytes to detect gzip content regardless of filename extension.
    # Return True if file starts with GZIP magic (0x1f 0x8b), else False.
    if not path or not os.path.exists(path):
        return False
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
            return magic == b"\x1f\x8b"
    except Exception:
        return False

# Parse semicolon-separated INFO string produced by generate.py
def parse_simple_info(info_str: str) -> Dict[str, str]:
    """
    Parse simple INFO string like "GENE=BRCA1;HGVSc=ENST...:c.123A>T;HGVSp=p.Val41Glu;REVEL=0.95;splice_ai_score=0.12;SEX=male"
    Returns dict with keys uppercased and values as strings.
    """
    out = {}
    if not info_str:
        return out
    parts = str(info_str).split(';')
    # DEBUG: print first few INFO strings
    if not hasattr(parse_simple_info, "_debug_count"):
        parse_simple_info._debug_count = 0
    if parse_simple_info._debug_count < 5:
        print(f"\nDEBUG parse_simple_info #{parse_simple_info._debug_count + 1}")
        print(f"  raw: {info_str[:200]}...")
        parse_simple_info._debug_count += 1
    for p in parts:
        if not p:
            continue
        if '=' in p:
            k, v = p.split('=', 1)
            k = k.strip()
            v = v.strip()
            if not k:
                continue
            out[k.upper()] = v
        else:
            # flag style INFO like 'SOMATIC' (unlikely here but keep)
            out[p.strip().upper()] = ""
    return out

# Parse SAMPLE field produced by generate.py (FORMAT GT:DP:GQ:AD:AF)
def parse_sample_field(sample_field: str, format_keys: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Parse a sample string like "0/1:35:99:20,15:0.300".
    If format_keys is None, assume GT:DP:GQ:AD:AF order used by generate.py.
    Returns dictionary with keys gt, dp, gq, ad (tuple), af (float)
    """
    out = {"gt": "./.", "dp": None, "gq": None, "ad": None, "af": None}
    if not sample_field:
        return out
    if format_keys is None:
        format_keys = ["GT", "DP", "GQ", "AD", "AF"]
    vals = sample_field.split(':')
    for idx, key in enumerate(format_keys):
        if idx >= len(vals):
            break
        v = vals[idx]
        k = key.upper()
        if k == "GT":
            out["gt"] = v or "./."
        elif k == "DP":
            try:
                out["dp"] = int(v) if v not in (".", "") else None
            except Exception:
                out["dp"] = None
        elif k == "GQ":
            try:
                out["gq"] = int(v) if v not in (".", "") else None
            except Exception:
                out["gq"] = None
        elif k == "AD":
            try:
                parts = v.split(',')
                ad = tuple(int(x) if x not in (".","") else 0 for x in parts)
                out["ad"] = ad
            except Exception:
                out["ad"] = None
        elif k == "AF":
            try:
                out["af"] = float(v) if v not in (".","") else None
            except Exception:
                out["af"] = None
        else:
            # store unknown keys as raw
            out[k] = v
    return out

# build VariantRecord from an extracted VCF-like line (one data line, not header)
def create_variant_from_extract_line(fields: List[str], sample_name: str = "proband") -> VariantRecord:
    """
    fields expected: CHROM, POS, ID (vid), REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
    This matches the SELECT ... produced by generate.py.
    
    IMPORTANT: sample_name is used as the sample identifier for this record.
    """
    # Basic columns
    chrom = fields[0]
    pos = int(fields[1]) if fields[1].isdigit() else int(float(fields[1]))
    vid = fields[2]
    ref = fields[3]
    alt = fields[4]
    info_raw = fields[7] if len(fields) > 7 else ""
    fmt = fields[8] if len(fields) > 8 else ""
    sample_field = fields[9] if len(fields) > 9 else ""

    # CRITICAL: Use the provided sample_name, not trying to infer from file
    vr = VariantRecord(
        chrom=normalize_chrom(chrom), 
        pos=pos, 
        ref=ref, 
        alt=alt, 
        sample=sample_name   # Explicitly set from parameter
    )

    try:
        vr.query_variant = f"{vr.chrom}-{vr.pos}-{vr.ref}-{vr.alt}"
    except Exception:
        vr.query_variant = f"{chrom}-{pos}-{ref}-{alt}"

    if not vr.sample or vr.sample.lower() in ("sample", "unknown", ""):
        vr.sample = "proband"
    # parse INFO -> vr._raw_ann
    ann = parse_simple_info(info_raw)

    vr._raw_ann = ann
    vr.tid = ann.get("TID") or ann.get("tid") or ann.get("TRANSCRIPT_ID") or vr.transcript or ""
    #  Extract sex from INFO field
    sample_sex = ann.get("SEX") or ann.get("Sex") or ann.get("sex") or ""
    if sample_sex:
        # Convert to expected format (Male/Female/Unknown)
        sex_lower = sample_sex.lower()
        if sex_lower in ("male", "м", "муж", "man", "m"):
            vr.sample_sex = "Male"
        elif sex_lower in ("female", "ж", "жен", "woman", "f"):
            vr.sample_sex = "Female"
        else:
            vr.sample_sex = "Unknown"
    else:
        vr.sample_sex = "Unknown"

    # Normalize Impact into vr._raw_ann and convenience attribute
    try:
        impact_val = (ann.get("IMPACT") or ann.get("Impact") or ann.get("impact") or "")
        impact_val = str(impact_val).strip() if impact_val is not None else ""
        vr._raw_ann["IMPACT"] = impact_val
        vr.impact = impact_val.upper() if impact_val else ""
    except Exception:
        vr.impact = ""

    # populate common fields from INFO (if present)
    # prefer keys exactly like generated: GENE, HGVSC, HGVSP, REVEL, ALPHAMISSENSE, SPLICE_AI_SCORE, TRANSCRIPT_SELECT, EXON
    vr.gene = ann.get("GENE") or ann.get("SYMBOL") or ann.get("GENE_SYMBOL") or ""
    vr.gene_raw = vr.gene
    vr.hgvsc = ann.get("HGVSC") or ann.get("CDNA") or ann.get("C_HGVS") or ""
    vr.hgvsp = ann.get("HGVSP") or ann.get("PROTEIN") or ann.get("P_HGVS") or ""
    # transcript: try MANE/transcript_select or parse from HGVSc if contains "ENST" or "NM_"
    tr = ann.get("TRANSCRIPT_SELECT") or ann.get("TRANSCRIPT") or ann.get("MANE_SELECT") or ""
    if tr:
        vr.transcript = tr
    else:
        # try to extract from HGVSc field like "ENST00000...:c.123"
        if vr.hgvsc and ":" in vr.hgvsc and vr.hgvsc.lower().startswith(("enst","nm_")):
            tx_part = vr.hgvsc.split(":", 1)[0]
            vr.transcript = tx_part

    # consequence (generate.py may not include VEP consequence). Try INFO key CONSEQUENCE or ANN/CSQ fallback
    vr.consequence = ann.get("CONSEQUENCE") or ann.get("CSQ") or ann.get("ANN") or ""

    # numeric scores
    try:
        if ann.get("REVEL") is not None:
            vr.revel = float(ann.get("REVEL"))
    except Exception:
        pass
    try:
        alpha_candidates = ["AlphaMissense", "ALPHAMISSENSE", "ALPHA_MISSENSE", "alpha_missense", 
                        "ALPHA-MISSENSE", "AlphaMissense_score", "ALPHAMISSENSE_SCORE"]
        
        for alpha_key in alpha_candidates:
            alpha_val = ann.get(alpha_key)
            if alpha_val is not None:
                alpha_str = str(alpha_val).strip().upper()  # Convert to uppercase for consistency
                if alpha_str and alpha_str.lower() not in ["", ".", "na", "none", "null"]:
                    # Try to parse as float first
                    try:
                        vr.alpha_missense = float(alpha_str)
                        vr.alpha_missense_category = None
                        logger.debug(f"Parsed AlphaMissense as float: {vr.alpha_missense}")
                        break
                    except ValueError:
                        # Not a float - must be categorical (P, B, A, etc.)
                        vr.alpha_missense_category = alpha_str
                        vr.alpha_missense = None
                        logger.debug(f"Parsed AlphaMissense as category: {vr.alpha_missense_category}")
                        break
    except Exception as e:
        logger.debug(f"Failed to parse AlphaMissense: {e}")
        pass
    try:
        if ann.get("SPLICE_AI_SCORE") is not None:
            vr.spliceai = float(ann.get("SPLICE_AI_SCORE"))
        elif ann.get("SPLICEAI") is not None:
            vr.spliceai = float(ann.get("SPLICEAI"))
    except Exception:
        pass
    # try:
    #     if ann.get("CADD") is not None:
    #         vr.cadd = float(ann.get("CADD"))
    # except Exception:
    #     pass

    # # Parse SQUIRLS (can be numeric score or category like 'N', 'D', 'P', 'B', 'T')
    # try:
        
    #     if ann.get("SQUIRLS") is not None:
    #         sq_val = ann.get("SQUIRLS")
    #         # Try to parse as float first
    #         try:
    #             vr.squirls_score = float(sq_val)
    #         except (ValueError, TypeError):
    #             # Not a float - must be category (N, D, P, B, T)
    #             vr.squirls_category = str(sq_val).strip().upper()
    #     else:
    #         # Only print for first few variants to avoid spam
    #         if not hasattr(create_variant_from_extract_line, "_debug_no_squirls_count"):
    #             create_variant_from_extract_line._debug_no_squirls_count = 0
    #         if create_variant_from_extract_line._debug_no_squirls_count < 5:
    #             create_variant_from_extract_line._debug_no_squirls_count += 1
    # except Exception as e:
    #     print(f"DEBUG: SQUIRLS error: {e}", flush=True)
    
    # try:
    #     if ann.get("SQUIRLS_SCORE") is not None:
    #         sq_score_val = ann.get("SQUIRLS_SCORE")
    #         try:
    #             vr.squirls_score = float(sq_score_val)
    #         except (ValueError, TypeError):
    #             print(f"DEBUG: SQUIRLS_SCORE not a float: {sq_score_val}", flush=True)
    # except Exception as e:
    #     print(f"DEBUG: SQUIRLS_SCORE error: {e}", flush=True)

    try:
        # Try multiple possible field names
        dscore_candidates = ["DSCORE", "Dscore", "D_SCORE", "D-SCORE", "DSPLICE_SCORE", "D_SPLICE"]
        for dkey in dscore_candidates:
            dscore_val = ann.get(dkey)
            if dscore_val is not None:
                try:
                    vr.dscore = float(dscore_val)
                    break
                except (ValueError, TypeError):
                    # If not float, try to extract numeric part
                    dscore_str = str(dscore_val).strip()
                    dscore_match = re.search(r'(\d+\.?\d*)', dscore_str)
                    if dscore_match:
                        vr.dscore = float(dscore_match.group(1))
                        break
                    pass
    except Exception as e:
        logger.debug(f"Failed to parse Dscore: {e}")

    # Parse gnomAD frequencies - take MINIMUM (ignoring -1 values)
    gnomad_candidates = ["GNOMAD_AF", "GNOMAD_GENOME", "MAX_FREQ", "gnom_genome", "max_freq", "GNOM_GENOME"]
    gnomad_vals = []
    for key in gnomad_candidates:
        val = ann.get(key)
        if val is not None:
            try:
                f = float(val)
                # CRITICAL: Ignore -1, -1.0, and any negative placeholder values
                if f >= 0:  # Only consider non-negative values (valid allele frequencies)
                    gnomad_vals.append(f)
            except (ValueError, TypeError):
                pass

    if gnomad_vals:
        vr.gnomad_af = min(gnomad_vals)  # Take MINIMUM of valid values only
        vr.gnomad_details = {"source": "VCF_INFO", "values": gnomad_vals, "type": "min"}
    else:
        vr.gnomad_af = 0.0  
        vr.gnomad_details = {"source": "VCF_INFO", "type": "no_valid_data"}

    # exon field if present
    vr.exon = ann.get("EXON") or ann.get("EXON_NUMBER") or ""

    # parse sample columns
    sample_parsed = parse_sample_field(sample_field)
    vr.proband_gt = sample_parsed.get("gt", "./.")
    vr.proband_dp = sample_parsed.get("dp")
    vr.proband_ad = None
    if sample_parsed.get("ad") is not None:
        ad = sample_parsed.get("ad")
        if isinstance(ad, (list, tuple)):
            if len(ad) >= 2:
                vr.proband_ad = (ad[0], ad[1])
            elif len(ad) == 1:
                vr.proband_ad = (ad[0], 0)
    vr.proband_af = sample_parsed.get("af")

    # parse AA ref, pos, alt from HGVSp using existing parse_hgvsp_three_letter function
    aa_ref, aa_pos, aa_alt = parse_hgvsp_extended(vr.hgvsp)
    vr.aa_ref = aa_ref
    vr.aa_pos = aa_pos
    vr.aa_alt = aa_alt
    
    return vr

def _is_het(gt: str) -> bool:
    return gt in ("0/1","1/0","0|1","1|0")

def _is_hom(gt: str) -> bool:
    return gt in ("1/1","1|1")

def _is_wt(gt: str) -> bool:
    return gt in ("0/0","0|0")

def clean_hgvsp(hgvsp: str) -> str:
    """Clean HGVSp string for comparison."""
    if not hgvsp:
        return ""
    # Remove everything before the last colon
    if ":" in hgvsp:
        parts = hgvsp.split(":")
        if len(parts) > 1:
            hgvsp = parts[-1]
    # Remove leading "p." if present
    hgvsp = re.sub(r'^\s*p\.', '', hgvsp.strip(), flags=re.IGNORECASE)
    # Remove parentheses if present
    hgvsp = re.sub(r'^\((.*)\)$', r'\1', hgvsp)
    return hgvsp.strip()

# Valid three-letter amino acid names (capitalized like "Pro")
VALID_THREE = {
    "Ala","Arg","Asn","Asp","Cys","Glu","Gln","Gly","His","Ile",
    "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val",
    "Sec","Pyl","Ter"
}

def _normalize_to_three_letter(aa: str) -> Optional[str]:
    """
    Normalize an amino-acid token to a three-letter code (capitalized, e.g. 'Pro').
    Accepts:
      - three-letter forms (any case): 'Pro', 'pro', 'PRO'
      - stop/termination tokens: '*', 'Ter', 'STOP' -> returns 'Ter'
    Returns None if token cannot be interpreted as a valid three-letter AA.
    """
    if not aa:
        return None
    a = aa.strip()
    if a == "*":
        return "Ter"
    a_norm = a.capitalize()
    if a_norm.lower() in ("ter", "stop"):
        return "Ter"
    if a_norm in VALID_THREE:
        return a_norm
    # try title case as a fallback (handles e.g. 'SEc' -> 'Sec')
    a_title = a.title()
    if a_title in VALID_THREE:
        return a_title
    return None

def parse_hgvsp_three_letter(hgvsp: str) -> Tuple[Optional[str], Optional[int], Optional[str]]:
    """
    Parse an HGVSp value and return (aa_ref_three, aa_pos, aa_alt_three).
    - aa_ref_three and aa_alt_three are three-letter codes (e.g. 'Pro', 'Ile', 'Ter') or None.
    - aa_pos is an integer or None.
    Supported simple forms:
      - p.Pro404Pro
      - p.Trp26Ter
      - p.Trp26*
      - p.Leu1199  (returns ref + pos, alt=None)
    Returns (None, None, None) for 'Not provided', '(=)', unsupported forms, or when parsing fails.
    """
    if not hgvsp:
        return (None, None, None)
    h = hgvsp.strip()
    # common "not provided" markers
    if h.lower() in ("not provided", "not_provided", ".", "-", ""):
        return (None, None, None)

    # Use existing clean_hgvsp if present in globals to preserve current cleaning behavior
    try:
        h_clean = clean_hgvsp(h) if 'clean_hgvsp' in globals() else h
    except Exception:
        h_clean = h

    # remove leading "p." and unwrap parentheses
    s = re.sub(r'^\s*p\.', '', h_clean, flags=re.IGNORECASE).strip()
    s = re.sub(r'^\((.*)\)$', r'\1', s)

    # ignore synonymous notation like p.(=)
    if s in ("(=)", "=", "?", "p.(=)"):
        return (None, None, None)

    # Case: REF POS ALT where REF and ALT are three-letter codes or '*' (e.g. "Pro404Pro")
    m = re.match(r'^(?P<ref>[A-Za-z]{3}|\*)\s*(?P<pos>\d+)\s*(?P<alt>[A-Za-z]{3}|\*)$', s, flags=re.IGNORECASE)
    if m:
        ref_raw = m.group("ref")
        pos_raw = m.group("pos")
        alt_raw = m.group("alt")
        try:
            pos = int(pos_raw)
        except Exception:
            pos = None
        ref_3 = _normalize_to_three_letter(ref_raw)
        alt_3 = _normalize_to_three_letter(alt_raw)
        return (ref_3, pos, alt_3)

    # Case: REFposAlt where Alt might be '*', 'Ter', 'Stop', 'X' (e.g. "Trp26*")
    m2 = re.match(r'^(?P<ref>[A-Za-z]{3})(?P<pos>\d+)(?P<alt>\*|Ter|Stop|X)$', s, flags=re.IGNORECASE)
    if m2:
        ref_raw = m2.group("ref")
        pos = int(m2.group("pos"))
        alt_raw = m2.group("alt")
        ref_3 = _normalize_to_three_letter(ref_raw)
        alt_3 = _normalize_to_three_letter(alt_raw)
        return (ref_3, pos, alt_3)

    # Case: REFpos (no ALT) e.g. "Leu1199"
    m3 = re.match(r'^(?P<ref>[A-Za-z]{3}|\*)(?P<pos>\d+)$', s, flags=re.IGNORECASE)
    if m3:
        ref_3 = _normalize_to_three_letter(m3.group("ref"))
        pos = int(m3.group("pos"))
        return (ref_3, pos, None)

    # Unsupported or complex HGVSp (indels, frameshifts, ranges, splice, etc.)
    return (None, None, None)

def parse_hgvsp_extended(hgvsp: str) -> Tuple[Optional[str], Optional[int], Optional[str]]:
    """
    Parse HGVSp string to extract AA_ref, AA_pos, AA_alt.
    Returns ONLY actual amino acid codes (three-letter) or Ter.
    For frameshift/deletions, AA_alt is None.
    """
    if not hgvsp or hgvsp.lower() in ["not provided", ".", "-", ""]:
        return None, None, None
    
    # Clean the HGVSp string
    hgvsp_clean = clean_hgvsp(hgvsp)
    
    # 1. Try standard three-letter parser first (missense/nonsense)
    ref, pos, alt = parse_hgvsp_three_letter(hgvsp_clean)
    if ref is not None or pos is not None:
        # For standard substitutions, return as-is
        return ref, pos, alt
    
    hgvsp_lower = hgvsp_clean.lower()
    
    # 2. Handle delins variants: p.Tyr2236_Arg2239delinsTer
    # Extract the INSERTED amino acid as AA_alt
    if "delins" in hgvsp_lower:
        # Pattern: p.Tyr2236_Arg2239delinsTer
        match = re.search(r'p\.([A-Za-z]{3})(\d+)_[A-Za-z]{3}\d+delins([A-Za-z\*]{1,3})', hgvsp_clean, re.IGNORECASE)
        if not match:
            # Pattern: p.Leu4delinsTer (single AA delins)
            match = re.search(r'p\.([A-Za-z]{3})(\d+)delins([A-Za-z\*]{1,3})', hgvsp_clean, re.IGNORECASE)
        
        if match:
            ref_three = _normalize_to_three_letter(match.group(1))
            try:
                pos = int(match.group(2))
            except:
                pos = None
            # Get the inserted amino acid
            ins_aa = match.group(3) if len(match.groups()) >= 3 else None
            alt_three = _normalize_to_three_letter(ins_aa) if ins_aa else None
            return ref_three, pos, alt_three
    
    # 3. For frameshifts/deletions/insertions: return AA_ref and pos, but AA_alt = None
    # The type will be determined from consequence field
    
    # Frameshift: p.Ile540fs
    if "fs" in hgvsp_lower and "delins" not in hgvsp_lower:
        match = re.search(r'p\.([A-Za-z]{3})(\d+)', hgvsp_clean, re.IGNORECASE)
        if match:
            ref_three = _normalize_to_three_letter(match.group(1))
            try:
                pos = int(match.group(2))
            except:
                pos = None
            return ref_three, pos, None  
    
    # Deletion: p.Lys1617del ot p.Leu4_Gln6del
    if "del" in hgvsp_lower and "delins" not in hgvsp_lower:
        # Get first AA and position
        match = re.search(r'p\.([A-Za-z]{3})(\d+)', hgvsp_clean, re.IGNORECASE)
        if match:
            ref_three = _normalize_to_three_letter(match.group(1))
            try:
                pos = int(match.group(2))
            except:
                pos = None
            return ref_three, pos, None  
    
    # Insertion: p.Lys1617_Leu1618insAla
    if "ins" in hgvsp_lower and "delins" not in hgvsp_lower:
        # Get the AA before insertion
        match = re.search(r'p\.([A-Za-z]{3})(\d+)_', hgvsp_clean, re.IGNORECASE)
        if match:
            ref_three = _normalize_to_three_letter(match.group(1))
            try:
                pos = int(match.group(2))
            except:
                pos = None
            return ref_three, pos, None
    
    # Duplication: p.Lys1617dup
    if "dup" in hgvsp_lower:
        match = re.search(r'p\.([A-Za-z]{3})(\d+)', hgvsp_clean, re.IGNORECASE)
        if match:
            ref_three = _normalize_to_three_letter(match.group(1))
            try:
                pos = int(match.group(2))
            except:
                pos = None
            return ref_three, pos, None  
    
    return None, None, None

def parse_protein_pos(hgvsp: str) -> Optional[int]:
    """
    Return amino-acid position parsed from HGVSp or None.
    Uses parse_hgvsp_three_letter to avoid duplicated logic.
    """
    if not hgvsp:
        return None
    ref, pos, alt = parse_hgvsp_three_letter(hgvsp)
    return pos

def extract_aa_position_from_hgvsp(hgvsp: str) -> Optional[int]:
    """
    Compatibility wrapper used in PM5/ClinVar matching.
    Returns amino-acid position parsed from HGVSp or None.
    """
    try:
        return parse_protein_pos(hgvsp)
    except Exception:
        return None

def parse_protein_change_three_letter(hgvsp: str) -> Tuple[Optional[str], Optional[int], Optional[str]]:
    """
    Parse HGVSp string to extract AA_ref, AA_pos, AA_alt in three-letter codes.
    Returns (aa_ref_three, aa_pos, aa_alt_three)
    """
    if not hgvsp:
        return None, None, None
    
    # Clean the HGVSp string
    hgvsp_clean = clean_hgvsp(hgvsp)
    
    # Skip common non-informative values
    if hgvsp_clean.lower() in ["not provided", ".", "-", "", "p.(=)", "p.?"]:
        return None, None, None
    
    # Use existing function that returns three-letter codes
    return parse_hgvsp_three_letter(hgvsp_clean)

# --- Parse cDNA pos from HGVSc and map cDNA position -> exon number using nmd_map / tx_map ---
def parse_cdna_from_hgvsc(hgvsc: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse cDNA coordinate from an HGVSc string.
    Examples:
      "ENST00000258415.9:c.646+1G>T" -> returns (646, +1)
      "c.123-2A>G" -> returns (123, -2)
      "c.876A>T" -> returns (876, None)
    Returns (cdna_pos, intronic_offset) where intronic_offset is None if absent,
    else a signed int (+N or -N).
    """
    if not hgvsc:
        return None, None
    s = str(hgvsc).strip()
    # If there is an ENST: prefix, split on ":"
    if ":" in s:
        parts = s.split(":")
        # find the part that starts with c.
        found = None
        for p in parts[::-1]:
            if p.startswith("c.") or p.startswith("C."):
                found = p
                break
        if not found:
            # fallback to entire string
            found = parts[-1]
    else:
        found = s

    # Parse patterns like c.646+1, c.876
    m = re.search(r'c\.(\d+)(?:([+-])(\d+))?', found, flags=re.I)
    if not m:
        return None, None
    try:
        cdna_pos = int(m.group(1))
    except Exception:
        cdna_pos = None
    intr = None
    if m.group(2) and m.group(3):
        try:
            sign = 1 if m.group(2) == '+' else -1
            intr = sign * int(m.group(3))
        except Exception:
            intr = None
    return cdna_pos, intr

def map_cdna_to_exon_number(transcript: str,
                            cdna_pos: Optional[int],
                            nmd_map: Dict[str, Dict[str, Any]],
                            tx_map: Optional[Dict[str, Any]] = None,
                            intronic_offset: Optional[int] = None
                            ) -> Optional[int]:
    """
    Map cDNA coordinate (cdna_pos) to an exon number using nmd_map.
    Accepts optional intronic_offset (signed int) from HGVSc (e.g., +1, -2).
    Improvements:
     - If cds_exon_length missing, try to infer exon cds length from genomic_start/genomic_end when available.
     - If intronic_offset provided, interpret c.(N)+k / c.(N)-k as adjacent exon boundaries:
         * For +k (intronic after cdna_pos) prefer the exon whose cumulative end == cdna_pos
         * For -k (intronic before cdna_pos) prefer exon whose cumulative start == cdna_pos+1
    Returns exon_number or None.
    """
    if not cdna_pos or cdna_pos <= 0:
        return None
    if not nmd_map:
        return None

    def try_map_on_entry(entry: Dict[str, Any]) -> Optional[int]:
        exs = entry.get("exons", []) or []
        cumulative = 0
        # build list of exon lengths, but if cds_exon_length missing use genomic coords fallback
        exon_entries = []
        for e in exs:
            # attempt to get cds_exon_length
            clen = None
            try:
                clen = e.get("cds_exon_length") or e.get("cds_length") or e.get("cds_length_exon")
                if clen is not None:
                    clen = int(clen)
                    if clen <= 0:
                        clen = None
            except Exception:
                clen = None

            # fallback: compute from genomic coords if both present
            if clen is None:
                try:
                    gs = e.get("genomic_start") or e.get("genomicStart") or e.get("start")
                    ge = e.get("genomic_end") or e.get("genomicEnd") or e.get("end")
                    if gs is not None and ge is not None:
                        gs_i = int(gs); ge_i = int(ge)
                        # approximate exon length (genomic) as coding length only if exon is coding — best-effort
                        clen = abs(ge_i - gs_i) + 1
                except Exception:
                    clen = None

            exon_entries.append({"exon_number": e.get("exon_number"), "cds_exon_length": clen,
                                 "genomic_start": e.get("genomic_start"), "genomic_end": e.get("genomic_end")})

        # If any exon still lacks length, we will attempt ordinal mapping by order
        # Now accumulate
        cumulative = 0
        for e in exon_entries:
            clen = e.get("cds_exon_length") or 0
            # if we do not have a reliable cds_exon_length for this exon, mark and break to fallback to order-based mapping
            if clen <= 0:
                cumulative = None
                break
            prev = cumulative
            cumulative += int(clen)
            # direct hit for a coding base inside this exon
            if prev < cdna_pos <= cumulative:
                try:
                    return int(e.get("exon_number")) if e.get("exon_number") is not None else None
                except Exception:
                    return None

        # If we couldn't map by cds_exon_length, attempt ordinal mapping using genomic coordinates if present:
        if cumulative is None:
            # try to find an exon whose genomic interval contains the variant genomic position —
            # but we do not have genomic variant position here; caller may handle that fallback.
            # Try to use approximate logic: if intronic_offset present, interpret cdna boundary:
            # For intronic offsets, prefer the exon at the boundary:
            if intronic_offset is not None:
                # For +N (intronic after cdna_pos) the exon with cumulative end == cdna_pos is the donor exon.
                # For -N (intronic before cdna_pos) the exon with cumulative start == cdna_pos is acceptor exon.
                # We'll attempt a best-effort by scanning exon_entries and using cumulative counts inferred from exon order if possible.
                running = 0
                for e in exon_entries:
                    clen = e.get("cds_exon_length")
                    # when length missing, try genomic length fallback
                    if not clen:
                        try:
                            gs = e.get("genomic_start"); ge = e.get("genomic_end")
                            if gs is not None and ge is not None:
                                clen = abs(int(ge) - int(gs)) + 1
                            else:
                                # cannot compute further
                                clen = None
                        except Exception:
                            clen = None
                    if not clen:
                        # skip if unable to derive
                        continue
                    prev = running
                    running += int(clen)
                    if intronic_offset > 0:
                        # donor side: if this exon ends at cdna_pos (running == cdna_pos) -> choose it
                        if running == cdna_pos:
                            try:
                                return int(e.get("exon_number")) if e.get("exon_number") is not None else None
                            except Exception:
                                return None
                    else:
                        # intronic_offset < 0: acceptor side — check if prev+1 == cdna_pos
                        if prev + 1 == cdna_pos:
                            try:
                                return int(e.get("exon_number")) if e.get("exon_number") is not None else None
                            except Exception:
                                return None
            # nothing reliable
            return None

        return None  

    # 1) direct transcript key attempts
    tx_candidates = []
    if transcript:
        tx_full = transcript.strip()
        tx_base = tx_full.split(".")[0] if "." in tx_full else tx_full
        tx_candidates.append(tx_full)
        tx_candidates.append(tx_base)

    for key in list(nmd_map.keys()):
        if tx_candidates and key in tx_candidates:
            res = try_map_on_entry(nmd_map[key])
            if res:
                return res

    # 2) if transcript is Ensembl (ENST), try to map to RefSeq via tx_map
    if tx_map and isinstance(tx_map, dict):
        enst_key = transcript.split(".")[0] if transcript else None
        if enst_key:
            try:
                refseqs = tx_map.get("ensembl_to_refseq", {}).get(enst_key) or tx_map.get("ensembl_to_refseq", {}).get(transcript)
                if refseqs:
                    for nm in refseqs:
                        if nm in nmd_map:
                            res = try_map_on_entry(nmd_map[nm])
                            if res:
                                return res
            except Exception:
                pass

    # 3) try to match by transcript base 
    tx_base = transcript.split(".")[0] if transcript else None
    if tx_base:
        for k in nmd_map.keys():
            if k.split(".")[0] == tx_base:
                res = try_map_on_entry(nmd_map[k])
                if res:
                    return res

    # Could not map
    return None


# --- Precise canonical splice-site check (+/-1-2 nt) ---
def _is_canonical_splice_site(vr: VariantRecord, 
                              nmd_map: Dict[str, Dict[str, Any]] = None) -> bool:
    """
    Determine if a splice variant is at a canonical splice site (+1/+2 for donor, -1/-2 for acceptor)
    using exon boundary information from nmd_map (RefSeq transcripts only).
    
    Canonical splice sites are:
      - Donor site (5' splice site): first 2 nucleotides of intron (position +1, +2 relative to exon end)
      - Acceptor site (3' splice site): last 2 nucleotides of intron (position -1, -2 relative to next exon start)
    
    Args:
        vr: VariantRecord with consequence, hgvsc, chrom, pos, transcript
        nmd_map: Pre-loaded NMD/exon table from --nmd-table (keyed by RefSeq transcript ID like NM_000138.5)
    
    Returns:
        True if variant is at canonical splice position, False otherwise
    """
    cons = (getattr(vr, "consequence", "") or "").lower()
    
    # Only consider splice-related consequences
    if not any(tok in cons for tok in ("splice_acceptor", "splice_donor", "splice_region")):
        return False
    
    # If consequence explicitly says splice_acceptor or splice_donor (without "region"),
    # and we have no exon table, assume canonical (conservative fallback)
    is_explicit_acceptor_donor = ("splice_acceptor" in cons or "splice_donor" in cons) and "region" not in cons
    
    # Step 1: Get RefSeq transcript ID for this variant
    transcript_id = getattr(vr, "transcript", "") or getattr(vr, "tid", "") or ""
    
    # If no transcript in variant record, try to extract from HGVSc
    if not transcript_id:
        hgvsc = str(getattr(vr, "hgvsc", "") or "")
        if hgvsc and ":" in hgvsc:
            parts = hgvsc.split(":", 1)
            # Look for RefSeq pattern NM_ or NR_
            if parts[0].startswith(("NM_", "NR_", "XM_", "XR_")):
                transcript_id = parts[0]
    
    if not transcript_id:
        logger.debug(f"No transcript ID for {vr.chrom}:{vr.pos}, cannot determine canonical splice status")
        return is_explicit_acceptor_donor
    
    # Normalize transcript ID (remove version suffix for matching)
    tx_base = transcript_id.split(".")[0] if "." in transcript_id else transcript_id
    
    # Step 2: Find exon structure for this transcript in nmd_map
    exon_structure = None
    best_transcript = None
    
    # Try exact match first, then base version
    if nmd_map:
        if transcript_id in nmd_map:
            exon_structure = nmd_map[transcript_id]
            best_transcript = transcript_id
        elif tx_base in nmd_map:
            exon_structure = nmd_map[tx_base]
            best_transcript = tx_base
    
    if not exon_structure:
        logger.debug(f"No exon structure found for transcript {transcript_id} in nmd_map")
        # Fallback: if VEP says splice_acceptor/donor, assume canonical
        return is_explicit_acceptor_donor
    
    # Step 3: Get strand and exons list
    strand = exon_structure.get("strand", "+")
    exons = exon_structure.get("exons", [])
    
    if not exons:
        logger.debug(f"No exons in structure for {best_transcript}")
        return is_explicit_acceptor_donor
    
    # Sort exons by genomic position (need to handle negative strand correctly)
    # For positive strand: exons sorted by genomic_start ascending
    # For negative strand: exons sorted by genomic_start descending (biological order)
    if strand == "+" or strand == "1":
        exons_sorted = sorted(exons, key=lambda x: x.get("genomic_start", 0))
    else:
        exons_sorted = sorted(exons, key=lambda x: x.get("genomic_start", 0), reverse=True)
    
    variant_pos = int(vr.pos)
    
    # Step 4: Check each exon boundary for canonical splice positions
    for i, exon in enumerate(exons_sorted):
        genomic_start = exon.get("genomic_start")
        genomic_end = exon.get("genomic_end")
        
        if genomic_start is None or genomic_end is None:
            continue
        
        gs = int(genomic_start)
        ge = int(genomic_end)
        
        if strand == "+" or strand == "1":
            # Positive strand
            # Donor site (5' splice site): positions after exon end: +1, +2
            if variant_pos == ge + 1 or variant_pos == ge + 2:
                logger.debug(f"Canonical donor site detected: {vr.chrom}:{vr.pos} at +{variant_pos - ge} from exon {exon.get('exon_number')}")
                return True
            
            # Acceptor site (3' splice site): positions before next exon start: -1, -2
            if i + 1 < len(exons_sorted):
                next_exon_start = int(exons_sorted[i + 1].get("genomic_start"))
                if variant_pos == next_exon_start - 1 or variant_pos == next_exon_start - 2:
                    logger.debug(f"Canonical acceptor site detected: {vr.chrom}:{vr.pos} at -{next_exon_start - variant_pos} from exon {exon.get('exon_number')}")
                    return True
        else:
            # Negative strand (transcribed from reverse strand)
            # Donor site: positions before exon start (since transcription goes opposite direction)
            if variant_pos == gs - 1 or variant_pos == gs - 2:
                logger.debug(f"Canonical donor site detected (negative strand): {vr.chrom}:{vr.pos}")
                return True
            
            # Acceptor site: positions after exon end
            if i + 1 < len(exons_sorted):
                next_exon_end = int(exons_sorted[i + 1].get("genomic_end"))
                if variant_pos == next_exon_end + 1 or variant_pos == next_exon_end + 2:
                    logger.debug(f"Canonical acceptor site detected (negative strand): {vr.chrom}:{vr.pos}")
                    return True
    
    # Step 5: Try using HGVSc offset as additional fallback
    hgvsc = str(getattr(vr, "hgvsc", "") or "")
    if hgvsc:
        cdna_pos, intronic_offset = _parse_cdna_from_hgvsc(hgvsc)
        if cdna_pos is not None and intronic_offset is not None:
            # Canonical offsets are 1 or 2 (absolute value)
            if abs(intronic_offset) in (1, 2):
                logger.debug(f"Canonical splice site inferred from HGVSc offset: {intronic_offset}")
                return True
    
    # Step 6: Not a canonical splice site
    if is_explicit_acceptor_donor:
        logger.warning(f"Variant {vr.chrom}:{vr.pos} has splice_acceptor/donor annotation but not at canonical position per exon boundaries")
    
    return False

def _is_noncanonical_splice_region_variant(vr: VariantRecord) -> bool:
    """
    Identify non-canonical splice region variants that should NOT get PP3.
    
    Returns True for variants that:
      1. Have 'splice_region' in consequence
      2. Do NOT have 'splice_acceptor' or 'splice_donor'
      3. Have intronic offset >= 3 from exon boundary (e.g., c.6313+5G>A)
    
    These variants should not receive PP3 because SpliceAI scores for
    positions beyond +/-2 are often unreliable for pathogenicity prediction.
    """
    cons = (getattr(vr, "consequence", "") or "").lower()
    
    # Must have splice_region but NOT canonical splice acceptor/donor
    if "splice_region" not in cons:
        return False
    
    if "splice_acceptor" in cons or "splice_donor" in cons:
        return False
    
    # Parse HGVSc to get intronic offset
    hgvsc = str(getattr(vr, "hgvsc", "") or "")
    if hgvsc:
        # Pattern: c.NNN+NN or c.NNN-NN
        match = re.search(r'c\.[^\s:+-]*([+-])(\d+)', hgvsc)
        if match:
            try:
                offset = int(match.group(2))
                # Canonical are only 1 or 2
                if offset in (1, 2):
                    return False
                # Offsets >= 3 are non-canonical
                if offset >= 3:
                    return True
            except Exception:
                pass
    
    # If we can't determine offset but it's splice_region without acceptor/donor,
    # treat as non-canonical (conservative)
    return True

def _parse_cdna_from_hgvsc(hgvsc: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse cDNA coordinate and intronic offset from HGVSc string.
    
    Examples:
      "c.646+1G>T" -> returns (646, 1)
      "c.123-2A>G" -> returns (123, -2)
      "c.876A>T" -> returns (876, None)
    
    Returns:
      (cdna_position, intronic_offset) where intronic_offset is None if not an intronic variant
    """
    if not hgvsc:
        return None, None
    
    s = str(hgvsc).strip()
    
    # Handle transcript prefix like "NM_000069.3:c.3716G>A"
    if ":" in s:
        parts = s.split(":")
        for p in parts:
            if p.startswith("c.") or p.startswith("C."):
                s = p
                break
    
    # Parse patterns: c.646+1G>T, c.123-2A>G, c.876A>T
    match = re.search(r'c\.(\d+)(?:([+-])(\d+))?', s, flags=re.I)
    if not match:
        return None, None
    
    try:
        cdna_pos = int(match.group(1))
    except (ValueError, TypeError):
        cdna_pos = None
    
    intronic_offset = None
    if match.group(2) and match.group(3):
        try:
            sign = 1 if match.group(2) == '+' else -1
            intronic_offset = sign * int(match.group(3))
        except (ValueError, TypeError):
            intronic_offset = None
    
    return cdna_pos, intronic_offset

# ---------------------------
# gnomAD AF aggregation utility
# ---------------------------
GNOMAD_KEY_RE = re.compile(r"(?i)gnomad|gnomad_r|gnomad\.|max_af|max_freq|af_popmax|gnom_genome|gnomad_genome|gnom3g|gnom3G")

def extract_gnomad_af_from_info(info: Dict[str, Any]) -> Tuple[Dict[str,Dict[str,float]], Dict[str,Any]]:
    """
    Extract gnomAD allele frequencies from INFO field.
    Returns per-version frequencies and GLOBAL MINIMUM AF (most conservative for PM2).
    Ignores placeholder values like -1, -1.0, 999, etc.
    """
    per_version = defaultdict(dict)
    global_min = {"min_af": 1.0, "pop": None, "version": None}
    all_afs = []
    
    for k, v in (info or {}).items():
        if not isinstance(k, str):
            continue
        if not GNOMAD_KEY_RE.search(k):
            continue
        kl = k.lower()
        version = "unknown"
        if re.search(r"\b(r|v)?2\b|2\.1", kl):
            version = "v2"
        elif re.search(r"\b(r|v)?3\b|3\.1", kl):
            version = "v3"
        elif re.search(r"\b(r|v)?4\b|4\.1", kl):
            version = "v4"
        pop = "global"
        mpop = re.search(r"(afr|amr|eas|sas|nfe|fin|asj|oth|popmax|pmax|global|adj|overall)", kl)
        if mpop:
            pop = mpop.group(1)
        
        vals = []
        try:
            if isinstance(v, (list,tuple)):
                for x in v:
                    try: 
                        f = float(x)
                        # CRITICAL: Ignore -1, -1.0, and any negative placeholder values
                        if f >= 0 and f <= 1:  # Valid allele frequency range
                            vals.append(f)
                    except Exception: 
                        pass
            else:
                s = str(v)
                for sep in ("/", ",", ";", "&", "|"):
                    if sep in s:
                        parts = s.split(sep)
                        break
                else:
                    parts = [s]
                for p in parts:
                    try: 
                        f = float(p)
                        # CRITICAL: Ignore -1, -1.0, and any negative placeholder values
                        if f >= 0 and f <= 1:  # Valid allele frequency range
                            vals.append(f)
                    except Exception: 
                        pass
        except Exception:
            vals = []
        
        if not vals:
            continue
        
        # Take MINIMUM value for this field (most conservative)
        val = min(vals)
        all_afs.append(val)
        
        # Store per-version per-population minimum
        if val < per_version[version].get(pop, 1.0):
            per_version[version][pop] = val
        
        # Update global minimum
        if val < global_min["min_af"]:
            global_min = {"min_af": val, "pop": pop, "version": version}
    
    # If no valid values found, set min_af to 1.0 (will be treated as no data)
    if not all_afs:
        global_min["min_af"] = 1.0
    
    logger.debug(f"gnomAD AF extraction: found {len(all_afs)} valid values, min={min(all_afs) if all_afs else 'None'}")
    
    return dict(per_version), global_min

# ---------------------------
# Gene validation and DatabaseManager (ClinVar/HGMD/Internal)
# ---------------------------
def validate_gene_match_db_evidence(variant_gene: str, db_gene_info: str, db_source: str) -> Tuple[bool, str]:
    variant_gene_norm = normalize_gene_name(variant_gene)
    if not variant_gene_norm:
        return False, f"Variant has no gene information"
    if not db_gene_info or db_gene_info.strip() in ["", ".", "Not provided", "not_provided", "not_specified"]:
        return True, f"No gene information in {db_source} - validation skipped"
    if db_source == "ClinVar":
        db_genes = []
        for gene_pair in str(db_gene_info).split('|'):
            if ':' in gene_pair:
                gene_raw = gene_pair.split(':')[0]
                gene_norm = normalize_gene_name(gene_raw)
                db_genes.append(gene_norm)
        if not db_genes:
            return True, f"No parsable gene info in ClinVar"
        if variant_gene_norm in db_genes:
            return True, f"Gene match in ClinVar: {variant_gene_norm}"
        else:
            return False, f"ClinVar gene mismatch: {variant_gene_norm} not in {db_genes}"
    elif db_source == "HGMD":
        db_gene_norm = normalize_gene_name(db_gene_info)
        if not db_gene_norm:
            return True, f"No gene info in HGMD - validation skipped"
        if variant_gene_norm == db_gene_norm:
            return True, f"Gene match in HGMD: {variant_gene_norm}"
        else:
            return False, f"HGMD gene mismatch: {variant_gene_norm} vs {db_gene_norm}"
    elif db_source == "Internal_DB":
        return True, f"Internal DB lacks gene information - validation skipped"
    else:
        return True, f"Unknown database {db_source} - validation skipped"

class DatabaseManager:
    """
    Centralized manager for external genomic databases (ClinVar, HGMD, gnomAD, internal DB).
    
    This class handles:
      - Loading and indexing of ClinVar VCF, HGMD CSV, internal DB CSV, and gnomAD VCF/TSV files.
      - Querying variant evidence (ClinVar significance, HGMD class, internal classifications).
      - Annotating gnomAD allele frequencies and population metrics.
      - Supporting protein-level queries for PS1/PM5 criteria (same-residue pathogenic variants).
    
    Attributes:
        clinvar_vcf (pysam.VariantFile): Opened ClinVar VCF file.
        hgmd_index (dict): Nested dictionary for HGMD variant lookups.
        gnomad_vcf_handles (dict): Opened gnomAD VCFs by version (v2, v3, v4).
        internal_db_index (dict): Index for internal laboratory database.
        clinvar_protein_index (dict): Pre-indexed ClinVar protein changes for PS1/PM5.
    """
    def __init__(self, db_paths: Dict[str,str], 
                 gnomad_v2: Optional[str]=None, gnomad_v3: Optional[str]=None, gnomad_v4: Optional[str]=None, 
                 internal_db_path: Optional[str]=None):
        self.paths = dict(db_paths or {})
        self.clinvar_vcf: Optional[pysam.VariantFile] = None
        self.hgmd_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.clinvar_protein_index: Dict[str, List[dict]] = defaultdict(list)
        self.clinvar_submission_index_by_vid = {}
        self.clinvar_protein_index_extended = {}
        self.internal_db_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.gnomad_vcf_handles: Dict[str, pysam.VariantFile] = {}
        self.gnomad_tsv_handles: Dict[str, Dict[str,Any]] = {}
        self._open_clinvar()
        self._load_hgmd_csv()
        self._open_gnomad_versions(gnomad_v2, gnomad_v3, gnomad_v4)
        if internal_db_path:
            self._load_internal_db(internal_db_path)
        exp_path = self.paths.get("CLINVAR_EXPANDED_CSV")
        if exp_path:
            try:
                self._load_clinvar_expanded_csv(exp_path)
            except Exception:
                logger.warning("Failed to load CLINVAR_EXPANDED_CSV: %s", exp_path)

    def _open_clinvar(self):
        p = self.paths.get("CLINVAR_VCF")
        if p and os.path.exists(p) and HAVE_PYSAM:
            try:
                self.clinvar_vcf = pysam.VariantFile(p)
                logger.info("Opened ClinVar VCF: %s", p)
                try: self._build_clinvar_indexes_from_vcf(p)
                except: pass
            except Exception as e:
                logger.warning("Failed to open ClinVar VCF: %s", e)

    def _detect_tsv_header(self, path: str) -> Tuple[List[str], bool]:
        openf = gzip.open if path.endswith((".gz", ".bgz")) else open
        cols = []
        header_hash = False
        try:
            with openf(path, "rt", encoding="utf-8", errors="replace") as fh:
                for _ in range(50):
                    ln = fh.readline()
                    if not ln:
                        break
                    ln = ln.rstrip("\n")
                    if ln.startswith("#"):
                        header_hash = True
                        ln2 = ln.lstrip("#").strip()
                        if ln2 and "\t" in ln2:
                            cols = [c.strip() for c in ln2.split("\t")]
                            break
                        else:
                            continue
                    else:
                        if "\t" in ln:
                            cols = [c.strip() for c in ln.split("\t")]
                        else:
                            cols = ln.split()
                        break
        except Exception as e:
            logger.debug("Failed to read TSV header from %s: %s", path, e)
        return cols, header_hash

    def _open_gnomad_versions(self, v2: Optional[str], v3: Optional[str], v4: Optional[str]):
        for ver, path in (("v2", v2), ("v3", v3), ("v4", v4)):
            if not path:
                continue
            if not os.path.exists(path):
                logger.warning("gnomAD %s path not found: %s", ver, path)
                continue
            opened = False
            if HAVE_PYSAM:
                try:
                    vf = pysam.VariantFile(path)
                    self.gnomad_vcf_handles[ver] = vf
                    logger.info("Opened gnomAD %s as VCF: %s", ver, path)
                    opened = True
                except Exception:
                    opened = False
            if opened:
                continue
            if HAVE_PYSAM:
                try:
                    tb = pysam.TabixFile(path)
                    cols, header_hash = self._detect_tsv_header(path)
                    if not cols:
                        cols = ["chrom","pos","ref","alt","af","ac","an"]
                    colmap = {c.lower(): i for i,c in enumerate(cols)}
                    self.gnomad_tsv_handles[ver] = {"tb": tb, "cols": colmap, "path": path, "header_hash": header_hash}
                    logger.info("Opened gnomAD %s as TSV (tabix): %s (cols: %s)", ver, path, list(colmap.keys())[:8])
                    opened = True
                except Exception as e:
                    logger.warning("Failed to open gnomAD %s as TSV: %s", ver, e)
            if not opened:
                logger.warning("Could not open gnomAD %s file %s as VCF or TSV", ver, path)

    def _parse_tsv_line_to_dict(self, ln: str, cols_map: Dict[str,int]) -> Dict[str,str]:
        parts = ln.rstrip("\n").split("\t")
        out = {}
        for name, idx in cols_map.items():
            try:
                out[name] = parts[idx] if idx < len(parts) else ""
            except Exception:
                out[name] = ""
        return out

    def _extract_af_from_tsv_record(self, rec_dict: Dict[str,str]) -> Tuple[Optional[float], Dict[str,Any]]:
        """
        Extract allele frequency from TSV record, taking MINIMUM value.
        """
        details = {}
        
        # Try common frequency field names
        for k in ("af","allele_freq","allele_frequency","freq","frequency","aaf","adj_af","af_global","af_raw"):
            if k in rec_dict and rec_dict[k] not in ("", "."):
                try:
                    # Parse as float
                    af = float(rec_dict[k])
                    details["source"] = k
                    # Return as is (will be compared for min later)
                    return af, details
                except Exception:
                    pass
        
        # Try population-specific fields and take MINIMUM
        pop_afs = []
        for k in list(rec_dict.keys()):
            if re.search(r"(afr|amr|eas|sas|nfe|fin|asj|oth)", k):
                v = rec_dict.get(k)
                if v and v not in (".", ""):
                    try:
                        af = float(v)
                        pop_afs.append(af)
                        details[f"pop_{k}"] = af
                    except Exception:
                        pass
        
        if pop_afs:
            min_af = min(pop_afs)
            details["source"] = "min_population"
            details["min_pop_af"] = min_af
            return min_af, details
        
        # Calculate from AC/AN
        ac = None
        an = None
        for k in ("ac","allele_count","AC"):
            if k.lower() in rec_dict and rec_dict.get(k.lower()) not in ("", "."):
                try:
                    ac = float(rec_dict.get(k.lower()))
                    break
                except Exception:
                    try:
                        ac = float(rec_dict.get(k.lower()).split(",")[0])
                        break
                    except Exception:
                        pass
        
        for k in ("an","allele_number","AN"):
            if k.lower() in rec_dict and rec_dict.get(k.lower()) not in ("", "."):
                try:
                    an = float(rec_dict.get(k.lower()))
                    break
                except Exception:
                    try:
                        an = float(rec_dict.get(k.lower()).split(",")[0])
                        break
                    except Exception:
                        pass
        
        if ac is not None and an and an != 0:
            af = ac / an
            details["source"] = "AC/AN"
            return af, details
        
        return None, details

    def annotate_gnomad_for_variant(self, vr: VariantRecord):
        """
        Annotate variant with gnomAD allele frequencies.
        Uses MINIMUM frequency across all versions and populations (most conservative for PM2).
        """
        # First try to get from VEP CSQ annotations in _raw_ann
        perv, gmin = extract_gnomad_af_from_info(getattr(vr, "_raw_ann", {}) or {})
        if gmin.get("min_af", 1.0) < 1.0:
            vr.gnomad_af = gmin["min_af"]
            vr.gnomad_details = {"version_pop": gmin, "source": "VEP_CSQ", "type": "min"}
            return
        
        chrom_clean = vr.chrom.replace("chr", "")
        chrom_vars = [chrom_clean, "chr" + chrom_clean]
        best = {"min_af": 1.0, "version": None, "pop": None, "orig": None}  # Start with 1.0, look for minimum
        
        # Query gnomAD VCFs
        for ver, vf in self.gnomad_vcf_handles.items():
            for c_query in chrom_vars:
                try:
                    for rec in vf.fetch(c_query, vr.pos - 1, vr.pos):
                        if int(rec.pos) != int(vr.pos): 
                            continue
                        if rec.ref != vr.ref: 
                            continue
                        if not (rec.alts and vr.alt in rec.alts): 
                            continue
                        
                        info_rec = dict(rec.info)
                        perv2, gmin2 = extract_gnomad_af_from_info(info_rec)
                        
                        # Take MINIMUM across all
                        if gmin2.get("min_af", 1.0) < best["min_af"]:
                            best.update({
                                "min_af": float(gmin2["min_af"]), 
                                "version": ver, 
                                "pop": gmin2.get("pop")
                            })
                        
                        # Also calculate from AC/AN if available
                        try:
                            ac = rec.info.get('AC')
                            an = rec.info.get('AN')
                            if ac and an and an > 0:
                                try:
                                    ac0 = float(ac[0]) if isinstance(ac, (list,tuple)) else float(ac)
                                    an0 = float(an[0]) if isinstance(an, (list,tuple)) else float(an)
                                    if ac0 >= 0 and an0 > 0:  # Ignore negative values
                                        afcalc = ac0 / an0
                                        if afcalc < best["min_af"]:
                                            best.update({"min_af": afcalc, "version": ver, "pop": "global_calc"})
                                except Exception:
                                    pass
                        except Exception:
                                    pass
                        
                    if best["min_af"] < 1.0: 
                        break
                except ValueError:
                    continue
                except Exception:
                    continue
            if best["min_af"] < 1.0:
                break
        
        # Query gnomAD TSV files
        for ver, obj in self.gnomad_tsv_handles.items():
            tb = obj.get("tb")
            colmap = obj.get("cols", {})
            if not tb: 
                continue
            for c_query in chrom_vars:
                try:
                    found_any = False
                    for ln in tb.fetch(c_query, vr.pos - 1, vr.pos):
                        found_any = True
                        recd = self._parse_tsv_line_to_dict(ln, colmap)
                        recd_l = {k.lower(): v for k, v in recd.items()}
                        
                        p_col = next((k for k in recd_l if k in ('pos','position')), None)
                        r_col = next((k for k in recd_l if k in ('ref','reference')), None)
                        a_col = next((k for k in recd_l if k in ('alt','alternate')), None)
                        
                        if p_col and int(recd_l[p_col]) != int(vr.pos): 
                            continue
                        if r_col and recd_l[r_col] != vr.ref: 
                            continue
                        if a_col and recd_l[a_col] != vr.alt: 
                            continue
                        
                        af_val, details = self._extract_af_from_tsv_record(recd_l)
                        
                        # Take MINIMUM frequency
                        if af_val is not None and af_val < best["min_af"]:
                            best.update({
                                "min_af": af_val, 
                                "version": ver, 
                                "pop": details.get("source")
                            })
                    
                    if found_any and best["min_af"] < 1.0: 
                        break
                except ValueError:
                    continue
                except Exception as e:
                    continue
        
        if best["min_af"] < 1.0:
            vr.gnomad_af = best["min_af"]
            vr.gnomad_details = {"version": best["version"], "pop": best["pop"], "type": "min"}
            logger.debug(f"gnomAD annotation for {vr.chrom}:{vr.pos}: min_af={vr.gnomad_af:.6f} from {best['version']}")
        else:
            # If no frequency found, set to 0.0 (not observed)
            vr.gnomad_af = 0.0
            vr.gnomad_details = {"type": "min", "source": "not_found"}

    def _load_clinvar_expanded_csv(self, path: Optional[str]) -> None:
        """
        Robust loader for clinvar_expanded CSV/TSV produced from submission_summary mapping.
        Builds:
        - self.clinvar_expanded_by_vid: {VariationID: [submission dicts]}
        - self.clinvar_expanded_by_qvar: {Query_Variant: [submission dicts]}
        Each submission dict minimally: {'phenotype','clnsig','review','submission_id','hgvsp','gene','stars'}
        This loader is tolerant to different column names and fallback parsing (pandas or csv).
        """
        if not path or not os.path.exists(path):
            logger.warning(f"ClinVar expanded CSV not found: {path}")
            return
        logger.info(f"Loading ClinVar expanded CSV: {path}")
        self.clinvar_expanded_by_vid = defaultdict(list)
        self.clinvar_expanded_by_qvar = defaultdict(list)

        def normalize_review_to_stars(review_text: str) -> int:
            if not review_text:
                return 0
            low = str(review_text).lower()
            if 'practice_guideline' in low or 'practice guideline' in low:
                return 4
            if 'expert_panel' in low or 'expert panel' in low:
                return 3
            # heuristics for 2 stars: criteria_provided + multiple_submitters + no_conflict
            if 'criteria_provided' in low and ('multiple_submitters' in low or 'multiple submitters' in low):
                if 'no_conflict' in low or 'no conflict' in low or 'no_conflicts' in low:
                    return 2
                return 1
            if 'criteria_provided' in low or 'criteria provided' in low:
                return 1
            return 0

        def safe_first_nonempty(row, candidates):
            """Return first non-empty value from row among candidate column names."""
            for cand in candidates:
                if not cand:
                    continue
                # pandas Series supports get(); csv.DictReader gives dict
                try:
                    if cand in row and row.get(cand) not in (None, "", ".", "NA", "None"):
                        return str(row.get(cand)).strip()
                except Exception:
                    try:
                        val = row.get(cand)
                        if val not in (None, "", ".", "NA", "None"):
                            return str(val).strip()
                    except Exception:
                        continue
            return ""

        # Try pandas first, but with robust candidate detection
        try:
            if HAVE_PANDAS:
                df = pd.read_csv(path, sep=None, engine='python', dtype=str, on_bad_lines='skip').fillna("")
                cols_lower = {c.lower().strip(): c for c in df.columns}
                # heuristics for columns
                def find_col(candidates):
                    for c in candidates:
                        if c and c.lower() in cols_lower:
                            return cols_lower[c.lower()]
                    # fuzzy substring
                    for k in cols_lower:
                        for cand in candidates:
                            if cand and cand in k:
                                return cols_lower[k]
                    return None

                qcol = find_col(['query_variant', 'queryvariant', 'qvar', 'query_var', 'Query_Variant'])
                vidcol = find_col(['variationid', 'variation_id', 'variation id', 'variationid', 'VariationID'])
                phencol = find_col(['reportedphenotypeinfo', 'reported_phenotype_info', 'phenotype', 'reportedphenotype', 'ReportedPhenotypeInfo'])
                # significance candidates (order of preference)
                sig_candidates = [
                    'ClinicalSignificance',      
                    'clinical_significance',
                    'clinicalsignificance',
                    'clinsig',
                    'clnsig',
                    'VCF_CLINICAL_SIGNIFICANCE',
                    'VCF_CLNSIG'                
                ]                
                # review/status and other columns
                revcol = find_col([
                    'ReviewStatus',              
                    'reviewstatus',
                    'review_stat',
                    'review',
                    'vcf_clnrevstat',            
                    'VCF_CLNREVSTAT'
                ])
                sidcol = find_col(['submissionid', 'submission_id', 'submission', 'SubmissionID'])
                hgvscol = find_col(['hgvs_protein', 'hgvs_prot', 'hgvsp', 'HGVS_Protein', 'protein'])
                genecol = find_col(['genesymbol', 'gene', 'gene_symbol', 'symbol', 'GeneSymbol'])

                for _, r in df.iterrows():
                    try:
                        # robust qvar extraction
                        qv = ""
                        if qcol:
                            qv = str(r.get(qcol, "") or "").strip()
                        if not qv:
                            # check several common keys
                            for cand in ('Query_Variant','query_variant','QueryVariant','qvar'):
                                if cand in r and r.get(cand):
                                    qv = str(r.get(cand)).strip(); break

                        vid = None
                        try:
                            vid_raw = r.get(vidcol, "") if vidcol else ""
                            vid = int(str(vid_raw).strip()) if str(vid_raw).strip() else None
                        except Exception:
                            vid = None

                        phen = ""
                        if phencol:
                            phen = str(r.get(phencol, "") or "").strip()
                        # clinical significance: try prioritized candidates, then scan keys
                        clnsig = safe_first_nonempty(r, sig_candidates)
                        if not clnsig:
                            # fallback: scan all columns for 'clnsig'/'significance'
                            for k in r.keys():
                                kl = str(k).lower()
                                if 'clnsig' in kl or 'significance' in kl:
                                    val = r.get(k)
                                    if val not in (None, "", ".", "NA", "None"):
                                        clnsig = str(val).strip()
                                        break

                        # Review/status
                        review = ""
                        if revcol:
                            review = str(r.get(revcol, "") or "").strip()
                        if not review:
                            # scan for likely review columns
                            for cand in ('VCF_CLNREVSTAT','ReviewStatus','reviewstatus','review_stat','review'):
                                if cand in r and r.get(cand):
                                    review = str(r.get(cand)).strip(); break

                        sid = str(r.get(sidcol, "") or "").strip() if sidcol else ""
                        hgvsp = ""
                        if hgvscol:
                            hgvsp = str(r.get(hgvscol, "") or "").strip()
                        gene = ""
                        if genecol:
                            gene = str(r.get(genecol, "") or "").strip()

                        stars = normalize_review_to_stars(review)

                        entry = {
                            "phenotype": phen,
                            "clnsig": clnsig,
                            "review": review,
                            "submission_id": sid,
                            "hgvsp": hgvsp,
                            "gene": gene,
                            "stars": stars
                        }

                        if vid is not None:
                            self.clinvar_expanded_by_vid[int(vid)].append(entry)
                        if qv:
                            self.clinvar_expanded_by_qvar[qv].append(entry)
                            # also index with/without chr- prefix
                            if qv.startswith("chr"):
                                alt = qv.replace("chr", "", 1)
                            else:
                                alt = "chr" + qv
                            if alt != qv:
                                self.clinvar_expanded_by_qvar[alt].append(entry)
                            # also index by gene-specific key
                            gene_norm = normalize_gene_name(gene) if gene else ""
                            if gene_norm:
                                key = f"{qv}|{gene_norm}"
                                self.clinvar_expanded_by_qvar[key].append(entry)
                                if alt != qv:
                                    self.clinvar_expanded_by_qvar[f"{alt}|{gene_norm}"].append(entry)
                    except Exception:
                        continue

                # Deduplicate stored lists (by vid/qvar) preserving order
                def dedupe_list(lst):
                    seen = set(); out = []
                    for s in lst:
                        try:
                            key = (str(s.get('submission_id','')), str(s.get('clnsig','')), str(s.get('phenotype','')), str(s.get('hgvsp','')))
                            if key in seen:
                                continue
                            seen.add(key)
                            out.append(s)
                        except Exception:
                            continue
                    return out

                for k in list(self.clinvar_expanded_by_vid.keys()):
                    self.clinvar_expanded_by_vid[k] = dedupe_list(self.clinvar_expanded_by_vid[k])
                for k in list(self.clinvar_expanded_by_qvar.keys()):
                    self.clinvar_expanded_by_qvar[k] = dedupe_list(self.clinvar_expanded_by_qvar[k])

                logger.info(f"✓ Loaded ClinVar expanded CSV with pandas: {len(self.clinvar_expanded_by_vid)} by VariationID, {len(self.clinvar_expanded_by_qvar)} by Query_Variant")
                return

        except Exception as e:
            logger.warning(f"Pandas parsing failed, trying fallback parser: {e}")

        # Fallback parser (csv.DictReader with semicolon or tab) - robust fallback
        try:
            import csv
            # try semicolon first, then tab/auto
            tried = False
            for delim in (';', '\t', ','):
                try:
                    with open(path, "r", encoding='utf-8', errors='replace') as fh:
                        reader = csv.DictReader(fh, delimiter=delim)
                        # if only one field name found with this delimiter, try next delimiter
                        if reader.fieldnames is None or len(reader.fieldnames) <= 1:
                            continue
                        for r in reader:
                            try:
                                # robust qv
                                qv = str(r.get('Query_Variant', "") or r.get('query_variant', "") or r.get('QueryVariant', "") or r.get(next(iter(r.keys())), "")).strip()
                                vid = None
                                try:
                                    vid_raw = r.get('VariationID') or r.get('variationid') or r.get('Variation_Id') or ""
                                    vid = int(str(vid_raw).strip()) if str(vid_raw).strip() else None
                                except Exception:
                                    vid = None
                                phen = str(r.get('ReportedPhenotypeInfo', "") or r.get('reportedphenotypeinfo', "") or r.get('phenotype', "")).strip()

                                # significance detection (robust)
                                clnsig = ""
                                for cand in ('VCF_CLNSIG', 'VCF_CLINICAL_SIGNIFICANCE', 'ClinicalSignificance', 'clinsig', 'clnsig', 'clinical_significance', 'clinicalsignificance'):
                                    if cand in r and r.get(cand) not in (None, "", ".", "NA", "None"):
                                        clnsig = str(r.get(cand)).strip()
                                        break
                                if not clnsig:
                                    # search keys for 'clnsig' or 'significance'
                                    for k in r.keys():
                                        if 'clnsig' in k.lower() or 'significance' in k.lower():
                                            val = r.get(k)
                                            clnsig = str(val).strip() if val not in (None, "") else ""
                                            if clnsig:
                                                break

                                review = ""
                                for cand in ('VCF_CLNREVSTAT', 'ReviewStatus', 'reviewstatus', 'review_stat', 'review'):
                                    if cand in r and r.get(cand):
                                        review = str(r.get(cand)).strip(); break
                                sid = str(r.get('SubmissionID', "") or r.get('submissionid', "") or "").strip()
                                hgvsp = str(r.get('HGVS_Protein', "") or r.get('hgvs_protein', "") or r.get('HGVSp', "") or "").strip()
                                gene = str(r.get('GeneSymbol', "") or r.get('genesymbol', "") or r.get('Gene', "") or "").strip()

                                stars = normalize_review_to_stars(review)

                                entry = {
                                    "phenotype": phen,
                                    "clnsig": clnsig,
                                    "review": review,
                                    "submission_id": sid,
                                    "hgvsp": hgvsp,
                                    "gene": gene,
                                    "stars": stars
                                }
                                if vid is not None:
                                    self.clinvar_expanded_by_vid[int(vid)].append(entry)
                                if qv:
                                    self.clinvar_expanded_by_qvar[qv].append(entry)
                                    if qv.startswith("chr"):
                                        alt_qv = qv.replace("chr","",1)
                                    else:
                                        alt_qv = "chr"+qv
                                    if alt_qv != qv:
                                        self.clinvar_expanded_by_qvar[alt_qv].append(entry)
                                    gene_norm = normalize_gene_name(gene) if gene else ""
                                    if gene_norm:
                                        key = f"{qv}|{gene_norm}"
                                        self.clinvar_expanded_by_qvar[key].append(entry)
                                        if alt_qv != qv:
                                            self.clinvar_expanded_by_qvar[f"{alt_qv}|{gene_norm}"].append(entry)
                            except Exception:
                                continue
                    # If we reached here without exception for this delimiter, break
                    tried = True
                    break
                except Exception:
                    continue

            # Deduplicate like above
            def dedupe_list(lst):
                seen = set(); out = []
                for s in lst:
                    try:
                        key = (str(s.get('submission_id','')), str(s.get('clnsig','')), str(s.get('phenotype','')), str(s.get('hgvsp','')))
                        if key in seen:
                            continue
                        seen.add(key)
                        out.append(s)
                    except Exception:
                        continue
                return out

            for k in list(self.clinvar_expanded_by_vid.keys()):
                self.clinvar_expanded_by_vid[k] = dedupe_list(self.clinvar_expanded_by_vid[k])
            for k in list(self.clinvar_expanded_by_qvar.keys()):
                self.clinvar_expanded_by_qvar[k] = dedupe_list(self.clinvar_expanded_by_qvar[k])

            logger.info(f"Loaded ClinVar expanded CSV (fallback parser): vids={len(self.clinvar_expanded_by_vid)} qvars={len(self.clinvar_expanded_by_qvar)}")
            return
        except Exception as e:
            logger.warning(f"Failed to load clinvar_expanded CSV with fallback parser: {e}", exc_info=True)
            return

    def _revstat_to_stars(self, revstat: str) -> int:
        """
        Convert ClinVar review status to star rating.
        
        Rules:
        - 4 stars: practice guideline
        - 3 stars: expert panel
        - 2 stars: criteria_provided, multiple_submitters, no_conflict
        - 1 star: criteria_provided (single submitter or with conflict)
        - 0 stars: no assertion criteria provided
        
        IMPORTANT: For conflicting classifications, we return 1 star (not 2)
        because the conflict prevents a "no_conflict" rating.
        """
        if not revstat:
            return 0
        low = str(revstat).lower()
        
        # Highest confidence: practice guideline
        if 'practice_guideline' in low or 'practice guideline' in low:
            return 4
        
        # Expert panel
        if 'expert_panel' in low or 'expert panel' in low:
            return 3
        
        # Check for conflict - this reduces confidence
        has_conflict = 'conflict' in low
        
        # criteria_provided with multiple submitters AND no conflict -> 2 stars
        if 'criteria_provided' in low and ('multiple_submitters' in low or 'multiple submitters' in low):
            if not has_conflict and ('no_conflict' in low or 'no conflict' in low):
                return 2
            # Has conflict -> return 1 star (downgraded)
            return 1
        
        # criteria_provided (single submitter or with conflict) -> 1 star
        if 'criteria_provided' in low or 'criteria provided' in low:
            return 1
        
        return 0

    def get_clinvar_vcf_aggregate(self, vr: VariantRecord) -> Tuple[Optional[str], int, str]:
        """
        Look up clinvar_vcf for variant and return (CLNSIG_normalized, stars, 'clinvar_vcf').
        
        CRITICAL: For conflicting classifications (CLNSIG=Conflicting_classifications_of_pathogenicity),
        this method returns the CONFLICT label, NOT the individual conflicting values.
        Stars are calculated based on CLNREVSTAT field.
        
        If no record found return (None,0,"")
        """
        # Get variant gene for filtering
        variant_gene_norm = normalize_gene_name(getattr(vr, "gene", "") or "")
        
        if not self.clinvar_vcf:
            # try to find clnsig in vr._raw_ann if generate.py included CLNSIG
            clnsig = (getattr(vr, "_raw_ann", {}) or {}).get("CLNSIG") or (getattr(vr, "clinvar_sig", "") or "")
            rev = (getattr(vr, "_raw_ann", {}) or {}).get("CLNREVSTAT") or ""
            if clnsig:
                cls = self._normalize_clinvar_label(clnsig)
                stars = self._revstat_to_stars(rev)
                return cls, int(stars), "clinvar_vcf_info"
            return None, 0, ""
        
        chrom_clean = vr.chrom
        chrom_candidates = [chrom_clean]
        if chrom_clean and chrom_clean.startswith("chr"):
            chrom_candidates.append(chrom_clean.replace("chr","",1))
        else:
            chrom_candidates.append("chr" + chrom_clean)
        
        for c in chrom_candidates:
            try:
                for rec in self.clinvar_vcf.fetch(c, vr.pos - 1, vr.pos):
                    if int(rec.pos) != int(vr.pos): 
                        continue
                    if rec.ref != vr.ref:
                        continue
                    if not (rec.alts and vr.alt in rec.alts):
                        continue
                    
                    # CHECK GENE MATCH before using VCF record
                    gene_info = rec.info.get('GENEINFO', '')
                    rec_genes = []
                    if gene_info:
                        for gene_pair in str(gene_info).split('|'):
                            if ':' in gene_pair:
                                gene_raw = gene_pair.split(':')[0]
                                rec_genes.append(normalize_gene_name(gene_raw))
                    
                    # If record has gene info and it doesn't match variant's gene, skip
                    if rec_genes and variant_gene_norm and variant_gene_norm not in rec_genes:
                        logger.debug(f"Skipping ClinVar VCF record for {vr.chrom}:{vr.pos}: gene mismatch "
                                f"(variant={variant_gene_norm}, record={rec_genes})")
                        continue
                    
                    info_rec = dict(rec.info)
                    
                    # Get CLNSIG - this is the AGGREGATE classification
                    clnsig = info_rec.get("CLNSIG")
                    if isinstance(clnsig, (list, tuple)):
                        clnsig_s = ";".join([str(x) for x in clnsig])
                    else:
                        clnsig_s = str(clnsig or "")
                    
                    # Get CLNREVSTAT for star rating
                    rev = info_rec.get("CLNREVSTAT") or ""
                    
                    # CRITICAL: Normalize the label correctly
                    # If CLNSIG contains "Conflicting", return the conflict label
                    # NOT the individual conflicting values from CLNSIGCONF
                    cls = self._normalize_clinvar_label(clnsig_s)
                    
                    # Calculate stars from review status
                    stars = self._revstat_to_stars(rev)
                    
                    # Debug logging for conflict detection
                    if "conflict" in clnsig_s.lower():
                        logger.debug(f"VCF conflict detected for {vr.chrom}:{vr.pos}: CLNSIG={clnsig_s}, stars={stars}")
                    
                    return cls, int(stars), "clinvar_vcf"
            except ValueError:
                continue
            except Exception:
                continue
        
        # last fallback: vr.clinvar_sig if set by previous logic
        existing_sig = getattr(vr, "clinvar_sig", "") or ""
        existing_rev = getattr(vr, "clinvar_review", "") or ""
        if existing_sig:
            cls = self._normalize_clinvar_label(existing_sig)
            stars = self._revstat_to_stars(existing_rev)
            return cls, int(stars), "variant_record_field"
        
        return None, 0, ""

    def get_clinvar_filtered_entries(self, vr: VariantRecord) -> List[Dict[str, Any]]:
        """Get all entries from filtered table for this variant."""
        if not hasattr(self, "clinvar_filtered_by_qvar"):
            # Try to load if not loaded
            filtered_path = self.paths.get("CLINVAR_FILTERED_CSV") or os.path.join("databases", "clinvar", "clinvar_conditions_summary.mane_filtered.csv")
            try:
                self._load_clinvar_filtered_csv(filtered_path)
            except Exception:
                return []
        
        variant_gene_norm = normalize_gene_name(getattr(vr, "gene", "") or "")
    
        entries = []
        try:
            qvar = f"{vr.chrom}-{vr.pos}-{vr.ref}-{vr.alt}"
            # Try with and without chr prefix
            for key in [qvar, qvar.replace("chr", ""), "chr" + qvar if not qvar.startswith("chr") else qvar[3:]]:
                for entry in self.clinvar_filtered_by_qvar.get(key, []):
                    # Filter by gene to prevent cross-gene contamination
                    entry_gene_norm = normalize_gene_name(entry.get("GeneSymbol", "") or "")
                    # If entry has a gene and it doesn't match variant's gene, skip
                    if entry_gene_norm and variant_gene_norm and entry_gene_norm != variant_gene_norm:
                        logger.debug(f"Skipping ClinVar entry for {vr.chrom}:{vr.pos}: gene mismatch "
                                f"({entry_gene_norm} vs {variant_gene_norm})")
                        continue
                    
                    entries.append(entry)
        except Exception:
            pass
        
        # Deduplicate
        seen = set()
        unique = []
        for e in entries:
            key = (e.get("VariationID"), e.get("ReportedPhenotypeInfo"))
            if key not in seen:
                seen.add(key)
                unique.append(e)
        return unique

    def get_clinvar_classification_with_fallback(self, vr: VariantRecord, acmg_disease: Optional[str] = None) -> Tuple[Optional[str], int, str, bool]:
        """
        Get ClinVar classification with proper phenotype-aware fallback.
        
        CRITICAL LOGIC (MUST FOLLOW EXACTLY):
        1. If acmg_disease is provided, search filtered_entries for entries whose phenotype MATCHES acmg_disease
        2. If ANY matching entries exist -> use consensus from those entries (high confidence)
        3. If NO matching entries exist -> IGNORE filtered_entries completely (do NOT use non-matching phenotypes!)
        4. Fallback to VCF aggregate (which gives true consensus across all submissions)
        
        This prevents using classifications from unrelated phenotypes (e.g., using "Likely benign" for 
        "Multiple self-healing squamous epithelioma" when ACMG disease is "Loeys-Dietz syndrome").
        """
        filtered_entries = self.get_clinvar_filtered_entries(vr)
    
        # ============================================================
        # FIX: Get variant gene for filtering in case entries slipped through
        # ============================================================
        variant_gene_norm = normalize_gene_name(getattr(vr, "gene", "") or "")
        
        # Case 1: We have ACMG disease to match against
        if acmg_disease and filtered_entries:
            matched_entries = []
            non_matched_phenotypes = []
            
            for entry in filtered_entries:
                entry_gene_norm = normalize_gene_name(entry.get("GeneSymbol", "") or "")
                if entry_gene_norm and variant_gene_norm and entry_gene_norm != variant_gene_norm:
                    logger.debug(f"Skipping entry in fallback: gene mismatch ({entry_gene_norm} vs {variant_gene_norm})")
                    continue
                
                phen = entry.get("ReportedPhenotypeInfo", "") or ""
                if phen:
                    matched, level, reason = check_disease_match(acmg_disease, phen)
                    if matched:
                        matched_entries.append(entry)
                        logger.debug(f"Phenotype match for {vr.chrom}:{vr.pos}: {acmg_disease} -> {phen}")
                    else:
                        non_matched_phenotypes.append(phen)
            
            # ONLY use filtered entries if we have phenotype matches
            if matched_entries:
                consensus_sig = self._consensus_classification_from_rows(matched_entries)
                stars = self._aggr_stars_for_rows(matched_entries)
                logger.debug(f"Using matched filtered entries for {vr.chrom}:{vr.pos}: sig={consensus_sig}, stars={stars}")
                return consensus_sig, stars, "clinvar_filtered_phenotype_match", True
            else:
                # No phenotype match - do NOT use filtered entries!
                logger.debug(f"No phenotype match for {vr.chrom}:{vr.pos}. ACMG disease={acmg_disease}, available phenotypes={non_matched_phenotypes[:3]}. Falling back to VCF.")
                # Continue to VCF fallback - do NOT return filtered entries
        
        # Case 2: No ACMG disease provided OR no phenotype matches found
        # CRITICAL: Do NOT use filtered_entries without phenotype match!
        # Fall back to VCF aggregate which gives the true consensus across all submissions
        cls, stars, src = self.get_clinvar_vcf_aggregate(vr)
        
        # Debug: log what VCF returned
        if cls:
            logger.debug(f"VCF fallback for {vr.chrom}:{vr.pos}: cls={cls}, stars={stars}, src={src}")
        
        # If VCF returned something meaningful, use it
        if cls and cls != "Not provided":
            return cls, stars, src or "clinvar_vcf", False
        
        # Ultimate fallback: if VCF also empty, return None (no evidence)
        return None, 0, "", False

    # -------------------------
    # ClinVar filtered table (from clinvar_merge.py)
    # -------------------------
    def _load_clinvar_filtered_csv(self, path: Optional[str]) -> None:
        """
        Load table produced by clinvar_merge.py / clinvar_conditions_summary.mane_filtered.csv
        Expected columns (semicolon-separated): Query_Variant;GeneSymbol;VariationID;ReportedPhenotypeInfo;Classifications;Stars;Evidence_Count
        Builds indices:
          - self.clinvar_filtered_by_qvar: { "chr1-...-A-T": [entry, ...] }
          - self.clinvar_filtered_by_vid: { 12345: [entry, ...] }
        Each entry is a dict with parsed fields.
        """
        self.clinvar_filtered_by_qvar = defaultdict(list)
        self.clinvar_filtered_by_vid = defaultdict(list)
        if not path or not os.path.exists(path):
            logger.info("ClinVar filtered CSV not found: %s", path)
            return
        logger.info("Loading ClinVar filtered CSV: %s", path)
        try:
            # prefer pandas if available (fast + tolerant)
            if HAVE_PANDAS:
                df = pd.read_csv(path, sep=';', dtype=str, encoding='utf-8', on_bad_lines='skip').fillna("")
                for _, r in df.iterrows():
                    qv = str(r.get('Query_Variant', "") or "").strip()
                    vid_raw = r.get('VariationID', "") or ""
                    try:
                        vid = int(str(vid_raw).strip()) if str(vid_raw).strip() else None
                    except Exception:
                        vid = None
                    entry = {
                        "Query_Variant": qv,
                        "GeneSymbol": str(r.get('GeneSymbol', "") or ""),
                        "VariationID": vid,
                        "ReportedPhenotypeInfo": str(r.get('ReportedPhenotypeInfo', "") or ""),
                        "Classifications": str(r.get('Classifications', "") or ""),
                        "Stars": int(r.get('Stars', 0) or 0),
                        "Evidence_Count": int(r.get('Evidence_Count', 0) or 0)
                    }
                    if qv:
                        self.clinvar_filtered_by_qvar[qv].append(entry)
                        # also store alt key without/with 'chr' prefix both ways
                        if qv.startswith("chr"):
                            alt = qv.replace("chr", "", 1)
                        else:
                            alt = "chr" + qv
                        if alt != qv:
                            self.clinvar_filtered_by_qvar[alt].append(entry)
                    if vid:
                        self.clinvar_filtered_by_vid[int(vid)].append(entry)
                logger.info("Loaded clinvar_filtered entries: qvars=%d vids=%d", len(self.clinvar_filtered_by_qvar), len(self.clinvar_filtered_by_vid))
                return
            # fallback simple csv parser
            import csv
            with open(path, "rt", encoding="utf-8", errors="replace") as fh:
                reader = csv.DictReader(fh, delimiter=';')
                for r in reader:
                    qv = str(r.get('Query_Variant', "") or "").strip()
                    vid_raw = r.get('VariationID', "") or ""
                    try:
                        vid = int(str(vid_raw).strip()) if str(vid_raw).strip() else None
                    except Exception:
                        vid = None
                    entry = {
                        "Query_Variant": qv,
                        "GeneSymbol": str(r.get('GeneSymbol', "") or ""),
                        "VariationID": vid,
                        "ReportedPhenotypeInfo": str(r.get('ReportedPhenotypeInfo', "") or ""),
                        "Classifications": str(r.get('Classifications', "") or ""),
                        "Stars": int(r.get('Stars', 0) or 0),
                        "Evidence_Count": int(r.get('Evidence_Count', 0) or 0)
                    }
                    if qv:
                        self.clinvar_filtered_by_qvar[qv].append(entry)
                        if qv.startswith("chr"):
                            alt = qv.replace("chr","",1)
                        else:
                            alt = "chr" + qv
                        if alt != qv:
                            self.clinvar_filtered_by_qvar[alt].append(entry)
                    if vid:
                        self.clinvar_filtered_by_vid[int(vid)].append(entry)
            logger.info("Loaded clinvar_filtered (fallback) entries: qvars=%d vids=%d", len(self.clinvar_filtered_by_qvar), len(self.clinvar_filtered_by_vid))
        except Exception as e:
            logger.warning("Failed to load clinvar_filtered CSV %s: %s", path, e, exc_info=True)
            self.clinvar_filtered_by_qvar = defaultdict(list)
            self.clinvar_filtered_by_vid = defaultdict(list)

    def _normalize_clinvar_label(self, s: str) -> str:
        if not s:
            return ""
        sl = str(s).lower()
        sl = sl.replace("_", " ")
        sl = sl.strip()
        
        # CRITICAL: Check conflict FIRST (before pathogenic)
        if "conflict" in sl:
            return "Conflicting classifications of pathogenicity"
        
        if "likely pathogenic" in sl:
            return "Likely pathogenic"
        
        if "pathogenic" in sl and "likely" not in sl:
            return "Pathogenic"
        
        if "likely benign" in sl:
            return "Likely benign"
        
        if "benign" in sl and "pathogen" not in sl:
            return "Benign"
        
        if "uncertain" in sl or "vus" in sl:
            return "Uncertain significance"
        
        return s

    def _aggr_stars_for_rows(self, rows: List[Dict[str,Any]]) -> int:
        """
        Compute aggregate stars following your clinvar_merge logic:
         - if any row has stars>=3 -> return max star
         - else if >=2 independent 1-star submissions -> promote to 2
         - else return max(stars)
        """
        if not rows:
            return 0
        stars_list = [int(r.get("Stars", 0) or 0) for r in rows]
        max_star = max(stars_list) if stars_list else 0
        if max_star >= 3:
            return max_star
        # count how many 1-star rows (evidence_count>0) - treat them as independent
        ones = sum(1 for s in stars_list if s == 1)
        if ones >= 2:
            return 2
        return max_star

    def _consensus_classification_from_rows(self, rows: List[Dict[str,Any]]) -> str:
        """
        Take several rows (each contains 'Classifications' computed by clinvar_merge) and apply your consensus rules:
        - Filter by authority (we assume rows are already per-phenotype aggregated)
        - Build vector presence and decide Conflict vs Pathogenic vs Benign vs VUS.
        This mirrors get_consensus_classification() in clinvar_merge.py
        """
        if not rows:
            return ""
        labels = [str(r.get("Classifications","") or "").strip() for r in rows if r.get("Classifications")]
        labels_low = [l.lower() for l in labels]
        has_path = any('pathogenic' in l and 'benign' not in l for l in labels_low)
        has_benign = any('benign' in l for l in labels_low)
        has_vus = any('uncertain' in l or 'vus' in l for l in labels_low)
        categories = sum([has_path, has_benign, has_vus])
        if categories > 1:
            return "Conflicting classifications of pathogenicity"
        if has_path:
            # if at least one strict Pathogenic (not likely) present -> Pathogenic; else Likely pathogenic
            if any('pathogenic' in l and 'likely' not in l for l in labels_low):
                return "Pathogenic"
            return "Likely pathogenic"
        if has_benign:
            if any('benign' in l and 'likely' not in l for l in labels_low):
                return "Benign"
            return "Likely benign"
        if has_vus:
            return "Uncertain significance"
        # fallback: join uniques
        uniq = sorted(set(labels))
        return "/".join(uniq)

    def get_clinvar_filtered_classification(self, vr: VariantRecord, acmg_disease: Optional[str] = None) -> Tuple[Optional[str], int, str]:
        """
        Get classification from filtered table ONLY if phenotype matches ACMG disease.
        
        IMPORTANT: This method should ONLY return a classification when there is a 
        direct phenotype match between the ACMG disease and the ReportedPhenotypeInfo.
        
        Returns (classification, stars, source_description) or (None, 0, "") if:
        - No filtered entries exist, OR
        - No phenotype match found (to prevent using wrong phenotype classifications)
        
        This method should NOT be used as a fallback without phenotype matching.
        """
        # ensure index loaded
        if not hasattr(self, "clinvar_filtered_by_qvar"):
            return None, 0, ""
        
        # build qvar keys to try
        try:
            qvar = None
            if getattr(vr, "chrom", None) and getattr(vr, "pos", None) and getattr(vr, "ref", None) is not None and getattr(vr, "alt", None) is not None:
                qvar = f"{vr.chrom}-{int(vr.pos)}-{vr.ref}-{vr.alt}"
        except Exception:
            qvar = None
        
        candidates = []
        if qvar:
            candidates.extend(self.clinvar_filtered_by_qvar.get(qvar, []) or [])
            # try alt form without/with chr prefix
            if qvar.startswith("chr"):
                alt = qvar.replace("chr","",1)
            else:
                alt = "chr" + qvar
            if alt != qvar:
                candidates.extend(self.clinvar_filtered_by_qvar.get(alt, []) or [])
        
        # by VariationID also
        vid = None
        try:
            vid = int(vr._raw_ann.get("VARIATIONID") or vr._raw_ann.get("VAR_ID") or vr._raw_ann.get("VARID") or vr._raw_ann.get("ALLELEID") or 0)
            if not vid:
                vid = None
        except Exception:
            vid = None
        
        if vid:
            candidates.extend(self.clinvar_filtered_by_vid.get(int(vid), []) or [])
        
        if not candidates:
            return None, 0, ""
        
        # CRITICAL: Only return classification if we have phenotype match
        if acmg_disease:
            matched_entries = []
            for r in candidates:
                phen = r.get("ReportedPhenotypeInfo", "") or ""
                if not phen:
                    continue
                try:
                    ok, level, reason = check_disease_match(acmg_disease, phen)
                    if ok:
                        matched_entries.append(r)
                except Exception:
                    continue
            
            # Only return if we found matching entries
            if matched_entries:
                cls = self._consensus_classification_from_rows(matched_entries)
                stars = self._aggr_stars_for_rows(matched_entries)
                return cls, int(stars), "clinvar_filtered_per_disease"
        
        # No phenotype match found - return None (caller should fall back to VCF)
        return None, 0, ""

    def _get_expanded_submissions_for_variant(self, variation_id: Optional[int], qvar: Optional[str]) -> List[Dict]:
        """
        Return deduplicated, ordered list of expanded ClinVar submission dicts for a variant.
        Combines data from:
        - self.clinvar_expanded_by_vid (VariationID -> [submissions])
        - self.clinvar_expanded_by_qvar (Query_Variant -> [submissions])
        Handles 'chr' prefix variance and deduplicates by (submission_id, clnsig, phenotype, hgvsp).
        """
        out = []
        try:
            # Gather by VariationID if available
            if variation_id is not None and hasattr(self, "clinvar_expanded_by_vid"):
                try:
                    out.extend(self.clinvar_expanded_by_vid.get(int(variation_id), []) or [])
                except Exception:
                    # defensive: ignore malformed vid
                    pass

            # Gather by Query_Variant if available
            if qvar and hasattr(self, "clinvar_expanded_by_qvar"):
                try:
                    # direct key
                    out.extend(self.clinvar_expanded_by_qvar.get(qvar, []) or [])
                    # try w/o/with chr prefix
                    if qvar.startswith("chr"):
                        alt = qvar.replace("chr", "", 1)
                    else:
                        alt = "chr" + qvar
                    if alt != qvar:
                        out.extend(self.clinvar_expanded_by_qvar.get(alt, []) or [])
                except Exception:
                    pass

            # If we still have no results but an alternative index exists (backwards compatibility),
            # try to use any clinvar_submissions_by_qvar-like attribute if present.
            # This keeps behavior robust if different loader sets a different attribute name.
            if not out:
                for attr in ("clinvar_submissions_by_qvar", "clinvar_expanded_by_qvar"):
                    if hasattr(self, attr):
                        try:
                            idx = getattr(self, attr)
                            if qvar in idx:
                                out.extend(idx.get(qvar, []) or [])
                            if qvar and qvar.startswith("chr") and qvar.replace("chr","") in idx:
                                out.extend(idx.get(qvar.replace("chr",""), []) or [])
                            if qvar and not qvar.startswith("chr") and ("chr" + qvar) in idx:
                                out.extend(idx.get("chr" + qvar, []) or [])
                        except Exception:
                            continue

            # Deduplicate preserving order
            seen = set()
            uniq = []
            for s in out:
                sid = str(s.get("submission_id", "") or "")
                sig = str(s.get("clnsig", "") or "")
                phen = str(s.get("phenotype", "") or "")
                hgvsp = str(s.get("hgvsp", "") or "")
                key = (sid, sig, phen, hgvsp)
                if key in seen:
                    continue
                seen.add(key)
                uniq.append(s)
            return uniq
        except Exception:
            # On unexpected errors, return empty list but do not raise to keep pipeline running.
            return []

    def get_submissions_for_query_variant(self, qvar: str, variant_gene: Optional[str] = None) -> List[Dict[str,str]]:
        """
        Return list of submission dicts for a query-variant string (robust against 'chr' prefix).
        If variant_gene is provided, prefer submissions annotated to that gene (indexed by "qvar|GENE").
        Falls back to non-gene-specific submissions if no gene-specific entries found.
        """
        out: List[Dict[str,str]] = []
        if not qvar:
            return out

        # 1) Try gene-specific key if variant_gene provided
        try:
            if variant_gene:
                gene_norm = normalize_gene_name(variant_gene)
                if gene_norm:
                    key = f"{qvar}|{gene_norm}"
                    if key in getattr(self, "clinvar_expanded_by_qvar", {}):
                        out.extend(self.clinvar_expanded_by_qvar.get(key, []) or [])
                        return out
                    # try alt with/without chr prefix
                    if qvar.startswith("chr"):
                        alt = qvar.replace("chr", "", 1)
                    else:
                        alt = f"chr{qvar}"
                    key_alt = f"{alt}|{gene_norm}"
                    if key_alt in getattr(self, "clinvar_expanded_by_qvar", {}):
                        out.extend(self.clinvar_expanded_by_qvar.get(key_alt, []) or [])
                        return out
        except Exception:
            # on failure, fall back to positional lookup below
            pass

        # 2) Fallback: positional queries (existing behavior)
        if hasattr(self, "clinvar_expanded_by_qvar"):
            try:
                out.extend(self.clinvar_expanded_by_qvar.get(qvar, []) or [])
                # try with/without chr prefix
                if qvar.startswith("chr"):
                    alt_qv = qvar.replace("chr", "", 1)
                else:
                    alt_qv = f"chr{qvar}"
                if alt_qv != qvar:
                    out.extend(self.clinvar_expanded_by_qvar.get(alt_qv, []) or [])
            except Exception:
                pass

        return out


    def _load_internal_db(self, path: str):
        if not os.path.exists(path):
            logger.warning(f"Internal DB path not found: {path}")
            return
        try:
            logger.info(f"Loading Internal DB from {path}...")
            df = pd.read_csv(path, dtype=str).fillna("")
            if 'vid' not in df.columns or 'anntype' not in df.columns:
                logger.warning("Internal DB CSV missing 'vid' or 'anntype' columns")
                return
            cls_map = {
                "P": "P", "Ps": "P", "Pc": "P",
                "LP": "LP", "LPs": "LP", "LPc": "LP",
                "B": "B",
                "LB": "LB",
                "VUS": "VUS"
            }
            count = 0
            for _, r in df.iterrows():
                vid = str(r.get('vid', '')).strip()
                if not vid: continue
                parts = vid.split(':')
                if len(parts) < 3: continue
                c_raw = parts[0]; p_raw = parts[1]; ref = parts[2]; alt = parts[3]
                raw_type = str(r.get('anntype', '')).strip()
                mapped_cls = cls_map.get(raw_type)
                if not mapped_cls:
                    continue
                chrom = normalize_chrom(c_raw)
                try: pos = int(p_raw)
                except ValueError: continue
                alts = [x.strip() for x in alt.split(',')]
                for a in alts:
                    if not a: continue
                    self.internal_db_index[chrom][pos][ref][a] = {
                        "cls": mapped_cls,
                        "date": str(r.get('date', ''))
                    }
                    count += 1
            logger.info(f"Internal DB loaded successfully. Indexed {count} entries.")
        except Exception as e:
            logger.error(f"Failed to load internal DB: {e}", exc_info=True)

    def _load_hgmd_csv(self) -> None:
        path = self.paths.get("HGMD_VCF")
        if not path or not os.path.exists(path):
            path = "/home/nik/share/ccu-ngs/ngs/bases/hgmd/dmsupport/hgmd_2025_vars_dmsupport.csv"
        if not os.path.exists(path):
            logger.warning("HGMD CSV not found: %s", path)
            return

        try:
            logger.info("Loading HGMD from CSV: %s", path)
            # Read with pandas if available 
            if HAVE_PANDAS:
                df = pd.read_csv(path, dtype=str, on_bad_lines='skip').fillna("")
                # Build lower->orig column map
                cols_map = {c.lower().strip(): c for c in df.columns}
                def get_col_val(row, candidates, default=""):
                    for cand in candidates:
                        if cand in cols_map:
                            return str(row.get(cols_map[cand], "")).strip()
                    return default

                # candidate names
                chrom_cands = ["chrom", "chr", "chromosome"]
                pos_cands = ["pos", "position", "start"]
                ref_cands = ["ref", "reference", "reference_allele"]
                alt_cands = ["alt", "alternate", "alt_allele"]
                gene_cands = ["gene", "symbol", "gene_symbol"]
                tag_cands = ["tag", "class", "hgmd_class"]
                pubs_cands = ["publications", "pubmed", "pubs", "npublications", "publication_count"]
                dmsup_cands = ["dmsupported", "dm_supported", "dm support", "dmsupport", "dm_supported_score", "dm_support"]

                count = 0
                for _, r in df.iterrows():
                    chrom_raw = get_col_val(r, chrom_cands, "")
                    if not chrom_raw:
                        continue
                    chrom = normalize_chrom(chrom_raw)
                    pos_raw = get_col_val(r, pos_cands, "")
                    try:
                        pos = int(float(pos_raw)) if pos_raw not in ("", "NA", "NULL") else 0
                    except Exception:
                        continue
                    ref = get_col_val(r, ref_cands, "")
                    alt = get_col_val(r, alt_cands, "")
                    if not ref or not alt:
                        # skip incomplete
                        continue

                    hgmd_gene = get_col_val(r, gene_cands, "")
                    tag = get_col_val(r, tag_cands, "")
                    pubs_val = get_col_val(r, pubs_cands, "")
                    try:
                        pub_count = int(float(pubs_val)) if pubs_val not in ("", "NA", "NULL") else 0
                    except Exception:
                        # try to extract number from text
                        nums = re.findall(r'\d+', pubs_val or "")
                        pub_count = int(nums[0]) if nums else 0

                    # dmsupported parsing
                    dmsupported_val = get_col_val(r, dmsup_cands, "")
                    dmsupported = 0
                    try:
                        if dmsupported_val and str(dmsupported_val).strip() not in ("", "NA", "NULL"):
                            dmsupported = int(float(re.sub(r'[^\d\.]', '', str(dmsupported_val))))
                    except Exception:
                        dmsupported = 0

                    rank = get_col_val(r, ["rankscore", "rank", "score"], "")
                    phen = get_col_val(r, ["disease", "phenotype", "phen"], "")

                    self.hgmd_index[chrom][pos][ref][alt] = {
                        "class": tag,
                        "phen": phen,
                        "pubs": pub_count,
                        "dmsupported": dmsupported,
                        "rankscore": rank,
                        "gene": hgmd_gene,
                    }
                    count += 1

                logger.info("Loaded HGMD index with %d variants", count)
                return
            else:
                # fallback simple parser if pandas not available
                with open(path, "rt", encoding="utf-8", errors="replace") as fh:
                    reader = csv.DictReader(fh)
                    count = 0
                    for r in reader:
                        chrom_raw = r.get("CHROM") or r.get("chrom") or r.get("chr")
                        if not chrom_raw:
                            continue
                        chrom = normalize_chrom(chrom_raw)
                        try:
                            pos = int(float(r.get("POS") or r.get("pos") or 0))
                        except Exception:
                            continue
                        ref = str(r.get("REF","")).strip()
                        alt = str(r.get("ALT","")).strip()
                        if not ref or not alt:
                            continue
                        hgmd_gene = r.get("GENE") or r.get("SYMBOL") or ""
                        tag = r.get("tag") or r.get("class") or ""
                        pubs_val = r.get("publications") or r.get("pubmed") or r.get("pubs") or ""
                        try:
                            pub_count = int(float(pubs_val)) if pubs_val not in ("", None) else 0
                        except Exception:
                            nums = re.findall(r'\d+', str(pubs_val) or "")
                            pub_count = int(nums[0]) if nums else 0
                        dmsupported_val = r.get("dmsupported") or r.get("dmsupport") or r.get("dm_supported") or ""
                        try:
                            dmsupported = int(float(dmsupported_val)) if dmsupported_val not in ("", None) else 0
                        except Exception:
                            dmsupported = 0
                        rank = r.get("rankscore") or ""
                        phen = r.get("disease") or ""
                        self.hgmd_index[chrom][pos][ref][alt] = {
                            "class": tag,
                            "phen": phen,
                            "pubs": pub_count,
                            "dmsupported": dmsupported,
                            "rankscore": rank,
                            "gene": hgmd_gene,
                        }
                        count += 1
                logger.info("Loaded HGMD index (fallback) with %d variants", count)
                return
        except Exception as e:
            logger.warning("HGMD CSV load failed: %s", e)
            return

    def get_variant_db_evidence(self, vr: VariantRecord) -> List[Tuple[str,str]]:
        results: List[Tuple[str,str]] = []
        vr.has_strong_pathogenic_evidence = False
        vr.has_strong_benign_evidence = False
        vr.is_hq_pathogenic_db = False
        vr.is_hq_benign_db = False
        vr.is_hq_conflict_db = False
        vr.clinvar_submissions = []
        
        # Initialize manual_reasons if needed
        if not hasattr(vr, "manual_reasons") or vr.manual_reasons is None:
            vr.manual_reasons = []
        
        acmg_disease = getattr(vr, "_acmg_disease", None)
        
        filtered_entries = self.get_clinvar_filtered_entries(vr)
        
        # Convert filtered entries to submission format for variant_row_with_disease
        for entry in filtered_entries:
            entry_gene = normalize_gene_name(entry.get("GeneSymbol", ""))
            variant_gene_norm = normalize_gene_name(vr.gene)
            
            # Skip entries for other genes at same coordinate
            if entry_gene and variant_gene_norm and entry_gene != variant_gene_norm:
                logger.debug(f"Skipping ClinVar entry for {vr.chrom}:{vr.pos}: gene mismatch "
                            f"({entry_gene} vs {variant_gene_norm})")
                continue 
            
            submission = {
                "submission_id": str(entry.get("VariationID", "")),
                "clnsig": entry.get("Classifications", ""),
                "phenotype": entry.get("ReportedPhenotypeInfo", ""),
                "stars": entry.get("Stars", 0),
                "review": "",  # Not available in filtered table
                "hgvsp": "",
                "gene": entry.get("GeneSymbol", ""),
            }
            # Deduplicate before adding
            already_exists = False
            for existing in vr.clinvar_submissions:
                if (existing.get("submission_id") == submission["submission_id"] and 
                    existing.get("phenotype") == submission["phenotype"]):
                    already_exists = True
                    break
            if not already_exists and submission["clnsig"] and submission["clnsig"] not in ("", "Not provided"):
                vr.clinvar_submissions.append(submission)
        
        # =========================================================
        # Now get classification with fallback (for scoring)
        # CRITICAL: Only use filtered entries if they match ACMG disease phenotype
        # Otherwise fall back to VCF aggregate (which may show conflicts)
        # =========================================================
        clinvar_cls, clinvar_stars, clinvar_source, phen_match = self.get_clinvar_classification_with_fallback(
            vr, acmg_disease=acmg_disease
        )
        
        if clinvar_cls:
            results.append((clinvar_cls, f"ClinVar_{clinvar_source}_{clinvar_stars}stars"))
            vr.clinvar_sig = clinvar_cls
            vr.clinvar_stars = str(clinvar_stars)
            vr._clinvar_phenotype_match = phen_match
            
            # CRITICAL: Store trait from VCF if available and no phenotype match
            if not phen_match and clinvar_source == "clinvar_vcf":
                # Try to get CLNDN from VCF record
                if self.clinvar_vcf:
                    chrom_clean = vr.chrom
                    chrom_candidates = [chrom_clean]
                    if chrom_clean and chrom_clean.startswith("chr"):
                        chrom_candidates.append(chrom_clean.replace("chr", "", 1))
                    else:
                        chrom_candidates.append("chr" + chrom_clean)
                    
                    # Get variant gene for filtering
                    variant_gene_norm = normalize_gene_name(getattr(vr, "gene", "") or "")
                    
                    for c in chrom_candidates:
                        try:
                            for rec in self.clinvar_vcf.fetch(c, vr.pos - 1, vr.pos):
                                if int(rec.pos) != int(vr.pos):
                                    continue
                                if rec.ref != vr.ref:
                                    continue
                                if not (rec.alts and vr.alt in rec.alts):
                                    continue
                                
                                # CHECK GENE MATCH before using CLNDN
                                gene_info = rec.info.get('GENEINFO', '')
                                rec_genes = []
                                if gene_info:
                                    for gene_pair in str(gene_info).split('|'):
                                        if ':' in gene_pair:
                                            gene_raw = gene_pair.split(':')[0]
                                            rec_genes.append(normalize_gene_name(gene_raw))
                                
                                # If record has gene info and it doesn't match variant's gene, skip
                                if rec_genes and variant_gene_norm and variant_gene_norm not in rec_genes:
                                    logger.debug(f"Skipping CLNDN: gene mismatch for {vr.chrom}:{vr.pos}")
                                    continue
                                
                                info_rec = dict(rec.info)
                                clndn = info_rec.get("CLNDN")
                                if clndn:
                                    if isinstance(clndn, (list, tuple)):
                                        vr.clinvar_trait = ";".join(str(x) for x in clndn)
                                    else:
                                        vr.clinvar_trait = str(clndn)
                                break
                        except:
                            continue

        # INTERNAL DB
        internal_evidence = []
        try:
            entry = self.internal_db_index.get(vr.chrom, {}).get(vr.pos, {}).get(vr.ref, {}).get(vr.alt)
            if entry:
                int_cls = str(entry.get("cls", "")) if isinstance(entry, dict) else str(entry)
                if int_cls in ["P", "LP"]:
                    classification = "Pathogenic" if int_cls == "P" else "Likely pathogenic"
                    internal_evidence.append((classification, "Internal_DB"))
                    vr.internal_db_sig = classification
                elif int_cls in ["B", "LB"]:
                    classification = "Benign" if int_cls == "B" else "Likely benign"
                    internal_evidence.append((classification, "Internal_DB"))
                    vr.internal_db_sig = classification
                else:
                    vr.internal_db_sig = int_cls
        except Exception as e:
            logger.debug(f"Internal DB lookup failed: {e}")
        results.extend(internal_evidence)

        # HGMD 
        hgmd_evidence = []
        try:
            hgmd_entry = self.hgmd_index.get(vr.chrom, {}).get(vr.pos, {}).get(vr.ref, {}).get(vr.alt)
            if hgmd_entry:
                tag = str(hgmd_entry.get('class', '')).strip()
                pubs_str = str(hgmd_entry.get('pubs', '0')).strip()
                hgmd_gene = str(hgmd_entry.get('gene', '')).strip()
                rank = str(hgmd_entry.get('rankscore', '')).strip()
                dmsupported = int(hgmd_entry.get('dmsupported', 0) or 0)

                gene_valid, gene_reason = validate_gene_match_db_evidence(vr.gene, hgmd_gene, "HGMD")
                if not gene_valid:
                    # do NOT treat this HGMD entry as evidence for this variant's gene
                    vr.manual_reasons.append(f"Skipped HGMD entry at {vr.chrom}:{vr.pos} {hgmd_gene}: {gene_reason}")
                    # keep hgmd meta defaults (do not set hgmd_class/hgmd_phen/pubs as DM)
                    vr.hgmd_class = ""
                    vr.hgmd_phen = ""
                    vr.hgmd_publication_count = 0
                    vr.hgmd_dmsupported = 0
                else:
                    # existing handling for valid gene
                    try:
                        pub_count = int(pubs_str)
                    except (ValueError, TypeError):
                        pub_count = 0
                    vr.hgmd_publication_count = pub_count
                    vr.hgmd_rankscore = rank if rank and rank != "NULL" else "Not provided"
                    vr.hgmd_phen = str(hgmd_entry.get('phen', '')).strip()
                    vr.hgmd_class = tag
                    vr.hgmd_dmsupported = dmsupported
                    is_well_supported = dmsupported >= 2
                    if tag == 'DM' and is_well_supported:
                        vr.has_strong_pathogenic_evidence = True
                        vr.is_hq_pathogenic_db = True
                        hgmd_evidence.append(("Pathogenic", f"HGMD_DM_supported{dmsupported}"))
                    elif tag == 'DM':
                        vr.has_weak_pathogenic_evidence = True
                        hgmd_evidence.append(("Pathogenic", "HGMD_DM_weak_support"))
                    elif tag:
                        hgmd_evidence.append((tag, "HGMD"))
            else:
                # Initialize HGMD fields to defaults if not found
                vr.hgmd_dmsupported = 0
                vr.hgmd_publication_count = 0
                vr.hgmd_rankscore = "Not provided"
                vr.hgmd_phen = ""
                vr.hgmd_class = ""
                
        except Exception as e:
            logger.debug(f"HGMD dictionary lookup error: {e}")
            vr.hgmd_dmsupported = 0
            vr.hgmd_publication_count = 0
            vr.hgmd_rankscore = "Not provided"
            vr.hgmd_phen = ""
            vr.hgmd_class = ""
        
        results.extend(hgmd_evidence)

        # --- ensure _clin_sigs_all and _clin_stars_all exist in scope (from earlier collection) ---
        # Fallback: if vr.clinvar_sig already has tokens, add them and try to align stars if vr.clinvar_stars is a ';' list
        existing_sig = getattr(vr, "clinvar_sig", "") or ""
        if existing_sig and str(existing_sig).strip() and str(existing_sig).strip().lower() not in ("not provided", ".", "null", "none"):
            try:
                parts = [p.strip() for p in re.split(r'[;,\|/]+', str(existing_sig)) if p.strip()]
                # try to split existing stars string, if it's a semicolon list
                existing_stars_raw = str(getattr(vr, "clinvar_stars", "") or "")
                if existing_stars_raw and ";" in existing_stars_raw:
                    stars_parts = [p.strip() for p in existing_stars_raw.split(";")]
                else:
                    stars_parts = []
                for idx_p, p2 in enumerate(parts):
                    _clin_sigs_all.append(p2)
                    if idx_p < len(stars_parts):
                        _clin_stars_all.append(stars_parts[idx_p])
                    else:
                        _clin_stars_all.append("")  # unknown alignment -> empty
            except Exception:
                pass

        # Build paired list of (sig, star) to keep alignment robust
        pairs = []
        for i, sig in enumerate(_clin_sigs_all):
            s = (str(sig or "")).strip()
            # guard against missing star entries
            star = _clin_stars_all[i] if i < len(_clin_stars_all) else ""
            star = "" if star in (None, "None") else str(star).strip()
            if not s:
                continue
            pairs.append((s, star))

        # Deduplicate preserving order (case-insensitive)
        uniq_pairs = []
        seen = set()
        for sig, star in pairs:
            key = sig.lower()
            if key in seen:
                continue
            seen.add(key)
            uniq_pairs.append((sig, star))

        # If we have no signatures at all, fallback to existing_sig / default
        if not uniq_pairs:
            if existing_sig and existing_sig.strip() and existing_sig.strip().lower() not in ("not provided", ".", "null", "none"):
                vr.clinvar_sig = existing_sig
            else:
                vr.clinvar_sig = "Not provided"
            if not getattr(vr, "clinvar_stars", ""):
                vr.clinvar_stars = "0"
        else:
            # decide whether definitive labels exist
            def _is_definitive_label(s):
                return bool(re.search(r'\b(likely\s+pathogenic|pathogenic|likely\s+benign|benign)\b', (s or "").lower()))
            has_definitive = any(_is_definitive_label(sig) for sig, _ in uniq_pairs)

            # Filter out conflict tokens only if definitive labels exist; keep star alignment
            final_pairs = []
            for sig, star in uniq_pairs:
                sig_l = sig.lower()
                if has_definitive and re.search(r'\b(conflict|conflicting)\b', sig_l):
                    # skip conflict token when other definitive tokens present
                    continue
                final_pairs.append((sig, star))

            # If filtering removed everything (only conflicts present), keep original uniq_pairs (so conflict preserved)
            if not final_pairs:
                final_pairs = uniq_pairs

            final_sigs = [p for p, _ in final_pairs]
            final_stars = [st for _, st in final_pairs if st not in ("", None)]

            vr.clinvar_sig = ";".join(final_sigs) if final_sigs else "Not provided"
            if final_stars:
                vr.clinvar_stars = ";".join(final_stars)
            else:
                vr.clinvar_stars = getattr(vr, "clinvar_stars", "") or "0"

        # CONFLICT DETECTION
        try:
            hq_sources = []
            for classification, source in results:
                is_hq = False
                if source.startswith("ClinVar_Submission_"):
                    try:
                        stars_match = re.search(r'(\d+)stars', source)
                        if stars_match:
                            stars = int(stars_match.group(1))
                            if stars >= 2:
                                is_hq = True
                    except Exception:
                        pass
                elif source.startswith("HGMD_DM_"):
                    is_hq = True
                elif source == "Internal_DB" and classification in ("Pathogenic", "Likely pathogenic"):
                    is_hq = True
                if is_hq:
                    hq_sources.append((classification, source))
            if len(hq_sources) >= 2:
                classifications = set(c[0].lower() for c in hq_sources)
                has_path = any('pathogen' in c for c in classifications)
                has_benign = any('benign' in c for c in classifications)
                if has_path and has_benign:
                    vr.is_hq_conflict_db = True
                    vr.conflict_details = f"Conflict between HQ sources: {', '.join([f'{c[0]} ({c[1]})' for c in hq_sources])}"
        except Exception as e:
            logger.debug(f"Conflict detection failed: {e}")
        
        vr.db_hits = results
        try:
            aggregate_clinvar_signatures(vr)
        except Exception:
            pass
        vr.db_hits = results
        return results

    # ClinVar protein index builder 
    def _build_clinvar_indexes_from_vcf(self, clinvar_path: str):
        if not clinvar_path or not os.path.exists(clinvar_path):
            logger.warning("ClinVar path invalid: %s", clinvar_path)
            return
        try:
            logger.info("Loading ClinVar from %s with gene validation...", clinvar_path)
            count = 0
            protein_count = 0
            with gzip.open(clinvar_path, "rt", encoding="utf-8", errors="replace") as fh:
                for ln in fh:
                    if ln.startswith("#"): 
                        continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 8: 
                        continue
                    chrom = parts[0]
                    try: 
                        pos = int(parts[1])
                    except: 
                        continue
                    info = parts[7]
                    
                    idx = info.find("GENEINFO=")
                    
                    rest = info[idx+9:]
                    end_idx = rest.find(";")
                    val = rest[:end_idx] if end_idx != -1 else rest
                    if not val:
                        continue
                    genes_raw = []
                    for gene_pair in val.split("|"):
                        if ":" in gene_pair:
                            gene_raw = gene_pair.split(":")[0]
                        else:
                            gene_raw = gene_pair
                        gene_norm = normalize_gene_name(gene_raw)
                        if gene_norm:
                            genes_raw.append((gene_raw, gene_norm))
                    
                    if not genes_raw:
                        continue
                                    
                    # Extract CLNSIG
                    clnsig = ""
                    if "CLNSIG=" in info:
                        m = re.search(r'CLNSIG=([^;]+)', info)
                        if m: 
                            clnsig = m.group(1)
                    
                    # Calculate stars from CLNREVSTAT
                    stars = 0
                    if "CLNREVSTAT=" in info:
                        low = info.lower()
                        if 'practice_guideline' in low: 
                            stars = 4
                        elif 'expert_panel' in low: 
                            stars = 3
                        elif 'criteria_provided' in low: 
                            stars = 1
                            if 'multiple_submitters' in low and 'no_conflict' in low:
                                stars = 2
                    
                    # Extract HGVS_Protein
                    hgvsp_matches = []
                    if "HGVS_Protein=" in info:
                        hgvs_match = re.findall(r'HGVS_Protein=([^;]+)', info)
                        for hgvsp in hgvs_match:
                            if hgvsp.startswith("p."):
                                hgvsp_matches.append(hgvsp)
                    if not hgvsp_matches and "p." in info:
                        p_matches = re.findall(r'p\.[A-Za-z*]+[0-9]+[A-Za-z*]*', info)
                        for m in p_matches:
                            if m.startswith("p."):
                                hgvsp_matches.append(m)
                    
                    # Also try to find p. pattern anywhere in INFO
                    if not hgvsp_matches and "p." in info:
                        old_matches = re.findall(r'(p\.[A-Za-z0-9\(\)\*]+)', info)
                        for mh in old_matches:
                            if mh not in hgvsp_matches:
                                hgvsp_matches.append(mh)
                    
                    ref = parts[3] if len(parts) > 3 else ""
                    alt = parts[4] if len(parts) > 4 else ""
                    
                    # Extract CLNDN (trait) and CLNSIGSCV (SCV ids) if present
                    clndn_matches = re.findall(r'CLNDN=([^;]+)', info)
                    clndn_val = ";".join(clndn_matches) if clndn_matches else ""
                    scv_matches = re.findall(r'CLNSIGSCV=([^;]+)', info)
                    scv_val = ";".join(scv_matches) if scv_matches else ""
                    
                    # Also extract VARIANT field for qvar
                    variant_val = ""
                    if "VARIANT=" in info:
                        m = re.search(r'VARIANT=([^;]+)', info)
                        if m:
                            variant_val = m.group(1)
                    
                    for gene_raw, gene_norm in genes_raw:
                        for hgvsp in hgvsp_matches:
                            entry = {
                                "hgvsp": hgvsp,
                                "clnsig": clnsig,
                                "stars": stars,
                                "chrom": chrom,
                                "pos": pos,
                                "ref": ref,
                                "alt": alt,
                                "gene_raw": gene_raw,
                                "gene_norm": gene_norm,
                                "clndn": clndn_val,
                                "clnsigscv": scv_val,
                                "variant": variant_val,
                                "qvar": f"{chrom}-{pos}-{ref}-{alt}"  # Add qvar for easier lookup
                            }
                            self.clinvar_protein_index[gene_norm].append(entry)
                            protein_count += 1
                    
                    count += 1
                    if count % 1000000 == 0:
                        logger.debug("Processed %d ClinVar variants, %d with protein", count, protein_count)
            
            logger.info("ClinVar loaded with gene validation. Protein variants: %d. Total variants: %d", protein_count, count)
            
            # Also index by qvar for faster lookup
            self.clinvar_by_qvar = defaultdict(list)
            for gene_norm, entries in self.clinvar_protein_index.items():
                for entry in entries:
                    qvar = entry.get("qvar")
                    if qvar:
                        self.clinvar_by_qvar[qvar].append(entry)
            
        except Exception as e:
            logger.warning("ClinVar load error: %s", e)
            import traceback
            traceback.print_exc()
        

    def find_protein_variant_by_hgvsp(self, gene: str, hgvsp: str, pos_aa: int, transcript: str,
                                        current_chrom: Optional[str] = None, current_pos: Optional[int] = None,
                                        current_ref: Optional[str] = None, current_alt: Optional[str] = None,
                                        acmg_disease: Optional[str] = None) -> List[Tuple[str, str, dict]]:
        """
        Find PS1 candidates: ANOTHER variant with SAME protein change.
        Skips Conflict entries.
        """
        if not hgvsp or not gene:
            return []
        
        norm_gene = normalize_gene_name(gene)
        results = []
        seen = set()
        current_hgvsp_clean = clean_hgvsp(hgvsp)
        
        def get_stars_value(stars_value) -> int:
            try:
                return int(stars_value)
            except (ValueError, TypeError):
                try:
                    stars_str = str(stars_value).strip()
                    if stars_str:
                        match = re.search(r'(\d+)', stars_str)
                        if match:
                            return int(match.group(1))
                except Exception:
                    pass
                return 0
        
        def is_conflict(clnsig: str) -> bool:
            """Check if ClinVar classification indicates conflict."""
            clnsig_lower = clnsig.lower()
            return "conflict" in clnsig_lower
        
        def get_ps1_strength(stars: int, clnsig: str) -> Tuple[Optional[str], Optional[int]]:
            """Determine PS1 strength based on stars and classification."""
            clnsig_lower = clnsig.lower()
            
            # Skip conflict entries completely
            if is_conflict(clnsig):
                return (None, None)
            
            is_pathogenic = "pathogenic" in clnsig_lower and "likely" not in clnsig_lower
            is_likely_pathogenic = "likely pathogenic" in clnsig_lower or ("likely" in clnsig_lower and "pathogen" in clnsig_lower)
            
            if stars >= 2:
                if is_pathogenic:
                    return ("PS1_Strong", 4)
                elif is_likely_pathogenic:
                    return ("PS1_Moderate", 2)
            elif stars >= 1:
                if is_pathogenic or is_likely_pathogenic:
                    return ("PS1_Supporting", 1)
            
            return (None, None)
        
        def is_same_genomic_variant(entry_chrom, entry_pos, entry_ref, entry_alt) -> bool:
            if current_chrom is None or current_pos is None:
                return False
            if entry_chrom is None or entry_pos is None:
                return False
            
            entry_chrom_norm = entry_chrom
            current_chrom_norm = current_chrom
            if not entry_chrom_norm.startswith("chr") and current_chrom_norm.startswith("chr"):
                entry_chrom_norm = "chr" + entry_chrom_norm
            elif entry_chrom_norm.startswith("chr") and not current_chrom_norm.startswith("chr"):
                current_chrom_norm = "chr" + current_chrom_norm
            
            try:
                if entry_chrom_norm != current_chrom_norm:
                    return False
                if int(entry_pos) != int(current_pos):
                    return False
                if str(entry_ref or "") != str(current_ref or ""):
                    return False
                if str(entry_alt or "") != str(current_alt or ""):
                    return False
                return True
            except Exception:
                return False
        
        qvar_key = None
        if current_chrom and current_pos and current_ref is not None and current_alt is not None:
            qvar_key = f"{current_chrom}-{current_pos}-{current_ref}-{current_alt}"
        
        # =========================================================
        # Search expanded ClinVar submissions
        # =========================================================
        expanded = []
        if qvar_key and hasattr(self, "clinvar_expanded_by_qvar"):
            expanded = self.get_submissions_for_query_variant(qvar_key, variant_gene=gene)
        
        for sub in expanded:
            try:
                sub_hgvsp = clean_hgvsp(sub.get("hgvsp", "") or "")
                if sub_hgvsp != current_hgvsp_clean:
                    continue
                
                sub_stars = get_stars_value(sub.get("stars", 0))
                if sub_stars < 1:
                    continue
                
                sub_sig = str(sub.get("clnsig", "") or "").strip()
                
                # CRITICAL: Skip Conflict entries
                if is_conflict(sub_sig):
                    logger.debug(f"Skipping PS1: conflict entry for {sub_hgvsp}")
                    continue
                
                sub_chrom = sub.get("chrom", "")
                sub_pos = sub.get("pos", None)
                sub_ref = sub.get("ref", "")
                sub_alt = sub.get("alt", "")
                
                # Skip the SAME genomic variant
                if is_same_genomic_variant(sub_chrom, sub_pos, sub_ref, sub_alt):
                    continue
                
                sub_gene = sub.get("gene", "") or sub.get("gene_raw", "")
                gene_valid, _ = validate_gene_match_db_evidence(gene, sub_gene, "ClinVar")
                if not gene_valid:
                    continue
                
                # Disease/phenotype matching
                if acmg_disease:
                    phen = sub.get("phenotype", "") or ""
                    if phen:
                        matched, _, _ = check_disease_match(acmg_disease, phen)
                        if not matched:
                            continue
                    else:
                        continue
                
                # Determine PS1 strength
                ps1_strength, ps1_points = get_ps1_strength(sub_stars, sub_sig)
                if ps1_strength is None:
                    continue
                
                key = f"{sub_sig}|{sub.get('submission_id', '')}|{ps1_strength}"
                if key not in seen:
                    seen.add(key)
                    results.append((sub_sig, ps1_strength, ps1_points, sub))
                    
            except Exception:
                continue
        
        # =========================================================
        # Fallback to protein_index
        # =========================================================
        if not results and hasattr(self, "clinvar_protein_index"):
            for entry in self.clinvar_protein_index.get(norm_gene, []):
                entry_hgvsp = clean_hgvsp(entry.get('hgvsp', ''))
                if not entry_hgvsp or entry_hgvsp != current_hgvsp_clean:
                    continue
                
                entry_stars = get_stars_value(entry.get('stars', 0))
                if entry_stars < 1:
                    continue
                
                clnsig_raw = str(entry.get('clnsig', '') or "").strip()
                
                # CRITICAL: Skip Conflict entries
                if is_conflict(clnsig_raw):
                    logger.debug(f"Skipping PS1 (protein_index): conflict entry for {entry_hgvsp}")
                    continue
                
                entry_chrom = entry.get('chrom', '')
                entry_pos = entry.get('pos', None)
                entry_ref = entry.get('ref', '')
                entry_alt = entry.get('alt', '')
                
                # Skip the SAME genomic variant
                if is_same_genomic_variant(entry_chrom, entry_pos, entry_ref, entry_alt):
                    continue
                
                # Determine PS1 strength
                ps1_strength, ps1_points = get_ps1_strength(entry_stars, clnsig_raw)
                if ps1_strength is None:
                    continue
                
                sub_dict = {
                    "hgvsp": entry_hgvsp,
                    "clnsig": clnsig_raw,
                    "stars": entry_stars,
                    "ps1_strength": ps1_strength,
                    "ps1_points": ps1_points,
                    "submission_id": f"protein_index_{entry.get('chrom', '')}_{entry.get('pos', '')}",
                    "chrom": entry_chrom,
                    "pos": entry_pos,
                    "ref": entry_ref,
                    "alt": entry_alt,
                    "gene": entry.get('gene_raw', '')
                }
                
                key = f"{clnsig_raw}|{sub_dict['submission_id']}|{ps1_strength}"
                if key not in seen:
                    seen.add(key)
                    results.append((clnsig_raw, ps1_strength, ps1_points, sub_dict))
        
        return results

    def check_pm5_variants(self, gene: str, pos_aa: int, current_hgvsp: str, 
                           current_chrom: str = None, current_pos: int = None,
                           current_ref: str = None, current_alt: str = None) -> List[Dict[str, Any]]:
        """
        Checks for other pathogenic/likely-pathogenic variants at the same amino-acid residue.
        Skips entries that represent the exact same genomic variant as current.
        """
        if not pos_aa or not gene:
            return []

        norm_gene = normalize_gene_name(gene)
        results = []

        # iterate ClinVar protein-index for the gene
        for entry in self.clinvar_protein_index.get(norm_gene, []):

            # skip if entry lacks protein annotation or gene mismatch
            entry_hgvsp = entry.get('hgvsp', '')
            if not entry_hgvsp:
                continue

            # Skip exact same protein change (PS1 case)
            try:
                if clean_hgvsp(entry_hgvsp) == clean_hgvsp(current_hgvsp):
                    continue
            except Exception:
                pass
            
            # Skip if this is the exact same genomic variant (same coordinates)
            try:
                entry_chrom = entry.get('chrom', '')
                entry_pos = entry.get('pos', None)
                entry_ref = entry.get('ref', '')
                entry_alt = entry.get('alt', '')
                if (current_chrom and current_pos and current_ref is not None and current_alt is not None):
                    if entry_chrom and entry_pos and entry_ref and entry_alt:
                        # Normalize chromosome
                        ec = entry_chrom
                        if not ec.startswith("chr") and current_chrom.startswith("chr"):
                            ec = "chr" + ec
                        if (ec == current_chrom and int(entry_pos) == int(current_pos) and 
                            entry_ref == current_ref and entry_alt == current_alt):
                            # Same genomic variant - skip for PM5
                            continue
            except Exception:
                pass

            # Parse position from entry HGVSp
            entry_pos_aa = extract_aa_position_from_hgvsp(entry_hgvsp)
            
            if not entry_pos_aa or entry_pos_aa != pos_aa:
                continue

            # check significance & stars
            clnsig_raw = str(entry.get('clnsig', '') or "").strip()
            stars = int(entry.get('stars', 0) or 0)
            if stars < 2:
                continue

            clnsig_norm = re.sub(r'[_\s]+', ' ', clnsig_raw.lower().strip())

            # skip explicit benign entries
            if re.search(r'\bbenign\b', clnsig_norm) and 'pathogen' not in clnsig_norm:
                continue

            # determine best_class with conflict-handling and combined token rule
            if re.search(r'\b(conflict|conflicting|conflicting_classifications)\b', clnsig_norm):
                # skip conflicting as pathogenic evidence for PM5
                continue

            has_likely = bool(re.search(r'\blikely\s+pathogenic\b', clnsig_norm))
            has_path = bool(re.search(r'\bpathogenic\b', clnsig_norm))

            if has_path and has_likely:
                best_class = "Pathogenic"
            elif has_path:
                best_class = "Pathogenic"
            elif has_likely:
                best_class = "Likely pathogenic"
            else:
                best_class = "Unknown"
            
            # Add detailed entry to results
            results.append({
                "hgvsp": entry_hgvsp,
                "clnsig": entry.get('clnsig', ''),
                "best": best_class,
                "stars": stars,
                "chrom": entry.get('chrom', ''),
                "pos": entry.get('pos', ''),
                "gene_raw": entry.get('gene_raw', '')
            })

        return results

# ---------------------------
# Exons table & NMD logic
# ---------------------------
def load_exons_table(path: str) -> Optional[pd.DataFrame]:
    if not path or not os.path.exists(path):
        logger.warning(f"Exons file not found: {path}. NMD prediction disabled.")
        return None
    if not HAVE_PANDAS:
        logger.warning("Pandas not installed. NMD prediction disabled.")
        return None
    try:
        logger.info(f"Loading exon structure from {path}...")
        df = pd.read_csv(path, dtype=str, sep=None, engine='python', on_bad_lines='skip').fillna("")
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')
        column_map = {}
        for col in ['transcript_id', 'transcript', 'transcriptid']:
            if col in df.columns:
                column_map['transcript'] = col; break
        for col in ['gene', 'gene_symbol', 'symbol', 'gene_name']:
            if col in df.columns:
                column_map['gene'] = col; break
        for col in ['exon', 'exon_number', 'exon_num']:
            if col in df.columns:
                column_map['exon'] = col; break
        for col in ['genomic_start', 'start', 'start_pos', 'position_start']:
            if col in df.columns:
                column_map['start'] = col; break
        for col in ['genomic_end', 'end', 'end_pos', 'position_end']:
            if col in df.columns:
                column_map['end'] = col; break
        required = ['transcript', 'exon', 'start', 'end']
        missing = [col for col in required if col not in column_map]
        if missing:
            logger.error(f"Exons file missing required columns: {missing}")
            logger.error(f"Available columns: {list(df.columns)}")
            return None
        df_clean = pd.DataFrame()
        df_clean['Transcript_ID'] = df[column_map['transcript']]
        df_clean['Gene'] = df[column_map['gene']] if 'gene' in column_map else ''
        df_clean['Exon'] = df[column_map['exon']]
        df_clean['Start'] = pd.to_numeric(df[column_map['start']], errors='coerce')
        df_clean['End'] = pd.to_numeric(df[column_map['end']], errors='coerce')
        df_clean = df_clean.dropna(subset=["Start", "End"])
        df_clean = df_clean[df_clean['Start'] > 0]
        df_clean = df_clean[df_clean['End'] > df_clean['Start']]
        def extract_exon_number(val):
            if pd.isna(val):
                return None
            s = str(val).lower()
            patterns = [r'exon(\d+)', r'exon_(\d+)', r'ex(\d+)', r'e(\d+)', r'(\d+)']
            for pattern in patterns:
                match = re.search(pattern, s)
                if match:
                    try:
                        return int(match.group(1))
                    except:
                        continue
            return None
        df_clean['Exon_Number'] = df_clean['Exon'].apply(extract_exon_number)
        missing_exons = df_clean['Exon_Number'].isna()
        if missing_exons.any():
            for tx_id in df_clean[missing_exons]['Transcript_ID'].unique():
                tx_mask = df_clean['Transcript_ID'] == tx_id
                tx_exons = df_clean[tx_mask].sort_values('Start')
                for i, (idx, _) in enumerate(tx_exons.iterrows(), 1):
                    df_clean.at[idx, 'Exon_Number'] = i
        df_clean = df_clean.dropna(subset=['Exon_Number'])
        df_clean['Exon_Number'] = df_clean['Exon_Number'].astype(int)
        transcript_last_exons = {}
        for tx_id, group in df_clean.groupby('Transcript_ID'):
            if not group['Exon_Number'].isna().all():
                last_exon_num = group['Exon_Number'].max()
                last_exons = group[group['Exon_Number'] == last_exon_num]
                if not last_exons.empty:
                    last_exon = last_exons.iloc[0]
                    transcript_last_exons[tx_id] = {
                        'start': int(last_exon['Start']),
                        'end': int(last_exon['End']),
                        'exon_number': int(last_exon['Exon_Number']),
                        'gene': last_exon.get('Gene', '')
                    }
        df_clean = df_clean.assign(_transcript_last_exons=transcript_last_exons)
        logger.info(f"Loaded exon data for {len(transcript_last_exons)} transcripts")
        return df_clean
    except Exception as e:
        logger.error(f"Failed to load exons file {path}: {e}", exc_info=True)
        return None

def compute_nmd_internal(vr: VariantRecord, exons_df: pd.DataFrame) -> str:
    relevant_consequences = {"stop_gained", "frameshift_variant", "frameshift", "splice_acceptor_variant", "splice_donor_variant"}
    if not any(t in (vr.consequence or "").lower() for t in relevant_consequences):
        return ""
    if exons_df is None or exons_df.empty:
        return ""
    tx_id = vr.transcript or vr._raw_ann.get("Feature") if vr._raw_ann else ""
    if not tx_id:
        return ""
    tx_base = tx_id.split('.')[0] if '.' in tx_id else tx_id
    if not hasattr(exons_df, '_transcript_last_exons'):
        logger.warning("exons_df doesn't have _transcript_last_exons attribute")
        return ""
    last_exon_info = None
    if tx_id in exons_df._transcript_last_exons:
        last_exon_info = exons_df._transcript_last_exons[tx_id]
    elif tx_base in exons_df._transcript_last_exons:
        last_exon_info = exons_df._transcript_last_exons[tx_base]
    if not last_exon_info:
        gene_exons = exons_df[exons_df['Gene'] == vr.gene_raw] if 'Gene' in exons_df.columns else pd.DataFrame()
        if not gene_exons.empty:
            gene_transcripts = gene_exons['Transcript_ID'].unique()
            for gtx in gene_transcripts:
                if tx_base in gtx or tx_id in gtx:
                    last_exon_info = exons_df._transcript_last_exons.get(gtx)
                    if last_exon_info:
                        break
        if not last_exon_info:
            logger.debug(f"No exon data found for transcript {tx_id} (gene {vr.gene})")
            return ""
    variant_pos = int(vr.pos)
    last_exon_start = last_exon_info['start']
    last_exon_end = last_exon_info['end']
    logger.debug(f"NMD check for {vr.gene}:{tx_id} at {variant_pos} (exon {last_exon_start}-{last_exon_end})")
    if last_exon_start <= variant_pos <= last_exon_end:
        return "Escaping_NMD_Last_Exon"
    else:
        return "Triggering_NMD"


def load_transcript_map(path: str) -> Dict[str, Any]:
    """
    Load mapping table between RefSeq (NM_) and Ensembl (ENST_) transcripts.

    Returns dict:
      {
        "ensembl_to_refseq": { "ENST...": ["NM_..."], ... },
        "refseq_to_ensembl": { "NM_...": ["ENST..."], ... },
        "mane": { "ENST...": "NM_..." }  # if indicated
      }
    Accepts CSV/TSV. Column fuzzy-detection applied.
    """
    out = {"ensembl_to_refseq": {}, "refseq_to_ensembl": {}, "mane": {}}
    if not path or not os.path.exists(path):
        logger.debug("Transcript map not provided or not found: %s", path)
        return out

    def _add_pair(enst, nm, is_mane=False):
        if not enst or not nm:
            return
        out["ensembl_to_refseq"].setdefault(enst, [])
        if nm not in out["ensembl_to_refseq"][enst]:
            out["ensembl_to_refseq"][enst].append(nm)
        out["refseq_to_ensembl"].setdefault(nm, [])
        if enst not in out["refseq_to_ensembl"][nm]:
            out["refseq_to_ensembl"][nm].append(enst)
        if is_mane:
            out["mane"][enst] = nm

    # Try pandas  
    if HAVE_PANDAS:
        try:
            df = pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")
            cols_lower = {c.lower(): c for c in df.columns}
            # possible column names
            ref_candidates = ["rna_nucleotide_accession.version", "transcript", "transcript_id", "refseq", "rna_nucleotide"]
            enst_candidates = ["ensembl_rna_identifier", "ensembl_rna", "enst", "ensembl"]
            mane_candidates = ["mane_select", "mane", "mane_select_y", "is_mane_select"]
            # find actual column names
            def find_col(cands):
                for k in cands:
                    if k in cols_lower:
                        return cols_lower[k]
                # fuzzy substring match
                for k in cols_lower:
                    for cand in cands:
                        if cand in k:
                            return cols_lower[k]
                return None
            ref_col = find_col(ref_candidates)
            enst_col = find_col(enst_candidates)
            mane_col = find_col(mane_candidates)
            for _, r in df.iterrows():
                nm = ""
                enst = ""
                try:
                    if ref_col:
                        nm = str(r.get(ref_col, "")).strip()
                    else:
                        # try common fallback columns
                        for c in ("RNA_nucleotide_accession.version","RNA_nucleotide"):
                            if c in r:
                                nm = str(r.get(c, "")).strip()
                                break
                    if enst_col:
                        enst = str(r.get(enst_col, "")).strip()
                    if not nm and not enst:
                        # try alternative columns
                        continue
                    is_mane = False
                    if mane_col:
                        mv = str(r.get(mane_col, "")).strip().lower()
                        if mv in ("true", "yes", "y", "1", "t"):
                            is_mane = True
                    # normalize empty strings
                    if nm:
                        nm = nm.strip()
                    if enst:
                        enst = enst.strip()
                    # Add both if present
                    if enst and nm:
                        _add_pair(enst, nm, is_mane=is_mane)
                except Exception:
                    continue
            logger.info("Loaded transcript map: %d Ensembl entries", len(out["ensembl_to_refseq"]))
            return out
        except Exception as e:
            logger.warning("Pandas read failed for transcript map %s: %s", path, e)

    # Fallback: manual CSV parse
    try:
        with open(path, "rt", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            # if fields seem not tab-separated, try comma
            if reader.fieldnames is None or len(reader.fieldnames) == 1:
                fh.seek(0)
                reader = csv.DictReader(fh)
            for row in reader:
                # try to find NM and ENST in row keys/values
                nm = ""
                enst = ""
                for k, v in row.items():
                    if not v:
                        continue
                    kn = k.lower()
                    if "nm_" in str(v).lower() or (kn and "refseq" in kn) or ("rna_nucleotide" in kn):
                        if not nm:
                            nm = str(v).strip()
                    if "enst" in str(v).lower() or (kn and "ensembl" in kn):
                        if not enst:
                            enst = str(v).strip()
                if enst and nm:
                    _add_pair(enst, nm, is_mane=False)
        logger.info("Loaded transcript map (fallback parser): %d Ensembl entries", len(out["ensembl_to_refseq"]))
        return out
    except Exception as e:
        logger.warning("Failed to parse transcript map %s: %s", path, e)
        return out


def load_nmd_table(path: str) -> Dict[str, Dict[str, Any]]:
    """
    Load a precomputed NMD/exon table (one-row-per-exon or aggregated) into a mapping keyed by RefSeq transcript ID.

    Expected minimal columns (case-insensitive): Transcript_ID (or transcript), number_cds_exon / Exon_Number, genomic_start/genomic_end optional,
    CDS_Exon_length, CDS_length, Protein_length, strand, Chrom, total_cds_exons, Total_Exons

    Returns:
      nmd_map = {
         "NM_0123.4": {
             "chrom": "chr1",
             "strand": "+",
             "total_exons": int,
             "total_cds_exons": int,
             "cds_length": int,
             "protein_length": int,
             "exons": [
                 {"exon_number": int, "cds_exon_length": int, "genomic_start": int or None, "genomic_end": int or None},
                 ...
             ],
             "last_exon": {"exon_number": n, "genomic_start": int, "genomic_end": int}  # optional
         }, ...
      }
    """
    nmd_map: Dict[str, Dict[str, Any]] = {}
    if not path or not os.path.exists(path):
        logger.debug("NMD table not provided: %s", path)
        return nmd_map

    # Try pandas 
    if HAVE_PANDAS:
        try:
            df = pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")
            cols = {c.lower().strip(): c for c in df.columns}
            # find columns
            tx_col = cols.get("transcript_id") or cols.get("transcript") or cols.get("name")
            chrom_col = cols.get("chrom") or cols.get("chr")
            strand_col = cols.get("strand")
            exon_num_col = (cols.get("exon_number") or cols.get("number_cds_exon") or cols.get("exon"))
            exon_len_col = cols.get("exon_length") or cols.get("cds_exon_length")
            cds_exon_len_col = cols.get("cds_exon_length")
            cds_length_col = cols.get("cds_length")
            protein_length_col = cols.get("protein_length")
            total_exons_col = cols.get("total_exons")
            total_cds_exons_col = cols.get("total_cds_exons")
            start_col = cols.get("genomic_start") or cols.get("start") or cols.get("genomicstart")
            end_col = cols.get("genomic_end") or cols.get("end") or cols.get("genomicend")
            # iterate rows
            for _, r in df.iterrows():
                tx = str(r.get(tx_col, "")).strip() if tx_col else ""
                if not tx:
                    continue
                entry = nmd_map.setdefault(tx, {"exons": []})
                if chrom_col and r.get(chrom_col):
                    entry["chrom"] = str(r.get(chrom_col)).strip()
                if strand_col and r.get(strand_col):
                    entry["strand"] = str(r.get(strand_col)).strip()
                try:
                    if total_exons_col and r.get(total_exons_col):
                        entry["total_exons"] = int(float(r.get(total_exons_col)))
                except Exception:
                    pass
                try:
                    if total_cds_exons_col and r.get(total_cds_exons_col):
                        entry["total_cds_exons"] = int(float(r.get(total_cds_exons_col)))
                except Exception:
                    pass
                try:
                    if cds_length_col and r.get(cds_length_col):
                        entry["cds_length"] = int(float(r.get(cds_length_col)))
                except Exception:
                    pass
                try:
                    if protein_length_col and r.get(protein_length_col):
                        entry["protein_length"] = int(float(r.get(protein_length_col)))
                except Exception:
                    pass
                exon_number = None
                try:
                    if exon_num_col and r.get(exon_num_col) not in ("", "."):
                        exon_number = int(float(r.get(exon_num_col)))
                except Exception:
                    exon_number = None
                cds_exon_len = None
                try:
                    if cds_exon_len_col and r.get(cds_exon_len_col) not in ("", "."):
                        cds_exon_len = int(float(r.get(cds_exon_len_col)))
                    elif exon_len_col and r.get(exon_len_col) not in ("", "."):
                        cds_exon_len = int(float(r.get(exon_len_col)))
                except Exception:
                    cds_exon_len = None
                genomic_start = None
                genomic_end = None
                try:
                    if start_col and r.get(start_col) not in ("", "."):
                        genomic_start = int(float(r.get(start_col)))
                    if end_col and r.get(end_col) not in ("", "."):
                        genomic_end = int(float(r.get(end_col)))
                except Exception:
                    genomic_start = genomic_start
                    genomic_end = genomic_end
                exon_entry = {
                    "exon_number": exon_number,
                    "cds_exon_length": cds_exon_len,
                    "genomic_start": genomic_start,
                    "genomic_end": genomic_end
                }
                entry["exons"].append(exon_entry)
            # postprocess: compute last_exon if possible
            for tx, info in nmd_map.items():
                exs = [e for e in info.get("exons", []) if e.get("exon_number") is not None]
                if exs:
                    exs_sorted = sorted(exs, key=lambda x: int(x["exon_number"]))
                    info["exons"] = exs_sorted
                    last = exs_sorted[-1]
                    info["last_exon"] = {
                        "exon_number": last.get("exon_number"),
                        "genomic_start": last.get("genomic_start"),
                        "genomic_end": last.get("genomic_end")
                    }
            logger.info("Loaded NMD table for %d transcripts from %s", len(nmd_map), path)
            return nmd_map
        except Exception as e:
            logger.warning("Pandas failed to parse NMD table %s: %s", path, e)

    # Fallback CSV parse if pandas not available
    try:
        with open(path, "rt", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or len(reader.fieldnames) == 1:
                fh.seek(0)
                reader = csv.DictReader(fh)
            for r in reader:
                tx = r.get("Transcript_ID") or r.get("transcript") or r.get("Transcript")
                if not tx:
                    continue
                tx = tx.strip()
                entry = nmd_map.setdefault(tx, {"exons": []})
                chrom = r.get("Chrom") or r.get("chrom") or r.get("chr")
                if chrom:
                    entry["chrom"] = chrom.strip()
                strand = r.get("strand") or r.get("STRAND")
                if strand:
                    entry["strand"] = strand.strip()
                try:
                    te = r.get("Total_Exons") or r.get("total_exons") or r.get("TotalExons")
                    if te:
                        entry["total_exons"] = int(float(te))
                except Exception:
                    pass
                try:
                    tce = r.get("total_cds_exons")
                    if tce:
                        entry["total_cds_exons"] = int(float(tce))
                except Exception:
                    pass
                try:
                    cds_len = r.get("CDS_length")
                    if cds_len:
                        entry["cds_length"] = int(float(cds_len))
                except Exception:
                    pass
                try:
                    prot_len = r.get("Protein_length")
                    if prot_len:
                        entry["protein_length"] = int(float(prot_len))
                except Exception:
                    pass
                exon_number = None
                try:
                    en = r.get("number_cds_exon") or r.get("Exon_Number") or r.get("exon")
                    if en:
                        exon_number = int(float(en))
                except Exception:
                    exon_number = None
                genomic_start = None; genomic_end = None
                try:
                    gs = r.get("Genomic_start") or r.get("genomic_start") or r.get("Start")
                    ge = r.get("Genomic_end") or r.get("genomic_end") or r.get("End")
                    if gs:
                        genomic_start = int(float(gs))
                    if ge:
                        genomic_end = int(float(ge))
                except Exception:
                    genomic_start = genomic_start; genomic_end = genomic_end
                cds_exon_len = None
                try:
                    cel = r.get("CDS_Exon_length") or r.get("CDS_Exon_length")
                    if cel:
                        cds_exon_len = int(float(cel))
                except Exception:
                    pass
                exon_entry = {"exon_number": exon_number, "cds_exon_length": cds_exon_len, "genomic_start": genomic_start, "genomic_end": genomic_end}
                entry["exons"].append(exon_entry)
        # postprocess last exon
        for tx, info in nmd_map.items():
            exs = [e for e in info.get("exons", []) if e.get("exon_number") is not None]
            if exs:
                exs_sorted = sorted(exs, key=lambda x: int(x["exon_number"]))
                info["exons"] = exs_sorted
                last = exs_sorted[-1]
                info["last_exon"] = {"exon_number": last.get("exon_number"), "genomic_start": last.get("genomic_start"), "genomic_end": last.get("genomic_end")}
        logger.info("Loaded NMD table (fallback parse) for %d transcripts", len(nmd_map))
        return nmd_map
    except Exception as e:
        logger.warning("Failed to parse NMD table %s: %s", path, e)
        return nmd_map


def _populate_genomic_exon_coords(nmd_map: Dict[str, Dict[str, Any]], exons_map: Any) -> None:
    """
    Fill genomic_start/genomic_end in nmd_map['exons'] from exons_map.
    exons_map may be:
      - pandas DataFrame returned by load_exons_table (has attribute _transcript_last_exons and rows)
      - dict keyed by RefSeq transcript with 'exons' list elements containing 'start'/'end'
    This function mutates nmd_map in-place.
    """
    if not nmd_map:
        return
    # Helper to fetch exons list for a refseq transcript
    def fetch_exons_for_tx(tx_key: str):
        # try pandas frame style
        if HAVE_PANDAS and isinstance(exons_map, pd.DataFrame):
            # parse prefix
            tx_full = tx_key
            tx_base = tx_key.split(".")[0] if "." in tx_key else tx_key
            # find rows
            df = exons_map
            candidates = []
            # check exact matches in Transcript_ID col if exists
            if "Transcript_ID" in df.columns:
                candidates = df[df["Transcript_ID"].astype(str).str.startswith(tx_base)]
            if candidates is None or len(candidates) == 0:
                # try matching gene etc - fallback
                return []
            # assemble exons list
            ex_list = []
            for _, row in candidates.sort_values("Start").iterrows():
                try:
                    s = int(row["Start"]); e = int(row["End"])
                except Exception:
                    continue
                ex_list.append({"start": s, "end": e, "exon_number": int(row.get("Exon_Number")) if "Exon_Number" in row and row.get("Exon_Number") else None})
            return ex_list
        # if dict-like
        if isinstance(exons_map, dict):
            # try exact key or base
            if tx_key in exons_map:
                info = exons_map[tx_key]
                return info.get("exons", [])
            base = tx_key.split(".")[0]
            if base in exons_map:
                return exons_map[base].get("exons", [])
        # otherwise no coordinates
        return []

    for tx, info in nmd_map.items():
        try:
            exs = info.get("exons", [])
            if not exs:
                continue
            # attempt to get genomic exons from exons_map
            ref_exs = fetch_exons_for_tx(tx)
            if not ref_exs:
                continue
            # try to match by exon_number if present
            if any(e.get("exon_number") for e in exs):
                # map ref_exs by exon index (1-based)
                for e in exs:
                    exnum = e.get("exon_number")
                    if exnum and 1 <= exnum <= len(ref_exs):
                        ref = ref_exs[exnum - 1]
                        if ref:
                            e["genomic_start"] = ref.get("start") or e.get("genomic_start")
                            e["genomic_end"] = ref.get("end") or e.get("genomic_end")
            else:
                # try to align by order
                for i, e in enumerate(exs):
                    if i < len(ref_exs):
                        ref = ref_exs[i]
                        e["genomic_start"] = ref.get("start") or e.get("genomic_start")
                        e["genomic_end"] = ref.get("end") or e.get("genomic_end")
            # update last_exon if missing
            le = info.get("last_exon")
            if not le:
                # find exon with highest exon_number
                with_coords = [e for e in exs if e.get("exon_number") is not None]
                if with_coords:
                    last = sorted(with_coords, key=lambda x: int(x["exon_number"]))[-1]
                    info["last_exon"] = {"exon_number": last.get("exon_number"), "genomic_start": last.get("genomic_start"), "genomic_end": last.get("genomic_end")}
                else:
                    # fallback: last element
                    last = exs[-1]
                    info["last_exon"] = {"exon_number": last.get("exon_number"), "genomic_start": last.get("genomic_start"), "genomic_end": last.get("genomic_end")}
        except Exception:
            continue


# ---------------------------
# Variant NMD predictor (updated)
# ---------------------------
def variant_triggers_nmd(vr: Any,
                         nmd_map: Dict[str, Dict[str, Any]],
                         exons_map: Any,
                         tx_map: Dict[str, Any],
                         penultimate_nt: int = 50
                         ) -> Tuple[bool, Dict[str, Any]]:
    """
    Predict whether a variant triggers nonsense-mediated decay (NMD).
    
    NMD prediction is critical for PVS1 strength assignment. This function uses:
      - Transcript-aware mapping (RefSeq/Ensembl) via tx_map.
      - Exon structure from nmd_map or exons_map for last-exon detection.
      - HGVSc cDNA coordinates for precise exon mapping.
      - Conservative logic: if any transcript suggests NMD escape, variant is considered NMD-escaped.
    
    Args:
        vr (VariantRecord): Variant to evaluate.
        nmd_map (Dict): Pre-loaded NMD/exon table keyed by transcript.
        exons_map (Any): Exon structure DataFrame or dict.
        tx_map (Dict): Transcript mapping between Ensembl and RefSeq.
        penultimate_nt (int): Terminal window size for penultimate exon (default 50nt).
    
    Returns:
        Tuple[bool, Dict]: (triggers_nmd, details_dict)
        - triggers_nmd (bool): True if variant predicted to trigger NMD.
        - details_dict (Dict): Contains per-transcript predictions and reasoning.
    """
    details: Dict[str, Any] = {"checked_transcripts": [], "reason": "", "variant_pos": None}
    try:
        chrom = getattr(vr, "chrom", None)
        pos = int(getattr(vr, "pos"))
        details["variant_pos"] = f"{chrom}:{pos}"
    except Exception:
        details["reason"] = "no_variant_position"
        return False, details

    try:
        exon_field = getattr(vr, "exon", "") or ""
        if exon_field:
            # typical formats: "42/42", "3/8", sometimes "exon 3/8" or "3/8 (something)"
            m = re.search(r'(\d+)\s*/\s*(\d+)', str(exon_field))
            if m:
                curr = int(m.group(1)); total = int(m.group(2))
                if total > 0 and curr == total:
                    details["reason"] = "variant_in_last_exon_by_EXON_field"
                    return False, details
    except Exception:
        # if this check fails for any reason, continue to normal logic
        pass

    # proceed with the rest of the transcript-aware prediction (unchanged, conservative combination)
    # Build list of RefSeq candidate transcripts for this variant.
    candidates: List[str] = []

    tx_raw = getattr(vr, "transcript", "") or ""
    if tx_raw.upper().startswith("NM_"):
        candidates.append(tx_raw)

    # Map Ensembl -> RefSeq if tx_map provided
    if tx_map and isinstance(tx_map, dict) and tx_raw:
        enst_full = tx_raw
        enst_base = enst_full.split(".")[0] if "." in enst_full else enst_full
        if "ensembl_to_refseq" in tx_map:
            m = tx_map["ensembl_to_refseq"].get(enst_full) or tx_map["ensembl_to_refseq"].get(enst_base)
            if m:
                for nm in m:
                    if nm not in candidates:
                        candidates.append(nm)
        if "mane" in tx_map and enst_full in tx_map["mane"]:
            mn = tx_map["mane"].get(enst_full)
            if mn and mn not in candidates:
                candidates.insert(0, mn)

    # Try to extract refseq candidates from nmd_map keys if available
    if not candidates and nmd_map:
        tx_base = tx_raw.split(".")[0] if tx_raw else None
        for key in list(nmd_map.keys()):
            if tx_base and key.startswith(tx_base):
                if key not in candidates:
                    candidates.append(key)

    # Fallback: if still no candidates, attempt to find transcripts in nmd_map by gene
    if not candidates and nmd_map:
        gene = getattr(vr, "gene", "") or getattr(vr, "gene_raw", "")
        if gene:
            gene_norm = gene.upper() if isinstance(gene, str) else gene
            for key, info in nmd_map.items():
                try:
                    ig = info.get("gene") or info.get("Gene") or ""
                    if ig and gene_norm == str(ig).upper():
                        if key not in candidates:
                            candidates.append(key)
                except Exception:
                    continue

    # NEW FALLBACK: if still no candidates, but exons_map (pandas) is available -> attempt to use it directly
    if not candidates and HAVE_PANDAS and isinstance(exons_map, pd.DataFrame) and hasattr(exons_map, "_transcript_last_exons"):
        tx_full = tx_raw
        tx_base = tx_full.split(".")[0] if tx_full else None
        tlast = None
        if tx_full and tx_full in exons_map._transcript_last_exons:
            tlast = tx_full
        elif tx_base and tx_base in exons_map._transcript_last_exons:
            tlast = tx_base
        if not tlast:
            gene_raw = getattr(vr, "gene_raw", "") or getattr(vr, "gene", "")
            if gene_raw:
                df = exons_map
                cand_rows = df[df["Gene"].astype(str).str.upper() == str(gene_raw).upper()]
                if not cand_rows.empty:
                    tlast = cand_rows.iloc[0]["Transcript_ID"]
        if tlast:
            candidates.append(tlast)

    # Deduplicate
    candidates = list(dict.fromkeys(candidates))

    if not candidates:
        details["reason"] = "no_refseq_candidates"
        return False, details  # conservative: do not predict NMD

    # Ensure genomic coords populated in nmd_map where possible
    try:
        _populate_genomic_exon_coords(nmd_map, exons_map)
    except Exception:
        pass

    per_tx_results = []
    for ref_tx in candidates:
        tx_result = {"transcript": ref_tx, "predicts_nmd": None, "details": {}}
        try:
            info = nmd_map.get(ref_tx) or nmd_map.get(ref_tx.split(".")[0])
            last_exon = None
            exons_list = []
            strand = "+"

            # If nmd_map has info - use it
            if info:
                last_exon = info.get("last_exon")
                exons_list = info.get("exons", []) or []
                strand = info.get("strand") or info.get("strand", "+")
            else:
                # Try to get last exon from exons_map (pandas frame) as fallback
                if HAVE_PANDAS and isinstance(exons_map, pd.DataFrame) and hasattr(exons_map, "_transcript_last_exons"):
                    base = ref_tx.split(".")[0]
                    last_ex = exons_map._transcript_last_exons.get(ref_tx) or exons_map._transcript_last_exons.get(base)
                    if last_ex:
                        last_exon = {"exon_number": last_ex.get("exon_number"), "genomic_start": last_ex.get("start"), "genomic_end": last_ex.get("end")}
                        df = exons_map
                        rows = df[df["Transcript_ID"].astype(str).str.startswith(base)]
                        exs = []
                        if not rows.empty:
                            for _, rr in rows.sort_values("Start").iterrows():
                                exs.append({"exon_number": int(rr["Exon_Number"]) if rr["Exon_Number"] else None, "genomic_start": int(rr["Start"]), "genomic_end": int(rr["End"])})
                            exons_list = exs
                elif isinstance(exons_map, dict):
                    em = exons_map.get(ref_tx) or exons_map.get(ref_tx.split(".")[0])
                    if em:
                        last_exon = em.get("last_exon")
                        exons_list = em.get("exons", [])
                        strand = em.get("strand", strand)

            # If still no last_exon coords -> unknown for this transcript
            if not last_exon or last_exon.get("genomic_start") is None or last_exon.get("genomic_end") is None:
                tx_result["predicts_nmd"] = None
                tx_result["details"]["reason"] = "no_last_exon_coords"
                per_tx_results.append(tx_result)
                continue

            last_start = int(last_exon.get("genomic_start"))
            last_end = int(last_exon.get("genomic_end"))

            # Variant in last exon -> escapes NMD
            if last_start <= pos <= last_end:
                tx_result["predicts_nmd"] = False
                tx_result["details"]["reason"] = "variant_in_last_exon"
                per_tx_results.append(tx_result)
                continue

            # ----------------------------
            # Splice-site aware exon-number inference using HGVSc c. position (preferred),
            # falling back to genomic-proximity inference if needed.
            # ----------------------------
            exon_num = None

            # 0) Try to parse cDNA coordinate from HGVSc and map to exon number using nmd_map (preferred)
            try:
                hgvsc_raw = getattr(vr, "hgvsc", "") or ""
                cdna_pos, intronic_offset = parse_cdna_from_hgvsc(hgvsc_raw)
                if cdna_pos is not None:
                    # Use nmd_map (if available) and tx_map (if provided) to map cdna_pos -> exon_number
                    try:
                        exon_from_cdna = None
                        if nmd_map:
                            # transcript id may be vr.transcript or info in _raw_ann Feature
                            tx_id = getattr(vr, "transcript", "") or (vr._raw_ann.get("Feature") if getattr(vr,"_raw_ann",None) else "")
                            exon_from_cdna = map_cdna_to_exon_number(tx_id, cdna_pos, nmd_map, tx_map, intronic_offset=intronic_offset)
                        if exon_from_cdna:
                            exon_num = exon_from_cdna
                        else:
                            logger.debug("cdna->exon mapping failed for %s (transcript=%s, cdna=%s, intr=%s). nmd_map keys: %d. exons_map present: %s",
                                        getattr(vr, "hgvsp", "?"), tx_id, cdna_pos, intronic_offset, len(nmd_map) if nmd_map else 0, isinstance(exons_map, pd.DataFrame))
                    except Exception:
                        exon_num = None
            except Exception:
                exon_num = None

            # 1) If exon_num still unknown, try EXON field like "3/8"
            if exon_num is None:
                try:
                    exon_field = getattr(vr, "exon", "") or ""
                    if exon_field:
                        m_ex = re.search(r'(\d+)\s*/\s*(\d+)', str(exon_field))
                        if m_ex:
                            try:
                                exon_num = int(m_ex.group(1))
                            except Exception:
                                exon_num = None
                except Exception:
                    exon_num = None

            # 2) Fallback: infer exon by genomic proximity using exons_list (if available)
            if exon_num is None and exons_list:
                for e in exons_list:
                    try:
                        en = e.get("exon_number")
                        gs = e.get("genomic_start")
                        ge = e.get("genomic_end")
                        if en is None or gs is None or ge is None:
                            continue
                        # donor site (intronic after exon end)
                        if pos > ge and pos <= ge + 10:
                            exon_num = int(en)
                            break
                        # acceptor site (intronic before exon start)
                        if pos < gs and pos >= gs - 10:
                            exon_num = int(en)
                            break
                        # inside exon
                        if gs <= pos <= ge:
                            exon_num = int(en)
                            break
                    except Exception:
                        continue

            # 3) If this is a splice consequence and we inferred exon_num, check last/penultimate exon
            cons_lower = (getattr(vr, "consequence", "") or "").lower()
            if any(tok in cons_lower for tok in ("splice_acceptor", "splice_donor", "splice_region")):
                try:
                    last_exon_number = int(last_exon.get("exon_number")) if last_exon and last_exon.get("exon_number") is not None else None
                    if exon_num is not None and last_exon_number is not None:
                        # treat splice-sites in last exon and penultimate exon as escaping NMD
                        if exon_num == last_exon_number or exon_num == (last_exon_number - 1):
                            tx_result["predicts_nmd"] = False
                            tx_result["details"]["reason"] = f"splice_site_in_last_or_penultimate_exon (exon={exon_num})"
                            per_tx_results.append(tx_result)
                            continue
                except Exception:
                    # if specialist logic fails, fall back to general logic below
                    pass

            # find penultimate CDS exon coordinates
            pen_ex = None
            last_num = None
            if last_exon.get("exon_number"):
                last_num = int(last_exon.get("exon_number"))
            if last_num and exons_list:
                for e in exons_list:
                    en = e.get("exon_number")
                    if en is not None and int(en) == last_num - 1:
                        pen_ex = e
                        break
            if not pen_ex and exons_list and len(exons_list) >= 2:
                pen_ex = exons_list[-2]
            if not pen_ex or pen_ex.get("genomic_start") is None or pen_ex.get("genomic_end") is None:
                tx_result["predicts_nmd"] = None
                tx_result["details"]["reason"] = "no_penultimate_coords"
                per_tx_results.append(tx_result)
                continue

            pen_start = int(pen_ex.get("genomic_start"))
            pen_end = int(pen_ex.get("genomic_end"))

            # compute terminal window depending on strand
            strand_chr = str(strand) if strand else "+"
            if strand_chr == "-" or strand_chr.lower().startswith("minus"):
                region_start = pen_start
                region_end = min(pen_end, pen_start + penultimate_nt - 1)
            else:
                region_start = max(pen_start, pen_end - penultimate_nt + 1)
                region_end = pen_end

            # If variant is inside terminal window -> escapes NMD
            if region_start <= pos <= region_end:
                tx_result["predicts_nmd"] = False
                tx_result["details"]["reason"] = f"variant_in_penultimate_terminal_{penultimate_nt}nt"
                per_tx_results.append(tx_result)
                continue

            # Else variant upstream of terminal window and not in last exon => predicts NMD
            tx_result["predicts_nmd"] = True
            tx_result["details"]["reason"] = "variant_upstream_of_penultimate_terminal_window"
            per_tx_results.append(tx_result)
        except Exception as e:
            tx_result["predicts_nmd"] = None
            tx_result["details"]["reason"] = f"error:{str(e)}"
            per_tx_results.append(tx_result)
            continue

    # Combine results conservatively
    details["per_transcript"] = per_tx_results
    found_true = any(r.get("predicts_nmd") is True for r in per_tx_results)
    found_false = any(r.get("predicts_nmd") is False for r in per_tx_results)
    found_unknown = any(r.get("predicts_nmd") is None for r in per_tx_results)

    # Rules:
    if found_false:
        details["reason"] = "transcript_indicates_escape"
        return False, details
    if per_tx_results and all(r.get("predicts_nmd") is True for r in per_tx_results):
        details["reason"] = "all_transcripts_predict_nmd"
        return True, details
    if found_true and found_unknown:
        details["reason"] = "mixed_or_unknown_transcripts_conservative_false"
        return False, details
    details["reason"] = "no_transcript_predicted_nmd"
    return False, details

# ---------------------------
# ACMG Engine
# ---------------------------
class ACMGEngine:
    """
    Core engine for applying ACMG/AMP guidelines to variant classification.
    
    This class implements the decision logic for assigning ACMG criteria (PVS1, PS1, PM1, etc.)
    based on variant annotations, database evidence, and gene-specific rules. It handles:
      - PVS1 strength modulation based on NMD prediction.
      - PS1/PM5 detection via ClinVar protein-level matches.
      - PM3 scoring for recessive genes with phasing-aware batching.
      - Computational evidence (PP3/BP4) integration with safeguards against double-counting.
      - Final classification using the Tavtigian et al. (2018) points-based model.
    
    Attributes:
        dbm (DatabaseManager): Database manager for evidence lookup.
        gene_rules (Dict[str, Dict]): ACMG secondary findings gene rules.
        exons_df (pd.DataFrame): Exon structure data for NMD prediction.
        spliceai_thr (float): SpliceAI threshold for PP3/BP4.
        require_mane (bool): Whether to require MANE transcripts for reporting.
    """
    def __init__(self, dbm: DatabaseManager, gene_rules: Dict[str, Dict[str,Any]],
                 exons_df: Optional[Any] = None, spliceai_thr: float = 0.20,
                 require_mane: bool = False,
                 nmd_map: Optional[Dict[str, Any]] = None,
                 exons_map: Optional[Any] = None,
                 tx_map: Optional[Dict[str, Any]] = None):
        self.dbm = dbm
        self.gene_rules = {normalize_gene_name(k): v for k,v in (gene_rules or {}).items()}
        self.exons_df = exons_df
        self.spliceai_thr = spliceai_thr
        self.require_mane = require_mane
        self.nmd_map = nmd_map or {}
        if exons_map is not None:
            self.exons_map = exons_map
        elif exons_df is not None:
            self.exons_map = exons_df
        else:
            self.exons_map = {}

    def _get_gene_rule(self, gene: str) -> Dict[str,Any]:
        return self.gene_rules.get((gene or "").upper(), {"reportable": False, "special_flags": {}, "raw_row": {}})

    # def _check_upstream_pathogenic_variants(self, vr: VariantRecord, upstream_distance: int = 1000) -> bool:
    #     if not hasattr(self, 'dbm') or self.dbm is None:
    #         logger.debug("No DatabaseManager available for upstream variant check")
    #         return False
    #     if self.dbm.clinvar_vcf is None:
    #         logger.debug("ClinVar VCF not available for upstream variant check")
    #         return False
    #     try:
    #         upstream_start = max(0, vr.pos - upstream_distance)
    #         upstream_end = vr.pos - 1
    #         if upstream_start >= upstream_end:
    #             logger.debug(f"Invalid upstream region for {vr.chrom}:{vr.pos}")
    #             return False
    #         gene_norm = normalize_gene_name(vr.gene)
    #         pathogenic_found = False
    #         pathogenic_details = []
    #         chrom_variants = [vr.chrom]
    #         if vr.chrom.startswith("chr"):
    #             chrom_variants.append(vr.chrom.replace("chr", ""))
    #         else:
    #             chrom_variants.append(f"chr{vr.chrom}")
    #         for chrom_query in chrom_variants:
    #             try:
    #                 for rec in self.dbm.clinvar_vcf.fetch(chrom_query, upstream_start, upstream_end):
    #                     gene_info = rec.info.get('GENEINFO', '')
    #                     rec_genes = []
    #                     if gene_info:
    #                         for gene_pair in str(gene_info).split('|'):
    #                             if ':' in gene_pair:
    #                                 clinvar_gene = gene_pair.split(':')[0]
    #                                 rec_genes.append(normalize_gene_name(clinvar_gene))
    #                     if not rec_genes or gene_norm not in rec_genes:
    #                         continue
    #                     hgvsp_info = rec.info.get('HGVS_Protein', '')
    #                     if not hgvsp_info:
    #                         continue
    #                     clnsig = rec.info.get('CLNSIG')
    #                     clnsig_str = ";".join(clnsig) if isinstance(clnsig, (list, tuple)) else str(clnsig)
    #                     rev = rec.info.get('CLNREVSTAT') or ''
    #                     rev_str = str(rev).lower()
    #                     stars = 0
    #                     if 'practice_guideline' in rev_str: stars = 4
    #                     elif 'expert_panel' in rev_str: stars = 3
    #                     elif 'criteria_provided' in rev_str and 'multiple_submitters' in rev_str and 'no_conflict' in rev_str: stars = 2
    #                     elif 'criteria_provided' in rev_str: stars = 1
    #                     if stars >= 2:
    #                         if ('pathogenic' in clnsig_str or 'likely pathogenic' in clnsig_str) and 'benign' not in clnsig_str:
    #                             pathogenic_found = True
    #                             variant_info = {
    #                                 'position': rec.pos,
    #                                 'significance': clnsig_str,
    #                                 'stars': stars,
    #                                 'hgvsp': hgvsp_info[0] if isinstance(hgvsp_info, (list, tuple)) else hgvsp_info
    #                             }
    #                             pathogenic_details.append(variant_info)
    #                             logger.debug(f"Found upstream pathogenic variant at {rec.chrom}:{rec.pos}: {clnsig_str} ({stars} stars)")
    #             except ValueError:
    #                 continue
    #             except Exception as e:
    #                 logger.debug(f"Error fetching upstream variants for {chrom_query}: {e}")
    #                 continue
    #         if pathogenic_found:
    #             vr.manual_reasons.append(f"Found {len(pathogenic_details)} pathogenic variant(s) upstream:")
    #             for i, var in enumerate(pathogenic_details, 1):
    #                 vr.manual_reasons.append(f"  {i}. {var['position']}: {var['significance']} ({var['stars']} stars) - {var['hgvsp']}")
    #         return pathogenic_found
    #     except Exception as e:
    #         logger.error(f"Error checking upstream pathogenic variants: {e}")
    #         return False

    def _find_nearest_inframe_start_codon(self, vr: VariantRecord, fasta_handle=None) -> Tuple[Optional[int], Optional[str]]:
        if fasta_handle is None or not HAVE_PYSAM:
            logger.debug("No FASTA handle available for start codon search")
            return None, None
        try:
            search_distance = 300
            search_end = vr.pos + search_distance
            sequence = fasta_handle.fetch(vr.chrom, vr.pos, search_end).upper()
            if len(sequence) < 3:
                return None, None
            for frame in range(3):
                for i in range(frame, len(sequence) - 2, 3):
                    codon = sequence[i:i+3]
                    if codon == "ATG":
                        atg_position = vr.pos + i
                        position_diff = atg_position - vr.pos
                        if position_diff % 3 == 0:
                            context_start = max(0, i - 9)
                            context_end = min(len(sequence), i + 12)
                            context = sequence[context_start:context_end]
                            logger.debug(f"Found in-frame ATG at position {atg_position} (frame {frame})")
                            logger.debug(f"Context: ...{context}...")
                            return atg_position, context
            logger.debug(f"No in-frame ATG found within {search_distance}bp downstream")
            return None, None
        except Exception as e:
            logger.error(f"Error finding nearest start codon: {e}")
            return None, None

    def _check_protein_truncation_percentage(self, vr: VariantRecord) -> Tuple[bool, Optional[float]]:
        """
        Check what percentage of protein is removed by truncation.
        
        Returns:
            Tuple[bool, Optional[float]]: (is_more_than_10_percent, percentage_removed)
        """
        if not vr.hgvsp:
            return False, None
        
        # Try to extract truncation position from HGVSp
        pos_match = re.search(r'p\.[A-Za-z*]{1,3}(\d+)', vr.hgvsp)
        if not pos_match:
            pos_match = re.search(r'p\.(\d+)', vr.hgvsp)
            if not pos_match:
                return False, None
        
        try:
            trunc_pos = int(pos_match.group(1))
            
            # Comprehensive protein lengths for ACMG SF v3.3 genes
            protein_lengths = {
                # Cancer genes
                "BRCA1": 1863, "BRCA2": 3418, "TP53": 393, "MLH1": 756,
                "MSH2": 934, "MSH6": 1360, "PMS2": 862, "APC": 2843,
                "RB1": 928, "VHL": 213, "MEN1": 610, "RET": 1114,
                "PTEN": 403, "NF2": 595, "MAX": 160, "STK11": 433,
                "SMAD4": 552, "PALB2": 1186, "MUTYH": 546, "SDHB": 280,
                "SDHC": 169, "SDHD": 159, "SDHAF2": 166, "TMEM127": 238,
                "TSC1": 1164, "TSC2": 1807, "WT1": 449,
                
                # Cardiovascular genes
                "MYH7": 1935, "MYBPC3": 1274, "TNNT2": 298, "TNNI3": 210,
                "KCNH2": 1159, "KCNQ1": 676, "SCN5A": 2016, "RYR2": 4967,
                "DSP": 2871, "PKP2": 881, "DSG2": 1118, "DSC2": 901,
                "LMNA": 664, "MYH11": 1936, "ACTA2": 377, "SMAD3": 425,
                "TGFBR1": 503, "TGFBR2": 592, "RYR1": 5037, "CACNA1S": 1873,
                "FLNC": 2725, "TTN": 34350, "PLN": 52, "DES": 470,
                "ACTN2": 894, "TPM1": 284, "MYL2": 166, "MYL3": 195,
                "ACTC1": 377, "BAG3": 575, "RBM20": 1227, "CASQ2": 399,
                "CALM1": 149, "CALM2": 149, "CALM3": 149, "TRDN": 729,
                "TTR": 147, "PRKAG2": 569, "TMEM43": 400,
                
                # Connective tissue / aortic
                "FBN1": 2871, "COL3A1": 1466, "ACTA2": 377, "SMAD3": 425,
                "TGFBR1": 503, "TGFBR2": 592, "MYH11": 1936,
                
                # Metabolic genes
                "ATP7B": 1465, "CFTR": 1480, "BTD": 543, "GAA": 952,
                "HFE": 343, "OTC": 354, "PMM2": 246, "CYP27A1": 531,
                "RPE65": 533, "ALDOB": 364, "PAH": 452, "SERPINA1": 418,
                "PCSK9": 692, "HNF1A": 631,
                
                # X-linked genes
                "ABCD1": 745, "GLA": 429, "DMD": 3685,
                
                # Other
                "ENG": 658, "ACVRL1": 503,
            }
            
            gene_key = vr.gene.upper()
            if gene_key in protein_lengths:
                protein_len = protein_lengths[gene_key]
                # Calculate percentage of protein RETAINED (from truncation point to end)
                # Or percentage REMOVED? For PVS1, we care about how much is removed.
                # If truncation removes >10% of protein -> PVS1_Strong
                # If removes <10% -> PVS1_Moderate
                
                # Percentage of protein RETAINED after truncation
                retained_percent = ((protein_len - trunc_pos + 1) / protein_len) * 100
                # Percentage REMOVED
                removed_percent = 100 - retained_percent
                
                logger.debug(f"Truncation check for {vr.gene}: trunc_pos={trunc_pos}, protein_len={protein_len}, "
                            f"retained={retained_percent:.1f}%, removed={removed_percent:.1f}%")
                
                # Return True if >10% of protein is removed
                return removed_percent > 10, removed_percent
            else:
                # Unknown gene - conservative: assume >10% if truncation occurs before position 100
                # (most proteins are >100 AA, so truncation at pos <100 likely removes >10%)
                is_gt_10_percent = trunc_pos < 100
                logger.debug(f"Unknown gene {vr.gene}, trunc_pos={trunc_pos}, is_gt_10_percent={is_gt_10_percent}")
                return is_gt_10_percent, None
                
        except Exception as e:
            logger.debug(f"Error checking protein truncation percentage for {vr.gene}: {e}")
            # Conservative: return False (treat as small truncation) to avoid over-calling
            return False, None

    def _is_critical_domain(self, vr: VariantRecord) -> bool:
        """
        Check if the truncated/altered region is critical to protein function.
        Based on ClinGen SVI guidance.
        
        Critical domains include:
        - DNA-binding domains
        - Active sites
        - Transmembrane domains
        - Known functional domains from literature
        - For ACMG SF v3.3 genes, many are LOF-sensitive
        """
        # Get amino acid position of truncation
        aa_pos = vr.aa_pos
        if aa_pos is None:
            return False
        
        # Critical domains for ACMG SF v3.3 genes (based on literature)
        # Format: {gene: {domain_name: (start, end)}}
        critical_domains = {
            # Cancer genes - highly LOF-sensitive
            "BRCA1": {"RING_finger": (1, 100), "BRCT": (1640, 1863)},
            "BRCA2": {"DNA_binding": (2450, 3200)},
            "TP53": {"DNA_binding": (100, 300), "TAD": (1, 40)},
            "MLH1": {"MutL": (1, 756)},
            "MSH2": {"MutS": (1, 934)},
            "MSH6": {"MutS": (1, 1360)},
            "PMS2": {"MutL": (1, 862)},
            "APC": {"ARM": (1, 2843)},
            "RB1": {"pocket": (373, 870)},
            "PTEN": {"phosphatase": (1, 403)},
            "VHL": {"elongin_binding": (1, 213)},
            "MEN1": {"menin": (1, 610)},
            "RET": {"tyrosine_kinase": (1, 1114)},
            "NF2": {"FERM": (1, 595)},
            "MAX": {"bHLHZ": (1, 160)},
            "STK11": {"kinase": (1, 433)},
            "SMAD4": {"MH1": (1, 150), "MH2": (300, 552)},
            "PALB2": {"WD40": (1, 1186)},
            "MUTYH": {"HhH": (1, 546)},
            
            # Cardiovascular genes
            "MYH7": {"motor_domain": (1, 800), "tail": (801, 1935)},
            "MYBPC3": {"Ig_C2": (1, 1274)},
            "KCNH2": {"pore": (550, 650), "voltage_sensor": (400, 500)},
            "KCNQ1": {"pore": (300, 350), "voltage_sensor": (150, 250)},
            "SCN5A": {"pore": (1, 2016), "voltage_sensor": (1, 2016)},
            "RYR2": {"ryanodine_receptor": (1, 4967)},
            "DSP": {"plakin": (1, 2871)},
            "PKP2": {"armadillo": (1, 881)},
            "DSG2": {"cadherin": (1, 1118)},
            "DSC2": {"cadherin": (1, 901)},
            "LMNA": {"rod": (1, 664)},
            "TNNT2": {"tropomyosin_binding": (1, 298)},
            "TNNI3": {"tropomyosin_binding": (1, 210)},
            "MYL2": {"EF_hand": (1, 166)},
            "MYL3": {"EF_hand": (1, 195)},
            "TPM1": {"tropomyosin": (1, 284)},
            "ACTC1": {"actin": (1, 377)},
            "ACTN2": {"spectrin": (1, 894)},
            "BAG3": {"BAG": (1, 575)},
            "DES": {"desmin": (1, 470)},
            "FLNC": {"Ig_like": (1, 2725)},
            "RBM20": {"RRM": (1, 1227)},
            "TTN": {"I-band": (1, 34350), "A-band": (1, 34350)},
            "PLN": {"transmembrane": (1, 52)},
            "TTR": {"thyroxine_binding": (1, 147)},
            "CASQ2": {"calsequestrin": (1, 399)},
            "CALM1": {"EF_hand": (1, 149)},
            "CALM2": {"EF_hand": (1, 149)},
            "CALM3": {"EF_hand": (1, 149)},
            "TRDN": {"triadin": (1, 729)},
            
            # Connective tissue / aortic genes
            "FBN1": {"cb_EGF": (1, 2871), "TGFB_binding": (1, 2871)},
            "COL3A1": {"collagen": (1, 1466)},
            "ACTA2": {"actin": (1, 377)},
            "SMAD3": {"MH1": (1, 145), "MH2": (250, 425)},
            "TGFBR1": {"kinase": (1, 503)},
            "TGFBR2": {"kinase": (1, 592)},
            "MYH11": {"motor": (1, 800)},
            
            # Metabolic genes
            "ATP7B": {"copper_binding": (1, 1465)},
            "CFTR": {"NBD": (1, 1480), "R_domain": (1, 1480)},
            "BTD": {"biotin": (1, 543)},
            "GAA": {"glycosyl_hydrolase": (1, 952)},
            "HFE": {"MHC_class_I": (1, 343)},
            "OTC": {"ornithine": (1, 354)},
            "PMM2": {"phosphomannomutase": (1, 246)},
            "CYP27A1": {"cytochrome_P450": (1, 531)},
            "RPE65": {"RPE65": (1, 533)},
            
            # X-linked genes
            "ABCD1": {"ABCD1": (1, 745)},
            "GLA": {"alpha_galactosidase": (1, 429)},
            "DMD": {"spectrin": (1, 3685)},
            
            # Other
            "PRKAG2": {"CBS": (1, 569)},
            "TMEM43": {"transmembrane": (1, 400)},
            "SDHB": {"succinate_dehydrogenase": (1, 280)},
            "SDHC": {"succinate_dehydrogenase": (1, 169)},
            "SDHD": {"succinate_dehydrogenase": (1, 159)},
            "SDHAF2": {"SDHAF2": (1, 166)},
            "ENG": {"ZP": (1, 658)},
            "ACVRL1": {"kinase": (1, 503)},
            "HNF1A": {"homeobox": (1, 631)},
            "TMEM127": {"transmembrane": (1, 238)},
            "PCSK9": {"proprotein_convertase": (1, 692)},
        }
        
        gene_domains = critical_domains.get(vr.gene.upper(), {})
        
        # Check if truncation occurs within or before any critical domain
        for domain_name, (domain_start, domain_end) in gene_domains.items():
            if aa_pos <= domain_end:
                return True
        
        # For highly LOF-sensitive genes without specific domain info,
        # any truncation is considered critical
        lof_sensitive_genes = {
            "TP53", "BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "PMS2", "APC", 
            "NF1", "NF2", "RB1", "VHL", "PTEN", "MEN1", "RET", "STK11", 
            "SMAD4", "PALB2", "MUTYH", "DMD", "CFTR", "ATP7B"
        }
        if vr.gene.upper() in lof_sensitive_genes:
            return True
        
        return False

    def _check_last_exon_status(self, vr: VariantRecord) -> bool:
        used_csv = False
        if self.exons_df is not None and not self.exons_df.empty:
            tx_id = vr.transcript
            if tx_id:
                if '.' in tx_id: tx_base = tx_id.split('.')[0]
                else: tx_base = tx_id
                if hasattr(self.exons_df, '_transcript_last_exons'):
                    last_exon_info = None
                    if tx_id in self.exons_df._transcript_last_exons:
                        last_exon_info = self.exons_df._transcript_last_exons[tx_id]
                    elif tx_base in self.exons_df._transcript_last_exons:
                        last_exon_info = self.exons_df._transcript_last_exons[tx_base]
                    if last_exon_info:
                        used_csv = True
                        variant_pos = int(vr.pos)
                        if last_exon_info['start'] <= variant_pos <= last_exon_info['end']:
                            return True
                        return False
        if not used_csv and vr.exon:
            try:
                if "/" in vr.exon:
                    parts = vr.exon.split("/")
                    curr = int(parts[0])
                    total = int(parts[1])
                    if curr == total and total > 0:
                        return True
            except Exception:
                pass
        return False

    def evaluate_variant(self, vr: VariantRecord, fasta_handle=None):
        """
        Evaluate ACMG criteria for a single VariantRecord (vr) and produce:
        - assigned: list of applied criteria (strings, e.g. "PVS1_VeryStrong")
        - pts: dict mapping criterion->points
        - auto_class: automated class ("Pathogenic"/"Likely pathogenic"/"VUS"/"Likely benign"/"Benign")
        - manual_flag (always False here; triage decides manual later)
        - manual_reasons: list of reasons to show to curators

        Also fills:
        - vr.criteria_explanations (dict): human-readable explanation per base criterion (PVS1, PS1, PM2, ...)
        - vr.ps1_matches / vr.pm5_matches for downstream reporting
        """
        vr.criteria_assigned = []
        vr.criteria_points = {}
        vr.criteria_explanations = {}
        vr.ps1_matches = []
        vr.pm5_matches = []
        vr.criteria_assigned_ext = []
        vr.criteria_points_ext = {}
        vr.total_ext_points = 0
        vr.extended_class = "VUS"
        vr._pp5ext_explanation = ""
        vr._clinvar_phenotype_match = False
        vr._hgmd_phenotype_match = False
        vr.automated_class = "VUS"
        vr.manual_review = False
        vr.manual_reasons = []
        vr._acmg_disease = self._get_gene_rule(vr.gene).get("disease", "")
        logger.debug("Starting evaluate_variant for %s:%s %s>%s (sample=%s)", getattr(vr,"chrom",None), getattr(vr,"pos",None), getattr(vr,"ref",None), getattr(vr,"alt",None), getattr(vr,"sample",None))
        
        pts: Dict[str,int] = {}
        assigned: List[str] = []
        auto_class = "VUS"
        manual_reasons: List[str] = []
        # Ensure explanation containers exist
        if not hasattr(vr, "criteria_explanations") or vr.criteria_explanations is None:
            vr.criteria_explanations = {}
        if not hasattr(vr, "ps1_matches"):
            vr.ps1_matches = []
        if not hasattr(vr, "pm5_matches"):
            vr.pm5_matches = []

        cons = (vr.consequence or "").lower()
        af = vr.gnomad_af if vr.gnomad_af is not None and vr.gnomad_af >= 0 else 0.0

        # -------------------------
        # 0. BA1 (stand-alone benign)
        # -------------------------
        is_ba1_exception = False
        exceptions = THRESH.get("BA1_EXCEPTIONS", {})
        try:
            if vr.gene in exceptions:
                for exc in exceptions[vr.gene]:
                    if exc and exc in (vr.hgvsp or ""):
                        is_ba1_exception = True
                        break
        except Exception:
            is_ba1_exception = False

        if af >= THRESH.get("BA1_AF", 0.05) and not is_ba1_exception:
            assigned.append("BA1")
            pts["BA1"] = POINTS_MAP.get("BA1", -8)
            vr.criteria_explanations["BA1"] = f"BA1: gnomAD max AF={af:.6g} >= BA1 threshold {THRESH.get('BA1_AF')}; not in exceptions"
            total = sum(pts.values())
            return assigned, pts, "Benign", False, manual_reasons

        # -------------------------
        # 1. BS1 (strong benign by freq)
        # -------------------------
        is_bs1_exception = False
        try:
            if vr.gene in exceptions:
                for exc in exceptions[vr.gene]:
                    if exc and exc in (vr.hgvsp or ""):
                        is_bs1_exception = True
                        break
        except Exception:
            is_bs1_exception = False

        if af >= THRESH.get("BS1_AF", 0.01) and not is_bs1_exception:
            assigned.append("BS1")
            pts["BS1"] = POINTS_MAP.get("BS1", -4)
            vr.criteria_explanations["BS1"] = f"BS1: gnomAD max AF={af:.6g} >= BS1 threshold {THRESH.get('BS1_AF')}; not in exceptions"

        # -------------------------
        # 2. PVS1 Decision tree
        # -------------------------
        special = self._get_gene_rule(vr.gene).get("special_flags", {}) if self.gene_rules else {}
        # support both older/short and VEP-style consequence tokens
        is_lof_molecular = any(x in cons for x in ["stop_gained", "frameshift", "frameshift_variant", "start_lost", "stop_lost"])
        # Determine canonical splice only for true canonical positions (+/-1-2 nt) using helper
        # Use nmd_map for accurate exon boundary detection (RefSeq transcripts only)
        is_canonical_splice = _is_canonical_splice_site(
            vr, 
            nmd_map=getattr(self, "nmd_map", {})
        )
        logger.debug(f"PVS1 check for {vr.gene}:{vr.hgvsp} - Consequence: {cons}")
        logger.debug(f"  is_lof_molecular: {is_lof_molecular}")
        logger.debug(f"  is_canonical_splice: {is_canonical_splice}")

        is_splice_acceptor_donor = any(x in cons for x in ["splice_acceptor_variant", "splice_donor_variant"])
        is_canonical_splice = _is_canonical_splice_site(vr, nmd_map=getattr(self, "nmd_map", {})) or is_splice_acceptor_donor

        if is_lof_molecular or is_canonical_splice:
            # Skip if gene explicitly marked LOF not reportable
            if not special.get("lof_not_reportable", False):
                triggers_nmd = False
                nmd_details = {}
                in_last_exon = None
                pred_error = None

                try:
                    # respect any internal NMD annotation first (case-insensitive)
                    internal_nmd_flag = str(getattr(vr, "nmd", "") or "").strip()
                    internal_nmd_low = internal_nmd_flag.lower()
                    if internal_nmd_flag and ("escape" in internal_nmd_low or "last_exon" in internal_nmd_low):
                        # explicit internal prediction indicates escape -> do not trigger NMD
                        triggers_nmd = False
                        nmd_details = {"reason": f"internal_nmd:{internal_nmd_flag}"}
                        logger.debug("Internal NMD predictor indicates escape -> skipping transcript-aware NMD for PVS1")
                    else:
                        # call transcript-aware predictor (conservative)
                        triggers_nmd, nmd_details = variant_triggers_nmd(
                            vr,
                            getattr(self, "nmd_map", {}),
                            getattr(self, "exons_map", {}),
                            getattr(self, "tx_map", {}),
                            penultimate_nt=50
                        )
                        logger.debug(f"  NMD prediction: triggers_nmd={triggers_nmd}")
                        logger.debug(f"  NMD details: {nmd_details}")
                except Exception as _e:
                    pred_error = str(_e)
                    try:
                        in_last_exon = self._check_last_exon_status(vr)
                    except Exception:
                        in_last_exon = False
                    triggers_nmd = not bool(in_last_exon)
                    nmd_details = {"reason": f"nmd_pred_error:{pred_error}", "fallback_in_last_exon": in_last_exon}
                    logger.debug(f"  NMD prediction failed, fallback to last-exon check: in_last_exon={in_last_exon}, triggers_nmd={triggers_nmd}; error: {pred_error}")

                # If internal predictor explicitly indicates escape at any time, honour it
                try:
                    internal_nmd_flag = str(getattr(vr, "nmd", "") or "").strip()
                    if internal_nmd_flag:
                        if "escape" in internal_nmd_flag.lower() or "last_exon" in internal_nmd_flag.lower():
                            triggers_nmd = False
                            nmd_details.setdefault("override", []).append(f"internal_nmd_override:{internal_nmd_flag}")
                            vr.criteria_explanations["PVS1_NMD_internal_override"] = internal_nmd_flag
                            logger.debug(f"Overriding NMD -> internal predictor indicates escape ({internal_nmd_flag})")
                except Exception:
                    # do not fail the pipeline for this non-critical step
                    pass

                # Persist a short explanation for audit
                vr.criteria_explanations["PVS1_NMD"] = nmd_details.get("reason", json.dumps(nmd_details))
                # Start-loss
                if "start_lost" in cons or "start_lost" in (vr.consequence or "").lower():
                    # Check for downstream in-frame ATG
                    atg_pos, context = self._find_nearest_inframe_start_codon(vr, fasta_handle=fasta_handle)
                    
                    if atg_pos is None:
                        # No in-frame downstream ATG -> PVS1_Strong (4 points)
                        assigned.append("PVS1_Strong")
                        pts["PVS1_Strong"] = POINTS_MAP.get("PVS1_Strong", 4)
                        expl = "PVS1_Strong: start-loss variant with no in-frame downstream ATG found."
                    else:
                        # Has in-frame downstream ATG -> PVS1_Moderate (4 points)
                        assigned.append("PVS1_Moderate")
                        pts["PVS1_Moderate"] = POINTS_MAP.get("PVS1_Moderate", 2)
                        expl = f"PVS1_Moderate: start-loss variant. Nearest in-frame ATG at pos {atg_pos} (maintains frame)."
                    # Check upstream pathogenic variants (strengthen caution)
                    # upstream_found = False
                    # try:
                    #     upstream_found = self._check_upstream_pathogenic_variants(vr)
                    # except Exception:
                    #     upstream_found = False
                    # if upstream_found:
                    #     expl += " Upstream pathogenic variants found -> manual inspection recommended."
                    vr.criteria_explanations["PVS1"] = expl

                elif triggers_nmd:
                    # Triggers NMD -> PVS1_VeryStrong
                    assigned.append("PVS1_VeryStrong")
                    pts["PVS1_VeryStrong"] = POINTS_MAP.get("PVS1_VeryStrong", 8)
                    vr.criteria_explanations["PVS1"] = ("PVS1_VeryStrong: Variant is a predicted loss-of-function event and is predicted to trigger NMD")
                else:
                    # Escapes NMD => DO NOT assign PVS1. For truncating variants in last exon, use PM4 instead.
                    # For canonical splice variants in last exon, treat conservatively as PM4_Moderate
                    if is_canonical_splice:
                        assigned.append("PM4_Moderate")
                        pts["PM4_Moderate"] = POINTS_MAP.get("PM4_Moderate", 2)
                        vr.criteria_explanations["PM4"] = ("PM4_Moderate: canonical splice variant located in last exon / predicted to escape NMD; treated as protein-altering in last exon")
                    else:
                        # Stop/Frameshift in last exon - use PVS1 per flowchart
                        is_gt_10_percent, removed_pct = self._check_protein_truncation_percentage(vr)
                        is_critical = self._is_critical_domain(vr)

                        if is_gt_10_percent or is_critical:
                            # PVS1_Strong (4 points) for >10% removal OR critical domain
                            assigned.append("PVS1_Strong")
                            pts["PVS1_Strong"] = POINTS_MAP.get("PVS1_Strong", 4)
                            if is_gt_10_percent and is_critical:
                                vr.criteria_explanations["PVS1"] = f"PVS1_Strong: truncating variant in last exon removing >10% of protein ({removed_pct}%) in critical domain"
                            elif is_gt_10_percent:
                                vr.criteria_explanations["PVS1"] = f"PVS1_Strong: truncating variant in last exon removing >10% of protein ({removed_pct}%)"
                            else:
                                vr.criteria_explanations["PVS1"] = "PVS1_Strong: truncating variant in last exon affecting critical protein domain"
                        else:
                            # PVS1_Moderate (2 points) for <10% removal outside critical domain
                            assigned.append("PVS1_Moderate")
                            pts["PVS1_Moderate"] = POINTS_MAP.get("PVS1_Moderate", 2)
                            if removed_pct is not None:
                                vr.criteria_explanations["PVS1"] = f"PVS1_Moderate: truncating variant in last exon removing {removed_pct}% of protein (<10%) outside critical domain"
                            else:
                                vr.criteria_explanations["PVS1"] = "PVS1_Moderate: truncating variant in last exon removing <10% of protein outside critical domain"

        # -------------------------
        # PM4 (protein length changes)
        # -------------------------
        if "stop_lost" in cons or "stop_retained_variant" in cons:
            assigned.append("PM4_Strong")
            pts["PM4_Strong"] = POINTS_MAP.get("PM4_Strong", 4)
            vr.criteria_explanations["PM4"] = "PM4_Strong: stop-loss predicted to alter C-terminal sequence significantly"

        if "inframe" in cons or "protein_altering_variant" in cons or "inframe_deletion" in cons or "inframe_insertion" in cons:
            # check size and critical domains
            small_noncrit = self._check_inframe_small_noncritical(vr)
            if small_noncrit:
                assigned.append("PM4_Supporting")
                pts["PM4_Supporting"] = POINTS_MAP.get("PM4_Supporting", 1)
                vr.criteria_explanations["PM4"] = "PM4_Supporting: small in-frame indel (<10 AA) outside annotated critical domains"
            else:
                assigned.append("PM4_Moderate")
                pts["PM4_Moderate"] = POINTS_MAP.get("PM4_Moderate", 2)
                vr.criteria_explanations["PM4"] = "PM4_Moderate: in-frame indel or protein length change of moderate size or in important region"

        # -------------------------
        # PS1 & PM5 (missense / protein changes)
        # -------------------------
        pos_aa = parse_protein_pos(vr.hgvsp)
        vr.aa_pos = pos_aa
        rule = self._get_gene_rule(vr.gene) if self.gene_rules else {}
        acmg_disease = str(rule.get("disease","") or "").strip()
        if "missense" in cons and pos_aa:
            
            # Helper function to determine PS1 strength based on stars and classification
            def get_ps1_strength(stars: int, clnsig: str) -> Tuple[Optional[str], Optional[int]]:
                """
                Determine PS1 strength based on stars and classification.
                Returns (criteria_name, points).
                
                Rules:
                - 2+ stars + Pathogenic -> PS1_Strong (4 points)
                - 2+ stars + Likely pathogenic -> PS1_Moderate (2 points)
                - 1 star + Pathogenic -> PS1_Supporting (1 point)
                - 1 star + Likely pathogenic -> PS1_Supporting (1 point)
                """
                clnsig_lower = clnsig.lower()
                is_pathogenic = "pathogenic" in clnsig_lower and "likely" not in clnsig_lower
                is_likely_pathogenic = "likely pathogenic" in clnsig_lower or ("likely" in clnsig_lower and "pathogen" in clnsig_lower)
                
                if stars >= 2:
                    if is_pathogenic:
                        return ("PS1_Strong", POINTS_MAP.get("PS1_Strong", 4))
                    elif is_likely_pathogenic:
                        return ("PS1_Moderate", POINTS_MAP.get("PS1_Moderate", 2))
                elif stars >= 1:
                    if is_pathogenic or is_likely_pathogenic:
                        return ("PS1_Supporting", POINTS_MAP.get("PS1_Supporting", 1))
                
                return (None, None)
            
            # PS1
            matches = []
            try:
                matches = self.dbm.find_protein_variant_by_hgvsp(
                    vr.gene, vr.hgvsp, pos_aa, vr.transcript,
                    current_chrom=vr.chrom, current_pos=vr.pos, 
                    current_ref=vr.ref, current_alt=vr.alt,
                    acmg_disease=acmg_disease
                )
                for m in matches:
                    logger.info(f"  PS1 match: {m}")
            except Exception as e:
                logger.error(f"PS1 lookup failed: {e}")
                import traceback
                traceback.print_exc()

            ps1_assigned = False
            for m in matches:
                clnsig = ""
                ps1_strength = None
                ps1_points = None
                info = {}
                
                try:
                    if isinstance(m, (list, tuple)):
                        clnsig = str(m[0] or "")
                        if len(m) > 1:
                            # Check if second element is a strength string (PS1_Strong, etc.)
                            if isinstance(m[1], str) and m[1].startswith("PS1_"):
                                ps1_strength = str(m[1])
                                if len(m) > 2:
                                    ps1_points = int(m[2] or 0)
                                if len(m) > 3 and isinstance(m[3], dict):
                                    info = m[3]
                            else:
                                # Old format: (clnsig, source)
                                src = str(m[1] or "")
                                if len(m) > 2 and isinstance(m[2], dict):
                                    info = m[2]
                                # Calculate strength from stars and clnsig
                                stars = int(info.get("stars", 0) or 0)
                                ps1_strength, ps1_points = get_ps1_strength(stars, clnsig)
                    elif isinstance(m, dict):
                        clnsig = str(m.get("clnsig",""))
                        src = str(m.get("source",""))
                        info = m
                        stars = int(info.get("stars", 0) or 0)
                        ps1_strength, ps1_points = get_ps1_strength(stars, clnsig)
                except Exception as e:
                    logger.debug(f"Error parsing PS1 match: {e}")
                    continue

                if ps1_strength is None or ps1_points is None:
                    continue
                
                cln_low = clnsig.lower()
                
                # CRITICAL: Skip conflict entries completely - do NOT add to ps1_matches
                if "conflict" in cln_low:
                    logger.info(f"Skipping PS1: conflict entry for {info.get('hgvsp', vr.hgvsp)}")
                    # Do NOT add to vr.ps1_matches - just skip
                    continue
                
                # Assign PS1 with appropriate strength
                assigned.append(ps1_strength)
                pts[ps1_strength] = ps1_points
                ps1_assigned = True
                vr.ps1_matches.append({
                    "hgvsp": info.get("hgvsp", vr.hgvsp),
                    "clnsig": clnsig,
                    "source": src if 'src' in locals() else "ClinVar",
                    "stars": info.get("stars", ""),
                    "ps1_strength": ps1_strength,
                    "ps1_points": ps1_points
                })
                logger.info(f"PS1 assigned: {ps1_strength} ({ps1_points} points) for {clnsig} (stars={info.get('stars', 0)})")
                break  # Only assign the strongest match

            if not ps1_assigned:
                # PM5: different AA change at same residue reported pathogenic/likely pathogenic
                try:
                    pm5_matches = self.dbm.check_pm5_variants(
                        vr.gene, pos_aa, vr.hgvsp,
                        current_chrom=vr.chrom, 
                        current_pos=vr.pos,
                        current_ref=vr.ref,
                        current_alt=vr.alt
                    )
                except Exception:
                    pm5_matches = []
                
                if pm5_matches:
                    # Check if any match is Pathogenic (not just Likely)
                    has_pathogenic = any(m.get("best") == "Pathogenic" for m in pm5_matches)
                    
                    if has_pathogenic:
                        assigned.append("PM5_Moderate")
                        pts["PM5_Moderate"] = POINTS_MAP.get("PM5_Moderate", 2)
                        vr.criteria_explanations["PM5"] = f"PM5_Moderate: other different AA at same residue reported Pathogenic in database (residue {pos_aa})"
                    else:
                        assigned.append("PM5_Supporting")
                        pts["PM5_Supporting"] = POINTS_MAP.get("PM5_Supporting", 1)
                        vr.criteria_explanations["PM5"] = f"PM5_Supporting: other different AA at same residue reported Likely pathogenic in database (residue {pos_aa})"
                    
                    # Store all matches for output
                    vr.pm5_matches.extend(pm5_matches)

            # PS1 explanation (if assigned)
            if any(c.startswith("PS1") for c in assigned):
                # describe matches
                try:
                    details = "; ".join([f"{m.get('hgvsp','?')} ({m.get('clnsig','')}) [{m.get('ps1_strength','')}]" for m in vr.ps1_matches])
                except Exception:
                    details = "Matched pathogenic protein change in DB"
                vr.criteria_explanations["PS1"] = f"PS1: exact same protein change observed in DB: {details}"

        # -------------------------
        # PM2 (population rarity) - MOI specific thresholds
        # -------------------------
        rule = self._get_gene_rule(vr.gene) if self.gene_rules else {}
        moi = str(rule.get("moi","")).upper() if rule else ""
        pm2_thr = THRESH.get("PM2_AD_AF", 0.0001)
        if "AR" in moi or "RECESSIVE" in moi:
            pm2_thr = THRESH.get("PM2_AR_AF", 0.005)
        elif "XLR" in moi:
            pm2_thr = THRESH.get("PM2_XLR_AF", 0.003)
        elif "XLD" in moi:
            pm2_thr = THRESH.get("PM2_XLD_AF", 0.0001)

        # CRITICAL: Only apply PM2 if we have valid frequency data (af >= 0, not -1)
        af = vr.gnomad_af if vr.gnomad_af is not None and vr.gnomad_af >= 0 else None

        if af is not None and af < pm2_thr:
            assigned.append("PM2_Moderate")
            pts["PM2_Moderate"] = POINTS_MAP.get("PM2_Moderate", 2)
            vr.criteria_explanations["PM2"] = f"PM2_Moderate: gnomAD max AF={af:.6g} < threshold {pm2_thr} for MOI {moi}"
        else:
            if af is None or af < 0:
                vr.criteria_explanations["PM2"] = f"PM2 not applied: no valid frequency data (AF={vr.gnomad_af})"

        # -------------------------
        # PP3 / BP4 (computational predictions)
        # -------------------------
        s_score = vr.spliceai if vr.spliceai is not None else 0.0
        r_score = vr.revel if vr.revel is not None else -1.0
        a_score = vr.alpha_missense if vr.alpha_missense is not None else -1.0
        a_category = getattr(vr, "alpha_missense_category", None)  # 'P', 'B', 'A', etc.
        sq_score = getattr(vr, "squirls_score", 0.0) 
        sq_category = getattr(vr, "squirls_category", None)

        # CADD threshold name compatibility
        # cadd_thr = THRESH.get("CADD_SUPPORTING", THRESH.get("CADD_MODERATE", 30))

        # Check for Moderate evidence (priority over Supporting)
        is_pp3_moderate = False
        if s_score is not None and s_score >= THRESH.get("SPLICEAI_MODERATE", 0.20):
            is_pp3_moderate = True
        elif r_score is not None and r_score >= THRESH.get("REVEL_MODERATE", 0.93):
            is_pp3_moderate = True
        # elif c_score := (vr.cadd if vr.cadd is not None else 0.0):
        #     if c_score >= cadd_thr:
        #         is_pp3_moderate = True

        # Check for Supporting evidence (only if Moderate not already met)
        # REVEL supporting threshold OR AlphaMissense numeric threshold OR AlphaMissense category 'P' OR SpliceAI
        is_pp3_supporting = False
        if not is_pp3_moderate:
            if r_score is not None and r_score >= THRESH.get("REVEL_SUPPORTING", 0.64):
                is_pp3_supporting = True
            elif a_score is not None and a_score >= THRESH.get("ALPHA_MISSENSE_SUPPORTING", 0.56):
                is_pp3_supporting = True
            elif a_category == "P":
                is_pp3_supporting = True
            elif s_score is not None and s_score >= THRESH.get("SPLICEAI_SUPPORTING", 0.20):
                is_pp3_supporting = True
            elif sq_score is not None and sq_score >= THRESH.get("SQUIRLS_SUPPORTING", 0.70):
                is_pp3_supporting = True
            elif sq_category in ("D", "P"):  # D = Deleterious, P = Pathogenic
                is_pp3_supporting = True

        # Do NOT apply PP3 if PVS1 was applied or the variant is an obvious LOF/canonical splice consequence.
        pvs1_present = any(x.startswith("PVS1") for x in assigned)
        cons_lof_like = any(x in cons for x in ["stop_gained", "frameshift", "frameshift_variant", "start_lost", "stop_lost"])
        # is_canonical_splice_local was computed earlier; keep conservative check for splice_acceptor/donor keywords
        is_canonical_splice_local = any(x in cons for x in ["splice_acceptor_variant", "splice_donor_variant","splice_acceptor","splice_donor"])

        # CRITICAL: Non-canonical splice region variants (offset >=3) should NOT receive PP3 by default
        is_noncanonical_splice = _is_noncanonical_splice_region_variant(vr)
        
        # =========================================================
        # Determine if PP3 should be blocked
        # =========================================================
        
        # Check if high SpliceAI overrides non-canonical status
        allow_pp3_for_noncanonical = False
        if is_noncanonical_splice:
            spliceai_score = getattr(vr, "spliceai", 0.0)
            if spliceai_score is not None and spliceai_score >= 0.7:
                allow_pp3_for_noncanonical = True
                spliceai_score_str = f"{spliceai_score:.3f}" if spliceai_score is not None else "NA"
                vr.criteria_explanations["PP3_override"] = (
                    f"Non-canonical splice variant but SpliceAI={spliceai_score_str} >= 0.7 - allowing PP3 evaluation"
                )
        
        # Block PP3 if:
        # 1. PVS1 present, OR
        # 2. LOF-like variant, OR  
        # 3. Canonical splice variant, OR
        # 4. Non-canonical splice variant WITHOUT high SpliceAI override (SpliceAI < 0.7)
        block_pp3 = (pvs1_present or cons_lof_like or is_canonical_splice_local or 
                     (is_noncanonical_splice and not allow_pp3_for_noncanonical))
        
        if block_pp3:
            # Do not assign PP3 in these situations
            if is_noncanonical_splice and not allow_pp3_for_noncanonical:
                s_score_str = f"{s_score:.3f}" if s_score is not None else "NA"
                vr.criteria_explanations["PP3"] = (
                    f"PP3 not applied: non-canonical splice region variant (SpliceAI={s_score_str} < 0.7)"
                )
            else:
                vr.criteria_explanations["PP3"] = (
                    "PP3 not applied: variant is LOF/canonical-splice or PVS1 was assigned; "
                    "computational evidence is not stacked with PVS1/LOF for this pipeline."
                )
        else:
            # Only consider PP3 when not blocked
            if is_pp3_moderate:
                assigned.append("PP3_Moderate")
                pts["PP3_Moderate"] = POINTS_MAP.get("PP3_Moderate", 2)
                s_score_str = f"{s_score:.3f}" if s_score is not None else "NA"
                vr.criteria_explanations["PP3"] = (f"PP3_Moderate: computational evidence (SpliceAI={s_score_str}, REVEL={r_score if r_score!=-1 else 'NA'}) met moderate thresholds")
            elif is_pp3_supporting:
                assigned.append("PP3_Supporting")
                pts["PP3_Supporting"] = POINTS_MAP.get("PP3_Supporting", 1)
                
                # Build explanation based on what triggered Supporting
                if r_score >= THRESH.get("REVEL_SUPPORTING", 0.64):
                    r_score_str = f"{r_score:.3f}" if r_score is not None else "NA"
                    vr.criteria_explanations["PP3"] = f"PP3_Supporting: REVEL={r_score_str} >= supporting threshold"
                elif a_score >= THRESH.get("ALPHA_MISSENSE_SUPPORTING", 0.56):
                    a_score_str = f"{a_score:.3f}" if a_score is not None else "NA"
                    vr.criteria_explanations["PP3"] = f"PP3_Supporting: AlphaMissense={a_score_str} >= supporting threshold"
                elif a_category == "P":
                    vr.criteria_explanations["PP3"] = "PP3_Supporting: AlphaMissense category 'P' (pathogenic prediction)"
                elif s_score >= THRESH.get("SPLICEAI_SUPPORTING", 0.20):
                    s_score_str = f"{s_score:.3f}" if s_score is not None else "NA"
                    vr.criteria_explanations["PP3"] = f"PP3_Supporting: SpliceAI={s_score_str} >= supporting threshold"

        # ---------------------------------------------------------------------
        # SUPPRESS PS2: do not allow PS2 (de novo) criteria in automated scoring.
        # ---------------------------------------------------------------------
        try:
            ps2_keys = [k for k in list(pts.keys()) if str(k).upper().startswith("PS2")]
            if ps2_keys:
                for k in ps2_keys:
                    # remove from points map
                    pts.pop(k, None)
                # remove from assigned list too
                assigned = [a for a in assigned if not str(a).upper().startswith("PS2")]
                # record an explanation for PS2 suppression so it appears in audit logs/explanations
                vr.criteria_explanations["PS2_suppressed"] = "PS2 suppressed: de novo evidence (PS2) is not considered in automated scoring."
        except Exception:
            # never fail classification because of suppression bookkeeping
            pass

        # -------------------------
        # Extended DB-weighted scoring (PP5ext/BP5ext) and final classification
        # -------------------------
        try:
            # initialize pts_ext from base pts (includes PVS1, PS1, PMx, PPx, PM3 if present)
            pts_ext = dict(pts)
            # remove legacy PP5/BP5 keys if present (we will use PP5ext/BP5ext only)
            for k in list(pts_ext.keys()):
                if k.startswith("PP5") or k.startswith("BP5"):
                    pts_ext.pop(k, None)
        except Exception:
            pts_ext = {}

        applied_ext_key = None
        applied_ext_points = 0
        clin_match_flag = False
        hgmd_match_flag = False
        ext_expl_str = ""

        # prepare base rule & disease once
        rule = self._get_gene_rule(vr.gene) if hasattr(self, "_get_gene_rule") else {}
        acmg_disease = str(rule.get("disease", "") or "").strip()

        # --- ClinVar block (Consensus-based evaluation) ---
        clin_block_raw = 0.0
        clin_match_flag = False

        CLIN_BASE_MAP = {
            "pathogenic": 3,
            "likely pathogenic": 2,
            "benign": -3,
            "likely benign": -2,
            "uncertain significance": -1,
            "vus": -1,
            "conflicting classifications of pathogenicity": -1,
        }
        CLIN_STAR_COEF = {5: 2.0, 4: 1.5, 3: 1.25, 2: 1.0, 1: 0.75, 0: 0.5}

        try:
            # Get matching submissions (phenotype matches ACMG disease)
            matching_subs_for_scoring = []
            all_clin_subs = getattr(vr, "clinvar_submissions", []) or []
            
            for sub in all_clin_subs:
                phen = str(sub.get("phenotype", "")).strip()
                if acmg_disease and phen:
                    try:
                        matched, _, _ = check_disease_match(acmg_disease, phen)
                        if matched:
                            matching_subs_for_scoring.append(sub)
                            clin_match_flag = True
                    except Exception:
                        pass
            
            if matching_subs_for_scoring:
                # Get consensus classification and star from matching submissions
                consensus_sig = get_consensus_classification_from_submissions(matching_subs_for_scoring)
                
                # For conflicting classifications, always use 1 star
                if consensus_sig and "conflict" in consensus_sig.lower():
                    # For conflicting classifications, always use 1 star (low confidence)
                    consensus_star = 1
                else:
                    consensus_star = calculate_consensus_star(matching_subs_for_scoring)
                
                # Map consensus classification to base score
                sig_lower = consensus_sig.lower()
                
                # CRITICAL: Split by common separators to handle multi-value entries
                # Examples: "Pathogenic/Likely pathogenic", "Pathogenic;Likely pathogenic"
                parts = re.split(r'[;/|]+', sig_lower)
                
                # Track the strongest (most pathogenic) classification
                has_pathogenic = False
                has_likely_pathogenic = False
                has_conflict = False
                has_benign = False
                has_likely_benign = False
                has_vus = False
                
                for part in parts:
                    part = part.strip()
                    if "conflict" in part:
                        has_conflict = True
                    elif "pathogenic" in part and "likely" not in part:
                        has_pathogenic = True
                    elif "likely pathogenic" in part:
                        has_likely_pathogenic = True
                    elif "likely benign" in part:
                        has_likely_benign = True
                    elif "benign" in part and "pathogen" not in part:
                        has_benign = True
                    elif "uncertain" in part or "vus" in part:
                        has_vus = True
                
                # Determine base score based on strongest evidence
                if has_pathogenic:
                    # If ANY part says "Pathogenic" (without likely), use Pathogenic (3 points)
                    sub_base = CLIN_BASE_MAP.get("pathogenic", 3)
                elif has_likely_pathogenic:
                    sub_base = CLIN_BASE_MAP.get("likely pathogenic", 2)
                elif has_conflict:
                    sub_base = CLIN_BASE_MAP.get("conflicting classifications of pathogenicity", -1)
                elif has_benign:
                    sub_base = CLIN_BASE_MAP.get("benign", -3)
                elif has_likely_benign:
                    sub_base = CLIN_BASE_MAP.get("likely benign", -2)
                elif has_vus:
                    sub_base = CLIN_BASE_MAP.get("uncertain significance", -1)
                else:
                    sub_base = 0
                
                # Apply star coefficient
                coef = CLIN_STAR_COEF.get(consensus_star, 0.5)
                clin_block_raw = float(sub_base) * float(coef)
               
                # Apply phenotype match penalty for ClinVar (same as HGMD)
                if not clin_match_flag:
                    if clin_block_raw > 1.0:
                        clin_block_raw = 1.0
                        vr._clinvar_capped = True
                    elif clin_block_raw < -1.0:
                        clin_block_raw = -1.0
                        vr._clinvar_capped = True
                
                # Store consensus info on variant for audit
                vr._clinvar_consensus_sig = consensus_sig
                vr._clinvar_consensus_star = consensus_star
                vr._clinvar_phenotype_match = clin_match_flag
                vr._clin_block_raw_display = f"{clin_block_raw:.3f}"
                
                # Store detailed calculation info for audit (NEW)
                clin_details = []
                clin_details.append(f"sig={consensus_sig}")
                clin_details.append(f"base_points={sub_base}")
                clin_details.append(f"stars={consensus_star}")
                clin_details.append(f"coef={coef:.2f}")
                clin_details.append(f"phenotype_match={clin_match_flag}")
                clin_details.append(f"clin={clin_block_raw:.3f}")
                
                vr._clin_block_details = "; ".join(clin_details)
                vr._clin_block_raw = clin_block_raw
                
            else:
                # No phenotype match - fallback to VCF aggregate
                cls, stars, src = self.dbm.get_clinvar_vcf_aggregate(vr)
                if cls:
                    sig_lower = cls.lower()
                    
                    # Split by common separators
                    parts = re.split(r'[;/|]+', sig_lower)
                    
                    # Track strongest evidence
                    has_pathogenic = False
                    has_likely_pathogenic = False
                    has_conflict = False
                    has_benign = False
                    has_likely_benign = False
                    
                    for part in parts:
                        part = part.strip()
                        if "conflict" in part:
                            has_conflict = True
                        elif "pathogenic" in part and "likely" not in part:
                            has_pathogenic = True
                        elif "likely pathogenic" in part:
                            has_likely_pathogenic = True
                        elif "likely benign" in part:
                            has_likely_benign = True
                        elif "benign" in part and "pathogen" not in part:
                            has_benign = True
                    
                    # Determine base score
                    if has_pathogenic:
                        sub_base = CLIN_BASE_MAP.get("pathogenic", 3)
                    elif has_likely_pathogenic:
                        sub_base = CLIN_BASE_MAP.get("likely pathogenic", 2)
                    elif has_conflict:
                        sub_base = CLIN_BASE_MAP.get("conflicting classifications of pathogenicity", -1)
                    elif has_benign:
                        sub_base = CLIN_BASE_MAP.get("benign", -3)
                    elif has_likely_benign:
                        sub_base = CLIN_BASE_MAP.get("likely benign", -2)
                    else:
                        sub_base = 0
                    
                    coef = CLIN_STAR_COEF.get(stars, 0.5)
                    clin_block_raw = float(sub_base) * float(coef)

                    # Apply phenotype match penalty for ClinVar (same as HGMD)
                    if not clin_match_flag:
                        if clin_block_raw > 1.0:
                            clin_block_raw = 1.0
                            vr._clinvar_capped = True
                        elif clin_block_raw < -1.0:
                            clin_block_raw = -1.0
                            vr._clinvar_capped = True

                    vr._clinvar_consensus_sig = cls
                    vr._clinvar_consensus_star = stars
                    
                    # Store detailed calculation info for audit (NEW for fallback)
                    clin_details = []
                    clin_details.append(f"sig={cls}")
                    clin_details.append(f"base_points={sub_base}")
                    clin_details.append(f"stars={stars}")
                    clin_details.append(f"coef={coef:.2f}")
                    clin_details.append(f"phenotype_match=False")
                    clin_details.append(f"clin={clin_block_raw:.3f}")
                    clin_details.append(f"source=fallback_VCF")
                    
                    vr._clin_block_details = "; ".join(clin_details)
                    vr._clin_block_raw = clin_block_raw
                else:
                    clin_block_raw = 0.0
                    vr._clinvar_consensus_sig = "Not provided"
                    vr._clinvar_consensus_star = 0
                    
                    # Store default details (NEW)
                    vr._clin_block_details = "sig=Not provided; clin=0.000; source=fallback_VCF_no_match"
                    vr._clin_block_raw = 0.0
                
                vr._clinvar_phenotype_match = False
                vr._clin_block_raw_display = f"{clin_block_raw:.3f}"

        except Exception as e:
            logger.debug(f"ClinVar block calculation failed: {e}")
            clin_block_raw = 0.0
            clin_match_flag = False
            vr._clinvar_consensus_sig = "Not provided"
            vr._clinvar_consensus_star = 0
            vr._clin_block_raw_display = "0.000"

        # --- HGMD block (per-variant) ---
        hgmd_class_raw = (vr.hgmd_class or "").strip()
        HGMD_CLASS_MAP = {"DM": 3, "DM?": 1, "DFP": -1, "DP": -2, "FP": -2, "R": 0, "": 0}
        hgmd_key = str(hgmd_class_raw).upper()
        if "DM" in hgmd_key and "DM?" not in hgmd_key:
            hgmd_key = "DM"
        elif "DM?" in hgmd_key:
            hgmd_key = "DM?"
        elif "DFP" in hgmd_key:
            hgmd_key = "DFP"
        elif "DP" in hgmd_key and "DFP" not in hgmd_key:
            hgmd_key = "DP"
        elif "FP" in hgmd_key and "DFP" not in hgmd_key:
            hgmd_key = "FP"

        hgmd_class_points = HGMD_CLASS_MAP.get(hgmd_key, 0)
        try:
            hgmd_pub_count = int(getattr(vr, "hgmd_publication_count", 0) or 0)
        except Exception:
            hgmd_pub_count = 0
        try:
            hgmd_dmsup = int(getattr(vr, "hgmd_dmsupported", 0) or 0)
        except Exception:
            hgmd_dmsup = 0

        hgmd_pub_score = float(hgmd_dmsup) / float(hgmd_pub_count) if hgmd_pub_count else 0.0

        # Map hgmd_pub_score to coefficient
        if hgmd_pub_score >= 0.9:
            hgmd_pub_coef = 2.0
        elif 0.7 <= hgmd_pub_score < 0.9:
            hgmd_pub_coef = 1.5
        elif 0.4 <= hgmd_pub_score < 0.7:
            hgmd_pub_coef = 1.0
        elif 0.2 <= hgmd_pub_score < 0.4:
            hgmd_pub_coef = 0.75
        else:
            hgmd_pub_coef = 0.5
        
        # Persist HGMD pub score + coefficient for downstream reporting/audit
        try:
            vr._hgmd_pub_score = float(hgmd_pub_score)
        except Exception:
            vr._hgmd_pub_score = None
        try:
            vr._hgmd_pub_coef = float(hgmd_pub_coef)
        except Exception:
            vr._hgmd_pub_coef = None

        hgmd_block_raw = float(hgmd_class_points) * float(hgmd_pub_coef)
        
        # Store detailed HGMD calculation info for audit (NEW)
        hgmd_details = []
        hgmd_details.append(f"class={hgmd_key}")
        hgmd_details.append(f"base_points={hgmd_class_points}")
        hgmd_details.append(f"publications={hgmd_pub_count}")
        hgmd_details.append(f"dmsupported={hgmd_dmsup}")
        hgmd_details.append(f"pub_score={hgmd_pub_score:.3f}")
        hgmd_details.append(f"coef={hgmd_pub_coef:.2f}")
        hgmd_details.append(f"phenotype_match={hgmd_match_flag}")
        hgmd_details.append(f"hgmd={hgmd_block_raw:.3f}")
        
        vr._hgmd_block_details = "; ".join(hgmd_details)
        vr._hgmd_block_raw = hgmd_block_raw

        # phenotype match for HGMD
        hgmd_match_flag = False
        hgmd_match_level = ""
        hgmd_match_reason = ""
        hgmd_trait_text = str(vr.hgmd_phen or "")
        if acmg_disease and hgmd_trait_text:
            try:
                matched_h, match_level_h, match_reason_h = check_disease_match(acmg_disease, hgmd_trait_text)
                hgmd_match_flag = bool(matched_h)
                hgmd_match_level = match_level_h or ""
                hgmd_match_reason = match_reason_h or ""
            except Exception:
                hgmd_match_flag = False
                hgmd_match_level = ""
                hgmd_match_reason = ""
        # persist HGMD phenotype-match meta on the variant record for downstream audit
        try:
            vr._hgmd_phenotype_match = bool(hgmd_match_flag)
            vr._hgmd_phenotype_match_level = str(hgmd_match_level)
            vr._hgmd_phenotype_match_reason = str(hgmd_match_reason)
        except Exception:
            pass

        # If no phenotype match in HGMD, cap HGMD block to [-1.0, +1.0]
        if not hgmd_match_flag:
            if hgmd_block_raw > 1.0:
                hgmd_block_raw = 1.0
            elif hgmd_block_raw < -1.0:
                hgmd_block_raw = -1.0

        # Sum blocks and map to PP5ext/BP5ext strength
        raw_total = clin_block_raw + hgmd_block_raw

        applied_ext_key = None
        applied_ext_points = 0
        # mapping ranges:
        if raw_total >= 8.0:
            applied_ext_key = "PP5ext_VeryStrong"
            applied_ext_points = POINTS_MAP.get("PP5ext_VeryStrong", 8)
        elif raw_total >= 4.0:
            applied_ext_key = "PP5ext_Strong"
            applied_ext_points = POINTS_MAP.get("PP5ext_Strong", 4)
        elif raw_total >= 2.0:
            applied_ext_key = "PP5ext_Moderate"
            applied_ext_points = POINTS_MAP.get("PP5ext_Moderate", 2)
        elif raw_total >= 0.5:
            applied_ext_key = "PP5ext_Supporting"
            applied_ext_points = POINTS_MAP.get("PP5ext_Supporting", 1)
        elif -1.9 <= raw_total <= -1:
            # Skip BP5ext for HIGH impact variants (frameshift, stop_gained, splice-site)
            impact = getattr(vr, "impact", "").upper()
            if impact == "HIGH":
                applied_ext_key = None
                applied_ext_points = 0
            else:
                applied_ext_key = "BP5ext_Supporting"
                applied_ext_points = POINTS_MAP.get("BP5ext_Supporting", -1)
        elif -3.9 <= raw_total <= -2.0:
            applied_ext_key = "BP5ext_Moderate" 
            applied_ext_points = POINTS_MAP.get("BP5ext_Moderate", -2)
        elif raw_total <= -4.0:
            applied_ext_key = "BP5ext_Strong"
            applied_ext_points = POINTS_MAP.get("BP5ext_Strong", -4)

        # Apply the ext criterion
        if applied_ext_key:
            pts_ext[applied_ext_key] = applied_ext_points
            if applied_ext_key not in assigned:
                assigned.append(applied_ext_key)

        vr.criteria_points_ext = dict(pts_ext)

        # compute totals (use integer points from POINTS_MAP)
        try:
            total_ext = int(sum(int(vp) for vp in pts_ext.values()))
        except Exception:
            # defensive fallback
            try:
                total_ext = int(round(sum(float(vp) for vp in pts_ext.values())))
            except Exception:
                total_ext = 0

        vr.total_ext_points = total_ext

        # map total_ext to class using CLASSIFICATION_THRESHOLDS exactly
        def _map_points_to_class(total_points: int) -> str:
            for label, min_thr, max_thr in CLASSIFICATION_THRESHOLDS:
                if min_thr is not None and max_thr is not None:
                    if min_thr <= total_points <= max_thr:
                        return label
                elif min_thr is not None and max_thr is None:
                    if total_points >= min_thr:
                        return label
                elif min_thr is None and max_thr is not None:
                    if total_points <= max_thr:
                        return label
            return "VUS"

        extended_class = _map_points_to_class(total_ext)
        vr.extended_class = extended_class

        # persist final scoring into variant record (prefer extended if non-VUS)
        if extended_class and extended_class != "VUS":
            vr.criteria_points_final = dict(pts_ext)
            vr.total_points_final = int(total_ext)
            vr.criteria_points = dict(pts_ext)
            vr.total_points = int(total_ext)
        else:
            try:
                base_total = int(sum(int(x) for x in (pts or {}).values())) if pts else 0
            except Exception:
                base_total = 0
            vr.criteria_points_final = dict(pts)
            vr.total_points_final = base_total
            vr.criteria_points = dict(pts)
            vr.total_points = base_total

        # audit info
        vr._clin_block_raw_display = f"{clin_block_raw:.3f}"
        vr._hgmd_block_raw_display = f"{hgmd_block_raw:.3f}"
        vr._pp5_raw_total = float(raw_total)
        vr._applied_ext_key = applied_ext_key or ""
        vr._clinvar_phenotype_match = bool(clin_match_flag)
        vr._hgmd_phenotype_match = bool(hgmd_match_flag)
        vr._pp5ext_explanation = f"clin={clin_block_raw:.3f}, hgmd={hgmd_block_raw:.3f}, raw_total={raw_total:.3f}, applied={applied_ext_key}"

        try:
            if hasattr(vr, "extended_class") and vr.extended_class:
                vr.automated_class = vr.extended_class
            else:
                vr.automated_class = auto_class
        except Exception:
            vr.automated_class = auto_class

        logger.debug("Variant %s:%s assigned pts=%s total_ext=%s extended_class=%s automated_class=%s",
                    vr.chrom, vr.pos, vr.criteria_points, getattr(vr, "total_ext_points", None), getattr(vr, "extended_class", None), getattr(vr, "automated_class", None))
        # make automated_class reflect extended_class (if non-VUS) so downstream triage sees ext scoring 
        try:
            # If the extended model produced a non-VUS class, prefer it as the automated class
            if extended_class and extended_class != "VUS":
                auto_class = extended_class
            # persist for other code that expects this attribute
            vr.automated_class = auto_class
            # also persist "final" points used for automated class (prefer pts_ext if it influenced class)
            # store both for debugging/audit
            vr.criteria_points_final = pts_ext if extended_class and extended_class != "VUS" else pts
            vr.total_points_final = sum(int(x) for x in (vr.criteria_points_final or {}).values())
            logger.debug("Variant scoring summary: pts=%s total=%s pts_ext=%s total_ext=%s automated_class=%s extended_class=%s",
                         pts, vr.total_points, pts_ext, vr.total_ext_points, auto_class, extended_class)
        except Exception as _e:
            logger.debug("Failed to persist automated/extended class info: %s", _e)

        # ---------------------------------------------------------------------
        # Additional rule: If exactly TWO pathogenic-type criteria are present,
        # the variant cannot be called "Pathogenic" (but may still be "Likely pathogenic").
        # This enforces an extra conservative guard: strict Pathogenic requires >2 such criteria.
        # ---------------------------------------------------------------------
        try:
            path_cnt = sum(1 for x in assigned if any(p in x for p in ["PVS", "PS", "PM", "PP"]))
            # If algorithmic/extended model produced "Pathogenic" but only 2 pathogenic criteria -> downgrade to Likely pathogenic
            if ((extended_class or "").startswith("Pathogenic") ) and path_cnt == 2:
                # downgrade Pathogenic -> Likely pathogenic
                extended_class = "Likely pathogenic"
                # ensure variant record reflects the downgrade for downstream code
                try:
                    vr.extended_class = "Likely pathogenic"
                    # CRITICAL: Mark that this is still reportable despite downgrade
                    vr._downgraded_from_pathogenic = True
                    vr._reportable_after_downgrade = True
                except Exception:
                    pass
                # add an audit/manual note explaining downgrade
                note = "Downgraded Pathogenic -> Likely pathogenic: only 2 pathogenic-type criteria present"
                if not hasattr(vr, "manual_reasons") or vr.manual_reasons is None:
                    vr.manual_reasons = []
                vr.manual_reasons.append(note)
                # keep criteria_points etc intact; only class is adjusted
        except Exception:
            # do not break classification pipeline if counting fails
            path_cnt = sum(1 for x in assigned if any(p in x for p in ["PVS", "PS", "PM", "PP"]))

        # Apply minimum-pathogenic-criteria guard once (use pathogenic-criterion count)
        path_cnt = sum(1 for x in assigned if any(p in x for p in ["PVS", "PS", "PM", "PP"]))
        if (extended_class.startswith("Path") or extended_class.startswith("Likely path")) and path_cnt < MIN_PATHOGENIC_CRITERIA:
            # downgrade to VUS and record reason
            vr.automated_class = "VUS"
            vr.extended_class = "VUS"
            vr.manual_reasons = getattr(vr, "manual_reasons", []) + [f"Downgraded: <{MIN_PATHOGENIC_CRITERIA} pathogenic criteria present"]
            extended_class = "VUS"

        # Ensure explanations exist for any assigned criterion not yet explained (single pass)
        for crit in assigned:
            base = crit.split("_")[0]
            if base not in vr.criteria_explanations:
                vr.criteria_explanations[base] = f"{base}: applied as {crit}; detailed reasoning not captured"

        # Final automated class for return: use vr.automated_class (which mirrors extended_class unless downgraded)
        auto_class = getattr(vr, "automated_class", extended_class or "VUS")

        # Helpful debug output for triage
        logger.debug(
            "evaluate_variant: %s:%s %s assigned=%s pts_ext=%s total_ext=%s auto_class=%s manual_reasons=%s",
            vr.chrom, vr.pos, getattr(vr, "gene", "?"),
            assigned, pts_ext, vr.total_points, vr.automated_class, (";".join(getattr(vr, "manual_reasons", [])) if getattr(vr, "manual_reasons", None) else "")
        )

        # =========================================================================
        # CRITICAL GUARD: Enforce minimum criteria count for Pathogenic classification
        # - Pathogenic requires >2 pathogenic criteria (>2, не >=2)
        # - If exactly 2 criteria present -> downgrade to Likely pathogenic
        # =========================================================================
        try:
            path_cnt = sum(1 for x in assigned if any(p in x for p in ["PVS", "PS", "PM", "PP"]))
            
            if vr.automated_class == "Pathogenic":
                if path_cnt <= 2:
                    vr.automated_class = "Likely pathogenic"
                    vr.extended_class = "Likely pathogenic"
                    note = f"Downgraded Pathogenic -> Likely pathogenic: only {path_cnt} pathogenic criteria present (requires >2)"
                    if note not in vr.manual_reasons:
                        vr.manual_reasons.append(note)
            
            elif vr.automated_class == "Likely pathogenic":
                if path_cnt < 1:
                    vr.automated_class = "VUS"
                    vr.extended_class = "VUS"
                    note = "Downgraded Likely pathogenic -> VUS: no strong pathogenic criteria"
                    if note not in vr.manual_reasons:
                        vr.manual_reasons.append(note)
                        
        except Exception as e:
            logger.debug(f"Failed to apply criteria count guard: {e}")

        # Final return now reflects extended points/class as authoritative
        return assigned, pts_ext, vr.automated_class, False, getattr(vr, "manual_reasons", [])

    def _check_canonical_splice_position(self, vr: VariantRecord) -> bool:
        return "splice_acceptor" in (vr.consequence or "") or "splice_donor" in (vr.consequence or "")

    def _check_inframe_small_noncritical(self, vr: VariantRecord) -> bool:
        hgvsp_clean = clean_hgvsp(vr.hgvsp)
        if not hgvsp_clean:
            return False
        size = 0
        del_match = re.search(r'del(\d+)[A-Za-z]*$', hgvsp_clean)
        if del_match:
            try:
                size = int(del_match.group(1))
            except:
                pass
        ins_match = re.search(r'ins(\d+)[A-Za-z]*$', hgvsp_clean)
        if ins_match:
            try:
                size = int(ins_match.group(1))
            except:
                pass
        single_del_match = re.search(r'p\.[A-Za-z]{1,3}\d+del$', hgvsp_clean)
        if single_del_match and size == 0:
            size = 1
        dup_match = re.search(r'dup(\d+)[A-Za-z]*$', hgvsp_clean)
        if dup_match:
            try:
                size = int(dup_match.group(1))
            except:
                pass
        if size == 0:
            return False
        return size < 10

    def _compare_revel_scores(self, variant_revel: Optional[float], reference_revel: Optional[float]) -> str:
        if variant_revel is None or reference_revel is None:
            return "insufficient_data"
        try:
            v_revel = float(variant_revel)
            r_revel = float(reference_revel)
            if v_revel >= r_revel:
                return "greater_or_equal"
            else:
                return "less"
        except (ValueError, TypeError):
            return "insufficient_data"

    def apply_pm3_batch(self, candidates: List[VariantRecord]) -> None:
        """
        PM3 batch scoring (reworked):

        - Treat every sample independently (mono-sample mode). This method still groups
          variants by sample+gene to compute PM3 contributions across variants in the same sample.
        - PM3 maximum strength is SUPPORTING (PM3_Supporting).
        - PM3 is applied ONLY when there is NO evidence that candidate variants are in CIS.
          If any confirmed cis evidence exists for a given variant relative to partner(s),
          PM3 is NOT applied for that variant and an audit note is added.
        - Phasing evidence (parental or phased GT) is respected to detect confirmed trans/cis.
        - Contributions logic (unphased vs confirmed in-trans) is preserved for scoring,
          but the final assigned strength is capped to MP3_SUPPORTING.
        """
        by_sample_gene = defaultdict(lambda: defaultdict(list))
        for v in candidates:
            by_sample_gene[v.sample][v.gene].append(v)

        for sample, genes in by_sample_gene.items():
            for gene, varlist in genes.items():
                rule = self._get_gene_rule(gene)
                moi = str(rule.get("moi", "")).upper()
                # Only consider PM3 for genes with recessive MOI
                if not ("AR" in moi or "RECESSIVE" in moi or "XLR" in moi):
                    continue

                for v in varlist:
                    pm3_score = 0.0
                    unphased_contributions = 0.0
                    criteria_details = []
                    cis_found_for_variant = False

                    # If variant is homozygous -> count as 0.5 (carrier/homozygous support) per original pipeline
                    if _is_hom(v.proband_gt):
                        pm3_score += 0.5
                        unphased_contributions += 0.5
                        criteria_details.append("homozygous: +0.5")

                    # compare v to other variants in same sample/gene
                    for other_v in varlist:
                        if other_v == v:
                            continue

                        # Determine simple phasing relationships using available genotype fields.
                        in_trans_confirmed = False
                        in_cis_confirmed = False

                        # parental-derived in-trans detection (if parents were provided as mono-samples,
                        # that does not serve as per-variant parental phasing here; however if parental
                        # GTs were present on the same VariantRecord they would appear in father_gt/mother_gt)
                        try:
                            if (_is_het(v.father_gt) and _is_wt(v.mother_gt) and
                                _is_het(other_v.mother_gt) and _is_wt(other_v.father_gt)):
                                in_trans_confirmed = True
                        except Exception:
                            pass
                        try:
                            if (_is_het(v.mother_gt) and _is_wt(v.father_gt) and
                                _is_het(other_v.father_gt) and _is_wt(other_v.mother_gt)):
                                in_trans_confirmed = True
                        except Exception:
                            pass

                        # simple same-parent het => cis evidence
                        try:
                            if (_is_het(v.father_gt) and _is_het(other_v.father_gt)) or (_is_het(v.mother_gt) and _is_het(other_v.mother_gt)):
                                in_cis_confirmed = True
                        except Exception:
                            pass

                        # Additionally, use fully-phased proband GT (|) if present to detect cis/trans
                        try:
                            def _phased_index(gt: str) -> Optional[int]:
                                if not gt or '|' not in gt:
                                    return None
                                parts = gt.split('|')
                                try:
                                    if parts.count('1') == 1:
                                        return parts.index('1')
                                except Exception:
                                    pass
                                for idx, p in enumerate(parts):
                                    if str(p) == '1':
                                        return idx
                                return None

                            idx_v = _phased_index(v.proband_gt)
                            idx_o = _phased_index(other_v.proband_gt)
                            if idx_v is not None and idx_o is not None:
                                # same index -> cis; different index -> trans
                                if idx_v == idx_o:
                                    in_cis_confirmed = True
                                else:
                                    in_trans_confirmed = True
                        except Exception:
                            pass

                        # External cooccurrence evidence (recs_map/cooccurrence) may be present in the broader pipeline,
                        # but apply_pm3_batch is conservative here and only uses explicit cis confirmation to disallow PM3.
                        # Record cis evidence at per-variant level:
                        if in_cis_confirmed:
                            cis_found_for_variant = True
                            criteria_details.append(f"confirmed in cis with {other_v.hgvsp} (no PM3)")
                            # No PM3 points contributed for cis pairs; continue to next partner
                            continue

                        # if confirmed in trans, award stronger contribution:
                        if in_trans_confirmed:
                            if other_v.automated_class in ("Pathogenic", "Likely pathogenic"):
                                pm3_score += 1.0
                                criteria_details.append(f"confirmed in trans with P/LP ({other_v.hgvsp}): +1.0")
                            elif other_v.automated_class == "VUS":
                                pm3_score += 0.25
                                criteria_details.append(f"confirmed in trans with VUS ({other_v.hgvsp}): +0.25")
                            continue

                        # unphased contribution logic (capped aggregate below)
                        contribution = 0.0
                        if other_v.automated_class == "Pathogenic":
                            contribution = 0.5
                            criteria_details.append(f"unphased with P ({other_v.hgvsp}): +0.5")
                        elif other_v.automated_class == "Likely pathogenic":
                            contribution = 0.5
                            criteria_details.append(f"unphased with LP ({other_v.hgvsp}): +0.5")

                        # cap unphased contributions at 1.0 per original design
                        if unphased_contributions + contribution <= 1.0:
                            unphased_contributions += contribution
                            pm3_score += contribution
                        elif unphased_contributions < 1.0:
                            remaining = 1.0 - unphased_contributions
                            unphased_contributions += remaining
                            pm3_score += remaining
                            criteria_details.append(f"unphased (capped): +{remaining:.2f}")

                    # ensure unphased_contributions does not exceed 1.0
                    if unphased_contributions > 1.0:
                        excess = unphased_contributions - 1.0
                        pm3_score -= excess
                        unphased_contributions = 1.0

                    # Determine PM3 final assignment with policy:
                    # - apply PM3 only if there is NO confirmed cis evidence for this variant
                    # - maximum strength is Supporting when pm3_score >= 0.5
                    pm3_level = None
                    if cis_found_for_variant:
                        # Do not apply PM3 if variant has any confirmed cis evidence
                        v.manual_reasons.append("PM3 not applied: confirmed cis evidence detected with partner variant(s)")
                    else:
                        # apply PM3 Supporting if pm3_score >= 0.5 (this is the pipeline maximum)
                        if pm3_score >= 0.5:
                            pm3_level = "PM3_Supporting"
                        else:
                            pm3_level = None

                    # Apply PM3_Supporting to variant record if appropriate and not already present
                    if pm3_level and pm3_level not in v.criteria_assigned:
                        v.criteria_assigned.append(pm3_level)
                        v.criteria_points[pm3_level] = POINTS_MAP.get(pm3_level, 1)
                        if criteria_details:
                            details_str = ", ".join(criteria_details)
                            v.manual_reasons.append(f"PM3 ({pm3_level}): total={pm3_score:.2f} [{details_str}]")
                        else:
                            v.manual_reasons.append(f"PM3 ({pm3_level}): total={pm3_score:.2f}")

                    # annotate phase field for downstream triage (trans/cis/mixed/unknown)
                    try:
                        trans_found = False; cis_found = False
                        for other_v in varlist:
                            if other_v is v: continue
                            if (_is_het(v.father_gt) and _is_het(other_v.mother_gt)) or (_is_het(v.mother_gt) and _is_het(other_v.father_gt)):
                                trans_found = True
                            if (_is_het(v.father_gt) and _is_het(other_v.father_gt)) or (_is_het(v.mother_gt) and _is_het(other_v.mother_gt)):
                                cis_found = True
                        if trans_found and not cis_found:
                            v.phase = "trans"
                        elif cis_found and not trans_found:
                            v.phase = "cis"
                        elif trans_found and cis_found:
                            v.phase = "mixed"
                        else:
                            v.phase = "unknown"
                    except Exception:
                        v.phase = getattr(v, "phase", "unknown")

# ---------------------------
# Deduplication helper
# ---------------------------
def deduplicate_candidates(candidates: List[VariantRecord]) -> List[VariantRecord]:
    keymap: Dict[Tuple[str,str,int,str,str,str], VariantRecord] = {}
    for v in candidates:
        # CRITICAL: Include sample in the key to keep variants from different individuals separate
        key = (
            v.sample or "unknown",           # Sample name is crucial
            v.chrom or "",                    # Chromosome
            int(v.pos),                        # Position
            v.ref or "",                       # Reference allele
            v.alt or "",                        # Alternative allele
            (v.gene or "").upper()              # Gene (normalized)
        )
        
        if key not in keymap:
            keymap[key] = v
            continue
            
        # If duplicate exists, keep the one with better annotations
        existing = keymap[key]
        
        # Prefer MANE transcript
        e_ann = getattr(existing, "_raw_ann", {}) or {}
        v_ann = getattr(v, "_raw_ann", {}) or {}
        e_can = str(e_ann.get("MANE_PLUS_CLINICAL") or e_ann.get("MANE_SELECT") or e_ann.get("CANONICAL") or "").upper()
        v_can = str(v_ann.get("MANE_PLUS_CLINICAL") or v_ann.get("MANE_SELECT") or v_ann.get("CANONICAL") or "").upper()
        
        if v_can in ("YES","Y","TRUE","1") and not (e_can in ("YES","Y","TRUE","1")):
            keymap[key] = v
            continue
            
        # Prefer higher total points
        try:
            ev_pts = int(getattr(existing, "total_points", 0) or 0)
        except Exception:
            ev_pts = 0
        try:
            v_pts = int(getattr(v, "total_points", 0) or 0)
        except Exception:
            v_pts = 0
            
        if v_pts > ev_pts:
            keymap[key] = v
            continue
            
    return list(keymap.values())

# ---------------------------
# Delins detection helper
# ---------------------------
def detect_delins(candidates: List[VariantRecord]) -> None:
    """
    Detect potential delins (complex adjacent indels) by clustering nearby variants.

    Revised behavior:
      - Build clusters as connected components of nearby variants (within ADJ_THRESHOLD bp).
      - Only mark is_delins_part = True for variants in a cluster if the cluster contains at least
        one Pathogenic/Likely-pathogenic variant that is in the SAME GENE as the candidate variant.
        (We require gene concordance to avoid flagging variants in a gene because of a P/LP in a different gene.)
      - For marking an individual variant as is_delins_part, compute the minimal genomic distance
        from that variant to any Pathogenic member of the same gene in the same cluster, and require
        that distance to be <= ADJ_THRESHOLD.
      - If a cluster contains P/LP variants only in other genes, do NOT flag any member as delins_part;
        instead add an explanatory audit note.
      - Benign variants are filtered out and not flagged as delins parts.
      - Phasing evidence (cis/trans) is recorded in manual_reasons as before.
    """
    def key_of(v: VariantRecord) -> str:
        return f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.sample}"

    def _is_benign_variant(v: VariantRecord) -> bool:
        """
        Conservative check whether a variant is benign enough to be ignored in delins detection.
        """
        # 1. Automated class
        try:
            if getattr(v, "automated_class", "") in ("Benign", "Likely benign"):
                return True
        except Exception:
            pass
        try:
            total_pts = int(getattr(v, "total_ext_points", getattr(v, "total_points", 0)) or 0)
            ext_cls = getattr(v, "extended_class", getattr(v, "automated_class", "VUS"))
            if ext_cls == "VUS" and total_pts < 3:
                return True  # Treat as benign for delins detection
        except Exception:
            pass

        # 2. ClinVar high-quality benign (>=2 stars)
        try:
            if getattr(v, "clinvar_stars", "") and str(v.clinvar_stars).isdigit() and int(v.clinvar_stars) >= 2:
                sig = str(getattr(v, "clinvar_sig", "") or "").lower()
                if ("benign" in sig and "pathogen" not in sig) or ("likely benign" in sig):
                    return True
        except Exception:
            pass

        # 3. HGMD benign classes with support
        try:
            hgmd_class = str(getattr(v, "hgmd_class", "") or "").upper()
            if hgmd_class in ("DP", "FP", "DFP"):
                pubs = int(getattr(v, "hgmd_publication_count", 0) or 0)
                if pubs >= 2:
                    return True
        except Exception:
            pass

        # 4. Strong negative points
        try:
            total_points = int(getattr(v, "total_points", 0) or 0)
            if total_points <= -4:
                return True
        except Exception:
            pass

        return False

    # Group by sample and chromosome
    by_sample_chr = defaultdict(lambda: defaultdict(list))
    for v in candidates:
        try:
            sam = v.sample or ""
            chrom = v.chrom or ""
            by_sample_chr[sam][chrom].append(v)
        except Exception:
            continue

    # Initialize delins metadata on all variants
    for v in candidates:
        if not hasattr(v, "manual_reasons") or v.manual_reasons is None:
            v.manual_reasons = []
        v._delins_flag = False
        v._delins_cluster_id = None
        v._delins_cluster_members = []
        v.is_delins_part = False

    cluster_counter = 0

    def _phased_hap_index(gt: str) -> Optional[int]:
        """
        Return haplotype index (0 or 1) where ALT allele '1' is located if GT is phased and exactly one '1' present.
        Otherwise return None.
        """
        if not gt or '|' not in gt:
            return None
        parts = gt.split('|')
        try:
            if parts.count('1') == 1:
                return parts.index('1')
            for idx, p in enumerate(parts):
                if str(p) == '1':
                    return idx
        except Exception:
            return None
        return None

    ADJ_THRESHOLD = 10  # adjacency threshold in base pairs

    for sample, chrom_map in by_sample_chr.items():
        for chrom, varlist in chrom_map.items():
            if len(varlist) < 2:
                continue

            # Compute approximate genomic span for a variant using ref length
            def _span(v: VariantRecord):
                try:
                    s = int(v.pos)
                    rlen = len(v.ref) if v.ref is not None else 1
                    e = s + max(1, int(rlen)) - 1
                    return s, e
                except Exception:
                    return int(v.pos), int(v.pos)

            # Sort and build adjacency graph: edges between variants within ADJ_THRESHOLD bp
            varlist_sorted = sorted(varlist, key=lambda x: int(x.pos))
            spans = { key_of(vv): _span(vv) for vv in varlist_sorted }
            neighbors = { key_of(vv): set() for vv in varlist_sorted }
            for i in range(len(varlist_sorted)):
                vi = varlist_sorted[i]; ki = key_of(vi)
                s_i, e_i = spans[ki]
                j = i + 1
                while j < len(varlist_sorted):
                    vj = varlist_sorted[j]; kj = key_of(vj)
                    s_j, e_j = spans[kj]
                    gap = 0 if s_j <= e_i else (s_j - e_i - 1)
                    if gap > ADJ_THRESHOLD:
                        break
                    neighbors[ki].add(kj)
                    neighbors[kj].add(ki)
                    j += 1

            # Find connected components (clusters)
            visited = set()
            key_to_var = { key_of(vv): vv for vv in varlist_sorted }
            for v_obj in varlist_sorted:
                start_k = key_of(v_obj)
                if start_k in visited:
                    continue
                if not neighbors.get(start_k):
                    visited.add(start_k)
                    continue

                comp_keys = []
                queue = [start_k]
                while queue:
                    cur = queue.pop()
                    if cur in visited:
                        continue
                    visited.add(cur)
                    comp_keys.append(cur)
                    for nb in neighbors.get(cur, set()):
                        if nb not in visited:
                            queue.append(nb)

                if len(comp_keys) < 2:
                    continue

                comp_vars = [key_to_var[kk] for kk in comp_keys if kk in key_to_var]
                if not comp_vars:
                    continue

                cluster_counter += 1
                member_str = ";".join(comp_keys)

                # Evaluate phasing inside component
                any_phased = False
                any_cis = False
                any_trans = False
                phasing_pairs = []
                for a_idx in range(len(comp_vars)):
                    for b_idx in range(a_idx + 1, len(comp_vars)):
                        a = comp_vars[a_idx]; b = comp_vars[b_idx]
                        ai = _phased_hap_index(getattr(a, "proband_gt", ""))
                        bi = _phased_hap_index(getattr(b, "proband_gt", ""))
                        if ai is not None and bi is not None:
                            any_phased = True
                            if ai == bi:
                                any_cis = True
                                phasing_pairs.append((a, b, "cis"))
                            else:
                                any_trans = True
                                phasing_pairs.append((a, b, "trans"))

                # Identify pathogenic members in this cluster (Pathogenic/LP or high-quality DB pathogenic)
                pathogenic_members = []
                for vv in comp_vars:
                    try:
                        if getattr(vv, "automated_class", "") in ("Pathogenic", "Likely pathogenic") or getattr(vv, "is_hq_pathogenic_db", False):
                            pathogenic_members.append(vv)
                    except Exception:
                        continue

                # If component has no pathogenics -> do not mark delins parts. Add audit notes.
                if not pathogenic_members:
                    for vv in comp_vars:
                        note = (f"Cluster id={cluster_counter} detected (members={member_str}) - no Pathogenic/LP present; "
                                "delins NOT flagged; requires P/LP to consider delins.")
                        if note not in vv.manual_reasons:
                            vv.manual_reasons.append(note)
                        vv._delins_cluster_id = cluster_counter
                        vv._delins_cluster_members = comp_keys.copy()
                        vv._delins_flag = False
                        vv.is_delins_part = False
                    continue

                # Build helper: for each candidate variant, select pathogenic members that belong to SAME GENE (normalized)
                path_genes = set()
                for p in pathogenic_members:
                    try:
                        path_genes.add(normalize_gene_name(getattr(p, "gene", "") or getattr(p, "gene_raw", "")))
                    except Exception:
                        pass

                # If there are pathogenic members but none share any gene with cluster members,
                # treat the cluster as informational only (do not flag). Record explanatory notes.
                # We'll decide per-variant using same-gene pathogenic subset.
                any_same_gene_pathogenic = False
                for vv in comp_vars:
                    try:
                        vv_gene_norm = normalize_gene_name(getattr(vv, "gene", "") or getattr(vv, "gene_raw", ""))
                        if vv_gene_norm and vv_gene_norm in path_genes:
                            any_same_gene_pathogenic = True
                            break
                    except Exception:
                        continue

                if not any_same_gene_pathogenic:
                    # There are pathogenic members but in different gene(s); do not flag as delins_part.
                    for vv in comp_vars:
                        note = (f"Cluster id={cluster_counter} contains Pathogenic variant(s) in other gene(s) "
                                f"{sorted(list(path_genes))}; no same-gene P/LP -> delins NOT flagged for members.")
                        if note not in vv.manual_reasons:
                            vv.manual_reasons.append(note)
                        vv._delins_cluster_id = cluster_counter
                        vv._delins_cluster_members = comp_keys.copy()
                        vv._delins_flag = False
                        vv.is_delins_part = False
                    continue

                # Compute minimal distance from each variant to pathogenic members in the SAME gene,
                # and mark as delins_part only when <= ADJ_THRESHOLD and variant is not benign.
                def _min_dist_to_path_same_gene(vv: VariantRecord, path_list: List[VariantRecord]) -> int:
                    try:
                        p = int(vv.pos)
                        dmin = None
                        for pm in path_list:
                            try:
                                if normalize_gene_name(pm.gene or pm.gene_raw) != normalize_gene_name(vv.gene or vv.gene_raw):
                                    continue
                                pp = int(pm.pos)
                                d = abs(pp - p)
                                if dmin is None or d < dmin:
                                    dmin = d
                            except Exception:
                                continue
                        return dmin if dmin is not None else 9999999
                    except Exception:
                        return 9999999

                # Now process each member: filter benign, compute min distance to same-gene pathogenics and flag accordingly
                for vv in comp_vars:
                    # Filter benign variants completely (not flagged)
                    if _is_benign_variant(vv):
                        vv._delins_flag = False
                        vv._delins_cluster_id = cluster_counter
                        vv._delins_cluster_members = comp_keys.copy()
                        vv.is_delins_part = False
                        vv._filtered = True
                        vv._exclude_from_reports = True
                        # ensure not included in auto/manual lists downstream
                        vv._in_auto = False
                        vv._in_manual = False
                        continue

                    # find same-gene pathogenic members for this variant
                    same_gene_paths = [p for p in pathogenic_members if normalize_gene_name(p.gene or p.gene_raw) == normalize_gene_name(vv.gene or vv.gene_raw)]
                    if not same_gene_paths:
                        # no same-gene pathogenic partner for this vv -> do not flag
                        vv._delins_flag = False
                        vv._delins_cluster_id = cluster_counter
                        vv._delins_cluster_members = comp_keys.copy()
                        vv.is_delins_part = False
                        note = (f"Cluster id={cluster_counter} contains P/LP but none in the same gene as this variant ({vv.gene or vv.gene_raw}) -> not flagged")
                        if note not in vv.manual_reasons:
                            vv.manual_reasons.append(note)
                        continue

                    dist = _min_dist_to_path_same_gene(vv, same_gene_paths)
                    if dist <= ADJ_THRESHOLD:
                        vv._delins_flag = True
                        vv._delins_cluster_id = cluster_counter
                        vv._delins_cluster_members = comp_keys.copy()
                        vv.is_delins_part = True
                        pair_notes = []
                        for (p1, p2, rel) in phasing_pairs:
                            pair_notes.append(f"{p1.chrom}:{p1.pos}-{p2.chrom}:{p2.pos}:{rel}")
                        note = (f"Delins cluster (cluster_id={cluster_counter}; members={member_str}); "
                                f"nearest same-gene Pathogenic distance={dist}bp; phasing evidence: {', '.join(pair_notes) if pair_notes else 'none'}")
                        if note not in vv.manual_reasons:
                            vv.manual_reasons.append(note)
                    else:
                        vv._delins_flag = False
                        vv._delins_cluster_id = cluster_counter
                        vv._delins_cluster_members = comp_keys.copy()
                        vv.is_delins_part = False
                        note = (f"Cluster id={cluster_counter} contains Pathogenic member(s) in same gene but this variant is {dist}bp away "
                                f"from nearest same-gene Pathogenic (> {ADJ_THRESHOLD}bp) -> not flagged as delins_part")
                        if note not in vv.manual_reasons:
                            vv.manual_reasons.append(note)

                # If confirmed CIS evidence exists, annotate non-benign delins parts for curator prioritization
                if any_cis:
                    for (p1, p2, rel) in phasing_pairs:
                        if rel == "cis":
                            for vv in (p1, p2):
                                if not _is_benign_variant(vv) and getattr(vv, "is_delins_part", False):
                                    cis_note = f"Phasing: confirmed CIS pair in cluster {cluster_counter} -> delins likely (cis evidence)"
                                    if cis_note not in vv.manual_reasons:
                                        vv.manual_reasons.append(cis_note)

# ---------------------------
# Position rules loader & checker (from recommendations table)
# ---------------------------
def load_recommendations_table(path: str) -> Dict[str, Dict[str,Any]]:
    """
    Parse ACMG SF recommendations table for gene-specific reporting guidance.
    
    This table contains additional constraints beyond the main ACMG SF list:
      - Requires_biallelic: Only report if biallelic variants are present.
      - Missense_only: Limit reporting to missense variants.
      - Hemizygous_only_male: Report only in males with hemizygous variants.
      - Positional_rules: Transcript-specific exon ranges or protein position limits.
      - Exclusions: Variants not recommended for return.
    
    Args:
        path (str): Path to CSV/TSV recommendations file.
    
    Returns:
        Dict[str, Dict]: Keyed by normalized gene name, containing rule dictionaries.
    """
    recs = defaultdict(dict)
    if not path or not os.path.exists(path):
        return recs
    if not HAVE_PANDAS:
        logger.warning("Pandas not installed; cannot parse recommendations table.")
        return recs
    try:
        df = pd.read_csv(path, sep=None, engine='python', dtype=str).fillna("")
    except Exception as e:
        logger.warning(f"Failed to read recommendations table {path}: {e}")
        return recs
    cols = {c.lower().strip(): c for c in df.columns}
    gene_col = cols.get("gene symbol") or cols.get("gene") or next(iter(df.columns), None)
    comment_col = cols.get("reporting guidance comment") or cols.get("reporting guidance") or cols.get("comment") or cols.get("reporting guidance comment")
    for _, r in df.iterrows():
        gene_raw = str(r.get(gene_col, "")).strip() if gene_col else ""
        if not gene_raw:
            continue
        g = normalize_gene_name(gene_raw)
        comment = str(r.get(comment_col, "")).strip() if comment_col else ""
        comment_l = comment.lower()
        rules = {
            "requires_biallelic": False,
            "missense_only": False,
            "hemizygous_only_male": False,
            "exclusions": [],
            "position_rules": {}
        }
        if "biallelic" in comment_l or "only reportable if there are two" in comment_l or "only reportable if there are biallelic" in comment_l:
            rules["requires_biallelic"] = True
        if "missense only" in comment_l or "missense variants only" in comment_l:
            rules["missense_only"] = True
        if "hemizygous p/lp variants in males are reportable" in comment_l and "females are not reportable" in comment_l:
            rules["hemizygous_only_male"] = True
        if "not recommended for return" in comment_l:
            rules["exclusions"].append(comment)
        pos_rules = {}
        # transcript-specific exons: look for "exons X-Y of NM_xxx" or "exons X-Y ... NM_..."
        m = re.search(r'(nm_[0-9]+(?:\.\d+)?)[^\n\r]{0,80}?exons?\s*(\d+)[\-\–—](\d+)', comment, flags=re.I)
        if not m:
            m = re.search(r'exons?\s*(\d+)[\-\–—](\d+)[^\n\r]{0,80}?of\s*(nm_[0-9]+(?:\.\d+)?)', comment, flags=re.I)
            if m:
                if m.group(3):
                    tx = m.group(3)
                    start = int(m.group(1))
                    end = int(m.group(2))
                    pos_rules.setdefault("transcript_exons", []).append({"transcript": tx, "start_exon": start, "end_exon": end})
        else:
            tx = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            pos_rules.setdefault("transcript_exons", []).append({"transcript": tx, "start_exon": start, "end_exon": end})
        # generic exon range
        if "transcript_exons" not in pos_rules:
            m2 = re.search(r'exons?\s*(\d+)[\-\–—](\d+)', comment, flags=re.I)
            if m2:
                pos_rules["exon_range"] = {"start_exon": int(m2.group(1)), "end_exon": int(m2.group(2))}
        # protein positions
        m3 = re.search(r'(?:positions?|aa positions?|amino acid positions?)\s*(\d+)[\-\–—](\d+)', comment, flags=re.I)
        if m3:
            pos_rules["protein_range"] = {"start_aa": int(m3.group(1)), "end_aa": int(m3.group(2))}
        if pos_rules:
            rules["position_rules"] = pos_rules
        rules["raw_comment"] = comment
        recs[g] = rules
    logger.info(f"Loaded reporting guidance (with positional rules) for {len(recs)} genes")
    return recs

def variant_matches_position_rules(v: VariantRecord, pos_rules: Dict[str,Any]) -> bool:
    """
    Returns True if variant v is allowed by pos_rules.
    pos_rules may contain:
      - transcript_exons: list of {transcript, start_exon, end_exon}
      - exon_range: {start_exon, end_exon}
      - protein_range: {start_aa, end_aa}
    """
    if not pos_rules:
        return True
    # 1) Transcript-specific exon ranges
    tx_rules = pos_rules.get("transcript_exons", [])
    if tx_rules:
        for tr in tx_rules:
            tx_req = tr.get("transcript")
            if not tx_req:
                continue
            tx_req_base = tx_req.split('.')[0].upper()
            v_tx = (v.transcript or "").split('.')[0].upper()
            if v_tx and v_tx == tx_req_base:
                exon_str = v.exon or ""
                try:
                    if "/" in exon_str:
                        curr = int(exon_str.split("/")[0])
                        if tr["start_exon"] <= curr <= tr["end_exon"]:
                            return True
                    else:
                        m = re.search(r'(\d+)', exon_str)
                        if m:
                            curr = int(m.group(1))
                            if tr["start_exon"] <= curr <= tr["end_exon"]:
                                return True
                except Exception:
                    pass
        return False
    # 2) Generic exon range
    exon_range = pos_rules.get("exon_range")
    if exon_range:
        exon_str = v.exon or ""
        try:
            if "/" in exon_str:
                curr = int(exon_str.split("/")[0])
                if exon_range["start_exon"] <= curr <= exon_range["end_exon"]:
                    return True
            else:
                m = re.search(r'(\d+)', exon_str)
                if m:
                    curr = int(m.group(1))
                    if exon_range["start_exon"] <= curr <= exon_range["end_exon"]:
                        return True
        except Exception:
            pass
        return False
    # 3) Protein position range
    pr = pos_rules.get("protein_range")
    if pr:
        aa = v.aa_pos or parse_protein_pos(v.hgvsp or "")
        if aa and pr["start_aa"] <= aa <= pr["end_aa"]:
            return True
        return False
    # Default allow
    return True

# ---------------------------
# Output helpers & triage
# ---------------------------

def select_moi_for_variant(rule: Dict[str, Any], variant: VariantRecord) -> str:
    """
    Choose the most relevant MOI string from rule['moi'] based on variant phenotype texts
    (variant.hgmd_phen and variant.clinvar_trait). If no clear mapping, return rule['moi'].
    """
    if not rule:
        return ""
    moi_field = str(rule.get("moi", "") or "").strip()
    disease_field = str(rule.get("disease", "") or "").strip()
    if not moi_field or not disease_field:
        return moi_field

    # Split disease and moi lists by common separators
    disease_list = [d.strip() for d in re.split(r'[;/\|]', disease_field) if d.strip()]
    moi_list = [m.strip() for m in re.split(r'[;/\|]', moi_field) if m.strip()]
    if not disease_list:
        return moi_field

    # Pad moi_list if shorter
    if len(moi_list) < len(disease_list):
        moi_list = (moi_list + [""] * len(disease_list))[: len(disease_list)]

    # Build variant phenotype string
    traits = []
    if getattr(variant, "hgmd_phen", "") and str(variant.hgmd_phen).lower() not in ("", "not provided", ".", "null"):
        traits.append(str(variant.hgmd_phen))
    if getattr(variant, "clinvar_trait", "") and str(variant.clinvar_trait).lower() not in ("", "not provided", ".", "null"):
        traits.append(str(variant.clinvar_trait))
    variant_trait = " | ".join(traits).strip()

    # Try to match each disease token to variant_trait
    for idx, disease_token in enumerate(disease_list):
        try:
            matched, _, _ = check_disease_match(disease_token, variant_trait, extra_recs_list=None, gene_rules=None)
            if matched:
                return moi_list[idx] if idx < len(moi_list) else moi_list[-1] if moi_list else moi_field
        except Exception:
            continue

    # fallback: if variant._disease_match True, return first non-empty MOI
    try:
        if getattr(variant, "_disease_match", False):
            for m in moi_list:
                if m:
                    return m
    except Exception:
        pass

    return moi_field

def choose_disease_for_gene(rule: Dict[str, Any], variants: List[VariantRecord]) -> Tuple[str, str, str]:
    """
    From a gene-level rule that may contain multiple diseases and MOIs (semicolon/pipe-separated),
    choose the most appropriate (disease, moi) tuple for this sample/gene based on the variants observed.

    Logic:
      - Parse `disease` and `moi` fields into parallel lists (split by ; / | / /).
      - Count the number of P/LP variants in `variants` (automated_class in Pathogenic/Likely pathogenic).
      - If exactly 1 P/LP variant -> prefer a disease whose MOI indicates dominant (contains 'AD' or 'DOMINANT').
        If none found, fall back to first disease.
      - If >= 2 P/LP heterozygous variants -> prefer a disease whose MOI indicates recessive (contains 'AR' or 'RECESSIVE').
        If none found, fall back to first disease.
      - If no P/LP found, return the first disease/moi (fallback).
    Returns:
      (selected_disease, selected_moi, reason)
    """
    if not rule:
        return "", "", "no_rule"

    disease_field = str(rule.get("disease", "") or "").strip()
    moi_field = str(rule.get("moi", "") or "").strip()

    # Split preserving order
    def split_list(s: str) -> List[str]:
        if not s:
            return []
        parts = [p.strip() for p in re.split(r'[;/\|]+', s) if p.strip()]
        return parts

    diseases = split_list(disease_field)
    mois = split_list(moi_field)

    # Pad mois to match diseases
    if len(mois) < len(diseases):
        mois = (mois + [""] * len(diseases))[:len(diseases)]
    if not diseases:
        return "", moi_field, "no_disease_list"

    # Count P/LP variants in this sample/gene
    plp_variants = [v for v in variants if getattr(v, "automated_class", "") in ("Pathogenic", "Likely pathogenic")]
    plp_count = len(plp_variants)

    # Helper to find index matching MOI tokens
    def find_index_with_moi_token(tokens: List[str]) -> Optional[int]:
        for i, m in enumerate(mois):
            ml = (m or "").upper()
            for t in tokens:
                if t in ml:
                    return i
        return None

    # Determine selection
    selected_idx = None
    reason = "fallback_first"

    if plp_count == 1:
        # prefer dominant disease
        idx = find_index_with_moi_token(["AD", "AUTOSOMAL_DOMINANT", "DOMINANT"])
        if idx is not None:
            selected_idx = idx
            reason = "single_plp_prefers_dominant"
    elif plp_count >= 2:
        # prefer recessive disease (AR)
        idx = find_index_with_moi_token(["AR", "AUTOSOMAL_RECESSIVE", "RECESSIVE"])
        if idx is not None:
            selected_idx = idx
            reason = "multi_plp_prefers_recessive"

    # fallback heuristics
    if selected_idx is None:
        # try to pick a disease with explicit AD/AR depending on count
        if plp_count == 1:
            # choose first disease that looks AD-like by moi tokens among mois (looser match)
            for i, m in enumerate(mois):
                if any(tok in (m or "").upper() for tok in ("AD", "DOMINANT")):
                    selected_idx = i; reason = "fallback_dominant_by_moi"
                    break
        elif plp_count >= 2:
            for i, m in enumerate(mois):
                if any(tok in (m or "").upper() for tok in ("AR", "RECESSIVE")):
                    selected_idx = i; reason = "fallback_recessive_by_moi"
                    break

    if selected_idx is None:
        selected_idx = 0
        reason = reason if reason else "final_fallback_first"

    selected_disease = diseases[selected_idx] if selected_idx < len(diseases) else diseases[0]
    selected_moi = mois[selected_idx] if selected_idx < len(mois) else (mois[0] if mois else "")

    return selected_disease, selected_moi, reason

def variant_row_with_disease(v: VariantRecord, rule: Dict[str,Any], family_id: str = "", recs_map: Dict = None, gene_rules: Dict = None) -> dict:
    """
    Build output row for CSV with proper ClinVar phenotype-aware logic.
    
    Key logic:
    - clinvar_phenotype_info: ALL phenotypes from ALL submissions
    - clinvar_trait: ONLY phenotypes matching ANY ACMG disease (split by ;)
    - clinvar_sig: Raw classifications from matching submissions (preserved as is)
    - clinvar_sig_consensus: Consensus classification from matching submissions (Pathogenic/Likely pathogenic/Benign/Conflicting)
    - clinvar_stars: Raw star ratings from matching submissions (preserved as is, e.g., "1;2")
    - clinvar_star_consensus: Aggregated star rating per clinvar_merge.py logic
    - clinvar_phenotype_match: Yes if any submission matches ANY ACMG disease
    """
    # Keep original raw string for display in CSV output
    acmg_disease_raw = rule.get("disease", "") if isinstance(rule, dict) else ""
    
    # Split multiple diseases (separated by semicolon, pipe, or slash) for matching
    acmg_diseases = []
    if acmg_disease_raw:
        for part in re.split(r'[;/\|]+', acmg_disease_raw):
            disease = part.strip()
            if disease:
                acmg_diseases.append(disease)
    
    # =========================================================
    # CLINVAR LOGIC - Phenotype-aware submission processing
    # =========================================================
    
    # Get all submissions for this variant
    all_submissions = getattr(v, "clinvar_submissions", []) or []
    
    # Filter out submissions skipped due to gene mismatch
    valid_submissions = []
    for sub in all_submissions:
        if isinstance(sub, dict) and str(sub.get("note", "")).lower().startswith("gene_mismatch"):
            continue
        valid_submissions.append(sub)
    
    # Split submissions into matching vs non-matching based on phenotype
    matching_subs = []
    non_matching_subs = []
    
    for sub in valid_submissions:
        phen = str(sub.get("phenotype", "")).strip()
        if acmg_diseases and phen:
            matched = False
            for acmg_disease in acmg_diseases:
                try:
                    is_match, _, _ = check_disease_match(acmg_disease, phen)
                    if is_match:
                        matched = True
                        break
                except Exception:
                    continue
            if matched:
                matching_subs.append(sub)
            else:
                non_matching_subs.append(sub)
        else:
            non_matching_subs.append(sub)
    
    # 1. clinvar_phenotype_info: ALL phenotypes from ALL submissions (preserve all, no dedup)
    all_phens = []
    for sub in valid_submissions:
        phen = str(sub.get("phenotype", "")).strip()
        if phen and phen.lower() not in ("not provided", ".", "null"):
            all_phens.append(re.sub(r'^[A-Z0-9]+:', '', phen))
    clinvar_phenotype_info = ";".join(all_phens) if all_phens else "Not provided"
    
    # =========================================================
    # CRITICAL: Check if we have phenotype match
    # =========================================================
    has_phenotype_match = len(matching_subs) > 0
    
    if has_phenotype_match:
        # MODE 1: Use CSV data (phenotype match)
        
        # 2. clinvar_trait: ONLY phenotypes from matching submissions
        match_phens = []
        for sub in matching_subs:
            phen = str(sub.get("phenotype", "")).strip()
            if phen:
                match_phens.append(re.sub(r'^[A-Z0-9]+:', '', phen))
        clinvar_trait = ";".join(match_phens) if match_phens else "Not provided"
        
        # 3. clinvar_sig: Raw classifications from matching submissions
        match_sigs_raw = []
        for sub in matching_subs:
            sig = str(sub.get("clnsig", "")).strip()
            if sig and sig.lower() not in ("not provided", ".", "null"):
                match_sigs_raw.append(sig)
        clinvar_sig_field = ";".join(match_sigs_raw) if match_sigs_raw else "Not provided"
        
        # 4. clinvar_sig_consensus: Consensus classification from matching submissions
        consensus_sig = get_consensus_classification_from_submissions(matching_subs)
        clinvar_sig_consensus = consensus_sig if consensus_sig else "Not provided"
        
        # 5. clinvar_stars: Raw star ratings from matching submissions
        match_stars_raw = []
        for sub in matching_subs:
            stars = str(sub.get("stars", 0))
            if stars and stars != "0":
                match_stars_raw.append(stars)
        clinvar_stars_field = ";".join(match_stars_raw) if match_stars_raw else "0"
        
        # 6. clinvar_star_consensus: Aggregated star rating
        # For conflicting classifications, always use 1 star
        if consensus_sig and "conflict" in consensus_sig.lower():
            star_consensus = 1
        else:
            star_consensus = calculate_consensus_star(matching_subs)
        
    else:
        # MODE 2: No phenotype match - use VCF data (fallback)
        
        # Get VCF aggregate data from variant record
        vcf_cls = getattr(v, "clinvar_sig", "")
        vcf_stars = getattr(v, "clinvar_stars", "0")
        vcf_trait = getattr(v, "clinvar_trait", "")
        
        # 2. clinvar_trait: phenotypes from VCF (CLNDN)
        if vcf_trait and vcf_trait.strip() and vcf_trait.lower() not in ("not provided", ".", "null", ""):
            clinvar_trait = vcf_trait
        else:
            clinvar_trait = "Not provided"
        
        # 3. clinvar_sig: classification from VCF
        if vcf_cls and vcf_cls.strip() and vcf_cls.lower() not in ("not provided", ".", "null", ""):
            clinvar_sig_field = vcf_cls
        else:
            clinvar_sig_field = "Not provided"
        
        # 4. clinvar_sig_consensus: EMPTY for VCF fallback (no consensus from multiple submissions)
        clinvar_sig_consensus = "Not provided"
        
        # 5. clinvar_stars: stars from VCF
        if vcf_stars and vcf_stars != "0":
            clinvar_stars_field = vcf_stars
        else:
            clinvar_stars_field = "0"
        
        # 6. clinvar_star_consensus: EMPTY for VCF fallback
        star_consensus = 0
    
    # 7. clinvar_phenotype_match: Yes if any submission matches ANY ACMG disease
    clinvar_phenotype_match_flag = "Yes" if has_phenotype_match else "No"
     
    # Special_Guidance should be "Yes" only when this gene has additional Rules
    has_guidance = "No"
    try:
        if isinstance(rule, dict) and str(rule.get("rules_text", "") or "").strip():
            has_guidance = "Yes"
    except Exception:
        pass

    # Prepare various textual fields for output
    reasons = "; ".join(v.manual_reasons) if v.manual_reasons else ""
    revel_val = f"{v.revel:.3f}" if v.revel is not None else "Not provided"
    dscore_val = f"{v.dscore:.2f}" if v.dscore is not None else "Not provided"
    splice_val = f"{v.spliceai:.3f}" if v.spliceai is not None else "Not provided"
    alpha_val = f"{v.alpha_missense:.3f}" if v.alpha_missense is not None else "Not provided"
    
    # Determine zygosity (Het/Hom/Hemi)
    zyg = "Het"
    if _is_hom(v.proband_gt):
        zyg = "Hom"
        if v.sample_sex == "Male" and ("X" in v.chrom or "Y" in v.chrom):
            zyg = "Hemi"

    # Build informative ClinVar string for display (includes signature, stars, trait)
    clin_parts = []
    if clinvar_sig_field and clinvar_sig_field != "Not provided":
        clin_parts.append(f"Sig={clinvar_sig_field}")
    if clinvar_stars_field and clinvar_stars_field != "0":
        clin_parts.append(f"Stars={clinvar_stars_field}")
    if clinvar_trait and clinvar_trait != "Not provided":
        clin_parts.append(f"Trait='{clinvar_trait}'")
    clin_block_display = getattr(v, "_clin_block_raw_display", "") or ""
    if clin_block_display:
        clin_parts.append(f"ClinBlock=({clin_block_display})")
    clinvar_info = "; ".join(clin_parts) if clin_parts else "Not provided"

    # Build informative HGMD string (class, publications, support, phenotype)
    hgmd_parts = []
    hgmd_class = getattr(v, "hgmd_class", "") or ""
    hgmd_pubs = getattr(v, "hgmd_publication_count", 0) or 0
    hgmd_dmsup = getattr(v, "hgmd_dmsupported", 0) or 0
    hgmd_rank = getattr(v, "hgmd_rankscore", "") or ""
    hgmd_phen = getattr(v, "hgmd_phen", "") or ""
    if hgmd_class:
        hgmd_parts.append(f"Class={hgmd_class}")
    if hgmd_pubs:
        hgmd_parts.append(f"Pubs={hgmd_pubs}")
    if hgmd_dmsup:
        hgmd_parts.append(f"dmsupported={hgmd_dmsup}")
    if hgmd_rank and str(hgmd_rank).strip().upper() not in ("", "NULL"):
        hgmd_parts.append(f"Rank={hgmd_rank}")
    if hgmd_phen:
        hgmd_parts.append(f"Phen='{hgmd_phen}'")
    hgmd_block_display = getattr(v, "_hgmd_block_raw_display", "") or ""
    if hgmd_block_display:
        hgmd_parts.append(f"HGMDBlock=({hgmd_block_display})")
    hgmd_info = "; ".join(hgmd_parts) if hgmd_parts else "Not provided"

    # Build Rules_applied string
    raw_rules = getattr(v, "_rules_applied", []) or []
    
    # CRITICAL FIX: Filter out false positive rules for genes that don't have them
    gene_name = getattr(v, "gene", "")
    if gene_name:
        gene_norm = normalize_gene_name(gene_name)
        # Get actual rules_text for this gene
        actual_rule = gene_rules.get(gene_norm, {}) if gene_rules else {}
        actual_rules_text = str(actual_rule.get("rules_text", "") or "").strip()
         
        # If gene has NO rules_text, clear _rules_applied
        if not actual_rules_text:
            raw_rules = []
            v._rules_applied = []
            logger.debug(f"Cleared _rules_applied for {gene_name} - no rules_text in ACMG table")
    
    unique_rules = list(dict.fromkeys(raw_rules))
    rules_applied_str = ";".join(unique_rules)

    # HGMD phenotype field with deduplication
    hgmd_ph = getattr(v, "hgmd_phen", "") or ""
    hgmd_phenotype_info = ";".join(_uniq_preserve([p.strip() for p in (hgmd_ph or "").split(";") if p.strip()])) if hgmd_ph else "Not provided"

    # Build database summary from ClinVar and HGMD evidence
    db_sum_parts = []
    if clinvar_sig_field and clinvar_sig_field != "Not provided":
        db_sum_parts.append(f"ClinVar: Sig={clinvar_sig_field}")
    if clinvar_trait and clinvar_trait != "Not provided":
        db_sum_parts.append(f"ClinVar: Trait={clinvar_trait}")
    if hgmd_ph and hgmd_ph != "":
        db_sum_parts.append(f"HGMD: Phen={hgmd_phenotype_info}")
    if v.internal_db_sig:
        db_sum_parts.append(f"Internal: {v.internal_db_sig}")
    db_summary = "; ".join(db_sum_parts) if db_sum_parts else "No database evidence"

    # gnomAD version and allele frequency info
    gnomad_ver = v.gnomad_details.get("version", "Unknown") if v.gnomad_details else "Unknown"
    in_auto = getattr(v, '_in_auto', False)
    in_manual = getattr(v, '_in_manual', False)
    filtered_flag = getattr(v, '_filtered', True)
    sample_af_val = f"{v.proband_af:.3f}" if v.proband_af is not None else "0.0"

    # Choose effective MOI for this specific variant (disease-specific mapping)
    try:
        # Pass the raw disease string (with possible multiple diseases)
        moi = select_moi_for_variant(rule if isinstance(rule, dict) else {}, v)
        if not moi:
            moi = rule.get("moi", "Unknown") if isinstance(rule, dict) else "Unknown"
    except Exception:
        moi = rule.get("moi", "Unknown") if isinstance(rule, dict) else "Unknown"
        
    # Format PS1 and PM5 match strings for output
    ps1_str = ";".join([f"{m.get('hgvsp','')}/{m.get('clnsig','')}/{m.get('source','')}" for m in v.ps1_matches]) if v.ps1_matches else ""
    pm5_str = ";".join([f"{m.get('hgvsp','')}/{m.get('best','')}" for m in v.pm5_matches]) if v.pm5_matches else ""

    # HGMD DM support value
    dmsupported_val = getattr(v, "hgmd_dmsupported", 0) or 0

    # In variant_row_with_disease, when building comments_field:
    comments_field = ""

    # Collect comments from various sources
    comment_parts = []

    # Add rule-related comments if present
    if hasattr(v, "_comments_parts") and v._comments_parts:
        comment_parts.extend(v._comments_parts)

    # Add any other comments from variant.comments
    if hasattr(v, "comments") and v.comments:
        comment_parts.append(v.comments)

    # Join all comment parts
    if comment_parts:
        comments_field = "; ".join(comment_parts)

    # Build the final output row dictionary
    output_row = {
        "Query_Variant": v.query_variant,
        "sample": v.sample,
        "Family_ID": family_id if family_id else v.sample,
        "chrom": v.chrom,
        "pos": v.pos,
        "ref": v.ref,
        "alt": v.alt,
        "gene": v.gene,
        "Inheritance": moi,
        "HGVSc": v.hgvsc or "Not provided",
        "HGVSp": v.hgvsp or "Not provided",
        "AA_ref": v.aa_ref or "",
        "AA_alt": v.aa_alt or "",
        "AA_pos": v.aa_pos or "",
        "tid": v.tid,
        "Exon": v.exon,
        "Consequence": v.consequence,
        "Impact": v._raw_ann.get("IMPACT", "") if getattr(v, "_raw_ann", None) else getattr(v, "impact", ""),
        "Zygosity": zyg,
        "DP": v.proband_dp if v.proband_dp is not None else "0",
        "Sample_AF": sample_af_val,
        "acmg_disease_association": acmg_disease_raw,
        "Special_Guidance": has_guidance,
        "Rules_applied": rules_applied_str,
        # ClinVar phenotype-aware fields
        "clinvar_phenotype_match": clinvar_phenotype_match_flag,
        "clinvar_phenotype_match_level": getattr(v, "_clinvar_phenotype_match_level", "") or "",
        "clinvar_phenotype_match_reason": getattr(v, "_clinvar_phenotype_match_reason", "") or "",
        "clinvar_phenotype_info": clinvar_phenotype_info,
        "hgmd_phenotype_info": hgmd_phenotype_info,
        "hgmd_phenotype_match": "Yes" if getattr(v, "_hgmd_phenotype_match", False) else "No",
        "hgmd_phenotype_match_level": getattr(v, "_hgmd_phenotype_match_level", "") or "",
        "hgmd_phenotype_match_reason": getattr(v, "_hgmd_phenotype_match_reason", "") or "",
        "clinvar_trait": clinvar_trait,  
        "clinvar_sig": clinvar_sig_field,
        "clinvar_sig_consensus": clinvar_sig_consensus, 
        "clinvar_stars": clinvar_stars_field,
        "clinvar_star_consensus": star_consensus,
        "hgmd_class": v.hgmd_class or "Not provided",
        "hgmd_phenotype": hgmd_phen or "Not provided",
        "hgmd_publication_count": v.hgmd_publication_count,
        "hgmd_dmsupported": dmsupported_val,
        "hgmd_rankscore": getattr(v, 'hgmd_rankscore', "Not provided"),
        "internal_db_sig": v.internal_db_sig or "Not provided",
        "database_summary": db_summary,
        "gnomAD_AF": f"{max(v.gnomad_af, 0.0):.6f}" if v.gnomad_af is not None and v.gnomad_af > 0 else "0.0",
        "gnomAD_version": gnomad_ver,
        "REVEL": revel_val,
        "DSCORE": dscore_val,
        "SpliceAI": splice_val,
        "SQUIRLS": getattr(v, "squirls_score", "") if getattr(v, "squirls_score", None) is not None else getattr(v, "squirls_category", ""),
        "AlphaMissense": alpha_val,
        "ps1_matches": ps1_str,
        "pm5_matches": pm5_str,
        "criteria_points_ext": json.dumps(getattr(v, "criteria_points_ext", {}) or {}, ensure_ascii=False),
        "total_ext_points": getattr(v, "total_ext_points", ""),
        "extended_class": getattr(v, "extended_class", ""),
        "PP5ext_BP5ext_expl": getattr(v, "_pp5ext_explanation", ""),
        "PP5_clin_block": getattr(v, "_clin_block_details", getattr(v, "_clin_block_raw_display", "")),
        "PP5_hgmd_block": getattr(v, "_hgmd_block_details", getattr(v, "_hgmd_block_raw_display", "")),
        "PP5_raw_total": f"{getattr(v, '_pp5_raw_total', 0.0):.3f}" if getattr(v, "_pp5_raw_total", None) is not None else "",
        "PP5_applied": getattr(v, "_applied_ext_key", ""),
        "phase": v.phase,
        "is_delins_part": "Yes" if v.is_delins_part else "No",
        "in_auto_report": "Yes" if in_auto else "No",
        "in_manual_review": "Yes" if in_manual else "No",
        "auto_report_reasons": ";".join(getattr(v, "auto_report_reasons", [])),
        "manual_reasons": reasons,
        "filtered_reasons": ";".join(getattr(v, "filtered_reasons", [])),
        "filtered_out": "Yes" if filtered_flag else "No",
        "comments": comments_field
    }

    # Add per-criterion explanation columns to avoid KeyError during CSV reindex
    for crit in ["PVS1","PS1","PM1","PM2","PM3","PM4","PM5","PP1","PP3","PP5","BA1","BS1","BS2","BP4","BP5","BP6","BP7"]:
        output_row[f"{crit}_expl"] = v.criteria_explanations.get(crit, "")

    return output_row

def get_consensus_classification_from_submissions(submissions: List[Dict]) -> str:
    """
    Calculate consensus classification from a list of submissions.
    Follows the same logic as clinvar_merge.py:
    - If conflicting vectors -> "Conflicting classifications of pathogenicity"
    - If pathogenic vector -> "Pathogenic" or "Likely pathogenic"
    - If benign vector -> "Benign" or "Likely benign"
    - If VUS -> "Uncertain significance"
    """
    if not submissions:
        return "Not provided"
    
    # Collect unique normalized classifications
    sigs = []
    for sub in submissions:
        sig = str(sub.get("clnsig", "")).strip()
        if sig and sig.lower() not in ("not provided", ".", "null"):
            sigs.append(normalize_clinvar_label(sig))
    
    if not sigs:
        return "Not provided"
    
    # First, check if any submission is explicitly CONFLICT
    # This must be checked BEFORE pathogenic/benign because "Conflicting classifications of pathogenicity" contains "pathogenic"
    has_conflict = any('conflict' in s.lower() for s in sigs)
    
    if has_conflict:
        # Check if there are also definitive pathogenic or benign labels
        has_path = any(s.lower() in ("pathogenic", "likely pathogenic") for s in sigs)
        has_benign = any(s.lower() in ("benign", "likely benign") for s in sigs)
        
        if has_path and has_benign:
            return "Conflicting classifications of pathogenicity"
        elif has_path:
            # Conflict + Pathogenic -> still Conflict (conservative)
            return "Conflicting classifications of pathogenicity"
        elif has_benign:
            # Conflict + Benign -> still Conflict (conservative)
            return "Conflicting classifications of pathogenicity"
        else:
            # Only conflict labels
            return "Conflicting classifications of pathogenicity"
    
    # Now check vectors for non-conflict submissions
    sigs_lower = [s.lower() for s in sigs]
    
    # Check vectors (safe because we already filtered out conflicts)
    has_path = any(s in ("pathogenic", "likely pathogenic") for s in sigs_lower)
    has_benign = any(s in ("benign", "likely benign") for s in sigs_lower)
    has_vus = any('uncertain' in s or 'vus' in s for s in sigs_lower)
    
    categories = sum([has_path, has_benign, has_vus])
    
    if categories > 1:
        return "Conflicting classifications of pathogenicity"
    
    if has_path:
        # If at least one strict Pathogenic (not likely) -> Pathogenic, else Likely pathogenic
        if any(s == "pathogenic" for s in sigs_lower):
            return "Pathogenic"
        return "Likely pathogenic"
    
    if has_benign:
        if any(s == "benign" for s in sigs_lower):
            return "Benign"
        return "Likely benign"
    
    if has_vus:
        return "Uncertain significance"
    
    # Fallback
    return sigs[0] if sigs else "Not provided"

def collapse_output_rows(rows: List[Dict[str,Any]]) -> List[Dict[str,Any]]:
    out = {}
    for r in rows:
        norm_gene = normalize_gene_name(r.get("gene") or r.get("normalized_gene","") or "")
        key = f"{r.get('sample','')}-{r.get('chrom','')}-{r.get('pos','')}-{r.get('ref','')}-{r.get('alt','')}-{norm_gene}"
        if key not in out:
            nr = dict(r)
            try:
                if isinstance(nr.get("criteria_points"), str) and nr.get("criteria_points"):
                    nr["criteria_points"] = json.loads(nr["criteria_points"])
            except Exception:
                pass
            out[key] = nr
            continue
        existing = out[key]
        mr_existing = set([x for x in str(existing.get("manual_reasons","")).split(";") if x])
        mr_new = set([x for x in str(r.get("manual_reasons","")).split(";") if x])
        merged_mr = ";".join(sorted(mr_existing.union(mr_new)))
        existing["manual_reasons"] = merged_mr
        ca_existing = set([x for x in str(existing.get("criteria_assigned","")).split(";") if x])
        ca_new = set([x for x in str(r.get("criteria_assigned","")).split(";") if x])
        merged_ca = ";".join(sorted(ca_existing.union(ca_new)))
        existing["criteria_assigned"] = merged_ca
        try:
            cp_existing = existing.get("criteria_points") or {}
            if isinstance(cp_existing, str) and cp_existing:
                cp_existing = json.loads(cp_existing)
        except Exception:
            cp_existing = {}
        try:
            cp_new = r.get("criteria_points") or {}
            if isinstance(cp_new, str) and cp_new:
                cp_new = json.loads(cp_new)
        except Exception:
            cp_new = {}
        for k, v in cp_new.items():
        # parse incoming value to int
            try:
                vi = int(v)
            except Exception:
                try:
                    vi = int(float(v))
                except Exception:
                    vi = 0

            if k not in cp_existing:
                cp_existing[k] = vi
            else:
                try:
                    existing_val = int(cp_existing.get(k, 0))
                except Exception:
                    try:
                        existing_val = int(float(cp_existing.get(k, 0)))
                    except Exception:
                        existing_val = 0

                # If either value is negative, keep the most negative (preserve penalty)
                if existing_val < 0 or vi < 0:
                    cp_existing[k] = min(existing_val, vi)
                else:
                    # both non-negative: keep the maximum (avoid double-counting of same criterion)
                    cp_existing[k] = max(existing_val, vi)

        existing["criteria_points"] = cp_existing
        try:
            existing["total_points"] = int(sum(int(x) for x in existing["criteria_points"].values()))
        except Exception:
            pass
        order = {"Pathogenic":4, "Likely pathogenic":3, "VUS":2, "Likely benign":1, "Benign":0}
        a1 = existing.get("automated_class","VUS"); a2 = r.get("automated_class","VUS")
        if order.get(a2,2) > order.get(a1,2):
            existing["automated_class"] = a2
    final = []
    for k,v in out.items():
        v2 = dict(v)
        try:
            v2["criteria_points"] = json.dumps(v2.get("criteria_points") or {}, ensure_ascii=False)
        except Exception:
            v2["criteria_points"] = "{}"
        final.append(v2)
    return final

# ---------------------------
# strict_triage_and_output (integrates positional checks)
# ---------------------------
def strict_triage_and_output(
    report_candidates: List[VariantRecord],
    audit_candidates: List[VariantRecord],
    gene_rules: Dict[str, Dict[str, Any]],
    recs_map: Dict[str, Dict[str, Any]],
    outdir: str,
    family_id_input: str = "Unknown",
    *,
    dbm: Optional["DatabaseManager"] = None,
    # cooccur_path: Optional[str] = None  # Handled by external script
) -> Tuple[str, str, str, int, int]:
    """
    Triage and output generation for ACMG secondary findings.
    
    This function implements the conservative auto/manual decision pipeline:
      - High-quality database-driven auto-reporting (ClinVar≥2★, HGMD DM with support).
      - Algorithmic auto-reporting for high-confidence pathogenic variants.
      - AR-gene carrier suppression (single het P/LP without partner).
      - AR pairing logic with phasing/co-occurrence evidence.
      - X-linked male hemizygosity checks and AF filtering.
      - Gene-specific rule evaluation (exclusions, missense-only, positional constraints).
      - Delins/cluster detection and escalation to manual review.
    
    Outputs:
      - auto_conclusions.csv: Variants recommended for automatic reporting.
      - manual_review_list.csv: Variants requiring curator review.
      - all_candidates.csv: Full annotated variant list with decision metadata.
    
    Args:
        report_candidates (List[VariantRecord]): Variants passing initial QC.
        audit_candidates (List[VariantRecord]): Variants failing QC but retained for audit.
        gene_rules (Dict): ACMG SF gene table rules.
        recs_map (Dict): Additional gene-specific recommendations (positional rules).
        outdir (str): Output directory for CSV files.
        family_id_input (str): Family identifier for reporting.
        dbm (DatabaseManager): Optional database manager for co-occurrence lookup.
        cooccur_path (str): Optional path to gnomAD co-occurrence TSV.
    
    Returns:
        Tuple[str, str, str, int, int]: Paths to output files and counts.
    """

    os.makedirs(outdir, exist_ok=True)
    family_id = family_id_input or "Unknown"
    candidates = list(report_candidates or []) + list(audit_candidates or [])

    # Diagnostic logging: show counts to help debugging when pipeline runs in batch mode.
    logger.info("strict_triage_and_output: report_candidates=%d, audit_candidates=%d, total_candidates=%d",
                len(report_candidates or []), len(audit_candidates or []), len(candidates))
    auto_ids = set() 
    manual_ids = set()

    # ensure DB manager available (prefer passed dbm; fallback to global DBM if present)
    dbm_local = dbm or globals().get("DBM") or globals().get("dbm") or None

    def get_key(v: VariantRecord) -> str:
        return f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.sample}"

    # Safe allele fraction extraction
    def _get_sample_af(v: VariantRecord) -> float:
        try:
            if getattr(v, "proband_af", None) is not None:
                return float(v.proband_af)
        except Exception:
            pass
        try:
            ad = getattr(v, "proband_ad", None)
            if ad and isinstance(ad, (list, tuple)) and len(ad) >= 2:
                ref = float(ad[0] or 0)
                alt = float(ad[1] or 0)
                tot = ref + alt
                if tot > 0:
                    return alt / tot
        except Exception:
            pass
        # fallback to 0.0
        return 0.0

    # Read cooccurrence TSV into dict for fast lookup (load once)
    coocc_dict = {}
    # if cooccur_path and os.path.exists(cooccur_path):
    #     try:
    #         with open(cooccur_path, "rt", encoding="utf-8") as fh:
    #             for ln in fh:
    #                 if not ln.strip() or ln.startswith("#"):
    #                     continue
    #                 parts = ln.rstrip("\n").split("\t")
    #                 if len(parts) < 3:
    #                     continue
    #                 a = parts[0].strip()
    #                 b = parts[1].strip()
    #                 cis_flag = parts[2].strip().lower()
    #                 evidence = parts[3].strip() if len(parts) > 3 else ""
    #                 key1 = f"{a}|{b}"
    #                 key2 = f"{b}|{a}"
    #                 cis_val = None
    #                 if cis_flag in ("yes", "y", "true", "1"):
    #                     cis_val = True
    #                 elif cis_flag in ("no", "n", "false", "0"):
    #                     cis_val = False
    #                 coocc_dict[key1] = {"cis": cis_val, "evidence": evidence, "source": "cooccurrence_tsv"}
    #                 coocc_dict[key2] = coocc_dict[key1]
    #     except Exception:
    #         # if loading fails, keep dict empty (function remains conservative)
    #         coocc_dict = {}
    coocc_dict = {}  # Always empty - cooccurrence handled by external script

    # Check cis/trans evidence from recs_map or cooccurrence dict
    def _check_pair_cis_evidence(v1: VariantRecord, v2: VariantRecord) -> Optional[dict]:
        vid1 = f"{v1.chrom}:{v1.pos}:{v1.ref}:{v1.alt}"
        vid2 = f"{v2.chrom}:{v2.pos}:{v2.ref}:{v2.alt}"
        # 1) recs_map per gene (if exists)
        try:
            rec = recs_map.get(v1.gene, {}) if recs_map else {}
            cis_map = rec.get("cis_pairs") or rec.get("cis_evidence") or {}
            if cis_map:
                k1 = f"{vid1}|{vid2}"
                k2 = f"{vid2}|{vid1}"
                if k1 in cis_map:
                    entry = cis_map[k1]
                    return {"cis": bool(entry.get("cis")), "source": "recs_map", "evidence": entry.get("evidence", "")}
                if k2 in cis_map:
                    entry = cis_map[k2]
                    return {"cis": bool(entry.get("cis")), "source": "recs_map", "evidence": entry.get("evidence", "")}
        except Exception:
            pass
        # 2) cooccurrence analysis handled by external script - return None
        # key = f"{vid1}|{vid2}"
        # if key in coocc_dict:
        #     return coocc_dict[key]
        # cooccurrence analysis handled by external script
        return None
    
    # --- Helpers to parse & evaluate gene Rules column for a variant ---
    def _split_rules_text(raw: Optional[str]) -> List[str]:
        """Split raw Rules cell text into non-empty rule lines."""
        if not raw:
            return []
        # Many cells may contain newlines or repeated identical lines; split by newline and keep non-empty
        lines = []
        for ln in re.split(r'[\r\n]+', str(raw)):
            s = ln.strip()
            if s:
                lines.append(s)
        return lines

    def _evaluate_single_rule(rule: str, vr: VariantRecord) -> bool:
        """
        Return True if this single rule matches the variant (meaning it triggers the rule).
        Supports patterns seen in ACMG_SF_v3.3.tsv:
        - Include/Exclude if Impact is "MODERATE" or "LOW" (or unquoted)
        - Exclude if Consequence contains "5_prime_UTR_variant"
        - Include/Exclude if HGVSp contains "p.Cys282Tyr"
        - Exclude if HGVSp contains any of: ["p.X","p.Y",...]
        """
        if not rule or not vr:
            return False
        r = rule.strip()

        # --- Extract/derive IMPACT ---
        try:
            impact_val = ""
            if hasattr(vr, "impact") and (getattr(vr, "impact") or "") != "":
                impact_val = str(getattr(vr, "impact"))
            else:
                ann = getattr(vr, "_raw_ann", {}) or {}
                impact_val = str(ann.get("IMPACT") or ann.get("Impact") or ann.get("impact") or "").strip()
            impact_val_up = impact_val.upper()
        except Exception:
            impact_val_up = ""

        # If IMPACT absent, derive from consequence (best-effort)
        if not impact_val_up:
            cons = str(getattr(vr, "consequence", "") or "").lower()
            if any(tok in cons for tok in ("stop_gained", "frameshift", "start_lost", "stop_lost", "splice_acceptor", "splice_donor")):
                impact_val_up = "HIGH"
            elif "missense" in cons or "inframe" in cons:
                impact_val_up = "MODERATE"
            elif "synonymous" in cons or "5_prime_utr_variant" in cons or "3_prime_utr_variant" in cons:
                impact_val_up = "LOW"
            else:
                impact_val_up = ""

        # --- 1) Impact rules - FIXED to handle "or" with spaces properly ---
        # Pattern matches: "MODERATE" or "LOW", "HIGH" or "MODERATE", etc.
        # Allows spaces around "or": "MODERATE or LOW", "MODERATE or  LOW", etc.
        
        # First try to match pattern with quoted values and optional "or"
        # This regex captures: first value, optional "or" with second value
        # Example: Impact is "MODERATE" or "LOW" -> group1="MODERATE", group2="LOW"
        m = re.search(r'impact\s+is\s+(?:"([^"]+)"|\'([^\']+)\'|([A-Za-z_]+))(?:\s+or\s+(?:"([^"]+)"|\'([^\']+)\'|([A-Za-z_]+)))?', r, flags=re.I)
        if m:
            toks = []
            # Collect first value (from any of the possible quote styles)
            for g in (m.group(1), m.group(2), m.group(3)):
                if g:
                    toks.append(str(g).strip().upper())
                    break
            # Collect second value if present (after "or")
            for g in (m.group(4), m.group(5), m.group(6)):
                if g:
                    toks.append(str(g).strip().upper())
                    break
            
            if toks and impact_val_up:
                # Check if impact matches ANY of the allowed values
                for t in toks:
                    if t and impact_val_up == t:  # Exact match, not substring
                        return True
            return False

        # Fallback: handle "Impact is MODERATE or LOW" without quotes and with spaces
        # This regex captures the entire list after "is" and splits by "or"
        m_alt = re.search(r'impact\s+is\s+([A-Za-z_]+(?:\s+or\s+[A-Za-z_]+)*)', r, flags=re.I)
        if m_alt:
            values_part = m_alt.group(1).strip()
            # Split by "or" (case insensitive)
            allowed_values = re.split(r'\s+or\s+', values_part, flags=re.I)
            allowed_values = [v.strip().upper() for v in allowed_values]
            if impact_val_up in allowed_values:
                return True
            return False

        # also support simple "Impact is MODERATE" without quotes (caught by above, but keep fallback)
        m2 = re.search(r'\bimpact\s+is\s+([A-Za-z_]+)\b', r, flags=re.I)
        if m2 and '"' not in r and "'" not in r:
            desired = m2.group(1).strip().upper()
            if desired and desired == impact_val_up:
                return True

        # --- 2) Consequence contains ... ---
        m3 = re.search(r'consequence\s+contains\s*["\']?([^"\']+)["\']?', r, flags=re.I)
        if m3:
            want = m3.group(1).strip().lower()
            cons = str(getattr(vr, "consequence", "") or "").lower()
            if want and want in cons:
                return True
            return False

        # --- 3) HGVSp contains "..." ---
        m4 = re.search(r'hgvsp\s+contains\s*["\']([^"\']+)["\']', r, flags=re.I)
        if m4:
            want = m4.group(1).strip()
            hgvsp = str(getattr(vr, "hgvsp", "") or "")
            if want and want in hgvsp:
                return True
            return False

        # --- 4) HGVSp contains any of: ["p.Arg...","p.X..."] ---
        m5 = re.search(r'hgvsp\s+contains\s+any\s+of\s*:\s*(\[[^\]]+\])', r, flags=re.I)
        if m5:
            list_text = m5.group(1)
            try:
                arr = ast.literal_eval(list_text)
                if isinstance(arr, (list, tuple)):
                    hgvsp = str(getattr(vr, "hgvsp", "") or "")
                    for item in arr:
                        if item and str(item) in hgvsp:
                            return True
            except Exception:
                # fallback: parse quoted substrings
                items = re.findall(r'["\']([^"\']+)["\']', list_text)
                hgvsp = str(getattr(vr, "hgvsp", "") or "")
                for it in items:
                    if it and it in hgvsp:
                        return True
            return False

        # No known pattern matched
        return False

    def evaluate_gene_rules_for_variant(rules_text: Optional[str], vr: VariantRecord) -> Tuple[str, List[str]]:
        """
        Evaluate gene-specific rule block against a variant and return decision + matched rule lines.
        - Returns ("exclude" | "include" | "no_match", matched_rules_list)
        """
        if not rules_text:
            return "no_match", []

        lines = _split_rules_text(rules_text)
        matched_excludes = []
        matched_includes = []
        include_lines_present = False
        include_lines_collected = []

        for line in lines:
            line_strip = line.strip()
            if not line_strip:
                continue

            is_excl = bool(re.match(r'^\s*Exclude\s+if', line_strip, flags=re.I))
            is_incl = bool(re.match(r'^\s*Include\s+if', line_strip, flags=re.I))

            if is_incl:
                include_lines_present = True
                include_lines_collected.append(line_strip)

            triggers = _evaluate_single_rule(line_strip, vr)
            if triggers:
                if is_excl:
                    matched_excludes.append(line_strip)
                elif is_incl:
                    matched_includes.append(line_strip)
                else:
                    if 'exclude' in line_strip.lower():
                        matched_excludes.append(line_strip)
                    elif 'include' in line_strip.lower():
                        matched_includes.append(line_strip)
                    else:
                        matched_includes.append(line_strip)

        # Exclude rules have precedence
        if matched_excludes:
            return "exclude", matched_excludes

        # If at least one include matched
        if matched_includes:
            return "include", matched_includes

        # If include lines present but none matched - do NOT filter
        if include_lines_present and not matched_includes:
            explanation = "No include rule matched"
            return "no_match", [explanation] + include_lines_collected

        return "no_match", []

    # Initialize bookkeeping fields on VariantRecord items
    for v in report_candidates:
        v._qc_passed = True
    qc_passed_keys = {get_key(v) for v in report_candidates}
    for v in audit_candidates:
        if get_key(v) in qc_passed_keys:
            v._qc_passed = True
        else:
            v._qc_passed = False
            if not hasattr(v, "filtered_reasons"):
                v.filtered_reasons = []
            if "Failed strict QC" not in v.filtered_reasons:
                v.filtered_reasons.append("Failed strict QC")

    # Group variants
    all_variants = report_candidates + audit_candidates
    # Make sure all variants have common internal attributes used later
    for v in all_variants:
        # boolean/flags
        if not hasattr(v, "_qc_passed"):
            v._qc_passed = True
        if not hasattr(v, "_in_auto"):
            v._in_auto = False
        if not hasattr(v, "_in_manual"):
            v._in_manual = False
        if not hasattr(v, "_filtered"):
            v._filtered = True
        # disease matching meta
        if not hasattr(v, "_disease_match"):
            v._disease_match = False
        if not hasattr(v, "_disease_match_level"):
            v._disease_match_level = "N/A"
        if not hasattr(v, "_disease_match_reason"):
            v._disease_match_reason = "N/A"
        # db evidence flag (used in auto-decision)
        if not hasattr(v, "_has_db_evidence"):
            v._has_db_evidence = False
        # ensure lists/dicts exist
        if not hasattr(v, "manual_reasons") or v.manual_reasons is None:
            v.manual_reasons = []
        if not hasattr(v, "auto_report_reasons") or v.auto_report_reasons is None:
            v.auto_report_reasons = []
        if not hasattr(v, "filtered_reasons") or v.filtered_reasons is None:
            v.filtered_reasons = []
        if not hasattr(v, "criteria_explanations") or v.criteria_explanations is None:
            v.criteria_explanations = {}
        if not hasattr(v, "ps1_matches") or v.ps1_matches is None:
            v.ps1_matches = []
        if not hasattr(v, "pm5_matches") or v.pm5_matches is None:
            v.pm5_matches = []
        if not hasattr(v, "comments") or v.comments is None:
          v.comments = ""
    for v in all_variants:
        try:
            if getattr(v, "is_delins_part", False):
                note = "Manual: variant flagged as part of delins cluster near pathogenic variant -> manual review required"
                if note not in v.manual_reasons:
                    v.manual_reasons.append(note)
                v._in_manual = True
                v._in_auto = False
                v._filtered = False
                manual_ids.add(f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.sample}")
        except Exception:
            # Do not fail triage on unexpected errors here; just continue
            continue
    acmg_reportable_genes = set()
    for gk, gv in (gene_rules or {}).items():
        try:
            gn = normalize_gene_name(gk)
            if isinstance(gv, dict) and gv.get("reportable", False):
                acmg_reportable_genes.add(gn)
        except Exception:
            continue
    by_sample_gene = defaultdict(lambda: defaultdict(list))
    for v in all_variants:
        # normalize gene name for lookup safety
        gene_key = normalize_gene_name(getattr(v, "gene", "") or "")
        # If a set of reportable genes is defined, skip variants whose genes are not in that set
        if acmg_reportable_genes and gene_key not in acmg_reportable_genes:
            # skip non-reportable genes
            continue
        # Ensure variant fields exist (avoid attribute errors later)
        if not hasattr(v, "manual_reasons"):
            v.manual_reasons = []
        if not hasattr(v, "auto_report_reasons"):
            v.auto_report_reasons = []
        if not hasattr(v, "filtered_reasons"):
            v.filtered_reasons = []
        if not hasattr(v, "criteria_explanations"):
            v.criteria_explanations = {}
        if not hasattr(v, "ps1_matches"):
            v.ps1_matches = []
        if not hasattr(v, "pm5_matches"):
            v.pm5_matches = []
        # Use normalized gene key as grouping key to match gene_rules map expectations
        by_sample_gene[v.sample][gene_key].append(v)
    
    # Collect rows for outputs
    auto_rows = []
    manual_rows = []
    all_rows = []

    # Criteria list for expl columns
    BASE_CRITERIA = ["PVS1","PS1","PM1","PM2","PM3","PM4","PM5","PP1","PP3","PP5",
                     "BA1","BS1","BS2","BP4","BP5","BP6","BP7"]

    # Iterate samples/genes
    for sample, genes in by_sample_gene.items():
        for gene, varlist in genes.items():
            filtered_varlist = []
            
            for v in varlist:
                try:
                    total_pts = int(getattr(v, "total_ext_points", getattr(v, "total_points", 0)) or 0)
                    ext_cls = getattr(v, "extended_class", getattr(v, "automated_class", "VUS"))
                    
                    # Filtered: VUS < 3 pts
                    if ext_cls == "VUS" and total_pts < 3:
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        v.is_delins_part = False
                        if not hasattr(v, "filtered_reasons") or v.filtered_reasons is None:
                            v.filtered_reasons = []
                        reason = f"Filtered: VUS with {total_pts} points (<3) - not eligible for reporting"
                        if reason not in v.filtered_reasons:
                            v.filtered_reasons.append(reason)
                        continue
                    
                    selected_moi = getattr(v, "_selected_moi", None)
                    if not selected_moi:
                        rule = gene_rules.get(gene, {})
                        selected_moi = str(rule.get("moi", "") or "").upper()
                    
                    #  # Filtered: AD VUS 
                    is_ad = any(x in selected_moi.upper() for x in ["AD", "DOMINANT", "AUTOSOMAL_DOMINANT", "XLD"])
                    
                    if ext_cls == "VUS" and is_ad and total_pts < 3:
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        if not hasattr(v, "filtered_reasons") or v.filtered_reasons is None:
                            v.filtered_reasons = []
                        reason = f"Filtered: AD inheritance with VUS <3 points (class={ext_cls}, points={total_pts})"
                        if reason not in v.filtered_reasons:
                            v.filtered_reasons.append(reason)
                        continue
                    
                    filtered_varlist.append(v)
                    
                except Exception as e:
                    logger.debug(f"Error filtering variant {getattr(v, 'query_variant', '?')}: {e}")
                    filtered_varlist.append(v)  
            
            varlist[:] = filtered_varlist
            
            if not varlist:
                continue
            
            # Load gene rule and recommendations for this gene
            rule = gene_rules.get(gene, {}) if gene_rules else {}
            spec_rules = recs_map.get(gene, {}) if recs_map else {}

             # --- Evaluate per-variant Rules text (record which rule lines matched) ---
            try:
                rules_text = ""
                if isinstance(rule, dict):
                    rules_text = str(rule.get("rules_text") or rule.get("Rules") or "").strip()
                for v in varlist:
                    try:
                        decision, matched_rules = evaluate_gene_rules_for_variant(rules_text, v)
                        v._rules_applied = matched_rules or []
                        v._rules_decision = decision or ""
                        if matched_rules:
                            note = f"Rules_applied: {', '.join(matched_rules)}"
                            if note not in (v.comments or ""):
                                if v.comments and v.comments.strip():
                                    v.comments = v.comments.strip() + "; " + note
                                else:
                                    v.comments = note
                    except Exception:
                        v._rules_applied = []
                        v._rules_decision = ""
            except Exception:
                for v in varlist:
                    v._rules_applied = []
                    v._rules_decision = ""

            # Determine sample sex
            sample_sex = varlist[0].sample_sex if varlist and hasattr(varlist[0], "sample_sex") else "Unknown"

            # Choose disease / MOI for this sample/gene based on variant counts
            selected_disease, selected_moi, disease_choice_reason = choose_disease_for_gene(rule, varlist)

            # Normalize selected MOI for later checks
            moi = str(selected_moi or "").upper()
            acmg_disease = selected_disease or str(rule.get("disease", "") or "")

            # Store choice on variants for provenance and downstream logic (used later for AR-blocking etc.)
            for v in varlist:
                try:
                    v._selected_acmg_disease = acmg_disease
                    v._selected_moi = moi
                    v._selected_disease_choice_reason = disease_choice_reason
                except Exception:
                    pass

            # --- Normalise effective MOI for this sample/gene and apply single-het carrier suppression ---
            try:
                eff_moi = (moi or "").upper()
            except Exception:
                eff_moi = str(rule.get("moi", "") or "").upper()

            # Check whether this gene/disease is recessive in effective MOI
            eff_is_recessive = ("AR" in eff_moi) or ("RECESSIVE" in eff_moi)
            # also treat explicit X-linked-recessive as recessive for carrier rules:
            eff_is_recessive = eff_is_recessive or ("XLR" in eff_moi or "X-LINKED_RECESSIVE" in eff_moi or "X_LINKED_RECESSIVE" in eff_moi)

            if eff_is_recessive:
                # collect P/LP and VUS>=3 in this sample/gene
                plp_vars = [vv for vv in varlist if getattr(vv, "automated_class", "") in ("Pathogenic", "Likely pathogenic")]
                plp_hets = [vv for vv in plp_vars if _is_het(getattr(vv, "proband_gt", "./."))]
                plp_homs_or_hemi = [vv for vv in plp_vars if _is_hom(getattr(vv, "proband_gt", "./.")) or (
                    (getattr(vv, "sample_sex", "").lower() == "male")
                    and (vv.chrom.upper().replace("CHR","") in ("X","Y"))
                    and (sum(1 for a in re.split(r'[\/|]', str(getattr(vv, "proband_gt", "") or "")) if str(a) == "1") == 1)
                )]
                vus5 = [vv for vv in varlist if getattr(vv, "automated_class", "") == "VUS" and (getattr(vv, "total_points", 0) >= 3) and _is_het(getattr(vv, "proband_gt", "./."))]

                # If there is exactly one heterozygous P/LP, NO homozygous/hemi P/LP, NO VUS>=3 partners and no other P/LP,
                # treat whole set as carrier -> FILTER ALL (do not escalate to manual nor auto).
                if len(plp_vars) == 1 and len(plp_hets) == 1 and len(plp_homs_or_hemi) == 0 and len(vus5) == 0:
                    note = "Filtered: Recessive heterozygous carrier (single P/LP het with no qualifying partner) - not returned as secondary finding"
                    for vv in varlist:
                        vv._filtered = True
                        vv._in_auto = False
                        vv._in_manual = False
                        if not hasattr(vv, "filtered_reasons") or vv.filtered_reasons is None:
                            vv.filtered_reasons = []
                        if note not in vv.filtered_reasons:
                            vv.filtered_reasons.append(note)
                    # skip further processing for this gene/sample
                    continue

            # --- AR-only handling with homo/hemizygote allowance and pairing/manual rules ---
            # Use the effective MOI determined per-sample (eff_moi) 
            moi_field_raw = (getattr(varlist[0], "_selected_moi", None) or moi or str(rule.get("moi", "") or "")).upper()
            is_ar_flag = ("AR" in moi_field_raw) or ("RECESSIVE" in moi_field_raw) or ("XLR" in moi_field_raw)
            has_dom_flag = any(tok in moi_field_raw for tok in ("AD", "DOMINANT", "AUTOSOMAL_DOMINANT", "XLD", "X_LINKED_DOMINANT"))
            is_ar_only = bool(is_ar_flag and not has_dom_flag)

            # Robust hemizygous detection (treat male X/Y single-alt as hemizygous)
            def _is_hemizygous_variant(v: VariantRecord) -> bool:
                gt = getattr(v, "proband_gt", "") or ""
                # clear homozygous first
                if _is_hom(gt):
                    return True
                # check genotypes like "1/." or "./1" or single-allele encodings
                if "|" in gt or "/" in gt:
                    parts = re.split(r'[\/|]', gt)
                    try:
                        # single alt allele present and other allele missing -> hemizygous-like
                        if len(parts) == 2 and (parts.count("1") == 1) and ("." in parts):
                            return True
                    except Exception:
                        pass
                # consider male on X/Y with a single alt allele as hemizygous
                try:
                    sex = getattr(varlist[0], "sample_sex", "").lower() if varlist else ""
                    chrom_norm = (varlist[0].chrom.upper().replace("CHR","")) if varlist else ""
                    if sex == "male" and chrom_norm in ("X", "Y"):
                        if "1" in gt and not _is_hom(gt):
                            return True
                except Exception:
                    pass
                return False

            # =============================================================
            # AR PAIR LOGIC (UPDATED)
            # =============================================================
            if is_ar_only:
                # Helper to get total points
                def _get_points_ar(v):
                    try:
                        pts = int(getattr(v, "total_ext_points", 0) or 0)
                        if pts == 0:
                            pts = int(getattr(v, "total_points", 0) or 0)
                        return pts
                    except Exception:
                        return 0
                
                # Classify variants in THIS SAMPLE + THIS GENE
                plp_vars = [v for v in varlist if v.automated_class in ("Pathogenic", "Likely pathogenic")]
                vus_med = [v for v in varlist if v.automated_class == "VUS" and _get_points_ar(v) >= 3]
                vus_low = [v for v in varlist if v.automated_class == "VUS" and _get_points_ar(v) < 3]
                
                total_plp = len(plp_vars)
                has_vus_ge3 = len(vus_med) > 0

                # Check for single heterozygous P/LP with no partner
                if total_plp == 1 and len(plp_hets) == 1 and len(plp_homs_or_hemi) == 0:
                    single_het = plp_hets[0]
                    
                    # Check if there's ANY VUS that could be a potential partner
                    potential_partners = [v for v in vus_med if _get_points_ar(v) >= 3]
                    
                    if not potential_partners:
                        # Single heterozygous P/LP with no VUS >=3 - FILTER (carrier)
                        note = f"Filtered: AR gene {gene} heterozygous carrier (single P/LP with no qualifying partner) - not returned as secondary finding"
                        k = get_key(single_het)
                        single_het._filtered = True
                        single_het._in_auto = False
                        single_het._in_manual = False
                        if note not in single_het.filtered_reasons:
                            single_het.filtered_reasons.append(note)
                        if k in manual_ids:
                            manual_ids.remove(k)
                        if k in auto_ids:
                            auto_ids.remove(k)
                        continue  
                                        
                # CASE 1: Exactly 2 P/LP variants
                if total_plp == 2:
                    v1, v2 = plp_vars[0], plp_vars[1]
                    
                    # Check if variants are from the SAME sample (proband)
                    are_same_sample = (v1.sample == v2.sample)
                    
                    # Check if there are any OTHER significant variants (VUS with >=3 points or P/LP)
                    other_significant_variants = [v for v in varlist 
                                                if v not in (v1, v2) 
                                                and (v.extended_class in ("Pathogenic", "Likely pathogenic") 
                                                        or (v.extended_class == "VUS" and _get_points_ar(v) >= 3))]
                    
                    has_extra_significant = len(other_significant_variants) > 0
                    
                    # Same sample with two P/LP variants -> MUST go to MANUAL review
                   
                    if are_same_sample:
                        # Same sample with two P/LP variants -> ALWAYS MANUAL
                        reason = f"Manual review required: AR gene {gene} with two P/LP variants in sample ({v1.sample}) - requires confirmation to verify compound heterozygosity"
                        
                        # Send BOTH variants to MANUAL
                        for pvar in (v1, v2):
                            k = get_key(pvar)
                            if k not in manual_ids:
                                if reason not in pvar.manual_reasons:
                                    pvar.manual_reasons.append(reason)
                                pvar._in_manual = True
                                pvar._in_auto = False
                                pvar._filtered = False
                                manual_ids.add(k)
                                if k in auto_ids:
                                    auto_ids.remove(k)
                        
                        # Filter out any other variants (VUS <3)
                        for pvar in varlist:
                            if pvar not in (v1, v2):
                                k = get_key(pvar)
                                pvar._filtered = True
                                pvar._in_auto = False
                                pvar._in_manual = False
                                note = f"Filtered: VUS with {_get_points_ar(pvar)} points (<3) in AR gene {gene}"
                                if note not in pvar.filtered_reasons:
                                    pvar.filtered_reasons.append(note)
                                if k in manual_ids:
                                    manual_ids.remove(k)
                                if k in auto_ids:
                                    auto_ids.remove(k)
                        
                        continue  # Skip further processing for this gene/sample
                    
                    # Different samples (e.g., proband and parent) -> treat as separate carriers
                    elif v1.sample != v2.sample:
                        logger.warning(f"AR Pair: variants from different samples {v1.sample} vs {v2.sample} - treating as separate carriers")
                        for pvar in (v1, v2):
                            k = get_key(pvar)
                            pvar._filtered = True
                            pvar._in_auto = False
                            pvar._in_manual = False
                            note = f"Filtered: Single P/LP in {gene} for sample {pvar.sample} - carrier"
                            if note not in pvar.filtered_reasons:
                                pvar.filtered_reasons.append(note)
                            if k in manual_ids:
                                manual_ids.remove(k)
                            if k in auto_ids:
                                auto_ids.remove(k)
                        continue
                    
                    # Same sample but has extra significant variants -> already handled above, but if there are extra:
                    elif has_extra_significant:
                        logger.info(f"AR gene {gene}: 2 P/LP + {len(other_significant_variants)} significant extra variant(s) in same sample -> MANUAL review")
                        
                        # Send ALL significant variants to MANUAL
                        all_significant = [v1, v2] + other_significant_variants
                        for pvar in all_significant:
                            k = get_key(pvar)
                            if k not in manual_ids:
                                pvar.manual_reasons.append(
                                    f"Manual review required: AR gene {gene} with 2 P/LP variants and {len(other_significant_variants)} additional significant variant(s) in same sample - complex case"
                                )
                                pvar._in_manual = True
                                pvar._in_auto = False
                                pvar._filtered = False
                                manual_ids.add(k)
                                if k in auto_ids:
                                    auto_ids.remove(k)
                        
                        # Filter out insignificant variants (VUS <3)
                        for pvar in varlist:
                            if pvar not in all_significant:
                                k = get_key(pvar)
                                pvar._filtered = True
                                pvar._in_auto = False
                                pvar._in_manual = False
                                note = f"Filtered: VUS with <3 points in AR gene {gene}"
                                if note not in pvar.filtered_reasons:
                                    pvar.filtered_reasons.append(note)
                                if k in manual_ids:
                                    manual_ids.remove(k)
                                if k in auto_ids:
                                    auto_ids.remove(k)
                        
                        continue
                
                # CASE 2: 1 P/LP + at least one VUS >=3 -> MANUAL (only P/LP and VUS >=3)
                if total_plp == 1 and has_vus_ge3:
                    # Collect ALL variants that should go to manual review
                    variants_to_send_to_manual = []
                    
                    # Add the P/LP variant
                    variants_to_send_to_manual.extend(plp_vars)
                    
                    # Add VUS with >=3 points
                    variants_to_send_to_manual.extend(vus_med)
                    
                    reason = f"Manual review required: AR gene {gene} with 1 P/LP + {len(vus_med)} VUS(>=3) - requires phasing/pairing review"
                    
                    # Send ALL collected variants to manual review
                    for v in variants_to_send_to_manual:
                        k = get_key(v)
                        if k not in manual_ids:
                            if reason not in v.manual_reasons:
                                v.manual_reasons.append(reason)
                            v._in_manual = True
                            v._in_auto = False
                            v._filtered = False
                            manual_ids.add(k)
                            if k in auto_ids:
                                auto_ids.remove(k)
                    
                    # Filter out VUS <3 (these stay filtered)
                    for v in vus_low:
                        k = get_key(v)
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        note = f"Filtered: VUS with {_get_points_ar(v)} points (<3) in AR gene {gene} - not part of manual review pair"
                        if note not in v.filtered_reasons:
                            v.filtered_reasons.append(note)
                        if k in manual_ids:
                            manual_ids.remove(k)
                        if k in auto_ids:
                            auto_ids.remove(k)
                    
                    continue  # Skip further processing for this gene/sample
                
                # CASE 3: 3+ P/LP variants -> MANUAL (all P/LP + VUS >=3 only)
                if total_plp >= 3:
                    all_candidates_for_manual = set()
                    
                    # Add ALL P/LP variants
                    for v in plp_vars:
                        all_candidates_for_manual.add(get_key(v))
                    
                    # Add VUS with >=3 points
                    for v in vus_med:
                        all_candidates_for_manual.add(get_key(v))
                    
                    reason = f"Manual review required: AR gene {gene} with {total_plp} P/LP variants (>2) - multiple candidates"
                    if vus_med:
                        reason += f" + {len(vus_med)} VUS(>=3)"
                    
                    for k in all_candidates_for_manual:
                        v = next((p for p in varlist if get_key(p) == k), None)
                        if v is None:
                            continue
                        if k not in manual_ids:
                            if reason not in v.manual_reasons:
                                v.manual_reasons.append(reason)
                            v._in_manual = True
                            v._in_auto = False
                            v._filtered = False
                            manual_ids.add(k)
                            if k in auto_ids:
                                auto_ids.remove(k)
                    
                    # Filter out VUS with 3-4 points
                    for v in vus_med:
                        k = get_key(v)
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        note = f"Filtered: VUS with 1-2 points in AR gene {gene} with 3+ P/LP - not part of manual review"
                        if note not in v.filtered_reasons:
                            v.filtered_reasons.append(note)
                        if k in manual_ids:
                            manual_ids.remove(k)
                        if k in auto_ids:
                            auto_ids.remove(k)
                    
                    # Filter out VUS <3
                    for v in vus_low:
                        k = get_key(v)
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        note = f"Filtered: VUS with <3 points in AR gene {gene}"
                        if note not in v.filtered_reasons:
                            v.filtered_reasons.append(note)
                        if k in manual_ids:
                            manual_ids.remove(k)
                        if k in auto_ids:
                            auto_ids.remove(k)
                    
                    continue  # Skip further processing for this gene/sample
                
                # CASE 4: ALL OTHER CASES -> FILTER EVERYTHING
                for v in varlist:
                    k = get_key(v)
                    v._filtered = True
                    v._in_auto = False
                    v._in_manual = False
                    pts = _get_points_ar(v)
                    
                    if total_plp == 1 and not has_vus_ge3:
                        note = f"Filtered: AR gene {gene} with single P/LP variant but no VUS(>=3) partner - carrier, not reported"
                    elif total_plp == 0 and has_vus_ge3:
                        note = f"Filtered: AR gene {gene} with VUS(>=3) but no P/LP - not reportable"
                    else:
                        note = f"Filtered: AR gene {gene} variant not reportable (class={v.automated_class}, points={pts})"
                    
                    if note not in v.filtered_reasons:
                        v.filtered_reasons.append(note)
                    
                    if k in manual_ids:
                        manual_ids.remove(k)
                    if k in auto_ids:
                        auto_ids.remove(k)
                
                continue  # Skip remaining logic for this gene/sample

            # HQ HGMD phenotype mismatch check with ClinVar fallback:
            # If a variant has high-quality HGMD evidence (is_hq_pathogenic_db) and an HGMD phenotype
            # that does not match any ACMG-listed disease for the gene, then check ClinVar phenotype.
            # Only if neither HGMD nor ClinVar match we escalate the variant to manual review.
            disease_tokens = [d.strip() for d in re.split(r'[;/\|]+', str(rule.get("disease","") or "")) if d.strip()]
            for v in varlist:
                try:
                    if not getattr(v, "is_hq_pathogenic_db", False):
                        continue
                    hgmd_phen = str(getattr(v, "hgmd_phen", "") or "").strip()
                    if not hgmd_phen:
                        continue

                    # 1) try HGMD -> ACMG disease match
                    matched_hgmd = False
                    for dtk in disease_tokens:
                        if not dtk:
                            continue
                        matched, level, reason = check_disease_match(dtk, hgmd_phen, extra_recs_list=None, gene_rules=None)
                        if matched:
                            matched_hgmd = True
                            break

                    if matched_hgmd:
                        # HGMD phenotype matches ACMG disease — no manual escalation based on HGMD mismatch
                        continue

                    # 2) HGMD did not match — try ClinVar phenotype (if present)
                    clinvar_trait = str(getattr(v, "clinvar_trait", "") or "").strip()
                    matched_clinvar = False
                    if clinvar_trait:
                        for dtk in disease_tokens:
                            if not dtk:
                                continue
                            matched, level, reason = check_disease_match(dtk, clinvar_trait, extra_recs_list=None, gene_rules=None)
                            if matched:
                                matched_clinvar = True
                                break

                    if matched_clinvar:
                        # ClinVar phenotype matches ACMG disease; add an info note but DO NOT force manual review
                        note = f"HGMD phenotype '{hgmd_phen}' mismatches ACMG disease but ClinVar trait '{clinvar_trait}' matches; no manual escalation."
                        if note not in v.manual_reasons:
                            v.manual_reasons.append(note)
                        continue

                    # 3) Neither HGMD nor ClinVar matched -> decide: escalate to manual or treat as carrier
                    # If this variant is a P/LP heterozygote in an effective recessive MOI with no partners -> treat as carrier (filter) instead of manual.
                    is_plp = getattr(v, "automated_class", "") in ("Pathogenic", "Likely pathogenic")
                    is_het = _is_het(getattr(v, "proband_gt", "./."))
                    # use effective MOI for this variant (fall back to gene-level eff_moi)
                    eff_moi_v = (getattr(v, "_selected_moi", None) or moi or str(rule.get("moi", "") or "")).upper()
                    is_recessive_v = ("AR" in eff_moi_v) or ("RECESSIVE" in eff_moi_v) or ("XLR" in eff_moi_v)

                    if is_plp and is_het and is_recessive_v:
                        # check if there are other qualifying partners in varlist (P/LP hom/hemi or VUS>=3 or other P/LP het)
                        other_plp_vars = [x for x in varlist if x is not v and getattr(x, "automated_class", "") in ("Pathogenic", "Likely pathogenic")]
                        other_plp_hom_hemi = [x for x in other_plp_vars if _is_hom(getattr(x, "proband_gt", "./.")) or (
                            (getattr(x, "sample_sex", "").lower() == "male") and (x.chrom.upper().replace("CHR","") in ("X","Y")) and (sum(1 for a in re.split(r'[\/|]', str(getattr(x, "proband_gt", "") or "")) if str(a) == "1") == 1)
                        )]
                        other_vus5 = [x for x in varlist if x is not v and getattr(x, "automated_class", "") == "VUS" and (getattr(x, "total_points", 0) >= 3) and _is_het(getattr(x, "proband_gt", "./."))]

                        if len(other_plp_vars) == 0 and len(other_plp_hom_hemi) == 0 and len(other_vus5) == 0:
                            # treat as carrier -> filter
                            carrier_note = "Filtered: Recessive heterozygous carrier (single P/LP het with no qualifying partner) - disease-mismatch does not escalate to manual"
                            if carrier_note not in v.filtered_reasons:
                                v.filtered_reasons.append(carrier_note)
                            v._filtered = True
                            v._in_auto = False
                            v._in_manual = False
                            # do not add to manual_ids
                        else:
                            # escalate to manual (complex case / potential true reportable mismatch)
                            if note not in v.manual_reasons:
                                v.manual_reasons.append(note)
                            v._in_manual = True
                            v._in_auto = False
                            v._filtered = False
                            manual_ids.add(get_key(v))
                    else:
                        # non-recessive or hom/germline cases: escalate as before
                        if note not in v.manual_reasons:
                            v.manual_reasons.append(note)
                        v._in_manual = True
                        v._in_auto = False
                        v._filtered = False
                        manual_ids.add(get_key(v))
                except Exception:
                    continue

            # inheritance flags
            is_ar = ("AR" in moi) or ("RECESSIVE" in moi)
            is_xlr = ("XLR" in moi) or ("X-LINKED RECESSIVE" in moi)
            is_xl_male = (("X" in moi) or ("X-LINKED" in moi) or is_xlr) and (sample_sex == "Male")


            biologically_valid: List[VariantRecord] = []
            for v in varlist:
                # Preserve any previously set manual/auto flags (e.g. from earlier HGMD/ClinVar checks).
                # Initialize flags only when they are not already present to avoid clearing decisions made above.
                v._in_auto = getattr(v, "_in_auto", False)
                v._in_manual = getattr(v, "_in_manual", False)
                # Default filtered state: preserve existing if set, otherwise True until proven otherwise
                v._filtered = getattr(v, "_filtered", True)
                # Ensure lists exist
                if not hasattr(v, "manual_reasons"):
                    v.manual_reasons = []
                if not hasattr(v, "auto_report_reasons"):
                    v.auto_report_reasons = []
                if not hasattr(v, "filtered_reasons"):
                    v.filtered_reasons = []
                if not hasattr(v, "criteria_explanations"):
                    v.criteria_explanations = {}
                if not hasattr(v, "ps1_matches"):
                    v.ps1_matches = []
                if not hasattr(v, "pm5_matches"):
                    v.pm5_matches = []
                if getattr(v, '_excluded_by_rules', False):
                    continue

                # 1. HQ benign filter (dominant)
                is_hq_benign = False
                try:
                    if getattr(v, "clinvar_stars", "") and str(v.clinvar_stars).isdigit() and int(v.clinvar_stars) >= 2:
                        sl = (v.clinvar_sig or "").lower()
                        if ("benign" in sl or "likely benign" in sl) and ("pathogen" not in sl) and ("conflict" not in sl):
                            is_hq_benign = True
                except Exception:
                    pass
                #try:
                 #   if getattr(v, "internal_db_sig", "") and ("benign" in v.internal_db_sig.lower() or "likely benign" in v.internal_db_sig.lower()):
                 #       is_hq_benign = True
               # except Exception:
                 #   pass
                try:
                    if getattr(v, "hgmd_class", "") in ("DP","FP","DFP"):
                        if int(getattr(v, "hgmd_publication_count", 0) or 0) >= 2:
                            is_hq_benign = True
                except Exception:
                    pass
                if is_hq_benign:
                    reason = "Filtered: HQ Benign (ClinVar≥2★ or HGMD benign or Internal DB benign)"
                    if reason not in v.filtered_reasons:
                        v.filtered_reasons.append(reason)
                    # remain filtered (do not add to biologically_valid)
                    continue

                # 2. Gene mismatch handling (if any manual_reasons indicate mismatch)
                gene_mismatch_detected = any("gene mismatch" in (r or "").lower() for r in v.manual_reasons)
                if gene_mismatch_detected:
                    # clear DB evidence coming from mismatched sources
                    if any("clinvar" in (r or "").lower() for r in v.manual_reasons):
                        v.clinvar_sig = ""; v.clinvar_stars = "0"; v.clinvar_trait = ""
                    if any("hgmd" in (r or "").lower() for r in v.manual_reasons):
                        v.hgmd_class = ""; v.hgmd_phen = ""; v.hgmd_publication_count = 0
                    # Add audit note but do NOT set manual flag
                    if "Gene mismatch detected -> DB evidence ignored for matching sources" not in v.manual_reasons:
                        v.manual_reasons.append("Gene mismatch detected -> DB evidence ignored for matching sources")

                # 3. disease match (collect DB traits and check)
                traits = []
                if getattr(v, "hgmd_phen", "") and v.hgmd_phen.lower() not in ("not provided", ".", "not_provided", "", "null"):
                    traits.append(v.hgmd_phen)
                
                # IMPORTANT: Use raw submissions, not v.clinvar_trait (which is set later)
                # Collect phenotypes from clinvar_submissions
                clinvar_phens_for_match = []
                for sub in getattr(v, "clinvar_submissions", []):
                    phen = str(sub.get("phenotype", "")).strip()
                    if phen and phen.lower() not in ("not provided", ".", "null"):
                        # Remove prefixes like "C1234567:"
                        phen_clean = re.sub(r'^[A-Z0-9]+:', '', phen)
                        clinvar_phens_for_match.append(phen_clean)
                if clinvar_phens_for_match:
                    traits.extend(clinvar_phens_for_match)

                full_trait = " | ".join(_uniq_preserve(traits)) if traits else ""

                # Always call check_disease_match and store results
                is_match, match_level, match_reason = check_disease_match(
                    acmg_disease,
                    full_trait,
                    [spec_rules] if spec_rules else None,
                    rule
                )

                # After constructing full_trait and calling check_disease_match:
                v._disease_match = bool(is_match)
                v._disease_match_level = str(match_level) if match_level else "N/A"
                v._disease_match_reason = str(match_reason) if match_reason else "N/A"

                # mark whether this variant has any database phenotype/evidence to consider
                v._has_db_evidence = (
                    (len(traits) > 0) or
                    bool((getattr(v, "clinvar_sig", "") or "").strip()) or
                    bool((getattr(v, "hgmd_class", "") or "").strip())
                )

                # --- detect conflict "VUS vs Benign" with ClinVar-first strategy ---
                # Prefer precise ClinVar CLNSIG parsing to detect a Benign vs VUS conflict.
                # If ClinVar is inconclusive, fall back to scanning v.db_hits (legacy).
                try:
                    def _norm_token(tok: str) -> str:
                        t = str(tok or "").lower().strip()
                        t = re.sub(r'[_\-\s]+', ' ', t)
                        return t

                    clin_raw = str(getattr(v, "clinvar_sig", "") or "")
                    clin_has_info = bool(clin_raw and clin_raw.strip() and clin_raw.strip() not in (".", "NA", "None", ""))
                    clin_benign = False
                    clin_vus = False

                    if clin_has_info:
                        # split on common separators used in CLNSIG fields
                        parts = re.split(r'[;,\|/]+', clin_raw)
                        for p in parts:
                            p2 = _norm_token(p)
                            if not p2:
                                continue
                            # treat 'likely benign' and 'benign' as benign
                            if "benign" in p2 and "pathogen" not in p2:
                                clin_benign = True
                            # VUS / Uncertain significance
                            if "uncertain" in p2 or "vus" in p2 or "uncertain significance" in p2:
                                clin_vus = True
                        # If ClinVar explicitly contains both benign and VUS -> conflict detected
                        if clin_benign and clin_vus:
                            if not getattr(v, "criteria_points", None):
                                v.criteria_points = {}
                            if "DB_CONFLICT_VUS_BENIGN" not in v.criteria_points:
                                v.criteria_points["DB_CONFLICT_VUS_BENIGN"] = -1
                            else:
                                v.criteria_points["DB_CONFLICT_VUS_BENIGN"] = -1
                            # recompute total_points conservatively (will be recomputed later as well)
                            try:
                                v.total_points = sum(int(x) for x in v.criteria_points.values())
                            except Exception:
                                v.total_points = getattr(v, "total_points", 0) - 1
                            conflict_msg = "Database conflict (ClinVar): Benign vs VUS detected -> conservative -1 point applied (info in comments)."
                            if conflict_msg not in v.manual_reasons:
                                v.manual_reasons.append(conflict_msg)
                            v._db_vus_benign_conflict = True
                        # skip fallback if ClinVar provided decisive info (either conflict or no conflict)
                        # If ClinVar had tokens but did not indicate both benign and VUS, do nothing here.
                    else:
                        # No usable ClinVar CLNSIG -> fallback to scanning v.db_hits (legacy behavior)
                        db_hits = getattr(v, "db_hits", []) or []
                        has_benign = False
                        has_vus = False
                        for h in db_hits:
                            sig = ""
                            try:
                                if isinstance(h, (list, tuple)) and len(h) > 0:
                                    sig = str(h[0] or "")
                                elif isinstance(h, dict):
                                    sig = str(h.get("clnsig") or h.get("significance") or "")
                                else:
                                    sig = str(h)
                            except Exception:
                                sig = str(h)
                            s = str(sig).lower()
                            if "benign" in s and "pathogen" not in s:
                                has_benign = True
                            if "vus" in s or "uncertain" in s or "uncertain_significance" in s:
                                has_vus = True
                        if has_benign and has_vus:
                            if not getattr(v, "criteria_points", None):
                                v.criteria_points = {}
                            if "DB_CONFLICT_VUS_BENIGN" not in v.criteria_points:
                                v.criteria_points["DB_CONFLICT_VUS_BENIGN"] = -1
                            else:
                                v.criteria_points["DB_CONFLICT_VUS_BENIGN"] = -1
                            try:
                                v.total_points = sum(int(x) for x in v.criteria_points.values())
                            except Exception:
                                v.total_points = getattr(v, "total_points", 0) - 1
                            conflict_msg = "Database conflict: Benign vs VUS detected (fallback) -> conservative -1 point applied (info in comments)."
                            if conflict_msg not in v.manual_reasons:
                                v.manual_reasons.append(conflict_msg)
                            v._db_vus_benign_conflict = True
                except Exception:
                    # conservative: if detection fails, do not apply penalty or force manual
                    pass

                # ----- Gene Rules (from ACMG_SF table 'Rules' column) -----
                # CRITICAL: rules_text must come ONLY from this gene's rule, not from any other gene
                try:
                    # Get rule for the current gene (normalized)
                    current_gene_norm = normalize_gene_name(gene)
                    rule_for_gene = gene_rules.get(current_gene_norm, {}) if gene_rules else {}
                    rules_text = str(rule_for_gene.get("rules_text", "") or "").strip()
                    
                    # ADD DIAGNOSTIC
                    if current_gene_norm == "MSH6":
                        logger.info(f"MSH6: rules_text from gene_rules = '{rules_text}'")
                        logger.info(f"MSH6: rule_for_gene keys = {list(rule_for_gene.keys())}")
                except Exception as e:
                    logger.debug(f"Failed to get rules_text for {gene}: {e}")
                    rules_text = ""

                if rules_text:
                    decision, matched_lines = evaluate_gene_rules_for_variant(rules_text, v)
                    # record applied rules on variant for output (preserve earlier entries)
                    if not hasattr(v, "_rules_applied"):
                        v._rules_applied = []
                    # If some lines matched, record them
                    if matched_lines:
                        v._rules_applied.extend(matched_lines)
                    else:
                        # No rule line matched -> show the raw rules text
                        raw_one_line = " | ".join([ln.strip() for ln in _split_rules_text(rules_text)])
                        v._rules_applied.append(f"RULES_PRESENT_NO_MATCH: {raw_one_line}")

                    if decision == "exclude":
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        reason = f"Filtered by gene 'Rules' (exclude): {'; '.join(matched_lines)}"
                        if reason not in v.filtered_reasons:
                            v.filtered_reasons.append(reason)
                        v._excluded_by_rules = True
                        continue
                                    
                    elif decision == "include":
                        # Check if variant is already filtered by classification
                        # Rules should ONLY apply to variants that are already reportable
                        auto_class = getattr(v, "automated_class", "VUS")
                        total_points = getattr(v, "total_ext_points", getattr(v, "total_points", 0))
                        
                        # Benign/Likely benign variants should NEVER go to manual review
                        # VUS with <3 points should also be filtered
                        is_benign = auto_class in ("Benign", "Likely benign") or total_points <= -1
                        is_low_vus = auto_class == "VUS" and total_points < 3
                        
                        if is_benign or is_low_vus:
                            # Filter out - rules do not apply to benign/low-quality variants
                            filter_reason = f"Filtered: Rules 'include' matched but variant is {auto_class} (points={total_points}) - not reportable"
                            if filter_reason not in v.filtered_reasons:
                                v.filtered_reasons.append(filter_reason)
                            continue  # Skip this variant entirely
                        
                        # Only reach here for Pathogenic/Likely pathogenic or VUS>=3
                        note = f"Included by gene 'Rules': {'; '.join(matched_lines)}"
                        if not hasattr(v, "_comments_parts"):
                            v._comments_parts = []
                        if note not in v._comments_parts:
                            v._comments_parts.append(note)
                        v._rules_included = True
                    # decision == "no_match" - do nothing, variant continues normal processing
                
                # 4. gene-specific positional filtering (if recommendations include position/exon restrictions)
                pos_rules = spec_rules.get("position_rules") if spec_rules else None
                if pos_rules:
                    try:
                        ok_pos = variant_matches_position_rules(v, pos_rules)
                    except Exception:
                        ok_pos = False
                    if not ok_pos:
                        reason = f"Filtered: outside recommended transcript/exon/protein region per recommendations: {pos_rules}"
                        if reason not in v.filtered_reasons:
                            v.filtered_reasons.append(reason)
                        continue

                # 5. gene-specific variant type restrictions
                if spec_rules.get("missense_only") and "missense" not in (v.consequence or "").lower():
                    reason = f"Filtered: gene rule 'missense only' but consequence is {v.consequence}"
                    if reason not in v.filtered_reasons:
                        v.filtered_reasons.append(reason)
                    continue
                if spec_rules.get("hemizygous_only_male") and sample_sex == "Female":
                    reason = f"Filtered: gene rule 'hemizygous only male' (female excluded)"
                    if reason not in v.filtered_reasons:
                        v.filtered_reasons.append(reason)
                    continue

                # If passed all above, mark as not filtered and add to biologically_valid
                v._filtered = False
                biologically_valid.append(v)

            # After building biologically_valid for the gene/sample, ensure filtered variants have at least one reason
            for vv in varlist:
                if getattr(vv, "_filtered", True):
                    if not getattr(vv, "filtered_reasons"):
                        vv.filtered_reasons.append("Filtered: did not pass gene-specific or QC/priority filters")

            # define drivers and passengers
            drivers = [v for v in biologically_valid if v._qc_passed and (getattr(v, "is_hq_pathogenic_db", False) or v.automated_class in ("Pathogenic", "Likely pathogenic"))]
            passengers = [v for v in biologically_valid if v.automated_class in ("Pathogenic", "Likely pathogenic", "VUS")]

            # detect delins / flags if helper exists
            try:
                # Build map cluster_id -> members (VariantRecord objects)
                cluster_map = defaultdict(list)
                for vv in biologically_valid:
                    cid = getattr(vv, "_delins_cluster_id", None)
                    if cid:
                        cluster_map[cid].append(vv)

                for cid, members in cluster_map.items():
                    # only consider clusters with >=2 variants
                    if not members or len(members) < 2:
                        continue
                    
                    passing_members = []
                    for m in members:
                        meets = False
                        try:
                            if getattr(m, "automated_class", "") in ("Pathogenic", "Likely pathogenic") and getattr(m, "total_points", 0) > 3:
                                meets = True
                            elif getattr(m, "total_points", 0) > 3:
                                meets = True
                        except Exception:
                            meets = False
                        
                        if meets:
                            passing_members.append(m)
                    
                    if len(passing_members) >= 2:
                        member_keys = [f"{x.chrom}:{x.pos}:{x.ref}:{x.alt}" for x in passing_members]
                        for m in passing_members:
                            k = get_key(m)
                            if k not in manual_ids:
                                m.manual_reasons.append(f"Manual: Delins/overlap cluster (id={cid}) meets criteria -> cluster manual review required. Members: {', '.join(member_keys)}")
                                m._in_manual = True
                                m._filtered = False
                                m._in_auto = False
                                manual_ids.add(k)
                                if k in auto_ids:
                                    auto_ids.remove(k)
                    
            except Exception:
                # conservative: if evaluation fails, do not change manual/auto status
                pass

            # =============================================================
            # X-chromosome variants -> manual review for ALL variants with >=3 points
            # =============================================================
            for v in passengers:
                if not v._qc_passed:
                    continue
                
                # Check if variant is on X chromosome
                chrom_norm = str(v.chrom or "").upper().replace("CHR", "")
                if chrom_norm != "X":
                    continue
                
                k = get_key(v)
                if k in auto_ids or k in manual_ids:
                    continue
                
                # Get points and class
                pts = int(getattr(v, "total_ext_points", getattr(v, "total_points", 0)) or 0)
                is_plp = v.automated_class in ("Pathogenic", "Likely pathogenic")
                is_vus = v.automated_class == "VUS"
                
                sample_sex = str(getattr(v, "sample_sex", "")).lower()
                sex_status = f"sex={sample_sex}" if sample_sex else "sex=Unknown"
                
                # CRITICAL: ALL X-chromosome variants with >=3 points go to MANUAL
                # Not just P/LP and not just VUS!
                if pts >= 3:
                    v.manual_reasons.append(
                        f"Manual review required: X-chromosome variant (points={pts} >= 3) - {sex_status}. "
                        f"Verify sex and carrier status (hemizygous/heterozygous)."
                    )
                    v._in_manual = True
                    v._in_auto = False
                    v._filtered = False
                    manual_ids.add(k)
                    if k in auto_ids:
                        auto_ids.remove(k)
                
                # VUS with <3 points on X -> filtered
                elif is_vus and pts < 3:
                    v._filtered = True
                    v._in_auto = False
                    v._in_manual = False
                    note = f"Filtered: X-chromosome VUS with {pts} points (<3) - below review threshold"
                    if note not in v.filtered_reasons:
                        v.filtered_reasons.append(note)

            # =============================================================
            # NON-CANONICAL SPLICE VARIANTS with high SpliceAI -> MANUAL review
            # =============================================================
            for v in passengers:
                if not v._qc_passed:
                    continue
                
                k = get_key(v)
                if k in auto_ids or k in manual_ids:
                    continue
                
                # Check if this is a non-canonical splice variant
                is_noncanonical = _is_noncanonical_splice_region_variant(v)
                
                if not is_noncanonical:
                    continue
                
                # Get SpliceAI score
                spliceai_score = getattr(v, "spliceai", None)
                if spliceai_score is None or spliceai_score < 0.7:
                    continue
                
                # Get total points
                pts = int(getattr(v, "total_ext_points", getattr(v, "total_points", 0)) or 0)
                
                # If points >= 3 and non-canonical splice with high SpliceAI -> MANUAL
                if pts >= 3:
                    spliceai_score_str = f"{spliceai_score:.3f}" if spliceai_score is not None else "NA"
                    v.manual_reasons.append(
                        f"Manual review required: Non-canonical splice variant with high SpliceAI={spliceai_score_str} (≥0.7) and {pts} points - "
                        f"requires experimental validation (minigene/RNA-seq)"
                    )
                    v._in_manual = True
                    v._in_auto = False
                    v._filtered = False
                    manual_ids.add(k)
                    if k in auto_ids:
                        auto_ids.remove(k)
                else:
                    # Low points - filter out
                    v._filtered = True
                    v._in_auto = False
                    v._in_manual = False
                    spliceai_score_str = f"{spliceai_score:.3f}" if spliceai_score is not None else "NA"
                    note = f"Filtered: Non-canonical splice variant with SpliceAI={spliceai_score_str} but only {pts} points (<3)"
                    if note not in v.filtered_reasons:
                        v.filtered_reasons.append(note)

            # =============================================================
            # NEW: DSCORE AND IMPACT-BASED MANUAL REVIEW FOR VUS (1-5 POINTS)
            # Rules:
            #   1. Impact = HIGH + VUS (1-5 points) -> ALWAYS manual (Dscore irrelevant)
            #   2. Dscore >= 0.7 + Impact = MODERATE + VUS (1-5 points) -> manual  
            #   3. Dscore == 1.0 + VUS (1-5 points) -> manual (any Impact)
            # =============================================================
            
            # Get thresholds from global config
            DSCORE_MANUAL_THRESHOLD = THRESH.get("DSCORE_MANUAL_THRESHOLD", 0.7)
            DSCORE_PERFECT = THRESH.get("DSCORE_PERFECT", 1.0)
            VUS_POINTS_MIN = THRESH.get("VUS_POINTS_MIN_FOR_MANUAL", 1)
            VUS_POINTS_MAX = THRESH.get("VUS_POINTS_MAX_FOR_MANUAL", 5)
            
            for v in passengers:
                if not v._qc_passed:
                    continue
                
                k = get_key(v)
                
                # Skip if already assigned to auto or manual
                if k in auto_ids or k in manual_ids:
                    continue
                
                # Check if this is a VUS variant with points in range 1-5
                auto_class = getattr(v, "automated_class", "VUS")
                total_points = getattr(v, "total_ext_points", getattr(v, "total_points", 0))
                
                # Only consider VUS variants with points between 1 and 5
                is_vus_range = (auto_class == "VUS" and 
                            VUS_POINTS_MIN <= total_points <= VUS_POINTS_MAX)
                
                if not is_vus_range:
                    continue
                
                # Get Impact and Dscore
                impact = getattr(v, "impact", "").upper()
                dscore = getattr(v, "dscore", None)
                
                trigger_manual = False
                review_reason = None
                
                # Rule 1: Impact = HIGH -> always manual (Dscore irrelevant)
                if impact == "HIGH":
                    trigger_manual = True
                    dscore_str = f"{dscore:.3f}" if dscore is not None else "NA"
                    review_reason = (f"Manual review required: HIGH impact variant with VUS ({total_points} points) - "
                                    f"requires curator assessment (Dscore={dscore_str})")
                    logger.info(f"HIGH impact VUS {total_points} points -> manual: {v.chrom}:{v.pos} {v.hgvsp}")
                
                # Rule 2: Dscore >= 0.7 with Impact = MODERATE -> manual
                elif dscore is not None and dscore >= DSCORE_MANUAL_THRESHOLD and impact == "MODERATE":
                    trigger_manual = True
                    dscore_str = f"{dscore:.3f}" if dscore is not None else "NA"
                    review_reason = (f"Manual review required: High Dscore={dscore_str} with MODERATE impact and VUS ({total_points} points) - "
                                f"requires splicing validation")
                
                # Rule 3: Perfect Dscore (1.0) -> manual regardless of Impact
                elif dscore is not None and abs(dscore - DSCORE_PERFECT) < 0.001:
                    trigger_manual = True
                    dscore_str = f"{dscore:.3f}" if dscore is not None else "NA"
                    review_reason = (f"Manual review required: Perfect Dscore={dscore_str} (==1.0) with VUS ({total_points} points) - "
                                f"requires curator assessment (Impact={impact})")
                
                # Apply manual review if triggered
                if trigger_manual:
                    if k not in manual_ids:
                        v.manual_reasons.append(review_reason)
                        v._in_manual = True
                        v._in_auto = False
                        v._filtered = False
                        manual_ids.add(k)
                        if k in auto_ids:
                            auto_ids.remove(k)


            # --- AUTOMATIC LOGIC ---
            for v in drivers:
                if v.automated_class in ("Benign", "Likely benign"):
                    continue

                k = get_key(v)
    
                # Skip auto-report for VUS (1-5 points) that triggered manual review by Dscore/Impact rules
                if k in manual_ids:
                    continue

                # CRITICAL: X-chromosome variants NEVER go to automatic report
                chrom_norm = str(v.chrom or "").upper().replace("CHR", "")
                if chrom_norm == "X":
                    # X-chromosome variants are handled in X-chromosome logic above
                    # Do not process here
                    continue

                # =========================================================
                # NEW: Non-canonical splice variants with high SpliceAI NEVER go to automatic report
                # =========================================================
                is_noncanonical = _is_noncanonical_splice_region_variant(v)
                if is_noncanonical:
                    spliceai_score = getattr(v, "spliceai", None)
                    if spliceai_score is not None and spliceai_score >= 0.7:
                        # Already handled in manual logic above, but skip here just in case
                        continue

                # Prefer immediate auto-report for homozygous / hemizygous P/LP
                is_plp = v.automated_class in ("Pathogenic", "Likely pathogenic")
                
                # Check if this gene has special rules that prevent auto-reporting of homozygotes
                rule = gene_rules.get(v.gene, {}) if gene_rules else {}
                special_flags = rule.get("special_flags", {})
                
                # For AD genes, homozygous variants should not auto-report without additional evidence
                is_ad_gene = any(tok in moi for tok in ("AD", "DOMINANT", "AUTOSOMAL_DOMINANT"))
                
                if is_plp and (_is_hom(v.proband_gt) or 
                    (getattr(v, "sample_sex", "").lower() == "male" and v.chrom.upper().replace("CHR","") in ("X","Y") and not _is_het(v.proband_gt))):
                    
                    # Skip auto-report for AD genes with homozygous variants (except X-linked)
                    if is_ad_gene and not ("X" in moi.upper()):
                        # Send to manual review instead
                        v.manual_reasons.append(
                            f"Manual review required: Homozygous {v.automated_class} variant in AD gene {gene} - requires curator assessment"
                        )
                        v._in_manual = True
                        v._in_auto = False
                        v._filtered = False
                        manual_ids.add(get_key(v))
                        continue
                    
                    # For AR genes or X-linked, proceed with auto-report
                    v._in_auto = True
                    v._filtered = False
                    v.auto_report_reasons.append("Homozygous/hemizygous P/LP -> automatic report")
                    auto_ids.add(get_key(v))
                    continue

                if get_key(v) in manual_ids:
                    continue

                # ensure lists exist
                if not hasattr(v, "auto_report_reasons"):
                    v.auto_report_reasons = []
                if not hasattr(v, "criteria_explanations"):
                    v.criteria_explanations = {}


                # --- AR blocking is now handled by AR PAIR LOGIC above ---
                # This block is disabled to avoid conflict with new AR logic
                pass

                # Case B: Algorithmic scoring (Tavtigian et al., (2018))
                # Use extended_class and total_ext_points as the single basis for automatic reporting decisions.
                try:
                    ext_cls = getattr(v, "extended_class", "") or getattr(v, "automated_class", "")
                    ext_pts = int(getattr(v, "total_ext_points", 0) or 0)
                except Exception:
                    ext_cls = getattr(v, "automated_class", "")
                    ext_pts = int(getattr(v, "total_points", 0) or 0)

                # Check if this variant was downgraded but still reportable
                is_downgraded_reportable = getattr(v, "_reportable_after_downgrade", False)
                
                if ((ext_cls == "Pathogenic" and ext_pts >= 10) or
                    (ext_cls == "Likely pathogenic" and ext_pts >= 6) or
                    (is_downgraded_reportable and ext_cls == "Likely pathogenic" and ext_pts >= 6)):
                        v._in_auto = True
                        v._filtered = False
                        if is_downgraded_reportable:
                            v.auto_report_reasons.append("Extended score (downgraded from Pathogenic)")
                        else:
                            v.auto_report_reasons.append("Extended score")
                        auto_ids.add(get_key(v))

            # --- MANUAL LOGIC ---
            for v in drivers:
                if get_key(v) in manual_ids:
                    continue

            # --- XL male heterozygous handling ---
            if is_xl_male:
                for v in drivers:
                    if _is_het(v.proband_gt):
                        af_val = _get_sample_af(v)
                        # treat hemizygous if AF >= 0.9 -> do NOT send to manual on this ground alone
                        if af_val <= 0.9:
                            if get_key(v) not in manual_ids:
                                af_val_str = f"{af_val:.3f}" if af_val is not None else "NA"
                                v.manual_reasons.append(f"XL male heterozygous with AF={af_val_str} -> Manual for hemizygosity check")
                                v._in_manual = True; v._filtered = False; v._in_auto = False
                                manual_ids.add(get_key(v))
                                if get_key(v) in auto_ids: auto_ids.remove(get_key(v))
                        else:
                            # treat as hemizygous; leave to auto/manual decisions based on other evidence
                            pass

            # --- Complex recessive: >=2 unique P/LP heterozygotes => manual review for all such P/LP hets ---
            ar_path_lp = []
            seen_keys = set()

            for x in drivers:
                # defensive: skip None entries
                if x is None:
                    continue

                # determine effective MOI for this variant (prefer per-variant selected MOI, then sample-level moi, then rule)
                try:
                    eff_moi_x = (getattr(x, "_selected_moi", None) or moi or str(rule.get("moi", "") or "")).upper()
                except Exception as e:
                    logger.debug("Failed to determine effective MOI for variant %s: %s", getattr(x, "hgvsp", "?"), e)
                    eff_moi_x = str(rule.get("moi", "") or "").upper()

                # defensive: compute unique key; if get_key may raise, catch and continue
                try:
                    k = get_key(x)
                except Exception as e:
                    logger.debug("Failed to compute key for variant object: %s", e)
                    continue

                # already seen this physical variant (dedup)
                if k in seen_keys:
                    continue

                try:
                    is_recessive_moi = any(tok in eff_moi_x for tok in ("AR", "RECESSIVE", "XLR", "X-LINKED_RECESSIVE"))
                    is_het_proband = _is_het(getattr(x, "proband_gt", "./."))
                    is_plp = x.automated_class in ("Pathogenic", "Likely pathogenic")
                except Exception as e:
                    logger.debug("Attribute access failed for variant %s: %s", k, e)
                    continue

                if is_recessive_moi and is_het_proband and is_plp:
                    ar_path_lp.append(x)
                    seen_keys.add(k)

            # Only escalate if there are at least two distinct P/LP heterozygotes
            if len(ar_path_lp) >= 2:
                for topv in ar_path_lp:
                    tk = get_key(topv)
                    if tk not in manual_ids:
                        topv.manual_reasons.append("Manual: Recessive gene with multiple heterozygous P/LP variants - phasing and pairing review required")
                        topv._in_manual = True; topv._filtered = False; topv._in_auto = False
                        manual_ids.add(tk)
                        if tk in auto_ids:
                            auto_ids.remove(tk)
            
            # =============================================================
            # AD VUS LOGIC 
            # Rule: VUS with >=3 points in AD genes -> MANUAL review
            # =============================================================
            
            # Get MOI for this gene 
            gene_moi = (getattr(varlist[0], "_selected_moi", None) if varlist else None) or moi or str(rule.get("moi", "") or "").upper()
            
            # Check if this gene has AD inheritance pattern
            is_ad_gene = any(tok in gene_moi for tok in ("AD", "DOMINANT", "AUTOSOMAL_DOMINANT", "XLD", "X_LINKED_DOMINANT"))
            
            if is_ad_gene:
                for v in passengers:
                    if not v._qc_passed:
                        continue
                    
                    # Skip if already assigned to auto or manual
                    k = get_key(v)
                    if k in auto_ids or k in manual_ids:
                        continue
                    
                    # Only consider VUS variants
                    if v.automated_class != "VUS":
                        continue
                    
                    # Get points for this variant
                    try:
                        pts = int(getattr(v, "total_ext_points", getattr(v, "total_points", 0)) or 0)
                    except Exception:
                        pts = 0
                    
                    # VUS with >=3 points -> manual review
                    if pts >= 3:
                        if k not in manual_ids:
                            v.manual_reasons.append(
                                f"Manual review required: AD inheritance with VUS (points={pts} >= 3)"
                            )
                            v._in_manual = True
                            v._in_auto = False
                            v._filtered = False
                            manual_ids.add(k)
                            if k in auto_ids:
                                auto_ids.remove(k)
                    else:
                        # VUS with <3 points in AD gene - filter out
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        note = f"Filtered: VUS with {pts} points (<3) in AD gene - below review threshold"
                        if note not in v.filtered_reasons:
                            v.filtered_reasons.append(note)

    # --- Build output rows (augment with explanations) ---
    # Define columns including per-criterion explanation columns
    COLUMNS = [
        "Query_Variant","sample","Family_ID","chrom","pos","ref","alt","gene","Inheritance",
        "HGVSc","HGVSp","AA_ref","AA_alt","AA_pos","tid","Exon","Consequence","Impact","Zygosity","DP","Sample_AF",
        "acmg_disease_association","Special_Guidance", "Rules_applied",
        "clinvar_phenotype_match","clinvar_phenotype_info",
        "hgmd_phenotype_match","hgmd_phenotype_info",
        "clinvar_trait","clinvar_sig","clinvar_sig_consensus", "clinvar_stars", "clinvar_star_consensus",
        "hgmd_class","hgmd_phenotype","hgmd_publication_count","hgmd_rankscore", "hgmd_dmsupported",
        "internal_db_sig","database_summary",
        "gnomAD_AF","gnomAD_version",
        "DSCORE","REVEL","SpliceAI","SQUIRLS","AlphaMissense",
        "ps1_matches","pm5_matches","phase","is_delins_part",
        "criteria_points_ext","total_ext_points","extended_class",
        "PP5ext_BP5ext_expl","PP5_clin_block", "PP5_hgmd_block", "PP5_raw_total", "PP5_applied", 
        "PVS1_expl","PS1_expl","PM1_expl","PM2_expl","PM3_expl","PM4_expl","PM5_expl",
        "PP1_expl","PP3_expl",
        "BA1_expl","BS1_expl","BS2_expl","BP4_expl","BP6_expl","BP7_expl",
        # reporting/audit
        "in_auto_report","in_manual_review","auto_report_reasons","manual_reasons","comments","filtered_reasons","filtered_out"
    ]

    processed = set()
    rows_auto = []
    rows_manual = []
    rows_all = []

    for v in all_variants:
        k = get_key(v)
        
        # Skip variants that were excluded by gene rules
        if getattr(v, '_excluded_by_rules', False):
            v._filtered = True
            v._in_auto = False
            v._in_manual = False
            manual_ids.discard(k)
            auto_ids.discard(k)
            # Do NOT continue - we still want this variant in all_candidates
            # but with filtered_out = "Yes"
        
        # CRITICAL: Filter out benign/LB variants
        auto_class = getattr(v, "automated_class", "VUS")
        total_points = getattr(v, "total_ext_points", getattr(v, "total_points", 0))
        
        # Benign/Likely benign variants must NEVER go to manual review
        if auto_class in ("Benign", "Likely benign") or total_points <= -1:
            v._filtered = True
            v._in_auto = False
            v._in_manual = False
            manual_ids.discard(k)
            auto_ids.discard(k)
            if not hasattr(v, "filtered_reasons"):
                v.filtered_reasons = []
            reason = f"Filtered: {auto_class} variant (points={total_points}) - not reportable"
            if reason not in v.filtered_reasons:
                v.filtered_reasons.append(reason)
            # Do NOT continue - we still want this variant in all_candidates
        
        # VUS with <3 points should not go to manual review
        elif auto_class == "VUS" and total_points < 3:
            v._filtered = True
            v._in_auto = False
            v._in_manual = False
            manual_ids.discard(k)
            auto_ids.discard(k)
            if not hasattr(v, "filtered_reasons"):
                v.filtered_reasons = []
            reason = f"Filtered: VUS with {total_points} points (<3) - below reporting threshold"
            if reason not in v.filtered_reasons:
                v.filtered_reasons.append(reason)
            # Do NOT continue - we still want this variant in all_candidates

        # Continue with adding to rows (same as before)
        if k in processed:
            continue
        processed.add(k)

        rule = gene_rules.get(v.gene, {}) if gene_rules else {}
        row = variant_row_with_disease(v, rule, family_id, recs_map, gene_rules)

        # Compute extended totals & class (prefer precomputed snapshot on variant).
        try:
            # Prefer an already computed extended total if available
            if getattr(v, "total_ext_points", None) is not None:
                total_ext = int(v.total_ext_points)
            else:
                # prefer criteria_points_ext (dict) if present
                pts_ext = getattr(v, "criteria_points_ext", None)
                if pts_ext is None:
                    # fallback: try to construct from criteria_points (base pts)
                    base_pts = getattr(v, "criteria_points", {}) or {}
                    # If base_pts is a JSON string, try to parse it
                    if isinstance(base_pts, str):
                        try:
                            base_pts = json.loads(base_pts)
                        except Exception:
                            base_pts = {}
                    # remove legacy PP5/BP5 keys if present so extended scoring remains canonical
                    pts_ext = dict(base_pts)
                    for kk in list(pts_ext.keys()):
                        if kk.startswith("PP5") or kk.startswith("BP5"):
                            pts_ext.pop(kk, None)
                    # attach to variant for later reuse
                    v.criteria_points_ext = pts_ext
                else:
                    # pts_ext present: ensure it's a dict (if string parse it)
                    if isinstance(pts_ext, str):
                        try:
                            pts_ext = json.loads(pts_ext)
                        except Exception:
                            pts_ext = {}
                        v.criteria_points_ext = pts_ext

                # sum up ext points
                try:
                    total_ext = sum(int(x) for x in (pts_ext or {}).values()) if pts_ext else 0
                except Exception:
                    # defensive fallback
                    total_ext = 0

            # compute extended_class using the Tavtigian points model function
            try:
                extended_class = compute_class_from_points(getattr(v, "criteria_points_ext", {}) or {})[0]
            except Exception:
                # fallback to existing automated_class if compute fails
                extended_class = getattr(v, "automated_class", "VUS") or "VUS"

        except Exception:
            total_ext = 0
            extended_class = getattr(v, "automated_class", "VUS") or "VUS"

        v.total_ext_points = int(getattr(v, "total_ext_points", 0) or 0)
        v.extended_class = getattr(v, "extended_class", "") or "VUS"
        v.automated_class = v.extended_class
        v.total_points = int(v.total_ext_points or 0)

        row["total_ext_points"] = v.total_ext_points
        row["extended_class"] = v.extended_class
        row["automated_class"] = v.automated_class

        for vv in candidates:
            try:
                # Ensure criteria_points_ext exists (preferred). If missing, construct from criteria_points
                pts_ext = getattr(vv, "criteria_points_ext", None)
                if pts_ext is None:
                    base_pts = getattr(vv, "criteria_points", {}) or {}
                    if isinstance(base_pts, str):
                        try:
                            base_pts = json.loads(base_pts)
                        except Exception:
                            base_pts = {}
                    pts_ext = dict(base_pts)
                    # Remove legacy PP5/BP5 keys if present so extended scoring remains canonical
                    for kk in list(pts_ext.keys()):
                        if kk.startswith("PP5") or kk.startswith("BP5"):
                            pts_ext.pop(kk, None)
                    vv.criteria_points_ext = pts_ext
                # Ensure numeric total_ext_points is present and consistent
                try:
                    vv.total_ext_points = int(getattr(vv, "total_ext_points", 0) or 0)
                except Exception:
                    try:
                        vv.total_ext_points = sum(int(x) for x in (vv.criteria_points_ext or {}).values()) if vv.criteria_points_ext else 0
                    except Exception:
                        vv.total_ext_points = int(getattr(vv, "total_points", 0) or 0)
                # Ensure extended_class is present; compute if missing
                try:
                    if not getattr(vv, "extended_class", None):
                        vv.extended_class = compute_class_from_points(vv.criteria_points_ext or {})[0]
                except Exception:
                    vv.extended_class = getattr(vv, "automated_class", "VUS") or "VUS"
                # Make automated_class mirror extended_class (single canonical classification)
                vv.automated_class = vv.extended_class
                # Keep legacy total_points synchronized (so older code still sees canonical value)
                vv.total_points = int(vv.total_ext_points or 0)
            except Exception:
                # Safe fallback to avoid crash; preserve existing values where possible
                vv.automated_class = getattr(vv, "automated_class", "VUS") or "VUS"
                vv.extended_class = vv.automated_class
                try:
                    vv.total_ext_points = int(getattr(vv, "total_ext_points", 0) or 0)
                except Exception:
                    vv.total_ext_points = int(getattr(vv, "total_points", 0) or 0)
                vv.total_points = int(vv.total_ext_points or 0)

        # --- extended totals & class (PP5ext/BP5ext) ---
        # Populate canonical extended scoring fields for CSV output.
        # Use criteria_points_ext / total_ext_points / extended_class as the single source of truth.
        try:
            pts_ext = getattr(v, "criteria_points_ext", {}) or {}
            if isinstance(pts_ext, str):
                try:
                    pts_ext = json.loads(pts_ext)
                except Exception:
                    pts_ext = {}
        except Exception:
            pts_ext = {}

        # Compute or prefer existing total_ext_points
        try:
            total_ext_val = int(getattr(v, "total_ext_points", 0) or 0)
        except Exception:
            try:
                total_ext_val = sum(int(x) for x in (pts_ext or {}).values()) if pts_ext else 0
            except Exception:
                total_ext_val = 0

        # Compute or prefer existing extended_class
        try:
            ext_cls = getattr(v, "extended_class", None)
            if not ext_cls:
                ext_cls = compute_class_from_points(pts_ext)[0] if pts_ext else "VUS"
        except Exception:
            ext_cls = getattr(v, "automated_class", "VUS") or "VUS"

        # Synchronize variant fields
        v.criteria_points_ext = pts_ext
        v.total_ext_points = int(total_ext_val)
        v.extended_class = ext_cls
        v.automated_class = v.extended_class  # keep mirrored for compatibility
        v.total_points = int(v.total_ext_points or 0)

        # Serialize extended points as JSON for output
        try:
            criteria_points_ext_json = json.dumps(pts_ext or {}, ensure_ascii=False)
        except Exception:
            criteria_points_ext_json = "{}"

        # Phenotype match flags and PP5ext explanation
        clin_match_flag = getattr(v, "_clinvar_phenotype_match", False)
        hgmd_match_flag = getattr(v, "_hgmd_phenotype_match", False)
        pp5ext_expl = getattr(v, "_pp5ext_explanation", "")

        # Fill output row with canonical extended fields
        row["criteria_assigned_ext"] = ";".join(getattr(v, "criteria_assigned_ext", []) or [])
        row["criteria_points_ext"] = criteria_points_ext_json
        row["total_ext_points"] = v.total_ext_points
        row["extended_class"] = v.extended_class
        # Keep automated_class column for backward compatibility (it will equal extended_class)
        row["automated_class"] = v.automated_class

        row["PP5ext_BP5ext_expl"] = pp5ext_expl

        # attach per-criterion explanations (already included in variant_row_with_disease but keep safe)
        for crit in ("PVS1","PS1","PM1","PM2","PM3","PM4","PM5","PP1","PP3","PP5",
                     "BA1","BS1","BS2","BP4","BP5","BP6","BP7"):
            expl_col = f"{crit}_expl"
            try:
                row[expl_col] = v.criteria_explanations.get(crit, "")
            except Exception:
                row[expl_col] = ""
        
        # attach Rules_applied string
        try:
            row["Rules_applied"] = ";".join(getattr(v, "_rules_applied", [])) if getattr(v, "_rules_applied", None) else ""
        except Exception:
            row["Rules_applied"] = ""

        # ensure list fields are stringified
        row["ps1_matches"] = json.dumps(getattr(v, "ps1_matches", []), ensure_ascii=False)
        row["pm5_matches"] = json.dumps(getattr(v, "pm5_matches", []), ensure_ascii=False)

        # auto_report_reasons remain as-is (strings)
        row["auto_report_reasons"] = ";".join(getattr(v, "auto_report_reasons", []))

        # COMMENTS vs MANUAL allocation:
        # - If variant is auto (_in_auto True) move informational manual_reasons into 'comments' field
        #   and leave 'manual_reasons' empty.
        # - If variant is manual (_in_manual True) keep manual_reasons as-is and comments empty.
        # - Otherwise (neither) keep manual_reasons as-is and comments empty.
        filtered_notes = ";".join(getattr(v, "filtered_reasons", []) or [])
        manual_notes = ";".join(getattr(v, "manual_reasons", []) or [])
        auto_notes = ";".join(getattr(v, "auto_report_reasons", []) or [])

        # synchronize sets -> flags 
        k = get_key(v)
        if k in auto_ids:
            v._in_auto = True; v._in_manual = False; v._filtered = False
        elif k in manual_ids:
            v._in_manual = True; v._in_auto = False; v._filtered = False

        from collections import OrderedDict

        def _dedupe_keep_order(items):
            # preserve order, remove exact duplicates and empty strings
            out = []
            seen = set()
            for it in items:
                if not it:
                    continue
                s = str(it).strip()
                if not s or s in seen:
                    continue
                seen.add(s)
                out.append(s)
            return out

        # partition manual_reasons into "manual-specific" vs "other" using keyword matching
        MANUAL_PAT = re.compile(r'(?i)\bmanual\b|delins|phasing|pairing|pairing review|phasing required|hq conflict|conflict|in trans|in cis|hemizygous|hemizygosity')

        def partition_manual_notes(notes_list):
            manual_specific = []
            other = []
            for note in _dedupe_keep_order(notes_list):
                if MANUAL_PAT.search(note):
                    manual_specific.append(note)
                else:
                    other.append(note)
            return manual_specific, other

        # collect raw lists from variant
        raw_manual = getattr(v, "manual_reasons", []) or []
        raw_filtered = getattr(v, "filtered_reasons", []) or []
        raw_auto = getattr(v, "auto_report_reasons", []) or []

        # dedupe each source while keeping order
        raw_manual = _dedupe_keep_order(raw_manual)
        raw_filtered = _dedupe_keep_order(raw_filtered)
        raw_auto = _dedupe_keep_order(raw_auto)

        # partition manual notes into manual-specific (explicit reasons for manual) and others
        manual_specific, manual_other = partition_manual_notes(raw_manual)

        # prepare combined other-notes array (non-manual info) preserving sensible order:
        # prefer: manual_other first (notes that came from manual list but aren't manual reasons),
        # then filtered notes, then auto notes — but dedupe overall
        combined_other = _dedupe_keep_order(manual_other + raw_filtered + raw_auto)

        # Now decide per triage state. Ensure consistent flags if id-sets used upstream:
        k = get_key(v)
        if k in auto_ids:
            v._in_auto = True; v._in_manual = False; v._filtered = False
        elif k in manual_ids:
            v._in_manual = True; v._in_auto = False; v._filtered = False
        # else keep existing flags

        # Fill row fields according to strict rules
        if getattr(v, "_filtered", False):
            # Filtered: filtered_reasons MUST be present
            row["filtered_reasons"] = ";".join(raw_filtered) if raw_filtered else "Filtered: reason not recorded"
            # Ensure manual/auto empty (per rules)
            row["manual_reasons"] = ""
            row["auto_report_reasons"] = ""
            # Put any other informational notes into comments (optional)
            row["comments"] = ";".join(combined_other) if combined_other else ""
        elif getattr(v, "_in_manual", False):
            # Manual: manual_reasons MUST be present (prefer manual_specific; fallback to raw_manual)
            if manual_specific:
                row["manual_reasons"] = ";".join(manual_specific)
            elif raw_manual:
                # if there were notes but none matched manual pattern, keep the deduped raw list
                row["manual_reasons"] = ";".join(raw_manual)
            else:
                row["manual_reasons"] = "Manual: reason not recorded"
            # filtered/auto empty
            row["filtered_reasons"] = ""
            row["auto_report_reasons"] = ""
            # comments: move other informational notes (e.g. PM3 downgrade, low-quality ClinVar, etc.)
            row["comments"] = ";".join(combined_other) if combined_other else ""
        elif getattr(v, "_in_auto", False):
            # Auto: auto_report_reasons MUST be present
            row["auto_report_reasons"] = ";".join(raw_auto) if raw_auto else "Auto: reason not recorded"
            row["manual_reasons"] = ""
            row["filtered_reasons"] = ""
            # comments: include any informative notes (both manual_other and filtered notes)
            row["comments"] = ";".join(combined_other) if combined_other else ""
        else:
            # Neutral: not filtered, not manual, not auto -> all triage fields empty; keep notes in comments
            row["filtered_reasons"] = ""
            row["manual_reasons"] = ""
            row["auto_report_reasons"] = ""
            row["comments"] = ";".join(_dedupe_keep_order(raw_manual + raw_filtered + raw_auto)) if (raw_manual or raw_filtered or raw_auto) else ""

        # if variant was selected for AUTO or MANUAL reporting -> it is not filtered out (No)
        # otherwise it is filtered out (Yes).
        if getattr(v, "_in_auto", False) or getattr(v, "_in_manual", False):
            row["filtered_out"] = "No"
        else:
            row["filtered_out"] = "Yes"

        rows_all.append(row)  

        if getattr(v, "_in_auto", False):
            rows_auto.append(row)
        if getattr(v, "_in_manual", False):
            rows_manual.append(row)

    # collapse duplicates
    auto_rows_c = collapse_output_rows(rows_auto)
    manual_rows_c = collapse_output_rows(rows_manual)
    all_rows_c = collapse_output_rows(rows_all)

    # Ensure we have rows for ALL samples, even if they have no reportable variants
    all_samples_in_data = set(v.sample for v in all_variants)
    samples_in_output = set(row.get("sample") for row in all_rows_c)

    missing_samples = all_samples_in_data - samples_in_output
    if missing_samples:
        logger.warning(f"Samples with no output rows: {missing_samples}")

    # write CSVs (pandas preferred)
    path_all = os.path.join(outdir, "all_candidates.csv")
    path_auto = os.path.join(outdir, "auto_conclusions.csv")
    path_manual = os.path.join(outdir, "manual_review_list.csv")

    if HAVE_PANDAS:
        df_all = pd.DataFrame(all_rows_c).reindex(columns=COLUMNS)
        df_auto = pd.DataFrame(auto_rows_c).reindex(columns=COLUMNS)
        df_manual = pd.DataFrame(manual_rows_c).reindex(columns=COLUMNS)
        df_all.to_csv(path_all, index=False)
        df_auto.to_csv(path_auto, index=False)
        df_manual.to_csv(path_manual, index=False)
    else:
        # fallback TSV writer
        def _write_tsv(path, rows, cols):
            if not rows:
                # create empty with header to be consistent
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write("\t".join(cols) + "\n")
                return
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("\t".join(cols) + "\n")
                for r in rows:
                    fh.write("\t".join(str(r.get(c,"")) for c in cols) + "\n")
        _write_tsv(path_all, all_rows_c, COLUMNS)
        _write_tsv(path_auto, auto_rows_c, COLUMNS)
        _write_tsv(path_manual, manual_rows_c, COLUMNS)

    def generate_ar_candidate_pairs_for_sample(varlist: List["VariantRecord"]) -> List[Tuple[str,str,Tuple]]:
        """
        From a list of VariantRecord objects for a sample/gene, return list of pairs (vid1, vid2, (vobj1_key,vobj2_key))
        vid format: chr:pos:ref:alt
        Follow AR pairing logic: take heterozygous P/LP as potential primary and choose partner candidate:
        - prefer another heterozygous P/LP
        - else strongest VUS with >=1 point
        Only pairs for these combos are returned.
        """
        pairs = []
        # build passengers as in triage (only heterozygous P/LP or VUS)
        passengers = [v for v in varlist if v.automated_class in ("Pathogenic","Likely pathogenic","VUS") and _is_het(v.proband_gt)]
        path_lp = [v for v in passengers if v.automated_class in ("Pathogenic","Likely pathogenic")]
        # iterate strong path variants
        for strongest in sorted(path_lp, key=lambda x: x.total_points, reverse=True):
            # find partner: other P/LP first, else VUS with >=1
            others = [p for p in passengers if p is not strongest]
            plp_partners = [p for p in others if p.automated_class in ("Pathogenic","Likely pathogenic")]
            if plp_partners:
                partner = sorted(plp_partners, key=lambda x: x.total_points, reverse=True)[0]
                pairs.append((f"{strongest.chrom}:{strongest.pos}:{strongest.ref}:{strongest.alt}",
                            f"{partner.chrom}:{partner.pos}:{partner.ref}:{partner.alt}",
                            (get_key(strongest), get_key(partner))))
                continue
            vus_partners = [p for p in others if p.automated_class == "VUS" and p.total_points >= 1]
            if vus_partners:
                partner = sorted(vus_partners, key=lambda x: x.total_points, reverse=True)[0]
                pairs.append((f"{strongest.chrom}:{strongest.pos}:{strongest.ref}:{strongest.alt}",
                            f"{partner.chrom}:{partner.pos}:{partner.ref}:{partner.alt}",
                            (get_key(strongest), get_key(partner))))
        return pairs

    return path_all, path_auto, path_manual, len(auto_rows_c), len(manual_rows_c)

# ---------------------------
# Build variant record from VCF
# ---------------------------
def build_variant_record_from_rec(rec, acmg_genes_normalized: set,
                                  dad_vf: Optional[pysam.VariantFile],
                                  mom_vf: Optional[pysam.VariantFile],
                                  exons_df: Optional[Any],
                                  proband_sample: Optional[str]=None,
                                  require_mane: bool=False,
                                  dbm: Optional[DatabaseManager]=None) -> Optional[VariantRecord]:
    """
    Robust builder for VariantRecord from cyvcf2 or pysam VariantRecord.

    Fixes applied:
      - corrected info_get calls (proper default argument placement)
      - defensive parsing of INFO/FORMAT values
      - ensures vr._raw_ann populated and vr.impact convenience attribute set
      - fetches parental GTs when pysam handles are provided
    """
    def _to_str(x):
        if x is None:
            return ""
        if isinstance(x, (bytes, bytearray)):
            try:
                return x.decode('utf-8', errors='ignore')
            except Exception:
                return str(x)
        try:
            if isinstance(x, (list, tuple)):
                if len(x) == 0:
                    return ""
                v = x[0]
                return _to_str(v)
        except Exception:
            pass
        return str(x)

    def _to_float(x):
        if x is None or (isinstance(x, str) and x == ""):
            return None
        try:
            if isinstance(x, (list, tuple)):
                if len(x) == 0:
                    return None
                x0 = x[0]
                return float(x0)
            if isinstance(x, (bytes, bytearray)):
                x = x.decode('utf-8', errors='ignore')
            return float(x)
        except Exception:
            try:
                s = _to_str(x)
                s = s.split(',')[0].strip()
                return float(s) if s != "" else None
            except Exception:
                return None

    def _to_int(x):
        if x is None or (isinstance(x, str) and x == ""):
            return None
        try:
            if isinstance(x, (list, tuple)):
                if len(x) == 0:
                    return None
                x0 = x[0]
                return int(x0)
            if isinstance(x, (bytes, bytearray)):
                x = x.decode('utf-8', errors='ignore')
            return int(x)
        except Exception:
            try:
                s = _to_str(x)
                if ',' in s:
                    s = s.split(',')[0]
                return int(float(s))
            except Exception:
                return None

    # Basic CHROM/POS/REF/ALT extraction
    try:
        if hasattr(rec, "CHROM"):  # cyvcf2
            chrom = rec.CHROM
            pos = int(rec.POS)
            ref = str(rec.REF)
            alts = list(rec.ALT) if rec.ALT is not None else []
            sample_names = [s for s in (rec.samples or [])]
        else:  # pysam
            chrom = rec.chrom
            pos = int(rec.pos)
            ref = str(rec.ref)
            alts = list(rec.alts) if rec.alts is not None else []
            sample_names = list(rec.samples.keys()) if hasattr(rec, "samples") else []
    except Exception:
        logger.debug("Failed to read basic variant coordinates from record", exc_info=True)
        return None

    if not alts:
        return None
    alt = str(alts[0])
    chrom_norm = normalize_chrom(chrom)

    # Extract INFO robustly
    info_raw = {}
    try:
        if hasattr(rec, "INFO"):
            try:
                info_raw = dict(rec.INFO)
            except Exception:
                info_raw = {}
                for k in rec.INFO.keys():
                    try:
                        info_raw[k] = rec.INFO.get(k)
                    except Exception:
                        info_raw[k] = None
        elif hasattr(rec, "info"):
            try:
                info_raw = dict(rec.info)
            except Exception:
                info_raw = {}
                try:
                    for k in rec.info.keys():
                        info_raw[k] = rec.info.get(k)
                except Exception:
                    info_raw = {}
    except Exception:
        info_raw = {}

    info_up = {}
    for k, v in (info_raw or {}).items():
        try:
            info_up[str(k).upper()] = v
        except Exception:
            try:
                info_up[repr(k).upper()] = v
            except Exception:
                continue

    def info_get(*keys, default=None):
        for k in keys:
            if k is None:
                continue
            v = info_up.get(k.upper())
            if v is not None and v != "":
                return v
        return default

    # Common annotation fields (with corrected calls)
    gene_raw = _to_str(info_get("GENE", "SYMBOL", "GENEINFO", "Gene", default="")) or ""
    hgvsc = _to_str(info_get("HGVSC", "HGV_SC", "HGV_S", "HGVSc", default="")) or ""
    hgvsp = _to_str(info_get("HGVSP", "HGV_SP", "HGV_P", "HGVSp", default="")) or ""
    transcript = _to_str(info_get("TRANSCRIPT", "TRANSCRIPT_SELECT", "MANE_SELECT", "CANONICAL", default="")) or ""
    
    tid = _to_str(info_get("tid", "TID", "TRANSCRIPT_ID", "transcript_id", "FEATURE", "feature", default="")) or ""
    
    if not tid and hgvsc:
        if ":" in hgvsc:
            tx_part = hgvsc.split(":", 1)[0]
            if tx_part.startswith(("ENST", "NM_", "XM_", "XR_", "NR_")):
                tid = tx_part
    consequence = _to_str(info_get("CONSEQUENCE", "CONSEQUENCE_VEP", "CONS", "Annotation", "Consequence", default="")) or ""

    revel = _to_float(info_get("REVEL_SCORE", "REVEL", default=None))
    dscore_value = None
    dscore_val = info_get("DSCORE", "Dscore", "D_SCORE", "D-SCORE", "DSPLICE_SCORE", "D_SPLICE", default=None)
    if dscore_val is not None:
        try:
            dscore_value = float(dscore_val)
        except (ValueError, TypeError):
            dscore_str = str(dscore_val).strip()
            dscore_match = re.search(r'(\d+\.?\d*)', dscore_str)
            if dscore_match:
                dscore_value = float(dscore_match.group(1))
    alpha_missense = None
    alpha_val = info_get("AlphaMissense", default=None)
    if alpha_val is not None and str(alpha_val).strip() and str(alpha_val).strip().lower() not in ["", ".", "na", "none", "null"]:
        try:
            alpha_missense = _to_float(alpha_val)
            logger.debug(f"Found AlphaMissense={alpha_missense} from key=AlphaMissense")
        except Exception:
            pass
    else:
        alpha_keys = ["ALPHAMISSENSE", "ALPHA_MISSENSE", "alpha_missense", 
                      "ALPHA-MISSENSE", "AlphaMissense_score", "ALPHAMISSENSE_SCORE"]
        for key in alpha_keys:
            val = info_get(key, default=None)
            if val is not None and str(val).strip() and str(val).strip().lower() not in ["", ".", "na", "none", "null"]:
                try:
                    alpha_missense = _to_float(val)
                    if alpha_missense is not None:
                        logger.debug(f"Found AlphaMissense={alpha_missense} from key={key}")
                        break
                except Exception:
                    continue
    
    spliceai = _to_float(info_get("SPLICE_AI_SCORE", "SPLICEAI", "splice_ai_score", "SpliceAI", default=None))
    # cadd = _to_float(info_get("CADD_PHRED", "CADD", "CADD_phred", default=None))

    squirls_score = None
    squirls_category = None
    
    # Try to get SQUIRLS from INFO (as category or score)
    squirls_val = info_get("SQUIRLS", default=None)
    if squirls_val is not None:
        try:
            squirls_score = float(squirls_val)
        except (ValueError, TypeError):
            squirls_category = str(squirls_val).strip().upper()
    
    # Try to get SQUIRLS_SCORE
    squirls_score_val = info_get("SQUIRLS_SCORE", default=None)
    if squirls_score_val is not None:
        try:
            squirls_score = float(squirls_score_val)
        except (ValueError, TypeError):
            pass

    gnomad_af = _to_float(info_get("GNOMAD_AF", "GNOMAD", "MAX_AF", "MAX_AF_POPS", "AF_GNOMAD", "max_freq", default=None))

    if gnomad_af is not None and gnomad_af < 0:
        gnomad_af = 0.0

    nmd_flag = _to_str(info_get("NMD", "NMDESC", "NMDesc", default="") or "")
    exon_field = _to_str(info_get("EXON", "EXON_NUMBER", "exon", "Exon", default="") or "")

    # Sample / FORMAT parsing
    sample_name = proband_sample
    if not sample_name:
        sample_name = (sample_names[0] if sample_names else "Proband")

    prob_gt = "./."
    prob_dp = None
    prob_ad = None
    prob_af = None
    prob_gq = None

    try:
        if hasattr(rec, "format") and callable(getattr(rec, "format")):
            # cyvcf2
            try:
                sample_list = []
                try:
                    sample_list = [s.split(':')[0] if isinstance(s, str) else str(s) for s in (rec.samples or [])]
                except Exception:
                    sample_list = []
                sample_idx = 0
                if sample_list and sample_name in sample_list:
                    sample_idx = sample_list.index(sample_name)
                try:
                    gt_arr = rec.format('GT')
                    if gt_arr is not None and len(gt_arr) > 0 and sample_idx < len(gt_arr):
                        gt_entry = gt_arr[sample_idx]
                        if isinstance(gt_entry, (list, tuple, bytes, bytearray)):
                            try:
                                alleles = []
                                for a in gt_entry:
                                    if a is None:
                                        alleles.append('.')
                                    else:
                                        alleles.append(str(int(a)))
                                prob_gt = f"{alleles[0]}/{alleles[1]}" if len(alleles) >= 2 else "/".join(alleles)
                            except Exception:
                                prob_gt = _to_str(gt_entry)
                        else:
                            prob_gt = _to_str(gt_entry)
                except Exception:
                    pass
                try:
                    dp_arr = rec.format('DP')
                    if dp_arr is not None and len(dp_arr) > 0 and sample_idx < len(dp_arr):
                        dp_val = dp_arr[sample_idx]
                        prob_dp = _to_int(dp_val)
                except Exception:
                    pass
                try:
                    gq_arr = rec.format('GQ')
                    if gq_arr is not None and len(gq_arr) > 0 and sample_idx < len(gq_arr):
                        gq_val = gq_arr[sample_idx]
                        prob_gq = _to_int(gq_val)
                except Exception:
                    pass
                try:
                    ad_arr = rec.format('AD')
                    if ad_arr is not None and len(ad_arr) > 0 and sample_idx < len(ad_arr):
                        ad_val = ad_arr[sample_idx]
                        if isinstance(ad_val, (list, tuple)):
                            try:
                                prob_ad = tuple(int(x) for x in ad_val)
                            except Exception:
                                prob_ad = tuple(_to_int(x) or 0 for x in ad_val)
                        else:
                            s = _to_str(ad_val)
                            if ',' in s:
                                parts = [p.strip() for p in s.split(',')]
                                try:
                                    prob_ad = (int(float(parts[0])), int(float(parts[1]))) if len(parts) >= 2 else None
                                except Exception:
                                    prob_ad = None
                except Exception:
                    pass
                try:
                    af_arr = rec.format('AF')
                    if af_arr is not None and len(af_arr) > 0 and sample_idx < len(af_arr):
                        af_val = af_arr[sample_idx]
                        prob_af = _to_float(af_val)
                except Exception:
                    pass
            except Exception:
                try:
                    samples_list = rec.samples or []
                    if samples_list:
                        idx = 0
                        if isinstance(samples_list[0], dict):
                            names = [s.get('name') if isinstance(s, dict) else str(s) for s in samples_list]
                            if sample_name in names:
                                idx = names.index(sample_name)
                            sdat = samples_list[idx]
                            prob_gt = _to_str(sdat.get('GT') or prob_gt)
                            prob_dp = _to_int(sdat.get('DP') or prob_dp)
                            ad = sdat.get('AD') or None
                            if ad:
                                if isinstance(ad, (list,tuple)):
                                    prob_ad = tuple(int(x) for x in ad)
                                else:
                                    s = _to_str(ad)
                                    if ',' in s:
                                        p = [int(float(x)) for x in s.split(',') if x.strip().isdigit()]
                                        prob_ad = tuple(p) if p else None
                except Exception:
                    pass
        else:
            # pysam
            try:
                sample_keys = list(rec.samples.keys()) if hasattr(rec, "samples") else []
                chosen = None
                if sample_name and sample_name in sample_keys:
                    chosen = sample_name
                elif sample_keys:
                    chosen = sample_keys[0]
                if chosen:
                    sdat = rec.samples[chosen]
                    try:
                        gt = sdat.get('GT')
                        if gt:
                            prob_gt = "/".join("." if x is None else str(x) for x in gt)
                    except Exception:
                        pass
                    try:
                        prob_dp = _to_int(sdat.get('DP'))
                    except Exception:
                        pass
                    try:
                        prob_gq = _to_int(sdat.get('GQ'))
                    except Exception:
                        pass
                    try:
                        ad = sdat.get('AD')
                        if ad is not None:
                            if isinstance(ad, (list, tuple)):
                                prob_ad = tuple(int(x) for x in ad)
                            else:
                                s = _to_str(ad)
                                if ',' in s:
                                    parts = [p.strip() for p in s.split(',')]
                                    try:
                                        prob_ad = (int(float(parts[0])), int(float(parts[1]))) if len(parts) >= 2 else None
                                    except Exception:
                                        prob_ad = None
                    except Exception:
                        pass
                    try:
                        afv = sdat.get('AF')
                        if afv is not None:
                            prob_af = _to_float(afv)
                    except Exception:
                        pass
            except Exception:
                pass
    except Exception:
        logger.debug("FORMAT parsing failed for record %s:%s", chrom, pos, exc_info=True)

    if prob_af is None:
        prob_af = _to_float(info_get("AF", "ALLELE_FREQ", "ALLELE_FREQUENCY", "AF_raw", "AF_1", "AF1"))

    if prob_dp is None:
        prob_dp = _to_int(info_get("DP", "DEPTH"))

    if prob_ad is None:
        raw_ad_info = info_get("AD", "ALLELIC_DEPTH")
        if raw_ad_info is not None:
            if isinstance(raw_ad_info, (list, tuple)):
                try:
                    prob_ad = tuple(int(x) for x in raw_ad_info)
                except Exception:
                    prob_ad = None
            else:
                s = _to_str(raw_ad_info)
                if ',' in s:
                    parts = [p.strip() for p in s.split(',')]
                    try:
                        prob_ad = (int(float(parts[0])), int(float(parts[1]))) if len(parts) >= 2 else None
                    except Exception:
                        prob_ad = None

    # Get SEX from INFO, handle empty values properly
    sample_sex_val = info_up.get("SEX") or info_up.get("Sex") or info_up.get("sex") or ""
    sample_sex_parsed = "Unknown"
    
    # Only process non-empty values (not empty string, not just semicolon, etc.)
    if sample_sex_val and str(sample_sex_val).strip() not in ("", ".", "NA", "null", "None", ";"):
        sex_lower = str(sample_sex_val).strip().lower()
        # Handle various representations including Russian (Cyrillic)
        # Male variants: English, Russian, abbreviations
        male_patterns = ("male", "м", "муж", "man", "m", "мужской", "мужчина", "boy", "man")
        # Female variants: English, Russian, abbreviations
        female_patterns = ("female", "ж", "жен", "woman", "f", "женский", "женщина", "girl", "woman")
        
        if sex_lower in male_patterns:
            sample_sex_parsed = "Male"
        elif sex_lower in female_patterns:
            sample_sex_parsed = "Female"
        elif any(pattern in sex_lower for pattern in ("male", "man", "м", "муж")):
            sample_sex_parsed = "Male"
        elif any(pattern in sex_lower for pattern in ("female", "woman", "ж", "жен")):
            sample_sex_parsed = "Female"
        else:
            logger.debug(f"Unrecognized SEX value: '{sample_sex_val}', keeping Unknown")
        logger.debug(f"Parsed SEX: '{sample_sex_val}' -> '{sample_sex_parsed}'")
    else:
        logger.debug(f"SEX field is empty or missing for sample {proband_sample}, will use inference later")

    vr = VariantRecord(
        chrom=chrom_norm,
        pos=int(pos),
        ref=str(ref),
        alt=str(alt),
        sample=sample_name,
        gene=normalize_gene_name(gene_raw) if gene_raw else "",
        gene_raw=gene_raw,
        transcript=transcript, 
        tid=tid,  
        consequence=consequence,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        exon=exon_field,
        is_last_exon=False,
        nmd=nmd_flag,
        sample_sex=sample_sex_parsed,
        gnomad_af=gnomad_af,
        gnomad_details={},
        clinvar_sig="",
        clinvar_trait="",
        clinvar_review="",
        clinvar_stars="0",
        internal_db_sig="",
        internal_db_details="",
        hgmd_id="",
        hgmd_class="",
        hgmd_phen="",
        hgmd_rankscore="Not provided",
        hgmd_support_positive=0,
        hgmd_support_total=0,
        pathogenic_criteria_count=0,
        benign_criteria_count=0,
        hgmd_publication_count=0,
        gnomad_hom=None,
        revel=revel,
        alpha_missense=alpha_missense,
        spliceai=spliceai,
        dscore=dscore_value,
        nmdesc_predictor=nmd_flag,
        proband_gt=prob_gt or "./.",
        father_gt="./.",
        mother_gt="./.",
        proband_dp=prob_dp,
        proband_ad=prob_ad,
        proband_af=prob_af,
    )

    try:
        vr.query_variant = f"{vr.chrom}-{vr.pos}-{vr.ref}-{vr.alt}"
    except Exception:
        vr.query_variant = f"{chrom}-{pos}-{ref}-{alt}"

    def _fetch_parental_gt(vf_handle, chrom_q, pos_q, ref_q):
        if vf_handle is None:
            return "./."
        try:
            candidates = []
            if chrom_q.startswith("chr"):
                candidates = [chrom_q, chrom_q.replace("chr", "")]
            else:
                candidates = [chrom_q, "chr" + chrom_q]
            for q in candidates:
                try:
                    for rr in vf_handle.fetch(q, pos_q - 1, pos_q):
                        if int(rr.pos) != int(pos_q):
                            continue
                        if getattr(rr, 'ref', None) != ref_q:
                            continue
                        smp_keys = list(rr.samples.keys()) if hasattr(rr, "samples") else []
                        if not smp_keys:
                            continue
                        sname = smp_keys[0]
                        sd = rr.samples[sname]
                        gt = sd.get("GT")
                        if gt:
                            return "/".join("." if x is None else str(x) for x in gt)
                        fmt_gt = sd.get("GT")
                        if fmt_gt:
                            return _to_str(fmt_gt)
                except Exception:
                    continue
        except Exception:
            return "./."
        return "./."

    try:
        vr.father_gt = _fetch_parental_gt(dad_vf, chrom_norm, vr.pos, vr.ref) if dad_vf else "./."
    except Exception:
        vr.father_gt = "./."
    try:
        vr.mother_gt = _fetch_parental_gt(mom_vf, chrom_norm, vr.pos, vr.ref) if mom_vf else "./."
    except Exception:
        vr.mother_gt = "./."

    raw_keys = ("GENE", "HGVSC", "HGVSP", "REVEL_SCORE", "REVEL", "ALPHA_MISSENSE",
                "ALPHAMISSENSE", "SPLICE_AI_SCORE", "SPLICEAI", "DSCORE", "NMD", "GNOMAD_AF", "MAX_AF")
    vr._raw_ann = {k: info_up.get(k) for k in raw_keys if k in info_up}

    try:
        if getattr(vr, "_raw_ann", None):
            for k in ("ALPHAMISSENSE", "ALPHA_MISSENSE", "AlphaMissense", "alpha_missense", "ALPHA-MISSENSE"):
                val = vr._raw_ann.get(k) if isinstance(vr._raw_ann, dict) else None
                if val in (None, "", "."):
                    continue
                try:
                    vr.alpha_missense = float(val)
                    break
                except Exception:
                    # sometimes value like "NA" or "Not_provided" — skip
                    try:
                        s = str(val).strip()
                        if s and s.lower() not in ("na", "not provided", "not_provided", ".", "null"):
                            vr.alpha_missense = float(s)
                            break
                    except Exception:
                        continue
    except Exception:
        vr.alpha_missense = getattr(vr, "alpha_missense", None)

    # Normalize Impact into vr._raw_ann and vr.impact convenience attribute
    try:
        impact_val = (info_up.get("IMPACT") or info_up.get("Impact") or info_up.get("impact") or "")
        impact_val = str(impact_val).strip() if impact_val is not None else ""
        vr._raw_ann["IMPACT"] = impact_val
        vr.impact = impact_val.upper() if impact_val else ""
    except Exception:
        vr.impact = ""

    try:
        aa_ref_3, aa_pos, aa_alt_3 = parse_hgvsp_three_letter(vr.hgvsp or "")
        if aa_ref_3:
            vr.aa_ref = aa_ref_3
        if aa_pos is not None:
            vr.aa_pos = aa_pos
        if aa_alt_3:
            vr.aa_alt = aa_alt_3
    except Exception:
        # preserve existing behavior on failure
        pass

    try:
        if dbm:
            dbm.get_variant_db_evidence(vr)
    except Exception:
        logger.debug("Database lookup failed for %s:%s", chrom_norm, pos, exc_info=True)

    return vr

# ---------------------------
# Pre-analysis - QC sex inference
# ---------------------------
def perform_pre_analysis_qc(vcf_path: str, proband_id: str) -> Dict[str, Any]:
    """
    Perform pre-analysis QC including sex inference.
    """
    warnings = []
    inferred_sex = "Unknown"
    
    try:
        if not HAVE_PYSAM:
            logger.warning("pysam not available for sex inference")
            return {"warnings": warnings, "sex": inferred_sex}
        
        vf = pysam.VariantFile(vcf_path)
        
        # Get the actual sample name from VCF header
        vcf_sample_names = list(vf.header.samples)
        if not vcf_sample_names:
            logger.warning(f"No samples found in VCF {vcf_path}")
            return {"warnings": warnings, "sex": inferred_sex}
        
        # Use the first sample from VCF (most VCFs have only one sample)
        sample_name = vcf_sample_names[0]
        logger.info(f"QC Sex Check: VCF sample name = '{sample_name}', proband_id = '{proband_id}'")
        
        # FIRST: Try to get sex from SEX field in INFO (non-empty)
        sex_from_info = None
        records_checked = 0
        
        for rec in vf.fetch():
            records_checked += 1
            
            # Try to get SEX from INFO
            sex_val = None
            for key in ['SEX', 'Sex', 'sex']:
                if key in rec.info:
                    val = rec.info.get(key)
                    if val is not None:
                        sex_val = str(val).strip()
                        break
            
            # Only accept non-empty values that actually contain sex information
            if sex_val and sex_val not in ('', '.', 'NA', 'null', 'None', ';'):
                sex_lower = sex_val.lower()
                logger.info(f"Found non-empty SEX='{sex_val}' at record {records_checked}")
                if sex_lower in ('male', 'm', 'man', 'муж', 'boy', 'мужской', 'мужчина'):
                    sex_from_info = "Male"
                    break
                elif sex_lower in ('female', 'f', 'woman', 'жен', 'girl', 'женский', 'женщина'):
                    sex_from_info = "Female"
                    break
            
            # Limit search
            if records_checked > 5000:
                break
        
        if sex_from_info:
            inferred_sex = sex_from_info
            logger.info(f"QC Sex Check: Sex from INFO field: {inferred_sex}")
            vf.close()
            return {"warnings": warnings, "sex": inferred_sex}
        
        # SECOND: Chromosome-based inference
        logger.info("QC Sex Check: No SEX field found, running chromosome inference...")
        
        # Check for Y chromosome variants using the correct sample name
        has_y = False
        try:
            for contig in vf.header.contigs:
                if contig in ('chrY', 'Y', '24'):
                    logger.info(f"Checking Y chromosome contig: {contig}")
                    for rec in vf.fetch(contig):
                        if sample_name in rec.samples:
                            gt = rec.samples[sample_name].get('GT')
                            if gt:
                                # Check if any ALT allele (not 0/0 or .)
                                alleles = [str(x) for x in gt if x is not None]
                                if any(a not in ('0', '.') for a in alleles):
                                    has_y = True
                                    logger.info(f"Found Y variant at {rec.contig}:{rec.pos} with GT={gt}")
                                    break
                        if has_y:
                            break
                    if has_y:
                        break
        except Exception as e:
            logger.debug(f"Y chromosome check failed: {e}")
        
        if has_y:
            inferred_sex = "Male"
            logger.info(f"QC Sex Check: Inferred Male (has Y chromosome variants)")
            vf.close()
            return {"warnings": warnings, "sex": inferred_sex}
        
        # Check X chromosome heterozygosity for Female inference
        het_count = 0
        total_x = 0
        
        try:
            for contig in vf.header.contigs:
                if contig in ('chrX', 'X', '23'):
                    logger.info(f"Checking X chromosome contig: {contig}")
                    for rec in vf.fetch(contig):
                        if sample_name in rec.samples:
                            gt = rec.samples[sample_name].get('GT')
                            if gt and len(gt) >= 2:
                                a0 = str(gt[0]) if gt[0] is not None else '0'
                                a1 = str(gt[1]) if gt[1] is not None else '0'
                                if a0 != a1:
                                    het_count += 1
                                total_x += 1
                        if total_x > 1000:  # Limit for performance
                            break
        except Exception as e:
            logger.debug(f"X chromosome check failed: {e}")
        
        if total_x > 0:
            ratio = het_count / total_x
            logger.info(f"X chromosome stats: total={total_x}, het={het_count}, ratio={ratio:.3f}")
            if ratio > 0.15:
                inferred_sex = "Female"
                logger.info(f"QC Sex Check: Inferred Female (high X heterozygosity: {ratio:.3f})")
                vf.close()
                return {"warnings": warnings, "sex": inferred_sex}
            elif ratio < 0.05 and total_x > 10:
                inferred_sex = "Male"
                logger.info(f"QC Sex Check: Inferred Male (low X heterozygosity: {ratio:.3f})")
                vf.close()
                return {"warnings": warnings, "sex": inferred_sex}
        
        # If we reach here, sex remains Unknown
        logger.warning(f"QC Sex Check: Could not infer sex, keeping Unknown")
        
    except Exception as e:
        logger.warning(f"QC Sex Check failed: {e}")
        import traceback
        traceback.print_exc()
    
    logger.info(f"QC Sex Check: Final sex = {inferred_sex}")
    return {"warnings": warnings, "sex": inferred_sex}

# ---------------------------
# Utilities for class computation
# ---------------------------
def compute_class_from_points(criteria_points: Dict[str,int]) -> Tuple[str,int,int]:
    total = sum(int(v) for v in (criteria_points or {}).values())
    path_cnt = sum(1 for k in (criteria_points or {}) if any(p in k for p in ["PVS","PS","PM","PP"]))
    auto_class = "VUS"
    if total >= 10:
        auto_class = "Pathogenic"
    elif total >= 6:
        auto_class = "Likely pathogenic"
    elif total <= -6:
        auto_class = "Benign"
    elif total <= -1:
        auto_class = "Likely benign"
    else:
        auto_class = "VUS"
    if auto_class == "Pathogenic" and path_cnt <= 2:
        auto_class = "Likely pathogenic"

    if (auto_class.startswith("Path") or auto_class.startswith("Likely path")) and path_cnt < MIN_PATHOGENIC_CRITERIA:
        auto_class = "VUS"
        
    return auto_class, total, path_cnt

# ---------------------------
# Process VCF 
# ---------------------------
def process_vcf(proband_vcf: str,
                father_vcf: Optional[str],
                mother_vcf: Optional[str],
                outdir: str,
                acmg_table: str,
                db_paths: Dict[str, str],
                exons_file: Optional[str],
                gene_bed: Optional[str] = None,
                proband_sample: Optional[str] = None,
                strict_qc: bool = True,
                gnomad_v2: Optional[str] = None,
                gnomad_v3: Optional[str] = None,
                acmg_recs: Optional[str] = None,
                internal_db_path: Optional[str] = None,
                gnomad_v4: Optional[str] = None,
                family_id: str = "Unknown",
                require_mane: bool = False,
                #cooccurrence_path: Optional[str] = None, 
                gnomad_api_endpoint: Optional[str] = None, 
                gnomad_dataset: Optional[str] = None, 
                transcript_map: Optional[str] = None, 
                nmd_table: Optional[str] = None):
    os.makedirs(outdir, exist_ok=True)
    input_vcf_for_tools = proband_vcf
    if not proband_vcf.endswith(".gz") and not proband_vcf.endswith(".bgz"):
        logger.info("Input VCF is plain text. Compressing and indexing for tool compatibility...")
        compressed_vcf = os.path.join(outdir, "input_prepped.vcf.gz")
        if shutil.which("bgzip"):
            subprocess.run(f"bgzip -c '{proband_vcf}' > '{compressed_vcf}'", shell=True, check=True)
        else:
            with open(proband_vcf, 'rb') as f_in, gzip.open(compressed_vcf, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        subprocess.run(["tabix", "-p", "vcf", compressed_vcf], check=False)
        input_vcf_for_tools = compressed_vcf
    elif not os.path.exists(proband_vcf + ".tbi") and not os.path.exists(proband_vcf + ".csi"):
        try: subprocess.run(["tabix", "-p", "vcf", proband_vcf], check=False)
        except: pass
    working_vcf = input_vcf_for_tools
    annotated_vcf = working_vcf
    recs_map = load_recommendations_table(acmg_recs) if acmg_recs else {}
    qc_res = perform_pre_analysis_qc(input_vcf_for_tools, proband_sample)
    inferred_sex = qc_res.get("sex", "Unknown")
    logger.info(f"Pre-analysis QC complete. Sex: {inferred_sex}")
    if gene_bed and os.path.exists(gene_bed):
        logger.info("Gene BED provided: %s. Expanding by 20kb and merging intervals...", gene_bed)
        expanded_bed = os.path.join(outdir, "expanded_gene_regions.bed")
        intervals = defaultdict(list)
        with open(gene_bed) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#"): continue
                parts = re.split(r"\s+", ln)
                if len(parts) < 4: continue
                chrom = parts[0]
                try:
                    s = max(0, int(parts[1]) - 20000); e = int(parts[2]) + 20000
                except Exception:
                    continue
                intervals[chrom].append((s,e))
        with open(expanded_bed, "w") as out:
            for chrom in sorted(intervals.keys()):
                lst = sorted(intervals[chrom])
                cur_s, cur_e = lst[0]; merged=[]
                for s,e in lst[1:]:
                    if s <= cur_e + 1: cur_e = max(cur_e,e)
                    else:
                        merged.append((cur_s,cur_e)); cur_s,cur_e = s,e
                merged.append((cur_s,cur_e))
                for s,e in merged:
                    out.write(f"{chrom}\t{s}\t{e}\n")
        filtered_vcf = os.path.join(outdir, os.path.basename(proband_vcf).replace(".vcf","").replace(".gz","") + ".acmg_only.vcf.gz")
        if shutil.which("bcftools") is None:
            raise RuntimeError("bcftools not found; gene-bed filtering requires bcftools")
        run_cmd(["bcftools", "view", "-R", expanded_bed, input_vcf_for_tools, "-O", "z", "-o", filtered_vcf])
        try:
            run_cmd(["tabix","-p","vcf", filtered_vcf], check=False)
        except Exception:
            logger.debug("tabix indexing failed")
        working_vcf = filtered_vcf; annotated_vcf = working_vcf
    working_vcf = annotated_vcf

    # Set path to FASTA
    fasta_path = "~/ngs-data/ngs/ref/hg38.fa" 
    fasta_handle = None
    if fasta_path and os.path.exists(fasta_path) and HAVE_PYSAM:
        try:
            fasta_handle = pysam.FastaFile(fasta_path)
        except Exception:
            fasta_handle = None
    gene_rules = load_acmg_table(acmg_table) if HAVE_PANDAS else {}
    # FORCE: Mark all genes as reportable (VCF is already filtered to ACMG genes)
    for gene in gene_rules:
        gene_rules[gene]["reportable"] = True
    logger.info(f"Forced reportable=True for all {len(gene_rules)} ACMG genes")
    acmg_genes_normalized = set(g for g in gene_rules.keys() if gene_rules.get(g,{}).get("reportable", False))
    dbm = DatabaseManager(db_paths, gnomad_v2=gnomad_v2, gnomad_v3=gnomad_v3, gnomad_v4=gnomad_v4, internal_db_path=internal_db_path)
    exons_df = load_exons_table(exons_file) if exons_file else None
    # load optional transcript mapping and NMD precomputed table (if provided)
    tx_map = {}
    if transcript_map:
        try:
            tx_map = load_transcript_map(transcript_map)
        except Exception as e:
            logger.warning("Failed to load transcript map %s: %s", transcript_map, e)
            tx_map = {}

        # Load NMD map (if provided)
    nmd_map = {}
    if nmd_table:
        try:
            nmd_map = load_nmd_table(nmd_table)
            # populate genomic exon coords in nmd_map from exons_df if possible
            try:
                _populate_genomic_exon_coords(nmd_map, exons_df if exons_df is not None else {})
            except Exception:
                pass
        except Exception as e:
            logger.warning("Failed to load NMD table %s: %s", nmd_table, e)
            nmd_map = {}

    # ---------------------------
    # VCF parsing: MONO-MODE (each sample independent)
    # ---------------------------
    candidates: List[VariantRecord] = []
    logger.info("Parsing VCFs in mono-mode (each sample independent)...")

    def _basename_no_vcf_ext(path):
        base = os.path.basename(path)
        for ext in (".vcf.gz", ".vcf.bgz", ".vcf", ".bgz", ".gz"):
            if base.lower().endswith(ext):
                return base[:-len(ext)]
        # fallback: strip single extension if any
        return os.path.splitext(base)[0]

    vcf_samples = []

    proband_name = proband_sample if proband_sample else _basename_no_vcf_ext(working_vcf)

    vcf_samples.append({
        "path": working_vcf,
        "name": proband_name,
        "role": "proband"
    })

    if father_vcf and os.path.exists(father_vcf):
        fname = get_sample_names_from_vcf_header(father_vcf)
        vcf_samples.append({
            "path": father_vcf,
            "name": fname[0] if fname else "Father",
            "role": "father"
        })

    if mother_vcf and os.path.exists(mother_vcf):
        mname = get_sample_names_from_vcf_header(mother_vcf)
        vcf_samples.append({
            "path": mother_vcf,
            "name": mname[0] if mname else "Mother",
            "role": "mother"
        })

    logger.info(f"Will process {len(vcf_samples)} sample(s): {[s['name'] for s in vcf_samples]}")

    all_candidates: List[VariantRecord] = []

    for vcf_info in vcf_samples:
        logger.info(f"Parsing {vcf_info['role']} ({vcf_info['name']}) from {vcf_info['path']}")
        
        vcf_reader = None
        if HAVE_PYSAM:
            try:
                vcf_reader = pysam.VariantFile(vcf_info["path"])
            except Exception as e:
                logger.warning(f"Failed to open {vcf_info['path']} with pysam: {e}")
        
        if vcf_reader is None and HAVE_CYVCF2:
            try:
                from cyvcf2 import VCF
                vcf_reader = VCF(vcf_info["path"])
            except Exception as e:
                logger.warning(f"Failed to open {vcf_info['path']} with cyvcf2: {e}")
        
        if vcf_reader is None:
            logger.error(f"Cannot open {vcf_info['path']} with any library")
            continue
        
        iterator = vcf_reader.fetch() if hasattr(vcf_reader, 'fetch') else vcf_reader
        
        record_count = 0
        for rec in iterator:
            record_count += 1
            try:
                vr = build_variant_record_from_rec(
                    rec=rec,
                    acmg_genes_normalized=acmg_genes_normalized,
                    dad_vf=None,
                    mom_vf=None,
                    exons_df=exons_df,
                    proband_sample=vcf_info["name"],
                    require_mane=require_mane,
                    dbm=dbm
                )
                if not vr:
                    continue
                
                # Only override sex if it wasn't set from INFO (i.e., still Unknown)
                # This preserves sex information from VCF when available
                if vr.sample_sex == "Unknown" or not vr.sample_sex:
                    vr.sample_sex = inferred_sex
                else:
                    # Log if there's a discrepancy between INFO and inference
                    if inferred_sex != "Unknown" and vr.sample_sex != inferred_sex:
                        logger.debug(f"Sex discrepancy for {vr.sample}: INFO says {vr.sample_sex}, inference says {inferred_sex} - keeping INFO value")
                
                vr.sample_role = vcf_info["role"]
                
                if (vr.gnomad_af is None or vr.gnomad_af == 0.0) and dbm is not None:
                    try:
                        dbm.annotate_gnomad_for_variant(vr)
                    except Exception:
                        pass
                
                # DB evidence
                try:
                    if dbm is not None:
                        vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
                except Exception:
                    pass
                
                all_candidates.append(vr)
                
            except Exception as e:
                logger.debug(f"Skipping record {record_count} in {vcf_info['name']}: {e}")
                continue
        
        if hasattr(vcf_reader, 'close'):
            try:
                vcf_reader.close()
            except:
                pass
        
        logger.info(f"  {vcf_info['name']}: parsed {record_count} records, {len([c for c in all_candidates if c.sample == vcf_info['name']])} candidates")

    candidates = deduplicate_candidates(all_candidates)
    logger.info(f"Total candidates after dedup: {len(candidates)}")

    for v in candidates:
        if not v.sample:
            v.sample = proband_sample or "Proband"
    
    if len(candidates) == 0:
        logger.error("No variants were parsed from VCF. Check file format.")
        return
    
    logger.info(f"Total candidates after parsing: {len(candidates)}")

    # Deduplicate, set sample names, attach engine resources
    candidates = deduplicate_candidates(candidates)
    if proband_sample:
        for v in candidates:
            v.sample = proband_sample

    # Correct ACMGEngine instantiation with named args
    engine = ACMGEngine(
        dbm=dbm,
        gene_rules=gene_rules,
        exons_df=exons_df,
        spliceai_thr=THRESH.get("SPLICEAI_MODERATE", 0.20),
        require_mane=require_mane,
        nmd_map=nmd_map
    )
    # Attach transcript/exon/NMD resources to engine for PVS1/NMD decision
    engine.tx_map = tx_map
    engine.exons_map = exons_df
    engine.nmd_map = nmd_map
    # pass gnomAD API options into engine/global space for on-the-fly cooccurrence queries if needed
    engine.gnomad_api_endpoint = gnomad_api_endpoint or "https://gnomad.broadinstitute.org/api"
    engine.gnomad_dataset = gnomad_dataset or "gnomad_r2_1"

    # Initial ACMG evaluation (without PM3)
    for v in candidates:
        assigned, pts, auto_class, manual_flag, manual_reasons = engine.evaluate_variant(v, fasta_handle=fasta_handle)
        v.criteria_assigned = assigned
        v.criteria_points = pts
        v.automated_class = auto_class
        v.manual_review = manual_flag
        v.manual_reasons = v.manual_reasons + manual_reasons
        v.total_points = sum(pts.values())
        v.pathogenic_criteria_count = sum(1 for x in assigned if any(p in x for p in ["PVS","PS","PM","PP"]))
        v.benign_criteria_count = sum(1 for x in assigned if any(p in x for p in ["BA","BS","BP"]))
        
        # Enforce minimum pathogenic criteria count (Tavtigian rule)
        if "BA1" not in v.criteria_assigned:
            if (v.automated_class.startswith("Path") or v.automated_class.startswith("Likely path")) and v.pathogenic_criteria_count < MIN_PATHOGENIC_CRITERIA:
                v.automated_class = "VUS"
                v.manual_reasons.append(f"Insufficient pathogenic criteria (<{MIN_PATHOGENIC_CRITERIA})")

    if fasta_handle: 
        fasta_handle.close()

    # Apply PM3 batch scoring (phasing-aware)
    engine.apply_pm3_batch(candidates)

    for v in candidates:
        # initialize criteria_points_ext from base criteria_points if missing,
        # while removing legacy PP5/BP5 to avoid double counting
        if not hasattr(v, "criteria_points_ext") or v.criteria_points_ext is None:
            try:
                base_ext = dict(v.criteria_points or {})
                for kk in list(base_ext.keys()):
                    if kk.startswith("PP5") or kk.startswith("BP5"):
                        base_ext.pop(kk, None)
                v.criteria_points_ext = base_ext
            except Exception:
                v.criteria_points_ext = dict(v.criteria_points or {})

        # copy PM3 entries into criteria_points_ext (PM3 scoring computed earlier)
        try:
            for k, val in (v.criteria_points or {}).items():
                if k.startswith("PM3"):
                    # ensure integer
                    try:
                        v.criteria_points_ext[k] = int(val)
                    except Exception:
                        try:
                            v.criteria_points_ext[k] = int(float(val))
                        except Exception:
                            v.criteria_points_ext[k] = POINTS_MAP.get(k, 0)
        except Exception:
            pass

        # Recalculate extended totals and class
        try:
            v.total_ext_points = sum(int(x) for x in (v.criteria_points_ext or {}).values()) if (v.criteria_points_ext) else 0
        except Exception:
            v.total_ext_points = 0
        try:
            v.extended_class = compute_class_from_points(v.criteria_points_ext)[0] if (v.criteria_points_ext) else "VUS"
        except Exception:
            v.extended_class = "VUS"

        # Adopt extended_class as the canonical automated classification downstream
        v.automated_class = v.extended_class
        # Keep legacy total_points synchronized (optional)
        v.total_points = int(v.total_ext_points or 0)
        v.pathogenic_criteria_count = sum(1 for k in (v.criteria_points_ext or {}) if any(p in k for p in ["PVS","PS","PM","PP"]))

    # Delins cluster detection
    detect_delins(candidates)

    # =========================================================================
    # Step 6: Strict QC filtering with chromosome-specific thresholds
    # =========================================================================
    
    filtered_for_report = []
    dropped_count = 0
    
    for v in candidates:
        if not STRICT_QC_ENABLED or not strict_qc:
            filtered_for_report.append(v)
            continue
        
        dp = int(v.proband_dp or 0)
        af = float(v.proband_af or 0.0)
        
        # Determine chromosome type for threshold adjustment
        chrom_norm = str(v.chrom or "").upper().replace("CHR", "")
        is_sex_chrom = chrom_norm in ("X", "Y")
        is_autosome = not is_sex_chrom and chrom_norm not in ("M", "MT")
        
        # Check zygosity and hemizygosity
        is_homo = _is_hom(getattr(v, "proband_gt", "./."))
        is_hemi = False
        sample_sex = str(getattr(v, "sample_sex", "")).lower()
        
        # For males on sex chromosomes, treat as hemizygous
        if sample_sex == "male" and is_sex_chrom:
            gt = getattr(v, "proband_gt", "") or ""
            # If any alt allele present, it's hemizygous
            if "1" in gt:
                is_hemi = True
                is_homo = False  # Override homozygous for sex chromosomes in males
        
        # Determine appropriate DP threshold
        if is_hemi:
            # Hemizygous variants (male sex chromosomes) need lower threshold
            min_dp = 5
            threshold_type = "hemizygous"
        elif is_homo:
            # Homozygous variants need lower threshold
            min_dp = 5
            threshold_type = "homozygous"
        else:
            # Heterozygous variants need standard threshold
            min_dp = STRICT_MIN_DP  # usually 10
            threshold_type = "heterozygous"
        
        # Special case: Y chromosome in males
        if chrom_norm == "Y" and sample_sex == "male":
            # Y chromosome variants are always hemizygous in males
            min_dp = 5
            threshold_type = "y_chromosome"

        if chrom_norm == "X" and sample_sex == "male":
            # X chromosome variants are always hemizygous in males
            min_dp = 5
            threshold_type = "x_chromosome"
        
        if dp >= min_dp and af >= STRICT_MIN_AF:
            filtered_for_report.append(v)
        else:
            reason = (f"Failed strict QC: {threshold_type} variant on {v.chrom} "
                     f"with DP={dp} (min={min_dp}), AF={af:.3f} (min={STRICT_MIN_AF})")
            if not hasattr(v, "filtered_reasons") or v.filtered_reasons is None:
                v.filtered_reasons = []
            if reason not in v.filtered_reasons:
                v.filtered_reasons.append(reason)
            dropped_count += 1
    
    logger.info(f"After strict QC filtering: {len(filtered_for_report)} kept for Report, {dropped_count} low quality (audit only)")


    all_path, auto_path, manual_path, cnt_auto, cnt_manual = strict_triage_and_output(
        report_candidates=filtered_for_report,
        audit_candidates=candidates,
        gene_rules=gene_rules,
        recs_map=recs_map,
        outdir=outdir,
        family_id_input=family_id,
        dbm=dbm,                        
        #cooccur_path=cooccurrence_path
    )
    run_info = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "proband_vcf": proband_vcf,
        "annotated_vcf": annotated_vcf,
        "candidates_pre_qc": len(candidates),
        "candidates_post_qc": len(filtered_for_report),
        "dropped_pre_qc": dropped_count,
        "outputs": {"all": all_path, "auto": auto_path, "manual": manual_path},
        "counts": {"auto": cnt_auto, "manual": cnt_manual}
    }
    with open(os.path.join(outdir, "run_info.json"), "w") as fh:
        json.dump(run_info, fh, indent=2)
    logger.info("Finished. Outputs: %s, %s, %s", all_path, auto_path, manual_path)
    print("\n" + "="*60)
    print("ACMG PIPELINE FINISHED SUCCESSFULLY")
    print(f"Total candidates processed (ACMG genes): {len(candidates)}")
    print("-" * 30)
    print(f"VARIANTS FOR AUTOMATIC REPORT: {cnt_auto}")
    print(f"VARIANTS FOR MANUAL REVIEW:    {cnt_manual}")
    print("-" * 30)
    print(f"Output directory: {outdir}")
    print("="*60 + "\n")

# ---------------------------
# CLI
# ---------------------------
def load_acmg_table(path: str) -> Dict[str, Dict[str,Any]]:
    """
    Loads the ACMG SF table (ACMG gene list). Aggressive column detection allows flexible headers.
    Returns normalized gene keys.
    """
    rules: Dict[str, Dict[str,Any]] = {}
    if not path or not os.path.exists(path):
        logger.error("ACMG table not found: %s", path)
        return rules
    if not HAVE_PANDAS:
        logger.error("pandas required to load ACMG table")
        return rules
    try:
        df = pd.read_csv(path, sep=None, engine='python', dtype=str, on_bad_lines='skip').fillna("")
    except Exception as e:
        logger.warning("Failed to read ACMG table %s: %s", path, e)
        return rules
    cols = {c.lower().strip(): c for c in df.columns}
    def find_col_fuzzy(keywords):
        for k in keywords:
            if k in cols: return cols[k]
        for k in keywords:
            for real_col_lower in cols.keys():
                if k in real_col_lower:
                    return cols[real_col_lower]
        return None
    gene_col = find_col_fuzzy(["gene symbol", "gene", "symbol", "hgnc"])
    disease_col = find_col_fuzzy(["disease/phentyope", "phentyope", "disease", "phenotype", "condition", "disorder"])
    moi_col = find_col_fuzzy(["inheritance", "mode of inheritance", "moi"])
    variants_col = find_col_fuzzy(["variants to report", "variants", "variants_note", "note", "comments"])
    logger.info(f"MAPPED COLUMNS: Gene='{gene_col}', Disease='{disease_col}', MOI='{moi_col}'")
    if gene_col is None:
        logger.error("CRITICAL: Could not find 'Gene' column in ACMG table.")
        return rules
    any_reportable = False
    
    processed_genes = set()
    for _, r in df.iterrows():
        gene_raw = str(r.get(gene_col, "")).strip()
        if not gene_raw:
            continue
        
        gene_norm = normalize_gene_name(gene_raw)
        
        if gene_norm in processed_genes:
            # MERGE duplicate gene entries - BUT be careful with rules_text!
            if gene_norm in rules:
                existing_disease = rules[gene_norm].get("disease", "")
                new_disease = str(r.get(disease_col, "")).strip() if disease_col else ""
                if new_disease and new_disease not in existing_disease:
                    rules[gene_norm]["disease"] = f"{existing_disease}; {new_disease}"
                
                existing_moi = rules[gene_norm].get("moi", "")
                new_moi = str(r.get(moi_col, "")).strip() if moi_col else ""
                if new_moi and new_moi not in existing_moi:
                    rules[gene_norm]["moi"] = f"{existing_moi}; {new_moi}"
                
                existing_variants_note = rules[gene_norm].get("variants_note", "")
                new_variants_note = str(r.get(variants_col, "")).strip() if variants_col else ""
                if new_variants_note and new_variants_note not in existing_variants_note:
                    rules[gene_norm]["variants_note"] = f"{existing_variants_note}; {new_variants_note}"
                
                # CRITICAL: Do NOT merge rules_text - each gene should have its own rules!
                # If merging, keep the first non-empty rules_text, don't concatenate
                existing_rules = rules[gene_norm].get("rules_text", "")
                new_rules = str(r.get(rules_col_name, "")).strip() if rules_col_name else ""
                if new_rules and not existing_rules:
                    rules[gene_norm]["rules_text"] = new_rules
                # If both have rules, this is suspicious - log warning
                elif new_rules and existing_rules and new_rules != existing_rules:
                    logger.warning(f"Gene {gene_norm} has conflicting rules_text from multiple rows. Keeping first: '{existing_rules}'")
            continue
        
        processed_genes.add(gene_norm)
        
        moi = str(r.get(moi_col, "")).strip() if moi_col else ""
        disease = str(r.get(disease_col, "")).strip() if disease_col else ""
        disease = disease.replace("\r", " ").replace("\n", " ").strip()
        variants_note = str(r.get(variants_col, "")).strip() if variants_col else ""
        reportable = False
        found_flag_col = False
        for potential_col in cols:
            if "sf list" in potential_col or "version" in potential_col or "report" in potential_col:
                val = str(r.get(cols[potential_col],"")).strip()
                if val and val not in ("0", "-", ""):
                     reportable = True
                     found_flag_col = True
                     break
        if not found_flag_col and variants_col:
             if variants_note and len(variants_note) > 1: reportable = True
        special_raw = variants_note or ""
        text_blob = (special_raw + " " + moi + " " + disease).lower()
        flags = {
            "require_biallelic": "2 variants" in text_blob or "biallelic" in text_blob or "recessive" in text_blob or "2 het" in text_blob or "p and lp (2 variants)" in text_blob,
            "report_truncating_only": "truncating variants only" in text_blob or ("truncating" in text_blob and "only" in text_blob),
            "report_missense_only": "missense only" in text_blob,
            "lof_not_reportable": "not be reported as sfs" in text_blob or "gain-of-function" in text_blob or "gof" in text_blob,
            "disease_requires_confirmation": bool(disease.strip())
        }
        # Try to detect a 'Rules' column (case-insensitive)
        rules_col_name = None
        for possible in ("rules", "rule", "Rules", "RULES"):
            if possible in df.columns:
                rules_col_name = possible; break
        # If not present, try fuzzy match
        if not rules_col_name:
            for c in df.columns:
                if 'rule' in c.lower():
                    rules_col_name = c
                    break

        rules_text = ""
        if rules_col_name:
            try:
                rules_text = str(r.get(rules_col_name, "")).strip()
            except Exception:
                rules_text = ""

        rules[gene_norm] = {
            "reportable": reportable,
            "moi": moi,
            "disease": disease,
            "variants_note": variants_note,
            "special_flags": flags,
            "rules_text": rules_text,
            "raw_row": {c: str(r.get(c, "")) for c in df.columns}
        }
        if reportable:
            any_reportable = True
    
    if processed_genes:
        logger.info(f"Unique genes processed: {len(processed_genes)}")
        logger.info(f"Original rows in the table: {len(df)}")
        logger.info(f"Duplicate entries merged: {len(df) - len(processed_genes)}")
        
        genes_with_multiple = [g for g in processed_genes if ";" in rules.get(g, {}).get("disease", "")]
        if genes_with_multiple:
            for gene in genes_with_multiple[:5]:  
                disease = rules[gene].get("disease", "")
    
    if (not any_reportable) and len(rules) > 0:
        logger.warning("No explicit reportable flags detected. Marking ALL genes in table as reportable.")
        for g in rules: rules[g]["reportable"] = True
    logger.info("Loaded ACMG rules for %d genes from %s", len(rules), path)


    return rules

def build_parser():
    p = argparse.ArgumentParser(prog="acmg_sf_classifier.py", description="ACMG SF classifier")
    p.add_argument("--proband", help="Proband VCF (bgz recommended)")
    p.add_argument("--father", help="Father VCF (optional)")
    p.add_argument("--mother", help="Mother VCF (optional)")
    p.add_argument("--outdir", default="acmg_sf_results")
    p.add_argument("--acmg-table", required=True, help="ACMG SF table CSV/TSV")
    p.add_argument("--db-paths-json", help="JSON file to override DB paths")
    p.add_argument("--exons-file", help="Transcript/exons CSV (optional)")
    p.add_argument("--require-mane", action="store_true", help="Require MANE or MANE Plus Clinical annotation")
    p.add_argument("--gene-bed", help="BED file with gene regions to prefilter (optional). Intervals expanded by 20kb and merged.")
    p.add_argument("--clinvar-tabix", help="tabix-indexed ClinVar TSV (optional)")
    p.add_argument("--proband-sample", help="If VCF is multi-sample, specify proband sample name (optional)")
    p.add_argument("--no-strict-qc", action="store_true", help="Disable strict QC filtering (not recommended)")
    p.add_argument("--gnomad-v2", help="gnomAD v2 site file (VCF or TSV bgz)")
    p.add_argument("--gnomad-v3", help="gnomAD v3 site file (VCF or TSV bgz)")
    p.add_argument("--gnomad-v4", help="gnomAD v4 site VCF (tabix'd)")
    p.add_argument("--batch-input-dir", help="Folder containing subfolders with VCF files. Runs pipeline for ALL samples and aggregates results.")
    p.add_argument("--acmg-recs", help="Secondary CSV with specific variant/gene recommendations")
    p.add_argument("--internal-db", help="Path to internal DB CSV (must be hg38 coordinates)")
    p.add_argument("--hgmd", help="Path to HGMD Pro CSV (hg38)")
    #p.add_argument("--gnomad-cooccurrence", dest="gnomad_cooccurrence", help="Path to TSV with precomputed variant cooccurrence/cis-trans evidence (optional). Format: var1_vid\\tvar2_vid\\tcis_flag\\tevidence")
    p.add_argument("--gnomad-api-endpoint", dest="gnomad_api_endpoint",
                   default="https://gnomad.broadinstitute.org/api",
                   help="gnomAD GraphQL API endpoint (default public endpoint used by site)")
    p.add_argument("--gnomad-dataset", dest="gnomad_dataset",
                   default="gnomad_r2_1",
                   help="gnomAD dataset id for GraphQL queries (default gnomad_r2_1)")
    p.add_argument("--transcript-map", dest="transcript_map",
                   help="CSV/TSV with Ensembl <-> RefSeq mapping (columns: Ensembl_rna_identifier, RNA_nucleotide_accession.version, mane_select)")
    p.add_argument("--nmd-table", dest="nmd_table",
                   help="Precomputed NMD/exons table (tab/tsv) with columns including Transcript_ID, total_cds_exons, CDS_length, Protein_length, number_cds_exon, CDS_Exon_length")
    return p

def process_single_sample_worker(payload):
    p_vcf, args_dict, db_paths = payload
    try:
        dirname = os.path.dirname(p_vcf)
        filename = os.path.basename(p_vcf)
        dna_id = filename.split('_')[0]
        for ext in ['.vcf.gz', '.vcf.bgz', '.vcf', '.bgz', '.gz']:
            if dna_id.endswith(ext):
                dna_id = dna_id[:-len(ext)]
                break
        
        map_id = os.path.basename(dirname)
        s_out = os.path.join(args_dict.outdir, "individual_results", dna_id)
        siblings = glob.glob(os.path.join(dirname, "*"))
        f_vcf = next((s for s in siblings if "_father" in s and (s.endswith(".vcf") or s.endswith(".vcf.gz"))), None)
        m_vcf = next((s for s in siblings if "_mother" in s and (s.endswith(".vcf") or s.endswith(".vcf.gz"))), None)
        process_vcf(
            proband_vcf=p_vcf,
            father_vcf=f_vcf,
            mother_vcf=m_vcf,
            outdir=s_out,
            acmg_table=args_dict.acmg_table,
            db_paths=db_paths,
            exons_file=args_dict.exons_file,
            gene_bed=args_dict.gene_bed,
            proband_sample=dna_id,
            strict_qc=not args_dict.no_strict_qc,
            gnomad_v2=args_dict.gnomad_v2,
            gnomad_v3=args_dict.gnomad_v3,
            gnomad_v4=args_dict.gnomad_v4,
            acmg_recs=args_dict.acmg_recs,
            internal_db_path=args_dict.internal_db,
            family_id=map_id,
            require_mane=args_dict.require_mane,
            #cooccurrence_path=getattr(args_dict, "gnomad_cooccurrence", None),
            gnomad_api_endpoint=getattr(args_dict, "gnomad_api_endpoint", None),
            gnomad_dataset=getattr(args_dict, "gnomad_dataset", None),
            transcript_map=getattr(args_dict, "transcript_map", None),
            nmd_table=getattr(args_dict, "nmd_table", None) 
        )
        return (dna_id, map_id, s_out)
    except Exception as e:
        print(f"ERROR in worker for {p_vcf}: {e}")
        traceback.print_exc()
        return None

def run_batch_mode(args, db_paths):
    if not HAVE_PANDAS:
        logger.error("Batch mode requires pandas.")
        sys.exit(1)
    logger.info("STARTING PARALLEL BATCH ANALYSIS")
    patterns = [
        os.path.join(args.batch_input_dir, "**", "*_proband*.vcf*"),
        os.path.join(args.batch_input_dir, "**", "*_sample*.vcf*"),
        os.path.join(args.batch_input_dir, "**", "*_proband*.vcf"),
        os.path.join(args.batch_input_dir, "**", "*_sample*.vcf"),
        os.path.join(args.batch_input_dir, "**", "*.vcf*"),

    ]
    # aggregate unique VCF files
    proband_files = []
    for p in patterns:
        proband_files.extend(glob.glob(p, recursive=True))
    proband_files = sorted(list(set(proband_files)))
    if not proband_files:
        logger.error("No VCF files found (searched patterns: %s)", patterns)
        sys.exit(1)
    logger.info(f"Found {len(proband_files)} VCF files. Using ProcessPoolExecutor.")

    tasks = []
    for f in proband_files:
        try:
            samples = get_sample_names_from_vcf_header(f)
            if samples:
                sample_name = samples[0]
            else:
                # fallback to filename (basename without extension)
                sample_name = os.path.splitext(os.path.basename(f))[0]
        except Exception:
            sample_name = os.path.splitext(os.path.basename(f))[0]
        tasks.append((f, args, db_paths))

    results_meta = []
    MAX_WORKERS = 10
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(process_single_sample_worker, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            try:
                res = future.result()
                if res:
                    results_meta.append(res)
                    logger.info(f"Finished sample: {res[0]}")
            except Exception as e:
                logger.error(f"Worker process failed: {e}")
    logger.info("Aggregating results into FINAL tables...")
    auto_dfs = []
    manual_dfs = []
    all_dfs = []
    for dna_id, map_id, s_out in results_meta:
        p_auto = os.path.join(s_out, "auto_conclusions.csv")
        p_manual = os.path.join(s_out, "manual_review_list.csv")
        p_all = os.path.join(s_out, "all_candidates.csv")
        def read_csv_safe(path):
            if os.path.exists(path) and os.path.getsize(path) > 1:
                try:
                    return pd.read_csv(path, dtype=str)
                except Exception:
                    return None
            return None
        df_a = read_csv_safe(p_auto)
        if df_a is not None and not df_a.empty:
            auto_dfs.append(df_a)
        df_m = read_csv_safe(p_manual)
        if df_m is not None and not df_m.empty:
            manual_dfs.append(df_m)
        df_all = read_csv_safe(p_all)
        if df_all is not None and not df_all.empty:
            all_dfs.append(df_all)
    os.makedirs(args.outdir, exist_ok=True)
    if auto_dfs:
        final_auto = pd.concat(auto_dfs, ignore_index=True)
        final_auto.to_csv(os.path.join(args.outdir, "FINAL_auto_conclusions.csv"), index=False)
        logger.info(f"Saved FINAL_auto_conclusions.csv ({len(final_auto)} rows)")
    if manual_dfs:
        final_manual = pd.concat(manual_dfs, ignore_index=True)
        final_manual.to_csv(os.path.join(args.outdir, "FINAL_manual_review.csv"), index=False)
        logger.info(f"Saved FINAL_manual_review.csv ({len(final_manual)} rows)")
    if all_dfs:
        final_all = pd.concat(all_dfs, ignore_index=True)
        final_all.to_csv(os.path.join(args.outdir, "FINAL_all_candidates.csv"), index=False)
        logger.info(f"Saved FINAL_all_candidates.csv ({len(final_all)} rows)")
    print(f"\nBATCH FINISHED. Aggregated results in {args.outdir}")
    print(f"Total processed samples: {len(proband_files)}")

def main():
    parser = build_parser()
    args = parser.parse_args()
    db_paths = dict(DB_PATHS_DEFAULT)
    if args.db_paths_json and os.path.exists(args.db_paths_json):
        with open(args.db_paths_json) as fh:
            db_paths.update(json.load(fh))
    if args.hgmd:
        db_paths["HGMD_VCF"] = args.hgmd
    if args.batch_input_dir:
        run_batch_mode(args, db_paths)
    elif args.proband:
        process_vcf(
            proband_vcf=args.proband,
            father_vcf=args.father,
            mother_vcf=args.mother,
            outdir=args.outdir,
            acmg_table=args.acmg_table,
            db_paths=db_paths,
            exons_file=args.exons_file,
            require_mane=args.require_mane,
            gene_bed=args.gene_bed,
            proband_sample=args.proband_sample,
            strict_qc=not args.no_strict_qc,
            gnomad_v2=args.gnomad_v2,
            gnomad_v3=args.gnomad_v3,
            acmg_recs=args.acmg_recs,
            internal_db_path=args.internal_db,
            gnomad_v4=args.gnomad_v4,
            #cooccurrence_path=getattr(args, "gnomad_cooccurrence", None)
        )
    else:
        parser.print_help()

if __name__ == "__main__":
    main()