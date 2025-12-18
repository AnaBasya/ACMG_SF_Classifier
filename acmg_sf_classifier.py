#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACMG Secondary Findings Classifier

This script:
- Annotates with local VEP if needed (user must have VEP + plugins (dbNSFP, SpliceAI) configured).
- Integrates Internal Database and HGMD
- Applies gene-specific rules from ACMG SF CSVs.
- Implements PVS1, PS1/PM5, PP3/BP4, PM3 batch, PS2/PM6, TTN meta-exon handling.

Inputs:
    - Annotated VCF (VEP + dbNSFP)
    - Optional father/mother VCFs
    - ACMG Gene List & Recommendations Table
    - Genome/Exome Aggregation Databases (gnomAD v2/v3/v4)

- Outputs:
1. In single mode:
    - all_candidates.csv  (all candidates, auditing)
    - auto_conclusions.csv (automated P/LP conclusions only)
    - manual_review_list.csv (P/LP variants that require manual review)
    - run_info.json
2. In batch mode:
    - FINAL_auto_conclusions.csv
    - FINAL_manual_review.csv
    - FINAL_all_candidates.csv

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
from dataclasses import dataclass, field
from typing import Tuple, List, Dict, Any, Optional, Set
from collections import defaultdict

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
DEFAULT_VEP_CACHE = "/home/anna/anna/ACMG_SF_Classifier/databases/vep/.vep"
DB_PATHS_DEFAULT = {
    "CLINVAR_VCF": "databases/clinvar/clinvar_with_protein.vcf.gz",
}

THRESH = {
    # Computational Evidence (PP3)
    "REVEL_SUPPORTING": 0.644,
    "REVEL_MODERATE": 0.932,
    "ALPHA_MISSENSE_SUPPORTING": 0.564,  
    "CADD_MODERATE": 30,              
    
    # Splicing Predictors
    "SPLICEAI_PVS1": 0.50,              
    "SPLICEAI_MODERATE": 0.20,          
    "SPLICEAI_PS1_THRESHOLD": 0.10,     
    
    # Population Frequency
    "BA1_AF": 0.05,                     # Stand-alone benign (MAF > 5%)
    "BS1_AF": 0.01,                     # Strong benign (MAF > 1%)
    
    # PM2 Thresholds - 
    "PM2_AD_AF": 0.0001,                # 0.01% for AD 
    "PM2_AR_AF": 0.005,                 # 0.5% for AR
    "PM2_XLD_AF": 0.0001,               # 0.01% for XLD
    "PM2_XLR_AF": 0.003,                # 0.3% for XLR
    
    # BA1 Exception List
    "BA1_EXCEPTIONS": {
        "HFE": ["p.Cys282Tyr"],
        "BTD": ["p.Asp444His"]
    }
}

# -----------------------------------------------------------
# Updated Points Map based on Tavtigian et al., 2018, 2020
# -----------------------------------------------------------
POINTS_MAP = {
    # Pathogenic criteria points
    "PVS1_VeryStrong": 8,     # Loss-of-function variants
    "PVS1_Strong": 4,
    "PVS1_Moderate": 2,
    "PVS1_Supporting": 1,
    
    "PS1_Strong": 4,          # Same amino acid change as established pathogenic variant
    "PS1_Moderate": 2,
    "PS1_Supporting": 1,
    
    "PS2_VeryStrong": 8,      # De novo
    "PS2_Moderate": 2,
    
    "PM1": 2,                 # Located in a mutational hot spot and/or critical functional domain
    
    "PM2_Moderate": 2,        # Absent or at very low frequency in population databases
    
    "PM3_VeryStrong": 8,      # For recessive disorders: detected in trans with a pathogenic variant
    "PM3_Strong": 4,
    "PM3_Moderate": 2,
    "PM3_Supporting": 1,
    
    "PM4_Strong": 4,          # Protein length changes
    "PM4_Moderate": 2,
    "PM4_Supporting": 1,
    
    "PM5_Moderate": 2,        # Different amino acid change at same residue as known pathogenic variant
    "PM5_Supporting": 1,
    
    "PP1_Strong": 4,          # Co-segregation with disease
    "PP1_Moderate": 2,
    "PP1_Supporting": 1,
    
    "PP3_Moderate": 2,        # Multiple computational evidence
    "PP3_Supporting": 1,
    
    "PP5_Supporting": 1,      # Reputable source recently reports pathogenic
    
    # Benign criteria points
    "BA1": -8,                # Stand-alone benign (MAF > 5%)
    "BS1": -4,                # Strong benign (MAF > 1%)
    "BS2": -4,                # Strong benign (confirmed in healthy adult)
    "BP4": -1,                # Supporting benign (computational evidence)
    "BP5": -1,                # Supporting benign (reputable source benign)
    "BP6": -1,                # Supporting benign (segregation data)
    "BP7": -1,                # Supporting benign (synonymous with no splicing impact)
}

# -----------------------------------------------------------
# Classification thresholds (Tavtigian et al., 2020)
# -----------------------------------------------------------
CLASSIFICATION_THRESHOLDS = [
    ("Pathogenic", 10, None),           
    ("Likely pathogenic", 6, 9),        
    ("VUS", 0, 5),                      
    ("Likely benign", -5, -1),          
    ("Benign", None, -6)                
]

# -----------------------------------------------------------
# Minimum criteria requirement (except BA1)
# -----------------------------------------------------------
MIN_PATHOGENIC_CRITERIA = 2    # Minimum 2 pathogenic criteria for P/LP
MIN_BENIGN_CRITERIA = 2        # Minimum 2 benign criteria for B/LB (except BA1)


LOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
CANONICAL_SPLICE = {"splice_acceptor_variant", "splice_donor_variant"}
STRICT_QC_ENABLED = True
STRICT_MIN_DP = 15
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

    gene: str = ""
    gene_raw: str = ""
    transcript: str = ""
    consequence: str = ""
    hgvsc: str = ""
    hgvsp: str = ""
    exon: str = ""
    is_last_exon: bool = False
    nmd: Optional[str] = ""
    
    # Sex info (Default Unknown)
    sample_sex: str = "Unknown" 

    gnomad_af: Optional[float] = None
    gnomad_details: Dict[str, Any] = field(default_factory=dict)
    
    clinvar_sig: str = ""       
    clinvar_trait: str = ""     
    clinvar_review: str = ""  
    clinvar_stars: str = "0"  

    internal_db_sig: str = ""      
    internal_db_details: str = ""

    hgmd_id: str = ""               
    hgmd_class: str = ""            
    hgmd_phen: str = ""             
    hgmd_rankscore: str = "Not provided"  # Updated field name  

    pathogenic_criteria_count: int = 0          
    benign_criteria_count: int = 0              
    hgmd_publication_count: int = 0         

    gnomad_hom: Optional[int] = None
    revel: Optional[float] = None
    alpha_missense: Optional[float] = None
    spliceai: Optional[float] = None
    cadd: Optional[float] = None
    nmdesc_predictor: Optional[str] = ""

    proband_gt: str = "./."
    father_gt: str = "./."
    mother_gt: str = "./."
    proband_dp: Optional[int] = None
    proband_ad: Optional[Tuple[int,int]] = None
    proband_af: Optional[float] = None

    db_hits: List[Tuple[str,str]] = field(default_factory=list)
    criteria_assigned: List[str] = field(default_factory=list)
    criteria_points: Dict[str,int] = field(default_factory=dict)
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
    _temp_is_path: bool = False         
    _temp_is_lp: bool = False 
    
    _disease_match: bool = False
    _triage_group: str = ""

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

def normalize_chrom(chrom: str) -> str:
    c = str(chrom)
    if c.startswith("chr"):
        return c
    if re.match(r'^\d+$', c) or c in ("X","Y","MT","M"):
        return "chr" + c
    return c

def check_disease_match(acmg_disease_string: str, variant_trait: str, 
                       extra_recs_list: List[Dict] = None,
                       gene_rules: Dict[str, Any] = None) -> Tuple[bool, str, str]:
    
    # 1. Normalize inputs
    vt_clean_raw = str(variant_trait).replace("_", " ").replace("|", " ").replace(";", " ").replace(",", " ")
    acmg_clean_raw = str(acmg_disease_string).replace("_", " ").replace("|", " ").replace(";", " ").replace(",", " ")
    
    vt_lower = vt_clean_raw.lower().strip()
    acmg_lower = acmg_clean_raw.lower().strip()
    
    # 2. Smart check for missing data (REVISED LOGIC)
    # List of keywords that mean "No Data"
    no_phenotype_keywords = {
        "not provided", "not specified", "not available", "not reported", 
        "conflicting classifications", "conflicting interpretations",
        "not_provided", "not_specified", "n/a", ".", "Not provided", "not_specified"
    }
    
    # Clean the trait string by removing punctuation and extra spaces
    temp_vt = re.sub(r'\s+', ' ', re.sub(r'[^\w\s]', '', vt_lower)).strip()
    
    # Check if the cleaned string IS one of the keywords or is empty
    if not temp_vt or temp_vt in no_phenotype_keywords:
         return False, "no_db_phenotype", "No phenotype in database"

    if not acmg_lower:
        return False, "no_acmg_disease", "No ACMG disease defined"

    # 3. Check specific recommendations whitelist
    if extra_recs_list:
        for rec in extra_recs_list:
            target = rec.get("disease", "").replace("_", " ").lower()
            if target and target in vt_lower:
                if rec.get("inclusion", True):
                    return True, "whitelist_match", f"Matched recommendation: {target}"
                else:
                    return False, "whitelist_exclude", f"Excluded by recommendation: {target}"

    # 4. Strict substring & Token Match
    def clean_text(t):
        return re.sub(r'[^\w\s]', ' ', t).strip()
    
    acmg_clean = clean_text(acmg_lower)
    vt_clean = clean_text(vt_lower)
    
    # Direct substring check
    if acmg_clean in vt_clean or vt_clean in acmg_clean:
         return True, "direct_match", "Direct substring match"
         
    # Token intersection
    acmg_tokens = set(acmg_clean.split())
    vt_tokens = set(vt_clean.split())
    
    common = acmg_tokens.intersection(vt_tokens)
    
    # Filter stopwords (very common words that don't imply disease match)
    stopwords = {} 
    significant_matches = [w for w in common if len(w) > 3 and w not in stopwords]
    
    if len(significant_matches) > 0:
        return True, "token_match", f"Shared terms: {','.join(significant_matches)}"

    return False, "mismatch", "No overlap found"

def apply_gene_specific_rules(v: VariantRecord, rule: Dict[str, Any], 
                             recs_list: List[Dict], varlist: List[VariantRecord]) -> Tuple[bool, List[str]]:
    """
    Apply gene-specific reporting rules from recommendations table.
    
    Args:
        v: VariantRecord to check
        rule: Gene rules from ACMG table
        recs_list: Additional recommendations for the gene
        varlist: All variants in this gene for this sample
    
    Returns:
        Tuple of (is_reportable: bool, reasons: List[str])
    """
    reasons = []
    is_reportable = True
    
    if not recs_list:
        return is_reportable, reasons
    
    gene_norm = v.gene.upper()
    moi = str(rule.get("moi", "")).upper()
    
    # Collect all rules from recommendations
    all_rules = {}
    for rec in recs_list:
        rec_rules = rec.get("rules", {})
        all_rules.update(rec_rules)
    
    # 1. Check biallelic requirements (e.g., ATP7B, BTD require biallelic P/LP variants)
    if all_rules.get("requires_biallelic", False):
        # Count P/LP variants in this gene for this sample
        path_vars = [var for var in varlist if var.automated_class in ["Pathogenic", "Likely pathogenic"]]
        
        if len(path_vars) < 2:
            is_reportable = False
            reasons.append(f"{gene_norm} requires biallelic P/LP variants")
        
        # Note: P/LP + VUS combinations require manual assessment
        vus_vars = [var for var in varlist if var.automated_class == "VUS"]
        if len(path_vars) == 1 and len(vus_vars) >= 1:
            # Don't set is_reportable to False here - just add informational reason
            # These will be handled by AR gene logic in strict_triage_and_output
            reasons.append(f"{gene_norm}: P/LP + VUS requires manual assessment")
    
    # 2. Check hemizygous requirement for males (e.g., ABCD1 in males)
    if all_rules.get("requires_hemizygous_male", False) and v.sample_sex == "Male":
        if not _is_hom(v.proband_gt) and not (_is_het(v.proband_gt) and "X" in v.chrom):
            is_reportable = False
            reasons.append("Hemizygous variant required in males for this gene")
    
    # 3. Check if single heterozygous variants are not reportable
    if all_rules.get("single_het_not_reportable", False):
        path_vars = [var for var in varlist if var.automated_class in ["Pathogenic", "Likely pathogenic"]]
        if len(path_vars) == 1 and _is_het(v.proband_gt):
            is_reportable = False
            reasons.append("Single heterozygous variant not reportable for this gene")
    
    # 4. Check specific variant type restrictions (e.g., APOB: missense only)
    if all_rules.get("specific_variant_types") == "missense_only":
        if "missense" not in v.consequence.lower():
            is_reportable = False
            reasons.append("Only missense variants are reportable for this gene")
    
    # 5. Store reasons in variant object for final reporting
    if hasattr(v, '_rule_reasons'):
        v._rule_reasons.extend(reasons)
    else:
        v._rule_reasons = reasons
    
    return is_reportable, reasons

def _text_match_exact(target: str, query: str) -> bool:
    """
    Determine if two disease strings represent the same condition.
    Uses exact string matching with basic normalization.
    """
    def normalize_text(text: str) -> str:
        """Basic text normalization for comparison."""
        if not text:
            return ""
        
        # Convert to lowercase
        text = text.lower()
        
        # Remove common punctuation and special characters
        text = re.sub(r'[^\w\s]', ' ', text)
        
        # Replace multiple spaces with single space
        text = re.sub(r'\s+', ' ', text).strip()
        
        return text
    
    norm_target = normalize_text(target)
    norm_query = normalize_text(query)
    
    if not norm_target or not norm_query:
        return False
    
    # Direct string equality
    if norm_target == norm_query:
        return True
    
    # Check if one is substring of another
    return norm_target in norm_query or norm_query in norm_target

def _text_match_partial(target: str, query: str) -> bool:
    """
    Determine if two disease strings share common terms.
    Simple implementation without keyword extraction.
    """
    def extract_terms(text: str) -> Set[str]:
        """Extract all terms from text."""
        if not text:
            return set()
        
        text = text.lower()
        # Remove punctuation and split
        text = re.sub(r'[^\w\s]', ' ', text)
        words = text.split()
        
        # Include all words longer than 2 characters
        return {w for w in words if len(w) > 2}
    
    target_terms = extract_terms(target)
    query_terms = extract_terms(query)
    
    if not target_terms:
        return False
    
    # Check for any overlapping terms
    overlapping_terms = target_terms.intersection(query_terms)
    
    # If there's at least one overlapping term, consider it a partial match
    return len(overlapping_terms) > 0


# Robust gene normalization
_RE_LEADING_TRAILING = re.compile(r'^[\s"\'`]+|[\s"\'`]+$')
_RE_NON_ALNUM = re.compile(r'[^A-Za-z0-9_]')

def normalize_gene_name(g: str) -> str:
    """
    Normalize gene name/string: strip whitespace and quotes, remove surrounding parens,
    replace non-alphanum with underscore, collapse underscores, uppercase.
    """
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

def _is_het(gt: str) -> bool:
    return gt in ("0/1","1/0","0|1","1|0")

def _is_hom(gt: str) -> bool:
    return gt in ("1/1","1|1")

def _is_wt(gt: str) -> bool:
    return gt in ("0/0","0|0")

def clean_hgvsp(hgvsp: str) -> str:
    if not hgvsp:
        return ""
    
    if ":" in hgvsp:
        parts = hgvsp.split(":")
        if len(parts) > 1:
            return parts[-1]  # p.Ser12Phe
    
    return hgvsp

def parse_protein_pos(hgvsp: str) -> Optional[int]:
    if not hgvsp:
        return None
    hgvsp_clean = clean_hgvsp(hgvsp)
    m = re.search(r'p\.[A-Za-z]{1,3}(\d+)', hgvsp_clean)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return None
    return None

# ---------------------------
# CSQ header and parser (robust)
# ---------------------------
def parse_csq_header_from_vcf(vcf_path: str) -> List[str]:
    """
    Robust extraction of CSQ/ANN format from VCF header.
    Critically important for mapping VEP fields correctly.
    """
    logger.info("Parsing CSQ header from %s", vcf_path)
    
    # 1. Try Pysam (Best method)
    if HAVE_PYSAM:
        try:
            vf = pysam.VariantFile(vcf_path)
            # Check CSQ or ANN
            for tag in ["CSQ", "ANN"]:
                if tag in vf.header.info:
                    # Get the Description string
                    rec = vf.header.info[tag]
                    # Description is usually: "Consequence ... Format: Allele|..."
                    # We need to extract the part after "Format: "
                    desc = str(rec.description)
                    m = re.search(r'Format:\s*([^"\'>]+)', desc)
                    if m:
                        raw_fmt = m.group(1).strip().strip('"').strip('>').strip("'")
                        header = raw_fmt.split("|")
                        logger.info("Found %s header via Pysam: %d fields", tag, len(header))
                        return header
        except Exception as e:
            logger.debug("Pysam header parse failed: %s", e)

    # 2. Try Text Parsing (Fallback)
    try:
        # Detect compression
        open_func = gzip.open if vcf_path.endswith((".gz", ".bgz")) else open
        
        with open_func(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break # Header ended
                
                # Look for the INFO definition line
                if line.startswith("##INFO=<ID=CSQ") or line.startswith("##INFO=<ID=ANN"):
                    # Extract the part after "Format: " using strict regex
                    # It handles cases where the description is quoted
                    match = re.search(r'Format:\s*([A-Za-z0-9_\|\.\-\+]+)', line)
                    if match:
                        raw_fmt = match.group(1)
                        header = raw_fmt.split("|")
                        logger.info("Found CSQ header via Text: %d fields", len(header))
                        return header
    except Exception as e:
        logger.warning("Text header parse failed: %s", e)

    # 3. Last Resort Fallback (Default VEP columns)
    logger.error("CRITICAL: Could not parse CSQ format from header. Using default VEP fallback.")
    return [
        "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
        "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
        "Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE",
        "STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL",
        "MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT",
        "TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS",
        "miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
        "AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
        "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
        "MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME",
        "MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS",
        "AlphaMissense_score","CADD_phred","REVEL_score","NMD","SpliceAI_pred_DP_AG",
        "SpliceAI_pred_DP_AL","SpliceAI_pred_DP_DG","SpliceAI_pred_DP_DL",
        "SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL","SpliceAI_pred_SYMBOL"
    ]

def parse_vep_csq_line(line: str, header: List[str]) -> Dict[str,str]:
    parts = line.split("|")
    out = {}
    if header:
        for i, key in enumerate(header):
            out[key] = parts[i] if i < len(parts) else ""
    else:
        for i, p in enumerate(parts):
            out[f"H{i}"] = p
    return out

# ---------------------------
# VEP & dbNSFP Helpers 
# ---------------------------
def detect_dbnsfp_header_fields(dbnsfp_path: str) -> List[str]:
    """Reads header columns from dbNSFP (gz/bgz) to avoid requesting missing fields."""
    if not dbnsfp_path or not os.path.exists(dbnsfp_path):
        return []
    try:
        # Open as gzip text
        with gzip.open(dbnsfp_path, "rt", encoding="utf-8", errors="replace") as fh:
            for _ in range(5): # Check first few lines
                ln = fh.readline()
                if not ln: break
                ln = ln.rstrip("\n")
                # Handle headers starting with # or just raw column names
                if ln.startswith("#") or (ln and "Chr" in ln and "\t" in ln):
                    clean_ln = ln.lstrip("#").strip()
                    if "\t" in clean_ln:
                        return clean_ln.split("\t")
    except Exception as e:
        logger.debug("Failed to detect dbNSFP header: %s", e)
    return []

def choose_dbnsfp_fields(dbnsfp_path: Optional[str], requested_fields: List[str]) -> Optional[str]:
    """Returns valid plugin string for VEP: 'path,Field1,Field2' or just 'path'."""
    if not dbnsfp_path or not os.path.exists(dbnsfp_path):
        return None
    
    available_fields = detect_dbnsfp_header_fields(dbnsfp_path)
    if not available_fields:
        logger.warning("Could not read dbNSFP header. Passing path only (all fields will be added).")
        return dbnsfp_path # Fallback: let VEP try to load everything
    
    # Filter requested fields against available ones
    valid_fields = []
    missing_fields = []
    
    # Normalize comparison 
    avail_set = set(available_fields)
    
    for req in requested_fields:
        if req in avail_set:
            valid_fields.append(req)
        else:
            missing_fields.append(req)
            
    if missing_fields:
        logger.warning("The following requested dbNSFP fields are missing in the DB file and will be skipped: %s", missing_fields)
        
    if not valid_fields:
        return dbnsfp_path # Return path only if list is empty or matches nothing
        
    return dbnsfp_path + "," + ",".join(valid_fields)

def run_vep_robust(in_vcf: str, out_vcf: str, vep_cmd: str, vep_cache: str, 
                   fasta: Optional[str], extra_args: List[str], 
                   dbnsfp_plugin_str: Optional[str]) -> str:
    """Runs VEP and validates that CSQ header exists in output."""
    
    # Prepare command
    cmd_base = [vep_cmd] if shutil.which(vep_cmd) else ["conda", "run", "-n", "vep", vep_cmd]
    
    args = cmd_base + [
        "--input_file", in_vcf,
        "--output_file", out_vcf,
        "--vcf",
        "--compress_output", "bgzip",
        "--offline",
        "--cache",
        "--fork", "4",
        "--force_overwrite",
        "--everything" # Ensures CSQ is generated with rich details
    ]
    
    if vep_cache:
        args += ["--dir_cache", vep_cache]
    if fasta:
        args += ["--fasta", fasta]
    
    # Always set assembly/species if not manually overridden
    if not any("assembly" in x for x in extra_args):
        args += ["--assembly", "GRCh38"]
    if not any("species" in x for x in extra_args):
        args += ["--species", "homo_sapiens"]

    if dbnsfp_plugin_str:
        args += ["--plugin", f"dbNSFP,{dbnsfp_plugin_str}"]
    
    args += ["--plugin", "NMD"]

    sp_snv = "/home/anna/anna/ACMG_SF_Classifier/databases/vep/data/spliceai_scores.raw.snv.hg38.vcf.gz"
    sp_indel = "/home/anna/anna/ACMG_SF_Classifier/databases/vep/data/spliceai_scores.raw.indel.hg38.vcf.gz"
    
    args += ["--plugin", f"SpliceAI,snv={sp_snv},indel={sp_indel}"]

    if extra_args:
        args += extra_args

    logger.info("Running VEP command: %s", " ".join(shlex.quote(x) for x in args))
    
    # Run VEP
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error("VEP failed (exit %d). stderr:\n%s", proc.returncode, proc.stderr)
        raise RuntimeError("VEP execution failed")
        
    # Index output
    try:
        subprocess.run(["tabix", "-p", "vcf", out_vcf], check=False)
    except Exception:
        logger.warning("Could not tabix index VEP output (optional but recommended)")

    # Validate Output Header for CSQ
    has_csq = False
    try:
        # Quick check using zgrep or bcftools
        check_cmd = ["bcftools", "view", "-h", out_vcf]
        if shutil.which("bcftools"):
            res = subprocess.run(check_cmd, capture_output=True, text=True)
            if "ID=CSQ" in res.stdout:
                has_csq = True
        else:
            # Fallback python read
            with gzip.open(out_vcf, "rt") as fh:
                for _ in range(500):
                    line = fh.readline()
                    if line.startswith("##INFO=<ID=CSQ"):
                        has_csq = True; break
                    if not line.startswith("#"): break
    except Exception as e:
        logger.warning("Could not validate CSQ header presence: %s", e)
        # Assume true if check failed but VEP returncode was 0
        has_csq = True 

    if not has_csq:
        logger.error("VEP finished but 'CSQ' field is missing in VCF header. Check VEP cache version or arguments.")
        raise RuntimeError("VEP output invalid: No CSQ header")
        
    return out_vcf

# ---------------------------
# gnomAD AF aggregation utility
# ---------------------------
GNOMAD_KEY_RE = re.compile(r"(?i)gnomad|gnomad_r|gnomad\." )

def extract_gnomad_af_from_info(info: Dict[str, Any]) -> Tuple[Dict[str,Dict[str,float]], Dict[str,Any]]:
    per_version = defaultdict(dict)
    global_max = {"max_af": 0.0, "pop": None, "version": None}
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
                    try: vals.append(float(x))
                    except Exception: pass
            else:
                s = str(v)
                for sep in ("/", ",", ";", "&", "|"):
                    if sep in s:
                        parts = s.split(sep); break
                else:
                    parts = [s]
                for p in parts:
                    try: vals.append(float(p))
                    except Exception: pass
        except Exception:
            vals = []
        if not vals:
            continue
        val = max(vals)
        per_version[version][pop] = max(per_version[version].get(pop, 0.0), val)
        if val > global_max["max_af"]:
            global_max = {"max_af": val, "pop": pop, "version": version}
    return dict(per_version), global_max

# ---------------------------
# Gene validation function
# ---------------------------
def validate_gene_match_db_evidence(variant_gene: str, db_gene_info: str, db_source: str) -> Tuple[bool, str]:
    """
    Strict validation that variant gene matches database gene.
    Handles different database formats.
    
    Args:
        variant_gene: Normalized gene name from variant
        db_gene_info: Gene information from database (string format varies by DB)
        db_source: Database name ('ClinVar', 'HGMD', 'Internal_DB')
    
    Returns:
        Tuple[is_valid: bool, reason: str]
    """
    variant_gene_norm = normalize_gene_name(variant_gene)
    
    if not variant_gene_norm:
        return False, f"Variant has no gene information"
    
    # No gene info in database - cannot validate
    if not db_gene_info or db_gene_info.strip() in ["", ".", "Not provided", "not_provided", "not_specified"]:
        return True, f"No gene information in {db_source} - validation skipped"
    
    if db_source == "ClinVar":
        # ClinVar format: "GENE1:123|GENE2:456"
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
        # HGMD has direct gene column
        db_gene_norm = normalize_gene_name(db_gene_info)
        
        if not db_gene_norm:
            return True, f"No gene info in HGMD - validation skipped"
        
        if variant_gene_norm == db_gene_norm:
            return True, f"Gene match in HGMD: {variant_gene_norm}"
        else:
            return False, f"HGMD gene mismatch: {variant_gene_norm} vs {db_gene_norm}"
    
    elif db_source == "Internal_DB":
        # Internal DB doesn't have gene info in your format
        # We'll be conservative and skip validation
        return True, f"Internal DB lacks gene information - validation skipped"
    
    else:
        return True, f"Unknown database {db_source} - validation skipped"

# ---------------------------
# Database manager with TSV support
# ---------------------------
class DatabaseManager:
    """
    Manages connections to variant databases: ClinVar, HGMD, gnomAD, and Internal DB.
    Handles data retrieval, normalization, and hierarchical evidence prioritization.
    """
    def __init__(self, db_paths: Dict[str,str], 
                 gnomad_v2: Optional[str]=None, gnomad_v3: Optional[str]=None, gnomad_v4: Optional[str]=None, 
                 internal_db_path: Optional[str]=None):
        
        self.paths = dict(db_paths or {})
        self.clinvar_vcf: Optional[pysam.VariantFile] = None
        
        # HGMD is now an index dictionary, not a VCF handle
        # Structure: [chrom][pos][ref][alt] -> {tag, disease, pubs, rankscore}
        self.hgmd_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

        self.clinvar_protein_index: Dict[str, List[dict]] = defaultdict(list)
        self.internal_db_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

        self.gnomad_vcf_handles: Dict[str, pysam.VariantFile] = {}
        self.gnomad_tsv_handles: Dict[str, Dict[str,Any]] = {}
        
        self._open_clinvar()
        self._load_hgmd_csv() 
        self._open_gnomad_versions(gnomad_v2, gnomad_v3, gnomad_v4)
        
        if internal_db_path:
            self._load_internal_db(internal_db_path)


    def _open_clinvar(self):
        p = self.paths.get("CLINVAR_VCF")
        if p and os.path.exists(p) and HAVE_PYSAM:
            try:
                self.clinvar_vcf = pysam.VariantFile(p)
                logger.info("Opened ClinVar VCF: %s", p)
                try: self._build_clinvar_indexes_from_vcf(p)
                except: pass
            except Exception as e: logger.warning("Failed to open ClinVar VCF: %s", e)

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
            # try open as VCF
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
            # fallback TabixFile for TSV
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
        details = {}
        for k in ("af","allele_freq","allele_frequency","freq","frequency","aaf","adj_af","af_global","af_raw"):
            if k in rec_dict and rec_dict[k] not in ("", "."):
                try:
                    af = float(rec_dict[k])
                    details["source"] = k
                    return af, details
                except Exception:
                    pass
        for k in list(rec_dict.keys()):
            if re.search(r"(afr|amr|eas|sas|nfe|fin|asj|oth)", k):
                v = rec_dict.get(k)
                if v and v not in (".",""):
                    try:
                        af = float(v)
                        details["source"] = k
                        return af, details
                    except Exception:
                        pass
        ac = None; an = None
        for k in ("ac","allele_count","AC"):
            if k.lower() in rec_dict and rec_dict.get(k.lower()) not in ("", "."):
                try:
                    ac = float(rec_dict.get(k.lower()))
                except Exception:
                    try:
                        ac = float(rec_dict.get(k.lower()).split(",")[0])
                    except Exception:
                        ac = None
        for k in ("an","allele_number","AN"):
            if k.lower() in rec_dict and rec_dict.get(k.lower()) not in ("", "."):
                try:
                    an = float(rec_dict.get(k.lower()))
                except Exception:
                    try:
                        an = float(rec_dict.get(k.lower()).split(",")[0])
                    except Exception:
                        an = None
        if ac is not None and an:
            try:
                af = ac / an if an != 0 else None
                details["source"] = "AC/AN"
                return af, details
            except Exception:
                pass
        return None, details

    def annotate_gnomad_for_variant(self, vr: VariantRecord):
        # 1. Try VEP annotations first (gnomADg_AF, etc from CSQ)
        info = getattr(vr, "_raw_ann", {}) or {}
        perv, gmax = extract_gnomad_af_from_info(info)
        if gmax.get("max_af", 0.0) > 0:
            vr.gnomad_af = gmax["max_af"]
            vr.gnomad_details = {"version_pop": gmax, "source": "VEP_CSQ"}
            return

        # Prepare chrom keys: try both "chr1" and "1"
        chrom_clean = vr.chrom.replace("chr", "")
        chrom_vars = [chrom_clean, "chr" + chrom_clean]
        
        best = {"max_af": 0.0, "version": None, "pop": None, "orig": None}

        # 2. Query VCF handles
        for ver, vf in self.gnomad_vcf_handles.items():
            for c_query in chrom_vars:
                try:
                    # fetch needs int positions
                    for rec in vf.fetch(c_query, vr.pos - 1, vr.pos):
                        # Verify Exact Match
                        if int(rec.pos) != int(vr.pos): continue
                        if rec.ref != vr.ref: continue
                        if not (rec.alts and vr.alt in rec.alts): continue
                        
                        # Extract AF
                        info_rec = dict(rec.info)
                        perv2, gmax2 = extract_gnomad_af_from_info(info_rec)
                        if gmax2.get("max_af", 0.0) > best["max_af"]:
                            best.update({"max_af": float(gmax2["max_af"]), "version": ver, "pop": gmax2.get("pop")})
                        
                        # Check AC/AN explicitly
                        try:
                            ac = rec.info.get('AC'); an = rec.info.get('AN')
                            if ac and an:
                                ac0 = ac[0] if isinstance(ac,(list,tuple)) else ac
                                an0 = an[0] if isinstance(an,(list,tuple)) else an
                                if an0 and float(an0) > 0:
                                    afcalc = float(ac0)/float(an0)
                                    if afcalc > best["max_af"]:
                                        best.update({"max_af": afcalc, "version": ver, "pop": "global_calc"})
                        except: pass
                    # If we found something in this chrom variant, stop checking other chrom variants
                    if best["max_af"] > 0: break 
                except ValueError:
                    # Happens if contig not in header
                    continue
                except Exception:
                    continue

        # 3. Query TSV handles (Tabix)
        for ver, obj in self.gnomad_tsv_handles.items():
            tb = obj.get("tb"); colmap = obj.get("cols", {})
            if not tb: continue
            
            for c_query in chrom_vars:
                try:
                    found_any = False
                    for ln in tb.fetch(c_query, vr.pos - 1, vr.pos):
                        found_any = True
                        recd = self._parse_tsv_line_to_dict(ln, colmap)
                        recd_l = {k.lower(): v for k,v in recd.items()}
                        
                        # Verify Pos/Ref/Alt from TSV columns
                        # (Adjust key names based on what is commonly in gnomad TSVs)
                        p_col = next((k for k in recd_l if k in ('pos','position')), None)
                        r_col = next((k for k in recd_l if k in ('ref','reference')), None)
                        a_col = next((k for k in recd_l if k in ('alt','alternate')), None)

                        if p_col and int(recd_l[p_col]) != int(vr.pos): continue
                        if r_col and recd_l[r_col] != vr.ref: continue
                        if a_col and recd_l[a_col] != vr.alt: continue

                        af_val, details = self._extract_af_from_tsv_record(recd_l)
                        if af_val is not None and af_val > best["max_af"]:
                            best.update({"max_af": af_val, "version": ver, "pop": details.get("source")})
                    
                    if found_any and best["max_af"] > 0: break
                except ValueError:
                    continue
                except Exception as e:
                    continue

        if best["max_af"] > 0:
            vr.gnomad_af = best["max_af"]
            vr.gnomad_details = {"version": best["version"], "pop": best["pop"]}

    def _load_internal_db(self, path: str):
        """
        Parses the Internal DB CSV with 'vid' column (chr:pos:ref:alt).
        Maps internal codes (P, Ps, Pc -> P; LP, LPs, LPc -> LP).
        """
        if not os.path.exists(path):
            logger.warning(f"Internal DB path not found: {path}")
            return
        
        try:
            logger.info(f"Loading Internal DB from {path}...")
            
            # Read CSV (pandas handles the index column automatically)
            df = pd.read_csv(path, dtype=str).fillna("")
            
            # Columns check
            if 'vid' not in df.columns or 'anntype' not in df.columns:
                logger.warning("Internal DB CSV missing 'vid' or 'anntype' columns")
                return
            
            # Mapping Logic based on your Enum
            # E (Error) and A (Other) are ignored
            cls_map = {
                "P": "P", "Ps": "P", "Pc": "P",
                "LP": "LP", "LPs": "LP", "LPc": "LP",
                "B": "B",
                "LB": "LB",
                "VUS": "VUS"
            }
            
            count = 0
            
            for _, r in df.iterrows():
                # Parse VID: chr8:24208812:C:T
                vid = str(r.get('vid', '')).strip()
                if not vid: continue
                
                parts = vid.split(':')
                if len(parts) < 4: continue
                
                c_raw = parts[0]
                p_raw = parts[1]
                ref = parts[2]
                alt = parts[3]
                
                # Check Annotation Type
                raw_type = str(r.get('anntype', '')).strip()
                mapped_cls = cls_map.get(raw_type)
                
                # Skip if type is E, A or unknown
                if not mapped_cls:
                    continue

                # Normalize and Store
                chrom = normalize_chrom(c_raw)
                try: 
                    pos = int(p_raw)
                except ValueError: 
                    continue
                
                # Handle possible multi-allelics in alt if they exist (rare in vid format but safe to check)
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
        """
        Loads HGMD data from the new CSV format (hgmd_2025_vars_dmsupport.csv).
        Expected columns: CHROM,POS,REF,ALT,tag,disease,publications,rankscore
        """
        path = self.paths.get("HGMD_VCF") 
        
        # Fallback if argument not passed correctly
        if not path or not os.path.exists(path):
            path = "/home/nik/share/ccu-ngs/ngs/bases/hgmd/dmsupport/hgmd_2025_vars_dmsupport.csv"
            
        if not os.path.exists(path):
            logger.warning("HGMD CSV not found: %s", path)
            return

        try:
            logger.info("Loading HGMD from CSV: %s", path)
            
            try:
                df = pd.read_csv(path, dtype=str, on_bad_lines='skip').fillna("")
            except TypeError:
                # Fallback for older pandas versions
                df = pd.read_csv(path, dtype=str, error_bad_lines=False).fillna("")
            
            count = 0
            for _, r in df.iterrows():
                chrom = normalize_chrom(r.get('CHROM'))
                try:
                    pos = int(r.get('POS'))
                except: continue
                
                ref = r.get('REF', '').strip()
                alt = r.get('ALT', '').strip()
                
                if not ref or not alt: continue
                
                hgmd_gene = r.get('GENE', '') or r.get('SYMBOL', '') or ''

                # Store data
                self.hgmd_index[chrom][pos][ref][alt] = {
                    "class": r.get('tag', ''),
                    "phen": r.get('disease', ''),
                    "pubs": r.get('publications', '0'),
                    "rankscore": r.get('rankscore', ''),
                    "gene": hgmd_gene
                }
                count += 1
                
            logger.info("Loaded HGMD index with %d variants", count)
            
        except Exception as e:
            logger.warning("Failed to load HGMD CSV: %s", e)

    def get_variant_db_evidence(self, vr: VariantRecord) -> List[Tuple[str,str]]:
        """
        Aggregates pathogenicity evidence from all databases with strict quality filtering.
        Only considers HIGH QUALITY evidence: ClinGen, ClinVar ≥2★, HGMD DM with ≥2 pubs, Internal DB P/LP.
        """
        results: List[Tuple[str,str]] = []
        
        # Reset flags
        vr.has_strong_pathogenic_evidence = False
        vr.has_strong_benign_evidence = False
        vr.is_hq_pathogenic_db = False
        vr.is_hq_benign_db = False
        vr.is_hq_conflict_db = False
        
        # =========================================================
        # 1. CLINVAR 
        # =========================================================
        clinvar_evidence = []
        if self.clinvar_vcf:
            try:
                c_clean = vr.chrom.replace("chr", "")
                queries = [c_clean, "chr" + c_clean]
                found_rec = None
                
                for c_q in queries:
                    try:
                        for rec in self.clinvar_vcf.fetch(c_q, vr.pos - 1, vr.pos):
                            if getattr(rec, 'ref', None) != vr.ref: continue
                            if not (rec.alts and vr.alt in rec.alts): continue
                            
                            # GENE VALIDATION
                            gene_info = rec.info.get('GENEINFO', '')
                            clinvar_genes = []
                            if gene_info:
                                for gene_pair in str(gene_info).split('|'):
                                    if ':' in gene_pair:
                                        clinvar_gene = gene_pair.split(':')[0]
                                        clinvar_genes.append(normalize_gene_name(clinvar_gene))
                            
                            if clinvar_genes and vr.gene not in clinvar_genes:
                                vr.manual_reasons.append(f"ClinVar gene mismatch: {clinvar_genes} vs {vr.gene}")
                                continue
                            
                            found_rec = rec
                            break
                    except: continue
                    if found_rec: break
                
                if found_rec:
                    # Extract fields
                    clnsig = found_rec.info.get('CLNSIG') or ''
                    clinvar_sig_str = ";".join(clnsig) if isinstance(clnsig, (list, tuple)) else str(clnsig)
                    vr.clinvar_sig = clinvar_sig_str
                    
                    # Trait
                    clndn = found_rec.info.get('CLNDN') or ''
                    vr.clinvar_trait = ";".join(clndn) if isinstance(clndn, (list, tuple)) else str(clndn)

                    # Stars
                    rev = found_rec.info.get('CLNREVSTAT') or ''
                    rev_str = str(rev).lower()
                    stars = 0
                    if 'practice_guideline' in rev_str: stars = 4
                    elif 'expert_panel' in rev_str: stars = 3
                    elif 'criteria_provided' in rev_str and 'multiple_submitters' in rev_str and 'no_conflict' in rev_str: stars = 2
                    elif 'criteria_provided' in rev_str: stars = 1
                    
                    vr.clinvar_stars = str(stars)
                    sig_lower = clinvar_sig_str.lower()
                    
                    # HQ Logic
                    if stars >= 2:
                        if 'benign' in sig_lower and 'pathogen' not in sig_lower:
                            vr.has_strong_benign_evidence = True
                            vr.is_hq_benign_db = True
                            clinvar_evidence.append(("Benign", f"ClinVar_{stars}stars"))
                        elif 'pathogenic' in sig_lower or 'likely pathogenic' in sig_lower:
                            vr.has_strong_pathogenic_evidence = True
                            vr.is_hq_pathogenic_db = True
                            clinvar_evidence.append((clinvar_sig_str, f"ClinVar_{stars}stars"))
                        elif 'conflict' in sig_lower:
                            vr.is_hq_conflict_db = True
                            vr.conflict_details = "ClinVar Conflict"
                            clinvar_evidence.append((clinvar_sig_str, f"ClinVar_{stars}stars_conflict"))
                    else:
                        vr.manual_reasons.append(f"Low quality ClinVar: {clinvar_sig_str} ({stars} stars)")
                        
            except Exception as e:
                logger.debug(f"ClinVar extraction failed: {e}")
        
        results.extend(clinvar_evidence)
        
        # =========================================================
        # 2. INTERNAL DATABASE
        # =========================================================
        internal_evidence = []
        try:
            entry = self.internal_db_index.get(vr.chrom, {}).get(vr.pos, {}).get(vr.ref, {}).get(vr.alt)
            if entry:
                int_cls = str(entry.get("cls", "")) if isinstance(entry, dict) else str(entry)
                if int_cls in ["P", "LP"]:
                    vr.has_strong_pathogenic_evidence = True
                    vr.is_hq_pathogenic_db = True
                    classification = "Pathogenic" if int_cls == "P" else "Likely pathogenic"
                    internal_evidence.append((classification, "Internal_DB"))
                    vr.internal_db_sig = classification
                elif int_cls in ["B", "LB"]:
                    vr.has_strong_benign_evidence = True
                    vr.is_hq_benign_db = True
                    classification = "Benign" if int_cls == "B" else "Likely benign"
                    internal_evidence.append((classification, "Internal_DB"))
                    vr.internal_db_sig = classification
        except Exception as e:
            logger.debug(f"Internal DB lookup failed: {e}")
        
        results.extend(internal_evidence)
        
        # =========================================================
        # 3. HGMD 
        # =========================================================
        hgmd_evidence = []
        try:
            # Lookup in index
            hgmd_entry = self.hgmd_index.get(vr.chrom, {}).get(vr.pos, {}).get(vr.ref, {}).get(vr.alt)
            
            if hgmd_entry:
                tag = hgmd_entry.get('class', '')
                pubs_str = hgmd_entry.get('pubs', '0')
                hgmd_gene = hgmd_entry.get('gene', '')
                rank = hgmd_entry.get('rankscore', '')
                
                try:
                    pub_count = int(pubs_str)
                except (ValueError, TypeError):
                    pub_count = 0
                
                vr.hgmd_publication_count = pub_count
                vr.hgmd_rankscore = rank if rank and rank != "NULL" else "Not provided"
                vr.hgmd_phen = hgmd_entry.get('phen', '')
                vr.hgmd_class = tag

                gene_valid, gene_reason = validate_gene_match_db_evidence(
                    vr.gene, hgmd_gene, "HGMD"
                )
                
                if not gene_valid:
                    vr.manual_reasons.append(gene_reason)
                    
                    vr.hgmd_class = f"Gene mismatch: {tag}"
                    
                else:
                    # HIGH QUALITY: DM with >= 2 publications
                    if tag == 'DM' and pub_count >= 2:
                        vr.has_strong_pathogenic_evidence = True
                        vr.is_hq_pathogenic_db = True
                        hgmd_evidence.append(("Pathogenic", f"HGMD_DM_{pub_count}pubs"))
                    elif tag == 'DM':
                        # DM but <2 pubs
                        vr.has_weak_pathogenic_evidence = True
                        hgmd_evidence.append(("Pathogenic", "HGMD_DM_low_pubs"))
                    else:
                        hgmd_evidence.append((tag, "HGMD"))
                    
        except Exception as e:
            logger.debug(f"HGMD dictionary lookup error: {e}")
        
        results.extend(hgmd_evidence)
        
        # Check for conflicts between HQ sources
        hq_sources = [e for e in results if any(q in e[1] for q in ["ClinVar_2", "ClinVar_3", "ClinVar_4", "Internal_DB", "HGMD_DM_"])]
        
        if len(hq_sources) >= 2:
            classifications = set(e[0].lower() for e in hq_sources)
            # Normalize to check P vs B
            has_path = any('pathogen' in c for c in classifications)
            has_benign = any('benign' in c for c in classifications)
            
            if has_path and has_benign:
                vr.is_hq_conflict_db = True
                vr.conflict_details = f"Conflict between HQ sources: {', '.join([f'{c[0]} ({c[1]})' for c in hq_sources])}"
        
        vr.db_hits = results
        return results
            

    def check_pm5_variants(self, gene: str, pos_aa: int, current_hgvsp: str) -> Optional[str]:
        """
        Checks for other pathogenic variants at the same amino acid residue.
        Returns the classification of the best match ('Pathogenic' or 'Likely pathogenic').
        
        PM5: Different amino acid change at same residue as known pathogenic variant.
        """
        if not pos_aa or not gene:
            return None
        
        norm_gene = normalize_gene_name(gene)
        best_class = None
        
        for entry in self.clinvar_protein_index.get(norm_gene, []):
            # Skip if same exact variant (that would be PS1)
            if entry['hgvsp'] == current_hgvsp:
                continue
            
            # Clean the HGVSp
            entry_hgvsp_clean = clean_hgvsp(entry['hgvsp'])
            
            # Parse position to ensure same residue
            pos_match = re.search(r'p\.[A-Za-z]{1,3}(\d+)', entry_hgvsp_clean)
            if not pos_match:
                continue
            
            try:
                entry_pos_aa = int(pos_match.group(1))
            except:
                continue
            
            # Check if same residue
            if entry_pos_aa != pos_aa:
                continue
            
            # Check significance
            clnsig = str(entry['clnsig']).lower()
            stars = entry.get('stars', 0)
            
            # Only consider variants with at least 1 star
            if stars < 1:
                continue
            
            # Skip benign variants
            if "benign" in clnsig and "pathogen" not in clnsig:
                continue
            
            # Check for pathogenic classification
            if "pathogenic" in clnsig and "benign" not in clnsig and "likely" not in clnsig:
                # Pathogenic variant at same residue
                return "Pathogenic"
            
            if "likely pathogenic" in clnsig or ("likely" in clnsig and "pathogen" in clnsig):
                # Likely pathogenic variant at same residue
                best_class = "Likely pathogenic"
        
        return best_class

    # ClinVar indexing 
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
                    if ln.startswith("#"): continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 8: continue
                    chrom = parts[0]
                    try: pos = int(parts[1])
                    except: continue
                    info = parts[7]
                    
                    # Extract GENEINFO with multiple gene support
                    idx = info.find("GENEINFO=")
                    if idx == -1: continue
                    
                    rest = info[idx+9:]
                    end_idx = rest.find(";")
                    val = rest[:end_idx] if end_idx != -1 else rest
                    
                    if not val: continue
                    
                    # Parse all genes in GENEINFO
                    genes_raw = []
                    for gene_pair in val.split("|"):
                        if ":" in gene_pair:
                            gene_raw = gene_pair.split(":")[0]
                            gene_norm = normalize_gene_name(gene_raw)
                            if gene_norm:
                                genes_raw.append((gene_raw, gene_norm))
                    
                    if not genes_raw:
                        continue

                    # Basic fields
                    clnsig = ""
                    if "CLNSIG=" in info:
                        m = re.search(r'CLNSIG=([^;]+)', info)
                        if m: clnsig = m.group(1)
                        
                    stars = 0
                    if "CLNREVSTAT=" in info:
                        low = info.lower()
                        if 'practice_guideline' in low: stars = 4
                        elif 'expert_panel' in low: stars = 3
                        elif 'criteria_provided' in low: stars = 1
                        if stars == 1 and 'multiple_submitters' in low and 'no_conflict' in low:
                            stars = 2
                    
                    # Extract protein annotations
                    hgvsp_matches = []
                    if "HGVS_Protein=" in info:
                        hgvs_match = re.findall(r'HGVS_Protein=([^;]+)', info)
                        for hgvsp in hgvs_match:
                            if hgvsp.startswith("p."):
                                hgvsp_matches.append(hgvsp)
                    
                    if "p." in info:
                        old_matches = re.findall(r'(p\.[A-Za-z0-9\(\)\*]+)', info)
                        for mh in old_matches:
                            if mh not in hgvsp_matches:
                                hgvsp_matches.append(mh)
                    
                    # Index for each gene in GENEINFO
                    for gene_raw, gene_norm in genes_raw:
                        for hgvsp in hgvsp_matches:
                            entry = {
                                "hgvsp": hgvsp, 
                                "clnsig": clnsig, 
                                "stars": stars, 
                                "chrom": chrom, 
                                "pos": pos,
                                "gene_raw": gene_raw,
                                "gene_norm": gene_norm
                            }
                            self.clinvar_protein_index[gene_norm].append(entry)
                            protein_count += 1
                    
                    count += 1
                    
                    if count % 1000000 == 0:
                        logger.debug("Processed %d ClinVar variants, %d with protein", count, protein_count)
            
            logger.info("ClinVar loaded with gene validation. Protein variants: %d. Total variants: %d", protein_count, count)
            logger.info("Genes with protein annotations: %d", len(self.clinvar_protein_index))

        except Exception as e:
            logger.warning("ClinVar load error: %s", e)

    def find_protein_variant_by_hgvsp(self, gene: str, hgvsp: str, pos_aa: int, transcript: str) -> List[Tuple[str, str]]:
        """
        Searches for exact amino acid matches in ClinVar with GENE VALIDATION.
        For PS1 criterion: same amino acid change as established pathogenic variant.
        """
        if not hgvsp or not gene:
            return []
        
        norm_gene = normalize_gene_name(gene)
        results = []
        
        for entry in self.clinvar_protein_index.get(norm_gene, []):
            # Critical: Validate gene matches exactly
            if entry.get('gene_norm', '') != norm_gene:
                continue  # Skip if gene doesn't match exactly
            
            # Need at least 1 star for consideration
            if entry.get('stars', 0) < 1:
                continue
            
            clnsig = str(entry['clnsig']).lower()
            
            # Skip benign variants
            if "benign" in clnsig and "pathogen" not in clnsig:
                continue
            
            # Check if same protein change
            entry_hgvsp_clean = clean_hgvsp(entry['hgvsp'])
            current_hgvsp_clean = clean_hgvsp(hgvsp)
            
            if entry_hgvsp_clean == current_hgvsp_clean:
                # Exact match found
                if "pathogenic" in clnsig and "benign" not in clnsig and "likely" not in clnsig:
                    results.append((entry['clnsig'], "ClinVar_Pathogenic"))
                elif "likely pathogenic" in clnsig or ("likely" in clnsig and "pathogen" in clnsig):
                    results.append((entry['clnsig'], "ClinVar_Likely_pathogenic"))
                else:
                    results.append((entry['clnsig'], "ClinVar"))
        
        return results

    
# ---------------------------
# Exons table & NMD logic 
# ---------------------------
def load_exons_table(path: str) -> Optional[pd.DataFrame]:
    """
    Load exon structure table for rapid last exon determination.
    Supports flexible column naming in exons.csv.
    
    Expected columns (case-insensitive):
    - Transcript_ID or Transcript: Transcript identifier
    - Gene: Gene symbol
    - Exon: Exon designation
    - Genomic_start or Start: Genomic start position (hg38)
    - Genomic_end or End: Genomic end position (hg38)
    
    Returns:
        DataFrame with added '_transcript_last_exons' attribute
    """
    if not path or not os.path.exists(path):
        logger.warning(f"Exons file not found: {path}. NMD prediction disabled.")
        return None
    
    if not HAVE_PANDAS:
        logger.warning("Pandas not installed. NMD prediction disabled.")
        return None

    try:
        logger.info(f"Loading exon structure from {path}...")
        
        # Read with flexible delimiter and encoding
        df = pd.read_csv(path, dtype=str, sep=None, engine='python', 
                        encoding='utf-8', on_bad_lines='skip').fillna("")
        
        # Standardize column names (lowercase, remove spaces)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')
        
        # Map expected columns
        column_map = {}
        
        # Transcript ID
        for col in ['transcript_id', 'transcript', 'transcriptid']:
            if col in df.columns:
                column_map['transcript'] = col
                break
        
        # Gene
        for col in ['gene', 'gene_symbol', 'symbol', 'gene_name']:
            if col in df.columns:
                column_map['gene'] = col
                break
        
        # Exon
        for col in ['exon', 'exon_number', 'exon_num']:
            if col in df.columns:
                column_map['exon'] = col
                break
        
        # Start position
        for col in ['genomic_start', 'start', 'start_pos', 'position_start']:
            if col in df.columns:
                column_map['start'] = col
                break
        
        # End position
        for col in ['genomic_end', 'end', 'end_pos', 'position_end']:
            if col in df.columns:
                column_map['end'] = col
                break
        
        # Check required columns
        required = ['transcript', 'exon', 'start', 'end']
        missing = [col for col in required if col not in column_map]
        
        if missing:
            logger.error(f"Exons file missing required columns: {missing}")
            logger.error(f"Available columns: {list(df.columns)}")
            return None
        
        # Extract columns with proper names
        df_clean = pd.DataFrame()
        df_clean['Transcript_ID'] = df[column_map['transcript']]
        df_clean['Gene'] = df[column_map['gene']] if 'gene' in column_map else ''
        df_clean['Exon'] = df[column_map['exon']]
        df_clean['Start'] = pd.to_numeric(df[column_map['start']], errors='coerce')
        df_clean['End'] = pd.to_numeric(df[column_map['end']], errors='coerce')
        
        # Drop rows with invalid coordinates
        df_clean = df_clean.dropna(subset=["Start", "End"])
        df_clean = df_clean[df_clean['Start'] > 0]
        df_clean = df_clean[df_clean['End'] > df_clean['Start']]
        
        # Extract exon number from 'Exon' column
        def extract_exon_number(val):
            if pd.isna(val):
                return None
            s = str(val).lower()
            # Try different patterns
            patterns = [
                r'exon(\d+)', r'exon_(\d+)', r'ex(\d+)', r'e(\d+)',
                r'(\d+)'
            ]
            for pattern in patterns:
                match = re.search(pattern, s)
                if match:
                    try:
                        return int(match.group(1))
                    except:
                        continue
            return None
        
        df_clean['Exon_Number'] = df_clean['Exon'].apply(extract_exon_number)
        
        # For transcripts with no exon number extraction, try to infer
        missing_exons = df_clean['Exon_Number'].isna()
        if missing_exons.any():
            # Try sorting by position and assigning sequential numbers
            for tx_id in df_clean[missing_exons]['Transcript_ID'].unique():
                tx_mask = df_clean['Transcript_ID'] == tx_id
                tx_exons = df_clean[tx_mask].sort_values('Start')
                for i, (idx, _) in enumerate(tx_exons.iterrows(), 1):
                    df_clean.at[idx, 'Exon_Number'] = i
        
        # Drop transcripts with no exon numbers
        df_clean = df_clean.dropna(subset=['Exon_Number'])
        df_clean['Exon_Number'] = df_clean['Exon_Number'].astype(int)
        
        # For each transcript, identify the last exon
        transcript_last_exons = {}
        for tx_id, group in df_clean.groupby('Transcript_ID'):
            if not group['Exon_Number'].isna().all():
                last_exon_num = group['Exon_Number'].max()
                # Get the last exon entry
                last_exons = group[group['Exon_Number'] == last_exon_num]
                if not last_exons.empty:
                    last_exon = last_exons.iloc[0]
                    transcript_last_exons[tx_id] = {
                        'start': int(last_exon['Start']),
                        'end': int(last_exon['End']),
                        'exon_number': int(last_exon['Exon_Number']),
                        'gene': last_exon.get('Gene', '')
                    }
        
        # Store as DataFrame attribute
        df_clean._transcript_last_exons = transcript_last_exons
        logger.info(f"Loaded exon data for {len(transcript_last_exons)} transcripts")
        
        # Log sample of loaded data
        sample_tx = list(transcript_last_exons.keys())[:5]
        for tx in sample_tx:
            info = transcript_last_exons[tx]
            logger.debug(f"  {tx}: last exon {info['exon_number']} ({info['start']}-{info['end']})")
        
        return df_clean
        
    except Exception as e:
        logger.error(f"Failed to load exons file {path}: {e}", exc_info=True)
        return None
    

def compute_nmd_internal(vr: VariantRecord, exons_df: pd.DataFrame) -> str:
    """
    Determine NMD prediction based on last exon position.
    Uses exons.csv table for accurate last exon determination.
    
    IMPORTANT: Only last exon escapes NMD. All other exons trigger NMD.
    
    Returns:
        "Escaping_NMD_Last_Exon": Variant in last exon (escapes NMD)
        "Triggering_NMD": Variant in any other exon (triggers NMD)
        "": Not applicable or unknown
    """
    # Only applicable to LoF variants
    relevant_consequences = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"}
    if not any(t in vr.consequence.lower() for t in relevant_consequences):
        return ""  # Not applicable to non-LoF variants
    
    if exons_df is None or exons_df.empty:
        return ""  # Exon data unavailable
    
    # Get canonical transcript ID from variant annotation
    tx_id = vr.transcript
    if not tx_id:
        # Try to extract from HGVSc or raw annotation
        if vr._raw_ann and vr._raw_ann.get("Feature"):
            tx_id = vr._raw_ann.get("Feature")
        else:
            return ""  # Transcript information missing
    
    # Clean transcript ID (remove version number if present)
    if '.' in tx_id:
        tx_base = tx_id.split('.')[0]
    else:
        tx_base = tx_id
    
    # Access pre-computed last exon information
    if not hasattr(exons_df, '_transcript_last_exons'):
        logger.warning("exons_df doesn't have _transcript_last_exons attribute")
        return ""  # Last exon data not loaded properly
    
    # Try exact match first, then base match
    last_exon_info = None
    if tx_id in exons_df._transcript_last_exons:
        last_exon_info = exons_df._transcript_last_exons[tx_id]
    elif tx_base in exons_df._transcript_last_exons:
        last_exon_info = exons_df._transcript_last_exons[tx_base]
    
    if not last_exon_info:
        # Try to find by gene if transcript not found
        gene_exons = exons_df[exons_df['Gene'] == vr.gene_raw]
        if not gene_exons.empty:
            # Get all transcripts for this gene
            gene_transcripts = gene_exons['Transcript_ID'].unique()
            for gtx in gene_transcripts:
                if tx_base in gtx or tx_id in gtx:
                    last_exon_info = exons_df._transcript_last_exons.get(gtx)
                    if last_exon_info:
                        break
        
        if not last_exon_info:
            logger.debug(f"No exon data found for transcript {tx_id} (gene {vr.gene})")
            return ""
    
    # Check if variant position falls within last exon coordinates
    variant_pos = int(vr.pos)
    last_exon_start = last_exon_info['start']
    last_exon_end = last_exon_info['end']
    
    # Debug logging
    logger.debug(f"NMD check for {vr.gene}:{tx_id} at {variant_pos} (exon {last_exon_start}-{last_exon_end})")
    
    if last_exon_start <= variant_pos <= last_exon_end:
        return "Escaping_NMD_Last_Exon"
    else:
        return "Triggering_NMD"

# ---------------------------
# DbNSFP cache 
# ---------------------------
class DbnsfpSqliteCache:
    def __init__(self, path: Optional[str]):
        self.path = path
        self.conn = None
        if not path:
            return
        d = os.path.dirname(os.path.abspath(path))
        if d and not os.path.exists(d):
            os.makedirs(d, exist_ok=True)
        self.conn = sqlite3.connect(path, timeout=30)
        cur = self.conn.cursor()
        cur.execute("PRAGMA journal_mode=WAL;")
        cur.execute("CREATE TABLE IF NOT EXISTS dbnsfp_cache (key TEXT PRIMARY KEY, value TEXT, updated INTEGER)")
        self.conn.commit()

    def get(self, key: str) -> Optional[dict]:
        if not self.conn:
            return None
        cur = self.conn.cursor()
        cur.execute("SELECT value FROM dbnsfp_cache WHERE key = ?", (key,))
        row = cur.fetchone()
        if not row:
            return None
        try:
            return json.loads(row[0])
        except Exception:
            return None

    def set(self, key: str, value: dict):
        if not self.conn:
            return
        cur = self.conn.cursor()
        cur.execute("INSERT OR REPLACE INTO dbnsfp_cache (key, value, updated) VALUES (?, ?, ?)", (key, json.dumps(value), int(time.time())))
        self.conn.commit()

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

# ---------------------------
# ACMG Engine (kept as in merged file, using normalized gene names)
# ---------------------------
class ACMGEngine:
    def __init__(self, dbm: DatabaseManager, gene_rules: Dict[str, Dict[str,Any]], 
                 exons_df: Optional[Any] = None, ttn_meta: Dict[str, IntervalTree] = None, spliceai_thr: float = 0.20):
        self.dbm = dbm
        self.gene_rules = gene_rules or {}
        self.exons_df = exons_df
        self.ttn_meta = ttn_meta or {}
        self.spliceai_thr = spliceai_thr

    def _get_gene_rule(self, gene: str) -> Dict[str,Any]:
        return self.gene_rules.get((gene or "").upper(), {"reportable": False, "special_flags": {}, "raw_row": {}})
    
    def _check_upstream_pathogenic_variants(self, vr: VariantRecord, upstream_distance: int = 1000) -> bool:
        """
        Check if there are pathogenic variants upstream of start-loss variant.
        Used in PVS1 decision tree for start-loss variants.
        
        Args:
            vr: VariantRecord with start-loss consequence
            upstream_distance: Distance to search upstream (in bp)
            
        Returns:
            True if pathogenic variants found upstream, False otherwise
        """
        if not hasattr(self, 'dbm') or self.dbm is None:
            logger.debug("No DatabaseManager available for upstream variant check")
            return False
        
        if self.dbm.clinvar_vcf is None:
            logger.debug("ClinVar VCF not available for upstream variant check")
            return False
        
        try:
            # Calculate upstream search region
            upstream_start = max(0, vr.pos - upstream_distance)
            upstream_end = vr.pos - 1  # Just before the variant
            
            if upstream_start >= upstream_end:
                logger.debug(f"Invalid upstream region for {vr.chrom}:{vr.pos}")
                return False
            
            # Search in ClinVar for pathogenic variants in this gene in upstream region
            gene_norm = normalize_gene_name(vr.gene)
            pathogenic_found = False
            pathogenic_details = []
            
            # Try different chromosome formats
            chrom_variants = [vr.chrom]
            if vr.chrom.startswith("chr"):
                chrom_variants.append(vr.chrom.replace("chr", ""))
            else:
                chrom_variants.append(f"chr{vr.chrom}")
            
            for chrom_query in chrom_variants:
                try:
                    # Fetch variants from ClinVar in upstream region
                    for rec in self.dbm.clinvar_vcf.fetch(chrom_query, upstream_start, upstream_end):
                        # Check if variant is in the same gene
                        gene_info = rec.info.get('GENEINFO', '')
                        rec_genes = []
                        if gene_info:
                            for gene_pair in str(gene_info).split('|'):
                                if ':' in gene_pair:
                                    clinvar_gene = gene_pair.split(':')[0]
                                    rec_genes.append(normalize_gene_name(clinvar_gene))
                        
                        # Skip if not same gene
                        if not rec_genes or gene_norm not in rec_genes:
                            continue
                        
                        # Check if variant is protein-affecting
                        hgvsp_info = rec.info.get('HGVS_Protein', '')
                        if not hgvsp_info:
                            continue
                        
                        # Check significance
                        clnsig = rec.info.get('CLNSIG')
                        clnsig_str = ";".join(clnsig) if isinstance(clnsig, (list, tuple)) else str(clnsig)
                        sig_lower = clnsig_str.lower()
                        
                        # Check review status
                        rev = rec.info.get('CLNREVSTAT') or ''
                        rev_str = str(rev).lower()
                        stars = 0
                        if 'practice_guideline' in rev_str: stars = 4
                        elif 'expert_panel' in rev_str: stars = 3
                        elif 'criteria_provided' in rev_str and 'multiple_submitters' in rev_str and 'no_conflict' in rev_str: stars = 2
                        elif 'criteria_provided' in rev_str: stars = 1
                        
                        # Consider variants with ≥1 star and pathogenic classification
                        if stars >= 1:
                            if ('pathogenic' in sig_lower or 'likely pathogenic' in sig_lower) and 'benign' not in sig_lower:
                                pathogenic_found = True
                                variant_info = {
                                    'position': rec.pos,
                                    'significance': clnsig_str,
                                    'stars': stars,
                                    'hgvsp': hgvsp_info[0] if isinstance(hgvsp_info, (list, tuple)) else hgvsp_info
                                }
                                pathogenic_details.append(variant_info)
                                logger.debug(f"Found upstream pathogenic variant at {rec.chrom}:{rec.pos}: {clnsig_str} ({stars} stars)")
                except ValueError:
                    # Chromosome not in ClinVar header
                    continue
                except Exception as e:
                    logger.debug(f"Error fetching upstream variants for {chrom_query}: {e}")
                    continue
            
            if pathogenic_found:
                vr.manual_reasons.append(f"Found {len(pathogenic_details)} pathogenic variant(s) upstream:")
                for i, var in enumerate(pathogenic_details, 1):
                    vr.manual_reasons.append(f"  {i}. {var['position']}: {var['significance']} ({var['stars']} stars) - {var['hgvsp']}")
            
            return pathogenic_found
            
        except Exception as e:
            logger.error(f"Error checking upstream pathogenic variants: {e}")
            return False  # Conservative: assume no pathogenic variants if error


    def _find_nearest_inframe_start_codon(self, vr: VariantRecord, fasta_handle=None) -> Tuple[Optional[int], Optional[str]]:
        """
        Find the nearest in-frame start codon (ATG) downstream of start-loss variant.
        
        Args:
            vr: VariantRecord with start-loss consequence
            fasta_handle: pysam.FastaFile handle for reference genome
            
        Returns:
            Tuple of (position, sequence) of nearest in-frame ATG, or (None, None) if not found
        """
        if fasta_handle is None or not HAVE_PYSAM:
            logger.debug("No FASTA handle available for start codon search")
            return None, None
        
        try:
            # Search for ATG in reading frame downstream
            # We'll search up to 300bp downstream (100 codons)
            search_distance = 300
            search_end = vr.pos + search_distance
            
            # Fetch sequence
            sequence = fasta_handle.fetch(vr.chrom, vr.pos, search_end).upper()
            
            if len(sequence) < 3:
                return None, None
            
            # Check all three reading frames
            for frame in range(3):
                for i in range(frame, len(sequence) - 2, 3):
                    codon = sequence[i:i+3]
                    if codon == "ATG":
                        # Found ATG in this reading frame
                        atg_position = vr.pos + i
                        # Check if it's in-frame with the original start
                        # For start-loss variants, we want an ATG that maintains the reading frame
                        position_diff = atg_position - vr.pos
                        if position_diff % 3 == 0:  # In-frame
                            # Get more context
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
    
    def _check_protein_length_for_nmd(self, vr: VariantRecord) -> bool:
        """
        Check if nonsense transcript is at least 90% of wild-type length.
        "Не меньше ли транскрипт с нонсенсом, чем у дикого типа, на 90%?"
        
        Returns True if transcript is ≥90% of wild-type length, False otherwise.
        """
        if not vr.hgvsp:
            return True  # Conservative: assume OK if no info
        
        # Parse protein position from HGVSp
        pos_match = re.search(r'p\.[A-Za-z*]{1,3}(\d+)', vr.hgvsp)
        if not pos_match:
            # Try other patterns
            pos_match = re.search(r'p\.(\d+)', vr.hgvsp)
            if not pos_match:
                return True  # Conservative
        
        try:
            trunc_pos = int(pos_match.group(1))
            
            # Known protein lengths for common ACMG genes
            protein_lengths = {
                "BRCA1": 1863, "BRCA2": 3418, "TP53": 393, "MLH1": 756,
                "MSH2": 934, "MSH6": 1360, "PMS2": 862, "APC": 2843,
                "MYH7": 1935, "MYBPC3": 1274, "TNNT2": 298, "TNNI3": 210,
                "KCNH2": 1159, "KCNQ1": 676, "SCN5A": 2016, "RYR2": 4967,
                "DSP": 2871, "PKP2": 881, "DSG2": 1118, "DSC2": 901,
                "LMNA": 664, "MYH11": 1936, "ACTA2": 377, "SMAD3": 425,
                "COL3A1": 1466, "FBN1": 2871, "TGFBR1": 503, "TGFBR2": 592,
                "CACNA1S": 1873, "RYR1": 5037, "BTD": 543, "HFE": 343,
                "GLA": 429, "GAA": 952, "PMM2": 246, "ALDOB": 364,
                "CFTR": 1480, "SERPINA1": 418, "ATP7B": 1465, "PAH": 452,
                "RB1": 928, "VHL": 213, "MEN1": 610, "RET": 1114,
                "SDHD": 159, "SDHAF2": 166, "SDHC": 169, "SDHB": 280
            }
            
            # Get protein length for this gene
            gene_key = vr.gene.upper()
            if gene_key in protein_lengths:
                protein_len = protein_lengths[gene_key]
                
                # Calculate percentage of protein before truncation
                percentage = (trunc_pos / protein_len) * 100
                
                # Check if ≥90%
                return percentage >= 90
            else:
                # Gene not in known lengths - use heuristic
                # If truncation is after position 100, assume it's ≥90%
                return trunc_pos > 100
                
        except Exception as e:
            logger.debug(f"Error checking protein length for NMD: {e}")
            return True  # Conservative
        
    def _check_last_exon_status(self, vr: VariantRecord) -> bool:
        """
        Check if variant is in the last exon using exons.csv table OR VEP annotation.
        """
        # 1. Попытка использовать точные координаты из exons.csv 
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
                        # Если нашли транскрипт в CSV, и вариант не попадает в координаты последнего,
                        # значит он точно не в последнем. Возвращаем False.
                        return False 

        # 2. Fallback: Если CSV недоступен или транскрипт там не найден, используем строку VEP "Exon A/B"
        if not used_csv and vr.exon:
            try:
                # Ожидаем формат "16/16" или "10/23"
                if "/" in vr.exon:
                    parts = vr.exon.split("/")
                    curr = int(parts[0])
                    total = int(parts[1])
                    if curr == total and total > 0:
                        return True
            except Exception:
                pass
        
        return False
    
    def _ttn_overlap(self, chrom: str, pos: int) -> Optional[dict]:
        chromn = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
        if chromn not in self.ttn_meta:
            return None
        hits = self.ttn_meta[chromn].at(pos - 1)
        if not hits:
            return None
        return max(hits, key=lambda it: float(it.data.get("psi_meta", 0) or 0)).data

    def evaluate_variant(self, vr: VariantRecord, fasta_handle=None):
        pts = {}
        assigned = []
        auto_class = "VUS"
        manual_review = False
        manual_reasons = []
        
        # ---------------------------------------------------------
        # 0. BA1: Stand-alone Benign
        # ---------------------------------------------------------
        af = vr.gnomad_af or 0.0
        is_ba1_exception = False
        exceptions = THRESH.get("BA1_EXCEPTIONS", {})
        if vr.gene in exceptions:
            for exc in exceptions[vr.gene]:
                if exc in vr.hgvsp: is_ba1_exception = True; break
        
        if af >= THRESH["BA1_AF"] and not is_ba1_exception:
             assigned.append("BA1"); pts["BA1"] = -8
             return assigned, pts, "Benign", False, ["BA1 (AF > 5%)"]
        
        # 1. BS1
        if af >= THRESH["BS1_AF"]: assigned.append("BS1"); pts["BS1"] = -4

        # Gene Rules
        rule = self.gene_rules.get(vr.gene, {})
        cons = vr.consequence.lower()
        
        # ---------------------------------------------------------
        # 2. PVS1 Decision Tree (Includes Canonical Splice)
        # ---------------------------------------------------------
        is_lof_molecular = any(x in cons for x in ["stop_gained", "frameshift", "start_lost"])
        is_canonical_splice = any(x in cons for x in ["splice_acceptor", "splice_donor"])
        
        if is_lof_molecular or is_canonical_splice:
            special = rule.get("special_flags", {})
            if not special.get("lof_not_reportable", False):
                # NMD Check (Decision Tree Node)
                # "Variant in last exon?"
                in_last_exon = self._check_last_exon_status(vr)
                
                # "Predicted NMD?" (If not last exon -> Yes)
                triggers_nmd = not in_last_exon 
                
                if "start_lost" in cons:
                    assigned.append("PVS1_Moderate"); pts["PVS1_Moderate"] = 2
                
                elif triggers_nmd:
                    # Triggers NMD -> PVS1 VeryStrong
                    assigned.append("PVS1_VeryStrong"); pts["PVS1_VeryStrong"] = 8
                
                else:
                    # Escapes NMD (Last exon or close to end)
                    if is_canonical_splice:
                        # Splice in last exon -> likely PVS1_Strong
                        assigned.append("PVS1_Strong"); pts["PVS1_Strong"] = 4
                    else:
                        # Stop/Frameshift in last exon -> check protein length
                        # "Removes >10% of protein?"
                        is_long_trunc = self._check_protein_truncation_percentage(vr)
                        if is_long_trunc:
                            assigned.append("PVS1_Strong"); pts["PVS1_Strong"] = 4
                        else:
                            assigned.append("PVS1_Moderate"); pts["PVS1_Moderate"] = 2

        # ---------------------------------------------------------
        # 3. PM4 (Protein Length)
        # ---------------------------------------------------------
        if "stop_lost" in cons:
            assigned.append("PM4_Strong"); pts["PM4_Strong"] = 4
        
        if "inframe" in cons:
            if self._check_inframe_small_noncritical(vr):
                assigned.append("PM4_Supporting"); pts["PM4_Supporting"] = 1
            else:
                assigned.append("PM4_Moderate"); pts["PM4_Moderate"] = 2

        # ---------------------------------------------------------
        # 4. PS1 & PM5 (Missense Analysis)
        # ---------------------------------------------------------
        pos_aa = parse_protein_pos(vr.hgvsp)
        if "missense" in cons and pos_aa:
            # PS1
            matches = self.dbm.find_protein_variant_by_hgvsp(vr.gene, vr.hgvsp, pos_aa, vr.transcript)
            ps1_assigned = False
            for m_cls, src in matches:
                # Strong check: must be P/LP and not conflicting
                if ("pathogenic" in m_cls.lower() and "benign" not in m_cls.lower() and "conflicting" not in m_cls.lower()):
                    assigned.append("PS1_Strong"); pts["PS1_Strong"] = 4
                    ps1_assigned = True
                    break
            
            # PM5 (If no PS1)
            if not ps1_assigned:
                pm5_match = self.dbm.check_pm5_variants(vr.gene, pos_aa, vr.hgvsp)
                if pm5_match:
                    if pm5_match == "Pathogenic":
                        assigned.append("PM5_Moderate"); pts["PM5_Moderate"] = 2
                    else:
                        assigned.append("PM5_Supporting"); pts["PM5_Supporting"] = 1

        # ---------------------------------------------------------
        # 5. PM2 (Pop Freq)
        # ---------------------------------------------------------
        moi = str(rule.get("moi","")).upper()
        pm2_thr = THRESH["PM2_AD_AF"]
        if "AR" in moi or "RECESSIVE" in moi: pm2_thr = THRESH["PM2_AR_AF"]
        elif "XLR" in moi: pm2_thr = THRESH["PM2_XLR_AF"]
        elif "XLD" in moi: pm2_thr = THRESH["PM2_XLD_AF"]
        
        if af < pm2_thr:
             assigned.append("PM2_Moderate"); pts["PM2_Moderate"] = 2

        # ---------------------------------------------------------
        # 6. PP3 (Computational)
        # ---------------------------------------------------------
        s_score = vr.spliceai if vr.spliceai is not None else 0.0
        r_score = vr.revel if vr.revel is not None else -1.0
        a_score = vr.alpha_missense if vr.alpha_missense is not None else -1.0
        c_score = vr.cadd if vr.cadd is not None else 0.0
        
        is_pp3_moderate = False
        if s_score >= THRESH["SPLICEAI_MODERATE"]: is_pp3_moderate = True
        elif r_score >= THRESH["REVEL_MODERATE"]: is_pp3_moderate = True
        elif c_score >= THRESH["CADD_MODERATE"]: is_pp3_moderate = True
        
        if is_pp3_moderate:
            assigned.append("PP3_Moderate"); pts["PP3_Moderate"] = 2
        elif r_score >= THRESH["REVEL_SUPPORTING"]:
            assigned.append("PP3_Supporting"); pts["PP3_Supporting"] = 1
        elif r_score >= THRESH["ALPHA_MISSENSE_SUPPORTING"]:
            assigned.append("PP3_Supporting"); pts["PP3_Supporting"] = 1
        

        # ---------------------------------------------------------
        # 7. PP5 (Database)
        # ---------------------------------------------------------
        # Calculate HQ locally to ensure availability
        is_hq = False
        if vr.clinvar_stars and vr.clinvar_stars.isdigit() and int(vr.clinvar_stars) >= 2:
             if "pathogenic" in vr.clinvar_sig.lower() and "benign" not in vr.clinvar_sig.lower():
                 is_hq = True
        if vr.hgmd_class == "DM" and vr.hgmd_publication_count >= 2:
             is_hq = True
        if vr.internal_db_sig and "pathogenic" in vr.internal_db_sig.lower():
             is_hq = True

        if is_hq and not getattr(vr, 'is_hq_conflict_db', False):
             assigned.append("PP5_Supporting"); pts["PP5_Supporting"] = 1

        # ---------------------------------------------------------
        # Scoring
        # ---------------------------------------------------------
        total = sum(pts.values())
        path_cnt = sum(1 for x in assigned if any(p in x for p in ["PVS","PS","PM","PP"]))
        
        if total >= 10: auto_class = "Pathogenic"
        elif total >= 6: auto_class = "Likely pathogenic"
        elif total <= -6: auto_class = "Benign"
        elif total <= -1: auto_class = "Likely benign"
        else: auto_class = "VUS"
        
        if (auto_class.startswith("Path") or auto_class.startswith("Likely path")) and path_cnt < 2:
            auto_class = "VUS"
            vr.manual_reasons.append("Downgraded: <2 pathogenic criteria")

        return assigned, pts, auto_class, False, manual_reasons

    def _check_canonical_splice_position(self, vr: VariantRecord) -> bool:
        """
        Check if variant is in canonical ±1,2 splice positions.
        Simplified implementation - in real pipeline would parse exact position.
        """
        return "splice_acceptor" in vr.consequence or "splice_donor" in vr.consequence

    def _check_protein_truncation_percentage(self, vr: VariantRecord) -> bool:
        """
        Check if truncation removes more than 10% of protein.
        Uses protein position from HGVSp annotation and known protein lengths.
        
        Returns:
            True if truncation removes >10% of protein
            False if truncation removes ≤10% of protein or unknown
        """
        if not vr.hgvsp:
            return True  # Conservative: assume >10% if no info
        
        # Parse protein position from HGVSp
        pos_match = re.search(r'p\.[A-Za-z*]{1,3}(\d+)', vr.hgvsp)
        if not pos_match:
            # Try other patterns
            pos_match = re.search(r'p\.(\d+)', vr.hgvsp)
            if not pos_match:
                return True  # Conservative
        
        try:
            trunc_pos = int(pos_match.group(1))
            
            # Known protein lengths for common ACMG genes (extend as needed)
            protein_lengths = {
                "BRCA1": 1863, "BRCA2": 3418, "TP53": 393, "MLH1": 756,
                "MSH2": 934, "MSH6": 1360, "PMS2": 862, "APC": 2843,
                "MYH7": 1935, "MYBPC3": 1274, "TNNT2": 298, "TNNI3": 210,
                "KCNH2": 1159, "KCNQ1": 676, "SCN5A": 2016, "RYR2": 4967,
                "DSP": 2871, "PKP2": 881, "DSG2": 1118, "DSC2": 901,
                "LMNA": 664, "MYH11": 1936, "ACTA2": 377, "SMAD3": 425,
                "COL3A1": 1466, "FBN1": 2871, "TGFBR1": 503, "TGFBR2": 592,
                "CACNA1S": 1873, "RYR1": 5037, "BTD": 543, "HFE": 343,
                "GLA": 429, "GAA": 952, "PMM2": 246, "ALDOB": 364,
                "CFTR": 1480, "SERPINA1": 418, "ATP7B": 1465, "PAH": 452,
                "RB1": 928, "VHL": 213, "MEN1": 610, "RET": 1114,
                "SDHD": 159, "SDHAF2": 166, "SDHC": 169, "SDHB": 280
            }
            
            # Get protein length for this gene
            gene_key = vr.gene.upper()
            if gene_key in protein_lengths:
                protein_len = protein_lengths[gene_key]
                
                # Calculate truncation percentage
                trunc_percentage = (trunc_pos / protein_len) * 100
                
                # Debug logging
                logger.debug(f"Truncation check for {vr.gene}: pos={trunc_pos}, len={protein_len}, %={trunc_percentage:.1f}%")
                
                return trunc_percentage > 10  # >10% returns True
            else:
                # Gene not in known lengths - use heuristic
                if trunc_pos < 100:
                    return True  # Early truncation likely >10%
                else:
                    return False  # Later truncation likely ≤10%
                    
        except Exception as e:
            logger.debug(f"Error calculating truncation percentage: {e}")
            return True  # Conservative

    def _check_inframe_small_noncritical(self, vr: VariantRecord) -> bool:
        """
        Check if in-frame variant is small (<10 AA).
        """
        hgvsp_clean = clean_hgvsp(vr.hgvsp)
        if not hgvsp_clean:
            return False
        
        # Determine variant size
        size = 0
        
        # Check for deletion size
        del_match = re.search(r'del(\d+)[A-Za-z]*$', hgvsp_clean)
        if del_match:
            try:
                size = int(del_match.group(1))
            except:
                pass
        
        # Check for insertion size
        ins_match = re.search(r'ins(\d+)[A-Za-z]*$', hgvsp_clean)
        if ins_match:
            try:
                size = int(ins_match.group(1))
            except:
                pass
        
        # Check for single AA deletion
        single_del_match = re.search(r'p\.[A-Za-z]{1,3}\d+del$', hgvsp_clean)
        if single_del_match and size == 0:
            size = 1
        
        # Check for duplication size
        dup_match = re.search(r'dup(\d+)[A-Za-z]*$', hgvsp_clean)
        if dup_match:
            try:
                size = int(dup_match.group(1))
            except:
                pass
        
        # Check if small (<10 AA)
        if size == 0:
            return False
        
        return size < 10

    def _compare_revel_scores(self, variant_revel: Optional[float], reference_revel: Optional[float]) -> str:
        """
        Compare REVEL scores for PM5 evaluation.
        
        Returns:
            "greater_or_equal": variant REVEL ≥ reference REVEL
            "less": variant REVEL < reference REVEL
            "insufficient_data": cannot compare (one or both missing)
            "equal_position_different_aa": same position but different AA change
        """
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
        Apply PM3 criteria for recessive (AR) and X-linked recessive (XLR) genes.
        FOLLOWS EXACT FLOWCHART SCORING SYSTEM with strict limitations for unphased variants.
        
        Scoring system from flowchart:
        - Confirmed in trans with P/LP variant: +1.0 point
        - Phase unknown with P variant: +0.5 point
        - Phase unknown with LP variant: +0.5 point
        - Homozygous case: +0.5 point (max 1.0 total for unphased)
        - With VUS variant in trans: +0.25 point
        
        Total score → PM3 level:
        - 0.5 → PM3_Supporting
        - 1.0 → PM3_Moderate
        - 2.0 → PM3_Strong
        - ≥4.0 → PM3_VeryStrong
        
        CRITICAL LIMITATION: Without parental phasing data, max PM3_Supporting.
        """
        # Group variants by sample and gene
        by_sample_gene = defaultdict(lambda: defaultdict(list))
        for v in candidates:
            by_sample_gene[v.sample][v.gene].append(v)
        
        for sample, genes in by_sample_gene.items():
            for gene, varlist in genes.items():
                rule = self._get_gene_rule(gene)
                moi = str(rule.get("moi", "")).upper()
                
                # Only apply PM3 to AR or XLR genes
                if not ("AR" in moi or "RECESSIVE" in moi or "XLR" in moi):
                    continue
                
                # Check if we have any parental data for this gene
                has_parental_data = any(
                    any(gt != "./." for gt in [v.father_gt, v.mother_gt])
                    for v in varlist
                )
                
                # Separate variants by classification
                path_vars = [v for v in varlist if v.automated_class in ["Pathogenic", "Likely pathogenic"]]
                vus_vars = [v for v in varlist if v.automated_class == "VUS"]
                
                # Process each variant
                for v in varlist:
                    pm3_score = 0.0
                    unphased_contributions = 0.0
                    
                    # Track which criteria contributed
                    criteria_details = []
                    
                    # 1. Homozygous case: +0.5 point (max 1.0 for unphased total)
                    if _is_hom(v.proband_gt):
                        pm3_score += 0.5
                        unphased_contributions += 0.5
                        criteria_details.append("homozygous: +0.5")
                    
                    # 2. Check compound heterozygous with other variants
                    for other_v in varlist:
                        if other_v == v:
                            continue
                        
                        # Check phasing relationship
                        in_trans_confirmed = False
                        
                        # Scenario 1: v from father, other from mother
                        if (_is_het(v.father_gt) and _is_wt(v.mother_gt) and
                            _is_het(other_v.mother_gt) and _is_wt(other_v.father_gt)):
                            in_trans_confirmed = True
                        
                        # Scenario 2: v from mother, other from father
                        elif (_is_het(v.mother_gt) and _is_wt(v.father_gt) and
                            _is_het(other_v.father_gt) and _is_wt(other_v.mother_gt)):
                            in_trans_confirmed = True
                        
                        # Apply scoring based on flowchart
                        if in_trans_confirmed:
                            # Confirmed in trans
                            if other_v.automated_class in ["Pathogenic", "Likely pathogenic"]:
                                pm3_score += 1.0
                                criteria_details.append(f"confirmed in trans with P/LP ({other_v.gene}): +1.0")
                            elif other_v.automated_class == "VUS":
                                pm3_score += 0.25
                                criteria_details.append(f"confirmed in trans with VUS ({other_v.gene}): +0.25")
                        
                        elif has_parental_data and not in_trans_confirmed:
                            # Parental data exists but not in trans → no points
                            pass
                        
                        else:
                            # No parental data (phase unknown)
                            contribution = 0.0
                            if other_v.automated_class == "Pathogenic":
                                contribution = 0.5
                                criteria_details.append(f"unphased with P ({other_v.gene}): +0.5")
                            elif other_v.automated_class == "Likely pathogenic":
                                contribution = 0.5
                                criteria_details.append(f"unphased with LP ({other_v.gene}): +0.5")
                            
                            # Limit unphased contributions to max 1.0
                            if unphased_contributions + contribution <= 1.0:
                                unphased_contributions += contribution
                                pm3_score += contribution
                            elif unphased_contributions < 1.0:
                                remaining = 1.0 - unphased_contributions
                                unphased_contributions += remaining
                                pm3_score += remaining
                                criteria_details.append(f"unphased (capped): +{remaining:.2f}")
                    
                    # Apply constraint: unphased contributions cannot exceed 1.0
                    if unphased_contributions > 1.0:
                        excess = unphased_contributions - 1.0
                        pm3_score -= excess
                        unphased_contributions = 1.0
                    
                    # Map score to PM3 evidence level
                    pm3_level = None
                    if pm3_score >= 4.0:
                        pm3_level = "PM3_VeryStrong"
                    elif pm3_score >= 2.0:
                        pm3_level = "PM3_Strong"
                    elif pm3_score >= 1.0:
                        pm3_level = "PM3_Moderate"
                    elif pm3_score >= 0.5:
                        pm3_level = "PM3_Supporting"
                    
                    # CRITICAL LIMITATION: Without parental data, max is PM3_Supporting
                    if pm3_level and not has_parental_data:
                        if pm3_level in ["PM3_Moderate", "PM3_Strong", "PM3_VeryStrong"]:
                            old_level = pm3_level
                            pm3_level = "PM3_Supporting"
                            v.manual_reasons.append(f"PM3 downgraded from {old_level} to {pm3_level} (no parental phasing)")
                    
                    # Apply PM3 criteria if applicable
                    if pm3_level and pm3_level not in v.criteria_assigned:
                        v.criteria_assigned.append(pm3_level)
                        v.criteria_points[pm3_level] = POINTS_MAP.get(pm3_level, 1)
                        
                        # Add detailed explanation
                        if criteria_details:
                            details_str = ", ".join(criteria_details)
                            v.manual_reasons.append(f"PM3 ({pm3_level}): total={pm3_score:.2f} [{details_str}]")
                        else:
                            v.manual_reasons.append(f"PM3 ({pm3_level}): total={pm3_score:.2f}")

# ---------------------------
# Deduplication helper to reduce duplicate rows in outputs
# ---------------------------
def collapse_output_rows(rows: List[Dict[str,Any]]) -> List[Dict[str,Any]]:
    """
    Collapse rows by unique key (sample,chrom,pos,ref,alt,normalized_gene).
    Merge manual_reasons and criteria_assigned/points if duplicates found.
    """
    out = {}
    for r in rows:
        norm_gene = normalize_gene_name(r.get("gene") or r.get("normalized_gene","") or "")
        key = f"{r.get('sample','')}-{r.get('chrom','')}-{r.get('pos','')}-{r.get('ref','')}-{r.get('alt','')}-{norm_gene}"
        if key not in out:
            nr = dict(r)
            # ensure criteria_points stored as dict
            try:
                if isinstance(nr.get("criteria_points"), str) and nr.get("criteria_points"):
                    nr["criteria_points"] = json.loads(nr["criteria_points"])
            except Exception:
                pass
            out[key] = nr
            continue
        existing = out[key]
        # merge manual_reasons
        mr_existing = set([x for x in str(existing.get("manual_reasons","")).split(";") if x])
        mr_new = set([x for x in str(r.get("manual_reasons","")).split(";") if x])
        merged_mr = ";".join(sorted(mr_existing.union(mr_new)))
        existing["manual_reasons"] = merged_mr
        # merge criteria_assigned
        ca_existing = set([x for x in str(existing.get("criteria_assigned","")).split(";") if x])
        ca_new = set([x for x in str(r.get("criteria_assigned","")).split(";") if x])
        merged_ca = ";".join(sorted(ca_existing.union(ca_new)))
        existing["criteria_assigned"] = merged_ca
        # merge criteria_points (dict merge, keep max)
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
        for k,v in cp_new.items():
            try:
                vi = int(v)
            except Exception:
                try: vi = int(float(v))
                except Exception: vi = 0
            if k not in cp_existing:
                cp_existing[k] = vi
            else:
                try:
                    cp_existing[k] = max(int(cp_existing.get(k,0)), vi)
                except Exception:
                    cp_existing[k] = vi
        existing["criteria_points"] = cp_existing
        # recalc total_points as sum
        try:
            existing["total_points"] = int(sum(int(x) for x in existing["criteria_points"].values()))
        except Exception:
            pass
        # prefer more severe automated_class
        order = {"Pathogenic":4, "Likely pathogenic":3, "VUS":2, "Likely benign":1, "Benign":0}
        a1 = existing.get("automated_class","VUS"); a2 = r.get("automated_class","VUS")
        if order.get(a2,2) > order.get(a1,2):
            existing["automated_class"] = a2
    # convert criteria_points back to JSON strings for stable CSV output
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
# Output helpers (strict_triage_and_output modified to call collapse_output_rows)
# ---------------------------

def variant_row_with_disease(v: VariantRecord, rule: Dict[str,Any], family_id: str = "", recs_map: Dict = None) -> dict:
    """
    Format output row with updated columns: Special Guidance, Enhanced DB Summary, Sample AF.
    """
    acmg_disease = rule.get("disease", "") if isinstance(rule, dict) else ""
    
    # Check if gene has special guidance
    has_guidance = "No"
    if recs_map and v.gene in recs_map:
        if recs_map[v.gene]: has_guidance = "Yes"

    # Combine reasons
    reasons = "; ".join(v.manual_reasons) if v.manual_reasons else ""
    
    # Format DB phenotypes
    clinvar_phen = v.clinvar_trait
    if not clinvar_phen or clinvar_phen.lower() in [".", "not specified", "not provided", ""]:
        clinvar_phen = "Not provided"
        
    hgmd_phen = v.hgmd_phen
    if not hgmd_phen or hgmd_phen.lower() in [".", "not specified", "not provided", ""]:
        hgmd_phen = "Not provided"

    # Format Scores
    revel_val = f"{v.revel:.3f}" if v.revel is not None else "Not provided"
    cadd_val = f"{v.cadd:.2f}" if v.cadd is not None else "Not provided"
    splice_val = f"{v.spliceai:.3f}" if v.spliceai is not None else "Not provided"
    alpha_val = f"{v.alpha_missense:.3f}" if v.alpha_missense is not None else "Not provided"

    # Zygosity & Hemizygous Check
    zyg = "Het"
    if _is_hom(v.proband_gt): 
        zyg = "Hom"
        if v.sample_sex == "Male" and ("X" in v.chrom or "Y" in v.chrom):
            zyg = "Hemi"
            
    # Database Summary
    db_sum_parts = []
    if v.clinvar_sig: 
        stars = v.clinvar_stars if v.clinvar_stars else "0"
        db_sum_parts.append(f"ClinVar: {v.clinvar_sig} ({stars}*)")
    if v.hgmd_class: 
        pubs = v.hgmd_publication_count if v.hgmd_publication_count else 0
        db_sum_parts.append(f"HGMD: {v.hgmd_class} ({pubs} pubs)")
    if v.internal_db_sig: 
        db_sum_parts.append(f"Internal: {v.internal_db_sig}")
        
    db_summary = "; ".join(db_sum_parts) if db_sum_parts else "No database evidence"
    gnomad_ver = v.gnomad_details.get("version", "Unknown") if v.gnomad_details else "Unknown"

    # Flags
    in_auto = getattr(v, '_in_auto', False)
    in_manual = getattr(v, '_in_manual', False)
    filtered = getattr(v, '_filtered', True)
    
    # Sample AF (Allele Fraction)
    sample_af_val = f"{v.proband_af:.3f}" if v.proband_af is not None else "0.0"

    # Get inheritance from gene rules
    moi = rule.get("moi", "Unknown") if isinstance(rule, dict) else "Unknown"

    output_row = {
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
        "Exon": v.exon,
        "Consequence": v.consequence,
        "Zygosity": zyg,
        "DP": v.proband_dp if v.proband_dp is not None else "0",
        "Sample_AF": sample_af_val, 
        
        "acmg_disease_association": acmg_disease,
        "Special_Guidance": has_guidance,
        
        "clinvar_trait": clinvar_phen,
        "clinvar_sig": v.clinvar_sig or "Not provided",
        "clinvar_stars": v.clinvar_stars,
        
        "disease_match": "Yes" if getattr(v, '_disease_match', False) else "No",
        "disease_match_level": getattr(v, '_disease_match_level', 'N/A'),
        "disease_match_reason": getattr(v, '_disease_match_reason', 'N/A'),
        
        "hgmd_class": v.hgmd_class or "Not provided",
        "hgmd_phenotype": hgmd_phen,
        "hgmd_publication_count": v.hgmd_publication_count,
        "hgmd_rankscore": getattr(v, 'hgmd_rankscore', "Not provided"),
        
        "internal_db_sig": v.internal_db_sig or "Not provided",
        "database_summary": db_summary,
        
        "gnomAD_AF": f"{v.gnomad_af:.6f}" if v.gnomad_af else "0.0",
        "gnomAD_version": gnomad_ver,
        
        "REVEL": revel_val,
        "CADD": cadd_val,
        "SpliceAI": splice_val,
        "AlphaMissense": alpha_val,
        
        "criteria_assigned": ";".join(v.criteria_assigned),
        "criteria_points": json.dumps(v.criteria_points, ensure_ascii=False),
        "total_points": v.total_points,
        "automated_class": v.automated_class,
        
        "in_auto_report": "Yes" if in_auto else "No",
        "in_manual_review": "Yes" if in_manual else "No",
        "manual_reasons": reasons,
        "filtered_out": "Yes" if filtered else "No"
    }
    
    return output_row

def strict_triage_and_output(
    report_candidates: List[VariantRecord], 
    audit_candidates: List[VariantRecord],
    gene_rules: Dict[str, Dict[str,Any]], 
    recs_map: Dict[str, Dict[str,Any]], 
    outdir: str,
    family_id_input: str = "Unknown"
):
    os.makedirs(outdir, exist_ok=True)
    family_id = family_id_input

    auto_rows = []
    manual_rows = []
    manual_ids = set()
    auto_ids = set()
    
    def get_key(v): return f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.sample}"
    
    # 1. Mark QC status
    for v in report_candidates:
        v._qc_passed = True
    qc_passed_set = {get_key(v) for v in report_candidates}

    for v in audit_candidates:
        if get_key(v) in qc_passed_set:
            v._qc_passed = True
        else:
            v._qc_passed = False
            if not any("QC" in r for r in v.manual_reasons):
                v.manual_reasons.append("Filtered: Failed Strict QC")

    # 2. Merge Lists for Grouping
    all_variants = report_candidates + audit_candidates
    
    by_sample_gene = defaultdict(lambda: defaultdict(list))
    for v in all_variants:
        by_sample_gene[v.sample][v.gene].append(v)
        
    for sample, genes in by_sample_gene.items():
        for gene, varlist in genes.items():
            if not varlist: continue

            rule = gene_rules.get(gene, {})
            moi = str(rule.get("moi","")).upper()
            acmg_disease = rule.get("disease", "")
            spec_rules = recs_map.get(gene, {})
            
            sample_sex = varlist[0].sample_sex
            # Determine inheritance type from MOI field
            is_ar = "AR" in moi or "RECESSIVE" in moi
            is_xlr = "XLR" in moi or "X-LINKED RECESSIVE" in moi
            is_xld = "XLD" in moi or "X-LINKED DOMINANT" in moi
            is_ad = "AD" in moi or "DOMINANT" in moi or "AUTOSOMAL DOMINANT" in moi or not (is_ar or is_xlr or is_xld)

            # For AD and XLD genes: DO NOT need two variants (single P/LP is sufficient)
            # For AR genes: always need two variants (unless homozygous)
            # For XLR in females: need two variants; in males: hemizygous (one variant)
            requires_biallelic = False
            if is_ar:
                requires_biallelic = True  # AR genes always require two variants
            elif is_xlr and sample_sex == "Female":
                requires_biallelic = True  # XLR in females requires two variants
            elif is_xlr and sample_sex == "Male":
                requires_biallelic = False  # XLR in males - hemizygous check
            # For AD and XLD genes: requires_biallelic = False (default)

            is_xl_male = ("X" in moi or "X-LINKED" in moi) and sample_sex == "Male"
            
            # --- Pre-calculation & Filter Benign ---
            biologically_valid = []
            for v in varlist:
                v._in_auto = False; v._in_manual = False; v._filtered = True
                
                # A. AGGRESSIVE HQ BENIGN FILTERING - ALWAYS APPLY FIRST
                # Filter variants with high-quality benign evidence REGARDLESS of gene mismatches
                is_clinvar_benign = False
                
                # Check ClinVar for high-quality benign evidence (≥2 stars, pure benign)
                if v.clinvar_stars and v.clinvar_stars.isdigit() and int(v.clinvar_stars) >= 2:
                    sig_lower = v.clinvar_sig.lower()
                    # Only pure benign classifications (no pathogenic mentions, no conflicts)
                    if ("benign" in sig_lower or "likely benign" in sig_lower) and \
                    "pathogenic" not in sig_lower and \
                    "conflict" not in sig_lower and \
                    "conflicting" not in sig_lower:
                        is_clinvar_benign = True

                # Check Internal Database for benign evidence
                if v.internal_db_sig:
                    internal_sig_lower = v.internal_db_sig.lower()
                    if ("benign" in internal_sig_lower or "likely benign" in internal_sig_lower):
                        is_clinvar_benign = True

                # Check HGMD for benign classifications (DP, FP, DFP) with sufficient publications
                if v.hgmd_class in ["DP", "FP", "DFP"]:  # HGMD benign classes
                    if v.hgmd_publication_count >= 2:  # Only with sufficient publications
                        is_clinvar_benign = True
                                
                # CRITICAL: Filter out variant if ANY high-quality benign evidence exists
                # This overrides ANY algorithm scores (even Pathogenic with high points)
                if is_clinvar_benign:
                    v.manual_reasons.append("Filtered: HQ Benign")
                    continue
                
                # 1. GENE MISMATCH HANDLING - AFTER benign filtering
                # Detect and handle cases where database gene annotations don't match variant gene
                gene_mismatch_detected = False
                gene_mismatch_sources = []

                for reason in v.manual_reasons:
                    reason_lower = reason.lower()
                    if "gene mismatch" in reason_lower:
                        gene_mismatch_detected = True
                        # Identify which database has the mismatch
                        if "clinvar" in reason_lower:
                            gene_mismatch_sources.append("ClinVar")
                        elif "hgmd" in reason_lower:
                            gene_mismatch_sources.append("HGMD")
                        # Internal DB doesn't store gene info, so no mismatch possible

                if gene_mismatch_detected:
                    # Clear data from databases with gene mismatches
                    # This ensures algorithm classification takes precedence
                    for source in gene_mismatch_sources:
                        if source == "ClinVar":
                            v.clinvar_sig = ""
                            v.clinvar_stars = "0"
                            v.clinvar_trait = ""
                            v.has_strong_pathogenic_evidence = False
                            v.is_hq_pathogenic_db = False
                        elif source == "HGMD":
                            v.hgmd_class = ""
                            v.hgmd_phen = ""
                            v.hgmd_publication_count = 0
                            v.has_strong_pathogenic_evidence = False
                            v.has_weak_pathogenic_evidence = False
                            v.is_hq_pathogenic_db = False
                    
                    # Keep only one clean warning
                    v.manual_reasons = [r for r in v.manual_reasons if "gene mismatch" not in r.lower()]
                    v.manual_reasons.append(f"Gene mismatch in {', '.join(gene_mismatch_sources)} - ignoring their evidence")
                
                # B. Calculate disease match AFTER gene mismatch handling
                traits = []
                if v.clinvar_trait and v.clinvar_trait.lower() not in ["not provided", ".", "not_provided"]:
                    traits.append(v.clinvar_trait)
                if v.hgmd_phen and v.hgmd_phen.lower() not in ["not provided", ".", "not_provided"]:
                    traits.append(v.hgmd_phen)
                full_trait = " | ".join(traits)
                
                is_match, lvl, reason = check_disease_match(acmg_disease, full_trait, [spec_rules] if spec_rules else None, rule)
                v._disease_match = is_match
                v._disease_match_level = lvl
                v._disease_match_reason = reason
                v._has_db_evidence = (len(traits) > 0)
                if not v._has_db_evidence and v.internal_db_sig: v._has_db_evidence = True

                # C. RE-CALCULATE HQ STATUS AFTER GENE MISMATCH HANDLING
                # This is critical - database classifications are ignored for mismatched genes
                is_hq = False
                # Only consider databases WITHOUT gene mismatches
                if not gene_mismatch_detected or "ClinVar" not in gene_mismatch_sources:
                    if (v.clinvar_stars and v.clinvar_stars.isdigit() and int(v.clinvar_stars) >= 2):
                        sig_low = v.clinvar_sig.lower()
                        if ("pathogenic" in sig_low or "likely pathogenic" in sig_low) and "benign" not in sig_low:
                            is_hq = True
                
                if not gene_mismatch_detected or "HGMD" not in gene_mismatch_sources:
                    if (v.hgmd_class == "DM" and v.hgmd_publication_count >= 2):
                        is_hq = True
                
                # Internal DB always considered (no gene info to mismatch)
                if v.internal_db_sig and ("pathogenic" in v.internal_db_sig.lower() or "likely pathogenic" in v.internal_db_sig.lower()):
                    is_hq = True
                
                v._is_high_quality_db = is_hq
                
                # D. DATABASE CLASSIFICATION - ONLY IF NO GENE MISMATCH
                db_classification = None
                db_source = None
                
                # Only use database classification if:
                # 1. No gene mismatch OR mismatch was in other database
                # 2. Variant is high-quality
                # 3. Disease matches
                if is_hq and (not gene_mismatch_detected or v.internal_db_sig):
                    # ClinVar (if no ClinVar gene mismatch)
                    if not gene_mismatch_detected or "ClinVar" not in gene_mismatch_sources:
                        if v.clinvar_stars and v.clinvar_stars.isdigit() and int(v.clinvar_stars) >= 2:
                            sig_lower = v.clinvar_sig.lower()
                            if "pathogenic" in sig_lower and "benign" not in sig_lower and "conflict" not in sig_lower:
                                db_classification = "Pathogenic"
                                db_source = "ClinVar"
                            elif "likely pathogenic" in sig_lower and "benign" not in sig_lower:
                                db_classification = "Likely pathogenic"
                                db_source = "ClinVar"

                    # HGMD (if no HGMD gene mismatch)
                    if not gene_mismatch_detected or "HGMD" not in gene_mismatch_sources:
                        if not db_classification and v.hgmd_class == "DM" and v.hgmd_publication_count >= 2:
                            db_classification = "Pathogenic"
                            db_source = "HGMD"

                    # Internal DB (always considered - no gene info)
                    if not db_classification and v.internal_db_sig:
                        internal_lower = v.internal_db_sig.lower()
                        if "pathogenic" in internal_lower:
                            db_classification = "Pathogenic"
                            db_source = "Internal_DB"
                        elif "likely pathogenic" in internal_lower:
                            db_classification = "Likely pathogenic"
                            db_source = "Internal_DB"

                # Apply database classification if available and disease matches
                if db_classification and v._disease_match:
                    v.automated_class = db_classification
                    v.manual_reasons.append(f"Using {db_source} classification: {db_classification}")
                
                # E. Gene Rules
                if spec_rules.get("missense_only") and "missense" not in v.consequence.lower():
                    v.manual_reasons.append("Filtered: Gene rule (missense only)")
                    continue
                if spec_rules.get("hemizygous_only_male") and sample_sex == "Female":
                    v.manual_reasons.append("Filtered: Gene rule (males only)")
                    continue
                
                biologically_valid.append(v)

            # --- Define Drivers and Passengers ---
            drivers = []
            for v in biologically_valid:
                if not v._qc_passed: 
                    continue  # Failed QC
                # Also check if variant has other filtering reasons
                if "Filtered:" in ";".join(v.manual_reasons):
                    continue  # Already filtered for other reasons
                if v._is_high_quality_db or v.automated_class in ["Pathogenic", "Likely pathogenic"]:
                    drivers.append(v)

            passengers = []
            for v in biologically_valid:
                if "Filtered:" in ";".join(v.manual_reasons):
                    continue
                if v.automated_class in ["Pathogenic", "Likely pathogenic", "VUS"]:
                    passengers.append(v)

            # --- AUTOMATIC LOGIC ---
            for v in drivers:
                if get_key(v) in manual_ids: continue
                
                can_be_automatic = True
                
                # ENHANCED CHECK:
                is_recessive_mode = requires_biallelic or ("AR" in moi) or ("RECESSIVE" in moi)
                
                # Handle X-linked in males (they behave like Hemizygous/Dominant -> single hit is enough)
                is_male_xl = is_xl_male
                if not is_male_xl and ("X" in moi) and (sample_sex == "Male"):
                    is_male_xl = True

                # If Recessive AND NOT an X-linked Male
                if is_recessive_mode and not is_male_xl:
                    
                    # 1. Homozygotes can be reported automatically (pathogenic diagnosis)
                    is_hom = _is_hom(v.proband_gt)
                    
                    # 2. If it is Heterozygous, we BLOCK Automatic Report.
                    # Rationale: For recessive diseases, we need to see the pair (compound het) or the full set of variants.
                    # Sending one variant to "Auto" breaks the triage logic in Manual Review.
                    # Therefore, ALL Het variants for AR genes must go to Manual Review.
                    if not is_hom:
                        can_be_automatic = False
                        
                        # Add explanation
                        reason_msg = "Auto Blocked: Recessive Het variant (requires pairing check in Manual)"
                        if not any(reason_msg in r for r in v.manual_reasons):
                            v.manual_reasons.append(reason_msg)

                # If automatic reporting is blocked, skip this variant.
                if not can_be_automatic:
                    continue

                # Case A: HQ DB (but check for gene mismatches first)
                if v._is_high_quality_db:
                    # Check if HQ status is based on mismatched databases
                    has_gene_mismatch = any("Gene mismatch" in r for r in v.manual_reasons)
                    
                    if not has_gene_mismatch:
                        match_ok = v._disease_match or (v._disease_match_reason == "no_db_phenotype")
                        
                        if match_ok and not v.is_hq_conflict_db:
                            v._in_auto = True; v._filtered = False
                            auto_ids.add(get_key(v))
                    # If there's a gene mismatch, we fall through to algorithm classification
                
                # Case B: Algo High Score (MAIN PATH for gene mismatch cases)
                if not v._in_auto:  # Only if not already auto from HQ DB
                    if ((v.automated_class == "Pathogenic" and v.total_points >= 10) or 
                        (v.automated_class == "Likely pathogenic" and v.total_points >= 6)):
                        # For gene mismatch cases, be more lenient with disease match
                        if gene_mismatch_detected:
                            # Algorithm classification takes precedence when databases have gene mismatches
                            match_ok = v._disease_match or (v._disease_match_reason == "no_db_phenotype") or (not v._has_db_evidence) or True
                        else:
                            match_ok = v._disease_match or (v._disease_match_reason == "no_db_phenotype") or (not v._has_db_evidence)
                        
                        if match_ok:
                            v._in_auto = True; v._filtered = False
                            auto_ids.add(get_key(v))

            # --- MANUAL LOGIC ---
            for v in drivers:
                # Skip if already in manual or auto
                if get_key(v) in manual_ids:
                    continue
                # 1. HQ Conflict
                if v._is_high_quality_db and v.is_hq_conflict_db:
                    if get_key(v) not in manual_ids:
                        v.manual_reasons.append(f"HQ Conflict: {v.conflict_details}")
                        v._in_manual = True; v._filtered = False; v._in_auto = False
                        manual_ids.add(get_key(v))
                        if get_key(v) in auto_ids: auto_ids.remove(get_key(v))

                # 2. HQ Disease Mismatch - REVISED LOGIC
                if v._is_high_quality_db and not v._disease_match:
                    if v._disease_match_reason == "no_db_phenotype":
                        if v.automated_class in ["Pathogenic", "Likely pathogenic"] and v.total_points >= 6:
                            v._in_auto = True
                            v._filtered = False
                            v.manual_reasons.append("Auto: HQ DB + good algorithmic score (no phenotype data)")
                            auto_ids.add(get_key(v))
                            continue
                        else:
                            continue
                    
                    elif v._disease_match_reason == "mismatch":
                        if get_key(v) not in manual_ids:
                            v.manual_reasons.append(f"HQ Disease Mismatch ({v._disease_match_reason})")
                            v._in_manual = True; v._filtered = False; v._in_auto = False
                            manual_ids.add(get_key(v))
                            if get_key(v) in auto_ids: auto_ids.remove(get_key(v))
                    
                    else:
                        pass

            # 3. XL Male Heterozygous
            if is_xl_male:
                for v in drivers:
                     if _is_het(v.proband_gt):
                         if get_key(v) not in manual_ids:
                             v.manual_reasons.append("XL Male Heterozygous (Pseudo-het)")
                             v._in_manual = True; v._filtered = False; v._in_auto = False
                             manual_ids.add(get_key(v))
                             if get_key(v) in auto_ids: auto_ids.remove(get_key(v))

            # 4. Complex Recessive Cases - ONLY for genes requiring two variants
            if requires_biallelic:
                # Count only P/LP variants (exclude VUS)
                path_lp_vars = [v for v in drivers if v.automated_class in ["Pathogenic", "Likely pathogenic"]]
                
                # If more than 2 P/LP variants in a recessive gene -> complex case
                if len(path_lp_vars) > 2:
                    # ADD ONLY THE TOP 2 HIGHEST-SCORING VARIANTS to manual review
                    # Sort by total_points descending
                    sorted_vars = sorted(path_lp_vars, key=lambda x: x.total_points, reverse=True)
                    
                    # Take only top 2 variants for manual review
                    for v in sorted_vars[:2]:
                        if get_key(v) not in manual_ids and get_key(v) not in auto_ids:
                            v.manual_reasons.append(f"Manual: {len(path_lp_vars)} P/LP variants in recessive gene, reviewing top 2")
                            v._in_manual = True; v._filtered = False; v._in_auto = False
                            manual_ids.add(get_key(v))
                    
                    # DO NOT add the rest - they should be filtered out
                    for v in sorted_vars[2:]:
                        if not any("Filtered" in r for r in v.manual_reasons):
                            v.manual_reasons.append(f"Filtered: Excluded due to >2 P/LP variants in recessive gene")

            # 5. AR Pair Rescue - ONLY for genes requiring two variants (AR/XLR in females)
            if requires_biallelic:
                if len(passengers) < 2: continue
                
                # Get P/LP variants
                path_lp_vars = [v for v in passengers if v.automated_class in ["Pathogenic", "Likely pathogenic"]]
                
                # If no P/LP variants, skip
                if not path_lp_vars:
                    continue
                
                # Sort by score to get strongest P/LP
                path_lp_vars_sorted = sorted(path_lp_vars, key=lambda x: x.total_points, reverse=True)
                strongest_path = path_lp_vars_sorted[0]
                
                # Look for potential partners (store only variants, not tuples)
                potential_partners_variants = []
                
                # First: check for other P/LP variants
                if len(path_lp_vars_sorted) > 1:
                    second_strongest = path_lp_vars_sorted[1]
                    potential_partners_variants.append(second_strongest)
                
                # Second: check for strongest VUS
                else:
                    vus_vars = [v for v in passengers if v.automated_class == "VUS"]
                    logger.debug(f"VUS candidates for AR pair: {len(vus_vars)}")
                    if vus_vars:
                        vus_sorted = sorted(vus_vars, key=lambda x: x.total_points, reverse=True)
                        strongest_vus = vus_sorted[0]
                        logger.debug(f"Strongest VUS: {strongest_vus.hgvsp}, points: {strongest_vus.total_points}")
                        
                        if strongest_vus.total_points >= 1:
                            potential_partners_variants.append(strongest_vus)
                
                # Process the selected pair
                if potential_partners_variants:
                    for partner_var in potential_partners_variants:
                        # Check phasing
                        is_trans = False
                        if (_is_het(strongest_path.father_gt) and _is_het(partner_var.mother_gt)) or \
                        (_is_het(strongest_path.mother_gt) and _is_het(partner_var.father_gt)):
                            is_trans = True
                        
                        is_unknown = False
                        has_parents = (strongest_path.father_gt not in ["./.", "."] or strongest_path.mother_gt not in ["./.", "."])
                        if not has_parents and not is_trans:
                            is_unknown = True
                        
                        if is_trans or is_unknown:
                            # Add ONLY this specific pair to manual review
                            for p_var in [strongest_path, partner_var]:
                                if get_key(p_var) not in manual_ids and get_key(p_var) not in auto_ids:
                                    other = partner_var if p_var == strongest_path else strongest_path
                                    prefix = "AR Pair in trans" if is_trans else "AR Pair (Phase Unknown)"
                                    qc_note = "" if p_var._qc_passed else "[Low QC Rescued]"
                                    
                                    # SAFELY get HGVSp
                                    other_hgvsp = getattr(other, 'hgvsp', 'Not provided')
                                    p_var.manual_reasons.append(f"{prefix} {qc_note} with {other_hgvsp}")
                                    p_var._in_manual = True; p_var._filtered = False; p_var._in_auto = False
                                    
                                    manual_ids.add(get_key(p_var))
                                    if get_key(p_var) in auto_ids: 
                                        auto_ids.remove(get_key(p_var))
                            
                            # BREAK after finding first valid pair
                            break
                else:
                    logger.debug(f"No partner found for {strongest_path.hgvsp}")
    
    # --- WRITE OUTPUTS ---
    COLUMNS = [
        "sample", "Family_ID", "chrom", "pos", "ref", "alt", "gene", "Inheritance",
        "HGVSc", "HGVSp", "Exon", "Consequence", "Zygosity", "DP", "Sample_AF", 
        "acmg_disease_association", "Special_Guidance", 
        "clinvar_trait", "clinvar_sig", "clinvar_stars", 
        "disease_match", "disease_match_level", "disease_match_reason", 
        "hgmd_class", "hgmd_phenotype", "hgmd_publication_count", "hgmd_rankscore",
        "internal_db_sig", "database_summary", 
        "gnomAD_AF", "gnomAD_version", 
        "REVEL", "CADD", "SpliceAI", "AlphaMissense", 
        "criteria_assigned", "criteria_points", "total_points", "automated_class", 
        "in_auto_report", "in_manual_review", "manual_reasons", "filtered_out"
    ]

    auto_rows_data = []
    manual_rows_data = []
    all_rows_data = []

    processed_keys = set()
    for v in all_variants:
        k = get_key(v)
        if k in processed_keys: continue
        processed_keys.add(k)
        
        rule = gene_rules.get(v.gene, {})
        row = variant_row_with_disease(v, rule, family_id, recs_map)
        
        if v._in_auto: auto_rows_data.append(row)
        if v._in_manual: manual_rows_data.append(row)
        
        all_rows_data.append(row)

    auto_rows_c = collapse_output_rows(auto_rows_data)
    manual_rows_c = collapse_output_rows(manual_rows_data)
    all_rows_c = collapse_output_rows(all_rows_data)

    if HAVE_PANDAS:
        pd.DataFrame(all_rows_c).reindex(columns=COLUMNS).to_csv(os.path.join(outdir, "all_candidates.csv"), index=False)
        pd.DataFrame(auto_rows_c).reindex(columns=COLUMNS).to_csv(os.path.join(outdir, "auto_conclusions.csv"), index=False)
        pd.DataFrame(manual_rows_c).reindex(columns=COLUMNS).to_csv(os.path.join(outdir, "manual_review_list.csv"), index=False)
    
    return os.path.join(outdir, "all_candidates.csv"), \
           os.path.join(outdir, "auto_conclusions.csv"), \
           os.path.join(outdir, "manual_review_list.csv"), \
           len(auto_rows_c), len(manual_rows_c)

# ---------------------------
# Build variant record from VCF record (modified to use normalized gene matching)
# ---------------------------

def build_variant_record_from_rec(rec, csq_header: List[str], acmg_genes_normalized: set, 
                                  dad_vf: Optional[pysam.VariantFile], 
                                  mom_vf: Optional[pysam.VariantFile], 
                                  exons_df: Optional[Any], 
                                  proband_sample: Optional[str]=None) -> Optional[VariantRecord]:
    try:
        # Basic VCF fields
        if hasattr(rec, "CHROM"): # cyvcf2
            chrom = rec.CHROM; pos = rec.POS; ref = rec.REF; alts = rec.ALT
        else: # pysam
            chrom = rec.chrom; pos = rec.pos; ref = rec.ref; alts = rec.alts
            
        alt = alts[0] if alts and len(alts) > 0 else ""
        if not alt: return None
        
        chrom_norm = normalize_chrom(chrom)
    except Exception: return None

    # Extract INFO
    info = {}
    if hasattr(rec, "INFO"): info = dict(rec.INFO)
    elif hasattr(rec, "info"): 
        try: info = dict(rec.info.items())
        except: info = dict(rec.info)
    info_l = {k.lower(): v for k, v in info.items()}

    # Transcript Parsing
    csq_lines = []
    csq_raw = info_l.get("csq") or info_l.get("ann")
    if csq_raw:
        s = csq_raw
        if isinstance(s, (list, tuple)): s = ",".join(str(x) for x in s)
        for part in str(s).split(","):
            if "|" in part: csq_lines.append(part)

    best_ann = None
    for ln in csq_lines:
        ann = parse_vep_csq_line(ln, csq_header)
        gene_raw = (ann.get("SYMBOL") or ann.get("Gene") or "").strip()
        norm_gene = normalize_gene_name(gene_raw)
        if norm_gene in acmg_genes_normalized:
            is_canon = ann.get("MANE_SELECT") or ann.get("MANE") or (ann.get("CANONICAL") == "YES")
            if is_canon:
                best_ann = ann; best_ann["_NORM_GENE"] = norm_gene; best_ann["_RAW_GENE"] = gene_raw
                break
            if best_ann is None:
                best_ann = ann; best_ann["_NORM_GENE"] = norm_gene; best_ann["_RAW_GENE"] = gene_raw

    if not best_ann: return None
    norm_gene = best_ann.get("_NORM_GENE")
    gene_raw = best_ann.get("_RAW_GENE")
    transcript_id = best_ann.get("Feature", "")

    # --- PROBAND GENOTYPE PARSING (ROBUST MANUAL FALLBACK) ---
    prob_gt = "./."; prob_dp = None; prob_ad = None; sample_name = ""
    prob_af = 0.0 
    
    # 1. Identify Sample Index
    idx = 0
    if hasattr(rec, "samples"):
        s_list = list(rec.samples)
        if proband_sample and proband_sample in s_list:
            idx = s_list.index(proband_sample)
            sample_name = proband_sample
        else:
            sample_name = s_list[0] if s_list else "Proband"

    # 2. Try Standard Library Parsing (cyvcf2/pysam)
    try:
        if hasattr(rec, "genotypes"): # cyvcf2
            gt_data = rec.genotypes[idx]
            if len(gt_data) >= 2:
                s1 = str(gt_data[0]) if gt_data[0] != -1 else "."; s2 = str(gt_data[1]) if gt_data[1] != -1 else "."
                prob_gt = f"{s1}/{s2}"
            
            # DP
            dps = rec.format('DP')
            if dps is not None: 
                val = dps[idx]
                prob_dp = int(val.item() if hasattr(val, "item") else val)
            
            # AF (Try library first)
            afs = rec.format('AF')
            if afs is not None:
                val = afs[idx]
                if hasattr(val, "__len__") and len(val)>0: prob_af = float(val[0])
                else: prob_af = float(val)

        else: # pysam
            s_dat = getattr(rec, "samples", {})[sample_name]
            g = s_dat.get("GT"); prob_gt = "/".join("." if x is None else str(x) for x in g) if g else "./."
            prob_dp = s_dat.get("DP")
            
            if "AF" in s_dat:
                val = s_dat["AF"]
                if isinstance(val, (list, tuple)): prob_af = float(val[0])
                else: prob_af = float(val)
    except: 
        pass

    # 3. MANUAL STRING PARSING FALLBACK (If AF is still 0.0)
    # This bypasses library type-checking issues (e.g. Header says Array but Data is Scalar)
    if prob_af == 0.0:
        try:
            # Convert record to string to get the raw VCF line
            raw_line = str(rec).strip()
            parts = raw_line.split('\t')
            
            # VCF columns: 0=CHROM ... 8=FORMAT, 9+=SAMPLES
            if len(parts) > 9:
                fmt_str = parts[8]
                # Adjust sample index (9 is the first sample)
                samp_str = parts[9 + idx]
                
                fmt_keys = fmt_str.split(':')
                samp_vals = samp_str.split(':')
                
                # Check consistency
                if len(fmt_keys) == len(samp_vals):
                    # Extract AF
                    if "AF" in fmt_keys:
                        i_af = fmt_keys.index("AF")
                        raw_val = samp_vals[i_af]
                        # Handle "0.33,0.01" case
                        if "," in raw_val: raw_val = raw_val.split(",")[0]
                        if raw_val and raw_val != ".":
                            prob_af = float(raw_val)
                    
                    # Extract DP if missing
                    if not prob_dp and "DP" in fmt_keys:
                        i_dp = fmt_keys.index("DP")
                        raw_val = samp_vals[i_dp]
                        if raw_val and raw_val != ".":
                            prob_dp = int(raw_val)
                    
                    # Extract GT if missing
                    if prob_gt == "./." and "GT" in fmt_keys:
                        i_gt = fmt_keys.index("GT")
                        prob_gt = samp_vals[i_gt]

                    # Fallback AD calculation if AF is missing entirely
                    if prob_af == 0.0 and "AD" in fmt_keys:
                         i_ad = fmt_keys.index("AD")
                         raw_ad = samp_vals[i_ad] # e.g. "14,14"
                         if "," in raw_ad:
                             ads = [int(x) if x.isdigit() else 0 for x in raw_ad.split(",")]
                             if len(ads) >= 2 and sum(ads) > 0:
                                 prob_af = ads[1] / sum(ads)
                                 if not prob_dp: prob_dp = sum(ads)
        except Exception: 
            pass

    # Parents Lookup
    def get_gt_from_handle(vf, c, p, r):
        if not vf: return "./."
        try:
            q_chr = c if c in vf.header.contigs else c.replace("chr","")
            candidates = [c, c.replace("chr", ""), "chr"+c.replace("chr","")]
            valid = [x for x in candidates if x in vf.header.contigs]
            if not valid: return "./."
            for q in valid:
                try:
                    for rec_p in vf.fetch(q, p-1, p):
                        if int(rec_p.pos) == int(p) and rec_p.ref == r: 
                            s = list(rec_p.samples.values())[0]
                            g = s.get("GT")
                            return "/".join("." if x is None else str(x) for x in g) if g else "./."
                except: continue
        except: pass
        return "./."

    father_gt = get_gt_from_handle(dad_vf, chrom_norm, pos, ref)
    mother_gt = get_gt_from_handle(mom_vf, chrom_norm, pos, ref)

    # Helper for scores
    def get_val_anywhere(keys, return_type=float):
        def _parse_raw(val):
            if isinstance(val, (list, tuple)): val = val[0]
            s = str(val)
            if not s or s == ".": return None
            if "&" in s:
                nums = []
                for p in s.split("&"):
                    if p and p != ".":
                        try: nums.append(float(p))
                        except: pass
                if not nums: return None 
                if return_type == float: return max(nums)
                return str(max(nums))
            try: return return_type(s)
            except: return None

        for k in keys:
            if k in best_ann and best_ann[k]:
                res = _parse_raw(best_ann[k])
                if res is not None: return res
            if k.lower() in info_l:
                res = _parse_raw(info_l[k.lower()])
                if res is not None: return res
        return None

    # NMD / Scores
    vep_nmdesc = ""
    if "NMD" in best_ann and best_ann["NMD"]: vep_nmdesc = best_ann["NMD"]
    if not vep_nmdesc: vep_nmdesc = get_val_anywhere(["NMD", "NMDEscPredictor", "nmdesc_predictor"], str)

    splice_keys = ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_DS_max", "SpliceAI"]
    splice_vals = []
    for key in splice_keys:
        val = get_val_anywhere([key])
        if val is not None: splice_vals.append(val)
    spliceai = max(splice_vals) if splice_vals else None

    gnomad_af = None; gdet = {}
    _, gm = extract_gnomad_af_from_info(info_l)
    if gm.get("max_af", 0) > 0: gnomad_af = gm["max_af"]; gdet = gm
    else:
        _, gm2 = extract_gnomad_af_from_info(best_ann)
        if gm2.get("max_af", 0) > 0: gnomad_af = gm2["max_af"]; gdet = gm2

    vr = VariantRecord(
        chrom=chrom_norm, pos=int(pos), ref=ref, alt=str(alt), sample=sample_name,
        gene=norm_gene, gene_raw=gene_raw,
        transcript=transcript_id, consequence=best_ann.get("Consequence",""),
        hgvsc=best_ann.get("HGVSc",""), hgvsp=best_ann.get("HGVSp",""),
        exon=best_ann.get("EXON",""),
        gnomad_af=gnomad_af, gnomad_details=gdet,
        revel=get_val_anywhere(["REVEL_score", "revel"]), 
        alpha_missense=get_val_anywhere(["AlphaMissense_score"]), 
        spliceai=spliceai, 
        cadd=get_val_anywhere(["CADD_phred"]), 
        nmdesc_predictor=vep_nmdesc, 
        proband_gt=prob_gt, father_gt=father_gt, mother_gt=mother_gt,
        proband_dp=prob_dp, proband_ad=prob_ad, proband_af=prob_af 
    )
    
    internal_nmd = ""
    if exons_df is not None: 
        internal_nmd = compute_nmd_internal(vr, exons_df)
    
    if internal_nmd:
        vr.nmd = internal_nmd
        vr.nmdesc_predictor = f"Internal: {internal_nmd}"
    else:
        vr.nmd = vep_nmdesc
    vr._raw_ann = best_ann
    return vr

# ---------------------------
# Pre-analysis - quality and biological sex control
# ---------------------------

def perform_pre_analysis_qc(vcf_path: str, proband_id: str) -> Dict[str, Any]:
    """
    Performs pre-analysis quality control.
    Infers biological sex based on Y-chromosome presence and X-chromosome heterozygosity.
    """
    warnings = []
    inferred_sex = "Unknown"
    
    if HAVE_PYSAM:
        try:
            vf = pysam.VariantFile(vcf_path)
            chroms = list(vf.header.contigs)
            
            # Identify Sex Chromosomes
            chr_x_name = next((c for c in chroms if c.lower() in ('chrx', 'x', '23')), None)
            chr_y_name = next((c for c in chroms if c.lower() in ('chry', 'y', '24')), None)
            
            # 1. Check Y Chromosome (Strongest indicator for Male)
            has_y_variants = False
            if chr_y_name:
                try:
                    # Check if there are ANY variants on Y
                    for rec in vf.fetch(chr_y_name):
                        has_y_variants = True
                        break # Found one, enough to suspect Male
                except: pass
            
            if has_y_variants:
                inferred_sex = "Male"
                logger.info(f"QC Sex Check: Inferred Male (Variants found on {chr_y_name})")
                return {"warnings": warnings, "sex": inferred_sex}

            # 2. Check X Heterozygosity (if no Y or Y check failed)
            if chr_x_name:
                het_count = 0
                hom_count = 0
                total_x = 0
                
                try:
                    for rec in vf.fetch(chr_x_name):
                        if len(rec.samples) == 0: break
                        
                        # Determine sample index/name
                        s_name = proband_id if proband_id in rec.samples else rec.samples[0].name
                        gt = rec.samples[s_name]['GT']
                        
                        # Check Genotype
                        if len(gt) >= 2:
                            if gt[0] != gt[1]: het_count += 1
                            else: hom_count += 1
                            total_x += 1
                        
                        # Limit to avoid slow processing on whole genomes
                        if total_x > 2000: break
                except Exception: pass
                
                if total_x > 0:
                    ratio = het_count / total_x
                    logger.info(f"QC Sex Check: ChrX Stats - Total: {total_x}, Het: {het_count}, Ratio: {ratio:.2f}")
                    
                    if ratio > 0.15: # Significant portion is Het -> Female
                        inferred_sex = "Female"
                    elif ratio < 0.10: # Mostly Hom/Hemi -> Male
                        inferred_sex = "Male"
                    else:
                        inferred_sex = "Unknown" # Ambiguous
                else:
                     logger.warning("QC Sex Check: No variants found on ChrX to infer sex.")
            
        except Exception as e:
            logger.warning(f"QC Sex Check failed: {e}")

    return {"warnings": warnings, "sex": inferred_sex}

# ---------------------------
# Main pipeline (glue) - includes added minimal tests/logging before parsing
# ---------------------------

def deduplicate_candidates(candidates: List[VariantRecord]) -> List[VariantRecord]:
    """
    Deduplicate VariantRecord objects by sample+coord+alleles+normalized gene.

    Prefer:
    - records annotated as MANE_SELECT / CANONICAL
    - records with higher total_points
    - records with stronger computational evidence (REVEL + scaled CADD)

    Returns a list with one representative per unique key.
    """
    keymap: Dict[Tuple[str,str,int,str,str,str], VariantRecord] = {}
    
    for v in candidates:
        key = (v.sample or "",
            v.chrom or "",
            int(v.pos),
            v.ref or "",
            v.alt or "",
            (v.gene or "").upper())
        
        if key not in keymap:
            keymap[key] = v
            continue
            
        existing = keymap[key]
        
        # Prefer MANE / CANONICAL (check raw annotation if present)
        e_ann = getattr(existing, "_raw_ann", {}) or {}
        v_ann = getattr(v, "_raw_ann", {}) or {}
        e_can = str(e_ann.get("MANE_SELECT") or e_ann.get("CANONICAL") or "").upper()
        v_can = str(v_ann.get("MANE_SELECT") or v_ann.get("CANONICAL") or "").upper()
        
        if v_can in ("YES","Y","TRUE") and not (e_can in ("YES","Y","TRUE")):
            keymap[key] = v
            continue
            
        # Prefer higher total_points (if available, though calculation happens later usually)
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
            
        # Tie-break by computational evidence (REVEL + scaled CADD)
        def safe_float(val):
            if val is None: return 0.0
            try: return float(val)
            except (ValueError, TypeError): return 0.0
        
        def score(rec: VariantRecord) -> float:
            r = safe_float(getattr(rec, "revel", None))
            c = safe_float(getattr(rec, "cadd", None))
            return r + (c / 25.0)
        
        v_score = score(v)
        e_score = score(existing)
        
        if v_score > e_score:
            keymap[key] = v
            continue
            
        # Keep existing otherwise
    
    return list(keymap.values())

def process_vcf(proband_vcf: str,
                father_vcf: Optional[str],
                mother_vcf: Optional[str],
                outdir: str,
                run_vep_dbnsfp: bool,
                dbnsfp_path: Optional[str],
                dbnsfp_fields: Optional[List[str]],
                vep_cmd: str,
                vep_cache: Optional[str],
                fasta: Optional[str],
                vep_extra: List[str],
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
                family_id: str = "Unknown"):
    os.makedirs(outdir, exist_ok=True)

    # 1. PREPARE INPUT VCF (Compress & Index FIRST)
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

    # 2. PERFORM QC
    recs_map = load_recommendations_table(acmg_recs) if acmg_recs else {}
    qc_res = perform_pre_analysis_qc(input_vcf_for_tools, proband_sample) 
    inferred_sex = qc_res.get("sex", "Unknown")
    logger.info(f"Pre-analysis QC complete. Sex: {inferred_sex}")

    # 3. GENE BED FILTERING
    if gene_bed and os.path.exists(gene_bed):
        logger.info("Gene BED provided: %s. Expanding by 20kb and merging intervals...", gene_bed)
        expanded_bed = os.path.join(outdir, "expanded_gene_regions.bed")
        intervals = defaultdict(list)
        with open(gene_bed) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#"): continue
                parts = re.split(r"\s+", ln)
                if len(parts) < 3: continue
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
    
    # Run VEP
    if run_vep_dbnsfp:
        if not dbnsfp_path or not os.path.exists(dbnsfp_path):
            raise FileNotFoundError(f"dbNSFP file not found: {dbnsfp_path}")
        
        out_v = os.path.join(outdir, os.path.basename(working_vcf).replace(".vcf","").replace(".gz","") + ".dbnsfp.vep.vcf.gz")
        dbnsfp_plugin_string = choose_dbnsfp_fields(dbnsfp_path, dbnsfp_fields)
        
        try:
            annotated_vcf = run_vep_robust(
                in_vcf=working_vcf,
                out_vcf=out_v,
                vep_cmd=vep_cmd,
                vep_cache=vep_cache,
                fasta=fasta,
                extra_args=vep_extra,
                dbnsfp_plugin_str=dbnsfp_plugin_string
            )
            working_vcf = annotated_vcf
        except RuntimeError as e:
            logger.error("Critical error during VEP annotation. Aborting pipeline.")
            sys.exit(1)

    # Load Resources
    fasta_handle = None
    if fasta and os.path.exists(fasta) and HAVE_PYSAM:
        try: fasta_handle = pysam.FastaFile(fasta)
        except: pass

    gene_rules = load_acmg_table(acmg_table) if HAVE_PANDAS else {}
    acmg_genes_normalized = set(g for g in gene_rules.keys() if gene_rules.get(g,{}).get("reportable", False))
    
    dbm = DatabaseManager(db_paths, gnomad_v2=gnomad_v2, gnomad_v3=gnomad_v3, gnomad_v4=gnomad_v4, internal_db_path=internal_db_path)
    exons_df = load_exons_table(exons_file) if exons_file else None
    csq_header = parse_csq_header_from_vcf(annotated_vcf)

    # Parse VCF
    vcf_reader = None
    if HAVE_CYVCF2: vcf_reader = VCF(annotated_vcf)
    else: vcf_reader = pysam.VariantFile(annotated_vcf)

    dad_vf = pysam.VariantFile(father_vcf) if father_vcf and os.path.exists(father_vcf) and HAVE_PYSAM else None
    mom_vf = pysam.VariantFile(mother_vcf) if mother_vcf and os.path.exists(mother_vcf) and HAVE_PYSAM else None

    candidates: List[VariantRecord] = []
    logger.info("Parsing VCF...")
    
    # Parsing loop
    if HAVE_CYVCF2 and isinstance(vcf_reader, VCF):
        for rec in vcf_reader:
            vr = build_variant_record_from_rec(rec, csq_header, acmg_genes_normalized, dad_vf, mom_vf, exons_df, proband_sample=proband_sample)
            if not vr: continue
            vr.sample_sex = inferred_sex
            if (vr.gnomad_af is None or vr.gnomad_af == 0.0): dbm.annotate_gnomad_for_variant(vr)
            try: vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
            except: pass
            candidates.append(vr)
    else:
        for rec in vcf_reader.fetch():
            # (Simplified fallback for pysam logic, relying on build_variant_record_from_rec logic mainly)
            vr = build_variant_record_from_rec(rec, csq_header, acmg_genes_normalized, dad_vf, mom_vf, exons_df, proband_sample=proband_sample)
            if not vr: continue
            vr.sample_sex = inferred_sex
            if (vr.gnomad_af is None or vr.gnomad_af == 0.0): dbm.annotate_gnomad_for_variant(vr)
            try: vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
            except: pass
            candidates.append(vr)

    candidates = deduplicate_candidates(candidates)
    if proband_sample:
        for v in candidates: v.sample = proband_sample

    # Scoring
    engine = ACMGEngine(dbm, gene_rules, exons_df, THRESH["SPLICEAI_MODERATE"])  
    for v in candidates:
        assigned, pts, auto_class, manual_flag, manual_reasons = engine.evaluate_variant(v, fasta_handle=fasta_handle)
        v.criteria_assigned = assigned
        v.criteria_points = pts
        v.automated_class = auto_class
        v.manual_review = manual_flag
        v.manual_reasons = manual_reasons
        v.total_points = sum(pts.values())
        
        # Recount criteria groups
        v.pathogenic_criteria_count = sum(1 for x in assigned if any(p in x for p in ["PVS","PS","PM","PP"]))
        v.benign_criteria_count = sum(1 for x in assigned if any(p in x for p in ["BA","BS","BP"]))
        
        if "BA1" not in v.criteria_assigned:
            if (v.automated_class.startswith("Path") or v.automated_class.startswith("Likely path")) and v.pathogenic_criteria_count < MIN_PATHOGENIC_CRITERIA:
                v.automated_class = "VUS"
                v.manual_reasons.append(f"Insufficient pathogenic criteria (<{MIN_PATHOGENIC_CRITERIA})")
    
    if fasta_handle: fasta_handle.close()
    
    engine.apply_pm3_batch(candidates) 

    # --- UPDATED STRICT QC (Ignores GT string, checks DP/AF only) ---
    filtered_for_report = [] 
    dropped_count = 0
    
    for v in candidates:
        if not STRICT_QC_ENABLED or not strict_qc:
            filtered_for_report.append(v); continue
        
        dp = v.proband_dp or 0
        af = v.proband_af or 0.0 # Using calculated AB/AF
        
        # QC Logic: DP >= 15 AND AF >= 0.2
        # We ignore GT (./. or 0/1) because if AF > 0.2, it's real enough.
        
        if dp >= STRICT_MIN_DP and af >= STRICT_MIN_AF:
            filtered_for_report.append(v)
        else:
            dropped_count += 1
            
    logger.info("After strict QC filtering: %d kept for Report, %d low quality (audit only)", len(filtered_for_report), dropped_count)

    # Output
    all_path, auto_path, manual_path, cnt_auto, cnt_manual = strict_triage_and_output(
        report_candidates=filtered_for_report,  
        audit_candidates=candidates,            
        gene_rules=gene_rules,
        recs_map=recs_map, 
        outdir=outdir,
        family_id_input=family_id
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

    # --- PRINT TERMINAL STATISTICS ---
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
# ACMG table loader (full, returns normalized gene keys)
# ---------------------------
def load_acmg_table(path: str) -> Dict[str, Dict[str,Any]]:
    """
    Loads the ACMG table with aggressive column detection to handle headers 
    with typos (e.g., 'Phentyope'), whitespace, or unusual separators.
    """
    rules: Dict[str, Dict[str,Any]] = {}
    if not path or not os.path.exists(path):
        logger.error("ACMG table not found: %s", path)
        return rules
    if not HAVE_PANDAS:
        logger.error("pandas required to load ACMG table")
        return rules
    
    try:
        # Load file. engine='python' and sep=None allow auto-detection of delimiters.
        df = pd.read_csv(path, sep=None, engine='python', dtype=str, on_bad_lines='skip').fillna("")
        
        # DEBUG: Log detected columns to help diagnosis
        logger.info("--- DEBUG: ACMG Table Columns Detected ---")
        logger.info(list(df.columns))
        logger.info("------------------------------------------")
        
    except Exception as e:
        logger.warning("Failed to read ACMG table %s: %s", path, e)
        return rules

    # Create a normalized map of columns (lowercase, stripped of whitespace)
    cols = {c.lower().strip(): c for c in df.columns}
    
    # Helper to find a column by partial matching (substring)
    def find_col_fuzzy(keywords):
        # 1. Try exact matches in normalized map
        for k in keywords:
            if k in cols: return cols[k]
        
        # 2. Try checking if keyword is a SUBSTRING of the column header
        for k in keywords:
            for real_col_lower in cols.keys():
                if k in real_col_lower:
                    return cols[real_col_lower]
        return None

    # Identify Columns using broad keywords
    gene_col = find_col_fuzzy(["gene symbol", "gene", "symbol", "hgnc"])
    
    # Explicitly looking for "phentyope" due to the typo in your file
    disease_col = find_col_fuzzy(["disease/phentyope", "phentyope", "disease", "phenotype", "condition", "disorder"])
    
    moi_col = find_col_fuzzy(["inheritance", "mode of inheritance", "moi"])
    variants_col = find_col_fuzzy(["variants to report", "variants", "variants_note", "note", "comments"])

    logger.info(f"MAPPED COLUMNS: Gene='{gene_col}', Disease='{disease_col}', MOI='{moi_col}'")

    if gene_col is None:
        logger.error("CRITICAL: Could not find 'Gene' column in ACMG table.")
        return rules

    truthy = set(["1","yes","true","y","t","included","include","report","reportable","x"])
    any_reportable = False

    for _, r in df.iterrows():
        gene_raw = str(r.get(gene_col, "")).strip()
        if not gene_raw:
            continue
        
        gene_norm = normalize_gene_name(gene_raw)
        moi = str(r.get(moi_col, "")).strip() if moi_col else ""
        
        # Extract disease and clean up internal newlines/tabs
        disease = str(r.get(disease_col, "")).strip() if disease_col else ""
        disease = disease.replace("\r", " ").replace("\n", " ").strip()
        
        variants_note = str(r.get(variants_col, "")).strip() if variants_col else ""
        
        # Determine reportable status
        reportable = False
        # Look for explicit "report" or "version" columns
        found_flag_col = False
        for potential_col in cols:
            if "sf list" in potential_col or "version" in potential_col or "report" in potential_col:
                val = str(r.get(cols[potential_col],"")).strip()
                if val and val not in ("0", "-", ""):
                     reportable = True
                     found_flag_col = True
                     break
        
        # Fallback: if 'Variants to report' column exists and has text, assume reportable
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

        rules[gene_norm] = {
            "reportable": reportable, 
            "moi": moi, 
            "disease": disease, 
            "variants_note": variants_note, 
            "special_flags": flags, 
            "raw_row": {c: str(r.get(c, "")) for c in df.columns}
        }
        if reportable:
            any_reportable = True

    if (not any_reportable) and len(rules) > 0:
        logger.warning("No explicit reportable flags detected. Marking ALL genes in table as reportable.")
        for g in rules: rules[g]["reportable"] = True

    logger.info("Loaded ACMG rules for %d genes from %s", len(rules), path)
    return rules

def load_recommendations_table(path: str) -> Dict[str, Dict[str,Any]]:
    """
    Parses 'Gene Symbol' and 'Reporting Guidance Comment' to extract rules.
    """
    recs = defaultdict(dict)
    if not path or not os.path.exists(path):
        return recs
    
    try:
        df = pd.read_csv(path, sep=None, engine='python', dtype=str).fillna("")
        
        # Detect columns
        cols = {c.lower().strip(): c for c in df.columns}
        gene_col = cols.get("gene symbol") or cols.get("gene")
        comment_col = cols.get("reporting guidance comment") or cols.get("comment") or cols.get("guidance")
        
        if not gene_col: return recs

        for _, r in df.iterrows():
            g = normalize_gene_name(r.get(gene_col, ""))
            if not g: continue
            
            comment = str(r.get(comment_col, "")).lower()
            
            # Extract Rules from text
            rules = {
                "requires_biallelic": False,
                "missense_only": False,
                "hemizygous_only_male": False,
                "exclusions": []
            }
            
            # Biallelic / Recessive logic
            if "biallelic" in comment or "only reportable if there are two" in comment:
                rules["requires_biallelic"] = True
            
            # Missense only (APOB)
            if "only p/lp missense variants" in comment or "missense only" in comment:
                rules["missense_only"] = True
                
            # Hemizygous males (ABCD1)
            if "hemizygous p/lp variants in males are reportable" in comment and "single heterozygous" in comment and "females are not reportable" in comment:
                rules["hemizygous_only_male"] = True # Report male hemi, ignore female het
            
            # Exclusions (APC -> GAPPS)
            if "not recommended for return" in comment:
                # Basic parsing for "not recommended" context
                # This is hard to generalize perfectly without NLP, but we can look for specific keywords
                if "gapps" in comment: rules["exclusions"].append("gastric adenocarcinoma")
            
            recs[g] = rules
            
        logger.info(f"Loaded reporting guidance for {len(recs)} genes")
    except Exception as e:
        logger.warning(f"Failed to load recommendations table: {e}")
    
    return recs


# ---------------------------
# CLI
# ---------------------------
def build_parser():
    p = argparse.ArgumentParser(prog="acmg_sf_classifier.py", description="ACMG SF classifier - final ready-to-run")
    p.add_argument("--proband", help="Proband VCF (bgz recommended)")
    p.add_argument("--father", help="Father VCF (optional)")
    p.add_argument("--mother", help="Mother VCF (optional)")
    p.add_argument("--outdir", default="acmg_sf_results")
    p.add_argument("--run-vep-dbnsfp", action="store_true", help="Run VEP with dbNSFP plugin")
    p.add_argument("--dbnsfp", help="dbNSFP archive for VEP plugin (optional)")
    p.add_argument("--dbnsfp-fields", default="REVEL_score,SpliceAI_pred_DS_max,CADD_phred,AlphaMissense,gnomAD_genomes_AF", help="Comma-separated dbNSFP field names")
    p.add_argument("--vep", default="vep", help="VEP executable")
    p.add_argument("--vep-cache", default=DEFAULT_VEP_CACHE, help="VEP cache dir")
    p.add_argument("--fasta", help="Reference FASTA for VEP")
    p.add_argument("--vep-extra", nargs="*", default=[], help="Extra args to pass to VEP")
    p.add_argument("--acmg-table", required=True, help="ACMG SF table CSV/TSV")
    p.add_argument("--db-paths-json", help="JSON file to override DB paths")
    p.add_argument("--ttn-meta", help="TTN meta-exon CSV (optional)")
    p.add_argument("--exons-file", help="Transcript/exons CSV (optional)")
    p.add_argument("--aggressive", action="store_true", help="Allow DB-driven overrides for VUS->(Likely)Pathogenic (not recommended)")
    p.add_argument("--gene-bed", help="BED file with gene regions to prefilter (optional). Intervals expanded by 20kb and merged.")
    p.add_argument("--cadd-tabix", help="tabix-indexed CADD TSV (optional)")
    p.add_argument("--clinvar-tabix", help="tabix-indexed ClinVar TSV (optional)")
    p.add_argument("--dbnsfp-variants-dir", help="Directory with per-chrom dbNSFP variant files (optional)")
    p.add_argument("--dbnsfp-cache-sqlite", help="Optional SQLite file path to cache dbNSFP lookups across runs")
    p.add_argument("--proband-sample", help="If VCF is multi-sample, specify proband sample name (optional)")
    p.add_argument("--no-strict-qc", action="store_true", help="Disable strict QC filtering (not recommended)")
    p.add_argument("--gnomad-v2", help="gnomAD v2 site file (VCF or TSV bgz)")
    p.add_argument("--gnomad-v3", help="gnomAD v3 site file (VCF or TSV bgz)")
    p.add_argument("--gnomad-v4", help="gnomAD v4 site VCF (tabix'd)")
    p.add_argument("--batch-input-dir", help="Folder containing subfolders with VCF files. Runs pipeline for ALL samples and aggregates results.")
    p.add_argument("--acmg-recs", help="Secondary CSV with specific variant/gene recommendations")
    p.add_argument("--internal-db", help="Path to internal DB CSV (must be hg38 coordinates)") 
    p.add_argument("--hgmd", help="Path to HGMD Pro VCF (hg38)")
    return p

def process_single_sample_worker(payload):
    """
    Worker function for parallel processing.
    """
    p_vcf, args_dict, db_paths = payload
    
    try:
        dirname = os.path.dirname(p_vcf)
        filename = os.path.basename(p_vcf)
        dna_id = filename.split('_')[0] 
        map_id = os.path.basename(dirname) # This is the Family ID (e.g. 10239_2024)
        
        s_out = os.path.join(args_dict.outdir, "individual_results", dna_id)
        
        siblings = glob.glob(os.path.join(dirname, "*.vcf"))
        f_vcf = next((s for s in siblings if "_father" in s), None)
        m_vcf = next((s for s in siblings if "_mother" in s), None)
        
        dbnsfp_fields = [x.strip() for x in args_dict.dbnsfp_fields.split(",")] if args_dict.dbnsfp else []

        process_vcf(
            proband_vcf=p_vcf, 
            father_vcf=f_vcf, 
            mother_vcf=m_vcf, 
            outdir=s_out,
            run_vep_dbnsfp=args_dict.run_vep_dbnsfp, 
            dbnsfp_path=args_dict.dbnsfp, 
            dbnsfp_fields=dbnsfp_fields,
            vep_cmd=args_dict.vep, 
            vep_cache=args_dict.vep_cache, 
            fasta=args_dict.fasta, 
            vep_extra=args_dict.vep_extra,
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
            family_id=map_id # Pass the Family ID correctly
        )
        return (dna_id, map_id, s_out)
    except Exception as e:
        print(f"ERROR in worker for {p_vcf}: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_batch_mode(args, db_paths):
    if not HAVE_PANDAS:
        logger.error("Batch mode requires pandas.")
        sys.exit(1)

    logger.info("STARTING PARALLEL BATCH ANALYSIS")
    
    # Recursive search for proband VCFs
    search_pattern = os.path.join(args.batch_input_dir, "**", "*_proband.vcf")
    proband_files = glob.glob(search_pattern, recursive=True)
    
    if not proband_files:
        logger.error("No *_proband.vcf files found")
        sys.exit(1)

    logger.info(f"Found {len(proband_files)} samples. Using ProcessPoolExecutor.")

    tasks = [(f, args, db_paths) for f in proband_files]
    
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
        # Paths to individual results
        p_auto = os.path.join(s_out, "auto_conclusions.csv")
        p_manual = os.path.join(s_out, "manual_review_list.csv")
        p_all = os.path.join(s_out, "all_candidates.csv")

        # Helper to read safely
        def read_csv_safe(path):
            if os.path.exists(path) and os.path.getsize(path) > 1:
                try:
                    # Read as string to preserve formatting (e.g. "01" genes)
                    return pd.read_csv(path, dtype=str)
                except Exception:
                    return None
            return None

        # Read and collect
        df_a = read_csv_safe(p_auto)
        if df_a is not None and not df_a.empty:
            auto_dfs.append(df_a)

        df_m = read_csv_safe(p_manual)
        if df_m is not None and not df_m.empty:
            manual_dfs.append(df_m)

        df_all = read_csv_safe(p_all)
        if df_all is not None and not df_all.empty:
            all_dfs.append(df_all)

    # Save Final Master Files
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
        dbnsfp_fields = [x.strip() for x in args.dbnsfp_fields.split(",")] if args.dbnsfp else []
        process_vcf(
            proband_vcf=args.proband,
            father_vcf=args.father,
            mother_vcf=args.mother,
            outdir=args.outdir,
            run_vep_dbnsfp=args.run_vep_dbnsfp,
            dbnsfp_path=args.dbnsfp,
            dbnsfp_fields=dbnsfp_fields,
            vep_cmd=args.vep,
            vep_cache=args.vep_cache,
            fasta=args.fasta,
            vep_extra=args.vep_extra,
            acmg_table=args.acmg_table,
            db_paths=db_paths,
            ttn_meta=args.ttn_meta,
            exons_file=args.exons_file,
            aggressive=args.aggressive,
            gene_bed=args.gene_bed,
            cadd_tabix=args.cadd_tabix,
            clinvar_tabix=args.clinvar_tabix,
            dbnsfp_variants_dir=args.dbnsfp_variants_dir,
            dbnsfp_cache_sqlite=args.dbnsfp_cache_sqlite,
            proband_sample=args.proband_sample,
            strict_qc=not args.no_strict_qc,
            gnomad_v2=args.gnomad_v2,
            gnomad_v3=args.gnomad_v3,
            gnomad_v4=args.gnomad_v4,
            acmg_recs=args.acmg_recs,
            internal_db_path=args.internal_db
        )
    else:
        parser.print_help()

if __name__ == "__main__":
    main()