#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACMG Secondary Findings pipeline 

Changes introduced in this version:
 - normalize_gene_name() added and applied everywhere:
     * load_acmg_table creates normalized gene keys
     * ClinVar protein/splice index keys are normalized
     * build_variant_record_from_rec uses normalized gene to match ACMG gene set,
       stores raw symbol in vr._raw_ann['SYMBOL_RAW'] and vr.gene set to normalized form
     * gene_rules lookups use normalized keys
 - DatabaseManager: supports gnomAD input as VCF (pysam.VariantFile) and TSV (pysam.TabixFile).
     - TSV header detection and heuristic AF extraction implemented
 - annotate_gnomad_for_variant uses both VCF and TSV handles
 - Minimal sanity checks / logging before parsing: counts for filtered VCF and sample of normalized gene list
 - Outputs include normalized_gene column and final CSV writing collapses exact duplicates
 - Reduce duplicate rows in outputs by collapsing rows by unique key (sample,chrom,pos,ref,alt,normalized_gene)
   and merging criteria/manual_reasons where appropriate
 - Preserves prior logic for PVS1/PM5/PM3 and conservative triage

 """
from __future__ import annotations
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
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any
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
logger = logging.getLogger("acmg_sf_pipeline_final")

# ---------------------------
# Defaults & thresholds
# ---------------------------
DEFAULT_VEP_CACHE = "/home/anna/anna/ACMG_SF_Classifier/databases/vep/.vep"
DB_PATHS_DEFAULT = {
    "CLINVAR_VCF": "databases/clinvar/clinvar.vcf.gz",
    "HGMD_VCF": "databases/hgmd/hgmd_pro.vcf.gz",
}
THRESH = {
    "REVEL_SUPPORT": 0.644,
    "REVEL_MODERATE": 0.932,
    "SPLICEAI": 0.40,       # For PVS1 
    "SPLICEAI_SUPPORT": 0.20, # For PP3
    "BA1_AF": 0.05,
    "BS1_AF": 0.01,
    "PM2_STRICT": 1e-5
}

NMD_DISTANCE_TO_LAST_JUNCTION = 50

# Points & classification map 
POINTS_MAP = {
    "PVS1_VeryStrong": 8, "PVS1_Strong": 4, "PVS1_Moderate": 2, "PVS1_Supporting": 1,
    "PS1_Strong": 4, "PS1_Moderate": 2, "PS1_Supporting": 1,
    "PS2_VeryStrong": 8, "PS2_Moderate": 2,
    "PM1": 2, "PM2_Moderate": 2, "PM2_Supporting": 1,
    "PM3_Supporting": 1, "PM3_Moderate": 2, "PM3_Strong": 4, "PM3_VeryStrong": 8,
    "PM4_Moderate": 2, "PM4_Supporting": 1, "PM5_Moderate": 2, "PM5_Supporting": 1,
    "PP1_Supporting": 1, "PP1_Moderate": 2, "PP1_Strong": 4,
    "PP3_Moderate": 2, "PP3_Supporting": 1,
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
STRICT_QC_ENABLED = True
STRICT_MIN_DP = 15
STRICT_MIN_AB = 0.20

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
    proband_ab: Optional[float] = None

    db_hits: List[Tuple[str,str]] = field(default_factory=list)
    criteria_assigned: List[str] = field(default_factory=list)
    criteria_points: Dict[str,int] = field(default_factory=dict)
    total_points: int = 0
    automated_class: str = "VUS"
    manual_review: bool = False
    manual_reasons: List[str] = field(default_factory=list)
    _raw_ann: Dict[str,Any] = field(default_factory=dict)

    # --- NEW FIELDS FOR STRICT LOGIC ---
    # True ONLY if ClinVar stars >= 2 AND (P or LP) OR HGMD DM >= 2 papers (if available)
    is_hq_pathogenic_db: bool = False
    # True if HQ DB entries conflict (e.g. 2-star P vs 2-star B)
    is_hq_conflict_db: bool = False
    # True if ClinVar trait matches ACMG disease keywords
    _disease_match: bool = False

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


def check_disease_match(acmg_disease_string: str, variant_trait: str, extra_recs_trait: Optional[str] = None) -> bool:
    """
    Strict checks if the variant's phenotype matches the ACMG designated disease.
    
    Logic:
    1. If 'extra_recs_trait' is provided (from the supplementary table), matching that is sufficient.
    2. Tokenize ACMG disease string and Variant Trait string.
    3. Require meaningful overlap of tokens.
    
    Args:
        acmg_disease_string: From ACMG main table (e.g., "Dilated cardiomyopathy").
        variant_trait: From ClinVar/HGMD (e.g., "Cardiomyopathy, dilated; Heart failure").
        extra_recs_trait: From optional recommendations table.

    Returns:
        True if matched or if ACMG field is empty (assumed catch-all).
        False if distinct mismatch.
    """
    # 1. Check Extra Recs Overrides first
    if extra_recs_trait:
        # If specific instructions exist for this gene/variant, we check against them
        # For simplicity, if the table exists, we assume the variant passed that filter,
        # OR we compare text again. Let's compare text.
        if _text_match(extra_recs_trait, variant_trait):
            return True

    # 2. Main ACMG Table Check
    if not acmg_disease_string:
        # If ACMG table lists no specific disease, we assume the gene association is sufficient.
        return True
    
    return _text_match(acmg_disease_string, variant_trait)

def _text_match(target: str, query: str) -> bool:
    """Helper to tokenize and overlap two strings."""
    if not query or query in [".", "not provided", "not specified", "unknown"]:
        # If DB has no trait description, we cannot confirm match.
        # Strict Mode: This is a MISMATCH (Manual Review triggered if variant is strong).
        return False

    def get_tokens(s):
        # Lowercase
        s = s.lower()
        # Remove generic stopwords for medical genetics
        stopwords = {"type", "type", "disease", "syndrome", "deficiency", "disorder", "susceptibility", "to", "of", "and", "with", "hereditary", "familial"}
        # Replace non-alphanumeric with space
        s = re.sub(r'[^a-z0-9]', ' ', s)
        tokens = set(w for w in s.split() if len(w) > 2 and w not in stopwords)
        return tokens

    target_tokens = get_tokens(target)
    query_tokens = get_tokens(query)

    if not target_tokens: return True # Target was just stopwords?

    # INTERSECTION check
    # We require at least one specific noun match (e.g. "Hypertrophic", "Fabry", "Breast", "Ovarian")
    overlap = target_tokens.intersection(query_tokens)
    
    return len(overlap) > 0

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
                    m = re.search(r'Format:\s*(\S.*)', desc)
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
    return ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
            "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
            "Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE",
            "STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","CANONICAL","MANE_SELECT",
            "MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL",
            "UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS","miRNA",
            "HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF",
            "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF",
            "gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF",
            "MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS",
            "HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","AlphaMissense_score",
            "CADD_phred","REVEL_score"]

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
                for sep in ("/", ",", ";"):
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
        self.hgmd_vcf: Optional[pysam.VariantFile] = None
        
        # In-memory indexes for fast lookups
        self.clinvar_protein_index: Dict[str, List[dict]] = defaultdict(list)
        self.clinvar_splice_index: Dict[str, List[dict]] = defaultdict(list)
        
        # Internal DB Index: Nested dictionary [Chrom][Pos][Ref][Alt] -> Classification
        # Using a nested structure allows O(1) lookup speed for specific variants.
        self.internal_db_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

        # gnomAD handles
        self.gnomad_vcf_handles: Dict[str, pysam.VariantFile] = {}
        self.gnomad_tsv_handles: Dict[str, Dict[str,Any]] = {}
        
        # Initialize connections
        self._open_clinvar()
        self._open_hgmd()
        self._open_gnomad_versions(gnomad_v2, gnomad_v3, gnomad_v4)
        
        # Load Internal Database if path is provided
        if internal_db_path:
            self._load_internal_db(internal_db_path)

    def _open_clinvar(self):
        p = self.paths.get("CLINVAR_VCF")
        if p and os.path.exists(p) and HAVE_PYSAM:
            try:
                self.clinvar_vcf = pysam.VariantFile(p)
                logger.info("Opened ClinVar VCF: %s", p)
                try:
                    self._build_clinvar_indexes_from_vcf(p)
                except Exception as e:
                    logger.debug("ClinVar indexing skipped/failed: %s", e)
            except Exception as e:
                logger.warning("Failed to open ClinVar VCF: %s", e)

    def _open_hgmd(self):
        p = self.paths.get("HGMD_VCF")
        if p and os.path.exists(p) and HAVE_PYSAM:
            try:
                self.hgmd_vcf = pysam.VariantFile(p)
                logger.info("Opened HGMD VCF: %s", p)
            except Exception as e:
                logger.warning("Failed to open HGMD VCF: %s", e)

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
        Parses the Internal DB CSV.
        Assumes the file has been converted to hg38 coordinates.
        Expected columns: CHROM, POS, REF, ALT, anntype (classification).
        Handles multi-allelic entries (e.g., ALT="A, G") by splitting them.
        """
        if not os.path.exists(path):
            logger.warning(f"Internal DB path not found: {path}")
            return
        
        try:
            logger.info(f"Loading Internal DB from {path}...")
            # Read CSV, forcing string type to avoid float conversion errors on chromosomes/alleles
            df = pd.read_csv(path, dtype=str).fillna("")
            count = 0
            
            for _, r in df.iterrows():
                # Normalize chromosome (e.g., "1" -> "chr1") to match pipeline standard
                chrom = normalize_chrom(r.get('CHROM', ''))
                try: 
                    pos = int(r.get('POS', 0))
                except ValueError: 
                    continue
                    
                ref = r.get('REF', '').strip()
                
                # Clean ALT field: remove quotes, whitespace, and split multiple alleles
                raw_alt = r.get('ALT', '').replace('"', '').replace("'", "")
                alts = [x.strip() for x in raw_alt.split(',')]
                
                # Classification code (e.g., P, B, VUS)
                cls = r.get('anntype', '').strip()
                
                for alt in alts:
                    if not alt: continue
                    # Store mapping: Chrom -> Pos -> Ref -> Alt -> Classification
                    self.internal_db_index[chrom][pos][ref][alt] = cls
                    count += 1
            
            logger.info(f"Internal DB loaded successfully. Indexed {count} entries.")
            
        except Exception as e:
            logger.error(f"Failed to load internal DB: {e}", exc_info=True)

    def get_variant_db_evidence(self, vr: VariantRecord) -> List[Tuple[str,str]]:
        """
        Aggregates pathogenicity evidence from all databases.
        
        Priority Hierarchy for "High Quality" determination:
        1. ClinGen / ClinVar Expert Panel (Stars 3-4)
        2. Internal Database
        3. ClinVar Standard (Stars 0-2)
        4. HGMD Professional
        
        Returns:
            A list of tuples (Classification, Source) for the final report.
            Also updates `vr.is_hq_pathogenic_db` and `vr.is_hq_conflict_db`.
        """
        results: List[Tuple[str,str]] = []
        
        # Variables to store ClinVar state for hierarchy logic
        clinvar_found = False
        clinvar_stars = 0
        clinvar_sig_str = ""

        # =========================================================
        # 1. ClinVar & ClinGen (Parse Step)
        # =========================================================
        # We parse ClinVar first to determine if it meets the "ClinGen" criteria (Stars >= 3).
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
                            found_rec = rec
                            break
                    except: continue
                    if found_rec: break
                
                if found_rec:
                    clinvar_found = True
                    
                    # Extract Significance
                    clnsig = found_rec.info.get('CLNSIG') or ''
                    clinvar_sig_str = ";".join(clnsig) if isinstance(clnsig, (list, tuple)) else str(clnsig)
                    vr.clinvar_sig = clinvar_sig_str

                    # Extract Traits
                    clndn = found_rec.info.get('CLNDN') or found_rec.info.get('CLNDISDB') or ''
                    vr.clinvar_trait = str(clndn).replace("_", " ").replace("|", "; ")

                    # Extract Review Status & Stars
                    rev = found_rec.info.get('CLNREVSTAT') or ''
                    rev_str = str(rev).lower()
                    vr.clinvar_review = rev_str

                    if 'practice_guideline' in rev_str: clinvar_stars = 4
                    elif 'expert_panel' in rev_str: clinvar_stars = 3
                    elif 'criteria_provided' in rev_str and 'multiple_submitters' in rev_str and 'no_conflict' in rev_str: clinvar_stars = 2
                    elif 'criteria_provided' in rev_str: clinvar_stars = 1
                    
                    vr.clinvar_stars = str(clinvar_stars)

                    # --- HIERARCHY TIER 1: ClinGen / Expert Panel ---
                    # If Stars >= 3, this is the highest level of evidence.
                    if clinvar_stars >= 3:
                        results.append((clinvar_sig_str, "ClinGen/Expert"))
                        
                        # Determine Pathogenicity based on Expert opinion
                        if "pathogen" in clinvar_sig_str.lower() and "benign" not in clinvar_sig_str.lower():
                            vr.is_hq_pathogenic_db = True
            except Exception as e:
                logger.debug(f"ClinVar extraction failed: {e}")

        # =========================================================
        # 2. Internal Database (Hierachy Tier 2)
        # =========================================================
        # Checked after ClinGen. If ClinGen is absent, this becomes the primary driver.
        try:
            # Look for exact match in loaded index
            int_cls = self.internal_db_index.get(vr.chrom, {}).get(vr.pos, {}).get(vr.ref, {}).get(vr.alt)
            
            if int_cls:
                # Map short codes to full ACMG terms for readability
                mapped_cls = int_cls
                if int_cls == "P": mapped_cls = "Pathogenic"
                elif int_cls == "LP": mapped_cls = "Likely pathogenic"
                elif int_cls == "B": mapped_cls = "Benign"
                elif int_cls == "LB": mapped_cls = "Likely benign"
                elif int_cls == "VUS": mapped_cls = "VUS"
                
                # Add to results list so it appears in the output table
                results.append((mapped_cls, "Internal_DB"))
                
                # Trust Internal DB as a High Quality source (unless overridden by ClinGen above)
                # Note: We use additive logic here. If Internal DB says Pathogenic, we flag it.
                if "Pathogenic" in mapped_cls or "P" == int_cls:
                    vr.is_hq_pathogenic_db = True
                    
        except Exception as e:
            logger.debug(f"Internal DB lookup failed: {e}")

        # =========================================================
        # 3. ClinVar Standard (Hierarchy Tier 3)
        # =========================================================
        # If we found a ClinVar entry but it wasn't Expert Panel (Stars < 3), add it now.
        if clinvar_found and clinvar_stars < 3:
            results.append((clinvar_sig_str, "ClinVar"))
            
            # Standard logic: Stars >= 2 is considered High Quality (though lower priority than Internal)
            sig_lower = clinvar_sig_str.lower()
            if clinvar_stars >= 2:
                if 'pathogen' in sig_lower and 'benign' not in sig_lower:
                    vr.is_hq_pathogenic_db = True
                elif 'conflict' in sig_lower:
                    vr.is_hq_conflict_db = True

        # =========================================================
        # 4. HGMD Professional (Hierarchy Tier 4)
        # =========================================================
        # Checked last. Used if no better evidence exists, or to support existing evidence.
        if self.hgmd_vcf:
            try:
                c_clean = vr.chrom.replace("chr", "")
                # HGMD usually uses simple chromosome names (1, 2, X)
                for rec in self.hgmd_vcf.fetch(c_clean, vr.pos - 1, vr.pos):
                    if getattr(rec, 'ref', None) != vr.ref: continue
                    if not (rec.alts and vr.alt in rec.alts): continue
                    
                    hgmd_class = rec.info.get('CLASS') or ""
                    
                    # Map 'DM' (Disease Causing Mutation) to Pathogenic
                    mapped_hgmd = hgmd_class
                    if str(hgmd_class).upper() == 'DM':
                        mapped_hgmd = "Pathogenic"
                        # Mark as HQ evidence
                        vr.is_hq_pathogenic_db = True
                    
                    results.append((mapped_hgmd, "HGMD_Pro"))
            except Exception: 
                pass
            
        return results

    # ClinVar indexing (normalize gene keys)
    def _build_clinvar_indexes_from_vcf(self, clinvar_path: str):
        if not clinvar_path or not os.path.exists(clinvar_path):
            logger.warning("ClinVar path invalid: %s", clinvar_path)
            return

        try:
            logger.info("Loading ClinVar from %s ...", clinvar_path)
            count = 0
            
            with gzip.open(clinvar_path, "rt", encoding="utf-8", errors="replace") as fh:
                for ln in fh:
                    if ln.startswith("#"): continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 8: continue
                    chrom = parts[0]
                    try: pos = int(parts[1])
                    except: continue
                    info = parts[7]
                    
                    # Manual parsing of GENEINFO
                    # Search for substring "GENEINFO="
                    idx = info.find("GENEINFO=")
                    if idx == -1: continue
                    
                    # Cut string starting from value
                    rest = info[idx+9:]
                    # Value ends at next semicolon or end of line
                    end_idx = rest.find(";")
                    val = rest[:end_idx] if end_idx != -1 else rest
                    
                    # GENEINFO format: "Sym:ID|Sym2:ID" or just "Sym:ID"
                    if not val: continue
                    
                    # Take first gene part "Sym:ID"
                    first_part = val.split("|")[0]
                    if ":" not in first_part: continue
                    
                    gene_raw = first_part.split(":")[0]
                    gene = normalize_gene_name(gene_raw)
                    if not gene: continue

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
                        # Adjust based on conflicts/multis
                        if stars == 1 and 'multiple_submitters' in low and 'no_conflict' in low:
                            stars = 2
                    
                    # Protein
                    if "p." in info:
                        # Extract p. codes
                        matches = re.findall(r'(p\.[A-Za-z0-9\(\)\*]+)', info)
                        for mh in matches:
                            entry = {"hgvsp": mh, "clnsig": clnsig, "stars": stars, "chrom": chrom, "pos": pos}
                            self.clinvar_protein_index[gene].append(entry)
                    
                    # Splice
                    splice_entry = {"pos": pos, "clnsig": clnsig, "stars": stars, "spliceai": None}
                    self.clinvar_splice_index[gene].append(splice_entry)
                    
                    count += 1
            
            logger.info("ClinVar loaded. Index size (genes): %d. Total variants processed: %d", len(self.clinvar_splice_index), count)

        except Exception as e:
            logger.warning("ClinVar load error: %s", e)

    def find_protein_variant_by_hgvsp(self, gene: str, hgvsp: Optional[str], pos: Optional[int] = None, transcript: Optional[str] = None) -> List[Tuple[str,str]]:
        """
        Return high-quality protein-level matches from ClinVar/internal DB.
        Normalizes gene input and checks clinvar_protein_index entries.
        Returns a list of tuples (classification_string, source), e.g. ("Pathogenic","ClinVar")
        """
        results: List[Tuple[str,str]] = []
        try:
            norm_gene = normalize_gene_name(gene)
            entries = self.clinvar_protein_index.get(norm_gene) or []
            # If provided exact HGVSp, prefer exact protein match (and require stars>=2)
            if hgvsp:
                hgvsp_norm = str(hgvsp).strip().lower()
                for e in entries:
                    e_hgvsp = e.get("hgvsp")
                    if not e_hgvsp:
                        continue
                    if str(e_hgvsp).strip().lower() == hgvsp_norm and (e.get("stars", 0) >= 2):
                        clnsig = e.get("clnsig","")
                        if "pathogen" in str(clnsig).lower():
                            results.append(("Pathogenic","ClinVar"))
                        elif "likely" in str(clnsig).lower():
                            results.append(("Likely pathogenic","ClinVar"))
            # If pos provided, match by residue position (p.X123)
            if pos is not None:
                for e in entries:
                    p = e.get("hgvsp")
                    if not p:
                        continue
                    m = re.search(r'p\.[A-Za-z]{1,3}(\d+)', str(p))
                    if not m:
                        continue
                    try:
                        ppos = int(m.group(1))
                    except Exception:
                        continue
                    if ppos == int(pos) and (e.get("stars", 0) >= 2):
                        clnsig = e.get("clnsig","")
                        if "pathogen" in str(clnsig).lower():
                            results.append(("Pathogenic","ClinVar"))
                        elif "likely" in str(clnsig).lower():
                            results.append(("Likely pathogenic","ClinVar"))
        except Exception:
            logger.debug("find_protein_variant_by_hgvsp failed for gene=%s hgvsp=%s pos=%s", gene, hgvsp, pos, exc_info=True)
        return results

    def find_splice_equivalent_variants(self, gene: str, pos: int, radius: int = 6) -> List[Tuple[str,str,float]]:
        """
        Return splice-equivalent ClinVar matches near given genomic position.
        Returns list of tuples (classification, source, spliceai_score_or_0.0).
        Matches entries within +/- radius and with stars>=2.
        """
        out: List[Tuple[str,str,float]] = []
        try:
            norm_gene = normalize_gene_name(gene)
            entries = self.clinvar_splice_index.get(norm_gene) or []
            for e in entries:
                p = e.get("pos")
                if p is None:
                    continue
                try:
                    if abs(int(p) - int(pos)) <= int(radius) and (e.get("stars", 0) >= 2):
                        cln = e.get("clnsig","")
                        sa = e.get("spliceai") or 0.0
                        if cln and "pathogen" in str(cln).lower():
                            out.append(("Pathogenic","ClinVar", float(sa or 0.0)))
                        elif cln and "likely" in str(cln).lower():
                            out.append(("Likely pathogenic","ClinVar", float(sa or 0.0)))
                except Exception:
                    continue
        except Exception:
            logger.debug("find_splice_equivalent_variants failed for gene=%s pos=%s", gene, pos, exc_info=True)
        return out       

# ---------------------------
# Exons table & NMD logic (unchanged)
# ---------------------------
def load_exons_table(path: str):
    if not path or not os.path.exists(path):
        return None
    if not HAVE_PANDAS:
        logger.warning("pandas required to load exons table")
        return None
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
        for col in ("Genomic_start","Genomic_end","CDS_start","CDS_end","Exon_rank"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
        return df
    except Exception as e:
        logger.warning("Failed to load exons table: %s", e)
        return None

def compute_nmd_for_variant(vr: VariantRecord, exons_df: Optional[Any]) -> Tuple[bool, Optional[str]]:
    if exons_df is None or exons_df.empty:
        return (False, None)
    tx = vr.transcript or ""
    if not tx:
        return (False, None)
    tx_df = exons_df[exons_df["Transcript_ID"] == tx]
    if tx_df.empty:
        return (False, None)
    pos = int(vr.pos)
    containing = tx_df[(tx_df["Genomic_start"] <= pos) & (tx_df["Genomic_end"] >= pos)]
    if containing.empty:
        last_exon_rank = int(tx_df["Exon_rank"].max()) if "Exon_rank" in tx_df.columns else None
        return (False, "variant_not_in_exon_assume_NMD")
    exon_row = containing.iloc[0]
    if "Exon_rank" in tx_df.columns:
        try:
            cur_rank = int(exon_row["Exon_rank"])
            last_rank = int(tx_df["Exon_rank"].max())
            if cur_rank == last_rank:
                return (True, "in_last_exon")
            penult = tx_df[tx_df["Exon_rank"] == last_rank - 1]
            if not penult.empty and "CDS_end" in penult.columns and penult["CDS_end"].notnull().any():
                penult_cds_end = int(penult.iloc[0]["CDS_end"])
                if pos <= penult_cds_end - NMD_DISTANCE_TO_LAST_JUNCTION:
                    return (False, "triggers_nmd")
                else:
                    return (False, "escapes_nmd_close_to_last_junction")
            else:
                return (False, "not_last_exon_no_penult_cds")
        except Exception:
            return (False, None)
    else:
        sorted_exons = tx_df.sort_values("Genomic_start")
        last_exon = sorted_exons.iloc[-1]
        if exon_row.equals(last_exon):
            return (True, "in_last_exon")
        else:
            return (False, "triggers_nmd_assumed")
    return (False, None)

# ---------------------------
# DbNSFP cache (kept unchanged)
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
    def __init__(self, dbm: DatabaseManager, gene_rules: Dict[str, Dict[str,Any]], ttn_meta: Dict[str, IntervalTree] = None, spliceai_thr: float = 0.20):
        self.dbm = dbm
        self.gene_rules = gene_rules or {}
        self.ttn_meta = ttn_meta or {}
        self.spliceai_thr = spliceai_thr

    def _get_gene_rule(self, gene: str) -> Dict[str,Any]:
        return self.gene_rules.get((gene or "").upper(), {"reportable": False, "special_flags": {}, "raw_row": {}})

    def _ttn_overlap(self, chrom: str, pos: int) -> Optional[dict]:
        chromn = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
        if chromn not in self.ttn_meta:
            return None
        hits = self.ttn_meta[chromn].at(pos - 1)
        if not hits:
            return None
        return max(hits, key=lambda it: float(it.data.get("psi_meta", 0) or 0)).data

    def evaluate_variant(self, vr: VariantRecord, fasta_handle: Optional[Any] = None) -> Tuple[List[str], Dict[str,int], str, bool, List[str]]:
        pts: Dict[str,int] = {}
        assigned: List[str] = []
        manual_reasons: List[str] = []
        auto_class = "VUS"
        manual_review = False

        # --- 1. BA1 / BS1 Exception List (Hardcoded per flowchart) ---
        # Specific exception for BTD and HFE common pathogenic variants with high AF.
        is_ba1_exception = False
        
        # BTD: p.Asp444His (c.1330G>C)
        if "BTD" in vr.gene:
            if "Asp444His" in vr.hgvsp or "c.1330G>C" in vr.hgvsc:
                is_ba1_exception = True
        # HFE: p.Cys282Tyr (c.845G>A)
        elif "HFE" in vr.gene:
            if "Cys282Tyr" in vr.hgvsp or "c.845G>A" in vr.hgvsc:
                is_ba1_exception = True

        af = vr.gnomad_af if vr.gnomad_af is not None else 0.0
        
        # Apply population filtering ONLY if not an exception
        if not is_ba1_exception:
            if af >= THRESH["BA1_AF"]:
                pts["BA1"] = POINTS_MAP["BA1"]; assigned.append("BA1")
                auto_class = "Benign"
                return assigned, pts, auto_class, False, ["BA1: high population AF"]
            if af >= THRESH["BS1_AF"]:
                pts["BS1"] = POINTS_MAP["BS1"]; assigned.append("BS1")
                manual_reasons.append(f"BS1: AF {af} >= {THRESH['BS1_AF']}")

        # DB Analysis 
        db_hits = vr.db_hits or []
        db_supports_path = False
        for cls, src in db_hits:
            if "pathogen" in str(cls).lower():
                db_supports_path = True
                if "Benign" not in auto_class: auto_class = "Pathogenic"

        cons = (vr.consequence or "").lower()

        # --- 2. PVS1 / LoF Logic (Strict Decision Tree) ---
        if any(x in cons for x in LOF_CONSEQUENCES):
            rule = self._get_gene_rule(vr.gene)
            special = rule.get("special_flags", {}) if isinstance(rule, dict) else {}
            gene_mech = self.dbm.get_gene_mechanism(vr.gene) if hasattr(self.dbm, "get_gene_mechanism") else "Unknown"
            
            # Check if LoF is a valid mechanism for this gene
            lof_allowed = True
            if special.get("lof_not_reportable") or (gene_mech != "Valid" and not rule.get("mechanism")):
                lof_allowed = False
            
            if lof_allowed:
                # BRANCH: Start-Loss
                if "start_lost" in cons:
                    # Flowchart: Is there an alternative start codon?
                    found_alt_start = False
                    if fasta_handle and HAVE_PYSAM:
                         # Look for 'ATG' within +90bp of variant
                         try:
                             seq = fasta_handle.fetch(vr.chrom, vr.pos, vr.pos+90).upper()
                             if "ATG" in seq: found_alt_start = True
                         except: pass
                    
                    if found_alt_start:
                        # Alternative start exists -> Moderate strength (safe assumption without detailed upstream analysis)
                        pts["PVS1_Moderate"] = POINTS_MAP["PVS1_Moderate"]; assigned.append("PVS1_Moderate")
                    else:
                        # No alternative start -> VeryStrong
                        pts["PVS1_VeryStrong"] = POINTS_MAP["PVS1_VeryStrong"]; assigned.append("PVS1_VeryStrong")
                
                # BRANCH: Splice Donor/Acceptor (Canonical)
                elif any(x in cons for x in CANONICAL_SPLICE):
                    # Flowchart implies SpliceAI check for PVS1 here, usually default to VeryStrong
                    # unless rescue or weak splice impact. Defaulting to VeryStrong for canonical.
                    pts["PVS1_VeryStrong"] = POINTS_MAP["PVS1_VeryStrong"]; assigned.append("PVS1_VeryStrong")

                # BRANCH: Nonsense / Frameshift
                else:
                    # Flowchart: "Is it a critical region?"
                    in_critical = special.get("in_critical_region", False) 
                    
                    # NMD Check hierarchy: NMDEscPredictor -> VEP NMD -> 50bp rule
                    nmd_pred_raw = (getattr(vr, "nmdesc_predictor", "") or "").lower()
                    nmd_predicted = False
                    
                    if "escape" in nmd_pred_raw or "not_predicted" in nmd_pred_raw:
                        nmd_predicted = False
                    elif "nmd" in nmd_pred_raw or "predicted" in nmd_pred_raw:
                        nmd_predicted = True
                    else:
                        # Fallback to standard 50bp rule flag
                        nmd_res = (vr.nmd or "").lower()
                        if "triggers_nmd" in nmd_res: nmd_predicted = True
                        
                    if nmd_predicted:
                        pts["PVS1_VeryStrong"] = POINTS_MAP["PVS1_VeryStrong"]; assigned.append("PVS1_VeryStrong")
                    else:
                        # No NMD (escapes). Check Critical Region
                        if in_critical:
                            pts["PVS1_Strong"] = POINTS_MAP["PVS1_Strong"]; assigned.append("PVS1_Strong")
                        else:
                            # Not NMD, Not Critical -> PVS1_Moderate (Removing >10% protein check as usually true)
                            pts["PVS1_Moderate"] = POINTS_MAP["PVS1_Moderate"]; assigned.append("PVS1_Moderate")

        # --- 3. PS1 / PM5 (Protein Matches) ---
        # Logic remains largely standard, checking for exact or codon matches in DB
        pos_aa = parse_protein_pos(vr.hgvsp)
        prot_matches = self.dbm.find_protein_variant_by_hgvsp(vr.gene, vr.hgvsp, pos_aa, vr.transcript)
        
        if prot_matches:
            has_path = any(cls and 'pathogen' in str(cls).lower() for cls,_ in prot_matches if cls)
            if has_path:
                 pts["PS1_Strong"] = 4; assigned.append("PS1_Strong") 
        
        # --- 4. Computational Evidence (PP3)  ---
        # Flowchart: SpliceAI >= 0.2 (Supporting), REVEL >= 0.644 (Supp), >= 0.932 (Mod)
        
        # Check SpliceAI
        splice_score = float(vr.spliceai) if vr.spliceai is not None else 0.0
        if splice_score >= THRESH["SPLICEAI_SUPPORT"]: 
            # Apply only if PVS1 logic didn't already consume splicing evidence
            if "PVS1" not in str(assigned):
                if "PP3_Supporting" not in assigned:
                    pts["PP3_Supporting"] = POINTS_MAP["PP3_Supporting"]; assigned.append("PP3_Supporting")

        # Check REVEL
        if ("missense" in cons or "inframe" in cons) and vr.revel is not None:
            rv = float(vr.revel)
            if rv >= THRESH["REVEL_MODERATE"]: # 0.932
                # Remove lower score if added by SpliceAI to avoid double counting
                if "PP3_Supporting" in pts: 
                     del pts["PP3_Supporting"]
                     assigned.remove("PP3_Supporting")
                pts["PP3_Moderate"] = POINTS_MAP["PP3_Moderate"]; assigned.append("PP3_Moderate")
            elif rv >= THRESH["REVEL_SUPPORT"]: # 0.644
                if "PP3_Moderate" not in pts and "PP3_Supporting" not in pts:
                    pts["PP3_Supporting"] = POINTS_MAP["PP3_Supporting"]; assigned.append("PP3_Supporting")

        # --- 5. PS2 / PM6 (De Novo / Trio) ---
        pg = vr.proband_gt; fg = vr.father_gt; mg = vr.mother_gt
        dp = vr.proband_dp or 0; ab = vr.proband_ab or 0.0
        if _is_het(pg):
            if dp >= STRICT_MIN_DP and ab >= STRICT_MIN_AB and _is_wt(fg) and _is_wt(mg):
                 pts["PS2_VeryStrong"] = POINTS_MAP["PS2_VeryStrong"]; assigned.append("PS2_VeryStrong")
        
        # --- 6. PM2 (Population) ---
        # Should also check if absent in gnomAD v4 (assumed available in AF data)
        if vr.gnomad_af is not None:
            afv = vr.gnomad_af
            if afv == 0.0:
                pts["PM2_Moderate"] = POINTS_MAP["PM2_Moderate"]; assigned.append("PM2_Moderate")
            elif 0.0 < afv < THRESH["PM2_STRICT"]:
                pts["PM2_Supporting"] = POINTS_MAP["PM2_Supporting"]; assigned.append("PM2_Supporting")

        # Sum Points
        total_pts = 0
        for v in pts.values():
            total_pts += int(v)
            
        # Determine Auto Class
        if total_pts >= 10 or db_supports_path: auto_class = "Pathogenic"
        elif total_pts >= 6: auto_class = "Likely pathogenic"
        else: auto_class = "VUS"

        return list(assigned), pts, auto_class, manual_review, manual_reasons

    def apply_pm3_batch(self, candidates: List[VariantRecord]):
        """
        Rewritten PM3 logic (Recessive/AR Compound Het) per flowchart.
        Uses additive scoring for the second allele:
           - In trans with P/LP -> +1.0
           - Phase Unknown with P -> +0.5
           - Phase Unknown with LP -> +0.5 (Max 1.0 cumulative for unphased)
           - Homozygous -> +0.5
           - In trans with VUS -> +0.25
           
        Final summation:
           0.5 -> PM3_Supporting
           1.0 -> PM3_Moderate
           2.0 -> PM3_Strong
           >=4.0 -> PM3_VeryStrong
        """
        by_sample_gene = defaultdict(lambda: defaultdict(list))
        for v in candidates:
            by_sample_gene[v.sample][v.gene].append(v)
            
        for sample, genes in by_sample_gene.items():
            for gene, varlist in genes.items():
                if len(varlist) == 0: continue
                rule = self._get_gene_rule(gene)
                moi = str(rule.get("moi","")).upper()
                if "AD" in moi and "AR" not in moi: continue # Only apply to AR context

                # Pre-calculate "pathogenicity status" for each variant based on current points/DB
                for v in varlist:
                    curr_pts = sum(v.criteria_points.values())
                    v._temp_is_path = (curr_pts >= 10) or any('pathogen' in str(c[0]).lower() for c in v.db_hits)
                    v._temp_is_lp = (6 <= curr_pts < 10) or any('likely' in str(c[0]).lower() for c in v.db_hits)
                
                # Pairwise evaluation
                for v1 in varlist:
                    pm3_score = 0.0
                    
                    # 1. Homozygous check
                    if _is_hom(v1.proband_gt):
                        pm3_score += 0.5 
                    
                    # 2. Compound Heterozygote check
                    unphased_points_acc = 0.0 # Limit max 1.0 for unphased contributors
                    
                    for v2 in varlist:
                        if v1 == v2: continue 
                        
                        # Determine phase via Trio data
                        in_trans = False
                        phase_unknown = True
                        
                        try:
                            if _is_het(v1.father_gt) and _is_het(v2.mother_gt) and _is_wt(v1.mother_gt) and _is_wt(v2.father_gt):
                                in_trans = True; phase_unknown = False
                            elif _is_het(v2.father_gt) and _is_het(v1.mother_gt) and _is_wt(v2.mother_gt) and _is_wt(v1.father_gt):
                                in_trans = True; phase_unknown = False
                        except:
                            phase_unknown = True
                        
                        # Add points
                        if in_trans:
                            if v2._temp_is_path or v2._temp_is_lp:
                                pm3_score += 1.0 
                            else: # VUS in trans
                                pm3_score += 0.25 
                        elif phase_unknown:
                            # If phase is unknown, we award partial points
                            add_val = 0.0
                            if v2._temp_is_path: add_val = 0.5 
                            elif v2._temp_is_lp: add_val = 0.5 
                            
                            # Cap unphased points
                            if unphased_points_acc + add_val <= 1.0:
                                unphased_points_acc += add_val
                                pm3_score += add_val
                            elif unphased_points_acc < 1.0:
                                # (Optional) add remaining fraction to reach 1.0 cap
                                pass
                    
                    # Map calculated score to ACMG Criteria
                    new_crit = None
                    if pm3_score >= 4.0: new_crit = "PM3_VeryStrong" 
                    elif pm3_score >= 2.0: new_crit = "PM3_Strong"     
                    elif pm3_score >= 1.0: new_crit = "PM3_Moderate"   
                    elif pm3_score >= 0.5: new_crit = "PM3_Supporting" 
                    
                    if new_crit and new_crit not in v1.criteria_assigned:
                        v1.criteria_assigned.append(new_crit)
                        v1.criteria_points[new_crit] = POINTS_MAP.get(new_crit, 1) 
                        

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

def variant_row_with_disease(v: VariantRecord, rule: Dict[str,Any]) -> dict:
    acmg_disease = (rule.get("disease") or "") if isinstance(rule, dict) else ""
    
    # Format ClinVar detailed info
    # Example output: "Pathogenic | Stars: 2 | criteria provided, multiple submitters, no conflicts"
    clinvar_full_desc = ""
    if v.clinvar_sig:
        # Clean up review string (remove underscores)
        status_clean = v.clinvar_review.replace("_", " ")
        clinvar_full_desc = f"{v.clinvar_sig} | Stars: {v.clinvar_stars} | {status_clean}"
    else:
        clinvar_full_desc = "Not in ClinVar"

    return {
        "sample": v.sample, "chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt,
        "gene": v.gene_raw or v.gene, 
        
        # Combined Disease Info
        "acmg_disease_association": acmg_disease,
        "clinvar_trait": v.clinvar_trait,
        
        "clinvar_details": clinvar_full_desc,
        
        "clinvar_sig": v.clinvar_sig,
        "clinvar_stars": v.clinvar_stars,
        
        "transcript": v.transcript, "consequence": v.consequence,
        "hgvsc": v.hgvsc, "hgvsp": v.hgvsp, "exon": v.exon,
        
        # External DB data
        "gnomad_af": v.gnomad_af, 
        "gnomad_details": json.dumps(v.gnomad_details, ensure_ascii=False),
        
        "revel": v.revel, "alpha_missense": v.alpha_missense, "spliceai": v.spliceai, "cadd": v.cadd,
        "nmdesc_predictor": v.nmdesc_predictor, # Added NMD to output
        
        "proband_gt": v.proband_gt, "proband_dp": v.proband_dp,
        "db_hits": ";".join([f"{c}|{s}" for (c,s) in v.db_hits]) if v.db_hits else "",
        
        "criteria_assigned": ";".join(v.criteria_assigned),
        "criteria_points": json.dumps(v.criteria_points, ensure_ascii=False),
        "total_points": v.total_points,
        "automated_class": v.automated_class,
        "manual_reasons": ";".join(v.manual_reasons)
    }

def strict_triage_and_output(
    report_candidates: List[VariantRecord], 
    audit_candidates: List[VariantRecord],
    gene_rules: Dict[str, Dict[str,Any]], 
    recs_map: Dict[str, List[Dict]],
    outdir: str
):
    """
    Implements STRICT filtering, reporting, and full auditing.
    
    1. Triage: Sorts candidates into Auto, Manual, or Filtered buckets.
       - Uses 'report_candidates' (strict QC passed variants).
       
    2. Audit: Generates a master list of ALL processed candidates.
       - Uses 'audit_candidates' (includes low QC variants).
    """
    os.makedirs(outdir, exist_ok=True)
    
    # Output containers
    auto_rows = []
    manual_rows = []
    
    # Tracking sets for final audit flags (Chrom:Pos:Ref:Alt:Gene)
    ids_in_auto = set()
    ids_in_manual = set()

    def get_key(v):
        return f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.gene}"

    # Group variants by Sample -> Gene using REPORT CANDIDATES (Strict QC)
    by_sample_gene = defaultdict(lambda: defaultdict(list))
    for v in report_candidates:
        by_sample_gene[v.sample][v.gene].append(v)

    # --- MAIN TRIAGE PROCESSING LOOP ---
    for sample, genes in by_sample_gene.items():
        for gene, varlist in genes.items():
            rule = gene_rules.get(gene, {})
            moi = str(rule.get("moi","")).upper()
            special = rule.get("special_flags", {})
            acmg_disease = rule.get("disease", "")
            
            gene_recs = recs_map.get(gene, [])

            # 1. Determine "Candidate" Status for Logic Grouping
            scored_vars = []
            for v in varlist:
                # A. White List Check
                is_whitelisted = False
                
                # B. High Quality Database Entry
                is_hq = v.is_hq_pathogenic_db
                
                # C. Algorithmic Score
                is_algo_strong = (v.total_points >= 6)
                
                # Assign grouping status:
                if is_hq or is_algo_strong or is_whitelisted:
                    v._triage_group = "CANDIDATE" 
                elif v.automated_class == "VUS":
                    v._triage_group = "VUS"
                else:
                    v._triage_group = "BENIGN"
                
                # D. Disease Match Check
                match_main = check_disease_match(acmg_disease, v.clinvar_trait)
                match_recs = False
                for rec in gene_recs:
                    if check_disease_match(rec.get("disease"), v.clinvar_trait):
                        match_recs = True; break
                        
                v._disease_match = (match_main or match_recs)
                scored_vars.append(v)
            
            # Active variants only: ignore Benign/LB to ensure High Specificity
            p_lp_list = [v for v in scored_vars if v._triage_group == "CANDIDATE"]
            vus_list = [v for v in scored_vars if v._triage_group == "VUS"]
            
            if not p_lp_list and not vus_list: 
                continue

            # =======================
            #  INHERITANCE LOGIC
            # =======================

            # CASE A: RECESSIVE (AR)
            if ("AR" in moi) or ("RECESSIVE" in moi) or special.get("require_biallelic"):
                
                # 1. Pair of Two P/LP Variants
                if len(p_lp_list) >= 2:
                    all_match = all(x._disease_match for x in p_lp_list)
                    any_conflict = any(x.is_hq_conflict_db for x in p_lp_list)
                    
                    if not all_match:
                        for x in p_lp_list: 
                            r = variant_row_with_disease(x, rule)
                            r["manual_reasons"] = "AR Pair: Disease Mismatch"
                            manual_rows.append(r); ids_in_manual.add(get_key(x))
                    elif any_conflict:
                        for x in p_lp_list:
                            r = variant_row_with_disease(x, rule)
                            r["manual_reasons"] = "AR Pair: DB Conflict"
                            manual_rows.append(r); ids_in_manual.add(get_key(x))
                    else:
                        for x in p_lp_list:
                            auto_rows.append(variant_row_with_disease(x, rule))
                            ids_in_auto.add(get_key(x))
                
                # 2. One P/LP + One VUS (Compound Heterozygote)
                elif len(p_lp_list) == 1 and len(vus_list) >= 1:
                    v_plp = p_lp_list[0]
                    in_trans_confirmed = False
                    partners = []
                    
                    for v_vus in vus_list:
                        is_trans = False
                        if _is_het(v_plp.father_gt) and _is_het(v_vus.mother_gt): is_trans = True
                        elif _is_het(v_plp.mother_gt) and _is_het(v_vus.father_gt): is_trans = True
                        
                        if is_trans:
                            partners.append(v_vus)
                            in_trans_confirmed = True

                    if in_trans_confirmed:
                        r_s = variant_row_with_disease(v_plp, rule)
                        r_s["manual_reasons"] = "AR CompHet: P/LP + VUS in Trans"
                        manual_rows.append(r_s); ids_in_manual.add(get_key(v_plp))
                        
                        for p in partners:
                            r_v = variant_row_with_disease(p, rule)
                            r_v["manual_reasons"] = "AR CompHet: VUS partner in Trans"
                            manual_rows.append(r_v); ids_in_manual.add(get_key(p))

                # 3. Homozygous P/LP
                elif len(p_lp_list) == 1:
                    v = p_lp_list[0]
                    if _is_hom(v.proband_gt):
                        if v.is_hq_conflict_db:
                            r = variant_row_with_disease(v, rule)
                            r["manual_reasons"] = "AR Hom P/LP: DB Conflict"
                            manual_rows.append(r); ids_in_manual.add(get_key(v))
                        elif not v._disease_match:
                            r = variant_row_with_disease(v, rule)
                            r["manual_reasons"] = "AR Hom P/LP: Disease Mismatch"
                            manual_rows.append(r); ids_in_manual.add(get_key(v))
                        else:
                            auto_rows.append(variant_row_with_disease(v, rule))
                            ids_in_auto.add(get_key(v))

            # CASE B: DOMINANT (AD) or X-LINKED (XL)
            else:
                is_xl = ("X" in v.chrom or "X" in moi)
                
                for v in p_lp_list:
                    # 1. Sex Check for XL
                    if is_xl and (v.sample_sex == "Male"):
                        if _is_het(v.proband_gt):
                            r = variant_row_with_disease(v, rule)
                            r["manual_reasons"] = "XL Male: Heterozygous P/LP call (Check Data)"
                            manual_rows.append(r); ids_in_manual.add(get_key(v))
                            continue

                    # 2. Conflict Check
                    if v.is_hq_conflict_db:
                        r = variant_row_with_disease(v, rule)
                        r["manual_reasons"] = "AD/XL P/LP: DB Conflict"
                        manual_rows.append(r); ids_in_manual.add(get_key(v))
                        continue

                    # 3. Disease Match Check
                    if not v._disease_match:
                        r = variant_row_with_disease(v, rule)
                        r["manual_reasons"] = f"AD/XL P/LP: Disease Mismatch ({v.clinvar_trait})"
                        manual_rows.append(r); ids_in_manual.add(get_key(v))
                        continue

                    # 4. Clean Auto
                    auto_rows.append(variant_row_with_disease(v, rule))
                    ids_in_auto.add(get_key(v))

    # --- 2. AUDIT GENERATION (Using audit_candidates - includes low QC) ---
    all_audit_rows = []
    for v in audit_candidates:
        rule = gene_rules.get((v.gene or "").upper(), {})
        row = variant_row_with_disease(v, rule)
        k = get_key(v)
        
        row['in_auto_report'] = (k in ids_in_auto)
        row['in_manual_review'] = (k in ids_in_manual)
        
        all_audit_rows.append(row)

    # --- 3. DEDUPLICATION AND WRITING ---
    
    all_rows_c = collapse_output_rows(all_audit_rows)
    auto_rows_c = collapse_output_rows(auto_rows)
    manual_rows_c = collapse_output_rows(manual_rows)

    if HAVE_PANDAS:
        pd.DataFrame(all_rows_c).to_csv(os.path.join(outdir, "all_candidates.csv"), index=False)
        pd.DataFrame(auto_rows_c).to_csv(os.path.join(outdir, "auto_conclusions.csv"), index=False)
        pd.DataFrame(manual_rows_c).to_csv(os.path.join(outdir, "manual_review_list.csv"), index=False)
    
    return os.path.join(outdir, "all_candidates.csv"), \
           os.path.join(outdir, "auto_conclusions.csv"), \
           os.path.join(outdir, "manual_review_list.csv"), \
           len(auto_rows_c), len(manual_rows_c)

# ---------------------------
# Build variant record from VCF record (modified to use normalized gene matching)
# ---------------------------
def build_variant_record_from_rec(rec, csq_header: List[str], acmg_genes_normalized: set, dad_vf: Optional[pysam.VariantFile], mom_vf: Optional[pysam.VariantFile], exons_df: Optional[Any], proband_sample: Optional[str]=None) -> Optional[VariantRecord]:
    try:
        # 1. Basic properties (Universal)
        chrom = getattr(rec, "CHROM", getattr(rec, "chrom", None)) or str(getattr(rec,"contig",""))
        pos = int(getattr(rec, "POS", getattr(rec, "pos", 0)))
        ref = getattr(rec, "REF", getattr(rec, "ref", ""))
        alts = getattr(rec, "ALT", getattr(rec, "alt", None)) or getattr(rec, "alts", None)
        alt = alts[0] if alts and len(alts) > 0 else ""
        if not alt: return None
    except Exception: return None

    # 2. Extract INFO (Universal)
    info = {}
    if hasattr(rec, "INFO"): # cyvcf2
        info = dict(rec.INFO)
    elif hasattr(rec, "info"): # pysam
        try: info = dict(rec.info.items())
        except: info = dict(rec.info)
    info_l = {k.lower(): v for k, v in info.items()}

    # 3. Find Best Transcript/Annotation
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

    # Fallback to INFO gene
    if not best_ann:
        gene_val = None
        for k in ("gene","symbol","hgnc","geneinfo"):
            if k in info_l:
                v = info_l[k]
                gene_val = str(v[0] if isinstance(v, (list,tuple)) else v).split(":")[0]
                break
        if gene_val:
            ng = normalize_gene_name(gene_val)
            if ng in acmg_genes_normalized:
                best_ann = {"SYMBOL": gene_val, "_NORM_GENE": ng, "_RAW_GENE": gene_val}
    if not best_ann: return None

    norm_gene = best_ann.get("_NORM_GENE")
    gene_raw = best_ann.get("_RAW_GENE")

    # 4. Extract Genotype Data (HANDLING CYVCF2 vs PYSAM)
    sample_name = ""
    prob_gt = "./."
    prob_dp = None
    prob_ad = None
    
    # Check object type by presence of specific methods
    is_cyvcf2 = hasattr(rec, "genotypes") and hasattr(rec, "format")
    
    if is_cyvcf2:
        # --- CYVCF2 Handling ---
        if hasattr(rec, "samples") and rec.samples:
            sample_name = proband_sample if (proband_sample and proband_sample in rec.samples) else rec.samples[0]
            
        # GT from rec.genotypes (list of [0, 1, bool])
        if rec.genotypes:
            # We take the first sample (index 0)
            gt_data = rec.genotypes[0] 
            # gt_data is e.g. [0, 1, False] for 0/1 phased=False
            # cyvcf2 uses -1 for dot
            if len(gt_data) >= 2:
                a1 = gt_data[0]; a2 = gt_data[1]
                s1 = str(a1) if a1 != -1 else "."
                s2 = str(a2) if a2 != -1 else "."
                sep = "|" if (len(gt_data)>2 and gt_data[2]) else "/"
                prob_gt = f"{s1}{sep}{s2}"
                
        # DP from rec.format('DP') -> numpy array
        try:
            dp_arr = rec.format('DP')
            if dp_arr is not None:
                prob_dp = int(dp_arr[0][0]) # Sample 0
        except: pass
        
        # AD from rec.format('AD')
        try:
            ad_arr = rec.format('AD')
            if ad_arr is not None:
                # Array is (n_samples, n_alleles) ?? Usually it's raw values
                # Actually cyvcf2 format('AD') returns (n_samples, X) usually REF,ALT
                row = ad_arr[0]
                if len(row) >= 2:
                    prob_ad = (int(row[0]), int(row[1]))
        except: pass

    else:
        # --- PYSAM Handling ---
        s_obj = getattr(rec, "samples", None)
        if s_obj:
            if proband_sample and proband_sample in s_obj:
                sample_name = proband_sample
            elif len(s_obj) > 0:
                sample_name = list(s_obj.keys())[0]
            
            if sample_name:
                s_dat = s_obj[sample_name]
                # GT
                if "GT" in s_dat:
                    g = s_dat["GT"]
                    prob_gt = "/".join("." if x is None else str(x) for x in g)
                # DP
                if "DP" in s_dat: prob_dp = s_dat["DP"]
                # AD
                if "AD" in s_dat:
                    a = s_dat["AD"]
                    if isinstance(a, (list, tuple)) and len(a)>=2: prob_ad = (a[0], a[1])

    # Fallback for DP/GT from INFO if empty
    if prob_dp is None and "dp" in info_l:
        try: prob_dp = int(info_l["dp"])
        except: pass

    # Calculate AB
    prob_ab = None
    if prob_ad and prob_dp and prob_dp > 0:
        prob_ab = prob_ad[1] / prob_dp
    elif prob_ad:
        t = prob_ad[0]+prob_ad[1]
        if t>0: prob_ab = prob_ad[1]/t
        if prob_dp is None: prob_dp = t

    # Helper to extract values from INFO or CSQ fields with robust parsing
    def get_val_anywhere(keys, return_type=float):
        
        def _parse_raw(val):
            # If the value is a list or tuple (common in pysam/cyvcf2), take the first element
            if isinstance(val, (list, tuple)):
                val = val[0]
            
            s = str(val)
            if not s or s == ".": 
                return None
            
            # If '&' is present, VEP returned multiple values (one per transcript).
            # Example: "0.178&0.289&0.331"
            if "&" in s:
                parts = s.split("&")
                nums = []
                for p in parts:
                    try: 
                        nums.append(float(p))
                    except: 
                        pass
                
                # For numerical scores (float), we take the MAXIMUM value.
                # This ensures we capture the most severe prediction across transcripts.
                if return_type == float and nums:
                    return max(nums)
                elif return_type == str and nums:
                    return str(max(nums))
                
                # If parsing failed or type is not float, return the raw string
                return s

            # Standard single value conversion
            try: 
                return return_type(s)
            except: 
                return None

        # 1. Search in VCF INFO fields (case-insensitive)
        for k in keys:
            if k.lower() in info_l:
                res = _parse_raw(info_l[k.lower()])
                if res is not None: 
                    return res
                
        # 2. Search in VEP CSQ annotations
        for k in keys:
            if k in best_ann and best_ann[k]:
                res = _parse_raw(best_ann[k])
                if res is not None: 
                    return res
                
        return None
    
    # 1. NMD: Extract Nonsense Mediated Decay prediction.
    # We check "NMD" (standard VEP plugin output) and "NMDEscPredictor" (older tools).
    nmdesc = get_val_anywhere(["NMDEscPredictor", "nmdesc_predictor", "NMD"], return_type=str)

    # 2. Pathogenicity Scores: REVEL, AlphaMissense, CADD.
    # We check multiple casing variations to ensure compatibility with different dbNSFP versions.
    revel = get_val_anywhere(["REVEL_score", "revel", "REVEL", "revel_score"])
    alpha = get_val_anywhere(["AlphaMissense_score", "AlphaMissense", "alphamissense_score"])
    cadd = get_val_anywhere(["CADD_phred", "CADD", "cadd_phred"])

    # 3. SpliceAI: Aggregate the maximum score from 4 specific delta scores.
    # The VEP plugin outputs separate scores for Acceptor/Donor Gain/Loss (AG, AL, DG, DL).
    # We take the maximum value to represent the most severe potential splicing impact.
    splice_keys = [
        "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", 
        "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
        "SpliceAI_pred_DS_max", "SpliceAI" # Fallback keys
    ]
    
    splice_vals = []
    for key in splice_keys:
        val = get_val_anywhere([key])
        if val is not None:
            splice_vals.append(val)
            
    # If any SpliceAI scores were found, use the highest one.
    spliceai = max(splice_vals) if splice_vals else None

    # 4. gnomAD AF: Global Allele Frequency extraction.
    # Logic: Try extracting from top-level INFO fields first. 
    # If not found, fall back to extracting from CSQ/ANN annotations.
    gnomad_af = None; gdet = {}
    _, gm = extract_gnomad_af_from_info(info_l)
    if gm.get("max_af", 0) > 0: 
        gnomad_af = gm["max_af"]; gdet = gm
    else:
        _, gm2 = extract_gnomad_af_from_info(best_ann)
        if gm2.get("max_af", 0) > 0: 
            gnomad_af = gm2["max_af"]; gdet = gm2

    vr = VariantRecord(
        chrom=normalize_chrom(chrom), pos=int(pos), ref=ref, alt=str(alt),
        sample=sample_name,
        gene=norm_gene, gene_raw=gene_raw,
        transcript=best_ann.get("Feature") or "", consequence=best_ann.get("Consequence") or "",
        hgvsc=best_ann.get("HGVSc") or "", hgvsp=best_ann.get("HGVSp") or "",
        exon=best_ann.get("EXON") or "",
        gnomad_af=gnomad_af, gnomad_details=gdet,
        revel=revel, alpha_missense=alpha, spliceai=spliceai, cadd=cadd, nmdesc_predictor=nmdesc,
        proband_gt=prob_gt, proband_dp=prob_dp, proband_ad=prob_ad, proband_ab=prob_ab
    )
    vr._raw_ann = best_ann
    return vr

# ---------------------------
# Pre-analysis - quality and biological sex control
# ---------------------------

def perform_pre_analysis_qc(vcf_path: str, proband_id: str) -> Dict[str, Any]:
    """
    Performs pre-analysis quality control.
    Specifically checks biological sex based on X-chromosome heterozygosity 
    to determine applicable criteria (PM2/Hemizygosity) later in the pipeline.
    """
    warnings = []
    inferred_sex = "Unknown"
    
    # 1. Sex Inference via X-chrom heterozygosity ratio
    if HAVE_PYSAM:
        try:
            vf = pysam.VariantFile(vcf_path)
            chroms = list(vf.header.contigs)
            # Find X chromosome designator
            chr_x = next((c for c in chroms if c.lower() in ('chrx', 'x')), None)
            
            if chr_x:
                het_count = 0
                total_x = 0
                # Iterate over a subset of variants on ChrX to save time
                for rec in vf.fetch(chr_x):
                    if len(rec.samples) == 0: break
                    sample = proband_id if proband_id in rec.samples else rec.samples[0]
                    gt = rec.samples[sample]['GT']
                    # Check for heterozygous genotype (e.g. 0/1)
                    if len(gt) >= 2 and gt[0] != gt[1]:
                        het_count += 1
                    total_x += 1
                    if total_x > 1000: break 
                
                ratio = het_count / total_x if total_x > 0 else 0
                
                # High het ratio -> Likely Female
                if ratio > 0.2:
                    inferred_sex = "Female"
                # Very low het ratio -> Likely Male (Hemizygous)
                elif ratio < 0.05 and total_x > 50:
                    inferred_sex = "Male"
                
                logger.info(f"QC Sex Check: Inferred {inferred_sex} (Het Ratio on X: {ratio:.2f})")
        except Exception as e:
            logger.warning(f"QC Sex Check failed: {e}")

    # 2. Relationship/Contamination checks would go here (requires external tools/logic)
    
    return {"warnings": warnings, "sex": inferred_sex}


# ---------------------------
# Main pipeline (glue) - includes added minimal tests/logging before parsing
# ---------------------------
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
                ttn_meta_csv: Optional[str],
                exons_file: Optional[str],
                aggressive: bool,
                gene_bed: Optional[str] = None,
                cadd_tabix: Optional[str] = None,
                clinvar_tabix: Optional[str] = None,
                dbnsfp_variants_dir: Optional[str] = None,
                dbnsfp_cache_sqlite: Optional[str] = None,
                proband_sample: Optional[str] = None,
                strict_qc: bool = True,
                gnomad_v2: Optional[str] = None,
                gnomad_v3: Optional[str] = None,
                acmg_recs: Optional[str] = None,
                internal_db_path: Optional[str] = None,
                gnomad_v4: Optional[str] = None):
    os.makedirs(outdir, exist_ok=True)

    # 1. PREPARE INPUT VCF (Compress & Index FIRST)
    # This must be done before QC because pysam.fetch requires an index
    input_vcf_for_tools = proband_vcf
    if not proband_vcf.endswith(".gz") and not proband_vcf.endswith(".bgz"):
        logger.info("Input VCF is plain text. Compressing and indexing for tool compatibility...")
        compressed_vcf = os.path.join(outdir, "input_prepped.vcf.gz")
        
        # Use bgzip if available (preferred for tabix), otherwise python gzip
        if shutil.which("bgzip"):
            subprocess.run(f"bgzip -c '{proband_vcf}' > '{compressed_vcf}'", shell=True, check=True)
        else:
            with open(proband_vcf, 'rb') as f_in, gzip.open(compressed_vcf, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Index the file
        subprocess.run(["tabix", "-p", "vcf", compressed_vcf], check=False)
        input_vcf_for_tools = compressed_vcf
    elif not os.path.exists(proband_vcf + ".tbi") and not os.path.exists(proband_vcf + ".csi"):
        # If compressed but missing index
        try: subprocess.run(["tabix", "-p", "vcf", proband_vcf], check=False)
        except: pass
    
    working_vcf = input_vcf_for_tools
    annotated_vcf = working_vcf

    # 2. PERFORM QC (Using the indexed file)
    # Load Extra Recommendations
    recs_map = load_recommendations_table(acmg_recs) if acmg_recs else {}
    
    # Pass input_vcf_for_tools (the .gz file), not raw proband_vcf
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
        logger.info("Filtered VCF written to %s", filtered_vcf)
    else:
        if gene_bed:
            logger.warning("gene-bed specified but not found: %s", gene_bed)

    # run VEP with dbNSFP 
    if run_vep_dbnsfp:
        if not dbnsfp_path or not os.path.exists(dbnsfp_path):
            raise FileNotFoundError(f"dbNSFP file not found: {dbnsfp_path}")
        
        out_v = os.path.join(outdir, os.path.basename(working_vcf).replace(".vcf","").replace(".gz","") + ".dbnsfp.vep.vcf.gz")
        
        # Smartly detect available fields to prevent VEP errors
        logger.info("Inspecting dbNSFP header to validate requested fields...")
        dbnsfp_plugin_string = choose_dbnsfp_fields(dbnsfp_path, dbnsfp_fields)
        
        # Run VEP using the robust wrapper
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

    # [FASTA STEP] Open FASTA handle for Start-Loss verification
    fasta_handle = None
    if fasta and os.path.exists(fasta) and HAVE_PYSAM:
        try:
            fasta_handle = pysam.FastaFile(fasta)
        except Exception as e:
            logger.warning(f"Could not open Fasta file: {e}")

    # load ACMG rules (use normalized gene keys)
    gene_rules = load_acmg_table(acmg_table) if HAVE_PANDAS else {}
    acmg_genes_normalized = set(g for g in gene_rules.keys() if gene_rules.get(g,{}).get("reportable", False))
    logger.info("Reportable genes (normalized): %d", len(acmg_genes_normalized))
    logger.debug("Sample of normalized genes: %s", list(sorted(acmg_genes_normalized))[:20])

    # init DB manager
    dbm = DatabaseManager(db_paths, gnomad_v2=gnomad_v2, gnomad_v3=gnomad_v3, gnomad_v4=gnomad_v4, internal_db_path=internal_db_path)

    exons_df = load_exons_table(exons_file) if exons_file else None

    csq_header = parse_csq_header_from_vcf(annotated_vcf)
    if not csq_header:
        logger.warning("Could not parse CSQ/ANN header from VCF; VEP annotations may be non-standard")

    # Minimal sanity checks/logging: count variants in filtered VCF
    try:
        if shutil.which("bcftools"):
            cnt = int(subprocess.run(["bcftools","view","-H", annotated_vcf], capture_output=True, text=True).stdout.count("\n"))
            logger.info("Sanity check: annotated VCF has approximately %d variant lines (header excluded)", cnt)
        else:
            logger.info("Sanity check: bcftools not available to count VCF lines")
    except Exception:
        logger.debug("Sanity check for VCF lines failed", exc_info=True)

    # VCF parsing
    vcf_reader = None
    if HAVE_CYVCF2:
        vcf_reader = VCF(annotated_vcf)
    else:
        if not HAVE_PYSAM:
            raise RuntimeError("Neither cyvcf2 nor pysam available for VCF parsing")
    dad_vf = pysam.VariantFile(father_vcf) if father_vcf and os.path.exists(father_vcf) and HAVE_PYSAM else None
    mom_vf = pysam.VariantFile(mother_vcf) if mother_vcf and os.path.exists(mother_vcf) and HAVE_PYSAM else None

    candidates: List[VariantRecord] = []
    logger.info("Parsing VCF and collecting candidate variants (ACMG genes)...")
    if vcf_reader:
        for rec in vcf_reader:
            vr = build_variant_record_from_rec(rec, csq_header, acmg_genes_normalized, dad_vf, mom_vf, exons_df, proband_sample=proband_sample)
            if not vr:
                continue
            # Assign inferred sex to record
            vr.sample_sex = inferred_sex
            
            if (vr.gnomad_af is None or vr.gnomad_af == 0.0) and hasattr(dbm, "annotate_gnomad_for_variant"):
                try:
                    dbm.annotate_gnomad_for_variant(vr)
                except Exception:
                    logger.debug("gnomad annotation failed for %s:%d", vr.chrom, vr.pos)
            try:
                vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
            except Exception:
                pass
            candidates.append(vr)
    else:
        try:
            vf = pysam.VariantFile(annotated_vcf)
            for rec in vf.fetch():
                info = dict(rec.info.items())
                csq = info.get("CSQ") or info.get("ANN")
                if not csq:
                    continue
                csq_str = csq[0] if isinstance(csq,(list,tuple)) else str(csq)
                best_ann = parse_vep_csq_line(csq_str.split(",")[0], csq_header) if csq_header else {}
                gene_raw = best_ann.get('SYMBOL') or best_ann.get('Gene') or ""
                norm_gene = normalize_gene_name(gene_raw)
                if norm_gene not in acmg_genes_normalized:
                    continue
                vr = VariantRecord(chrom=normalize_chrom(rec.chrom), pos=int(rec.pos), ref=rec.ref, alt=','.join(rec.alts) if rec.alts else "", sample=(list(rec.samples.keys())[0] if rec.samples else ""))
                vr.gene = norm_gene; vr.gene_raw = gene_raw; vr.hgvsp = best_ann.get('HGVSp','') or ''
                if (vr.gnomad_af is None or vr.gnomad_af == 0.0) and hasattr(dbm, "annotate_gnomad_for_variant"):
                    dbm.annotate_gnomad_for_variant(vr)
                try:
                    vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
                except Exception:
                    pass
                candidates.append(vr)
        except Exception as e:
            logger.error("Failed to parse VCF with pysam fallback: %s", e)
            raise

    logger.info("Collected %d candidate variants (pre-QC)", len(candidates))
    
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
            # Prefer higher total_points
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
            def score(rec: VariantRecord) -> float:
                r = getattr(rec, "revel", None) or 0.0
                c = getattr(rec, "cadd", None) or 0.0
                return float(r) + float(c) / 25.0
            if score(v) > score(existing):
                keymap[key] = v
                continue
            # Keep existing otherwise
        return list(keymap.values())

    candidates = deduplicate_candidates(candidates)
    logger.info("After deduplication: %d candidates", len(candidates))

    # Force sample name override
    if proband_sample:
        for v in candidates:
            v.sample = proband_sample

    # --- 1. SCORING ALL CANDIDATES FIRST ---
    
    # Pass fasta_handle to evaluation
    engine = ACMGEngine(dbm, gene_rules, {}, THRESH["SPLICEAI"]) 
    
    for v in candidates:
        assigned, pts, auto_class, manual_flag, manual_reasons = engine.evaluate_variant(v, fasta_handle=fasta_handle)
        v.criteria_assigned = assigned
        v.criteria_points = pts
        v.automated_class = auto_class
        v.manual_review = manual_flag
        v.manual_reasons = manual_reasons
        total = 0
        for val in v.criteria_points.values():
            try: total += int(val)
            except Exception: pass
        v.total_points = int(total)
    
    if fasta_handle: fasta_handle.close()

    # PM3 batch
    engine.apply_pm3_batch(candidates) 

    # Recompute totals and finalize (Standardized Strings)
    for v in candidates:
        total = 0
        for val in v.criteria_points.values():
            try: total += int(val)
            except Exception: pass
        v.total_points = int(total)
        
        if v.total_points >= 10: v.automated_class = "Pathogenic"
        elif v.total_points >= 6: v.automated_class = "Likely Pathogenic"
        elif v.total_points <= -6: v.automated_class = "Benign"
        elif v.total_points <= -1: v.automated_class = "Likely Benign"
        else:
            v.automated_class = "VUS"

    # --- 2. NOW APPLY STRICT QC (Filtering for Report) ---
    filtered_for_report = [] 
    dropped_count = 0
    
    for v in candidates:
        if not STRICT_QC_ENABLED or not strict_qc:
            filtered_for_report.append(v); continue
        if v.db_hits:
            filtered_for_report.append(v); continue
            
        dp = v.proband_dp or 0; ab = v.proband_ab or 0.0; gt = v.proband_gt or "./."
        
        if _is_het(gt):
            if dp >= STRICT_MIN_DP and ab >= STRICT_MIN_AB:
                filtered_for_report.append(v)
            else:
                dropped_count += 1
        elif _is_hom(gt):
            if dp >= STRICT_MIN_DP:
                filtered_for_report.append(v)
            else:
                dropped_count += 1
        else:
            dropped_count += 1
            
    logger.info("After strict QC filtering: %d kept for Report, %d low quality (audit only)", len(filtered_for_report), dropped_count)

    # --- 3. FINAL OUTPUT ---
    # 1. filtered_for_report -> used for Auto/Manual logic
    # 2. candidates -> used for All Candidates table
    all_path, auto_path, manual_path, cnt_auto, cnt_manual = strict_triage_and_output(
        report_candidates=filtered_for_report,  
        audit_candidates=candidates,            
        gene_rules=gene_rules,
        recs_map=recs_map, 
        outdir=outdir
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

# --- INSERT AFTER load_acmg_table ---

def load_recommendations_table(path: str) -> Dict[str, Dict[str,Any]]:
    """
    Loads the secondary CSV with specific variant/gene recommendations.
    Expected columns: Gene, Variant (optional), Disease, Action.
    Returns dictionary keyed by Normalized Gene.
    """
    recs = defaultdict(list)
    if not path or not os.path.exists(path):
        return recs
    
    try:
        df = pd.read_csv(path, sep=None, engine='python', dtype=str).fillna("")
        for _, r in df.iterrows():
            g = normalize_gene_name(r.get("Gene Symbol", "") or r.get("Gene", ""))
            if not g: continue
            rec_data = {
                "variant_id": r.get("Variant", "").strip(), # specific variant whitelist
                "disease": r.get("Disease", "").strip(),
                "inclusion": r.get("Inclusion", "True").lower() == "true"
            }
            recs[g].append(rec_data)
        logger.info(f"Loaded {len(recs)} specific recommendations from {path}")
    except Exception as e:
        logger.warning(f"Failed to load recommendations table: {e}")
    
    return recs

# ---------------------------
# CLI
# ---------------------------
def build_parser():
    p = argparse.ArgumentParser(prog="acmg_sf_pipeline_final.py", description="ACMG SF classifier - final ready-to-run")
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
    return p

def run_batch_mode(args, db_paths):
    if not HAVE_PANDAS:
        logger.error("Batch mode requires pandas.")
        sys.exit(1)

    logger.info("STARTING BATCH MODE ANALYSIS")
    
    search_pattern = os.path.join(args.batch_input_dir, "**", "*_proband.vcf")
    proband_files = glob.glob(search_pattern, recursive=True)
    
    if not proband_files:
        logger.error("No *_proband.vcf files found in %s", args.batch_input_dir)
        sys.exit(1)

    logger.info("Found %d samples.", len(proband_files))
    
    # Lists to hold dataframes for final concatenation
    auto_results_list = []
    manual_results_list = []
    all_candidates_list = [] 

    for p_vcf in proband_files:
        dirname = os.path.dirname(p_vcf)
        filename = os.path.basename(p_vcf)
        
        dna_id = filename.split('_')[0] 
        map_id = os.path.basename(dirname)

        logger.info(f">>> Processing: {dna_id} (Family: {map_id})")

        siblings = glob.glob(os.path.join(dirname, "*.vcf"))
        f_vcf = next((s for s in siblings if "_father" in s), None)
        m_vcf = next((s for s in siblings if "_mother" in s), None)

        # Output folder for this specific sample
        s_out = os.path.join(args.outdir, "individual_results", dna_id)
        
        try:
            dbnsfp_fields = [x.strip() for x in args.dbnsfp_fields.split(",")] if args.dbnsfp else []
            
            # Run pipeline (creates 3 CSVs in s_out)
            process_vcf(
                proband_vcf=p_vcf, father_vcf=f_vcf, mother_vcf=m_vcf, outdir=s_out,
                run_vep_dbnsfp=args.run_vep_dbnsfp, dbnsfp_path=args.dbnsfp, dbnsfp_fields=dbnsfp_fields,
                vep_cmd=args.vep, vep_cache=args.vep_cache, fasta=args.fasta, vep_extra=args.vep_extra,
                acmg_table=args.acmg_table, db_paths=db_paths, ttn_meta_csv=args.ttn_meta,
                exons_file=args.exons_file, aggressive=args.aggressive, gene_bed=args.gene_bed,
                cadd_tabix=args.cadd_tabix, clinvar_tabix=args.clinvar_tabix,
                dbnsfp_variants_dir=args.dbnsfp_variants_dir, dbnsfp_cache_sqlite=args.dbnsfp_cache_sqlite,
                proband_sample=dna_id, strict_qc=not args.no_strict_qc,
                gnomad_v2=args.gnomad_v2, gnomad_v3=args.gnomad_v3, gnomad_v4=args.gnomad_v4,
                acmg_recs=args.acmg_recs
            )

            # Paths to the newly created files
            path_auto = os.path.join(s_out, "auto_conclusions.csv")
            path_manual = os.path.join(s_out, "manual_review_list.csv")
            path_all = os.path.join(s_out, "all_candidates.csv")

            # --- Helpers to identify unique rows ---
            def make_keys(df):
                if df is None or df.empty: return set()
                # Create a set of tuples: (chrom, pos, ref, alt, gene)
                return set(zip(df['chrom'].astype(str), df['pos'].astype(str), 
                               df['ref'].astype(str), df['alt'].astype(str), 
                               df['gene'].astype(str)))

            df_auto = pd.read_csv(path_auto) if (os.path.exists(path_auto) and os.path.getsize(path_auto)>1) else pd.DataFrame()
            df_manual = pd.read_csv(path_manual) if (os.path.exists(path_manual) and os.path.getsize(path_manual)>1) else pd.DataFrame()
            df_all = pd.read_csv(path_all) if (os.path.exists(path_all) and os.path.getsize(path_all)>1) else pd.DataFrame()

            # Aggregate Auto
            if not df_auto.empty:
                df_auto['sample'] = dna_id
                df_auto.insert(1, 'Family_ID', map_id)
                auto_results_list.append(df_auto)

            # Aggregate Manual
            if not df_manual.empty:
                df_manual['sample'] = dna_id
                df_manual.insert(1, 'Family_ID', map_id)
                manual_results_list.append(df_manual)

            # Aggregate All + Check inclusion
            if not df_all.empty:
                df_all['sample'] = dna_id
                df_all.insert(1, 'Family_ID', map_id)
                
                # Get sets of keys present in auto and manual
                keys_auto = make_keys(df_auto)
                keys_manual = make_keys(df_manual)
                
                # Apply check logic row by row
                in_auto_flags = []
                in_manual_flags = []
                
                for _, r in df_all.iterrows():
                    k = (str(r['chrom']), str(r['pos']), str(r['ref']), str(r['alt']), str(r['gene']))
                    
                    if k in keys_auto: in_auto_flags.append(True)
                    else: in_auto_flags.append(False)
                        
                    if k in keys_manual: in_manual_flags.append(True)
                    else: in_manual_flags.append(False)

                df_all['in_auto_report'] = in_auto_flags
                df_all['in_manual_review'] = in_manual_flags
                
                all_candidates_list.append(df_all)

        except Exception as e:
            logger.error(f"Error processing {dna_id}: {e}", exc_info=True)

    # Save Master Files
    os.makedirs(args.outdir, exist_ok=True)
    
    if auto_results_list:
        final_auto = pd.concat(auto_results_list, ignore_index=True)
        final_auto.to_csv(os.path.join(args.outdir, "FINAL_auto_conclusions.csv"), index=False)
    
    if manual_results_list:
        final_manual = pd.concat(manual_results_list, ignore_index=True)
        final_manual.to_csv(os.path.join(args.outdir, "FINAL_manual_review.csv"), index=False)

    if all_candidates_list:
        final_all = pd.concat(all_candidates_list, ignore_index=True)
        final_all.to_csv(os.path.join(args.outdir, "FINAL_all_candidates.csv"), index=False)
        
    print(f"\nBATCH FINISHED. Aggregated results in {args.outdir}")
    print(f"Total processed samples: {len(proband_files)}")

def main():
    parser = build_parser()
    args = parser.parse_args()
    db_paths = dict(DB_PATHS_DEFAULT)
    if args.db_paths_json and os.path.exists(args.db_paths_json):
        with open(args.db_paths_json) as fh:
            db_paths.update(json.load(fh))
            
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
            ttn_meta_csv=args.ttn_meta,
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