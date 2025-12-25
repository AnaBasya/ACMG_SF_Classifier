#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACMG Secondary Findings Classifier - Ready-to-run

This script:
- Annotates (optionally) with VEP + dbNSFP
- Integrates ClinVar, HGMD and Internal DB evidence
- Applies gene-specific rules from ACMG SF table (with positional rules parsing)
- Implements PVS1 decision tree, PS1/PM5, PP3/BP4, PM3 batch, TTN meta-exon support
- Adds delins adjacency detection and transcript-position enforcement per recommendations table
- CADD_SUPPORTING threshold set to 30

Dependencies:
  pandas, pysam, cyvcf2 (optional), intervaltree (optional), tabix/bcftools/bgzip for some operations

Usage: see CLI at bottom (run with --help)

Author: adapted for user's dataset and requested positional checks
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
from dataclasses import dataclass, field
from typing import Tuple, List, Dict, Any, Optional, Set
from collections import defaultdict
from cooccurrence_helper import query_pair_gnomad, fetch_cooccurrence_for_pairs

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
    "CADD_SUPPORTING": 30,               

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
    }
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

    "PP5_Supporting": 1,

    "BA1": -8,
    "BS1": -4,
    "BS2": -4,
    "BP4": -1,
    "BP5": -1,
    "BP6": -1,
    "BP7": -1,
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
    hgmd_rankscore: str = "Not provided"
    hgmd_support_positive: int = 0
    hgmd_support_total: int = 0

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

def check_disease_match(acmg_disease_string: str, variant_trait: str, 
                       extra_recs_list: List[Dict] = None,
                       gene_rules: Dict[str, Any] = None) -> Tuple[bool, str, str]:
    vt_raw = str(variant_trait or "").strip()
    acmg_raw = str(acmg_disease_string or "").strip()
    # if no ACMG disease -> cannot match
    if not acmg_raw:
        return False, "no_acmg_disease", "No ACMG disease defined"

    # Clean & normalize
    def norm(s):
        s2 = str(s).lower()
        # replace underscores and pipes with space, collapse punctuation to spaces
        s2 = s2.replace("_", " ").replace("|", " ").replace("/", " ")
        # replace various punctuation with space
        s2 = re.sub(r'[^\w\s]', ' ', s2)
        # multiple whitespace -> single space
        s2 = re.sub(r'\s+', ' ', s2).strip()
        # remove leading/trailing non alphanumeric chars (extra safety)
        s2 = re.sub(r'^[^a-z0-9]+|[^a-z0-9]+$', '', s2)
        return s2

    acmg_norm = norm(acmg_raw)
    vt_norm = norm(vt_raw)

    # If DB trait empty -> mark as unknown DB phenotype (do not treat as explicit absence)
    if not vt_norm:
        return False, "unknown_db_phenotype", "No phenotype text in DB (unknown)"

    # Whitelist/exclusion matching (if provided)
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

    # Direct substring / token intersection (prefer whole-word/token matches)
    if acmg_norm in vt_norm or vt_norm in acmg_norm:
        return True, "direct_match", "Direct substring match"

    # Token intersection ignoring short/stopwords
    stopwords = {}
    acmg_tokens = {t for t in re.split(r'\s+', acmg_norm) if len(t) > 3 and t not in stopwords}
    vt_tokens = {t for t in re.split(r'\s+', vt_norm) if len(t) > 3 and t not in stopwords}

    # exact token overlap
    common = acmg_tokens.intersection(vt_tokens)
    if common:
        return True, "token_match", f"Shared terms: {', '.join(sorted(common))}"

    # try fuzzy token matching (allow small differences: remove common endings)
    def stem_token(tok):
        return re.sub(r'(itis|opathy|osis|ic|al|us|a)$', '', tok)
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
            return parts[-1]
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
    logger.info("Parsing CSQ header from %s", vcf_path)
    if HAVE_PYSAM:
        try:
            vf = pysam.VariantFile(vcf_path)
            for tag in ["CSQ", "ANN"]:
                if tag in vf.header.info:
                    rec = vf.header.info[tag]
                    desc = str(rec.description)
                    m = re.search(r'Format:\s*("([^"]+)"|\'([^\']+)\'|([^">]+))', desc)
                    if m:
                        raw_fmt = (m.group(2) or m.group(3) or m.group(4) or "").strip()
                        header = raw_fmt.split("|")
                        logger.info("Found %s header via Pysam: %d fields", tag, len(header))
                        return header
        except Exception as e:
            logger.debug("Pysam header parse failed: %s", e)

    try:
        open_func = gzip.open if vcf_path.endswith((".gz", ".bgz")) else open
        with open_func(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if line.startswith("##INFO=<ID=CSQ") or line.startswith("##INFO=<ID=ANN"):
                    match = re.search(r'Format:\s*("([^"]+)"|\'([^\']+)\'|([^>]+))', line)
                    if match:
                        raw_fmt = (match.group(2) or match.group(3) or match.group(4) or "").strip().rstrip('>')
                        header = raw_fmt.split("|")
                        logger.info("Found CSQ header via Text: %d fields", len(header))
                        return header
    except Exception as e:
        logger.warning("Text header parse failed: %s", e)

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
# VEP & dbNSFP Helpers (unchanged implementations)
# ---------------------------

def detect_dbnsfp_header_fields(dbnsfp_path: str) -> List[str]:
    if not dbnsfp_path or not os.path.exists(dbnsfp_path):
        return []
    try:
        with gzip.open(dbnsfp_path, "rt", encoding="utf-8", errors="replace") as fh:
            for _ in range(5):
                ln = fh.readline()
                if not ln: break
                ln = ln.rstrip("\n")
                if ln.startswith("#") or (ln and "Chr" in ln and "\t" in ln):
                    clean_ln = ln.lstrip("#").strip()
                    if "\t" in clean_ln:
                        return clean_ln.split("\t")
    except Exception as e:
        logger.debug("Failed to detect dbNSFP header: %s", e)
    return []

def choose_dbnsfp_fields(dbnsfp_path: Optional[str], requested_fields: List[str]) -> Optional[str]:
    if not dbnsfp_path or not os.path.exists(dbnsfp_path):
        return None
    available_fields = detect_dbnsfp_header_fields(dbnsfp_path)
    if not available_fields:
        logger.warning("Could not read dbNSFP header. Passing path only (all fields will be added).")
        return dbnsfp_path
    valid_fields = []
    missing_fields = []
    avail_set = set(available_fields)
    for req in requested_fields:
        if req in avail_set:
            valid_fields.append(req)
        else:
            missing_fields.append(req)
    if missing_fields:
        logger.warning("The following requested dbNSFP fields are missing in the DB file and will be skipped: %s", missing_fields)
    if not valid_fields:
        return dbnsfp_path
    return dbnsfp_path + "," + ",".join(valid_fields)

def run_vep_robust(in_vcf: str, out_vcf: str, vep_cmd: str, vep_cache: str, 
                   fasta: Optional[str], extra_args: List[str], 
                   dbnsfp_plugin_str: Optional[str]) -> str:
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
        "--everything"
    ]
    if vep_cache:
        args += ["--dir_cache", vep_cache]
    if fasta:
        args += ["--fasta", fasta]
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
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error("VEP failed (exit %d). stderr:\n%s", proc.returncode, proc.stderr)
        raise RuntimeError("VEP execution failed")
    try:
        subprocess.run(["tabix", "-p", "vcf", out_vcf], check=False)
    except Exception:
        logger.warning("Could not tabix index VEP output (optional but recommended)")
    has_csq = False
    try:
        check_cmd = ["bcftools", "view", "-h", out_vcf]
        if shutil.which("bcftools"):
            res = subprocess.run(check_cmd, capture_output=True, text=True)
            if "ID=CSQ" in res.stdout:
                has_csq = True
        else:
            with gzip.open(out_vcf, "rt") as fh:
                for _ in range(500):
                    line = fh.readline()
                    if line.startswith("##INFO=<ID=CSQ"):
                        has_csq = True; break
                    if not line.startswith("#"): break
    except Exception as e:
        logger.warning("Could not validate CSQ header presence: %s", e)
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
    def __init__(self, db_paths: Dict[str,str], 
                 gnomad_v2: Optional[str]=None, gnomad_v3: Optional[str]=None, gnomad_v4: Optional[str]=None, 
                 internal_db_path: Optional[str]=None):
        self.paths = dict(db_paths or {})
        self.clinvar_vcf: Optional[pysam.VariantFile] = None
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
        perv, gmax = extract_gnomad_af_from_info(getattr(vr, "_raw_ann", {}) or {})
        if gmax.get("max_af", 0.0) > 0:
            vr.gnomad_af = gmax["max_af"]
            vr.gnomad_details = {"version_pop": gmax, "source": "VEP_CSQ"}
            return
        chrom_clean = vr.chrom.replace("chr", "")
        chrom_vars = [chrom_clean, "chr" + chrom_clean]
        best = {"max_af": 0.0, "version": None, "pop": None, "orig": None}
        for ver, vf in self.gnomad_vcf_handles.items():
            for c_query in chrom_vars:
                try:
                    for rec in vf.fetch(c_query, vr.pos - 1, vr.pos):
                        if int(rec.pos) != int(vr.pos): continue
                        if rec.ref != vr.ref: continue
                        if not (rec.alts and vr.alt in rec.alts): continue
                        info_rec = dict(rec.info)
                        perv2, gmax2 = extract_gnomad_af_from_info(info_rec)
                        if gmax2.get("max_af", 0.0) > best["max_af"]:
                            best.update({"max_af": float(gmax2["max_af"]), "version": ver, "pop": gmax2.get("pop")})
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
                    if best["max_af"] > 0: break 
                except ValueError:
                    continue
                except Exception:
                    continue
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
                if len(parts) < 4: continue
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
            # Read with pandas if available for robustness
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
                support_cands = ["support", "support_field", "support_confirmed", "support_total"]

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

                    # support confirmed/total parsing
                    support_field = get_col_val(r, support_cands, "")
                    support_confirmed_count = 0
                    support_total_count = 0
                    if support_field:
                        try:
                            s = str(support_field).strip().replace("+", "")
                            parts = re.split(r'[;,\/\s]+', s)
                            for part in parts:
                                if part.strip():
                                    try:
                                        val = int(float(re.sub(r'[^\d\-]', '', part.strip())))
                                        support_total_count += 1
                                        if val == 1:
                                            support_confirmed_count += 1
                                    except Exception:
                                        # try parse patterns like "1/0/1"
                                        sub_nums = re.findall(r'\d+', part)
                                        if sub_nums:
                                            for sn in sub_nums:
                                                support_total_count += 1
                                                if int(sn) == 1:
                                                    support_confirmed_count += 1
                        except Exception:
                            pass

                    rank = get_col_val(r, ["rankscore", "rank", "score"], "")
                    phen = get_col_val(r, ["disease", "phenotype", "phen"], "")

                    self.hgmd_index[chrom][pos][ref][alt] = {
                        "class": tag,
                        "phen": phen,
                        "pubs": pub_count,
                        "dmsupported": dmsupported,
                        "support_confirmed": support_confirmed_count,
                        "support_total": support_total_count,
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
                        support_field = r.get("support") or ""
                        support_confirmed_count = 0
                        support_total_count = 0
                        if support_field:
                            parts = re.split(r'[;,\/\s]+', str(support_field))
                            for part in parts:
                                if part.strip():
                                    try:
                                        val = int(float(part.strip()))
                                        support_total_count += 1
                                        if val == 1:
                                            support_confirmed_count += 1
                                    except:
                                        pass
                        rank = r.get("rankscore") or ""
                        phen = r.get("disease") or ""
                        self.hgmd_index[chrom][pos][ref][alt] = {
                            "class": tag,
                            "phen": phen,
                            "pubs": pub_count,
                            "dmsupported": dmsupported,
                            "support_confirmed": support_confirmed_count,
                            "support_total": support_total_count,
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
        # CLINVAR
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
                    clnsig = found_rec.info.get('CLNSIG') or ''
                    clinvar_sig_str = ";".join(clnsig) if isinstance(clnsig, (list, tuple)) else str(clnsig)
                    vr.clinvar_sig = clinvar_sig_str
                    clndn = found_rec.info.get('CLNDN') or ''
                    vr.clinvar_trait = ";".join(clndn) if isinstance(clndn, (list, tuple)) else str(clndn)
                    rev = found_rec.info.get('CLNREVSTAT') or ''
                    rev_str = str(rev).lower()
                    stars = 0
                    if 'practice_guideline' in rev_str: stars = 4
                    elif 'expert_panel' in rev_str: stars = 3
                    elif 'criteria_provided' in rev_str and 'multiple_submitters' in rev_str and 'no_conflict' in rev_str: stars = 2
                    elif 'criteria_provided' in rev_str: stars = 1
                    vr.clinvar_stars = str(stars)
                    sig_lower = clinvar_sig_str.lower()
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
        # INTERNAL DB
        internal_evidence = []
        try:
            entry = self.internal_db_index.get(vr.chrom, {}).get(vr.pos, {}).get(vr.ref, {}).get(vr.alt)
            if entry:
                int_cls = str(entry.get("cls", "")) if isinstance(entry, dict) else str(entry)
                # Только записываем информацию для отображения в таблице
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
                tag = hgmd_entry.get('class', '')
                pubs_str = hgmd_entry.get('pubs', '0')
                hgmd_gene = hgmd_entry.get('gene', '')
                rank = hgmd_entry.get('rankscore', '')
                dmsupported = int(hgmd_entry.get('dmsupported', 0) or 0)
                support_confirmed = int(hgmd_entry.get('support_confirmed', 0) or 0)
                
                # Store values on variant record
                vr.hgmd_dmsupported = dmsupported
                vr.hgmd_support_confirmed = support_confirmed
                
                try:
                    pub_count = int(pubs_str)
                except (ValueError, TypeError):
                    pub_count = 0
                vr.hgmd_publication_count = pub_count
                vr.hgmd_rankscore = rank if rank and rank != "NULL" else "Not provided"
                vr.hgmd_phen = hgmd_entry.get('phen', '')
                vr.hgmd_class = tag
                
                # Check pathogenicity support
                is_well_supported = False
                if dmsupported >= 2:  # At least 2 points of support
                    is_well_supported = True
                elif support_confirmed >= 2:  # At least 2 confirming publications
                    is_well_supported = True
                
                # Gene validation
                gene_valid, gene_reason = validate_gene_match_db_evidence(
                    vr.gene, hgmd_gene, "HGMD"
                )
                
                if not gene_valid:
                    vr.manual_reasons.append(gene_reason)
                    vr.hgmd_class = f"Gene mismatch: {tag}"
                else:
                    if tag == 'DM' and is_well_supported:
                        vr.has_strong_pathogenic_evidence = True
                        vr.is_hq_pathogenic_db = True
                        hgmd_evidence.append(("Pathogenic", f"HGMD_DM_supported{dmsupported}"))
                    elif tag == 'DM':
                        vr.has_weak_pathogenic_evidence = True
                        hgmd_evidence.append(("Pathogenic", "HGMD_DM_weak_support"))
                    else:
                        hgmd_evidence.append((tag, "HGMD"))
        except Exception as e:
            logger.debug(f"HGMD dictionary lookup error: {e}")
        results.extend(hgmd_evidence)
        # Check conflicts
        hq_sources = [e for e in results if any(q in e[1] for q in ["ClinVar_","Internal_DB","HGMD_DM_"])]
        if len(hq_sources) >= 2:
            classifications = set(e[0].lower() for e in hq_sources)
            has_path = any('pathogen' in c for c in classifications)
            has_benign = any('benign' in c for c in classifications)
            if has_path and has_benign:
                vr.is_hq_conflict_db = True
                vr.conflict_details = f"Conflict between HQ sources: {', '.join([f'{c[0]} ({c[1]})' for c in hq_sources])}"
        vr.db_hits = results
        return results

    # ClinVar protein index builder (same as earlier)
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
                    idx = info.find("GENEINFO=")
                    if idx == -1: continue
                    rest = info[idx+9:]
                    end_idx = rest.find(";")
                    val = rest[:end_idx] if end_idx != -1 else rest
                    if not val: continue
                    genes_raw = []
                    for gene_pair in val.split("|"):
                        if ":" in gene_pair:
                            gene_raw = gene_pair.split(":")[0]
                            gene_norm = normalize_gene_name(gene_raw)
                            if gene_norm:
                                genes_raw.append((gene_raw, gene_norm))
                    if not genes_raw:
                        continue
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

    def find_protein_variant_by_hgvsp(self, gene: str, hgvsp: str, pos_aa: int, transcript: str) -> List[Tuple[str, str, dict]]:
        if not hgvsp or not gene:
            return []
        norm_gene = normalize_gene_name(gene)
        results = []
        for entry in self.clinvar_protein_index.get(norm_gene, []):
            if entry.get('gene_norm', '') != norm_gene:
                continue
            if entry.get('stars', 0) < 2:
                continue
            clnsig = str(entry['clnsig']).lower()
            if "benign" in clnsig and "pathogen" not in clnsig:
                continue
            entry_hgvsp_clean = clean_hgvsp(entry['hgvsp'])
            current_hgvsp_clean = clean_hgvsp(hgvsp)
            if entry_hgvsp_clean == current_hgvsp_clean:
                if "pathogenic" in clnsig and "benign" not in clnsig and "likely" not in clnsig:
                    results.append((entry['clnsig'], "ClinVar_Pathogenic", entry))
                elif "likely pathogenic" in clnsig or ("likely" in clnsig and "pathogen" in clnsig):
                    results.append((entry['clnsig'], "ClinVar_Likely_pathogenic", entry))
                else:
                    results.append((entry['clnsig'], "ClinVar", entry))
        return results
    
    def check_pm5_variants(self, gene: str, pos_aa: int, current_hgvsp: str) -> Optional[str]:
        """
        Checks for other pathogenic/likely-pathogenic variants at the same amino-acid residue.
        Returns:
        - "Pathogenic" if a Pathogenic (not Likely) variant at same residue (≥1★) is found;
        - "Likely pathogenic" if only Likely pathogenic matches found;
        - None if no suitable matches.
        Logic:
        - Uses self.clinvar_protein_index[normalize_gene]
        - Skips exact same protein change (this would be PS1)
        - Requires at least 1-star in ClinVar index (stars >= 2)
        - Ignores benign-only annotations
        """
        if not pos_aa or not gene:
            return None

        norm_gene = normalize_gene_name(gene)
        best_class = None

        # iterate ClinVar protein-index for the gene
        for entry in self.clinvar_protein_index.get(norm_gene, []):
            # skip if entry lacks protein annotation or gene mismatch
            entry_hgvsp = entry.get('hgvsp', '')
            if not entry_hgvsp:
                continue

            # skip exact same protein change (PS1 case)
            try:
                if clean_hgvsp(entry_hgvsp) == clean_hgvsp(current_hgvsp):
                    continue
            except Exception:
                pass

            # parse AA position from entry
            m = re.search(r'p\.[A-Za-z]{1,3}(\d+)', clean_hgvsp(entry_hgvsp) or "")
            if not m:
                # fallback: try generic digits
                m2 = re.search(r'(\d+)', clean_hgvsp(entry_hgvsp) or "")
                if not m2:
                    continue
                try:
                    entry_pos = int(m2.group(1))
                except:
                    continue
            else:
                try:
                    entry_pos = int(m.group(1))
                except:
                    continue

            if entry_pos != pos_aa:
                continue

            # check significance & stars
            clnsig = str(entry.get('clnsig', '')).lower()
            stars = int(entry.get('stars', 0) or 0)

            # require at least 1 star
            if stars < 2:
                continue

            # skip benign entries
            if "benign" in clnsig and "pathogen" not in clnsig:
                continue

            # if exact Pathogenic (not likely) — return strongest
            if "pathogenic" in clnsig and "likely" not in clnsig:
                return "Pathogenic"

            # if LP found, keep as best candidate
            if "likely pathogenic" in clnsig or ("likely" in clnsig and "pathogen" in clnsig):
                best_class = "Likely pathogenic"

        return best_class

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
        df_clean._transcript_last_exons = transcript_last_exons
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

    # Try pandas first for robustness
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

    # Try pandas for robust parsing
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
    Predict if variant `vr` triggers NMD.

    Improvements:
     - If variant EXON field (from VEP/CSQ) explicitly indicates last exon (e.g. "42/42"),
       treat as escaping NMD (conservative) and return False immediately.
     - Otherwise fallback to nmd_map / exons_map / transcript mapping logic as before.
    """
    details: Dict[str, Any] = {"checked_transcripts": [], "reason": "", "variant_pos": None}
    try:
        chrom = getattr(vr, "chrom", None)
        pos = int(getattr(vr, "pos"))
        details["variant_pos"] = f"{chrom}:{pos}"
    except Exception:
        details["reason"] = "no_variant_position"
        return False, details

    # QUICK: honor explicit VEP/CSQ exon annotation "curr/total" - if curr == total -> last exon -> escapes NMD
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

            # variant in last exon -> escapes NMD
            if last_start <= pos <= last_end:
                tx_result["predicts_nmd"] = False
                tx_result["details"]["reason"] = "variant_in_last_exon"
                per_tx_results.append(tx_result)
                continue

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
# ACMG Engine
# ---------------------------
class ACMGEngine:
    def __init__(self, dbm: DatabaseManager, gene_rules: Dict[str, Dict[str,Any]], 
                 exons_df: Optional[Any] = None, ttn_meta: Dict[str, IntervalTree] = None, spliceai_thr: float = 0.20,
                 require_mane: bool = False):
        self.dbm = dbm
        self.gene_rules = {normalize_gene_name(k): v for k,v in (gene_rules or {}).items()}
        self.exons_df = exons_df
        self.ttn_meta = ttn_meta or {}
        self.spliceai_thr = spliceai_thr
        self.require_mane = require_mane

    def _get_gene_rule(self, gene: str) -> Dict[str,Any]:
        return self.gene_rules.get((gene or "").upper(), {"reportable": False, "special_flags": {}, "raw_row": {}})

    def _check_upstream_pathogenic_variants(self, vr: VariantRecord, upstream_distance: int = 1000) -> bool:
        if not hasattr(self, 'dbm') or self.dbm is None:
            logger.debug("No DatabaseManager available for upstream variant check")
            return False
        if self.dbm.clinvar_vcf is None:
            logger.debug("ClinVar VCF not available for upstream variant check")
            return False
        try:
            upstream_start = max(0, vr.pos - upstream_distance)
            upstream_end = vr.pos - 1
            if upstream_start >= upstream_end:
                logger.debug(f"Invalid upstream region for {vr.chrom}:{vr.pos}")
                return False
            gene_norm = normalize_gene_name(vr.gene)
            pathogenic_found = False
            pathogenic_details = []
            chrom_variants = [vr.chrom]
            if vr.chrom.startswith("chr"):
                chrom_variants.append(vr.chrom.replace("chr", ""))
            else:
                chrom_variants.append(f"chr{vr.chrom}")
            for chrom_query in chrom_variants:
                try:
                    for rec in self.dbm.clinvar_vcf.fetch(chrom_query, upstream_start, upstream_end):
                        gene_info = rec.info.get('GENEINFO', '')
                        rec_genes = []
                        if gene_info:
                            for gene_pair in str(gene_info).split('|'):
                                if ':' in gene_pair:
                                    clinvar_gene = gene_pair.split(':')[0]
                                    rec_genes.append(normalize_gene_name(clinvar_gene))
                        if not rec_genes or gene_norm not in rec_genes:
                            continue
                        hgvsp_info = rec.info.get('HGVS_Protein', '')
                        if not hgvsp_info:
                            continue
                        clnsig = rec.info.get('CLNSIG')
                        clnsig_str = ";".join(clnsig) if isinstance(clnsig, (list, tuple)) else str(clnsig)
                        rev = rec.info.get('CLNREVSTAT') or ''
                        rev_str = str(rev).lower()
                        stars = 0
                        if 'practice_guideline' in rev_str: stars = 4
                        elif 'expert_panel' in rev_str: stars = 3
                        elif 'criteria_provided' in rev_str and 'multiple_submitters' in rev_str and 'no_conflict' in rev_str: stars = 2
                        elif 'criteria_provided' in rev_str: stars = 1
                        if stars >= 2:
                            if ('pathogenic' in clnsig_str or 'likely pathogenic' in clnsig_str) and 'benign' not in clnsig_str:
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
            return False

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

    def _check_protein_truncation_percentage(self, vr: VariantRecord) -> bool:
        if not vr.hgvsp:
            return True
        pos_match = re.search(r'p\.[A-Za-z*]{1,3}(\d+)', vr.hgvsp)
        if not pos_match:
            pos_match = re.search(r'p\.(\d+)', vr.hgvsp)
            if not pos_match:
                return True
        try:
            trunc_pos = int(pos_match.group(1))
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
            gene_key = vr.gene.upper()
            if gene_key in protein_lengths:
                protein_len = protein_lengths[gene_key]
                percentage = (trunc_pos / protein_len) * 100
                logger.debug(f"Truncation check for {vr.gene}: pos={trunc_pos}, len={protein_len}, %={percentage:.1f}%")
                return percentage > 10
            else:
                return trunc_pos > 100
        except Exception as e:
            logger.debug(f"Error checking protein length for NMD: {e}")
            return True

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

    def _ttn_overlap(self, chrom: str, pos: int) -> Optional[dict]:
        chromn = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
        if chromn not in self.ttn_meta:
            return None
        hits = self.ttn_meta[chromn].at(pos - 1)
        if not hits:
            return None
        return max(hits, key=lambda it: float(it.data.get("psi_meta", 0) or 0)).data

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
        af = vr.gnomad_af or 0.0

        # -------------------------
        # 0. BA1 (stand-alone benign)
        # -------------------------
        is_ba1_exception = False
        exceptions = THRESH.get("BA1_EXCEPTIONS", {})
        # vr.hgvsp may be empty; guard
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
            # BA1 stand-alone — we can short-circuit and return Benign
            total = sum(pts.values())
            return assigned, pts, "Benign", False, manual_reasons

        # -------------------------
        # 1. BS1 (strong benign by freq)
        # -------------------------
        if af >= THRESH.get("BS1_AF", 0.01):
            assigned.append("BS1")
            pts["BS1"] = POINTS_MAP.get("BS1", -4)
            vr.criteria_explanations["BS1"] = f"BS1: gnomAD max AF={af:.6g} >= BS1 threshold {THRESH.get('BS1_AF')}"

        # -------------------------
        # 2. PVS1 Decision tree
        # -------------------------
        special = self._get_gene_rule(vr.gene).get("special_flags", {}) if self.gene_rules else {}
        # support both older/short and VEP-style consequence tokens
        is_lof_molecular = any(x in cons for x in ["stop_gained", "frameshift", "frameshift_variant", "start_lost", "stop_lost"])
        is_canonical_splice = any(x in cons for x in ["splice_acceptor_variant", "splice_donor_variant", "splice_acceptor", "splice_donor"])
        logger.debug(f"PVS1 check for {vr.gene}:{vr.hgvsp} - Consequence: {cons}")
        logger.debug(f"  is_lof_molecular: {is_lof_molecular}")
        logger.debug(f"  is_canonical_splice: {is_canonical_splice}")
        logger.debug(f"  lof_not_reportable flag: {special.get('lof_not_reportable', False)}")

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

                # Final safeguard: if internal predictor explicitly indicates escape at any time, honour it
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
                    # Default: moderate, but check for downstream in-frame ATG and upstream pathogenic evidence
                    assigned.append("PVS1_Moderate")
                    pts["PVS1_Moderate"] = POINTS_MAP.get("PVS1_Moderate", 2)
                    # Explain decision
                    atg_pos, context = self._find_nearest_inframe_start_codon(vr, fasta_handle=fasta_handle)
                    expl = "PVS1_Moderate: start-loss variant."
                    if atg_pos:
                        expl += f" Nearest in-frame ATG at pos {atg_pos} (maintains frame) -> supports moderate rather than very strong."
                    else:
                        expl += " No in-frame downstream ATG found within search window."
                    # Check upstream pathogenic variants (strengthen caution)
                    upstream_found = False
                    try:
                        upstream_found = self._check_upstream_pathogenic_variants(vr)
                    except Exception:
                        upstream_found = False
                    if upstream_found:
                        expl += " Upstream pathogenic variants found -> manual inspection recommended."
                    vr.criteria_explanations["PVS1"] = expl

                elif triggers_nmd:
                    # Triggers NMD -> PVS1_VeryStrong
                    assigned.append("PVS1_VeryStrong")
                    pts["PVS1_VeryStrong"] = POINTS_MAP.get("PVS1_VeryStrong", 8)
                    vr.criteria_explanations["PVS1"] = ("PVS1_VeryStrong: ... predicted NMD")
                else:
                    # Escapes NMD => DO NOT assign PVS1. For truncating variants in last exon, use PM4 instead.
                    # For canonical splice variants in last exon, treat conservatively as PM4_Moderate
                    if is_canonical_splice:
                        assigned.append("PM4_Moderate")
                        pts["PM4_Moderate"] = POINTS_MAP.get("PM4_Moderate", 2)
                        vr.criteria_explanations["PM4"] = ("PM4_Moderate: canonical splice variant located in last exon / predicted to escape NMD; treated as protein-altering in last exon")
                    else:
                        # Stop/Frameshift in last exon -> prefer PM4_Moderate (unless very large deletion -> PM4_Strong)
                        long_trunc = self._check_protein_truncation_percentage(vr)
                        if long_trunc:
                            # If truncation removes >10% — still prefer PM4_Moderate 
                            assigned.append("PM4_Moderate")
                            pts["PM4_Moderate"] = POINTS_MAP.get("PM4_Moderate", 2)
                            vr.criteria_explanations["PM4"] = ("PM4_Moderate: truncating variant in last exon removing >10% of protein length; do not apply PVS1 for last-exon truncations")
                        else:
                            assigned.append("PM4_Moderate")
                            pts["PM4_Moderate"] = POINTS_MAP.get("PM4_Moderate", 2)
                            vr.criteria_explanations["PM4"] = ("PM4_Moderate: truncating variant in last exon removing <=10% of protein (do not apply PVS1)")

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
        if "missense" in cons and pos_aa:
            # PS1: exact same AA change reported Pathogenic/LP
            matches = []
            try:
                matches = self.dbm.find_protein_variant_by_hgvsp(vr.gene, vr.hgvsp, pos_aa, vr.transcript) or []
            except Exception:
                matches = []
            ps1_assigned = False
            # matches may be in various formats; allow tuples/lists/dicts
            for m in matches:
                # m might be (clnsig, source) or (clnsig, source, entry_dict)
                clnsig = ""
                src = ""
                info = {}
                try:
                    if isinstance(m, (list, tuple)):
                        clnsig = str(m[0] or "")
                        src = str(m[1] or "")
                        if len(m) > 2 and isinstance(m[2], dict):
                            info = m[2]
                    elif isinstance(m, dict):
                        clnsig = str(m.get("clnsig",""))
                        src = str(m.get("source",""))
                        info = m
                except Exception:
                    continue

                cln_low = clnsig.lower()
                if "pathogenic" in cln_low and "likely" not in cln_low:
                    # Strong PS1
                    assigned.append("PS1_Strong")
                    pts["PS1_Strong"] = POINTS_MAP.get("PS1_Strong", 4)
                    ps1_assigned = True
                    vr.ps1_matches.append({"hgvsp": info.get("hgvsp", vr.hgvsp), "clnsig": clnsig, "source": src, "stars": info.get("stars", "")})
                    break
                elif "likely pathogenic" in cln_low or ("likely" in cln_low and "pathogen" in cln_low):
                    # Moderate or Supporting depending on source; conservative: Supporting
                    assigned.append("PS1_Supporting")
                    pts["PS1_Supporting"] = POINTS_MAP.get("PS1_Supporting", 1)
                    ps1_assigned = True
                    vr.ps1_matches.append({"hgvsp": info.get("hgvsp", vr.hgvsp), "clnsig": clnsig, "source": src, "stars": info.get("stars", "")})
                    # do not break — we prefer a real Pathogenic if present, but conservatively break to avoid duplicating
                    break

            if not ps1_assigned:
                # PM5: different AA change at same residue reported pathogenic/likely pathogenic
                try:
                    pm5_match = self.dbm.check_pm5_variants(vr.gene, pos_aa, vr.hgvsp)
                except Exception:
                    pm5_match = None
                if pm5_match:
                    if pm5_match == "Pathogenic":
                        assigned.append("PM5_Moderate")
                        pts["PM5_Moderate"] = POINTS_MAP.get("PM5_Moderate", 2)
                        vr.pm5_matches.append({"best": "Pathogenic"})
                        vr.criteria_explanations["PM5"] = f"PM5_Moderate: other different AA at same residue reported Pathogenic in database (residue {pos_aa})"
                    else:
                        assigned.append("PM5_Supporting")
                        pts["PM5_Supporting"] = POINTS_MAP.get("PM5_Supporting", 1)
                        vr.pm5_matches.append({"best": pm5_match})
                        vr.criteria_explanations["PM5"] = f"PM5_Supporting: other different AA at same residue reported {pm5_match} in database"

            # PS1 explanation (if assigned)
            if any(c.startswith("PS1") for c in assigned):
                # describe matches
                try:
                    details = "; ".join([f"{m.get('hgvsp','?')} ({m.get('clnsig','')})" for m in vr.ps1_matches])
                except Exception:
                    details = "Matched pathogenic protein change in DB"
                vr.criteria_explanations["PS1"] = f"PS1: exact same protein change observed in DB: {details}"

        # -------------------------
        # PM2 (population rarity) - MOI specific thresholds
        # -------------------------
        rule = self._get_gene_rule(vr.gene) if self.gene_rules else {}
        moi = str(rule.get("moi","")).upper() if rule else ""
        pm2_thr = THRESH.get("PM2_AD_AF", 1e-4)
        if "AR" in moi or "RECESSIVE" in moi:
            pm2_thr = THRESH.get("PM2_AR_AF", THRESH.get("PM2_AD_AF", 1e-4))
        elif "XLR" in moi:
            pm2_thr = THRESH.get("PM2_XLR_AF", THRESH.get("PM2_AD_AF", 1e-4))
        elif "XLD" in moi:
            pm2_thr = THRESH.get("PM2_XLD_AF", THRESH.get("PM2_AD_AF", 1e-4))

        if af < pm2_thr:
            assigned.append("PM2_Moderate")
            pts["PM2_Moderate"] = POINTS_MAP.get("PM2_Moderate", 2)
            vr.criteria_explanations["PM2"] = f"PM2_Moderate: gnomAD max AF={af:.6g} < threshold {pm2_thr} for MOI {moi}"

        # -------------------------
        # PP3 / BP4 (computational predictions)
        # -------------------------
        s_score = vr.spliceai if vr.spliceai is not None else 0.0
        r_score = vr.revel if vr.revel is not None else -1.0
        a_score = vr.alpha_missense if vr.alpha_missense is not None else -1.0
        # CADD threshold name compatibility
        cadd_thr = THRESH.get("CADD_SUPPORTING", THRESH.get("CADD_MODERATE", 30))

        is_pp3_moderate = False
        if s_score >= THRESH.get("SPLICEAI_MODERATE", 0.20):
            is_pp3_moderate = True
        elif r_score >= THRESH.get("REVEL_MODERATE", 0.932):
            is_pp3_moderate = True
        elif c_score := (vr.cadd if vr.cadd is not None else 0.0):
            if c_score >= cadd_thr:
                # treat CADD as supportive/moderate depending on threshold chosen; here we use is_pp3_moderate true
                is_pp3_moderate = True

        # Do NOT apply PP3 if PVS1 was applied or the variant is an obvious LOF/canonical splice consequence.
        pvs1_present = any(x.startswith("PVS1") for x in assigned)
        cons_lof_like = any(x in cons for x in ["stop_gained", "frameshift", "frameshift_variant", "start_lost", "stop_lost"])
        # is_canonical_splice was computed earlier; keep conservative check for splice_acceptor/donor keywords
        is_canonical_splice_local = any(x in cons for x in ["splice_acceptor_variant", "splice_donor_variant","splice_acceptor","splice_donor"])

        if pvs1_present or cons_lof_like or is_canonical_splice_local:
            # Do not assign PP3 in these situations. Add explanation so curators/audit knows why PP3 skipped.
            vr.criteria_explanations["PP3"] = (
                "PP3 not applied: variant is LOF/canonical-splice or PVS1 was assigned; "
                "computational evidence is not stacked with PVS1/LOF for this pipeline."
            )
        else:
            # Only consider PP3 when PVS1 is not present and variant is not an obvious LOF/splice case
            if is_pp3_moderate:
                assigned.append("PP3_Moderate")
                pts["PP3_Moderate"] = POINTS_MAP.get("PP3_Moderate", 2)
                vr.criteria_explanations["PP3"] = (f"PP3_Moderate: computational evidence (SpliceAI={s_score:.3f}, REVEL={r_score if r_score!=-1 else 'NA'}, "
                                                f"CADD={vr.cadd if vr.cadd is not None else 'NA'}) met moderate thresholds")
            elif r_score >= THRESH.get("REVEL_SUPPORTING", 0.644) or a_score >= THRESH.get("ALPHA_MISSENSE_SUPPORTING", 0.564):
                assigned.append("PP3_Supporting")
                pts["PP3_Supporting"] = POINTS_MAP.get("PP3_Supporting", 1)
                vr.criteria_explanations["PP3"] = f"PP3_Supporting: REVEL={r_score:.3f} or AlphaMissense={a_score:.3f} >= supporting thresholds"
        # -------------------------
        # PP5 (reputable source) - high-quality db evidence used conservatively
        # -------------------------
        is_hq = False
        try:
            if vr.clinvar_stars and str(vr.clinvar_stars).isdigit() and int(vr.clinvar_stars) >= 2:
                sig_low = (vr.clinvar_sig or "").lower()
                if ("pathogenic" in sig_low or "likely pathogenic" in sig_low) and "benign" not in sig_low:
                    is_hq = True
        except Exception:
            is_hq = is_hq

        try:
            if getattr(vr, "hgmd_class", "") == "DM" and int(getattr(vr, "hgmd_publication_count", 0) or 0) >= 2:
                is_hq = True
        except Exception:
            pass

        #try:
         #   if getattr(vr, "internal_db_sig", "") and "pathogen" in (vr.internal_db_sig or "").lower():
                # internal DB considered high-quality for PP5 if signal indicates Pathogenic/LP
          #      is_hq = True
        #except Exception:
          #  pass

        if is_hq and not getattr(vr, "is_hq_conflict_db", False):
            assigned.append("PP5_Supporting")
            pts["PP5_Supporting"] = POINTS_MAP.get("PP5_Supporting", 1)
            vr.criteria_explanations["PP5"] = "PP5_Supporting: high-quality database annotation (ClinVar≥2★ or HGMD DM≥2 pubs or Internal DB P/LP)"

        # -------------------------
        # Final scoring & class assignment (Tavtigian points model)
        # -------------------------
        total = sum(pts.values())
        path_cnt = sum(1 for x in assigned if any(p in x for p in ["PVS","PS","PM","PP"]))

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

        # Minimum pathogenic criteria count guard
        if (auto_class.startswith("Path") or auto_class.startswith("Likely path")) and path_cnt < MIN_PATHOGENIC_CRITERIA:
            auto_class = "VUS"
            manual_reasons.append(f"Downgraded: <{MIN_PATHOGENIC_CRITERIA} pathogenic criteria present")

        # Ensure that any criteria present but missing textual explanation get a default text
        for crit in assigned:
            base = crit.split("_")[0]  # e.g. PVS1_VeryStrong -> PVS1
            if base not in vr.criteria_explanations:
                vr.criteria_explanations[base] = f"{base}: applied as {crit}; detailed reasoning not captured"

        return assigned, pts, auto_class, False, manual_reasons

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
        by_sample_gene = defaultdict(lambda: defaultdict(list))
        for v in candidates:
            by_sample_gene[v.sample][v.gene].append(v)
        for sample, genes in by_sample_gene.items():
            for gene, varlist in genes.items():
                rule = self._get_gene_rule(gene)
                moi = str(rule.get("moi", "")).upper()
                if not ("AR" in moi or "RECESSIVE" in moi or "XLR" in moi):
                    continue
                for v in varlist:
                    pm3_score = 0.0
                    unphased_contributions = 0.0
                    criteria_details = []
                    if _is_hom(v.proband_gt):
                        pm3_score += 0.5
                        unphased_contributions += 0.5
                        criteria_details.append("homozygous: +0.5")
                    for other_v in varlist:
                        if other_v == v:
                            continue
                        in_trans_confirmed = False
                        if (_is_het(v.father_gt) and _is_wt(v.mother_gt) and
                            _is_het(other_v.mother_gt) and _is_wt(other_v.father_gt)):
                            in_trans_confirmed = True
                        elif (_is_het(v.mother_gt) and _is_wt(v.father_gt) and
                            _is_het(other_v.father_gt) and _is_wt(other_v.mother_gt)):
                            in_trans_confirmed = True
                        in_cis_confirmed = False
                        if (_is_het(v.father_gt) and _is_het(other_v.father_gt)) or (_is_het(v.mother_gt) and _is_het(other_v.mother_gt)):
                            in_cis_confirmed = True
                        if in_trans_confirmed:
                            if other_v.automated_class in ["Pathogenic", "Likely pathogenic"]:
                                pm3_score += 1.0
                                criteria_details.append(f"confirmed in trans with P/LP ({other_v.hgvsp}): +1.0")
                            elif other_v.automated_class == "VUS":
                                pm3_score += 0.25
                                criteria_details.append(f"confirmed in trans with VUS ({other_v.hgvsp}): +0.25")
                        elif in_cis_confirmed:
                            criteria_details.append(f"confirmed in cis with {other_v.hgvsp} (no PM3 points)")
                        else:
                            contribution = 0.0
                            if other_v.automated_class == "Pathogenic":
                                contribution = 0.5
                                criteria_details.append(f"unphased with P ({other_v.hgvsp}): +0.5")
                            elif other_v.automated_class == "Likely pathogenic":
                                contribution = 0.5
                                criteria_details.append(f"unphased with LP ({other_v.hgvsp}): +0.5")
                            if unphased_contributions + contribution <= 1.0:
                                unphased_contributions += contribution
                                pm3_score += contribution
                            elif unphased_contributions < 1.0:
                                remaining = 1.0 - unphased_contributions
                                unphased_contributions += remaining
                                pm3_score += remaining
                                criteria_details.append(f"unphased (capped): +{remaining:.2f}")
                    if unphased_contributions > 1.0:
                        excess = unphased_contributions - 1.0
                        pm3_score -= excess
                        unphased_contributions = 1.0
                    pm3_level = None
                    if pm3_score >= 4.0:
                        pm3_level = "PM3_VeryStrong"
                    elif pm3_score >= 2.0:
                        pm3_level = "PM3_Strong"
                    elif pm3_score >= 1.0:
                        pm3_level = "PM3_Moderate"
                    elif pm3_score >= 0.5:
                        pm3_level = "PM3_Supporting"
                    has_parental_data = any(
                        any(gt != "./." for gt in [v.father_gt, v.mother_gt])
                        for v in varlist
                    )
                    if pm3_level and not has_parental_data:
                        if pm3_level in ["PM3_Moderate", "PM3_Strong", "PM3_VeryStrong"]:
                            old_level = pm3_level
                            pm3_level = "PM3_Supporting"
                            v.manual_reasons.append(f"PM3 downgraded from {old_level} to {pm3_level} (no parental phasing)")
                    if pm3_level and pm3_level not in v.criteria_assigned:
                        v.criteria_assigned.append(pm3_level)
                        v.criteria_points[pm3_level] = POINTS_MAP.get(pm3_level, 1)
                        if criteria_details:
                            details_str = ", ".join(criteria_details)
                            v.manual_reasons.append(f"PM3 ({pm3_level}): total={pm3_score:.2f} [{details_str}]")
                        else:
                            v.manual_reasons.append(f"PM3 ({pm3_level}): total={pm3_score:.2f}")
                    if pm3_score > 0 and has_parental_data:
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

# ---------------------------
# Deduplication helper
# ---------------------------
def deduplicate_candidates(candidates: List[VariantRecord]) -> List[VariantRecord]:
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
        e_ann = getattr(existing, "_raw_ann", {}) or {}
        v_ann = getattr(v, "_raw_ann", {}) or {}
        e_can = str(e_ann.get("MANE_PLUS_CLINICAL") or e_ann.get("MANE_SELECT") or e_ann.get("CANONICAL") or "").upper()
        v_can = str(v_ann.get("MANE_PLUS_CLINICAL") or v_ann.get("MANE_SELECT") or v_ann.get("CANONICAL") or "").upper()
        if v_can in ("YES","Y","TRUE","1") and not (e_can in ("YES","Y","TRUE","1")):
            keymap[key] = v
            continue
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
    return list(keymap.values())

# ---------------------------
# Delins detection helper
# ---------------------------
def detect_delins(candidates: List[VariantRecord]) -> None:
    """
    Improved delins detection by interval overlap / clustering.

    For each sample and chromosome:
      - compute genomic span for variant as [pos, pos + max(len(ref), len(alt)) - 1]
      - consider two variants overlapping if spans intersect
      - build connected components (clusters) of overlapping variants
      - for each cluster of size >= 2, mark variants with _delins_flag, assign _delins_cluster_id
        and store cluster member keys (for later triage decision).
      - add an informational manual_reasons note describing the cluster membership (so it will
        appear in 'comments' for auto variants).
    """
    # helper: build key
    def key_of(v: VariantRecord) -> str:
        return f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.sample}"

    # organize by sample and chromosome to avoid cross-sample/cross-chrom checks
    by_sample_chr = defaultdict(lambda: defaultdict(list))
    for v in candidates:
        try:
            sam = v.sample or ""
            chrom = v.chrom or ""
            by_sample_chr[sam][chrom].append(v)
        except Exception:
            continue

    cluster_counter = 0
    # clear previous cluster metadata if present
    for v in candidates:
        if hasattr(v, "_delins_flag"):
            delattr(v, "_delins_flag")
        v._delins_flag = False
        v._delins_cluster_id = None
        v._delins_cluster_members = []

    for sample, chrom_map in by_sample_chr.items():
        for chrom, varlist in chrom_map.items():
            n = len(varlist)
            if n < 2:
                continue
            # compute spans
            spans = []
            for v in varlist:
                try:
                    s = int(v.pos)
                except Exception:
                    s = 0
                ref_len = len(v.ref or "")
                alt_len = len(v.alt or "")
                span_len = max(ref_len, alt_len, 1)
                start = s
                end = s + span_len - 1
                spans.append((v, start, end))
            # build adjacency graph where edges link overlapping spans
            neighbors = {key_of(v): set() for v,_,_ in spans}
            for i in range(len(spans)):
                vi, si, ei = spans[i]
                ki = key_of(vi)
                for j in range(i+1, len(spans)):
                    vj, sj, ej = spans[j]
                    kj = key_of(vj)
                    # overlap test: intervals intersect if (si <= ej and sj <= ei)
                    overlap = (si <= ej and sj <= ei)
                    if overlap:
                        neighbors[ki].add(kj)
                        neighbors[kj].add(ki)
            # find connected components (clusters) in adjacency graph
            visited = set()
            for v_obj, s_v, e_v in spans:
                k = key_of(v_obj)
                if k in visited:
                    continue
                if not neighbors.get(k):
                    visited.add(k)
                    continue
                # BFS to collect component
                comp = []
                queue = [k]
                while queue:
                    cur = queue.pop()
                    if cur in visited:
                        continue
                    visited.add(cur)
                    comp.append(cur)
                    for nb in neighbors.get(cur, set()):
                        if nb not in visited:
                            queue.append(nb)
                if len(comp) < 2:
                    continue
                # Assign a cluster id and annotate members
                cluster_counter += 1
                # Resolve keys -> VariantRecord objects for the component
                comp_vars = []
                key_to_var = {key_of(vv): vv for vv in varlist}
                for kk in comp:
                    vv = key_to_var.get(kk)
                    if vv:
                        comp_vars.append(vv)
                # annotate each variant in cluster
                member_str = ";".join(comp)
                for vv in comp_vars:
                    vv._delins_flag = True
                    vv._delins_cluster_id = cluster_counter
                    vv._delins_cluster_members = comp[:]  # list of keys
                    note = f"Delins cluster detected (cluster_id={cluster_counter}, members={member_str}) - possible complex indel / delins adjacency"
                    if note not in getattr(vv, "manual_reasons", []):
                        if not hasattr(vv, "manual_reasons") or vv.manual_reasons is None:
                            vv.manual_reasons = []
                        vv.manual_reasons.append(note)

# ---------------------------
# Position rules loader & checker (from recommendations table)
# ---------------------------
def load_recommendations_table(path: str) -> Dict[str, Dict[str,Any]]:
    """
    Parses the ACMG recommendations table (annotations) and extracts:
      - requires_biallelic, missense_only, hemizygous_only_male
      - positional rules:
          * transcript-specific exon ranges (e.g., "exons 1-8 of NM_001256021.1")
          * generic exon ranges (e.g., "exons 2-5")
          * protein position ranges ("positions 100-200")
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

def variant_row_with_disease(v: VariantRecord, rule: Dict[str,Any], family_id: str = "", recs_map: Dict = None) -> dict:
    acmg_disease = rule.get("disease", "") if isinstance(rule, dict) else ""
    has_guidance = "No"
    if recs_map and v.gene in recs_map:
        if recs_map[v.gene]: has_guidance = "Yes"
    # prepare various textual fields
    reasons = "; ".join(v.manual_reasons) if v.manual_reasons else ""
    clinvar_phen = v.clinvar_trait
    if not clinvar_phen or clinvar_phen.lower() in [".", "not specified", "not provided", ""]:
        clinvar_phen = "Not provided"
    hgmd_phen = v.hgmd_phen
    if not hgmd_phen or hgmd_phen.lower() in [".", "not specified", "not provided", ""]:
        hgmd_phen = "Not provided"
    revel_val = f"{v.revel:.3f}" if v.revel is not None else "Not provided"
    cadd_val = f"{v.cadd:.2f}" if v.cadd is not None else "Not provided"
    splice_val = f"{v.spliceai:.3f}" if v.spliceai is not None else "Not provided"
    alpha_val = f"{v.alpha_missense:.3f}" if v.alpha_missense is not None else "Not provided"
    zyg = "Het"
    if _is_hom(v.proband_gt):
        zyg = "Hom"
        if v.sample_sex == "Male" and ("X" in v.chrom or "Y" in v.chrom):
            zyg = "Hemi"

    # database summary
    db_sum_parts = []
    if v.clinvar_sig:
        stars = v.clinvar_stars if v.clinvar_stars else "0"
        db_sum_parts.append(f"ClinVar: {v.clinvar_sig} ({stars}*)")
    if v.hgmd_class:
        pubs = v.hgmd_publication_count if v.hgmd_publication_count else 0
        # try to fetch dmsupported/support_confirmed if present
        dmsup = getattr(v, "hgmd_dmsupported", 0) or 0
        support_c = getattr(v, "hgmd_support_confirmed", 0) or 0
        db_sum_parts.append(f"HGMD: {v.hgmd_class} ({pubs} pubs, dmsupported={dmsup}, confirm={support_c})")
    if v.internal_db_sig:
        db_sum_parts.append(f"Internal: {v.internal_db_sig}")
    db_summary = "; ".join(db_sum_parts) if db_sum_parts else "No database evidence"

    gnomad_ver = v.gnomad_details.get("version", "Unknown") if v.gnomad_details else "Unknown"
    in_auto = getattr(v, '_in_auto', False)
    in_manual = getattr(v, '_in_manual', False)
    filtered_flag = getattr(v, '_filtered', True)
    sample_af_val = f"{v.proband_af:.3f}" if v.proband_af is not None else "0.0"
    # choose effective MOI for this specific variant (prefer disease-specific mapping)
    try:
        moi = select_moi_for_variant(rule if isinstance(rule, dict) else {}, v)
        if not moi:
            moi = rule.get("moi", "Unknown") if isinstance(rule, dict) else "Unknown"
    except Exception:
        moi = rule.get("moi", "Unknown") if isinstance(rule, dict) else "Unknown"
    ps1_str = ";".join([f"{m.get('hgvsp','')}/{m.get('clnsig','')}/{m.get('source','')}" for m in v.ps1_matches]) if v.ps1_matches else ""
    pm5_str = ";".join([f"{m.get('hgvsp','')}/{m.get('best','')}" for m in v.pm5_matches]) if v.pm5_matches else ""

    # dmsupported / support_confirmed values
    dmsupported_val = getattr(v, "hgmd_dmsupported", 0) or 0
    support_confirmed_val = getattr(v, "hgmd_support_confirmed", 0) or 0

    # disease match meta
    disease_match = getattr(v, '_disease_match', False)
    disease_match_level = getattr(v, '_disease_match_level', 'N/A')
    disease_match_reason = getattr(v, '_disease_match_reason', 'N/A')

    # comments: default empty; for auto variants we will move informational manual_reasons here later (in triage)
    comments_field = ""

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
        "AA_ref": v.aa_ref or "",
        "AA_alt": v.aa_alt or "",
        "AA_pos": v.aa_pos or "",
        "Exon": v.exon,
        "Consequence": v.consequence,
        "Zygosity": zyg,
        "DP": v.proband_dp if v.proband_dp is not None else "0",
        "Sample_AF": sample_af_val,
        "acmg_disease_association": acmg_disease,
        "Special_Guidance": has_guidance,
        "disease_match": "Yes" if disease_match else "No",
        "disease_match_level": disease_match_level,
        "disease_match_reason": disease_match_reason,
        "clinvar_trait": clinvar_phen,
        "clinvar_sig": v.clinvar_sig or "Not provided",
        "clinvar_stars": v.clinvar_stars,
        "hgmd_class": v.hgmd_class or "Not provided",
        "hgmd_phenotype": hgmd_phen,
        "hgmd_publication_count": v.hgmd_publication_count,
        "hgmd_dmsupported": dmsupported_val,
        "hgmd_support_confirmed": support_confirmed_val,
        "hgmd_rankscore": getattr(v, 'hgmd_rankscore', "Not provided"),
        "internal_db_sig": v.internal_db_sig or "Not provided",
        "database_summary": db_summary,
        "gnomAD_AF": f"{v.gnomad_af:.6f}" if v.gnomad_af else "0.0",
        "gnomAD_version": gnomad_ver,
        "REVEL": revel_val,
        "CADD": cadd_val,
        "SpliceAI": splice_val,
        "AlphaMissense": alpha_val,
        "ps1_matches": ps1_str,
        "pm5_matches": pm5_str,
        "criteria_assigned": ";".join(getattr(v, "criteria_assigned", [])),
        "criteria_points": json.dumps(getattr(v, "criteria_points", {}) or {}, ensure_ascii=False),
        "total_points": v.total_points,
        "automated_class": v.automated_class,
        "phase": v.phase,
        "is_delins_part": "Yes" if v.is_delins_part else "No",
        "in_auto_report": "Yes" if in_auto else "No",
        "in_manual_review": "Yes" if in_manual else "No",
        # placeholder fields - will be filled/overwritten in triage output section
        "auto_report_reasons": ";".join(getattr(v, "auto_report_reasons", [])),
        "manual_reasons": reasons,
        "filtered_reasons": ";".join(getattr(v, "filtered_reasons", [])),
        "filtered_out": "Yes" if filtered_flag else "No",
        "comments": comments_field
    }

    # ensure per-criterion explanation keys present (avoid KeyError on reindex)
    for crit in ["PVS1","PS1","PM1","PM2","PM3","PM4","PM5","PP1","PP3","PP5","BA1","BS1","BS2","BP4","BP5","BP6","BP7"]:
        output_row[f"{crit}_expl"] = v.criteria_explanations.get(crit, "")

    return output_row

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
    cooccur_path: Optional[str] = None
) -> Tuple[str, str, str, int, int]:
    """
    Strict triage and output.

    Replaces the original strict_triage_and_output implementation with:
     - Conservative auto/manual rules described in the spec (HQ DB-driven auto, algorithmic auto,
       AR blocking, AR pair rescue with phasing/cooccurrence checks, XL-male hemizygosity AF check).
     - Uses dbm if provided for cis/trans co-occurrence lookup; if not provided, remains conservative.
     - Augments output rows with per-criterion explanation columns (e.g. PVS1_expl, PS1_expl, ...).
     - Accepts cooccurrence TSV path (cooccur_path) to check population evidence of cis/trans for pairs.

    Notes / Integration:
     - Place this function in the main module (acmg_sf_classifier.py) replacing the previous strict_triage_and_output.
     - Callers should pass dbm and cooccur_path if available, e.g. from process_vcf:
           strict_triage_and_output(..., dbm=dbm, cooccur_path=args.gnomad_cooccurrence)
       If not passed, the function still works but will be conservative for unphased AR pairs.
     - This function expects VariantRecord to have (or will set) these attributes/lists:
         - v.criteria_explanations: Dict[str,str]
         - v.auto_report_reasons: List[str]
         - v.filtered_reasons: List[str]
         - v.ps1_matches, v.pm5_matches : lists
    """

    os.makedirs(outdir, exist_ok=True)
    family_id = family_id_input or "Unknown"
    auto_ids = set() 
    manual_ids = set()

    # ensure DB manager available (prefer passed dbm; fallback to global DBM if present)
    dbm_local = dbm or globals().get("DBM") or globals().get("dbm") or None

    def _select_moi_for_variant(rule: Dict[str, Any], variant: VariantRecord) -> str:
        """
        Select disease-specific MOI for given variant using rule dict and variant DB phenotype fields.
        - rule: gene rule dict (as from load_acmg_table), may contain 'disease' and 'moi' (both possibly with ';' separators)
        - variant: VariantRecord with fields hgmd_phen, clinvar_trait (used to match disease text)
        Returns the matching MOI string (raw), or fallback to rule.get('moi','') if no per-disease match.
        """
        if not rule:
            return ""
        moi_field = str(rule.get("moi", "") or "").strip()
        disease_field = str(rule.get("disease", "") or "").strip()
        if not moi_field or not disease_field:
            return moi_field

        # Split disease and moi lists by common separators preserving order
        disease_list = [d.strip() for d in re.split(r'[;/\|]', disease_field) if d.strip()]
        moi_list = [m.strip() for m in re.split(r'[;/\|]', moi_field) if m.strip()]

        # Pad moi_list to length of disease_list if needed
        if len(moi_list) < len(disease_list):
            moi_list = (moi_list + [""] * len(disease_list))[: len(disease_list)]

        # Build variant phenotype string (HGMD + ClinVar)
        traits = []
        if getattr(variant, "hgmd_phen", "") and str(variant.hgmd_phen).lower() not in ("", "not provided", ".", "null"):
            traits.append(str(variant.hgmd_phen))
        if getattr(variant, "clinvar_trait", "") and str(variant.clinvar_trait).lower() not in ("", "not provided", ".", "null"):
            traits.append(str(variant.clinvar_trait))
        variant_trait = " | ".join(traits).strip()

        # Try to match each disease token individually against the variant trait
        for idx, disease_token in enumerate(disease_list):
            try:
                matched, level, reason = check_disease_match(disease_token, variant_trait, extra_recs_list=None, gene_rules=None)
                if matched:
                    # return corresponding MOI if available
                    return moi_list[idx] if idx < len(moi_list) else moi_list[-1] if moi_list else moi_field
            except Exception:
                continue

        # Fallbacks:
        # 1) if variant._disease_match True (was matched earlier against full rule.disease) — prefer first MOI token that looks non-empty
        try:
            if getattr(variant, "_disease_match", False):
                for m in moi_list:
                    if m:
                        return m
        except Exception:
            pass

        # 2) return full MOI field as fallback
        return moi_field

    # Helper: unique key for a variant
    def get_key(v: VariantRecord) -> str:
        return f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}:{v.sample}"

    # Helper: safe allele fraction extraction
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

    # Helper: read cooccurrence TSV into dict for fast lookup (load once)
    coocc_dict = {}
    if cooccur_path and os.path.exists(cooccur_path):
        try:
            with open(cooccur_path, "rt", encoding="utf-8") as fh:
                for ln in fh:
                    if not ln.strip() or ln.startswith("#"):
                        continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    a = parts[0].strip()
                    b = parts[1].strip()
                    cis_flag = parts[2].strip().lower()
                    evidence = parts[3].strip() if len(parts) > 3 else ""
                    key1 = f"{a}|{b}"
                    key2 = f"{b}|{a}"
                    cis_val = None
                    if cis_flag in ("yes", "y", "true", "1"):
                        cis_val = True
                    elif cis_flag in ("no", "n", "false", "0"):
                        cis_val = False
                    coocc_dict[key1] = {"cis": cis_val, "evidence": evidence, "source": "cooccurrence_tsv"}
                    coocc_dict[key2] = coocc_dict[key1]
        except Exception:
            # if loading fails, keep dict empty (function remains conservative)
            coocc_dict = {}

    # Helper: check cis/trans evidence from recs_map or cooccurrence dict
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
        # 2) cooccurrence TSV dict
        key = f"{vid1}|{vid2}"
        if key in coocc_dict:
            return coocc_dict[key]
        # 3) if DatabaseManager provides such lookup (schema dependent) try conservative attempt
        try:
            if dbm_local and hasattr(dbm_local, "internal_db_index") and dbm_local.internal_db_index:
                # No standard way; skip (conservative)
                pass
        except Exception:
            pass
        return None
    
    def _db_high_quality_consensus(v: VariantRecord) -> Tuple[Optional[str], List[Dict[str,str]], bool, Optional[str]]:
        """
        Inspect variant v and return:
        (consensus_class, sources_list, consensus_flag, conflict_detail)

        - consensus_class: 'Pathogenic' or 'Likely pathogenic' (or None)
        - sources_list: list of dicts {"source":..., "raw":..., "norm":...}
        - consensus_flag: True if we may auto-assign (consensus or mild P/LP conflict)
        - conflict_detail: None, or text describing a substantive conflict (if consensus_flag False),
                        or a short note describing a mild P/LP disagreement (if consensus_flag True but mixed P/LP)

        High-quality sources considered:
        - ClinVar with stars >= 2 (use clinvar_sig)
        - HGMD DM with pubs >= 2 or dmsupported >=2 (treated as Pathogenic)
        - Internal DB annotated as Pathogenic/Likely pathogenic
        """
        sources: List[Dict[str,str]] = []

        # Helper to normalise significance strings
        def _norm_sig(s: Optional[str]) -> str:
            if not s:
                return "Unknown"
            s2 = str(s).lower()
            if "pathogenic" in s2 and "likely" not in s2:
                return "Pathogenic"
            if "likely pathogenic" in s2 or ("likely" in s2 and "pathogen" in s2):
                return "Likely pathogenic"
            if "benign" in s2:
                return "Benign"
            if "conflict" in s2 or "conflicting" in s2:
                return "Conflict"
            return s.strip() or "Unknown"

        # ClinVar high-quality entries (>=2 stars)
        try:
            if getattr(v, "clinvar_stars", "") and str(v.clinvar_stars).isdigit() and int(v.clinvar_stars) >= 2:
                raw = getattr(v, "clinvar_sig", "") or ""
                norm = _norm_sig(raw)
                sources.append({"source": "ClinVar", "raw": raw, "norm": norm})
                # if clinvar reports explicit 'conflict' mark it
                if "conflict" in (raw or "").lower():
                    # treat as internal conflict indicator
                    pass
        except Exception:
            pass

        # HGMD DM with evidence
        try:
            hgmd_tag = (getattr(v, "hgmd_class", "") or "").strip()
            if hgmd_tag:
                pubs = int(getattr(v, "hgmd_publication_count", 0) or 0)
                dmsup = int(getattr(v, "hgmd_dmsupported", 0) or 0)
                if hgmd_tag.upper() == "DM" and (pubs >= 2 or dmsup >= 2):
                    sources.append({"source": "HGMD", "raw": hgmd_tag, "norm": "Pathogenic"})
                elif hgmd_tag:
                    # still record lower-support HGMD entries as non-HQ (but include in sources list)
                    sources.append({"source": "HGMD", "raw": hgmd_tag, "norm": _norm_sig(hgmd_tag)})
        except Exception:
            pass

        # Internal DB
        try:
            internal_raw = (getattr(v, "internal_db_sig", "") or "").strip()
            if internal_raw:
                sources.append({"source": "Internal_DB_display", "raw": internal_raw, "norm": _norm_sig(internal_raw)})
        except Exception:
            pass

        if not sources:
            return None, [], False, None

        # Build set of normalized labels for HQ sources only
        # Consider as HQ those entries that are ClinVar(>=2★), HGMD DM with pubs/dmsup>=2, Internal P/LP
        hq_norms = []
        hq_sources = []
        for s in sources:
            src = s["source"]
            norm = s["norm"]
            raw = s["raw"]
            # decide HQ status
            is_hq = False
            if src == "ClinVar":
                # we already ensured only clinvar >=2 stars were added as ClinVar entries above
                is_hq = True
            elif src == "HGMD":
                # if we recorded as Pathogenic from HGMD it's HQ; else not HQ
                if norm == "Pathogenic":
                    is_hq = True
            elif src == "Internal_DB":
                # if norm in ("Pathogenic", "Likely pathogenic"):
                #     is_hq = True
                is_hq = False 
            if is_hq:
                hq_sources.append(s)
                hq_norms.append(norm)

        # If no HQ sources, nothing to auto-decide here
        if not hq_sources:
            # not HQ consensus; return full sources for logging but not consensus
            return None, sources, False, None

        unique = set(hq_norms)
        # Strict unanimous consensus
        if len(unique) == 1:
            return list(unique)[0], sources, True, None

        # Mild P vs LP disagreement -> treat as mild conflict but allow auto (choose conservative: Pathogenic)
        if unique == {"Pathogenic", "Likely pathogenic"} or unique == {"Likely pathogenic", "Pathogenic"}:
            detail = "Mild disagreement among HQ sources (Pathogenic vs Likely pathogenic); assigning conservative class 'Pathogenic' and adding note."
            # conservative assignment: Pathogenic (stronger)
            return "Pathogenic", sources, True, detail

        # Otherwise this is a substantive conflict (e.g. Pathogenic vs Benign or explicit ClinVar conflict)
        # Build detailed conflict string
        parts = []
        for s in hq_sources:
            parts.append(f"{s['source']}='{s['raw']}'")
        conflict_detail = "Conflict among high-quality sources: " + "; ".join(parts)
        # Also include any explicit ClinVar conflict flag if present
        try:
            clin_raw = getattr(v, "clinvar_sig", "") or ""
            if "conflict" in clin_raw.lower():
                conflict_detail += " (ClinVar reports conflict)"
        except Exception:
            pass

        return None, sources, False, conflict_detail

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
            if not varlist:
                continue

            rule = gene_rules.get(gene, {}) if gene_rules else {}
            spec_rules = recs_map.get(gene, {}) if recs_map else {}
            moi = str(rule.get("moi","")).upper() if rule else ""
            acmg_disease = rule.get("disease", "") if rule else ""
            sample_sex = varlist[0].sample_sex if varlist and hasattr(varlist[0], "sample_sex") else "Unknown"

            # inheritance flags
            is_ar = ("AR" in moi) or ("RECESSIVE" in moi)
            is_xlr = ("XLR" in moi) or ("X-LINKED RECESSIVE" in moi)
            is_xl_male = (("X" in moi) or ("X-LINKED" in moi) or is_xlr) and (sample_sex == "Male")

            # per-variant pre-filtering and checks

            biologically_valid: List[VariantRecord] = []
            for v in varlist:
                v._in_auto = False
                v._in_manual = False
                v._filtered = True  # default filtered until explicitly unfiltered
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
                if getattr(v, "clinvar_trait", "") and v.clinvar_trait.lower() not in ("not provided", ".", "not_provided", "", "null"):
                    traits.append(v.clinvar_trait)

                full_trait = " | ".join(traits) if traits else ""

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

                # --- Special case: High-quality ClinVar VUS should be treated conservatively ---
                # If ClinVar has >=2 stars and significance is Uncertain_significance / VUS,
                # treat as HQ-VUS: assign automated_class="VUS" and DO NOT put into auto/manual lists.
                try:
                    clin_stars_val = getattr(v, "clinvar_stars", "")
                    clin_sig_val = str(getattr(v, "clinvar_sig", "") or "").lower()
                    if str(clin_stars_val).isdigit() and int(clin_stars_val) >= 2 and ("uncertain" in clin_sig_val or "vus" in clin_sig_val):
                        # set class to VUS and keep variant filtered (not auto, not manual)
                        v.automated_class = "VUS"
                        v._in_auto = False
                        v._in_manual = False
                        v._filtered = True
                        note = "Filtered: ClinVar ≥2★ with Uncertain_significance (HQ-VUS) — conservative non-report"
                        if note not in v.filtered_reasons:
                            v.filtered_reasons.append(note)
                        # do not add to biologically_valid set for reporting
                        continue
                except Exception:
                    # conservative: if parsing fails, proceed as normal
                    pass

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
                    # check per-variant criterion: each variant must be HQ P/LP OR total_points > 5
                    all_meet = True
                    for m in members:
                        meets = False
                        try:
                            if getattr(m, "is_hq_pathogenic_db", False):
                                meets = True
                            elif getattr(m, "automated_class", "") in ("Pathogenic", "Likely pathogenic") and getattr(m, "total_points", 0) > 5:
                                # If automated_class is P/LP and total_points > 5 -> meets
                                meets = True
                            elif getattr(m, "total_points", 0) > 5:
                                meets = True
                        except Exception:
                            meets = False
                        if not meets:
                            all_meet = False
                            break
                    if all_meet:
                        # Move all cluster members to manual review
                        member_keys = [f"{x.chrom}:{x.pos}:{x.ref}:{x.alt}" for x in members]
                        for m in members:
                            k = get_key(m)
                            if k not in manual_ids:
                                m.manual_reasons.append(f"Manual: Delins/overlap cluster (id={cid}) meets criteria -> cluster manual review required. Members: {', '.join(member_keys)}")
                                m._in_manual = True
                                m._filtered = False
                                m._in_auto = False
                                manual_ids.add(k)
                                if k in auto_ids:
                                    auto_ids.remove(k)
                    else:
                        # Add an informational note (will appear in comments for auto variants)
                        for m in members:
                            note = f"Delins/overlap cluster (id={cid}) detected but not all members meet P/LP or >5pts -> not auto-escalated (possible carrier)."
                            if note not in m.manual_reasons:
                                m.manual_reasons.append(note)
            except Exception:
                # conservative: if evaluation fails, do not change manual/auto status
                pass
            # --- AUTOMATIC LOGIC ---
            for v in drivers:
                if get_key(v) in manual_ids:
                    continue

                # ensure lists exist
                if not hasattr(v, "auto_report_reasons"):
                    v.auto_report_reasons = []
                if not hasattr(v, "criteria_explanations"):
                    v.criteria_explanations = {}

                # --- Conservative AR blocking block with disease-specific MOI check ---
                # Determine disease-specific MOI for this variant
                eff_moi = _select_moi_for_variant(rule, v).upper()
                eff_is_ar = ("AR" in eff_moi) or ("RECESSIVE" in eff_moi)
                eff_is_xlr = ("XLR" in eff_moi) or ("X-LINKED RECESSIVE" in eff_moi)
                eff_is_xl_male = (("X" in eff_moi) or ("X-LINKED" in eff_moi) or eff_is_xlr) and (sample_sex == "Male")

                # If the variant maps to a recessive disease (by disease match), then block heterozygotes until pairing confirmed.
                # If the variant maps to a dominant disease (e.g. breast cancer AD) eff_is_ar will be False and it will not be auto-blocked.
                if eff_is_ar and not eff_is_xl_male:
                    if not _is_hom(v.proband_gt):
                        reason = "Filtered: Recessive heterozygous variant — requires pairing/phasing before reporting"
                        if reason not in v.filtered_reasons:
                            v.filtered_reasons.append(reason)
                        v._filtered = True
                        v._in_auto = False
                        v._in_manual = False
                        auto_ids.discard(get_key(v))
                        continue

                # Case A: HQ DB-driven automatic (reworked to ignore single-DB gene-mismatch and to defer AR-het manualization)
                if getattr(v, "is_hq_pathogenic_db", False):
                    # Detect per-source gene-mismatch flags recorded earlier in v.manual_reasons
                    clinvar_mismatch = any("clinvar gene mismatch" in (r or "").lower() for r in v.manual_reasons)
                    hgmd_mismatch = any("hgmd gene mismatch" in (r or "").lower() for r in v.manual_reasons)
                    internal_mismatch = any("internal db" in (r or "").lower() and "gene mismatch" in (r or "").lower() for r in v.manual_reasons)

                    # Determine HQ consensus and raw sources (original)
                    db_class, db_sources, consensus_flag, conflict_detail = _db_high_quality_consensus(v)

                    # If there were gene-mismatch notes, filter out corresponding DB sources and re-evaluate consensus
                    if (clinvar_mismatch or hgmd_mismatch or internal_mismatch) and db_sources:
                        filtered_sources = []
                        for s in db_sources:
                            sname = (s.get("source") or "").lower()
                            # skip the source if it was flagged as gene-mismatch
                            if clinvar_mismatch and "clinvar" in sname:
                                continue
                            if hgmd_mismatch and "hgmd" in sname:
                                continue
                            if internal_mismatch and ("internal" in sname or "internal_db" in sname or "internal_db_display" in sname):
                                continue
                            filtered_sources.append(s)

                        # Recompute HQ consensus from filtered_sources using same HQ logic as _db_high_quality_consensus
                        # Only consider ClinVar (already HQ-filtered by stars) and HGMD entries with norm == "Pathogenic" as HQ
                        hq_norms = []
                        hq_sources_filtered = []
                        for s in filtered_sources:
                            src = s.get("source", "")
                            norm = s.get("norm", "")
                            # Mark HQ: ClinVar entries included, and HGMD entries normalized as Pathogenic
                            if "clinvar" in (src or "").lower():
                                hq_sources_filtered.append(s); hq_norms.append(norm)
                            elif "hgmd" in (src or "").lower() and str(norm).lower().startswith("pathogen"):
                                hq_sources_filtered.append(s); hq_norms.append("Pathogenic")
                            # Internal DB is not considered HQ in original logic unless explicitly desired -> skip as HQ

                        # derive new consensus info
                        if not hq_sources_filtered:
                            # no HQ sources remain after excluding mismatched DB -> do not auto-consensus, fall back to algorithmic scoring
                            consensus_flag = False
                            db_class = None
                            conflict_detail = None
                        else:
                            unique = set([str(x).strip() for x in hq_norms if x])
                            if len(unique) == 1:
                                db_class = list(unique)[0]
                                consensus_flag = True
                                conflict_detail = None
                            elif unique == {"Pathogenic", "Likely pathogenic"} or unique == {"Likely pathogenic", "Pathogenic"}:
                                # mild disagreement -> conservative assignment to Pathogenic (as earlier)
                                db_class = "Pathogenic"
                                consensus_flag = True
                                conflict_detail = "Mild disagreement among HQ sources (Pathogenic vs Likely pathogenic) after excluding mismatched DB"
                            else:
                                # substantive conflict among remaining HQ sources -> manual
                                db_class = None
                                consensus_flag = False
                                conflict_detail = "Conflict among remaining high-quality sources after excluding mismatched DB: " + "; ".join([f"{s.get('source')}='{s.get('raw')}'" for s in hq_sources_filtered])

                    # After possible filtering & recomputation: act only if consensus_flag True
                    if consensus_flag:
                        # Determine effective MOI for gene-level rule (use rule-level MOI as fallback)
                        v_moi = str(rule.get("moi", "") or "").upper()

                        # If MOI suggests recessive and variant is heterozygous, DO NOT immediately send to manual:
                        # first defer to phasing/pairing checks (will be handled later in AR-pairing stage).
                        # Add informative note and proceed (do not set _in_manual here).
                        if ("AR" in v_moi or "RECESSIVE" in v_moi) and not _is_hom(v.proband_gt):
                            note = "Recessive inheritance suggested; heterozygous P/LP — pairing/phasing review required before reporting."
                            if note not in v.manual_reasons:
                                v.manual_reasons.append(note)
                            # leave v._in_auto False for now; allow algorithmic path and AR pairing checks to consider this variant
                            # do not continue / do not force manual assignment here
                        else:
                            # For dominant / non-recessive MOI (or homozygotes), auto-assign per DB consensus
                            v._in_auto = True
                            v._filtered = False
                            v.auto_report_reasons.append("HQ DB consensus")
                            if conflict_detail:
                                v.auto_report_reasons.append(f"HQ note: {conflict_detail}")
                            auto_ids.add(get_key(v))
                            # set automated_class to DB consensus if available
                            if db_class:
                                v.automated_class = db_class
                            # Ensure PP5_Supporting for HQ evidence (if not already present)
                            if not getattr(v, "criteria_points", None):
                                v.criteria_points = {}
                            if "PP5_Supporting" not in v.criteria_points:
                                v.criteria_points["PP5_Supporting"] = POINTS_MAP.get("PP5_Supporting", 1)
                            try:
                                v.total_points = sum(int(x) for x in v.criteria_points.values())
                            except Exception:
                                v.total_points = getattr(v, "total_points", 0)


                # Case B: Algorithmic scoring (Tavtigian)
                if ((v.automated_class == "Pathogenic" and v.total_points >= 10) or
                    (v.automated_class == "Likely pathogenic" and v.total_points >= 6)):
                    if getattr(v, "is_hq_conflict_db", False):
                        # conflicts must be manual-reviewed
                        continue
                    # require disease match or no DB phenotype or no DB evidence
                    if v._disease_match or (v._disease_match_reason == "no_db_phenotype") or (not v._has_db_evidence):
                        v._in_auto = True
                        v._filtered = False
                        v.auto_report_reasons.append("Algorithmic score")
                        auto_ids.add(get_key(v))
                        continue

            # --- MANUAL LOGIC ---
            for v in drivers:
                if get_key(v) in manual_ids:
                    continue

                # 1) HQ conflict between HQ sources -> manual
                if getattr(v, "is_hq_pathogenic_db", False) and getattr(v, "is_hq_conflict_db", False):
                    v.manual_reasons.append(f"HQ Conflict: {getattr(v,'conflict_details','')}")
                    v._in_manual = True; v._filtered = False; v._in_auto = False
                    manual_ids.add(get_key(v))
                    if get_key(v) in auto_ids:
                        auto_ids.remove(get_key(v))
                    continue

                # 2) HQ DB pathogenic but explicit disease mismatch -> manual
                if getattr(v, "is_hq_pathogenic_db", False) and not v._disease_match:
                    if v._disease_match_reason == "mismatch":
                        v.manual_reasons.append("HQ DB pathogenic but disease annotation mismatches ACMG disease -> Manual")
                        v._in_manual = True; v._filtered = False; v._in_auto = False
                        manual_ids.add(get_key(v))
                        if get_key(v) in auto_ids: auto_ids.remove(get_key(v))
                        continue
                    # if reason == no_db_phenotype -> leave for algorithmic path (not forced manual)

            # --- XL male heterozygous handling ---
            if is_xl_male:
                for v in drivers:
                    if _is_het(v.proband_gt):
                        af_val = _get_sample_af(v)
                        # treat hemizygous if AF >= 0.9 -> do NOT send to manual on this ground alone
                        if af_val <= 0.9:
                            if get_key(v) not in manual_ids:
                                v.manual_reasons.append(f"XL male heterozygous with AF={af_val:.3f} -> Manual for hemizygosity check")
                                v._in_manual = True; v._filtered = False; v._in_auto = False
                                manual_ids.add(get_key(v))
                                if get_key(v) in auto_ids: auto_ids.remove(get_key(v))
                        else:
                            # treat as hemizygous; leave to auto/manual decisions based on other evidence
                            pass

            # --- Complex recessive: >=2 P/LP variants => manual review for all P/LP hets
            ar_path_lp = []
            for x in drivers:
                # only consider variants that match a recessive disease according to per-variant MOI selection
                eff_moi_x = _select_moi_for_variant(rule, x).upper()
                if (("AR" in eff_moi_x) or ("RECESSIVE" in eff_moi_x)) and _is_het(x.proband_gt) and x.automated_class in ("Pathogenic", "Likely pathogenic"):
                    ar_path_lp.append(x)
            if len(ar_path_lp) >= 2:
                for topv in ar_path_lp:
                        if get_key(topv) not in manual_ids:
                            topv.manual_reasons.append("Manual: Recessive gene with multiple heterozygous P/LP variants - phasing and pairing review required")
                            topv._in_manual = True; topv._filtered = False; topv._in_auto = False
                            manual_ids.add(get_key(topv))
                            if get_key(topv) in auto_ids:
                                auto_ids.remove(get_key(topv))

            # --- AR Pair Rescue ---
            # Build list of heterozygous variants in this gene/sample that specifically map to a recessive disease MOI
            passengers_hets = [p for p in passengers if _is_het(p.proband_gt) and (("AR" in _select_moi_for_variant(rule, p).upper()) or ("RECESSIVE" in _select_moi_for_variant(rule, p).upper()))]

            # If none, skip AR pairing logic
            if passengers_hets:
                # P/LP heterozygotes as primary candidates
                plp_hets = [p for p in passengers_hets if p.automated_class in ("Pathogenic", "Likely pathogenic")]
                # For each P/LP, evaluate potential VUS partners
                def _phased_trans(gt1: str, gt2: str) -> Optional[bool]:
                    """
                    Return True if phased and alleles indicate trans (alts on different haplotypes),
                    False if phased and indicate cis, None if not phased / cannot decide.
                    Uses left-side allele comparison: if left allele differs -> trans.
                    Assumes simple diploid GTs like '0|1' or '1|0'.
                    """
                    try:
                        if not gt1 or not gt2:
                            return None
                        if "|" not in gt1 or "|" not in gt2:
                            return None
                        a1_left = gt1.split("|")[0].strip()
                        a2_left = gt2.split("|")[0].strip()
                        if a1_left == "." or a2_left == ".":
                            return None
                        # If left alleles equal => cis (alts on same haplotype), else trans
                        return a1_left != a2_left
                    except Exception:
                        return None

                # helper to confirm trans using phased GT, parental evidence, or cooccurrence
                def _confirm_trans(v1: VariantRecord, v2: VariantRecord) -> bool:
                    # 1) GT phased check
                    ph = _phased_trans(v1.proband_gt, v2.proband_gt)
                    if ph is True:
                        return True
                    if ph is False:
                        return False
                    # 2) parental phasing check (existing pattern)
                    try:
                        # inherited from different parents -> in trans
                        if (_is_het(v1.father_gt) and _is_wt(v1.mother_gt) and
                            _is_het(v2.mother_gt) and _is_wt(v2.father_gt)):
                            return True
                        if (_is_het(v1.mother_gt) and _is_wt(v1.father_gt) and
                            _is_het(v2.father_gt) and _is_wt(v2.mother_gt)):
                            return True
                    except Exception:
                        pass
                    # 3) cooccurrence TSV / recs_map evidence
                    try:
                        coinfo = _check_pair_cis_evidence(v1, v2)
                        if coinfo and "cis" in coinfo:
                            # cis==False means observed in trans
                            if coinfo.get("cis") is False:
                                return True
                            # cis==True means observed in cis => not trans
                            if coinfo.get("cis") is True:
                                return False
                    except Exception:
                        pass
                    # could not confirm trans
                    return False

                # If there are no P/LP heterozygotes, nothing to do
                if not plp_hets:
                    # no primary P/LP; nothing to rescue
                    pass
                else:
                    # For each P/LP variant, find candidate VUS partners (hets)
                    for primary in sorted(plp_hets, key=lambda x: x.total_points, reverse=True):
                        # Gather heterozygous partners excluding primary
                        partners = [p for p in passengers_hets if p is not primary]
                        if not partners:
                            continue

                        # Partition partners into VUS with >=5, VUS <=4, and other P/LP
                        vus5 = [p for p in partners if p.automated_class == "VUS" and (getattr(p, "total_points", 0) >= 5)]
                        vus_le4 = [p for p in partners if p.automated_class == "VUS" and (getattr(p, "total_points", 0) <= 4)]
                        other_plp_partners = [p for p in partners if p.automated_class in ("Pathogenic", "Likely pathogenic")]

                        # CASE: Multiple VUS with >=5 points -> send all related variants (primary + all vus5) to manual.
                        if len(vus5) >= 2:
                            for pvar in ([primary] + vus5):
                                if get_key(pvar) not in manual_ids:
                                    pvar.manual_reasons.append("Manual: Recessive gene with multiple heterozygous VUS (>=5 points) together with P/LP -> phasing/pairing review required")
                                    pvar._in_manual = True; pvar._filtered = False; pvar._in_auto = False
                                    manual_ids.add(get_key(pvar))
                                    if get_key(pvar) in auto_ids:
                                        auto_ids.remove(get_key(pvar))
                            # we handled this primary, move to next
                            continue

                        # CASE: Single VUS with >=5 points: require confirmed trans to send to manual
                        if len(vus5) == 1:
                            partner = vus5[0]
                            in_trans_confirmed = _confirm_trans(primary, partner)
                            if in_trans_confirmed:
                                # send both to manual for pairing review
                                for pvar in (primary, partner):
                                    if get_key(pvar) not in manual_ids:
                                        pvar.manual_reasons.append("Manual: Recessive pair confirmed in trans: P/LP + VUS(>=5) -> phasing/pairing review required")
                                        pvar._in_manual = True; pvar._filtered = False; pvar._in_auto = False
                                        manual_ids.add(get_key(pvar))
                                        if get_key(pvar) in auto_ids:
                                            auto_ids.remove(get_key(pvar))
                            else:
                                # not confirmed trans -> conservative: do not report (carrier); do not send anywhere
                                # add an informational note (will end up in comments for auto variants)
                                note = "AR pairing candidate: P/LP + single VUS(>=5) found but not confirmed in trans -> conservative: not reported (carrier)."
                                if note not in primary.manual_reasons:
                                    primary.manual_reasons.append(note)
                                if note not in partner.manual_reasons:
                                    partner.manual_reasons.append(note)
                            continue

                        # CASE: No VUS>=5 present:
                        # If only VUS <=4 exist alongside P/LP, treat as carrier - do not send to manual.
                        if vus_le4 and not vus5:
                            # do nothing (carrier). Add explanatory note to primary (and optionally to partners)
                            note = "AR pairing candidate: P/LP with VUS(s) <=4 points -> treated as carrier (not returned for secondary findings)."
                            if note not in primary.manual_reasons:
                                primary.manual_reasons.append(note)
                            for pv in vus_le4:
                                if note not in pv.manual_reasons:
                                    pv.manual_reasons.append(note)
                            # ensure not moved to manual
                            continue

                        # CASE: Other combinations (e.g., other P/LP partners, mixed classes)
                        # Fall back to previous conservative behavior: attempt to confirm trans and if confirmed send to manual for pairing review.
                        for partner in partners:
                            # skip if partner is the same as primary
                            if partner is primary:
                                continue
                            in_trans_confirmed = _confirm_trans(primary, partner)
                            if in_trans_confirmed:
                                # If partner is VUS with points >=1 we may still want manual review
                                # Conservative rule: confirmed trans with P/LP and any other (VUS>=1) -> manual pairing review
                                if partner.automated_class in ("Pathogenic", "Likely pathogenic"):
                                    # two P/LP in trans -> manual
                                    for pvar in (primary, partner):
                                        if get_key(pvar) not in manual_ids:
                                            pvar.manual_reasons.append("Manual: Two P/LP variants confirmed in trans -> pairing review required")
                                            pvar._in_manual = True; pvar._filtered = False; pvar._in_auto = False
                                            manual_ids.add(get_key(pvar))
                                            if get_key(pvar) in auto_ids:
                                                auto_ids.remove(get_key(pvar))
                                else:
                                    # partner is VUS or other; require at least 1 point and then send to manual
                                    if getattr(partner, "total_points", 0) >= 1:
                                        for pvar in (primary, partner):
                                            if get_key(pvar) not in manual_ids:
                                                pvar.manual_reasons.append("Manual: Recessive pair confirmed in trans -> phasing/pairing review required")
                                                pvar._in_manual = True; pvar._filtered = False; pvar._in_auto = False
                                                manual_ids.add(get_key(pvar))
                                                if get_key(pvar) in auto_ids:
                                                    auto_ids.remove(get_key(pvar))
                                # break after handling this partner for this primary
                                break
                            # if not confirmed -> continue trying other partners

    # --- Build output rows (augment with explanations) ---
    # Define columns including per-criterion explanation columns
    COLUMNS = [
        "sample","Family_ID","chrom","pos","ref","alt","gene","Inheritance",
        "HGVSc","HGVSp","AA_ref","AA_alt","AA_pos","Exon","Consequence","Zygosity","DP","Sample_AF",
        "acmg_disease_association","Special_Guidance",
        "disease_match","disease_match_level","disease_match_reason",
        "clinvar_trait","clinvar_sig","clinvar_stars",
        "hgmd_class","hgmd_phenotype","hgmd_publication_count","hgmd_rankscore", "hgmd_dmsupported",
        "hgmd_support_confirmed",
        "internal_db_sig","database_summary",
        "gnomAD_AF","gnomAD_version",
        "REVEL","CADD","SpliceAI","AlphaMissense",
        "ps1_matches","pm5_matches","phase","is_delins_part",
        "criteria_assigned","criteria_points","total_points","automated_class",
        # per-criterion explanation columns
        "PVS1_expl","PS1_expl","PM1_expl","PM2_expl","PM3_expl","PM4_expl","PM5_expl",
        "PP1_expl","PP3_expl","PP5_expl",
        "BA1_expl","BS1_expl","BS2_expl","BP4_expl","BP5_expl","BP6_expl","BP7_expl",
        # reporting/audit
        "in_auto_report","in_manual_review","auto_report_reasons","manual_reasons","comments","filtered_reasons","filtered_out"
    ]

    processed = set()
    rows_auto = []
    rows_manual = []
    rows_all = []

    for v in all_variants:
        k = get_key(v)
        if k in processed:
            continue
        processed.add(k)

        rule = gene_rules.get(v.gene, {}) if gene_rules else {}
        row = variant_row_with_disease(v, rule, family_id, recs_map)

        # attach per-criterion explanations (already included in variant_row_with_disease but keep safe)
        for crit in ("PVS1","PS1","PM1","PM2","PM3","PM4","PM5","PP1","PP3","PP5",
                     "BA1","BS1","BS2","BP4","BP5","BP6","BP7"):
            expl_col = f"{crit}_expl"
            try:
                row[expl_col] = v.criteria_explanations.get(crit, "")
            except Exception:
                row[expl_col] = ""

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
        if getattr(v, "_in_auto", False):
            # Auto: informational notes go into comments; manual_reasons cleared
            row["comments"] = ";".join(getattr(v, "manual_reasons", [])) if getattr(v, "manual_reasons", None) else ""
            row["manual_reasons"] = ""
        elif getattr(v, "_in_manual", False):
            # Manual: keep manual_reasons in the manual_reasons column, comments empty
            row["comments"] = ""
            row["manual_reasons"] = ";".join(getattr(v, "manual_reasons", []))
        else:
            # Neither auto nor manual: move notes into comments (manual_reasons empty)
            row["comments"] = ";".join(getattr(v, "manual_reasons", [])) if getattr(v, "manual_reasons", None) else ""
            row["manual_reasons"] = ""

        # filtered reasons and filtered_out flag
        row["filtered_reasons"] = ";".join(getattr(v, "filtered_reasons", []))
        # final filtered_out: if variant was selected for AUTO or MANUAL reporting -> it is not filtered out (No)
        # otherwise it is filtered out (Yes). Это более интуитивное и устойчивое поведение, чем использование промежуточного _filtered.
        if getattr(v, "_in_auto", False) or getattr(v, "_in_manual", False):
            row["filtered_out"] = "No"
        else:
            row["filtered_out"] = "Yes"

        if getattr(v, "_in_auto", False):
            rows_auto.append(row)
        if getattr(v, "_in_manual", False):
            rows_manual.append(row)
        rows_all.append(row)

    # collapse duplicates
    auto_rows_c = collapse_output_rows(rows_auto)
    manual_rows_c = collapse_output_rows(rows_manual)
    all_rows_c = collapse_output_rows(rows_all)

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
        We follow your AR pairing logic: take heterozygous P/LP as potential primary and choose partner candidate:
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
# Build variant record from VCF record (robust)
# ---------------------------
def build_variant_record_from_rec(rec, csq_header: List[str], acmg_genes_normalized: set,
                                  dad_vf: Optional[pysam.VariantFile],
                                  mom_vf: Optional[pysam.VariantFile],
                                  exons_df: Optional[Any],
                                  proband_sample: Optional[str]=None,
                                  require_mane: bool=False,
                                  dbm: Optional[DatabaseManager]=None) -> Optional[VariantRecord]:
    try:
        if hasattr(rec, "CHROM"):
            chrom = rec.CHROM; pos = rec.POS; ref = rec.REF; alts = rec.ALT
        else:
            chrom = rec.chrom; pos = rec.pos; ref = rec.ref; alts = rec.alts
        alt = alts[0] if alts and len(alts) > 0 else ""
        if not alt: return None
        chrom_norm = normalize_chrom(chrom)
    except Exception:
        return None
    info = {}
    if hasattr(rec, "INFO"): info = dict(rec.INFO)
    elif hasattr(rec, "info"):
        try: info = dict(rec.info.items())
        except: info = dict(rec.info)
    info_l = {k.lower(): v for k, v in info.items()}
    csq_lines = []
    csq_raw = info_l.get("csq") or info_l.get("ann")
    if csq_raw:
        s = csq_raw
        if isinstance(s, (list, tuple)): s = ",".join(str(x) for x in s)
        for part in str(s).split(","):
            if "|" in part: csq_lines.append(part)
    best_ann = None
    best_rank = -1
    for ln in csq_lines:
        ann = parse_vep_csq_line(ln, csq_header)
        gene_raw = (ann.get("SYMBOL") or ann.get("Gene") or "").strip()
        norm_gene = normalize_gene_name(gene_raw)
        mpc = str(ann.get("MANE_PLUS_CLINICAL") or "").upper()
        ms = str(ann.get("MANE_SELECT") or "").upper()
        canon = str(ann.get("CANONICAL") or "").upper()
        rank = 0
        if mpc in ("YES","Y","TRUE","1"): rank += 100
        if ms in ("YES","Y","TRUE","1"): rank += 50
        if canon in ("YES","Y","TRUE","1"): rank += 10
        if norm_gene in acmg_genes_normalized: rank += 5
        if require_mane and not (mpc in ("YES","Y","TRUE","1") or ms in ("YES","Y","TRUE","1")):
            continue
        if rank > best_rank:
            best_rank = rank
            best_ann = ann
            best_ann["_NORM_GENE"] = norm_gene
            best_ann["_RAW_GENE"] = gene_raw
    if not best_ann: return None
    norm_gene = best_ann.get("_NORM_GENE")
    gene_raw = best_ann.get("_RAW_GENE")
    # ENSURE we only build records for ACMG SF genes (acmg_genes_normalized passed from caller)
    if acmg_genes_normalized and norm_gene not in acmg_genes_normalized:
        return None
    transcript_id = best_ann.get("Feature", "")
    prob_gt = "./."; prob_dp = None; prob_ad = None; sample_name = ""
    prob_af = 0.0
    idx = 0
    if hasattr(rec, "samples"):
        s_list = list(rec.samples)
        if proband_sample and proband_sample in s_list:
            idx = s_list.index(proband_sample)
            sample_name = proband_sample
        else:
            sample_name = s_list[0] if s_list else "Proband"
    try:
        if hasattr(rec, "genotypes"):
            gt_data = rec.genotypes[idx]
            if len(gt_data) >= 2:
                s1 = str(gt_data[0]) if gt_data[0] != -1 else "."; s2 = str(gt_data[1]) if gt_data[1] != -1 else "."
                prob_gt = f"{s1}/{s2}"
            dps = rec.format('DP')
            if dps is not None:
                val = dps[idx]
                try: prob_dp = int(val.item() if hasattr(val, "item") else val)
                except:
                    try: prob_dp = int(val[0])
                    except: prob_dp = None
            afs = rec.format('AF')
            if afs is not None:
                val = afs[idx]
                if hasattr(val, "__len__") and len(val)>0: prob_af = float(val[0])
                else: prob_af = float(val)
        else:
            s_dat = getattr(rec, "samples", {})[sample_name]
            g = s_dat.get("GT"); prob_gt = "/".join("." if x is None else str(x) for x in g) if g else "./."
            prob_dp = s_dat.get("DP")
            if "AF" in s_dat:
                val = s_dat.get("AF")
                if isinstance(val, (list, tuple)): prob_af = float(val[0])
                else: prob_af = float(val)
    except:
        pass
    if prob_af == 0.0:
        try:
            raw_line = str(rec).strip()
            parts = raw_line.split('\t')
            if len(parts) > 9:
                fmt_str = parts[8]
                samp_str = parts[9 + idx]
                fmt_keys = fmt_str.split(':')
                samp_vals = samp_str.split(':')
                if len(fmt_keys) == len(samp_vals):
                    if "AF" in fmt_keys:
                        i_af = fmt_keys.index("AF")
                        raw_val = samp_vals[i_af]
                        if "," in raw_val: raw_val = raw_val.split(",")[0]
                        if raw_val and raw_val != ".":
                            prob_af = float(raw_val)
                    if not prob_dp and "DP" in fmt_keys:
                        i_dp = fmt_keys.index("DP")
                        raw_val = samp_vals[i_dp]
                        if raw_val and raw_val != ".":
                            prob_dp = int(raw_val)
                    if prob_gt == "./." and "GT" in fmt_keys:
                        i_gt = fmt_keys.index("GT")
                        prob_gt = samp_vals[i_gt]
                    if prob_af == 0.0 and "AD" in fmt_keys:
                         i_ad = fmt_keys.index("AD")
                         raw_ad = samp_vals[i_ad]
                         if "," in raw_ad:
                             ads = [int(x) if x.isdigit() else 0 for x in raw_ad.split(",")]
                             if len(ads) >= 2 and sum(ads) > 0:
                                 prob_af = ads[1] / sum(ads)
                                 if not prob_dp: prob_dp = sum(ads)
        except Exception:
            pass
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
                        except Exception: pass
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
    vep_nmdesc = ""
    if "NMD" in best_ann and best_ann["NMD"]: vep_nmdesc = best_ann["NMD"]
    if not vep_nmdesc: vep_nmdesc = get_val_anywhere(["NMD", "NMDEscPredictor", "nmdesc_predictor"], str)
    splice_keys = ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_DS_max", "SpliceAI"]
    splice_vals = []
    for key in splice_keys:
        val = None
        if key in best_ann and best_ann[key]:
            try:
                val = float(best_ann[key])
            except:
                try:
                    if isinstance(best_ann[key], (list,tuple)):
                        val = max([float(x) for x in best_ann[key] if x not in ("", ".")])
                except:
                    val = None
        if val is None and key.lower() in info_l:
            try:
                val = float(info_l[key.lower()])
            except:
                val = None
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
    try:
        h = clean_hgvsp(vr.hgvsp or "")
        m = re.search(r'p\.([A-Za-z]{1,3}|[A-Za-z])(\d+)([A-Za-z\*]{1,3})', h)
        if m:
            vr.aa_ref = m.group(1); vr.aa_pos = int(m.group(2)); vr.aa_alt = m.group(3)
    except:
        pass
    try:
        if dbm:
            dbm.get_variant_db_evidence(vr)
    except Exception:
        pass
    return vr

# ---------------------------
# Pre-analysis - QC sex inference
# ---------------------------
def perform_pre_analysis_qc(vcf_path: str, proband_id: str) -> Dict[str, Any]:
    warnings = []
    inferred_sex = "Unknown"
    if HAVE_PYSAM:
        try:
            vf = pysam.VariantFile(vcf_path)
            chroms = list(vf.header.contigs)
            chr_x_name = next((c for c in chroms if c.lower() in ('chrx', 'x', '23')), None)
            chr_y_name = next((c for c in chroms if c.lower() in ('chry', 'y', '24')), None)
            has_y_variants = False
            if chr_y_name:
                try:
                    for rec in vf.fetch(chr_y_name):
                        has_y_variants = True
                        break
                except: pass
            if has_y_variants:
                inferred_sex = "Male"
                logger.info(f"QC Sex Check: Inferred Male (Variants found on {chr_y_name})")
                return {"warnings": warnings, "sex": inferred_sex}
            if chr_x_name:
                het_count = 0
                hom_count = 0
                total_x = 0
                try:
                    for rec in vf.fetch(chr_x_name):
                        if len(rec.samples) == 0: break
                        s_name = proband_id if proband_id in rec.samples else rec.samples[0].name
                        gt = rec.samples[s_name]['GT']
                        if len(gt) >= 2:
                            if gt[0] != gt[1]: het_count += 1
                            else: hom_count += 1
                            total_x += 1
                        if total_x > 2000: break
                except Exception: pass
                if total_x > 0:
                    ratio = het_count / total_x
                    logger.info(f"QC Sex Check: ChrX Stats - Total: {total_x}, Het: {het_count}, Ratio: {ratio:.2f}")
                    if ratio > 0.15:
                        inferred_sex = "Female"
                    elif ratio < 0.10:
                        inferred_sex = "Male"
                    else:
                        inferred_sex = "Unknown"
                else:
                     logger.warning("QC Sex Check: No variants found on ChrX to infer sex.")
        except Exception as e:
            logger.warning(f"QC Sex Check failed: {e}")
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
    if (auto_class.startswith("Path") or auto_class.startswith("Likely path")) and path_cnt < MIN_PATHOGENIC_CRITERIA:
        auto_class = "VUS"
    return auto_class, total, path_cnt

# ---------------------------
# Process VCF (main pipeline)
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
                ttn_meta: Optional[str] = None,
                require_mane: bool = False,
                cooccurrence_path: Optional[str] = None, 
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
    if run_vep_dbnsfp:
        if not dbnsfp_path or not os.path.exists(dbnsfp_path):
            raise FileNotFoundError(f"dbNSFP file not found: {dbnsfp_path}")
        out_v = os.path.join(outdir, os.path.basename(working_vcf).replace(".vcf","").replace(".gz","") + ".dbnsfp.vep.vcf.gz")
        dbnsfp_plugin_string = choose_dbnsfp_fields(dbnsfp_path, dbnsfp_fields or [])
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
    fasta_handle = None
    if fasta and os.path.exists(fasta) and HAVE_PYSAM:
        try: fasta_handle = pysam.FastaFile(fasta)
        except: pass
    gene_rules = load_acmg_table(acmg_table) if HAVE_PANDAS else {}
    acmg_genes_normalized = set(g for g in gene_rules.keys() if gene_rules.get(g,{}).get("reportable", False))
    dbm = DatabaseManager(db_paths, gnomad_v2=gnomad_v2, gnomad_v3=gnomad_v3, gnomad_v4=gnomad_v4, internal_db_path=internal_db_path)
    exons_df = load_exons_table(exons_file) if exons_file else None
    # load optional transcript mapping and NMD precomputed table (if provided)
    tx_map = {}
    if transcript_map:
        try:
            from typing import Any
            tx_map = load_transcript_map(transcript_map)
        except Exception as e:
            logger.warning("Failed to load transcript map %s: %s", transcript_map, e)
            tx_map = {}

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
    csq_header = parse_csq_header_from_vcf(annotated_vcf)
    vcf_reader = None
    if HAVE_CYVCF2: vcf_reader = VCF(annotated_vcf)
    else: vcf_reader = pysam.VariantFile(annotated_vcf)
    dad_vf = pysam.VariantFile(father_vcf) if father_vcf and os.path.exists(father_vcf) and HAVE_PYSAM else None
    mom_vf = pysam.VariantFile(mother_vcf) if mother_vcf and os.path.exists(mother_vcf) and HAVE_PYSAM else None
    candidates: List[VariantRecord] = []
    logger.info("Parsing VCF...")
    if HAVE_CYVCF2 and isinstance(vcf_reader, VCF):
        for rec in vcf_reader:
            vr = build_variant_record_from_rec(rec, csq_header, acmg_genes_normalized, dad_vf, mom_vf, exons_df, proband_sample=proband_sample, require_mane=require_mane, dbm=dbm)
            if not vr: continue
            vr.sample_sex = inferred_sex
            if (vr.gnomad_af is None or vr.gnomad_af == 0.0): dbm.annotate_gnomad_for_variant(vr)
            try: vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
            except: pass
            candidates.append(vr)
    else:
        for rec in vcf_reader.fetch():
            vr = build_variant_record_from_rec(rec, csq_header, acmg_genes_normalized, dad_vf, mom_vf, exons_df, proband_sample=proband_sample, require_mane=require_mane, dbm=dbm)
            if not vr: continue
            vr.sample_sex = inferred_sex
            if (vr.gnomad_af is None or vr.gnomad_af == 0.0): dbm.annotate_gnomad_for_variant(vr)
            try: vr.db_hits.extend(dbm.get_variant_db_evidence(vr))
            except: pass
            candidates.append(vr)
    candidates = deduplicate_candidates(candidates)
    if proband_sample:
        for v in candidates: v.sample = proband_sample
    engine = ACMGEngine(dbm, gene_rules, exons_df, {}, THRESH["SPLICEAI_MODERATE"], require_mane=require_mane)
    # Attach transcript/exon/NMD resources to engine for PVS1/NMD decision
    engine.tx_map = tx_map
    engine.exons_map = exons_df
    engine.nmd_map = nmd_map
    # pass gnomAD API options into engine/global space for on-the-fly cooccurrence queries if needed
    engine.gnomad_api_endpoint = gnomad_api_endpoint or "https://gnomad.broadinstitute.org/api"
    engine.gnomad_dataset = gnomad_dataset or "gnomad_r2_1"
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
        if "BA1" not in v.criteria_assigned:
            if (v.automated_class.startswith("Path") or v.automated_class.startswith("Likely path")) and v.pathogenic_criteria_count < MIN_PATHOGENIC_CRITERIA:
                v.automated_class = "VUS"
                v.manual_reasons.append(f"Insufficient pathogenic criteria (<{MIN_PATHOGENIC_CRITERIA})")
    if fasta_handle: fasta_handle.close()
    engine.apply_pm3_batch(candidates)
    detect_delins(candidates)
    for v in candidates:
        v.total_points = sum(int(x) for x in v.criteria_points.values()) if v.criteria_points else 0
        auto_class, total, path_cnt = compute_class_from_points(v.criteria_points)
        v.automated_class = auto_class
        v.total_points = total
        v.pathogenic_criteria_count = path_cnt
    filtered_for_report = []
    dropped_count = 0
    for v in candidates:
        if not STRICT_QC_ENABLED or not strict_qc:
            filtered_for_report.append(v); continue
        dp = v.proband_dp or 0
        af = v.proband_af or 0.0
        if dp >= STRICT_MIN_DP and af >= STRICT_MIN_AF:
            filtered_for_report.append(v)
        else:
            v.filtered_reasons.append("Failed strict QC (DP/AF)")
            dropped_count += 1
    logger.info("After strict QC filtering: %d kept for Report, %d low quality (audit only)", len(filtered_for_report), dropped_count)
    all_path, auto_path, manual_path, cnt_auto, cnt_manual = strict_triage_and_output(
        report_candidates=filtered_for_report,
        audit_candidates=candidates,
        gene_rules=gene_rules,
        recs_map=recs_map,
        outdir=outdir,
        family_id_input=family_id,
        dbm=dbm,                        
        cooccur_path=cooccurrence_path
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
        logger.info("--- DEBUG: ACMG Table Columns Detected ---")
        logger.info(list(df.columns))
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
    
    if processed_genes:
        logger.info(f"Обработано {len(processed_genes)} уникальных генов")
        logger.info(f"Исходных строк в таблице: {len(df)}")
        logger.info(f"Объединено {len(df) - len(processed_genes)} дублирующих записей")
        
        genes_with_multiple = [g for g in processed_genes if ";" in rules.get(g, {}).get("disease", "")]
        if genes_with_multiple:
            logger.info(f"Гены с объединенными заболеваниями ({len(genes_with_multiple)}):")
            for gene in genes_with_multiple[:5]:  
                disease = rules[gene].get("disease", "")
                logger.info(f"  {gene}: {disease}")
    
    if (not any_reportable) and len(rules) > 0:
        logger.warning("No explicit reportable flags detected. Marking ALL genes in table as reportable.")
        for g in rules: rules[g]["reportable"] = True
    logger.info("Loaded ACMG rules for %d genes from %s", len(rules), path)
    return rules

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
    p.add_argument("--require-mane", action="store_true", help="Require MANE or MANE Plus Clinical annotation")
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
    p.add_argument("--hgmd", help="Path to HGMD Pro CSV (hg38)")
    p.add_argument("--gnomad-cooccurrence", dest="gnomad_cooccurrence", help="Path to TSV with precomputed variant cooccurrence/cis-trans evidence (optional). Format: var1_vid\\tvar2_vid\\tcis_flag\\tevidence")
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
        map_id = os.path.basename(dirname)
        s_out = os.path.join(args_dict.outdir, "individual_results", dna_id)
        siblings = glob.glob(os.path.join(dirname, "*"))
        f_vcf = next((s for s in siblings if "_father" in s and (s.endswith(".vcf") or s.endswith(".vcf.gz"))), None)
        m_vcf = next((s for s in siblings if "_mother" in s and (s.endswith(".vcf") or s.endswith(".vcf.gz"))), None)
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
            family_id=map_id,
            ttn_meta=args_dict.ttn_meta,
            require_mane=args_dict.require_mane,
            cooccurrence_path=getattr(args_dict, "gnomad_cooccurrence", None),
            gnomad_api_endpoint=getattr(args_dict, "gnomad_api_endpoint", None),
            gnomad_dataset=getattr(args_dict, "gnomad_dataset", None),
            transcript_map=getattr(args_dict, "transcript_map", None),
            nmd_table=getattr(args_dict, "nmd_table", None) 
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
    patterns = [
        os.path.join(args.batch_input_dir, "**", "*_proband*.vcf*"),
        os.path.join(args.batch_input_dir, "**", "*_sample*.vcf*"),
        os.path.join(args.batch_input_dir, "**", "*_proband*.vcf"),
        os.path.join(args.batch_input_dir, "**", "*_sample*.vcf"),
    ]
    # aggregate unique
    proband_files = []
    for p in patterns:
        proband_files.extend(glob.glob(p, recursive=True))
    proband_files = sorted(list(set(proband_files)))
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
            require_mane=args.require_mane,
            gene_bed=args.gene_bed,
            proband_sample=args.proband_sample,
            strict_qc=not args.no_strict_qc,
            gnomad_v2=args.gnomad_v2,
            gnomad_v3=args.gnomad_v3,
            acmg_recs=args.acmg_recs,
            internal_db_path=args.internal_db,
            gnomad_v4=args.gnomad_v4,
            cooccurrence_path=getattr(args, "gnomad_cooccurrence", None)
        )
    else:
        parser.print_help()

if __name__ == "__main__":
    main()