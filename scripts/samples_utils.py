#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
samples_utils.py

Tool to match an incoming sample/variant list (TSV/CSV) to a samples.json inventory,
augment rows with sample metadata (matched sample name, db path, meta path, sex),
and optionally generate extraction commands (extract_all_variants ...) grouped by
'Карта' (family/map) so that members of the same Karte go into one folder (trio handling).

Enhancements borrowed from generate.py:
 - normalize roles (proband/mother/father/sibling)
 - safe folder naming based on Karte
 - attempt to locate on-disk variants.sqlite using STATUS/SCI CSVs and a BASE_DIR
 - create run_jobs.sh with extract_all_variants calls (one per sample/member)
 - optional --write-vcfs: directly write VCFs from variants.sqlite (Python implementation)
 - optional --passport-csv: passport sex mapping and sex mismatch column
"""
from __future__ import annotations

import os
import sys
import json
import csv
import argparse
import logging
import re
import sqlite3
from typing import List, Dict, Any, Optional, Tuple, Set

# optional pandas / pyyaml
try:
    import pandas as pd  # type: ignore
    HAVE_PANDAS = True
except Exception:
    pd = None
    HAVE_PANDAS = False

try:
    import yaml  # type: ignore
    HAVE_YAML = True
except Exception:
    yaml = None
    HAVE_YAML = False

logger = logging.getLogger("samples_utils")
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Default header file used when writing VCFs
HEADER_VCF = "header.vcf"


# ----------------------
# Utility helpers (from generate.py)
# ----------------------
def clean_id(val: Any) -> str:
    """Remove leading zeros and special chars from DNA ID (best-effort)."""
    if val is None:
        return ""
    s = str(val).strip().split('.')[0]
    # if purely digits -> remove leading zeros
    if s.isdigit():
        return s.lstrip('0') or "0"
    # otherwise return trimmed
    return s


def safe_filename(name: str) -> str:
    """Make safe directory name (replace slashes/backslashes and strip)."""
    if name is None:
        return ""
    return str(name).replace('/', '_').replace('\\', '_').strip()


def normalize_role(role_raw: Any) -> str:
    """Map role text (possibly Russian) to one of: proband, mother, father, sibling, sample."""
    if role_raw is None:
        return "sample"
    r = str(role_raw).lower()
    if any(x in r for x in ("proband", "prob", "пробанд", "обратившийся", "обратившийся")):
        return "proband"
    if any(x in r for x in ("mother", "mom", "мать", "мама")):
        return "mother"
    if any(x in r for x in ("father", "dad", "отец", "папа")):
        return "father"
    if any(x in r for x in ("sibling", "sib", "сиб", "брат", "сестра")):
        return "sibling"
    return "sample"


def get_folder_type(study_name: Optional[str]) -> str:
    s = str(study_name).lower() if study_name is not None else ""
    return "panels" if "панел" in s or "панель" in s else "smk"


def detect_delimiter(sample_path: str, bytes_to_read: int = 4096) -> str:
    try:
        with open(sample_path, "rt", encoding="utf-8", errors="replace") as fh:
            sample = fh.read(bytes_to_read)
            sn = csv.Sniffer()
            dialect = sn.sniff(sample)
            return dialect.delimiter
    except Exception:
        # fallback heuristics
        try:
            with open(sample_path, "rt", encoding="utf-8", errors="replace") as fh:
                head = fh.readline()
                if "\t" in head:
                    return "\t"
                if ";" in head:
                    return ";"
                if "," in head:
                    return ","
        except Exception:
            pass
    return "\t"


# ----------------------
# File loaders / parsers
# ----------------------
def load_samples(samples_json: str) -> List[Dict[str, Any]]:
    if not os.path.exists(samples_json):
        raise FileNotFoundError(f"samples.json not found: {samples_json}")
    with open(samples_json, "rt", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        raise ValueError("samples.json must contain a top-level array of sample objects")
    return data


def load_filter_list(path: str, id_cols: Optional[List[str]] = None) -> Tuple[List[Dict[str, str]], List[str]]:
    """
    Load a filter CSV/TSV and return rows (list of dicts) and headers list.
    id_cols not used here except in main to select ids.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    delim = detect_delimiter(path)
    rows = []
    headers: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=delim)
        try:
            first = next(reader)
        except StopIteration:
            return [], []
        # detect header-like row
        header_like = any(re.search(r'[A-Za-zА-Яа-яЁё]', str(x or "")) for x in first)
        if header_like:
            headers = [h.strip() for h in first]
        else:
            # no header -> create numeric headers
            headers = [str(i) for i in range(len(first))]
            rows.append({headers[i]: first[i].strip() for i in range(len(first))})
        for r in reader:
            if not any(cell.strip() for cell in r):
                continue
            if len(r) < len(headers):
                r = r + [""] * (len(headers) - len(r))
            row = {headers[i]: r[i].strip() for i in range(len(headers))}
            rows.append(row)
    return rows, headers


# ----------------------
# Indexing / matching helpers
# ----------------------
def build_name_index(samples: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    idx: Dict[str, Dict[str, Any]] = {}
    for s in samples:
        name = str(s.get("name", "")).strip()
        if not name:
            continue
        idx[name] = s
    return idx


def build_value_indexes(samples: List[Dict[str, Any]]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Returns:
      - value_to_names_cs: exact value -> [sample_name,...]
      - value_to_names_ci: lowercase value -> [sample_name,...]
    Values are taken from all string scalar fields in each sample record.
    """
    name_index = build_name_index(samples)
    v_cs: Dict[str, List[str]] = {}
    v_ci: Dict[str, List[str]] = {}
    for s in samples:
        sname = str(s.get("name", "")).strip()
        for k, v in s.items():
            if v is None:
                continue
            if isinstance(v, (list, dict)):
                continue
            vv = str(v).strip()
            if vv == "":
                continue
            v_cs.setdefault(vv, []).append(sname)
            v_ci.setdefault(vv.lower(), []).append(sname)
    return v_cs, v_ci


def read_sex_from_meta(sample: Dict[str, Any], meta_subdir: Optional[str] = None) -> Optional[str]:
    """
    Try to read sex from meta YAML under sample['work_dir']/meta_subdir/<sample['name']>.yaml.
    Uses PyYAML if available; otherwise does simple regex search for 'sex:' line.
    """
    genome = sample.get("genome", "hg38")
    meta_subdir = meta_subdir or f"meta_{genome}"
    wd = sample.get("work_dir", "")
    name = sample.get("name", "")
    meta_fn = os.path.join(wd, meta_subdir, f"{name}.yaml")
    if not os.path.exists(meta_fn):
        return None
    # PyYAML preferred
    if HAVE_YAML:
        try:
            with open(meta_fn, "rt", encoding="utf-8") as fh:
                data = yaml.safe_load(fh)
            if isinstance(data, dict):
                for k, v in data.items():
                    if str(k).lower() == "sex":
                        return str(v).upper() if v is not None else None
                # nested search
                for val in data.values():
                    if isinstance(val, dict):
                        for k2, v2 in val.items():
                            if str(k2).lower() == "sex":
                                return str(v2).upper() if v2 is not None else None
        except Exception:
            pass
    # fallback simple text search
    try:
        with open(meta_fn, "rt", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                m = re.match(r'^\s*sex\s*:\s*["\']?(?P<v>[^"\']+)', ln, flags=re.I)
                if m:
                    return m.group("v").strip().upper()
    except Exception:
        pass
    return None


# ----------------------
# Disk search helpers (adapted/simplified from generate.py)
# ----------------------
def validate_and_fix_path(base_dir: str, dir_name: str, sub_type: str, dna_padded: str) -> Optional[str]:
    """
    Try to build plausible path to variants.sqlite for given dir_name and dna_padded.
    Will try common genome directories (hg38, hg19).
    """
    base_dir = os.path.expanduser(base_dir)
    alt = "panels" if sub_type == "smk" else "smk"
    # Try typical genomes
    for genome in ("hg38", "hg19"):
        path_primary = os.path.join(base_dir, dir_name, sub_type, dna_padded, genome, "db", "variants.sqlite")
        if os.path.isfile(path_primary):
            return path_primary
        path_alt = os.path.join(base_dir, dir_name, alt, dna_padded, genome, "db", "variants.sqlite")
        if os.path.isfile(path_alt):
            return path_alt
    return None


def load_csv_table(path: str) -> Optional[Any]:
    """Load status/sci CSV using pandas if available, else return None for simple CSV parse."""
    if path is None:
        return None
    if not os.path.exists(path):
        logger.warning("CSV path not found: %s", path)
        return None
    if HAVE_PANDAS:
        try:
            df = pd.read_csv(path, dtype=str, sep=None, engine="python").fillna("")
            return df
        except Exception as e:
            logger.warning("Pandas failed to read %s: %s", path, e)
            return None
    else:
        # fallback simple reader -> list of dicts
        try:
            delim = detect_delimiter(path)
            with open(path, "rt", encoding="utf-8", errors="replace") as fh:
                reader = csv.DictReader(fh, delimiter=delim)
                rows = [r for r in reader]
            return rows
        except Exception:
            return None


def find_col(df_or_rows, candidates: List[str]) -> Optional[str]:
    """
    Find a column name in DataFrame or list-of-dicts by case-insensitive candidate lookup.
    Returns actual column name (for DataFrame) or key (for dict rows) or None.
    """
    if df_or_rows is None:
        return None
    # pandas DataFrame
    if HAVE_PANDAS and isinstance(df_or_rows, pd.DataFrame):
        cols_lower = {c.lower(): c for c in df_or_rows.columns}
        for cand in candidates:
            if cand.lower() in cols_lower:
                return cols_lower[cand.lower()]
        # fuzzy: return a column that contains substring
        for cand in candidates:
            for real in cols_lower:
                if cand.lower() in real:
                    return cols_lower[real]
        return None
    else:
        # assume list of dicts -> check keys in first row
        try:
            first = df_or_rows[0] if df_or_rows else {}
            keys_lower = {k.lower(): k for k in first.keys()}
            for cand in candidates:
                if cand.lower() in keys_lower:
                    return keys_lower[cand.lower()]
            for cand in candidates:
                for real in keys_lower:
                    if cand.lower() in real:
                        return keys_lower[real]
            return None
        except Exception:
            return None


def search_status_and_sci_for_group(
    group_rows: List[Dict[str, str]],
    status_table: Any,
    sci_table: Any,
    base_dir: Optional[str],
    found_ids_tracker: Set[str],
) -> List[Tuple[str, str, str]]:
    """
    For a group of rows (same 'Карта') try to locate disk DBs and produce extract commands.
    Returns list of tuples: (db_path, out_param, desc)
    Uses STATUS table first, then SCI if not found.
    This is simplified/adapted from generate.py.
    """
    if status_table is None and sci_table is None:
        return []
    commands = []
    # Attempt to extract DNK ids from group
    dnk_list = [clean_id(r.get("ДНК") or r.get("ДНК ", "") or "") for r in group_rows if r.get("ДНК") or r.get("ДНК ")]
    dnk_list = [d for d in dnk_list if d]
    # We'll try status first
    def try_table(table, is_status: bool):
        if table is None:
            return []
        found = []
        id_col = find_col(table, ['ID', 'DNA', 'Id', 'Id_sample', 'sample'])
        dir_col = find_col(table, ['dir_name', 'directory', 'path', 'dir'])
        fam_col = find_col(table, ['fam', 'map', 'карта', 'Карта'])
        rel_col = find_col(table, ['rel', 'role', 'отношение'])
        if not id_col or not dir_col:
            return []
        # Pandas
        rows = []
        if HAVE_PANDAS and isinstance(table, pd.DataFrame):
            rows = table.to_dict(orient="records")
        else:
            rows = table
        for r in rows:
            try:
                r_id = clean_id(r.get(id_col, "") if isinstance(r, dict) else "")
            except Exception:
                r_id = ""
            if not r_id:
                continue
            if r_id in dnk_list:
                d_name = str(r.get(dir_col, "")).strip()
                if not d_name or d_name.lower() in ("не получен", "nan", ""):
                    continue
                # pad dna to 12 digits to try typical path formats
                dna_pad = r_id.zfill(12)
                # choose a folder type from incoming group's first study
                stype = get_folder_type(group_rows[0].get("НазваниеИсследования", ""))
                db_path = None
                if base_dir:
                    db_path = validate_and_fix_path(base_dir, d_name, stype, dna_pad)
                if db_path:
                    # decide role
                    role = normalize_role(next((gr.get("Отношение") for gr in group_rows if clean_id(gr.get("ДНК") or "") == r_id), "sample"))
                    # folder_name based on Karte (we will handle outer folder) -> but here return out_param suffix
                    # out_param will be filled by caller using Karte folder
                    desc = ("STATUS" if is_status else "SCI") + " | " + (str(r.get(fam_col, "") or "") if fam_col else "")
                    found.append((r_id, db_path, role, desc))
        return found

    res = try_table(status_table, True)
    if not res:
        res = try_table(sci_table, False)
    # convert to extract command tuples
    cmds = []
    for r_id, dbp, role, desc in res:
        out_param = f"{r_id}_{role}"
        cmds.append((dbp, out_param, desc))
        found_ids_tracker.add(r_id)
    return cmds


# ----------------------
# VCF writer (Python implementation to match generate.py output)
# ----------------------
def _normalize_colname_lookup(cols: List[Tuple[int, str]]) -> Dict[str, str]:
    """Helper to build case-insensitive lookup dict from PRAGMA table_info results."""
    mapping = {}
    for _, name in cols:
        mapping[name.lower()] = name
    return mapping


def make_vcfs_from_sqlite(db_path: str, out_vcf_path: str, header_vcf_path: str = HEADER_VCF) -> Tuple[bool, str]:
    """
    Read table vcf from SQLite DB and write a VCF file similar to generate.py output.
    Best-effort column detection. Returns (success, message).
    """
    if not os.path.exists(db_path):
        return False, "DB_NOT_FOUND"
    try:
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("PRAGMA table_info(vcf);")
        cols = cur.fetchall()
        if not cols:
            # fallback: try to find any table that looks like vcf
            cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [r[0] for r in cur.fetchall()]
            possible = next((t for t in tables if 'vcf' in t.lower()), None)
            if possible:
                cur.execute(f"PRAGMA table_info({possible});")
                cols = cur.fetchall()
            if not cols:
                conn.close()
                return False, "NO_VCF_TABLE"
        colnames = [c[1] for c in cols]
        colmap = _normalize_colname_lookup([(c[0], c[1]) for c in cols])

        # required coordinate cols
        required = ['chrom', 'pos', 'ref', 'alt']
        lower_cols = set([c.lower() for c in colnames])
        if not all(rc in lower_cols for rc in required):
            conn.close()
            return False, "MISSING_COORD_COLUMNS"

        # helper for choosing column by candidate substrings
        def pick(cands: List[str]) -> Optional[str]:
            for cand in cands:
                if cand.lower() in colmap:
                    return colmap[cand.lower()]
            # substring match
            for cand in cands:
                for k in colmap:
                    if cand.lower() in k:
                        return colmap[k]
            return None

        gene_col = pick(['gene', 'gene_symbol', 'symbol', 'geneinfo'])
        exon_col = pick(['exon', 'exon_number'])
        hgvsc_col = pick(['hgvsc', 'c_dna', 'cdna'])
        hgvsp_col = pick(['hgvsp', 'protein', 'aa_change'])
        impact_col = pick(['impact'])
        ann_col = pick(['consequence', 'ann', 'csq', 'annotation'])
        revel_col = pick(['revel'])
        alpha_col = pick(['alpha_missense', 'alphamissense'])
        splice_col = pick(['splice_ai', 'splice_ai_score', 'spliceai'])

        gt_col = pick(['gt', 'genotype'])
        dp_col = pick(['dp', 'depth'])
        gq_col = pick(['gq', 'genotype_quality'])
        ad_col = pick(['ad', 'allelic_depth'])
        af_col = pick(['af', 'allele_freq', 'max_af'])

        sample_expr_parts = []
        sample_expr_parts.append(f"COALESCE({gt_col}, './.')") if gt_col else sample_expr_parts.append("'./.'")
        sample_expr_parts.append(f"COALESCE(CAST({dp_col} AS TEXT), '.')") if dp_col else sample_expr_parts.append("'.'")
        sample_expr_parts.append(f"COALESCE(CAST({gq_col} AS TEXT), '.')") if gq_col else sample_expr_parts.append("'.'")
        sample_expr_parts.append(f"COALESCE(CAST({ad_col} AS TEXT), '.,.')") if ad_col else sample_expr_parts.append("'.,.'")
        sample_expr_parts.append(f"COALESCE(CAST({af_col} AS TEXT), '0.0')") if af_col else sample_expr_parts.append("'0.0'")
        sample_expr = " || ':' || ".join(sample_expr_parts)

        info_parts = []
        if gene_col:
            info_parts.append(f"'GENE=' || COALESCE({gene_col}, '')")
        else:
            info_parts.append("'GENE=.'")
        if exon_col:
            info_parts.append(f"'Exon=' || COALESCE({exon_col}, '')")
        if hgvsc_col:
            info_parts.append(f"'HGVSc=' || COALESCE({hgvsc_col}, '')")
        if hgvsp_col:
            info_parts.append(f"'HGVSp=' || COALESCE({hgvsp_col}, '')")
        if impact_col:
            info_parts.append(f"'Impact=' || COALESCE({impact_col}, '')")
        if ann_col:
            info_parts.append(f"'Consequence=' || COALESCE({ann_col}, '')")
        if revel_col:
            info_parts.append(f"'REVEL=' || COALESCE(CAST({revel_col} AS TEXT),'')")
        if alpha_col:
            info_parts.append(f"'AlphaMissense=' || COALESCE(CAST({alpha_col} AS TEXT),'')")
        if splice_col:
            info_parts.append(f"'splice_ai_score=' || COALESCE(CAST({splice_col} AS TEXT),'')")

        info_expr = " || ';' || ".join(info_parts)

        sql = f"SELECT CHROM, POS, vid, REF, ALT, COALESCE(QUAL,50) as QUAL, COALESCE(FILTER,'PASS') as FILTER, ({info_expr}) as INFO, 'GT:DP:GQ:AD:AF' as FORMAT, ({sample_expr}) as SAMPLE FROM vcf ORDER BY rowid ASC;"
        try:
            cur.execute(sql)
            rows = cur.fetchall()
        except Exception as e:
            conn.close()
            return False, f"SQL_FAIL:{e}"

        os.makedirs(os.path.dirname(out_vcf_path) or ".", exist_ok=True)
        with open(out_vcf_path, "wt", encoding="utf-8") as out:
            if os.path.exists(header_vcf_path):
                with open(header_vcf_path, "rt", encoding="utf-8") as hh:
                    out.write(hh.read())
            else:
                # minimal header
                out.write('##fileformat=VCFv4.2\n')
                out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\n')
                out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
                out.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">\n')
                out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">\n')
                out.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Fraction">\n')
                out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n")
            for r in rows:
                chrom, pos, vid, ref, alt, qual, filt, info, fmt, sample = r
                out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{fmt}\t{sample}\n")
        conn.close()
        return True, f"WROTE:{out_vcf_path}"
    except Exception as e:
        return False, f"ERROR:{e}"


# ----------------------
# Passport CSV loader & sex normalization
# ----------------------
def load_passport_sex_from_csv(csv_path: str, sample_name_col_candidates=("sample", "name", "Sample", "sample_name", "dir_name")) -> Dict[str, str]:
    """
    Return dict sample_name -> passport_sex (raw string).
    Uses pandas if available for convenience.
    """
    out: Dict[str, str] = {}
    if not csv_path or not os.path.exists(csv_path):
        return out
    if HAVE_PANDAS:
        try:
            df = pd.read_csv(csv_path, dtype=str).fillna("")
            cols = [c for c in df.columns]
            # find sample & sex columns
            sample_col = None
            for cand in sample_name_col_candidates:
                for c in cols:
                    if cand.lower() == c.lower():
                        sample_col = c
                        break
                if sample_col:
                    break
            sex_col = None
            for c in cols:
                if "sex" in c.lower() or "gender" in c.lower():
                    sex_col = c
                    break
            if sample_col and sex_col:
                for _, r in df.iterrows():
                    name = str(r[sample_col]).strip()
                    sex = str(r[sex_col]).strip()
                    if name:
                        out[name] = sex
        except Exception:
            # fallback to simple reader
            pass
    # simple reader fallback
    if not out:
        try:
            delim = detect_delimiter(csv_path)
            with open(csv_path, "rt", encoding="utf-8", errors="replace") as fh:
                dr = csv.DictReader(fh, delimiter=delim)
                header = dr.fieldnames or []
                sample_col = None
                sex_col = None
                for cand in sample_name_col_candidates:
                    for h in header:
                        if cand.lower() == h.lower():
                            sample_col = h; break
                    if sample_col:
                        break
                for h in header:
                    if "sex" in h.lower() or "gender" in h.lower():
                        sex_col = h; break
                if sample_col and sex_col:
                    for r in dr:
                        name = r.get(sample_col, "").strip()
                        sex = r.get(sex_col, "").strip()
                        if name:
                            out[name] = sex
        except Exception:
            pass
    return out


def normalize_sex_token(s: Optional[str]) -> Optional[str]:
    """Normalize sex tokens into canonical tokens for comparison (XX / XY / XXY / M / F)."""
    if s is None:
        return None
    st = str(s).strip()
    if not st:
        return None
    up = st.upper()
    # Common canonicalizations
    if up in ("XX", "F", "FEMALE", "Ж", "ЖЕН", "ЖЕНСКИЙ", "ЖЕНЩИНА"):
        return "XX"
    if up in ("XY", "M", "MALE", "М", "МУЖ", "МУЖСКОЙ", "МУЖЧИНА"):
        return "XY"
    if "XX" in up and "Y" in up:
        return "XXY"
    if up in ("UNKNOWN", "NA", "N/A", "."):
        return None
    # fallback: return uppercase token trimmed
    return up


# ----------------------
# Core matching + output generation
# ----------------------
def match_ids_to_samples(samples: List[Dict[str, Any]], ids: List[str], meta_subdir: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
    """
    Given list of identifier strings, attempt match against samples.json.
    Returns mapping ident -> {found, matched_sample_name, sample(dict), db_path, meta_path, sex}
    """
    idx = build_name_index(samples)
    v_cs, v_ci = build_value_indexes(samples)
    out: Dict[str, Dict[str, Any]] = {}
    for ident in ids:
        ident_s = str(ident).strip()
        rec = {"found": False, "matched_sample_name": None, "sample": None, "db_path": None, "meta_path": None, "sex": None}
        if ident_s in idx:
            s = idx[ident_s]
            rec["found"] = True
            rec["matched_sample_name"] = s.get("name")
            rec["sample"] = s
            rec["db_path"] = os.path.join(s.get("work_dir", ""), s.get("name", ""), s.get("genome", "hg38"), "db", "variants.sqlite")
            rec["meta_path"] = None
            if s.get("work_dir"):
                meta_dir = meta_subdir or f"meta_{s.get('genome','hg38')}"
                meta_fn = os.path.join(s.get("work_dir", ""), meta_dir, f"{s.get('name')}.yaml")
                if os.path.exists(meta_fn):
                    rec["meta_path"] = meta_fn
            rec["sex"] = read_sex_from_meta(s, meta_subdir=meta_subdir)
            out[ident_s] = rec
            continue
        if ident_s in v_cs:
            nm = v_cs[ident_s][0]
            s = idx.get(nm)
            if s:
                rec["found"] = True
                rec["matched_sample_name"] = s.get("name")
                rec["sample"] = s
                rec["db_path"] = os.path.join(s.get("work_dir", ""), s.get("name", ""), s.get("genome", "hg38"), "db", "variants.sqlite")
                rec["meta_path"] = None
                if s.get("work_dir"):
                    meta_dir = meta_subdir or f"meta_{s.get('genome','hg38')}"
                    meta_fn = os.path.join(s.get("work_dir", ""), meta_dir, f"{s.get('name')}.yaml")
                    if os.path.exists(meta_fn):
                        rec["meta_path"] = meta_fn
                rec["sex"] = read_sex_from_meta(s, meta_subdir=meta_subdir)
                out[ident_s] = rec
                continue
        low = ident_s.lower()
        if low in v_ci:
            nm = v_ci[low][0]
            s = idx.get(nm)
            if s:
                rec["found"] = True
                rec["matched_sample_name"] = s.get("name")
                rec["sample"] = s
                rec["db_path"] = os.path.join(s.get("work_dir", ""), s.get("name", ""), s.get("genome", "hg38"), "db", "variants.sqlite")
                rec["meta_path"] = None
                if s.get("work_dir"):
                    meta_dir = meta_subdir or f"meta_{s.get('genome','hg38')}"
                    meta_fn = os.path.join(s.get("work_dir", ""), meta_dir, f"{s.get('name')}.yaml")
                    if os.path.exists(meta_fn):
                        rec["meta_path"] = meta_fn
                rec["sex"] = read_sex_from_meta(s, meta_subdir=meta_subdir)
                out[ident_s] = rec
                continue
        out[ident_s] = rec
    return out


def enrich_rows_and_build_commands(
    rows: List[Dict[str, str]],
    headers: List[str],
    samples: List[Dict[str, Any]],
    id_cols: List[str],
    group_by_karta: bool,
    status_table: Any,
    sci_table: Any,
    base_dir: Optional[str],
    meta_subdir: Optional[str],
    out_script_path: Optional[str] = None
) -> Tuple[List[Dict[str, Any]], List[Tuple[str, str, str]]]:
    """
    Main worker:
      - groups rows by 'Карта' (if present) or each row as its own group
      - for each group, attempts to match identifiers to samples.json
      - if missing, attempts to locate disk DB via status/sci tables and base_dir
      - returns enriched rows and list of extract commands (db, out_param, desc)
    """
    # prepare id list to match
    ids_to_try = []
    for r in rows:
        # prefer provided id_cols in order
        val = None
        for c in id_cols:
            if c in r and str(r[c]).strip():
                val = str(r[c]).strip()
                break
        if not val:
            # fallback: first header
            val = r.get(headers[0], "").strip() if headers else ""
        if val:
            ids_to_try.append(val)
    matches = match_ids_to_samples(samples, ids_to_try, meta_subdir=meta_subdir)

    # grouping
    grouped: Dict[str, List[Dict[str, str]]] = {}
    if group_by_karta:
        for r in rows:
            k = str(r.get("Карта", "")).strip() if "Карта" in r else ""
            if not k:
                k = "UNKNOWN"
            grouped.setdefault(k, []).append(r)
    else:
        # each row its own group (by DNK)
        for r in rows:
            k = str(r.get("ДНК") or r.get("ДНК ", "") or "").strip()
            if not k:
                k = f"ROW_{id(rows)}"
            grouped.setdefault(k, []).append(r)

    enriched: List[Dict[str, Any]] = []
    commands: List[Tuple[str, str, str]] = []
    found_ids_tracker: Set[str] = set()

    for karta, grp in grouped.items():
        folder_name = safe_filename(karta) if karta and karta != "UNKNOWN" else None
        # try to find available dbs for group by checking matches for members first
        group_cmds: List[Tuple[str, str, str]] = []
        # iterate group rows: enrich with matched sample info when available
        for r in grp:
            # pick identifier for this row
            ident = None
            for c in id_cols:
                if c in r and str(r[c]).strip():
                    ident = str(r[c]).strip()
                    break
            if not ident:
                ident = r.get(headers[0], "").strip() if headers else ""
            match = matches.get(ident, {"found": False})
            matched = bool(match and match.get("found"))
            rowout: Dict[str, Any] = dict(r)  # copy original
            rowout["matched_sample_name"] = match.get("matched_sample_name") or ""
            rowout["found"] = "Yes" if matched else "No"
            rowout["db_path"] = match.get("db_path") or ""
            rowout["meta_path"] = match.get("meta_path") or ""
            rowout["genetic_sex"] = match.get("sex") or ""
            # store role normalized
            rowout["_role_norm"] = normalize_role(r.get("Отношение", "") or r.get("Отношение ", ""))
            # folder_name fallback: if no Karte, use family/map from matched sample
            if folder_name is None:
                # Try to derive folder_name from matched sample (if present)
                if matched and match.get("sample"):
                    folder_name = safe_filename(str(match["sample"].get("name") or ident))
            enriched.append(rowout)

        # Attempt to produce extract commands for members missing db or when user requested run script
        if base_dir:
            scmds = search_status_and_sci_for_group(grp, status_table, sci_table, base_dir, found_ids_tracker)
            for dbp, out_param_suffix, desc in scmds:
                out_param = f"{folder_name}/{out_param_suffix}" if folder_name else out_param_suffix
                group_cmds.append((dbp, out_param, desc))

        # Additionally, add commands for matched samples from samples.json (if file exists on disk)
        for r in grp:
            ident = None
            for c in id_cols:
                if c in r and str(r[c]).strip():
                    ident = str(r[c]).strip(); break
            if not ident:
                ident = r.get(headers[0], "").strip() if headers else ""
            match = matches.get(ident, {"found": False})
            if match and match.get("found"):
                s = match.get("sample")
                dbp = os.path.join(s.get("work_dir", ""), s.get("name", ""), s.get("genome", "hg38"), "db", "variants.sqlite")
                if os.path.isfile(dbp):
                    role = normalize_role(r.get("Отношение", ""))
                    out_param = f"{folder_name}/{s.get('name')}_{role}" if folder_name else f"{s.get('name')}_{role}"
                    desc = f"SAMPLES_JSON | {folder_name or s.get('name')}"
                    group_cmds.append((dbp, out_param, desc))
                    found_ids_tracker.add(clean_id(ident))

        # dedupe group_cmds by db & out_param
        seen = set()
        for dbp, outp, desc in group_cmds:
            key = (dbp, outp)
            if key in seen:
                continue
            seen.add(key)
            commands.append((dbp, outp, desc))

    # Optionally write a simple run script with extract_all_variants calls
    if out_script_path:
        with open(out_script_path, "wt", encoding="utf-8") as fh:
            fh.write("#!/bin/bash\nset -e\n\n")
            for dbp, outp, desc in commands:
                fh.write(f'extract_all_variants "{dbp}" "{outp}" "{desc}"\n')
        try:
            st = os.stat(out_script_path)
            os.chmod(out_script_path, st.st_mode | 0o111)
        except Exception:
            pass
        logger.info("Wrote run script: %s (%d commands)", out_script_path, len(commands))

    return enriched, commands


# ----------------------
# CLI
# ----------------------
def build_parser():
    p = argparse.ArgumentParser(prog="samples_utils.py", description="Match filter list to samples.json and optionally generate run script grouped by 'Карта'")
    p.add_argument("-s", "--samples-json", required=True, help="samples.json path")
    p.add_argument("-f", "--filter-list", required=True, help="input filter TSV/CSV with columns like 'ДНК','Карта','ФИО.Код','Отношение', etc.")
    p.add_argument("--filter-id-col", help="Comma-separated column names or 0-based indices to use as identifier (default: try common names or first column)")
    p.add_argument("--no-group-by-karta", dest="group_by_karta", action="store_false", default=True, help="Do not group rows by 'Карта' (default groups by карта)")
    p.add_argument("--status-csv", help="Optional STATUS CSV (to find disk paths like generate.py)")
    p.add_argument("--sci-csv", help="Optional SCI CSV (to find disk paths)")
    p.add_argument("--base-dir", help="Base directory prefix used by validate_and_fix_path to locate variants.sqlite (optional)")
    p.add_argument("--meta-subdir", help="meta subdir name inside work_dir (default meta_<genome>)")
    p.add_argument("--out-json", help="Output enriched JSON")
    p.add_argument("--out-csv", help="Output enriched CSV")
    p.add_argument("--out-script", help="Output run script (simple list of extract_all_variants calls)")
    p.add_argument("--write-vcfs", action="store_true", help="If set, write VCFs directly from found variants.sqlite using embedded SQL")
    p.add_argument("--vcf-outdir", default="vcf_out", help="Directory to write VCFs when --write-vcfs is used")
    p.add_argument("--passport-csv", help="CSV with passport sex mapping (optional). Will add passport_sex and sex_mismatch columns")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv=None):
    p = build_parser()
    args = p.parse_args(argv)
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    samples = load_samples(args.samples_json)
    rows, headers = load_filter_list(args.filter_list)
    if not rows:
        logger.error("No rows parsed from filter list.")
        sys.exit(1)

    # interpret filter-id-col
    id_cols: List[str] = []
    if args.filter_id_col:
        for spec in [s.strip() for s in args.filter_id_col.split(",") if s.strip()]:
            if spec.isdigit():
                idx = int(spec)
                if idx < 0 or idx >= len(headers):
                    logger.error("filter-id-col index %d out of range (0..%d)", idx, len(headers) - 1)
                    sys.exit(1)
                id_cols.append(headers[idx])
            else:
                # try case-insensitive header match
                match = next((h for h in headers if h.lower() == spec.lower()), None)
                if match:
                    id_cols.append(match)
                else:
                    # if not found, still append raw spec and let later code handle missing keys
                    id_cols.append(spec)
    else:
        # defaults try common names
        candidates = ["ФИО.Код", "ДНК", "Variant_ID", "Variant ID", headers[0]]
        for cand in candidates:
            found = next((h for h in headers if h.lower() == str(cand).lower()), None)
            if found:
                id_cols = [found]; break
        if not id_cols:
            id_cols = [headers[0]]

    # load status/sci tables if provided
    status_table = load_csv_table(args.status_csv) if args.status_csv else None
    sci_table = load_csv_table(args.sci_csv) if args.sci_csv else None

    # load passport CSV map if requested
    passport_map = {}
    if args.passport_csv:
        passport_map = load_passport_sex_from_csv(args.passport_csv)
        logger.info("Loaded passport sexes for %d samples from %s", len(passport_map), args.passport_csv)

    enriched, commands = enrich_rows_and_build_commands(
        rows=rows,
        headers=headers,
        samples=samples,
        id_cols=id_cols,
        group_by_karta=args.group_by_karta,
        status_table=status_table,
        sci_table=sci_table,
        base_dir=args.base_dir,
        meta_subdir=args.meta_subdir,
        out_script_path=args.out_script
    )

    # Post-process enriched rows: add passport_sex, normalize, compute mismatch, optionally write VCFs
    os.makedirs(args.vcf_outdir, exist_ok=True) if args.write_vcfs else None
    for r in enriched:
        matched_name = r.get("matched_sample_name") or ""
        passport_raw = passport_map.get(matched_name, "") if matched_name else ""
        r["passport_sex"] = passport_raw
        genetic_raw = r.get("genetic_sex") or ""
        r["genetic_sex"] = genetic_raw
        r["passport_sex_norm"] = normalize_sex_token(passport_raw)
        r["genetic_sex_norm"] = normalize_sex_token(genetic_raw)
        r["sex_mismatch"] = False
        if r["passport_sex_norm"] and r["genetic_sex_norm"] and r["passport_sex_norm"] != r["genetic_sex_norm"]:
            r["sex_mismatch"] = True

        # Optionally write VCF directly from found db_path
        if args.write_vcfs and r.get("found") == "Yes":
            dbp = r.get("db_path") or ""
            if dbp and os.path.isfile(dbp):
                out_vcf = os.path.join(args.vcf_outdir, f"{r.get('matched_sample_name')}.vcf")
                ok, msg = make_vcfs_from_sqlite(dbp, out_vcf, HEADER_VCF)
                r["vcf_written"] = "Yes" if ok else "No"
                r["vcf_msg"] = msg
            else:
                r["vcf_written"] = "No"
                r["vcf_msg"] = "DB_NOT_FOUND"

    # Write outputs
    if args.out_json:
        with open(args.out_json, "wt", encoding="utf-8") as fh:
            json.dump({"enriched": enriched, "commands": commands}, fh, ensure_ascii=False, indent=2)
        logger.info("Wrote JSON output: %s", args.out_json)

    if args.out_csv:
        # columns: original headers + new meta columns
        addcols = ["matched_sample_name", "found", "db_path", "meta_path", "genetic_sex", "passport_sex", "passport_sex_norm", "genetic_sex_norm", "sex_mismatch", "_role_norm"]
        out_headers = list(headers) + addcols
        with open(args.out_csv, "wt", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(out_headers)
            for r in enriched:
                rowvals = [r.get(h, "") for h in headers] + [r.get(c, "") for c in addcols]
                writer.writerow(rowvals)
        logger.info("Wrote CSV output: %s", args.out_csv)

    # summary
    total = len(rows)
    matched = sum(1 for r in enriched if r.get("found") == "Yes")
    logger.info("Processed %d rows. Matched to samples.json: %d. Commands generated: %d", total, matched, len(commands))
    if args.out_script:
        print(f"Run script written: {args.out_script} with {len(commands)} commands.")
    else:
        if commands:
            print("Commands generated (not written to file):")
            for dbp, outp, desc in commands:
                print(f'extract_all_variants "{dbp}" "{outp}" "{desc}"')


if __name__ == "__main__":
    main()