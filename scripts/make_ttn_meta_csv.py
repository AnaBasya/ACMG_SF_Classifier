#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create TTN meta-exon CSV used by ACMG SF classifier

Purpose
-------
This utility takes an input file describing TTN exon coordinates and PSI (percent spliced in)
and normalizes it to a CSV with the following columns:
  - chr
  - start
  - end
  - exon_id (optional, empty string if not provided)
  - psi_meta (integer 0-100)
  - lof_allowed (1 if psi_meta >= 70 else 0)

Input formats
-------------
- TSV/CSV with header containing at least chr/start/end and a psi column (psi, psi_meta, or similar).
- PSI values can be in 0-1 (fraction) or 0-100 (percent) and are normalized.

Output
------
- CSV suitable for use as --ttn-meta in the classifier.

Example
-------
python make_ttn_meta_csv.py --input ttn_input.tsv --output ttn_meta.csv

"""
from __future__ import annotations
import argparse
import pandas as pd
import os
import sys

def normalize_psi(val):
    try:
        v = float(val)
        if v <= 1.0:
            v = v * 100.0
        return float(v)
    except Exception:
        return None

def make_meta_df(input_path: str) -> pd.DataFrame:
    # Read file autodetecting delimiter
    try:
        df = pd.read_csv(input_path, sep=None, engine='python', dtype=str)
    except Exception as e:
        raise RuntimeError(f"Failed to read input {input_path}: {e}")

    cols = {c.lower(): c for c in df.columns}
    chrcol = cols.get('chr') or cols.get('chrom') or cols.get('chromosome')
    startcol = cols.get('start') or cols.get('exon_start') or cols.get('start_pos')
    endcol = cols.get('end') or cols.get('exon_end') or cols.get('end_pos')
    psicol = None
    for k, v in cols.items():
        if 'psi' in k:
            psicol = v
            break
    exoncol = cols.get('exon') or cols.get('exon_id') or cols.get('exon_name') or ""

    if not chrcol or not startcol or not endcol or not psicol:
        raise RuntimeError("Input must contain chr, start, end and a psi column (headers case-insensitive).")

    rows = []
    for _, r in df.iterrows():
        chrom = r.get(chrcol)
        try:
            start = int(float(r.get(startcol)))
            end = int(float(r.get(endcol)))
        except Exception:
            continue
        exon_id = r.get(exoncol) if exoncol else ""
        psi_raw = r.get(psicol)
        psi = normalize_psi(psi_raw)
        if psi is None:
            continue
        psi_int = int(round(psi))
        lof_allowed = 1 if psi_int >= 70 else 0
        rows.append({"chr": chrom, "start": start, "end": end, "exon_id": exon_id, "psi_meta": psi_int, "lof_allowed": lof_allowed})
    out = pd.DataFrame(rows)
    out = out.sort_values(['chr','start'])
    return out

def main():
    p = argparse.ArgumentParser(description="Create TTN meta-exon CSV")
    p.add_argument("--input", required=True, help="Input TSV/CSV with chr,start,end and psi column")
    p.add_argument("--output", required=True, help="Output CSV path")
    args = p.parse_args()

    out_df = make_meta_df(args.input)
    out_df.to_csv(args.output, index=False)
    print(f"Wrote {len(out_df)} rows to {args.output}")

if __name__ == "__main__":
    main()
