#!/usr/bin/env python3
"""
Compute per-chromosome MAX_AF restricted to regions and preserve AF column names + values.

Outputs TSV with header:
CHROM  POS  REF  ALT  MAX_AF  AF_COLS  AF_VALUES

AF_COLS: comma-separated AF column names (in source order)
AF_VALUES: semicolon-separated values corresponding to AF_COLS (empty string for missing)
"""
import argparse, gzip, bisect, sys

def load_regions_for_chrom(regions_bed, chrom):
    regs=[]
    with open(regions_bed) as fh:
        for ln in fh:
            if not ln.strip(): continue
            parts=ln.split()
            c=parts[0].replace('chr','') if parts[0].startswith('chr') else parts[0]
            if c!=chrom: continue
            try:
                s=int(parts[1]); e=int(parts[2])
            except Exception:
                continue
            regs.append((s,e))
    regs.sort(); merged=[]
    for s,e in regs:
        if not merged or s>merged[-1][1]:
            merged.append([s,e])
        else:
            if e>merged[-1][1]: merged[-1][1]=e
    return merged

def pos_in_regions(pos0, regs):
    i=bisect.bisect_right(regs,[pos0,10**18])-1
    return i>=0 and regs[i][0]<=pos0<regs[i][1]

def open_text(path):
    if path.endswith('.gz') or path.endswith('.bgz'):
        return gzip.open(path,'rt',encoding='utf-8',errors='replace')
    return open(path,'rt',encoding='utf-8',errors='replace')

def try_float(x):
    try:
        if x is None or x=='' or x=='.' or (isinstance(x,str) and x.upper()=='NA'):
            return None
        return float(x)
    except Exception:
        return None

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tsv', required=True)
    p.add_argument('--regions', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--out', required=True)
    args=p.parse_args()

    regs=load_regions_for_chrom(args.regions, args.chrom)
    if not regs:
        with open(args.out,'w') as oh:
            oh.write("CHROM\tPOS\tREF\tALT\tMAX_AF\tAF_COLS\tAF_VALUES\n")
        return

    fh=open_text(args.tsv)
    header_line = fh.readline()
    if not header_line:
        with open(args.out,'w') as oh:
            oh.write("CHROM\tPOS\tREF\tALT\tMAX_AF\tAF_COLS\tAF_VALUES\n")
        return
    header = header_line.rstrip('\n').split('\t')
    hl = [h.lower() for h in header]

    # Standard columns
    try:
        i_ch = hl.index('chrom')
        i_pos = hl.index('pos')
        i_ref = hl.index('ref')
        i_alt = hl.index('alt')
    except ValueError:
        # fallback to first four columns
        i_ch,i_pos,i_ref,i_alt = 0,1,2,3

    # Find AF-like columns: include 'popmax' and any containing 'af' but exclude FAF and 'fraction'
    af_cols = []
    for idx, name in enumerate(header):
        n = name.lower()
        if 'popmax' in n or ('af' in n and 'faf' not in n and 'fraction' not in n):
            af_cols.append((idx, name))
    af_idx = [t[0] for t in af_cols]
    af_names = [t[1] for t in af_cols]

    out = open(args.out, 'wt')
    out.write("CHROM\tPOS\tREF\tALT\tMAX_AF\tAF_COLS\tAF_VALUES\n")

    for ln in fh:
        parts = ln.rstrip('\n').split('\t')
        chrom_raw = parts[i_ch] if i_ch < len(parts) else ''
        chrom = chrom_raw.replace('chr','') if chrom_raw.startswith('chr') else chrom_raw
        if chrom != args.chrom:
            continue
        try:
            pos = int(parts[i_pos])
        except Exception:
            continue
        # pos is 1-based in VCF-derived TSV; convert to 0-based for checking
        if not pos_in_regions(pos-1, regs):
            continue
        ref = parts[i_ref] if i_ref < len(parts) else ''
        alt = parts[i_alt] if i_alt < len(parts) else ''
        vals = []
        nums = []
        for j in af_idx:
            v = parts[j] if j < len(parts) else ''
            fv = try_float(v)
            if fv is None:
                vals.append('')
            else:
                vals.append(str(fv))
                nums.append(fv)
        maxaf = max(nums) if nums else 0.0
        af_cols_str = ",".join(af_names) if af_names else ""
        af_vals_str = ";".join(vals) if vals else ""
        out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{maxaf}\t{af_cols_str}\t{af_vals_str}\n")
    out.close()

if __name__ == '__main__':
    main()