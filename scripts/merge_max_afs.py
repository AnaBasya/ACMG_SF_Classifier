#!/usr/bin/env python3
"""
Merge per-release chr*.max_af.tsv files.

Output columns:
VAR    MAX_AF    RELEASE_OF_MAX    AF_<pop1>    AF_<pop1>_RELEASE    AF_<pop2>    AF_<pop2>_RELEASE ...

Where VAR is chrom:pos:ref:alt and AF_<pop> is the maximal AF observed for that pop across releases,
and AF_<pop>_RELEASE is the release that provided that maximum.
"""
import argparse, glob, os, csv
from collections import defaultdict

p=argparse.ArgumentParser()
p.add_argument('--input-dir', default='gnomad', help='dir containing per-release subdirs with *.max_af.tsv')
p.add_argument('--out', default='gnomad_all_releases_max_af_per_pop.tsv')
args = p.parse_args()

pattern = os.path.join(args.input_dir, '*', '*.max_af.tsv')
files = sorted(glob.glob(pattern))
if not files:
    print("No per-release *.max_af.tsv files found under", args.input_dir)
    raise SystemExit(0)

overall_max = {}
release_of_max = {}
per_pop_max = defaultdict(dict)            # key -> {popname: value}
per_pop_release = defaultdict(dict)        # key -> {popname: release}
pop_order = []                             # keep discovery order
pop_seen = set()

for fn in files:
    rel = os.path.basename(os.path.dirname(fn))
    with open(fn) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            chrom = row.get('CHROM') or row.get('chrom') or ''
            pos = row.get('POS') or row.get('pos') or ''
            ref = row.get('REF') or row.get('ref') or ''
            alt = row.get('ALT') or row.get('alt') or ''
            if not (chrom and pos and ref and alt):
                continue
            key = f"{chrom}:{pos}:{ref}:{alt}"

            # release-specific overall maxaf
            try:
                maxaf_val = float(row.get('MAX_AF') or row.get('Max_AF') or row.get('max_af') or 0.0)
            except Exception:
                maxaf_val = 0.0
            if key not in overall_max or maxaf_val > overall_max[key]:
                overall_max[key] = maxaf_val
                release_of_max[key] = rel

            # parse AF_COLS and AF_VALUES
            af_cols_raw = row.get('AF_COLS') or row.get('af_cols') or ''
            af_vals_raw = row.get('AF_VALUES') or row.get('af_values') or ''
            if af_cols_raw and af_vals_raw:
                cols = [c.strip() for c in af_cols_raw.split(',') if c.strip()]
                vals = [v.strip() for v in af_vals_raw.split(';')]
                for i, col in enumerate(cols):
                    v = vals[i] if i < len(vals) else ''
                    try:
                        fv = float(v) if v not in ('', '.', 'NA', None) else None
                    except Exception:
                        fv = None
                    pname = col.strip()
                    if pname and pname not in pop_seen:
                        pop_seen.add(pname); pop_order.append(pname)
                    if fv is not None:
                        prev = per_pop_max[key].get(pname, None)
                        if (prev is None) or (fv > prev):
                            per_pop_max[key][pname] = fv
                            per_pop_release[key][pname] = rel

# Build header
pop_list = pop_order
out_header = ["VAR", "MAX_AF", "RELEASE_OF_MAX"]
for pnm in pop_list:
    safe = pnm.replace(" ", "_")
    out_header.append("AF_"+safe)
    out_header.append("AF_"+safe+"_RELEASE")

with open(args.out, 'wt') as out:
    out.write("\t".join(out_header) + "\n")
    for k in sorted(overall_max.keys()):
        row = [k, str(overall_max[k]), release_of_max.get(k, "")]
        pops = per_pop_max.get(k, {})
        prs = per_pop_release.get(k, {})
        for pnm in pop_list:
            v = pops.get(pnm, "")
            rel = prs.get(pnm, "")
            row.append(str(v) if v != "" else "")
            row.append(rel if rel else "")
        out.write("\t".join(row) + "\n")

print("Wrote", args.out, "variants:", len(overall_max))