#!/usr/bin/env bash
# Download and index dbNSFP variant files (dbNSFP4.1a_variant.chr{CHR}.gz) from Zenodo
# Usage:
#   ./download_and_index_dbnsfp.sh /path/to/output_dir  # uses default chromosome list
#   ./download_and_index_dbnsfp.sh /path/to/outdir 1 2 3 X
#
# Requirements: curl, bgzip, tabix, python3, gunzip (gzip)
set -euo pipefail

OUTDIR="${1:-./}"
shift || true

# Default list requested by you:
DEFAULT_CHRS=(1 2 3 5 6 7 9 10 11 12 13 14 15 16 17 18 19 22 X)

if [ "$#" -ge 1 ]; then
    CHRS=("$@")
else
    CHRS=("${DEFAULT_CHRS[@]}")
fi

BASE_URL="https://zenodo.org/records/4323592/files"
mkdir -p "$OUTDIR"

# check prerequisites
for cmd in curl python3 tabix bgzip gunzip; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 2
    fi
done

echo "Target directory: $OUTDIR"
echo "Chromosomes to download: ${CHRS[*]}"
echo "Base URL: $BASE_URL"
echo ""

download_file() {
    local url="$1"
    local out="$2"
    if [ -s "$out" ]; then
        echo "Skip download (exists): $out"
        return 0
    fi
    echo "Downloading: $url -> $out"
    # use curl with retries
    curl -L --retry 6 --retry-delay 3 --retry-max-time 60 -o "$out" "$url"
    if [ ! -s "$out" ]; then
        echo "Download failed or file empty: $out" >&2
        return 1
    fi
    return 0
}

detect_columns_and_index() {
    local file="$1"
    echo "Detecting CHR/POS columns and indexing: $file"

    # get header (first non-empty line)
    header=$(gunzip -c "$file" 2>/dev/null | sed -n '1,1p' || true)
    if [ -z "$header" ]; then
        echo "Failed to read header from $file" >&2
        return 2
    fi

    # Use python to detect column indices (1-based)
    read -r chrom_idx pos_idx <<< "$(python3 - "$header" <<'PY' 2>/dev/null
import sys
hdr = sys.stdin.read().strip().split('\t')
# normalize
low = [h.lower() for h in hdr]
chrom_idx = None
pos_idx = None
# common names for chrom
for i, v in enumerate(low, start=1):
    if v in ("chrom","chr","chromosome","#chrom","#chr"):
        chrom_idx = i
        break
# fallback: column that contains 'chr' in name (rare)
if chrom_idx is None:
    for i,v in enumerate(low, start=1):
        if v.startswith("chr") or "chrom" in v:
            chrom_idx = i; break
# find position column
for i, v in enumerate(low, start=1):
    if v in ("pos","position","start","coordinate"):
        pos_idx = i
        break
if pos_idx is None:
    # sometimes 'pos' is second column if file is variant-level
    if len(low) >= 2:
        pos_idx = 2
    else:
        pos_idx = 2
if chrom_idx is None:
    chrom_idx = 1
# output
print(chrom_idx, pos_idx)
PY
)"

    echo "  -> using CHR column: $chrom_idx, POS column: $pos_idx"

    # Try tabix indexing, if fails, repack to bgzip then index
    set +e
    tabix -f -s "$chrom_idx" -b "$pos_idx" -e "$pos_idx" "$file" >/dev/null 2>&1
    rc=$?
    set -e
    if [ $rc -eq 0 ]; then
        echo "Indexed OK: ${file}.tbi"
        return 0
    fi

    # Repack to bgzip (in-place)
    echo "Tabix indexing failed (likely not BGZF). Repacking to BGZF and indexing..."
    mv "$file" "${file}.orig"
    gunzip -c "${file}.orig" | bgzip -c > "$file"
    tabix -f -s "$chrom_idx" -b "$pos_idx" -e "$pos_idx" "$file"
    echo "Indexed after repacking: ${file}.tbi"
    return 0
}

for CHR in "${CHRS[@]}"; do
    FNAME="dbNSFP4.1a_variant.chr${CHR}.gz"
    URL="${BASE_URL}/${FNAME}"
    OUT="${OUTDIR}/${FNAME}"
    if download_file "$URL" "$OUT"; then
        # ensure file is non-empty
        if [ ! -s "$OUT" ]; then
            echo "ERROR: downloaded file is empty: $OUT" >&2
            continue
        fi
        # index
        if detect_columns_and_index "$OUT"; then
            echo "OK: $OUT indexed"
        else
            echo "WARNING: failed to index $OUT" >&2
        fi
    else
        echo "ERROR: download failed for $URL" >&2
    fi
    echo ""
done

echo "All done. Files saved to: $OUTDIR"
echo "Next: run VEP with --plugin dbNSFP,<path-to-file>,FIELD1,FIELD2,... using the exact field names from the file header."
