#!/usr/bin/env bash
# Быстро скачивает указанные хромосомы gnomAD v4.1 (joint) и их .tbi
# Использует gsutil -m cp (предпочтительно). fallback: aria2c по HTTPS.
# Пример:
#   ./download_gnomad_fast.sh 1,2,17,X  # скачает chr1,chr2,chr17,chrX

set -euo pipefail

CHRS_CSV="${1:-}"   # список хромосом через запятую, обязателен
OUTDIR="${2:-./databases/gnomad/3.1/vcf}"
THREADS="${3:-8}"   # число параллельных задач (для aria2c/gsutil)
ARIA_CONN_PER_FILE="${4:-8}" # aria2c: соединений на файл

if [[ -z "$CHRS_CSV" ]]; then
  echo "Usage: $0 CHRS_CSV [OUTDIR] [THREADS] [ARIA_CONN_PER_FILE]"
  echo "Example: $0 1,2,17,X ./data 8 8"
  exit 2
fi

mkdir -p "$OUTDIR"
IFS=',' read -ra CHRS <<< "$CHRS_CSV"

GS_PREFIX="gs://gcp-public-data--gnomad/release/3.1/vcf/genomes"
HTTPS_PREFIX="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes"
TEMPLATE="gnomad.genomes.v3.1.sites.chr%s.vcf.bgz"

files=()
for c in "${CHRS[@]}"; do
  fname=$(printf "$TEMPLATE" "$c")
  files+=("$fname")
done

# prefer gsutil if present
if command -v gsutil >/dev/null 2>&1; then
  echo "Using gsutil -m cp (recommended). Parallelism: $THREADS"
  # build list of gs:// URIs
  uris=()
  for f in "${files[@]}"; do
    uris+=("${GS_PREFIX}/${f}")
    uris+=("${GS_PREFIX}/${f}.tbi")
  done
  # shell brace expansion isn't needed here because we build list
  echo "Copying ${#uris[@]} files to $OUTDIR ..."
  # gsutil -m cp supports multiple sources; it parallelizes transfers
  gsutil -m cp -r "${uris[@]}" "$OUTDIR"/
  echo "Done."
  exit 0
fi

# fallback: use aria2c if available (fast multi-connection download over HTTPS)
if command -v aria2c >/dev/null 2>&1; then
  echo "gsutil not found — using aria2c (HTTPS multi-connection) fallback."
  listfile="$(mktemp)"
  for f in "${files[@]}"; do
    echo "${HTTPS_PREFIX}/${f}" >> "$listfile"
    echo "${HTTPS_PREFIX}/${f}.tbi" >> "$listfile"
  done
  echo "Starting aria2c with $THREADS parallel downloads and $ARIA_CONN_PER_FILE connections per file..."
  aria2c -i "$listfile" -j "$THREADS" -x "$ARIA_CONN_PER_FILE" -s "$ARIA_CONN_PER_FILE" -d "$OUTDIR" --allow-overwrite=true
  rm -f "$listfile"
  echo "Done."
  exit 0
fi

# last fallback: wget in parallel (less efficient)
if command -v parallel >/dev/null 2>&1 && command -v wget >/dev/null 2>&1; then
  echo "Using GNU parallel + wget as last resort."
  jobs=()
  for f in "${files[@]}"; do
    jobs+=("${HTTPS_PREFIX}/${f}")
    jobs+=("${HTTPS_PREFIX}/${f}.tbi")
  done
  printf "%s\n" "${jobs[@]}" | parallel -j "$THREADS" wget -c -q -P "$OUTDIR" {}
  echo "Done."
  exit 0
fi

echo "No suitable downloader found. Install gsutil (recommended) or aria2c or GNU parallel + wget."
exit 3
