#!/usr/bin/env bash
# run_all.sh -- download gnomAD releases (2.1.1, 3.1, 4.1), restrict to ACMG regions,
# compute per-release MAX_AF including population AF columns, and merge across releases.
#
# Usage:
#   chmod +x run_all.sh
#   ./run_all.sh
#
# Requirements: gsutil, bcftools, tabix, bgzip, python3 (with gzip, bisect; merge script uses only stdlib)
set -euo pipefail

#######################
## CONFIG - edit if needed
#######################
MANE_GTF="${MANE_GTF:-./data/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz}"
ACMG_CSV="${ACMG_CSV:-ACMG_SF_v3.3_full.csv}"
REF_FA="${REF_FA:-/home/anna/ngs-data/ngs/ref/hg38.fa}"
CLINVAR_VCF="${CLINVAR_VCF:-databases/clinvar/clinvar.vcf.gz}"

WORKDIR="${WORKDIR:-./preprocessing}"
MANE_OUT="${MANE_OUT:-./data/mane_out}"
GNOMAD_DIR="${GNOMAD_DIR:-./databases/gnomad}"

# releases to process
RELEASES=("2.1.1" "3.1" "4.1")

# gsutil bucket
GNB="gs://gcp-public-data--gnomad/release"

MANE_FLANK="${MANE_FLANK:-20000}"
PARALLEL_JOBS="${PARALLEL_JOBS:-8}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PY_COMPUTE="${SCRIPT_DIR}/compute_max_af_per_chrom.py"
PY_MERGE="${SCRIPT_DIR}/merge_max_afs.py"
MANE_FINAL_BED="${MANE_OUT}/mane_final.with_ttn_meta.bed"

# HFE position (user specified): GRCh38 6:26092913 with +/-20kb flank added in BED
HFE_POS="6:26092913"

#######################
## END CONFIG
#######################

mkdir -p "$WORKDIR" "$MANE_OUT" "$GNOMAD_DIR"
echo "WORKDIR: $WORKDIR"
echo "MANE_OUT: $MANE_OUT"
echo "GNOMAD_DIR: $GNOMAD_DIR"


echo "2) Create merged ACMG regions bed (0-based) -> acmg_regions.bed"
cut -f1-3 "$MANE_FINAL_BED" | sort -k1,1 -k2,2n > "${MANE_OUT}/mane_intervals.sorted.bed"
bedtools merge -i "${MANE_OUT}/mane_intervals.sorted.bed" > "${MANE_OUT}/acmg_regions.bed"
ACMG_REGIONS="${MANE_OUT}/acmg_regions.bed"

cut -f1 "$ACMG_REGIONS" | sed 's/^chr//' | sort -u > "${WORKDIR}/chroms.txt"
CHROMS_FILE="${WORKDIR}/chroms.txt"
echo "Chromosomes to process:"
cat "$CHROMS_FILE"

echo "3) Ensure helper Python scripts exist and are executable"
# compute_max_af_per_chrom.py and merge_max_afs.py should be in the same directory as this script
if [ ! -f "$PY_COMPUTE" ]; then
  echo "ERROR: missing $PY_COMPUTE"
  exit 1
fi
if [ ! -f "$PY_MERGE" ]; then
  echo "ERROR: missing $PY_MERGE"
  exit 1
fi
chmod +x "$PY_COMPUTE" "$PY_MERGE"

echo "4) Download gnomAD per-release per-chrom VCFs (or TSVs if present)"
# helper: list matches for a gsutil glob pattern and copy matches into DSTDIR
gs_find_and_copy() {
  SRC_PATTERN="$1"; DSTDIR="$2"
  if ! command -v gsutil >/dev/null 2>&1; then
    return 1
  fi
  # list matches (gsutil ls accepts globs)
  matches=$(gsutil ls "${SRC_PATTERN}" 2>/dev/null || true)
  if [ -z "$matches" ]; then
    return 1
  fi
  echo "Found matches for pattern: ${SRC_PATTERN}"
  # copy each match (avoid eval issues)
  echo "$matches" | while IFS= read -r obj; do
    # skip empty lines
    [ -z "$obj" ] && continue
    echo "Copying $obj -> $DSTDIR/"
    gsutil -m cp "${obj}" "${DSTDIR}/"
  done
  return 0
}

# For each chromosome and release, try release-specific patterns (based on your gsutil ls output)
while read CHR; do
  echo "=== CHR $CHR ==="
  for REL in "${RELEASES[@]}"; do
    mkdir -p "${GNOMAD_DIR}/${REL}"
    DST="${GNOMAD_DIR}/${REL}"
    GOT=0

    # Per-release specific candidate patterns (prefer VCFs; TSVs rarely available publicly for these releases)
    CAND_PATTERNS=()
    if [ "$REL" = "2.1.1" ]; then
      # release 2.1.1 naming: gnomad.genomes.r2.1.1.sites.<n>.vcf.bgz and exomes similar
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/gnomad.genomes.r${REL}.sites.${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/gnomad.genomes.r${REL}.sites.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/exomes/gnomad.exomes.r${REL}.sites.${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/exomes/gnomad.exomes.r${REL}.sites.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/*sites*chr${CHR}*.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/exomes/*sites*chr${CHR}*.vcf.bgz")
    elif [ "$REL" = "3.1" ]; then
      # release 3.1 naming (from your ls): gnomad.genomes.v3.1.sites.chr*.vcf.bgz and hgdp_1kg_subset variants
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/gnomad.genomes.v3.1.sites.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/gnomad.genomes.v3.1.hgdp_1kg_subset.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/*chr${CHR}*.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/*/*chr${CHR}*.vcf.bgz")
    elif [ "$REL" = "4.1" ]; then
      # release 4.1 naming: joint/genomes/exomes under vcf/joint vcf/genomes vcf/exomes
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/joint/gnomad.joint.v4.1.sites.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/genomes/gnomad.genomes.v4.1.sites.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/exomes/gnomad.exomes.v4.1.sites.chr${CHR}.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/*/*chr${CHR}*.vcf.bgz")
    else
      # generic fallback
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/*/*chr${CHR}*.vcf.bgz")
      CAND_PATTERNS+=("${GNB}/${REL}/vcf/*/*sites*${CHR}*.vcf.bgz")
    fi

    # Try patterns in order
    for P in "${CAND_PATTERNS[@]}"; do
      if gs_find_and_copy "$P" "${DST}"; then GOT=1; fi
    done

    if [ "$GOT" -eq 0 ]; then
      echo "WARNING: no VCF found for release $REL chr $CHR (patterns tried)."
      continue
    fi

    # After copying, find a VCF file in DST for this chr (pick first)
    vfile=$(ls "${DST}"/*chr${CHR}*.vcf.bgz "${DST}"/*.${CHR}.vcf.bgz "${DST}"/*sites.*.${CHR}*.vcf.bgz "${DST}"/*.vcf.bgz 2>/dev/null | head -n1 || true)
    if [ -z "$vfile" ]; then
      echo "ERROR: expected VCF in $DST but none found (after copy). Skipping."
      continue
    fi
    # Ensure index exists (copy may have .tbi file or we can tabix)
    if [ ! -f "${vfile}.tbi" ] && command -v tabix >/dev/null 2>&1; then
      echo "Indexing $vfile with tabix"
      tabix -p vcf "$vfile" || true
    fi

    # Convert VCF -> TSV (AF fields), create TSV_OUT name
    TSV_OUT="${DST}/$(basename "${vfile%.*}").af.tsv.gz"
    echo "Converting VCF -> TSV (AF fields) -> $TSV_OUT"
    if command -v bcftools >/dev/null 2>&1; then
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_popmax\t%INFO/AF\t%INFO/AF_joint\t%INFO/AF_genomes\t%INFO/AF_exomes\t%INFO/AF_afr\t%INFO/AF_amr\t%INFO/AF_asj\t%INFO/AF_eas\t%INFO/AF_fin\t%INFO/AF_nfe\t%INFO/AF_sas\t%INFO/AF_oth\n' "$vfile" | bgzip -c > "${TSV_OUT}.tmp"
      printf 'CHROM\tPOS\tREF\tALT\tAF_popmax\tAF\tAF_joint\tAF_genomes\tAF_exomes\tAF_afr\tAF_amr\tAF_asj\tAF_eas\tAF_fin\tAF_nfe\tAF_sas\tAF_oth\n' | cat - <(zcat "${TSV_OUT}.tmp") | bgzip -c > "$TSV_OUT"
      rm -f "${TSV_OUT}.tmp"
    else
      echo "bcftools missing; cannot convert VCF to TSV automatically for $vfile"
      continue
    fi

  done
done < "$CHROMS_FILE"

echo "5) Compute per-chrom MAX_AF restricted to ACMG regions"
for REL in "${RELEASES[@]}"; do
  for f in "${GNOMAD_DIR}/${REL}/"*.tsv.bgz "${GNOMAD_DIR}/${REL}"/*.af.tsv.gz; do
    [ -f "$f" ] || continue
    fname=$(basename "$f")
    chrom=""
    if [[ "$fname" =~ chr([0-9XYM]+) ]]; then chrom="${BASH_REMATCH[1]}"; fi
    if [ -z "$chrom" ]; then chrom=$(zcat "$f" | awk 'NR>1{print $1; exit}'); fi
    [ -n "$chrom" ] || { echo "Cannot determine chrom for $f; skipping"; continue; }
    out="${GNOMAD_DIR}/${REL}/chr${chrom}.max_af.tsv"
    echo "Processing release $REL file $f -> $out (chrom: $chrom)"
    python3 "$PY_COMPUTE" --tsv "$f" --regions "$ACMG_REGIONS" --chrom "$chrom" --out "$out"
  done
done

echo "6) Merge per-release results (compute global MAX_AF and per-population maxima + release sources)"
python3 "$PY_MERGE" --input-dir "$GNOMAD_DIR" --out "${WORKDIR}/gnomad_all_releases_max_af_per_pop.tsv"

echo "Done. Final merged: ${WORKDIR}/gnomad_all_releases_max_af_per_pop.tsv"