#!/usr/bin/env bash
#
# build_mane_bed.sh (updated)
# Same as before but genome file prefix normalized to match BED (auto add/remove "chr").
#
set -euo pipefail

if [ $# -lt 3 ]; then
  echo "Usage: $0 MANE_GTF.gz ACMG_GENES_FILE OUT_DIR [chr|nochr] [flank_bp] [REF_FA] [HFE_POS] [CLINVAR_VCF]" >&2
  exit 2
fi

MANE_GTF="$1"
ACMG_GENES_FILE="$2"
OUTDIR="$3"
PREFIX_MODE="${4:-chr}"
FLANK="${5:-20000}"
REF_FA="${6:-}"
DEFAULT_HFE="6:26092913"
HFE_POS_ARG="${7:-$DEFAULT_HFE}"
CLINVAR_VCF="${8:-}"

mkdir -p "$OUTDIR"

RAW_BED="${OUTDIR}/mane_raw.bed"
NORM_BED="${OUTDIR}/mane_norm.bed"
SORTED_BED="${OUTDIR}/mane_sorted.bed"
ACMG_BED="${OUTDIR}/mane_acmg_only.bed"
FINAL_BED="${OUTDIR}/mane_final.sorted.bed"
GENES_TXT="${OUTDIR}/acmg_genes.txt"

echo "build_mane_bed.sh starting"
echo " MANE_GTF: $MANE_GTF"
echo " ACMG_GENES_FILE: $ACMG_GENES_FILE"
echo " OUTDIR: $OUTDIR"
echo " PREFIX_MODE: $PREFIX_MODE"
echo " FLANK: $FLANK"
echo " REF_FA: ${REF_FA:-<none>}"
echo " HFE_POS_ARG: ${HFE_POS_ARG:-<none>}"
echo " CLINVAR_VCF: ${CLINVAR_VCF:-<none>}"
echo

if [ ! -f "$MANE_GTF" ]; then
  echo "ERROR: MANE GTF file not found: $MANE_GTF" >&2
  exit 3
fi
if [ ! -f "$ACMG_GENES_FILE" ]; then
  echo "ERROR: ACMG genes file not found: $ACMG_GENES_FILE" >&2
  exit 4
fi

# Extract gene symbols
awk -F'[;, \t]' 'NR==1{ h=$1; if(tolower(h)~/gene/) { next } }
{ g=$1; gsub(/["'\''"]/,"",g); gsub(/^[ \t]+|[ \t]+$/,"",g); if(g!="") print g }' "$ACMG_GENES_FILE" | sort -u > "$GENES_TXT"
if [ ! -s "$GENES_TXT" ]; then
  cut -f1 "$ACMG_GENES_FILE" | sed '1d' | sed 's/"//g' | sed '/^\s*$/d' | sort -u > "$GENES_TXT"
fi
echo "  extracted $(wc -l < "$GENES_TXT") genes"

# Extract MANE transcripts
echo "Extracting MANE transcripts -> $RAW_BED"
zcat -f "$MANE_GTF" | awk -F'\t' '
$3=="transcript" && ($0 ~ /MANE_Select/ || $0 ~ /MANE_Plus_Clinical/) {
  attr = $9; gene=""; tx="";
  if (match(attr, /gene_name "[^"]+"/)) { s = substr(attr, RSTART, RLENGTH); sub(/^gene_name "/,"",s); sub(/"$/,"",s); gene=s }
  if (match(attr, /transcript_id "[^"]+"/)) { t = substr(attr, RSTART, RLENGTH); sub(/^transcript_id "/,"",t); sub(/"$/,"",t); tx=t }
  if (gene=="" || tx=="") next;
  chrom=$1; start=$4; end=$5; bedstart = start-1;
  printf("%s\t%d\t%d\t%s|%s\n", chrom, bedstart, end, gene, tx)
}' > "${RAW_BED}"

if [ ! -s "${RAW_BED}" ]; then
  echo "ERROR: No MANE transcript lines extracted. Please verify MANE GTF and tags." >&2
  exit 5
fi
echo "  extracted $(wc -l < ${RAW_BED}) MANE transcript lines"

# Normalize chromosome names in BED
echo "Normalizing chromosome names -> $NORM_BED (mode: $PREFIX_MODE)"
if [ "$PREFIX_MODE" = "chr" ]; then
  awk 'BEGIN{OFS="\t"}{ if ($1 !~ /^chr/) $1="chr"$1; print $0 }' "${RAW_BED}" > "${NORM_BED}"
else
  awk 'BEGIN{OFS="\t"}{ sub(/^chr/,"",$1); print $0 }' "${RAW_BED}" > "${NORM_BED}"
fi

# Sort + unique
echo "Sorting and unique -> $SORTED_BED"
sort -k1,1 -k2,2n "${NORM_BED}" | uniq > "${SORTED_BED}"

# Filter to ACMG genes
echo "Filtering to ACMG genes -> $ACMG_BED"
awk -v genes_file="$GENES_TXT" 'BEGIN{ while((getline g<genes_file)>0) genes[g]=1 }
{ if($4=="") next; n=$4; split(n,a,"|"); gene=a[1]; if(gene in genes) print $0 }' "${SORTED_BED}" > "${ACMG_BED}"
echo "  ACMG transcripts kept: $(wc -l < "${ACMG_BED}")"

# Add HFE position (explicit)
HFE_ADDED=0
if [ -n "${HFE_POS_ARG}" ]; then
  if [[ "${HFE_POS_ARG}" =~ : ]]; then
    CHR=${HFE_POS_ARG%%:*}; POS=${HFE_POS_ARG##*:}; START=$((POS-1)); END=$POS
    echo -e "${CHR}\t${START}\t${END}\tHFE|HFE_C282Y" >> "${ACMG_BED}"
    HFE_ADDED=1
    echo "Added HFE position ${HFE_POS_ARG} to ACMG BED"
  fi
fi

# Prepare genome file from REF_FA (if provided)
if [ -n "$REF_FA" ] && command -v samtools >/dev/null 2>&1; then
  FAI="${REF_FA}.fai"
  if [ ! -f "$FAI" ]; then
    echo "Indexing REF_FA with samtools faidx..."
    samtools faidx "$REF_FA"
  fi
  GENOME_FILE="${OUTDIR}/genome_from_fai.txt"
  cut -f1,2 "$FAI" > "${GENOME_FILE}.orig"

  # Detect prefix mismatch between BED and genome file and normalize genome file accordingly
  # Get first contig from BED and from genome orig
  BED_FIRST=$(awk 'NR==1{print $1; exit}' "$ACMG_BED" | sed 's/:.*//')
  GEN_FIRST=$(awk 'NR==1{print $1; exit}' "${GENOME_FILE}.orig")

  # If BED has "chr" and genome doesn't -> add "chr"
  if [[ "$BED_FIRST" == chr* && "$GEN_FIRST" != chr* ]]; then
    echo "Adding 'chr' prefix to genome file contig names to match BED..."
    awk '{print "chr"$1"\t"$2}' "${GENOME_FILE}.orig" > "${GENOME_FILE}"
  # If BED lacks "chr" and genome has -> remove "chr"
  elif [[ "$BED_FIRST" != chr* && "$GEN_FIRST" == chr* ]]; then
    echo "Removing 'chr' prefix from genome file contig names to match BED..."
    awk '{g=$1; sub(/^chr/,"",g); print g"\t"$2}' "${GENOME_FILE}.orig" > "${GENOME_FILE}"
  else
    # no mismatch
    mv "${GENOME_FILE}.orig" "${GENOME_FILE}"
  fi
else
  GENOME_FILE=""
fi

# Ensure ACMG_BED has the same chromosome naming (add/remove chr) as chosen by PREFIX_MODE
echo "Ensuring ACMG_BED chromosome names are consistent with PREFIX_MODE ($PREFIX_MODE)"
if [ "$PREFIX_MODE" = "chr" ]; then
  awk 'BEGIN{OFS="\t"}{ if ($1 !~ /^chr/) $1="chr"$1; print $0 }' "${ACMG_BED}" > "${ACMG_BED}.norm"
else
  awk 'BEGIN{OFS="\t"}{ sub(/^chr/,"",$1); print $0 }' "${ACMG_BED}" > "${ACMG_BED}.norm"
fi
mv "${ACMG_BED}.norm" "${ACMG_BED}"

# Apply flank using bedtools (with the normalized genome file if present)
if [ -n "$GENOME_FILE" ] && command -v bedtools >/dev/null 2>&1; then
  echo "Applying bedtools slop with genome file ${GENOME_FILE} (flank ${FLANK})..."
  TMP_SLOP="${OUTDIR}/mane_acmg_slop.bed"
  bedtools slop -i "${ACMG_BED}" -g "${GENOME_FILE}" -b "${FLANK}" > "${TMP_SLOP}"
  mv "${TMP_SLOP}" "${FINAL_BED}"
else
  echo "bedtools or genome file not available; applying heuristic flank (clamped at 0) of ${FLANK} bp"
  awk -v b="${FLANK}" 'BEGIN{OFS="\t"}{s=$2-b; if(s<0) s=0; print $1,s,$3+b,$4}' "${ACMG_BED}" > "${FINAL_BED}"
fi

# Final sort & uniq
sort -k1,1 -k2,2n "${FINAL_BED}" | uniq > "${FINAL_BED}.tmp" && mv "${FINAL_BED}.tmp" "${FINAL_BED}"

echo
echo "MANE ACMG BED generated: ${FINAL_BED}"
echo "Lines: $(wc -l < ${FINAL_BED})"
head -n 20 "${FINAL_BED}"
if [ "${HFE_ADDED}" -eq 1 ]; then
  echo "HFE position added."
else
  echo "HFE position NOT added."
fi
echo "Done."
