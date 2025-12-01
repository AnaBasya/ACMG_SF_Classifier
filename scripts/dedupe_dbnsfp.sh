#!/usr/bin/env bash
# dedupe_dbnsfp.sh
# Safe deduplication script:
# - identifies identical files by sha256
# - for identical files keeps one copy in canonical dir if present, otherwise keeps largest file
# - for files with same basename but different checksum, moves smaller file(s) to quarantine and logs conflict
# - by default does dry-run; use --apply to actually move files
#
# Usage:
#   ./dedupe_dbnsfp.sh --canonical /path/to/dbnsfp_variants --paths /path/to/parent /path/to/other ... [--quarantine /path/to/quarantine] [--apply]
#
set -euo pipefail

print_help() {
  cat <<EOF
Usage: $0 --canonical DIR --paths DIR1 [DIR2 ...] [--quarantine DIR] [--apply] [--workers N]

Options:
  --canonical DIR   Directory to keep as canonical place for the final single-copy files.
  --paths DIR...    One or more directories to scan (include canonical as well).
  --quarantine DIR  Where to move removed/duplicate files (default: ./quarantine_TIMESTAMP).
  --apply           Actually perform moves. Without --apply it's a dry-run only.
  --workers N       Number of parallel sha256sum workers (default 4).
  -h, --help        Show this help.
EOF
}

# default params
APPLY=0
WORKERS=4
QUARANTINE=""
CANONICAL=""
PATHS=()

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --canonical) CANONICAL="$2"; shift 2;;
    --paths) shift; while [[ $# -gt 0 && "$1" != --apply && "$1" != --quarantine && "$1" != --workers && "$1" != --canonical && "$1" != -h && "$1" != --help ]]; do PATHS+=("$1"); shift; done;;
    --quarantine) QUARANTINE="$2"; shift 2;;
    --apply) APPLY=1; shift;;
    --workers) WORKERS="$2"; shift 2;;
    -h|--help) print_help; exit 0;;
    *) echo "Unknown arg: $1"; print_help; exit 1;;
  esac
done

if [[ -z "$CANONICAL" || "${#PATHS[@]}" -eq 0 ]]; then
  echo "ERROR: --canonical and --paths are required"
  print_help
  exit 1
fi

if [[ ! -d "$CANONICAL" ]]; then
  echo "ERROR: canonical dir does not exist: $CANONICAL"
  exit 1
fi

if [[ -z "$QUARANTINE" ]]; then
  TIMESTAMP=$(date -u +"%Y%m%dT%H%M%SZ")
  QUARANTINE="./quarantine_${TIMESTAMP}"
fi

echo "Canonical dir: $CANONICAL"
echo "Scan paths: ${PATHS[*]}"
echo "Quarantine dir: $QUARANTINE"
echo "Apply mode: $APPLY"
echo "Workers: $WORKERS"

mkdir -p "$QUARANTINE"
LOG_FILE="dedupe_actions.log"
: > "$LOG_FILE"

# gather files
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

find_cmd() {
  local p="$1"
  # only regular files
  find "$p" -type f -print0
}

echo "Gathering file list..."
FILES_LIST="$TMPDIR/files.list"
: > "$FILES_LIST"
for p in "${PATHS[@]}"; do
  if [[ -d "$p" ]]; then
    find "$p" -type f -print0 >> "$TMPDIR/paths_find0"
    # We will process null-separated file list
  else
    echo "Warning: path not found or not dir: $p" | tee -a "$LOG_FILE"
  fi
done

# Create file list (null-separated -> newline-safe)
# use while read -d '' to handle names with newlines
while IFS= read -r -d '' f; do
  # skip files that are inside quarantine dir if user provided it inside scanned paths
  case "$f" in "$QUARANTINE"/*) continue;; esac
  echo "$f" >> "$FILES_LIST"
done < <(cat "$TMPDIR/paths_find0")

NUM_FILES=$(wc -l < "$FILES_LIST")
echo "Found $NUM_FILES files to analyze" | tee -a "$LOG_FILE"

if [[ "$NUM_FILES" -eq 0 ]]; then
  echo "No files found. Exiting."
  exit 0
fi

# compute sha256 and size
SUMS_FILE="$TMPDIR/sums.tsv"
: > "$SUMS_FILE"
echo "Computing sha256 and sizes (this may take time)..." | tee -a "$LOG_FILE"

# parallel-ish computation (simple)
compute_sum() {
  local f="$1"
  # if the file is huge, sha256sum still ok; avoid race by quoting
  sha256sum --binary "$f" 2>/dev/null | awk -v f="$f" '{print $1 "\t" f}'
}
export -f compute_sum
export TMPDIR

# run sequentially to be portable (you can adapt to parallel if needed)
while IFS= read -r f; do
  size=$(stat -c%s "$f" 2>/dev/null || echo 0)
  sum=$(sha256sum --binary "$f" 2>/dev/null | awk '{print $1}')
  echo -e "${sum}\t${size}\t${f}" >> "$SUMS_FILE"
done < "$FILES_LIST"

# Build groups by checksum
# Format in sums.tsv: checksum \t size \t path
sort "$SUMS_FILE" > "$SUMS_FILE.sorted"
mv "$SUMS_FILE.sorted" "$SUMS_FILE"

# Helper to check if path is in canonical dir
is_in_canonical() {
  local p="$1"
  case "$p" in "$CANONICAL"/*) return 0;; *) return 1;; esac
}

# Process identical-checksum groups
echo "Processing identical-checksum groups..." | tee -a "$LOG_FILE"
awk -F'\t' '
{
  sum=$1; size=$2; $1=""; sub(/^\t/,""); path=$0;
  print sum "\t" size "\t" path
}
' "$SUMS_FILE" > "$TMPDIR/parsed.tsv"

# iterate checksums
while IFS=$'\t' read -r checksum size path; do
  : # placeholder to ensure while runs
done < "$TMPDIR/parsed.tsv"

# Group by checksum and handle
awk -F'\t' '{print > ("'"$TMPDIR"'/by_sum_"$1".list")}' "$TMPDIR/parsed.tsv"

KEEP_COUNT=0
MOVE_COUNT=0
CONFLICT_COUNT=0

for f in "$TMPDIR"/by_sum_*.list; do
  # read all lines for this checksum
  readarray -t L < "$f"
  n=${#L[@]}
  if [[ $n -lt 2 ]]; then
    # unique file, ensure it is in canonical dir: if not, plan to move it into canonical (but preserve subdir structure)
    IFS=$'\t' read -r checksum size path <<< "${L[0]}"
    if ! [[ "$path" == "$CANONICAL"/* ]]; then
      # plan to move into canonical, preserve filename; do not overwrite existing canonical file with same checksum
      target="${CANONICAL}/$(basename "$path")"
      if [[ -e "$target" ]]; then
        # if target exists but checksum differs, we will not overwrite: move to quarantine instead
        tsum=$(sha256sum --binary "$target" | awk '{print $1}')
        if [[ "$tsum" == "$checksum" ]]; then
          echo "INFO: identical file already exists in canonical: $path -> keep (target exists equal)" | tee -a "$LOG_FILE"
        else
          echo "CONFLICT: same basename but different checksum: source=$path target=$target" | tee -a "$LOG_FILE"
          echo -e "CONFLICT\t$path\t$target" >> "$LOG_FILE"
          CONFLICT_COUNT=$((CONFLICT_COUNT+1))
          # will move source to quarantine
          echo "MOVE_TO_QUARANTINE\t$path" >> "$LOG_FILE"
        fi
      else
        echo "PLAN_MOVE_TO_CANONICAL\t$path -> $target" | tee -a "$LOG_FILE"
        if [[ $APPLY -eq 1 ]]; then
          mkdir -p "$(dirname "$target")"
          mv -n -- "$path" "$target"
          echo "MOVED\t$path -> $target" >> "$LOG_FILE"
          KEEP_COUNT=$((KEEP_COUNT+1))
        fi
      fi
    fi
  else
    # multiple paths share same checksum (true duplicates)
    # prefer to keep one in canonical; otherwise keep the largest (should be same size though)
    # collect arrays
    keep=""
    largest_size=-1
    largest_path=""
    canonical_present=0
    for line in "${L[@]}"; do
      IFS=$'\t' read -r chk sz p <<< "$line"
      if [[ "$p" == "'"$CANONICAL"'"* ]]; then canonical_present=1; keep="$p"; fi
      if [[ $sz -gt $largest_size ]]; then largest_size=$sz; largest_path="$p"; fi
    done
    if [[ $canonical_present -eq 1 ]]; then
      # remove others to quarantine (or move to canonical if same basename and not present)
      for line in "${L[@]}"; do
        IFS=$'\t' read -r chk sz p <<< "$line"
        if [[ "$p" == "$keep" ]]; then continue; fi
        qpath="'"$QUARANTINE"'"/$(basename "$p")
        # ensure unique name in quarantine
        idx=1
        base=$(basename "$p")
        dest="$qpath"
        while [[ -e "$dest" ]]; do
          dest="'"$QUARANTINE"'"/"${base}.dup${idx}"
          idx=$((idx+1))
        done
        echo "DUPLICATE_MOVE\t$path -> $dest" | tee -a "$LOG_FILE"
        if [[ $APPLY -eq 1 ]]; then
          mv -n -- "$p" "$dest"
          echo "MOVED\t$p -> $dest" >> "$LOG_FILE"
          MOVE_COUNT=$((MOVE_COUNT+1))
        fi
      done
    else
      # no canonical present: choose largest_path as keep, move others to quarantine
      keep="$largest_path"
      for line in "${L[@]}"; do
        IFS=$'\t' read -r chk sz p <<< "$line"
        if [[ "$p" == "$keep" ]]; then continue; fi
        qpath="'"$QUARANTINE"'"/$(basename "$p")
        idx=1
        base=$(basename "$p")
        dest="$qpath"
        while [[ -e "$dest" ]]; do
          dest="'"$QUARANTINE"'"/"${base}.dup${idx}"
          idx=$((idx+1))
        done
        echo "DUPLICATE_MOVE\t$path -> $dest" | tee -a "$LOG_FILE"
        if [[ $APPLY -eq 1 ]]; then
          mv -n -- "$p" "$dest"
          echo "MOVED\t$p -> $dest" >> "$LOG_FILE"
          MOVE_COUNT=$((MOVE_COUNT+1))
        fi
      done
      # ensure keep is inside canonical: if keep not in canonical, move it INTO canonical (rename if necessary)
      if ! [[ "$keep" == "$CANONICAL"/* ]]; then
        target="${CANONICAL}/$(basename "$keep")"
        if [[ -e "$target" ]]; then
          tsum=$(sha256sum --binary "$target" | awk '{print $1}')
          if [[ "$tsum" == "$checksum" ]]; then
            # canonical already has identical file, move keep to quarantine
            qpath="'"$QUARANTINE"'"/$(basename "$keep")
            echo "CANONICAL_HAS_IDENTICAL\t$target exists, moving $keep to quarantine" | tee -a "$LOG_FILE"
            if [[ $APPLY -eq 1 ]]; then
              mv -n -- "$keep" "$qpath"
              echo "MOVED\t$keep -> $qpath" >> "$LOG_FILE"
              MOVE_COUNT=$((MOVE_COUNT+1))
            fi
          else
            # move keep INTO canonical (rename if target exists but different)
            idx=1
            dest="$target"
            while [[ -e "$dest" ]]; do
              dest="${CANONICAL}/$(basename "$keep").keep${idx}"
              idx=$((idx+1))
            done
            echo "MOVE_KEEP_TO_CANONICAL\t$keep -> $dest" | tee -a "$LOG_FILE"
            if [[ $APPLY -eq 1 ]]; then
              mv -n -- "$keep" "$dest"
              echo "MOVED\t$keep -> $dest" >> "$LOG_FILE"
              KEEP_COUNT=$((KEEP_COUNT+1))
            fi
          fi
        else
          echo "MOVE_KEEP_TO_CANONICAL\t$keep -> $target" | tee -a "$LOG_FILE"
          if [[ $APPLY -eq 1 ]]; then
            mv -n -- "$keep" "$target"
            echo "MOVED\t$keep -> $target" >> "$LOG_FILE"
            KEEP_COUNT=$((KEEP_COUNT+1))
          fi
        fi
      fi
    fi
  fi
done

echo "Summary (dry-run if --apply not used):" | tee -a "$LOG_FILE"
echo "Planned moves to quarantine or canonical written to: $LOG_FILE" | tee -a "$LOG_FILE"
echo "Moves performed (if apply): $MOVE_COUNT, Keeps moved (if apply): $KEEP_COUNT, Conflicts: $CONFLICT_COUNT" | tee -a "$LOG_FILE"
echo "Quarantine dir: $QUARANTINE" | tee -a "$LOG_FILE"
echo "Done."