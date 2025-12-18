import pandas as pd
import sys
import os

# === SETTINGS ===
INPUT_TABLE_PATH = "./data/input_list.tsv"       
SCI_CSV_PATH = "./data/sci-fullset.csv"   
STATUS_CSV_PATH = "./data/status.csv"     
# Using os.path.expanduser to correctly handle the tilde symbol in Python
BASE_DIR = os.path.expanduser("~/ngs-data/smk/GENOMES/BTK/") 
OUTPUT_SCRIPT = "run_jobs.sh"

# --- Functions ---

def clean_col_names(df):
    """Strips whitespace from column names."""
    if not df.empty: df.columns = df.columns.str.strip()
    return df

def find_col(df, candidates):
    """Finds a column name case-insensitively from a list of candidates."""
    if df.empty: return None
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower: return cols_lower[cand.lower()]
    return None

def normalize_role(role_raw):
    """Normalizes family role to English (proband, mother, father, sibling)."""
    if pd.isna(role_raw): return "sample"
    r = str(role_raw).lower().strip()
    if any(x in r for x in ["пробанд", "обратившийся", "proband", "prob"]): return "proband"
    if any(x in r for x in ["мать", "mother", "mom"]): return "mother"
    if any(x in r for x in ["отец", "father", "dad"]): return "father"
    if any(x in r for x in ["сибс", "брат", "сестра", "sib"]): return "sibling"
    return "sample"

def get_folder_type(study_name):
    """Determines subdirectory type: 'panels' or 'smk'."""
    s = str(study_name).lower() if pd.notna(study_name) else ""
    return "panels" if "панель" in s else "smk"

def clean_id(val):
    """Removes leading zeros and special chars from DNA ID."""
    s = str(val).strip().split('.')[0]
    return s.lstrip('0') if s.isdigit() else s

def safe_filename(name):
    """Replaces slashes with underscores for directory safety."""
    return str(name).replace('/', '_').replace('\\', '_').strip()

def validate_and_fix_path(base_dir, dir_name, sub_type, dna_padded):
    """
    Checks for the physical existence of the DB file.
    If not found, tries switching between 'smk' and 'panels'.
    Returns the valid absolute path or None.
    """
    # Construct primary path
    # base_dir should be expanded (no ~)
    path_primary = f"{base_dir}{dir_name}/{sub_type}/{dna_padded}/hg38/db/variants.sqlite"
    
    if os.path.isfile(path_primary):
        return path_primary
    
    # If not found, try the alternate sub_type
    path_alt = None
    if sub_type == "smk":
        path_alt = f"{base_dir}{dir_name}/panels/{dna_padded}/hg38/db/variants.sqlite"
    elif sub_type == "panels":
        path_alt = f"{base_dir}{dir_name}/smk/{dna_padded}/hg38/db/variants.sqlite"
        
    if path_alt and os.path.isfile(path_alt):
        return path_alt
        
    return None

def load_data():
    """Loads CSV/TSV data into DataFrames."""
    try:
        # fillna("") is important for groupby operations later
        df_inp = pd.read_csv(INPUT_TABLE_PATH, sep='\t', dtype=str).fillna("")
    except:
        print(f"CRITICAL: File not found: {INPUT_TABLE_PATH}"); sys.exit(1)

    try:
        df_st = clean_col_names(pd.read_csv(STATUS_CSV_PATH, dtype=str))
    except:
        print(f"WARN: File not found: {STATUS_CSV_PATH}"); df_st = pd.DataFrame()

    try:
        df_sci = clean_col_names(pd.read_csv(SCI_CSV_PATH, dtype=str))
    except:
        print(f"WARN: File not found: {SCI_CSV_PATH}"); df_sci = pd.DataFrame()

    return df_inp, df_sci, df_st

# --- SCI Table Search Logic ---
def search_sci(fam_id, group, df_sci, found_ids_tracker):
    commands = []
    id_col = find_col(df_sci, ['ID', 'sample'])
    dir_col = find_col(df_sci, ['dir_name', 'directory', 'path'])
    if not (id_col and dir_col): return []

    for _, row in group.iterrows():
        dna = clean_id(row['ДНК'])
        if not dna: continue

        # Filter by DNA ID (clean match)
        mask = df_sci[id_col].apply(clean_id) == dna
        match = df_sci[mask]
        if match.empty: continue
        
        rec = match.iloc[0]
        d_name = str(rec[dir_col]).strip()
        if d_name.lower() in ['не получен', 'nan', '', 'none'] or len(d_name) < 2: continue

        dna_pad = dna.zfill(12)
        role = normalize_role(row['Отношение'])
        stype = get_folder_type(row['НазваниеИсследования'])
        
        # VERIFY FILE ON DISK
        valid_path = validate_and_fix_path(BASE_DIR, d_name, stype, dna_pad)
        if not valid_path:
            print(f"[MISSING FILE] ID: {dna} (Fam: {fam_id}) in {d_name}")
            continue

        f_id_real = fam_id if fam_id not in ["UNKNOWN", ""] else f"SCI_{dna}"
        folder_name = safe_filename(f_id_real)
        
        # Use short DNA ID for filename
        out_param = f"{folder_name}/{dna}_{role}"
        
        commands.append(f'extract_all_variants "{valid_path}" "{out_param}" "SCI | {f_id_real}"')
        found_ids_tracker.add(dna)

    return commands

# --- STATUS Table Search Logic ---
def search_status(fam_id, group, df_status, found_ids_tracker):
    commands = []
    fam_col = find_col(df_status, ['fam', 'map', 'семья', 'карта'])
    id_col = find_col(df_status, ['ID', 'DNA'])
    dir_col = find_col(df_status, ['dir_name', 'dir'])
    rel_col = find_col(df_status, ['rel', 'role', 'отношение'])
    stat_col = find_col(df_status, ['deseased2', 'status', 'type'])

    if not (fam_col and id_col and dir_col): return []

    # 1. Determine list of families to process
    families_to_process = set()
    
    if fam_id and fam_id != "UNKNOWN":
        fam_str = str(fam_id).strip()
        if not df_status[df_status[fam_col].astype(str).str.strip() == fam_str].empty:
            families_to_process.add(fam_str)
    
    # If no families found via Map ID, search via DNA IDs
    if not families_to_process:
        inp_ids = [clean_id(x) for x in group['ДНК'].tolist() if x]
        mask = df_status[id_col].apply(clean_id).isin(inp_ids)
        found_rows = df_status[mask]
        if not found_rows.empty:
            uniq_fams = found_rows[fam_col].astype(str).str.strip().unique()
            for f in uniq_fams: families_to_process.add(f)

    if not families_to_process: return []

    # 2. Process each found family
    for real_fam in families_to_process:
        fam_rows = df_status[df_status[fam_col].astype(str).str.strip() == real_fam]
        if fam_rows.empty: continue

        # Is it a trio/duo/quadro?
        is_trio = False
        if stat_col:
            val = str(fam_rows.iloc[0][stat_col]).lower()
            if any(x in val for x in ['трио', 'trio', 'дуо', 'duo', 'квадро']):
                is_trio = True
        
        # If Trio, fetch everyone. If Mono, only fetch matching DNA IDs from input.
        targets = fam_rows if is_trio else fam_rows[fam_rows[id_col].apply(clean_id).isin([clean_id(x) for x in group['ДНК'].tolist()])]

        for _, row in targets.iterrows():
            d_raw = clean_id(row[id_col]) # Short ID
            d_pad = d_raw.zfill(12)       # Padded ID for path
            d_name = str(row[dir_col]).strip()
            
            if d_name.lower() in ['не получен', 'nan', ''] or len(d_name) < 2: continue

            stype = get_folder_type(group.iloc[0]['НазваниеИсследования'])
            
            # VERIFY FILE ON DISK
            valid_path = validate_and_fix_path(BASE_DIR, d_name, stype, d_pad)
            if not valid_path:
                print(f"[MISSING FILE] ID: {d_raw} (Map: {real_fam}) in {d_name}")
                continue

            # Determine role
            if is_trio and rel_col:
                role = normalize_role(row[rel_col])
            else:
                m = group[group['ДНК'].apply(clean_id) == d_raw]
                role = normalize_role(m.iloc[0]['Отношение']) if not m.empty else normalize_role(row.get(rel_col, ''))

            # Define output names
            out_fam_label = fam_id if (fam_id and fam_id != "UNKNOWN") else real_fam
            folder_name = safe_filename(out_fam_label)
            out_param = f"{folder_name}/{d_raw}_{role}"
            
            commands.append(f'extract_all_variants "{valid_path}" "{out_param}" "STATUS | {out_fam_label}"')
            
            # Track if we found a requested sample
            if d_raw in [clean_id(x) for x in group['ДНК'].tolist()]:
                found_ids_tracker.add(d_raw)

    return commands

def main():
    inp, sci, status = load_data()
    inp['Карта'] = inp['Карта'].replace('', 'UNKNOWN')
    
    found_ids_set = set()
    all_cmds = []

    try:
        grouped = inp.groupby(['Карта', 'НазваниеИсследования'])
    except KeyError as e:
        print(f"CRITICAL: Missing columns in input table: {e}"); sys.exit(1)

    print("Starting file verification and path generation...")

    for (fam, study), group in grouped:
        cmds = []
        is_sci = str(study).strip() == "Научное исследование"
        
        # Search strategy
        if fam == "UNKNOWN":
            c = search_status(fam, group, status, found_ids_set)
            if not c: c = search_sci(fam, group, sci, found_ids_set)
            cmds = c
        elif is_sci:
            c = search_sci(fam, group, sci, found_ids_set)
            if not c: c = search_status(fam, group, status, found_ids_set)
            cmds = c
        else:
            c = search_status(fam, group, status, found_ids_set)
            if not c: c = search_sci(fam, group, sci, found_ids_set)
            cmds = c
        
        if cmds: all_cmds.extend(cmds)

    # Remove duplicates
    unique_cmds = list(dict.fromkeys(all_cmds))
    
    # Write Bash script
    with open(OUTPUT_SCRIPT, 'w') as f:
        f.write(BASH_HEADER)
        f.write(f"\n# Verified Jobs: {len(unique_cmds)}\n")
        f.write("\n".join(unique_cmds))
        f.write("\n\necho 'All jobs finished. Check folders.'")

    input_ids_all = set([clean_id(x) for x in inp['ДНК'].tolist() if str(x).strip()])
    print("-" * 40)
    print(f"Samples requested: {len(input_ids_all)}")
    print(f"Verified paths found: {len(unique_cmds)} (includes family members if trio)")
    print(f"Script created: {OUTPUT_SCRIPT}")

# === BASH ===
BASH_HEADER = r"""#!/bin/bash
set -e

# Создаем заголовок VCF
cat > header.vcf << 'EOF'
##fileformat=VCFv4.2
##source=SQLiteExport
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=HGVSc,Number=1,Type=String,Description="HGVS cDNA notation">
##INFO=<ID=HGVSp,Number=1,Type=String,Description="HGVS protein notation">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="Reference call">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=LowGQ,Description="Low genotype quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Fraction">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample
EOF

extract_all_variants() {
    local db="$1"
    local output_rel_path="$2"
    local desc="$3"
    
    local out_vcf="${output_rel_path}.vcf"
    local dir_name=$(dirname "$out_vcf")
    mkdir -p "$dir_name"

    if [ -f "$out_vcf" ]; then echo "SKIP: $out_vcf"; return 0; fi
    
    if [ ! -f "$db" ]; then 
        echo "   [ERR] DB not found during execution: $db"; return 0; 
    fi
    
    # === AUTO-DETECT COLUMNS ===
    local cols=$(sqlite3 "$db" "PRAGMA table_info(vcf);")    
    # 1. Info: Gene
    local c_gene="''"
    if echo "$cols" | grep -q "|Gene|"; then c_gene="Gene";
    elif echo "$cols" | grep -q "|gene|"; then c_gene="gene"; fi
    
    # 2. Info: HGVSc
    local c_c="''"
    if echo "$cols" | grep -q "|HGVSc|"; then c_c="HGVSc";
    elif echo "$cols" | grep -q "|hgvsc|"; then c_c="hgvsc"; fi
    
    # 3. Info: HGVSp
    local c_p="''"
    if echo "$cols" | grep -q "|HGVSp|"; then c_p="HGVSp";
    elif echo "$cols" | grep -q "|hgvsp|"; then c_p="hgvsp"; fi

    # --- ОПРЕДЕЛЕНИЕ РЕАЛЬНЫХ ДАННЫХ ГЕНОТИПА ---
    
    # 4. GT (Genotype)
    local c_gt="'./.'"
    if echo "$cols" | grep -q "|GT|"; then c_gt="GT"; fi

    # 5. DP (Depth)
    local c_dp="'.'"
    if echo "$cols" | grep -q "|DP|"; then c_dp="DP"; fi

    # 6. GQ (Genotype Quality)
    local c_gq="'.'"
    if echo "$cols" | grep -q "|GQ|"; then c_gq="GQ"; fi
    
    # 7. AF
    local c_ad="'.,.'"
    local c_af="'.'"
    
    if echo "$cols" | grep -q "|AD|"; then 
        c_ad="AD"
        # If AD exists, try to calculate AF from it if not present? 
        # SQLite split is hard, assume caller provided AF if they provided AD.
    elif echo "$cols" | grep -q "|RO|" && echo "$cols" | grep -q "|AO|"; then
        # Build AD from RO,AO
        c_ad="RO || ',' || AO"
        
        # Calculate AF = AO / (RO + AO)
        # Check for division by zero
        c_af="CASE WHEN (RO + AO) > 0 THEN CAST(AO AS FLOAT) / (RO + AO) ELSE 0.0 END"
    fi
    
    # If explicit AF column exists, use it override
    if echo "$cols" | grep -q "|AF|"; then c_af="AF"; fi

    echo ">> $desc -> $out_vcf (detected Info: $c_gene, $c_c; Genotype: $c_gt, $c_dp, $c_ad, $c_af)"
    
    local temp_vcf="${out_vcf}_temp"
    
    # Added AF to FORMAT string and value list
    sqlite3 -separator $'\t' "$db" "
SELECT 
    CHROM, POS, vid as ID, REF, ALT, 
    COALESCE(QUAL, 50), COALESCE(FILTER, 'PASS'), 
    'GENE=' || COALESCE($c_gene, '') || ';HGVSc=' || COALESCE($c_c, '') || ';HGVSp=' || COALESCE($c_p, '') as INFO,
    'GT:DP:GQ:AD:AF' as FORMAT,
    COALESCE($c_gt, './.') || ':' || COALESCE($c_dp, '.') || ':' || COALESCE($c_gq, '.') || ':' || COALESCE($c_ad, '.,.') || ':' || COALESCE($c_af, '0.0') as Sample
FROM vcf;" > "$temp_vcf"
    
    cat header.vcf > "$out_vcf"
    sort -k1,1V -k2,2n "$temp_vcf" >> "$out_vcf"
    rm -f "$temp_vcf"
}
"""

if __name__ == "__main__":
    main()
