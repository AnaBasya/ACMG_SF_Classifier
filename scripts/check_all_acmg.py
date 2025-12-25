import pandas as pd

# Load your ACMG table
df = pd.read_csv('ACMG_SF_v3.3_full.csv', sep=';', dtype=str)
print(f"Total ACMG genes: {len(df)}")

# Known coordinates for common ACMG genes (simplified)
acmg_coords = {
    'BRCA1': ('chr17', 43044294, 43125483),
    'BRCA2': ('chr13', 32315473, 32400266),
    'TP53': ('chr17', 7668421, 7687550),
    'MLH1': ('chr3', 37034841, 37107662),
    'MSH2': ('chr2', 47641410, 47709935),
    'MSH6': ('chr2', 48010187, 48034092),
    'APC': ('chr5', 112043202, 112181936),
    'MYH7': ('chr14', 23416000, 23436000),
    'MYBPC3': ('chr11', 47340000, 47370000),
    'KCNQ1': ('chr11', 2594000, 2910000),
    'KCNH2': ('chr7', 150600000, 150700000),
    'SCN5A': ('chr3', 38560000, 38630000),
    'HFE': ('chr6', 26091141, 26107124),
    'BTD': ('chr3', 128198313, 128242590),
    'CFTR': ('chr7', 117292062, 117384629),
}

print("\nChecking if your variants overlap any ACMG genes:")
for gene, (chrom, start, end) in acmg_coords.items():
    cmd = f"bcftools view -H -r {chrom}:{start}-{end} ./data/result.vcf.gz 2>/dev/null | wc -l"
    import subprocess
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    count = int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0
    if count > 0:
        print(f"  ✓ {gene}: {count} variants")
    else:
        print(f"  ✗ {gene}: 0 variants")
