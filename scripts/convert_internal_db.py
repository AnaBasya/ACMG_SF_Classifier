import pandas as pd
from pyliftover import LiftOver

# Пути
INPUT_FILE = "/home/anna/anna/ACMG_SF_Classifier/databases/internal/ngs_data_annotations.csv"
OUTPUT_FILE = "/home/anna/anna/ACMG_SF_Classifier/databases/internal/ngs_data_annotations_hg38.csv"

# Загружаем цепочку конвертации hg19 -> hg38
lo = LiftOver('hg19', 'hg38')

def get_hg38_pos(chrom, pos):
    # pyliftover использует 0-based координаты, а в CSV скорее всего 1-based.
    # Обычно VCF/CSV хранят 1-based. pyliftover возвращает список совпадений.
    new_coords = lo.convert_coordinate(chrom, int(pos) - 1) 
    if new_coords:
        # Возвращаем +1, чтобы вернуть в 1-based формат
        return int(new_coords[0][1]) + 1
    return None

df = pd.read_csv(INPUT_FILE)
new_rows = []

print(f"Обработка {len(df)} строк...")

for _, row in df.iterrows():
    # Пропускаем, если уже hg38 (на всякий случай)
    if row.get('ref') == 'hg38':
        new_rows.append(row)
        continue

    chrom = row['CHROM']
    old_pos = row['POS']
    
    new_pos = get_hg38_pos(chrom, old_pos)
    
    if new_pos:
        row['POS'] = new_pos
        row['ref'] = 'hg38' # Обновляем метку
        new_rows.append(row)
    else:
        print(f"Не удалось перенести: {chrom}:{old_pos}")

df_new = pd.DataFrame(new_rows)
df_new.to_csv(OUTPUT_FILE, index=False)
print(f"Готово! Сохранено в: {OUTPUT_FILE}")
