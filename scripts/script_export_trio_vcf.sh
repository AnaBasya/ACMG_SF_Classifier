#!/bin/bash

set -e  # Прерывать выполнение при ошибках

# Создаем ПОЛНЫЙ VCF заголовок со всеми определениями
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

# Функция для извлечения ВСЕХ вариантов
extract_all_variants() {
    local db_path="$1"
    local output_file="$2"
    local sample_name="$3"
    
    echo "Извлекаем ВСЕ варианты для $sample_name..."
    echo "Путь к базе: $db_path"
    
    if [ ! -f "$db_path" ]; then
        echo "Ошибка: База данных не найдена: $db_path" >&2
        return 1
    fi
    
    # Извлекаем ВСЕ варианты с реалистичным распределением генотипов
    sqlite3 -separator $'\t' "$db_path" "
SELECT 
    CHROM,
    POS, 
    vid as ID, 
    REF, 
    ALT, 
    COALESCE(QUAL, 50) as QUAL, 
    COALESCE(FILTER, 'PASS') as FILTER, 
    'GENE=' || COALESCE(Gene, '') || ';HGVSc=' || COALESCE(HGVSc, '') || ';HGVSp=' || COALESCE(HGVSp, '') as INFO,
    'GT:DP:GQ:AD' as FORMAT,
    -- Реалистичное распределение генотипов:
    CASE (abs(random()) % 100)
        WHEN 0 THEN '0/0:35:99:35,0'    -- референсный гомозиготный
        WHEN 1 THEN '0/0:40:99:40,0'    -- референсный гомозиготный
        WHEN 2 THEN '1/1:25:80:0,25'    -- альтернативный гомозиготный
        WHEN 3 THEN '1/1:30:85:0,30'    -- альтернативный гомозиготный
        ELSE '0/1:28:90:14,14'          -- гетерозиготный (большинство)
    END as Sample
FROM vcf;" > "${output_file}_data.vcf"
    
    # Объединяем заголовок и данные
    cat header.vcf "${output_file}_data.vcf" > "${output_file}.vcf"
    
    # Проверяем результат
    line_count=$(wc -l < "${output_file}.vcf")
    variant_count=$((line_count - 38))  # 37 header строк + 1 заголовок
    echo "  ✓ Создан: ${output_file}.vcf ($variant_count вариантов)"
    
    # Покажем распределение FILTER значений
    echo "  FILTER значения:"
    grep -v '^#' "${output_file}.vcf" | cut -f7 | sort | uniq -c | while read count filter; do
        echo "    $filter: $count вариантов"
    done
}

# Используем полные пути
BASE_DIR="$HOME/ngs-data/smk/GENOMES/BTK"

echo "Создание VCF файлов с ПОЛНЫМ заголовком (включая FILTER)..."
echo "=========================================================="

# Извлекаем данные
extract_all_variants \
    "$BASE_DIR/2025/02/2/panels/000007118610/hg38/db/variants.sqlite" \
    "proband_test_no_sort" \
    "пробанда"

extract_all_variants \
    "$BASE_DIR/2025/02/3/panels/000007068710/hg38/db/variants.sqlite" \
    "father_test_no_sort" \
    "отца"

extract_all_variants \
    "$BASE_DIR/2025/02/3/panels/000007068700/hg38/db/variants.sqlite" \
    "mother_test_no_sort" \
    "матери"

vcf-sort proband_test_no_sort.vcf > proband_test.vcf
vcf-sort father_test_no_sort.vcf > father_test.vcf
vcf-sort mother_test_no_sort.vcf > mother_test.vcf

# Очистка временных файлов
rm -f header.vcf *data.vcf *no_sort.vcf

echo ""
echo "=========================================="
echo "Готово! Созданы файлы с полным заголовком:"
echo "  - proband_test.vcf"
echo "  - father_test.vcf" 
echo "  - mother_test.vcf"
echo ""
echo "Проверка структуры:"
echo "-------------------"

# Проверяем наличие всех определений
for file in proband_test.vcf father_test.vcf mother_test.vcf; do
    if [ -f "$file" ]; then
        echo "✓ $file:"
        echo "  CONTIG определений: $(grep -c '^##contig' "$file")"
        echo "  FILTER определений: $(grep -c '^##FILTER' "$file")"
        echo "  Всего вариантов: $(grep -v '^#' "$file" | wc -l)"
        echo "  FILTER значения: $(grep -v '^#' "$file" | cut -f7 | sort | uniq | tr '\n' ' ')"
    fi
    echo ""
done
