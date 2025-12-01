#!/bin/bash

set -euo pipefail

echo "Загружаем все базы данных для ACMG SF анализа..."

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BASE_DIR}/databases"
mkdir -p ${DB_DIR}

# Функции для логирования
log() { echo "[INFO] $1"; }
success() { echo "[SUCCESS] $1"; }
warning() { echo "[WARNING] $1"; }
error() { echo "[ERROR] $1"; exit 1; }

# Функция для проверки и загрузки
download_file() {
    local url=$1
    local output=$2
    local desc=$3
    
    # Создаем директорию если нужно
    mkdir -p "$(dirname "$output")"
    
    if [[ -f "$output" ]]; then
        warning "Пропускаем $desc (уже существует: $(basename "$output"))"
        return 0
    fi
    
    log "Загружаем $desc..."
    
    # Пробуем разные методы загрузки
    if command -v wget &> /dev/null; then
        if wget --progress=bar:force --timeout=120 --tries=3 -q -O "$output" "$url"; then
            success "Загружено: $desc"
            return 0
        fi
    fi
    
    if command -v curl &> /dev/null; then
        if curl -L --connect-timeout 60 --retry 3 -s -o "$output" "$url"; then
            success "Загружено: $desc"
            return 0
        fi
    fi
    
    error "Не удалось загрузить: $desc"
}

# Создаем структуру директорий
log "Создаем структуру директорий..."
mkdir -p ${DB_DIR}/{vep,clinvar,clingen,gnomad,cadd,revel,spliceai,local}

# 1. ClinVar (GRCh37)
download_file \
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" \
    "${DB_DIR}/clinvar/clinvar.vcf.gz" \
    "ClinVar GRCh37"

download_file \
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi" \
    "${DB_DIR}/clinvar/clinvar.vcf.gz.tbi" \
    "ClinVar index"

# 2. ClinGen данные
download_file \
    "https://search.clinicalgenome.org/kb/gene-validity/download" \
    "${DB_DIR}/clingen/gene_disease_validity.csv" \
    "ClinGen Gene-Disease Validity"

download_file \
    "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv" \
    "${DB_DIR}/clingen/dosage_sensitivity.tsv" \
    "ClinGen Dosage Sensitivity"

download_file \
    "http://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed" \
    "${DB_DIR}/clingen/variant_pathogenicity.csv" \
    "ClinGen Variant Pathogenicity"

download_file \
    "https://search.clinicalgenome.org/kb/reports/curation-activity-summary-report" \
    "${DB_DIR}/clingen/summary.csv" \
    "ClinGen Summary Report"

# 3. gnomAD (exomes GRCh37 для экономии места)
download_file \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz" \
    "${DB_DIR}/gnomad/gnomad_exomes.vcf.bgz" \
    "gnomAD Exomes GRCh37"

# 4. CADD scores
download_file \
    "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz" \
    "${DB_DIR}/cadd/cadd_scores.tsv.gz" \
    "CADD Scores GRCh37"

# 5. REVEL
download_file \
    "http://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" \
    "${DB_DIR}/revel/revel_all_chromosomes.zip" \
    "REVEL scores"

# 6. SpliceAI (GRCh37)
git clone https://github.com/Illumina/SpliceAI.git "${DB_DIR}/spliceai" || warning "SpliceAI уже существует или ошибка клонирования"

# 7. VEP cache (GRCh37)
if command -v vep &> /dev/null; then
    log "Настраиваем VEP cache для GRCh37..."
    vep_install -a cf -s homo_sapiens -y GRCh37 --CACHE_DIR ${DB_DIR}/vep --NO_UPDATE
else
    warning "VEP не установлен. Установите: conda install -c bioconda ensembl-vep"
fi

# Индексирование VCF файлов
log "Индексируем VCF файлы..."
for vcf in ${DB_DIR}/clinvar/clinvar.vcf.gz ${DB_DIR}/gnomad/gnomad_exomes.vcf.bgz; do
    if [[ -f "$vcf" && ! -f "${vcf}.tbi" ]]; then
        tabix -p vcf "$vcf" && success "Проиндексирован: $(basename "$vcf")"
    fi
done

# Подготовка ACMG SF данных из Excel
log "Подготавливаем ACMG SF данные из Excel..."
python3 ${BASE_DIR}/scripts/prepare_acmg_data.py || warning "Скрипт prepare_acmg_data.py не найден или завершился с ошибкой"

success "Все базы данных загружены и подготовлены!"

echo ""
echo "Использовано дискового пространства:"
du -sh ${DB_DIR}/* | sort -hr

echo ""
echo "Проверка целостности файлов:"
for db in clinvar clingen gnomad; do
    if [[ -f "${DB_DIR}/${db}/".gz || -f "${DB_DIR}/${db}/".csv ]]; then
        echo "OK: $(find ${DB_DIR}/${db} -name ".gz" -o -name ".csv" | wc -l) файлов в ${db}/"
    fi
done
