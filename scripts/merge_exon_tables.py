#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Создание объединенной таблицы экзонов для NMD-предсказания
ИСПРАВЛЕННАЯ ВЕРСИЯ - с правильным сохранением колонки Gene
"""

import pandas as pd
import numpy as np
import re
import argparse
import os
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_exons_csv(exons_csv_path: str) -> pd.DataFrame:
    """Загрузка таблицы exons.csv с геномными координатами"""
    logger.info(f"Загружаем exons.csv: {exons_csv_path}")
    
    # Читаем CSV
    df = pd.read_csv(exons_csv_path)
    logger.info(f"Загружено {len(df)} строк из exons.csv")
    logger.info(f"Колонки: {list(df.columns)}")
    
    # Извлекаем номер экзона из колонки Exon
    def extract_exon_number(exon_str):
        if pd.isna(exon_str):
            return None
        s = str(exon_str).lower()
        # Ищем номер экзона: exon2 -> 2, exon3 -> 3
        match = re.search(r'exon(\d+)', s)
        if match:
            return int(match.group(1))
        # Если просто число
        try:
            return int(s)
        except:
            return None
    
    df['Exon_Number'] = df['Exon'].apply(extract_exon_number)
    
    # Создаем колонку Gene из Gene.refGene
    if 'Gene.refGene' in df.columns:
        df['Gene'] = df['Gene.refGene']
        logger.info(f"Создана колонка 'Gene' из 'Gene.refGene'")
    elif 'Gene' in df.columns:
        df['Gene'] = df['Gene']
    else:
        df['Gene'] = ''
        logger.warning("Нет колонки с названием гена в exons.csv")
    
    # Создаем новый DataFrame с нужными колонками
    result = pd.DataFrame()
    
    # Обязательные колонки
    required_cols = ['Transcript_ID', 'Exon_Number', 'Gene']
    for col in required_cols:
        if col in df.columns:
            result[col] = df[col]
        else:
            logger.error(f"Отсутствует обязательная колонка: {col}")
            result[col] = ''
    
    # Координаты
    if 'Genomic_start' in df.columns:
        result['Genomic_start'] = df['Genomic_start']
    elif 'Start' in df.columns:
        result['Genomic_start'] = df['Start']
    else:
        logger.error("Нет колонки с начальной координатой экзона")
        result['Genomic_start'] = 0
    
    if 'Genomic_end' in df.columns:
        result['Genomic_end'] = df['Genomic_end']
    elif 'End' in df.columns:
        result['Genomic_end'] = df['End']
    else:
        logger.error("Нет колонки с конечной координатой экзона")
        result['Genomic_end'] = 0
    
    # Дополнительные колонки
    optional_cols = ['CDS_start', 'CDS_end', 'Exon_Length']
    for col in optional_cols:
        if col in df.columns:
            result[col] = df[col]
    
    # Сортируем
    if not result.empty:
        result = result.sort_values(['Gene', 'Transcript_ID', 'Exon_Number']).reset_index(drop=True)
    
    logger.info(f"Обработано exons.csv: {len(result)} строк")
    logger.info(f"Уникальных генов: {result['Gene'].nunique()}")
    logger.info(f"Пример генов: {list(result['Gene'].unique())[:10]}")
    
    return result

def load_exon_cds_tsv(cds_tsv_path: str) -> pd.DataFrame:
    """Загрузка таблицы exon_df_cds.tsv с CDS информацией"""
    logger.info(f"Загружаем exon_df_cds.tsv: {cds_tsv_path}")
    
    # Читаем TSV
    try:
        df = pd.read_csv(cds_tsv_path, sep='\t')
    except Exception as e:
        logger.error(f"Ошибка при чтении {cds_tsv_path}: {e}")
        # Попробуем с другим разделителем
        df = pd.read_csv(cds_tsv_path, sep='\t', engine='python')
    
    logger.info(f"Загружено {len(df)} строк из exon_df_cds.tsv")
    logger.info(f"Колонки: {list(df.columns)}")
    
    # Приводим названия колонок к стандартному виду
    df.columns = [col.strip() for col in df.columns]
    
    # Переименовываем для единообразия
    rename_dict = {}
    
    # Transcript ID
    if 'Transcript_ID' in df.columns:
        rename_dict['Transcript_ID'] = 'Transcript_ID'
    elif 'transcript' in df.columns:
        rename_dict['transcript'] = 'Transcript_ID'
    
    # Exon number
    if 'Exon_Number' in df.columns:
        rename_dict['Exon_Number'] = 'Exon_Number'
    elif 'number_cds_exon' in df.columns:
        rename_dict['number_cds_exon'] = 'Exon_Number'
    
    # Остальные колонки
    col_mapping = {
        'Exon_length': 'Exon_length',
        'CDS_Exon_length': 'CDS_Exon_length',
        'Total_Exons': 'Total_Exons',
        'total_cds_exons': 'Total_Exons',
        'CDS_length': 'CDS_length',
        'Protein_length': 'Protein_length',
        'Chrom': 'Chrom',
        'strand': 'strand'
    }
    
    for source, target in col_mapping.items():
        if source in df.columns:
            rename_dict[source] = target
    
    # Переименовываем
    if rename_dict:
        df = df.rename(columns=rename_dict)
        logger.info(f"Переименованы колонки: {rename_dict}")
    
    # Оставляем только нужные колонки
    keep_cols = ['Transcript_ID', 'Exon_Number', 'Exon_length', 'CDS_Exon_length',
                 'Total_Exons', 'CDS_length', 'Protein_length', 'Chrom', 'strand']
    
    available_cols = [col for col in keep_cols if col in df.columns]
    df = df[available_cols].copy()
    
    # Убедимся, что номера экзонов - целые числа
    if 'Exon_Number' in df.columns:
        df['Exon_Number'] = pd.to_numeric(df['Exon_Number'], errors='coerce')
    
    # Добавляем базовую версию Transcript ID для слияния
    df['Transcript_ID_base'] = df['Transcript_ID'].astype(str).str.split('.').str[0]
    
    logger.info(f"Обработано exon_df_cds.tsv: {len(df)} строк")
    logger.info(f"Уникальных транскриптов: {df['Transcript_ID'].nunique()}")
    
    return df

def load_transcript_map(transcript_map_path: str) -> pd.DataFrame:
    """Загрузка таблицы соответствия RefSeq-Ensembl"""
    logger.info(f"Загружаем refseq2enst_mane.tsv: {transcript_map_path}")
    
    try:
        df = pd.read_csv(transcript_map_path, sep='\t')
    except Exception as e:
        logger.error(f"Ошибка при чтении {transcript_map_path}: {e}")
        df = pd.read_csv(transcript_map_path, sep='\t', engine='python')
    
    logger.info(f"Загружено {len(df)} строк из transcript map")
    logger.info(f"Колонки: {list(df.columns)}")
    
    # Приводим названия колонок к стандартному виду
    df.columns = [col.strip() for col in df.columns]
    
    # Переименовываем для единообразия
    rename_dict = {}
    
    # RefSeq ID
    if 'RNA_nucleotide_accession.version' in df.columns:
        rename_dict['RNA_nucleotide_accession.version'] = 'RefSeq_ID'
    elif 'transcript' in df.columns:
        rename_dict['transcript'] = 'RefSeq_ID'
    
    # Ensembl ID
    if 'Ensembl_rna_identifier' in df.columns:
        rename_dict['Ensembl_rna_identifier'] = 'Ensembl_ID'
    elif 'Ensembl_rna' in df.columns:
        rename_dict['Ensembl_rna'] = 'Ensembl_ID'
    
    # Gene
    if 'gene' in df.columns:
        rename_dict['gene'] = 'Gene'
    elif 'GeneID' in df.columns:
        rename_dict['GeneID'] = 'Gene'
    elif 'Gene' in df.columns:
        rename_dict['Gene'] = 'Gene'
    
    # MANE select
    if 'mane_select' in df.columns:
        rename_dict['mane_select'] = 'MANE_Select'
    
    if rename_dict:
        df = df.rename(columns=rename_dict)
    
    # Оставляем только нужные колонки
    keep_cols = ['RefSeq_ID', 'Ensembl_ID', 'Gene', 'MANE_Select']
    available_cols = [col for col in keep_cols if col in df.columns]
    
    if 'Gene' not in available_cols:
        logger.warning("В transcript map нет информации о генах")
        available_cols = [col for col in keep_cols if col != 'Gene']
    
    df = df[available_cols].copy()
    
    # Очищаем значения
    for col in df.columns:
        df[col] = df[col].astype(str).str.strip()
    
    # Для RefSeq_ID создаем базовую версию
    df['RefSeq_ID_base'] = df['RefSeq_ID'].astype(str).str.split('.').str[0]
    
    logger.info(f"Обработано transcript map: {len(df)} строк")
    if 'Gene' in df.columns:
        logger.info(f"Уникальных генов: {df['Gene'].nunique()}")
    
    return df

def merge_datasets_simple(exons_df: pd.DataFrame, cds_df: pd.DataFrame, transcript_map_df: pd.DataFrame) -> pd.DataFrame:
    """Простое объединение всех трёх таблиц"""
    logger.info("Объединяем таблицы (простой метод)...")
    
    # 1. Сначала подготовим exons_df
    exons_df = exons_df.copy()
    exons_df['Transcript_ID_base'] = exons_df['Transcript_ID'].astype(str).str.split('.').str[0]
    
    # 2. Подготовим cds_df
    cds_df = cds_df.copy()
    cds_df['Transcript_ID_base'] = cds_df['Transcript_ID'].astype(str).str.split('.').str[0]
    
    # 3. Объединяем exons и cds по Transcript_ID_base и Exon_Number
    logger.info(f"Объединяем exons ({len(exons_df)} строк) и cds ({len(cds_df)} строк)...")
    
    # Используем merge с how='left' чтобы сохранить все строки из exons
    merged = pd.merge(
        exons_df,
        cds_df.drop(columns=['Transcript_ID']),
        left_on=['Transcript_ID_base', 'Exon_Number'],
        right_on=['Transcript_ID_base', 'Exon_Number'],
        how='left',
        suffixes=('', '_cds')
    )
    
    logger.info(f"После слияния exons+cds: {len(merged)} строк")
    
    # 4. Добавляем информацию из transcript map
    if not transcript_map_df.empty and 'RefSeq_ID' in transcript_map_df.columns:
        logger.info("Добавляем информацию о транскриптах...")
        
        # Объединяем по базовой версии RefSeq ID
        merged = pd.merge(
            merged,
            transcript_map_df.drop(columns=['RefSeq_ID']),
            left_on='Transcript_ID_base',
            right_on='RefSeq_ID_base',
            how='left'
        )
        
        # Если Gene_x (из exons) пустой, заполняем из Gene_y (из transcript_map)
        if 'Gene_x' in merged.columns and 'Gene_y' in merged.columns:
            mask = (merged['Gene_x'].isna() | (merged['Gene_x'] == ''))
            merged.loc[mask, 'Gene_x'] = merged.loc[mask, 'Gene_y']
            # Переименовываем обратно в Gene
            merged['Gene'] = merged['Gene_x']
            merged = merged.drop(columns=['Gene_x', 'Gene_y'])
        elif 'Gene' in merged.columns:
            # Уже есть колонка Gene
            pass
        elif 'Gene_y' in merged.columns:
            merged['Gene'] = merged['Gene_y']
            merged = merged.drop(columns=['Gene_y'])
    
    # 5. Удаляем временные колонки
    temp_cols = ['Transcript_ID_base', 'RefSeq_ID_base']
    for col in temp_cols:
        if col in merged.columns:
            merged = merged.drop(columns=[col])
    
    # 6. Сортируем
    sort_cols = []
    if 'Gene' in merged.columns:
        sort_cols.append('Gene')
    if 'Transcript_ID' in merged.columns:
        sort_cols.append('Transcript_ID')
    if 'Exon_Number' in merged.columns:
        sort_cols.append('Exon_Number')
    
    if sort_cols:
        merged = merged.sort_values(sort_cols).reset_index(drop=True)
    
    logger.info(f"Итоговая таблица: {len(merged)} строк, {len(merged.columns)} колонок")
    logger.info(f"Колонки: {list(merged.columns)}")
    
    if 'Gene' in merged.columns:
        logger.info(f"Уникальных генов: {merged['Gene'].nunique()}")
        logger.info(f"Пример генов: {list(merged['Gene'].dropna().unique())[:10]}")
    
    return merged

def create_simple_exon_table(merged_df: pd.DataFrame, output_path: str):
    """Создание простой таблицы экзонов в формате для load_exons_table"""
    logger.info(f"Создаем простую таблицу экзонов: {output_path}")
    
    # Минимальный набор колонок для load_exons_table
    required_cols = ['Transcript_ID', 'Exon_Number', 'Genomic_start', 'Genomic_end']
    
    # Проверяем наличие колонок
    missing_cols = [col for col in required_cols if col not in merged_df.columns]
    if missing_cols:
        logger.error(f"Отсутствуют обязательные колонки: {missing_cols}")
        return None
    
    # Создаем простую таблицу
    simple_df = pd.DataFrame()
    simple_df['Transcript_ID'] = merged_df['Transcript_ID']
    simple_df['Exon'] = merged_df['Exon_Number'].astype(str)  # load_exons_table ожидает строку или число
    
    # Проверяем координаты
    if 'Genomic_start' in merged_df.columns:
        simple_df['Start'] = pd.to_numeric(merged_df['Genomic_start'], errors='coerce')
    else:
        logger.error("Нет колонки Genomic_start")
        return None
    
    if 'Genomic_end' in merged_df.columns:
        simple_df['End'] = pd.to_numeric(merged_df['Genomic_end'], errors='coerce')
    else:
        logger.error("Нет колонки Genomic_end")
        return None
    
    # Добавляем Gene если есть
    if 'Gene' in merged_df.columns:
        simple_df['Gene'] = merged_df['Gene']
    
    # Удаляем строки с пустыми координатами
    initial_len = len(simple_df)
    simple_df = simple_df.dropna(subset=['Start', 'End'])
    simple_df = simple_df[(simple_df['Start'] > 0) & (simple_df['End'] > simple_df['Start'])]
    
    if len(simple_df) < initial_len:
        logger.info(f"Удалено {initial_len - len(simple_df)} строк с некорректными координатами")
    
    # Сохраняем
    simple_df.to_csv(output_path, index=False)
    logger.info(f"Сохранено: {len(simple_df)} строк в {output_path}")
    
    return simple_df

def create_nmd_table(merged_df: pd.DataFrame, output_path: str):
    """Создание таблицы в формате для NMD-предсказания"""
    logger.info(f"Создаем NMD-таблицу: {output_path}")
    
    # Проверяем необходимые колонки
    required_cols = ['Transcript_ID', 'Exon_Number', 'Total_Exons', 'CDS_Exon_length', 
                    'CDS_length', 'Protein_length', 'Chrom', 'strand']
    
    missing_cols = [col for col in required_cols if col not in merged_df.columns]
    if missing_cols:
        logger.warning(f"Отсутствуют колонки для NMD таблицы: {missing_cols}")
        # Пробуем создать с тем, что есть
        available_cols = [col for col in required_cols if col in merged_df.columns]
    else:
        available_cols = required_cols
    
    # Добавляем координаты если есть
    if 'Genomic_start' in merged_df.columns and 'Genomic_end' in merged_df.columns:
        available_cols.extend(['Genomic_start', 'Genomic_end'])
    
    if 'Gene' in merged_df.columns:
        available_cols.append('Gene')
    
    # Создаем NMD таблицу
    nmd_df = merged_df[available_cols].copy()
    
    # Переименовываем для совместимости с load_nmd_table
    rename_map = {
        'Exon_Number': 'number_cds_exon',
        'Total_Exons': 'total_cds_exons',
        'CDS_Exon_length': 'CDS_Exon_length',
        'CDS_length': 'CDS_length',
        'Protein_length': 'Protein_length',
        'Chrom': 'chrom',
        'strand': 'strand',
        'Genomic_start': 'genomic_start',
        'Genomic_end': 'genomic_end'
    }
    
    # Применяем переименование только для существующих колонок
    actual_rename = {k: v for k, v in rename_map.items() if k in nmd_df.columns}
    if actual_rename:
        nmd_df = nmd_df.rename(columns=actual_rename)
    
    # Сохраняем
    nmd_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Сохранено: {len(nmd_df)} строк в {output_path}")
    
    return nmd_df

def main():
    parser = argparse.ArgumentParser(description='Объединение таблиц экзонов для NMD-предсказания')
    parser.add_argument('--exons-csv', required=True, help='Таблица exons.csv с геномными координатами')
    parser.add_argument('--cds-tsv', required=True, help='Таблица exon_df_cds.tsv с CDS информацией')
    parser.add_argument('--transcript-map', required=True, help='Таблица refseq2enst_mane.tsv с соответствием транскриптов')
    parser.add_argument('--output-dir', default='nmd_data', help='Выходная директория')
    parser.add_argument('--prefix', default='acmg', help='Префикс для выходных файлов')
    
    args = parser.parse_args()
    
    # Создаем выходную директорию
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Загружаем таблицы
    exons_df = load_exons_csv(args.exons_csv)
    cds_df = load_exon_cds_tsv(args.cds_tsv)
    transcript_map_df = load_transcript_map(args.transcript_map)
    
    # Проверяем данные перед объединением
    logger.info("\n" + "="*50)
    logger.info("ПРОВЕРКА ДАННЫХ:")
    logger.info(f"exons_df колонки: {list(exons_df.columns)}")
    logger.info(f"cds_df колонки: {list(cds_df.columns)[:10]}...")
    logger.info(f"transcript_map_df колонки: {list(transcript_map_df.columns)}")
    logger.info("="*50 + "\n")
    
    # Объединяем
    merged_df = merge_datasets_simple(exons_df, cds_df, transcript_map_df)
    
    if merged_df is None or merged_df.empty:
        logger.error("Не удалось объединить таблицы")
        return
    
    # Сохраняем объединенную таблицу
    merged_output = os.path.join(args.output_dir, f'{args.prefix}_merged_exons.tsv')
    merged_df.to_csv(merged_output, sep='\t', index=False)
    logger.info(f"Объединенная таблица сохранена: {merged_output}")
    
    # Создаем таблицу для NMD-предсказания
    nmd_output = os.path.join(args.output_dir, f'{args.prefix}_nmd_table.tsv')
    create_nmd_table(merged_df, nmd_output)
    
    # Создаем простую таблицу экзонов
    simple_output = os.path.join(args.output_dir, f'{args.prefix}_simple_exons.csv')
    create_simple_exon_table(merged_df, simple_output)
    
    # Выводим статистику
    logger.info("\n" + "="*50)
    logger.info("СТАТИСТИКА:")
    logger.info(f"Всего транскриптов: {merged_df['Transcript_ID'].nunique()}")
    logger.info(f"Всего экзонов: {len(merged_df)}")
    
    if 'Gene' in merged_df.columns:
        gene_count = merged_df['Gene'].nunique()
        logger.info(f"Всего генов: {gene_count}")
        if gene_count < 10:
            logger.info(f"Список генов: {sorted(merged_df['Gene'].dropna().unique())}")
    
    logger.info(f"\nВыходные файлы:")
    logger.info(f"1. Объединенная таблица: {merged_output}")
    logger.info(f"2. NMD таблица: {nmd_output}")
    logger.info(f"3. Простая таблица экзонов: {simple_output}")
    logger.info("="*50)

if __name__ == '__main__':
    main()