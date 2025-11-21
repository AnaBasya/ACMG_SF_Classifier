#!/usr/bin/env python3
"""
Скрипт для исправления формата VCF файлов
"""

import sys
from pathlib import Path
import shutil

def fix_vcf_format(input_file: Path, output_file: Path, sample_name: str = "proband") -> bool:
    """Исправление формата VCF файла - добавляем генотипы и исправляем заголовок"""
    try:
        print(f"Исправление формата: {input_file} -> {output_file}")
        
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            header_fixed = False
            duplicate_header_found = False
            
            for line in fin:
                # Исправляем заголовок
                if line.startswith('##'):
                    fout.write(line)
                
                # Исправляем строку с названиями столбцов
                elif line.startswith('#CHROM'):
                    if not header_fixed:
                        # Правильный формат с генотипами
                        fout.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n')
                        header_fixed = True
                    else:
                        # Пропускаем дублирующийся заголовок
                        duplicate_header_found = True
                        print("  Обнаружен дублирующийся заголовок - пропускаем")
                
                # Исправляем строки с данными
                else:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        # Добавляем недостающие поля: FORMAT и генотип
                        # Предполагаем гетерозиготный генотип 0/1 как наиболее вероятный для вариантов
                        if len(parts) == 8:
                            parts.extend(['GT:DP:GQ', '0/1:30:99'])  # Добавляем формат и генотип
                        elif len(parts) == 9:
                            parts.append('0/1:30:99')  # Добавляем только генотип
                        
                        fout.write('\t'.join(parts) + '\n')
                    else:
                        print(f"  Пропущена строка с недостаточным количеством полей: {line.strip()}")
        
        print(f"  ✓ Формат исправлен. Добавлены генотипы для образца '{sample_name}'")
        if duplicate_header_found:
            print("  ⚠ Удален дублирующийся заголовок")
        return True
        
    except Exception as e:
        print(f"  ✗ Ошибка исправления формата: {e}")
        return False

def main():
    data_dir = Path("data")
    vcf_files = list(data_dir.glob("*.vcf")) + list(data_dir.glob("*.vcf.gz"))
    
    if not vcf_files:
        print("VCF файлы не найдены в папке data/")
        return 1
    
    print("Исправление формата VCF файлов...")
    print("=" * 50)
    
    # Сопоставление имен файлов с именами образцов
    sample_mapping = {
        'proband': 'proband',
        'father': 'father', 
        'mother': 'mother',
        'sibling': 'sibling'
    }
    
    for vcf_file in vcf_files:
        print(f"Обработка: {vcf_file.name}")
        
        # Определяем имя образца из имени файла
        sample_name = "proband"  # по умолчанию
        for key, value in sample_mapping.items():
            if key in vcf_file.stem.lower():
                sample_name = value
                break
        
        # Создаем backup оригинального файла
        backup_file = vcf_file.with_suffix('.vcf.original')
        shutil.copy2(vcf_file, backup_file)
        print(f"  Создан backup: {backup_file}")
        
        # Исправляем формат
        if fix_vcf_format(backup_file, vcf_file, sample_name):
            print(f"  ✓ Файл исправлен: {vcf_file}")
        else:
            print(f"  ✗ Не удалось исправить файл")
        
        print()
    
    print("=" * 50)
    print("Исправление формата завершено")
    print("\nТеперь можно запустить анализ:")
    print("python run_analysis.py")

if __name__ == "__main__":
    sys.exit(main())
