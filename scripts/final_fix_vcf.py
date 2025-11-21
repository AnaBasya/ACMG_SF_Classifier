#!/usr/bin/env python3
"""
Финальное исправление VCF файлов - удаляем дублирующий заголовок
"""

import sys
from pathlib import Path
import shutil

def final_fix_vcf(input_file: Path, output_file: Path) -> bool:
    """Финальное исправление - удаляем дублирующий заголовок"""
    try:
        print(f"Финальное исправление: {input_file}")
        
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            header_written = False
            duplicate_found = False
            
            for line_num, line in enumerate(fin, 1):
                line = line.strip()
                
                # Пропускаем дублирующую строку заголовка
                if line.startswith('CHROM\tPOS') or line.startswith('CHROM  POS'):
                    if not duplicate_found:
                        print(f"  Обнаружен дублирующий заголовок в строке {line_num} - пропускаем")
                        duplicate_found = True
                    continue
                
                # Все остальные строки пишем как есть
                fout.write(line + '\n')
                header_written = True
        
        if duplicate_found:
            print(f"  ✓ Дублирующий заголовок удален")
        else:
            print(f"  ✓ Файл уже в правильном формате")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Ошибка: {e}")
        return False

def main():
    data_dir = Path("data")
    vcf_files = list(data_dir.glob("*.vcf"))
    
    if not vcf_files:
        print("VCF файлы не найдены")
        return 1
    
    print("Финальное исправление VCF файлов...")
    print("=" * 50)
    
    for vcf_file in vcf_files:
        print(f"Обработка: {vcf_file.name}")
        
        # Создаем временный файл
        temp_file = vcf_file.with_suffix('.temp.vcf')
        
        if final_fix_vcf(vcf_file, temp_file):
            # Заменяем оригинальный файл исправленным
            backup_file = vcf_file.with_suffix('.vcf.backup2')
            shutil.copy2(vcf_file, backup_file)
            shutil.move(temp_file, vcf_file)
            
            print(f"  ✓ Файл исправлен")
            print(f"  Создан backup: {backup_file}")
            
            # Покажем первые строки исправленного файла
            print("  Первые строки исправленного файла:")
            with open(vcf_file, 'r') as f:
                for i, line in enumerate(f):
                    if i < 3:
                        print(f"    {line.strip()}")
                    if i >= 5:
                        break
        else:
            print(f"  ✗ Не удалось исправить файл")
        
        print()
    
    print("=" * 50)
    print("Исправление завершено")

if __name__ == "__main__":
    sys.exit(main())
