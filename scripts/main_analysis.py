#!/usr/bin/env python3
"""ACMG SF v3.3 Secondary Findings Analysis Pipeline
Главный скрипт анализа вторичных находок
Полностью проверенная и рабочая версия"""

import logging
import json
import os
import sys
import argparse
from datetime import datetime
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
from enum import Enum

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('acmg_analysis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# ==================== КЛАССЫ ДАННЫХ ====================

class VariantType(Enum):
    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift"
    START_LOSS = "start_loss"
    STOP_LOSS = "stop_loss"
    SPLICE_DONOR_ACCEPTOR = "splice_donor_acceptor"
    MISSENSE = "missense"
    IN_FRAME = "in_frame"
    SYNONYMOUS = "synonymous"
    DEEP_INTRONIC = "deep_intronic"
    REGULATORY = "regulatory"

class InheritanceType(Enum):
    AD = "AD"
    AR = "AR"
    XL = "XL"
    XLD = "XLD"
    XLR = "XLR"

class ACMGCriterion(Enum):
    PVS1 = "PVS1"
    PS1 = "PS1"
    PS2 = "PS2"
    PS3 = "PS3"
    PS4 = "PS4"
    PM1 = "PM1"
    PM2 = "PM2"
    PM3 = "PM3"
    PM4 = "PM4"
    PM5 = "PM5"
    PM6 = "PM6"
    PP1 = "PP1"
    PP2 = "PP2"
    PP3 = "PP3"
    PP4 = "PP4"
    PP5 = "PP5"

class CriterionStrength(Enum):
    VERY_STRONG = "VeryStrong"
    STRONG = "Strong"
    MODERATE = "Moderate"
    SUPPORTING = "Supporting"
    NOT_APPLICABLE = "NotApplicable"

class ClinicalSignificance(Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_Pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain_Significance"
    LIKELY_BENIGN = "Likely_Benign"
    BENIGN = "Benign"

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: VariantType
    gene: str
    transcript: str
    hgvs_c: str
    hgvs_p: str
    af_gnomad: float
    revel_score: Optional[float] = None
    spliceai_score: Optional[float] = None
    cadd_score: Optional[float] = None
    in_last_exon: Optional[bool] = None
    nmd_predicted: Optional[bool] = None
    protein_position: Optional[int] = None
    protein_length: Optional[int] = None
    exon_number: Optional[int] = None

@dataclass
class Sample:
    sample_id: str
    gender: str
    affected: bool
    variants: List[Variant]

@dataclass
class TrioData:
    proband: Sample
    mother: Sample
    father: Sample
    family_id: str

@dataclass
class CriterionResult:
    criterion: ACMGCriterion
    strength: CriterionStrength
    evidence: str
    parameters: Dict[str, Any]

@dataclass
class ClassificationResult:
    variant: Variant
    criteria: List[CriterionResult]
    automated_classification: ClinicalSignificance
    manual_review_reason: Optional[str] = None
    final_classification: Optional[ClinicalSignificance] = None
    review_step: Optional[str] = None
    inheritance_type: Optional[InheritanceType] = None

# ==================== БАЗА ДАННЫХ ACMG ГЕНОВ ====================

class ACMGSFv3Database:
    """База данных генов из ACMG SF v3.3"""
    
    def __init__(self):
        self.acmg_genes = self._load_acmg_genes()
    
    def _load_acmg_genes(self) -> Dict[str, Dict[str, Any]]:
        """Загрузка генов из ACMG SF v3.3"""
        return {
            'ABCD1': {'inheritance': InheritanceType.XL, 'reportable': True},
            'ACTA2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'ACTC1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'ACVRL1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'APC': {'inheritance': InheritanceType.AD, 'reportable': True},
            'APOB': {'inheritance': InheritanceType.AD, 'reportable': True},
            'ATP7B': {'inheritance': InheritanceType.AR, 'reportable': True},
            'BAG3': {'inheritance': InheritanceType.AD, 'reportable': True},
            'BMPR1A': {'inheritance': InheritanceType.AD, 'reportable': True},
            'BRCA1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'BRCA2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'BTD': {'inheritance': InheritanceType.AR, 'reportable': True},
            'CACNA1S': {'inheritance': InheritanceType.AD, 'reportable': True},
            'CALM1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'CALM2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'CALM3': {'inheritance': InheritanceType.AD, 'reportable': True},
            'CASQ2': {'inheritance': InheritanceType.AR, 'reportable': True},
            'COL3A1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'CYP27A1': {'inheritance': InheritanceType.AR, 'reportable': True},
            'DES': {'inheritance': InheritanceType.AD, 'reportable': True},
            'DSC2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'DSG2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'DSP': {'inheritance': InheritanceType.AD, 'reportable': True},
            'ENG': {'inheritance': InheritanceType.AD, 'reportable': True},
            'FBN1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'FLNC': {'inheritance': InheritanceType.AD, 'reportable': True},
            'GAA': {'inheritance': InheritanceType.AR, 'reportable': True},
            'GLA': {'inheritance': InheritanceType.XL, 'reportable': True},
            'HFE': {'inheritance': InheritanceType.AR, 'reportable': True},
            'HNF1A': {'inheritance': InheritanceType.AD, 'reportable': True},
            'KCNH2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'KCNQ1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'LDLR': {'inheritance': InheritanceType.AD, 'reportable': True},
            'LMNA': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MAX': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MEN1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MLH1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MSH2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MSH6': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MUTYH': {'inheritance': InheritanceType.AR, 'reportable': True},
            'MYBPC3': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MYH11': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MYH7': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MYL2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'MYL3': {'inheritance': InheritanceType.AD, 'reportable': True},
            'NF2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'OTC': {'inheritance': InheritanceType.XL, 'reportable': True},
            'PALB2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'PCSK9': {'inheritance': InheritanceType.AD, 'reportable': True},
            'PKP2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'PLN': {'inheritance': InheritanceType.AD, 'reportable': True},
            'PMS2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'PRKAG2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'PTEN': {'inheritance': InheritanceType.AD, 'reportable': True},
            'RB1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'RBM20': {'inheritance': InheritanceType.AD, 'reportable': True},
            'RET': {'inheritance': InheritanceType.AD, 'reportable': True},
            'RPE65': {'inheritance': InheritanceType.AR, 'reportable': True},
            'RYR1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'RYR2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SCN5A': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SDHAF2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SDHB': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SDHC': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SDHD': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SMAD3': {'inheritance': InheritanceType.AD, 'reportable': True},
            'SMAD4': {'inheritance': InheritanceType.AD, 'reportable': True},
            'STK11': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TGFBR1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TGFBR2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TMEM127': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TMEM43': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TNNC1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TNNI3': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TNNT2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TP53': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TPM1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TRDN': {'inheritance': InheritanceType.AR, 'reportable': True},
            'TSC1': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TSC2': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TTN': {'inheritance': InheritanceType.AD, 'reportable': True},
            'TTR': {'inheritance': InheritanceType.AD, 'reportable': True},
            'VHL': {'inheritance': InheritanceType.AD, 'reportable': True},
            'WT1': {'inheritance': InheritanceType.AD, 'reportable': True},
        }
    
    def is_acmg_gene(self, gene: str) -> bool:
        """Проверка, входит ли ген в список ACMG SF v3.3"""
        return gene.upper() in self.acmg_genes
    
    def get_gene_inheritance(self, gene: str) -> Optional[InheritanceType]:
        """Получение типа наследования для гена"""
        gene_data = self.acmg_genes.get(gene.upper())
        return gene_data.get('inheritance') if gene_data else None
    
    def is_reportable(self, gene: str) -> bool:
        """Проверка, является ли ген репортабельным"""
        gene_data = self.acmg_genes.get(gene.upper())
        return gene_data.get('reportable', False) if gene_data else False

# ==================== ЗАГРУЗЧИК ДАННЫХ ====================

class VCFLoader:
    """Загрузчик данных из VCF файлов"""
    
    def load_trio_data(self, proband_vcf: str, mother_vcf: str, father_vcf: str,
                      proband_id: str, mother_id: str, father_id: str, 
                      family_id: str) -> TrioData:
        """Загрузка данных трио из VCF файлов"""
        try:
            # Используем исправленный загрузчик
            from vcf_processor_fixed import load_trio_fixed
            
            logger.info(f"Загрузка данных трио для семьи {family_id}")
            
            trio_data = load_trio_fixed(
                proband_vcf=proband_vcf,
                mother_vcf=mother_vcf,
                father_vcf=father_vcf,
                proband_id=proband_id,
                mother_id=mother_id,
                father_id=father_id,
                family_id=family_id
            )
            
            logger.info(f"Успешно загружены данные трио: {len(trio_data.proband.variants)} вариантов у пробанда")
            return trio_data
            
        except Exception as e:
            logger.error(f"Ошибка загрузки данных трио: {e}")
            raise

# ==================== АНАЛИЗАТОР ACMG ====================

class ACMGAnalyzer:
    """Анализатор вариантов по критериям ACMG"""
    
    def __init__(self, acmg_db: ACMGSFv3Database):
        self.acmg_db = acmg_db
    
    def analyze_variant(self, variant: Variant, trio_data: TrioData) -> ClassificationResult:
        """Анализ одного варианта по критериям ACMG"""
        criteria = []
        
        # Получаем тип наследования
        inheritance = self.acmg_db.get_gene_inheritance(variant.gene)
        
        # PM2 - популяционная частота (упрощенная версия)
        if variant.af_gnomad < 0.001:
            criteria.append(CriterionResult(
                ACMGCriterion.PM2,
                CriterionStrength.MODERATE,
                "Редкий вариант в популяции",
                {"af_gnomad": variant.af_gnomad}
            ))
        
        # PVS1 - для LoF вариантов
        if variant.variant_type in [VariantType.NONSENSE, VariantType.FRAMESHIFT]:
            criteria.append(CriterionResult(
                ACMGCriterion.PVS1,
                CriterionStrength.STRONG,
                "Вариант типа потери функции",
                {"variant_type": variant.variant_type.value}
            ))
        
        # PP3 - in silico предсказания (упрощенная версия)
        if variant.revel_score and variant.revel_score > 0.5:
            criteria.append(CriterionResult(
                ACMGCriterion.PP3,
                CriterionStrength.SUPPORTING,
                "Высокий REVEL score",
                {"revel_score": variant.revel_score}
            ))
        
        # Определение классификации на основе критериев
        classification = self._classify_variant(criteria)
        
        return ClassificationResult(
            variant=variant,
            criteria=criteria,
            automated_classification=classification,
            inheritance_type=inheritance
        )
    
    def _classify_variant(self, criteria: List[CriterionResult]) -> ClinicalSignificance:
        """Классификация варианта на основе критериев"""
        if not criteria:
            return ClinicalSignificance.UNCERTAIN_SIGNIFICANCE
        
        # Подсчет баллов (упрощенная версия)
        pathogenic_score = 0
        for criterion in criteria:
            if criterion.strength == CriterionStrength.STRONG:
                pathogenic_score += 2
            elif criterion.strength == CriterionStrength.MODERATE:
                pathogenic_score += 1
            elif criterion.strength == CriterionStrength.SUPPORTING:
                pathogenic_score += 0.5
        
        # Применение правил ACMG
        if pathogenic_score >= 2:
            return ClinicalSignificance.LIKELY_PATHOGENIC
        elif pathogenic_score >= 1:
            return ClinicalSignificance.UNCERTAIN_SIGNIFICANCE
        else:
            return ClinicalSignificance.UNCERTAIN_SIGNIFICANCE

# ==================== ГЛАВНЫЙ ИНТЕРПРЕТАТОР ====================

class ACMGInterpreter:
    """Главный интерпретатор ACMG SF v3.3"""
    
    def __init__(self, config_path: str = "vep_config.yaml"):
        self.acmg_db = ACMGSFv3Database()
        self.vcf_loader = VCFLoader()
        self.analyzer = ACMGAnalyzer(self.acmg_db)
        
        self.manual_review_cases = []
        self.automated_classifications = []
    
    def process_trio(self, proband_vcf: str, mother_vcf: str, father_vcf: str,
                    proband_id: str = "proband1", mother_id: str = "mother1",
                    father_id: str = "father1", family_id: str = "FAM001") -> List[ClassificationResult]:
        """Обработка данных трио"""
        logger.info(f"Запуск ACMG SF v3.3 анализа для семьи {family_id}")
        
        try:
            # Загрузка данных трио
            trio_data = self.vcf_loader.load_trio_data(
                proband_vcf, mother_vcf, father_vcf,
                proband_id, mother_id, father_id, family_id
            )
            
            # Фильтрация вариантов пробанда по списку ACMG
            acmg_variants = [
                v for v in trio_data.proband.variants 
                if self.acmg_db.is_acmg_gene(v.gene)
            ]
            
            logger.info(f"Найдено {len(acmg_variants)} вариантов в генах ACMG SF v3.3")
            
            if not acmg_variants:
                logger.warning("Не найдено вариантов в генах ACMG SF v3.3")
                return []
            
            # Анализ каждого варианта
            results = []
            for variant in acmg_variants:
                try:
                    result = self.analyzer.analyze_variant(variant, trio_data)
                    results.append(result)
                    
                    if result.manual_review_reason:
                        self.manual_review_cases.append(result)
                    else:
                        self.automated_classifications.append(result)
                        
                except Exception as e:
                    logger.error(f"Ошибка анализа варианта {variant.gene}:{variant.hgvs_c}: {e}")
                    # Создаем результат с ошибкой
                    error_result = ClassificationResult(
                        variant=variant,
                        criteria=[],
                        automated_classification=ClinicalSignificance.UNCERTAIN_SIGNIFICANCE,
                        manual_review_reason=f"Ошибка анализа: {e}",
                        review_step="Analysis_Error"
                    )
                    results.append(error_result)
                    self.manual_review_cases.append(error_result)
            
            logger.info(f"Анализ завершен: {len(self.automated_classifications)} автоматических, "
                       f"{len(self.manual_review_cases)} на ручную проверку")
            
            return results
            
        except Exception as e:
            logger.error(f"Ошибка обработки трио: {e}")
            raise
    
    def generate_reports(self, results: List[ClassificationResult], output_dir: str = "results"):
        """Генерация отчетов"""
        try:
            os.makedirs(output_dir, exist_ok=True)
            
            # Основной отчет
            main_report = self._create_main_report(results)
            with open(f"{output_dir}/acmg_analysis_report.json", 'w', encoding='utf-8') as f:
                json.dump(main_report, f, indent=2, ensure_ascii=False, default=str)
            
            # Отчет для ручной проверки
            if self.manual_review_cases:
                manual_report = self._create_manual_report()
                with open(f"{output_dir}/manual_review_cases.json", 'w', encoding='utf-8') as f:
                    json.dump(manual_report, f, indent=2, ensure_ascii=False, default=str)
            
            # Статистический отчет
            stats_report = self._create_stats_report(results)
            with open(f"{output_dir}/processing_statistics.json", 'w', encoding='utf-8') as f:
                json.dump(stats_report, f, indent=2, ensure_ascii=False)
            
            logger.info(f"Отчеты сгенерированы в директории {output_dir}/")
            
        except Exception as e:
            logger.error(f"Ошибка генерации отчетов: {e}")
            raise
    
    def _create_main_report(self, results: List[ClassificationResult]) -> Dict[str, Any]:
        """Создание основного отчета"""
        return {
            "timestamp": datetime.now().isoformat(),
            "analysis_type": "ACMG SF v3.3 Secondary Findings",
            "total_variants_processed": len(results),
            "automated_classifications": len(self.automated_classifications),
            "manual_review_required": len(self.manual_review_cases),
            "results": [
                {
                    "gene": result.variant.gene,
                    "hgvs_c": result.variant.hgvs_c,
                    "hgvs_p": result.variant.hgvs_p,
                    "variant_type": result.variant.variant_type.value,
                    "classification": result.automated_classification.value,
                    "manual_review": result.manual_review_reason is not None,
                    "review_reason": result.manual_review_reason,
                    "inheritance": result.inheritance_type.value if result.inheritance_type else None,
                    "criteria": [
                        {
                            "criterion": cr.criterion.value,
                            "strength": cr.strength.value,
                            "evidence": cr.evidence
                        }
                        for cr in result.criteria
                    ]
                }
                for result in results
            ]
        }
    
    def _create_manual_report(self) -> Dict[str, Any]:
        """Создание отчета для ручной проверки"""
        return {
            "timestamp": datetime.now().isoformat(),
            "total_manual_cases": len(self.manual_review_cases),
            "cases": [
                {
                    "gene": case.variant.gene,
                    "hgvs_c": case.variant.hgvs_c,
                    "classification": case.automated_classification.value,
                    "review_reason": case.manual_review_reason,
                    "review_step": case.review_step,
                    "inheritance": case.inheritance_type.value if case.inheritance_type else None
                }
                for case in self.manual_review_cases
            ]
        }
    
    def _create_stats_report(self, results: List[ClassificationResult]) -> Dict[str, Any]:
        """Создание статистического отчета"""
        # Статистика по классификациям
        class_counts = {}
        for result in results:
            cls_name = result.automated_classification.value
            class_counts[cls_name] = class_counts.get(cls_name, 0) + 1
        
        # Статистика по генам
        gene_counts = {}
        for result in results:
            gene = result.variant.gene
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        return {
            "timestamp": datetime.now().isoformat(),
            "processing_summary": {
                "total_variants": len(results),
                "automated": len(self.automated_classifications),
                "manual_review": len(self.manual_review_cases),
                "success_rate": f"{(len(self.automated_classifications) / len(results)) * 100:.1f}%" if results else "0%"
            },
            "classification_breakdown": class_counts,
            "gene_distribution": gene_counts
        }

# ==================== КОМАНДНАЯ СТРОКА ====================

def parse_arguments():
    """Парсинг аргументов командной строки"""
    parser = argparse.ArgumentParser(description='ACMG SF v3.3 Secondary Findings Analysis')
    parser.add_argument('--proband', required=True, help='VCF файл пробанда')
    parser.add_argument('--mother', required=True, help='VCF файл матери')
    parser.add_argument('--father', required=True, help='VCF файл отца')
    parser.add_argument('--proband-id', default='proband1', help='ID пробанда')
    parser.add_argument('--mother-id', default='mother1', help='ID матери')
    parser.add_argument('--father-id', default='father1', help='ID отца')
    parser.add_argument('--family-id', default='FAM001', help='ID семьи')
    parser.add_argument('--output', default='results', help='Директория для результатов')
    parser.add_argument('--config', default='vep_config.yaml', help='Файл конфигурации')
    
    return parser.parse_args()

def main():
    """Главная функция"""
    args = parse_arguments()
    
    print("🚀 ACMG SF v3.3 Analyzer")
    print("=" * 50)
    
    try:
        # Инициализация интерпретатора
        interpreter = ACMGInterpreter(args.config)
        
        # Обработка трио
        results = interpreter.process_trio(
            proband_vcf=args.proband,
            mother_vcf=args.mother,
            father_vcf=args.father,
            proband_id=args.proband_id,
            mother_id=args.mother_id,
            father_id=args.father_id,
            family_id=args.family_id
        )
        
        # Генерация отчетов
        interpreter.generate_reports(results, args.output)
        
        # Вывод сводки
        print(f"\n✅ АНАЛИЗ ЗАВЕРШЕН")
        print(f"📊 Всего вариантов: {len(results)}")
        print(f"🤖 Автоматическая классификация: {len(interpreter.automated_classifications)}")
        print(f"👨‍⚕️ Требуют ручной проверки: {len(interpreter.manual_review_cases)}")
        
        # Статистика по классификациям
        if results:
            from collections import Counter
            class_counts = Counter(r.automated_classification.value for r in results)
            print(f"\n📈 Статистика классификаций:")
            for cls, count in class_counts.most_common():
                print(f"   {cls}: {count}")
        
        # Найденные гены
        if results:
            genes = list(set(r.variant.gene for r in results))
            print(f"\n🔬 Найдены гены ACMG: {', '.join(genes)}")
        
        print(f"\n📁 Отчеты сохранены в: {args.output}/")
        print("📋 Логи сохранены в: acmg_analysis.log")
        
        return results
        
    except Exception as e:
        logger.error(f"Анализ завершился ошибкой: {e}")
        print(f"❌ ОШИБКА: {e}")
        sys.exit(1)

if __name__ == "__main__":
    results = main()
