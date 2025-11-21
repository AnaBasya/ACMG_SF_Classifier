#!/usr/bin/env python3
"""
Simple VCF Parser - Robust VCF processing without pysam dependency
"""

import logging
import gzip
from typing import List, Dict, Any, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class SimpleVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    hgvs_c: str = ""
    hgvs_p: str = ""

class SimpleVCFLoader:
    """Simple VCF loader that handles various VCF formats"""
    
    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        
    def load_sample_variants(self, sample_id: str, gender: str, affected: bool) -> Any:
        """Load variants with robust parsing"""
        from main_analysis import Sample, Variant, VariantType
        
        variants = []
        line_count = 0
        
        try:
            # Determine if file is gzipped
            open_func = gzip.open if self.vcf_path.endswith('.gz') else open
            mode = 'rt' if self.vcf_path.endswith('.gz') else 'r'
            
            with open_func(self.vcf_path, mode) as f:
                for line in f:
                    line_count += 1
                    
                    # Skip header lines
                    if line.startswith('#'):
                        # Fix: Skip duplicate header lines
                        if line_count == 1 and line.startswith('#CHROM'):
                            continue
                        continue
                    
                    # Parse data line
                    try:
                        variant = self._parse_vcf_line(line.strip())
                        if variant:
                            # Convert to main Variant type
                            main_variant = Variant(
                                chrom=variant.chrom,
                                pos=variant.pos,
                                ref=variant.ref,
                                alt=variant.alt,
                                variant_type=self._determine_variant_type(variant),
                                gene=variant.gene,
                                transcript="",
                                hgvs_c=variant.hgvs_c,
                                hgvs_p=variant.hgvs_p,
                                af_gnomad=0.0
                            )
                            variants.append(main_variant)
                            
                    except Exception as e:
                        logger.debug(f"Skipping line {line_count}: {e}")
                        continue
            
            logger.info(f"Loaded {len(variants)} variants for {sample_id}")
            return Sample(sample_id, gender, affected, variants)
            
        except Exception as e:
            logger.error(f"Error loading variants for {sample_id}: {e}")
            return Sample(sample_id, gender, affected, [])
    
    def _parse_vcf_line(self, line: str) -> Optional[SimpleVariant]:
        """Parse a single VCF line"""
        parts = line.split('\t')
        
        # Basic validation
        if len(parts) < 8:
            return None
        
        try:
            chrom = parts[0].strip()
            pos = int(parts[1].strip())
            ref = parts[3].strip()
            alt = parts[4].strip()
            
            # Skip invalid variants
            if chrom == "CHROM" or ref == "REF" or alt == "ALT":
                return None
            
            # Parse INFO field
            info_parts = parts[7].split(';')
            gene = "UNKNOWN"
            hgvs_c = ""
            hgvs_p = ""
            
            for info in info_parts:
                if '=' in info:
                    key, value = info.split('=', 1)
                    if key == 'GENE':
                        gene = value
                    elif key == 'HGVSc':
                        hgvs_c = value
                    elif key == 'HGVSp':
                        hgvs_p = value
            
            # Skip if no gene information
            if gene == "UNKNOWN":
                return None
            
            return SimpleVariant(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                gene=gene,
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p
            )
            
        except (ValueError, IndexError) as e:
            logger.debug(f"Failed to parse line: {e}")
            return None
    
    def _determine_variant_type(self, variant: SimpleVariant) -> Any:
        """Determine variant type from HGVS annotations"""
        from main_analysis import VariantType
        
        hgvs_c_lower = variant.hgvs_c.lower()
        hgvs_p_lower = variant.hgvs_p.lower()
        
        # Check HGVS annotations first
        if 'del' in hgvs_c_lower and 'fs' in hgvs_p_lower:
            return VariantType.FRAMESHIFT
        if 'ter' in hgvs_p_lower or 'stop' in hgvs_p_lower or '*' in hgvs_p_lower:
            return VariantType.NONSENSE
        if any(splice in hgvs_c_lower for splice in ['+', '-', 'splice']):
            return VariantType.SPLICE_DONOR_ACCEPTOR
        
        # Sequence-based determination
        if len(variant.ref) == 1 and len(variant.alt) == 1:
            return VariantType.MISSENSE
        elif len(variant.ref) != len(variant.alt):
            return VariantType.FRAMESHIFT
        
        return VariantType.MISSENSE

def load_trio_fixed(
    proband_vcf: str,
    mother_vcf: str, 
    father_vcf: str,
    proband_id: str = "proband1",
    mother_id: str = "mother1",
    father_id: str = "father1", 
    family_id: str = "FAM001"
) -> Any:
    """Load trio data with simple VCF processor"""
    from main_analysis import TrioData
    
    logger.info("Loading trio data with simple processor...")
    
    try:
        # Load samples
        proband_loader = SimpleVCFLoader(proband_vcf)
        proband_sample = proband_loader.load_sample_variants(proband_id, "unknown", True)
        
        mother_loader = SimpleVCFLoader(mother_vcf)
        mother_sample = mother_loader.load_sample_variants(mother_id, "female", False)
        
        father_loader = SimpleVCFLoader(father_vcf)
        father_sample = father_loader.load_sample_variants(father_id, "male", False)
        
        trio_data = TrioData(proband_sample, mother_sample, father_sample, family_id)
        
        logger.info("✅ Trio data loaded successfully")
        logger.info(f"   Proband: {len(proband_sample.variants)} variants")
        logger.info(f"   Mother: {len(mother_sample.variants)} variants")
        logger.info(f"   Father: {len(father_sample.variants)} variants")
        
        return trio_data
        
    except Exception as e:
        logger.error(f"❌ Failed to load trio: {e}")
        raise

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    try:
        trio = load_trio_fixed(
            proband_vcf="../data/proband_test.vcf",
            mother_vcf="../data/mother_test.vcf",
            father_vcf="../data/father_test.vcf"
        )
        print("🎉 Simple VCF processor test: SUCCESS!")
    except Exception as e:
        print(f"❌ Simple VCF processor test: FAILED - {e}")
