"""
VCF writer for MitoSAlt events
Outputs structural variants in VCF 4.3 format
"""

from __future__ import annotations
from typing import Optional, TextIO, List
from datetime import datetime
from pathlib import Path
import logging
import pandas as pd

logger = logging.getLogger(__name__)


class VCFWriter:
    """Write mitochondrial events in VCF format"""
    
    def __init__(
        self,
        reference_name: str = "chrM",
        sample_name: Optional[str] = None
    ) -> None:
        """
        Initialize VCF writer
        
        Args:
            reference_name: Chromosome/contig name (default: "chrM")
            sample_name: Sample identifier (extracted from data if None)
        """
        self.reference_name: str = reference_name
        self.sample_name: Optional[str] = sample_name
    
    def write(
        self,
        events_df: pd.DataFrame,
        output_file: str,
        genome_length: int = 16569
    ) -> None:
        """
        Write events to VCF file
        
        Args:
            events_df: DataFrame with events (must have group column from classify)
            output_file: Path to output VCF file
            genome_length: Length of mitochondrial genome (default: 16569)
        """
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Extract sample name if not provided
        if self.sample_name is None:
            if 'sample' in events_df.columns:
                self.sample_name = str(events_df['sample'].iloc[0])
            else:
                self.sample_name = "SAMPLE"
        
        # Open file and write
        with open(output_file, 'w') as f:
            # Write header
            self._write_header(f, genome_length)
            
            # Write events
            for _, event in events_df.iterrows():
                vcf_line = self._event_to_vcf(event)
                f.write(vcf_line + "\n")

        logger.info(f"VCF output saved to {output_file}")
        logger.info(f"Wrote {len(events_df)} events in VCF format")

    def _write_header(self, f: TextIO, genome_length: int) -> None:
        """
        Write VCF header
        
        Args:
            f: File handle to write to
            genome_length: Length of mitochondrial genome
        """
        # File format
        f.write("##fileformat=VCFv4.3\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=SaltShaker_v0.1.0\n")
        f.write("##reference=mitochondrial_genome\n")
        
        # Contig
        f.write(f"##contig=<ID={self.reference_name},length={genome_length}>\n")
        
        # INFO fields
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT">\n')
        f.write('##INFO=<ID=HF,Number=1,Type=Float,Description="Heteroplasmy fraction (0-1)">\n')
        f.write('##INFO=<ID=GROUP,Number=1,Type=String,Description="Spatial group identifier">\n')
        f.write('##INFO=<ID=CLUSTER,Number=1,Type=String,Description="Cluster identifier from MitoSAlt">\n')
        f.write('##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="Variant crosses D-loop region">\n')
        f.write('##INFO=<ID=BLCROSS,Number=0,Type=Flag,Description="Variant crosses blacklist region">\n')
        
        # FORMAT fields
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (ALT,REF)">\n')
        
        # ALT descriptions
        f.write('##ALT=<ID=DEL,Description="Deletion">\n')
        f.write('##ALT=<ID=DUP,Description="Duplication">\n')
        
        # Column header
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}\n")
    
    def _event_to_vcf(self, event: pd.Series) -> str:
        """
        Convert single event to VCF line
        
        Args:
            event: Series with event data
            
        Returns:
            Formatted VCF line string
        """
        # Basic fields
        chrom: str = self.reference_name
        pos: int = int(event['del_start_median'])
        id_field: str = "."
        ref: str = "N"  # Symbolic
        
        # ALT based on event type
        svtype: str = "DEL" if event['final_event'] == 'del' else "DUP"
        alt: str = f"<{svtype}>"
        
        qual: str = "."
        filter_field: str = "PASS"
        
        # INFO field
        info_parts: List[str] = [
            f"SVTYPE={svtype}",
            f"END={int(event['del_end_median'])}",
            f"SVLEN={int(event['delsize'])}",
            f"HF={event['perc']/100:.4f}",
            f"GROUP={event['group']}",
            f"CLUSTER={event.get('cluster', 'NA')}"
        ]
        
        # Add flags
        if event.get('dloop') == 'yes':
            info_parts.append("DLOOP")
        bl_cross = event.get('blacklist_crossing', False)
        if bl_cross is True or bl_cross == 'yes' or str(bl_cross).lower() == 'true':
            info_parts.append("BLCROSS")
        
        info: str = ";".join(info_parts)
        
        # FORMAT field
        format_field: str = "GT:AD"
        
        # Sample field
        gt: str = "0/1"
        alt_reads: int = int(event.get('nread', 0))
        ref_reads: int = int(event.get('tread', 0)) - alt_reads
        sample_field: str = f"{gt}:{alt_reads},{ref_reads}"
        
        # Combine
        vcf_line: str = "\t".join([
            chrom, str(pos), id_field, ref, alt, qual, filter_field,
            info, format_field, sample_field
        ])
        
        return vcf_line


def write_vcf(
    events_df: pd.DataFrame,
    output_file: str,
    genome_length: int,
    reference_name: str = "chrM",
    sample_name: Optional[str] = None
) -> None:
    """
    Convenience function to write VCF
    
    Args:
        events_df: DataFrame with mitochondrial events (must have group column)
        output_file: Output VCF file path
        genome_length: Mitochondrial genome length
        reference_name: Chromosome/contig name (default: "chrM")
        sample_name: Sample identifier (default: None, extracted from data)
    """
    writer = VCFWriter(reference_name=reference_name, sample_name=sample_name)
    writer.write(events_df, output_file, genome_length)