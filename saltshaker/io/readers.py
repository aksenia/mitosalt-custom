"""
I/O Readers

Handles reading of various input files.
"""

from __future__ import annotations
from typing import List, Tuple, Optional
from importlib.resources import path
import pandas as pd
from pathlib import Path
import logging

from ..types import BlacklistRegion, GeneAnnotation

logger = logging.getLogger(__name__)


class BlacklistReader:
    """Reads and parses blacklist BED files"""
    
    @staticmethod
    def load_blacklist_regions(blacklist_file: str) -> List[BlacklistRegion]:
        """
        Robustly load blacklist regions from BED file
        
        Tries multiple parsing strategies to handle various BED formats.
        
        Args:
            blacklist_file: Path to BED file with regions to exclude
            
        Returns:
            List of dicts with 'chr', 'start', 'end' keys
        """
        blacklist_regions: List[BlacklistRegion] = []
        if not blacklist_file or not Path(blacklist_file).exists():
            return blacklist_regions
        
        try:
            blacklist: Optional[pd.DataFrame] = None
            # Try different separators
            for sep in ['\t', None, ' ', '\\s+']:
                try:
                    if sep == '\\s+':
                        blacklist = pd.read_csv(blacklist_file, sep=sep, header=None, engine='python')
                    elif sep is None:
                        blacklist = pd.read_csv(blacklist_file, sep=sep, header=None, engine='python')
                    else:
                        blacklist = pd.read_csv(blacklist_file, sep=sep, header=None)
                    
                    if blacklist.shape[1] >= 3:
                        break
                    else:
                        blacklist = None
                except Exception:
                    continue
            
            # Fallback: manual parsing
            if blacklist is None or blacklist.shape[1] < 3:
                with open(blacklist_file, 'r') as f:
                    lines = f.readlines()
                parsed_lines = []
                for line in lines:
                    line = line.strip()
                    if line:
                        parts = line.split()
                        if len(parts) >= 3:
                            parsed_lines.append([parts[0], parts[1], parts[2]])
                if parsed_lines:
                    blacklist = pd.DataFrame(parsed_lines, columns=['chr', 'start', 'end'])
            
            if blacklist is not None:
                blacklist = blacklist.iloc[:, :3].copy()
                blacklist.columns = ['chr', 'start', 'end']
                blacklist['start'] = pd.to_numeric(blacklist['start'], errors='coerce')
                blacklist['end'] = pd.to_numeric(blacklist['end'], errors='coerce')
                blacklist = blacklist.dropna()
                blacklist_regions = blacklist.to_dict('records')  # type: ignore
        
        except Exception as e:
            logger.warning(f"Could not load blacklist file: {e}")
        
        return blacklist_regions
    

class IntermediateReader:
    """Reads events in internal format"""
    
    @staticmethod
    def read(filepath: str) -> Tuple[pd.DataFrame, int]:
        """
        Read events with metadata
        
        Args:
            filepath: Path to intermediate TSV file
        
        Returns:
            Tuple of (events_df, genome_length)
        """
        genome_length: Optional[int] = None
        
        # Read metadata
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('# genome_length='):
                    genome_length = int(line.strip().split('=')[1])
                    break
        
        # Read full dataframe
        events: pd.DataFrame = pd.read_csv(filepath, sep='\t', comment='#')
        
        if genome_length is None:
            raise ValueError(f"No genome_length metadata found in {filepath}")
        
        return events, genome_length


def read_intermediate(filepath: str) -> Tuple[pd.DataFrame, int]:
    """
    Convenience function to read intermediate format
    
    Args:
        filepath: Path to intermediate TSV file
        
    Returns:
        Tuple of (events_df, genome_length)
    """
    return IntermediateReader.read(filepath)


class GeneAnnotationReader:
    """Read mitochondrial gene annotations from BED file"""
    
    @staticmethod
    def load_gene_annotations(bed_file: str) -> List[GeneAnnotation]:
        """
        Load gene annotations from BED file
        
        Expected format:
        chr  start  end  name  score  strand  thickStart  thickEnd  itemRgb
        chrM 576    647  MT-TF 0      +       576         647       255,255,0
        
        Args:
            bed_file: Path to BED file with gene annotations
        
        Returns:
            List of dicts with gene info including chr, start, end, name, color
        """
        annotations: List[GeneAnnotation] = []
        
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Parse RGB color
                rgb: Tuple[float, float, float] = tuple(
                    int(x)/255.0 for x in fields[8].split(',')
                )  # type: ignore
                
                gene_annotation: GeneAnnotation = {
                    'chr': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'name': fields[3],
                    'strand': fields[5] if len(fields) > 5 else '+',
                    'color': rgb
                }
                
                annotations.append(gene_annotation)
        
        return annotations