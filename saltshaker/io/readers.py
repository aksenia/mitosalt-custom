"""
I/O Readers

Handles reading of various input files.
"""

from importlib.resources import path
import pandas as pd
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class BlacklistReader:
    """Reads and parses blacklist BED files"""
    
    @staticmethod
    def load_blacklist_regions(blacklist_file):
        """
        Robustly load blacklist regions from BED file
        
        Tries multiple parsing strategies to handle various BED formats.
        
        Args:
            blacklist_file: Path to BED file with regions to exclude
            
        Returns:
            List of dicts with 'chr', 'start', 'end' keys
        """
        blacklist_regions = []
        if not blacklist_file or not Path(blacklist_file).exists():
            return blacklist_regions
        
        try:
            blacklist = None
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
                blacklist_regions = blacklist.to_dict('records')
        
        except Exception as e:
            logger.warning(f"Could not load blacklist file: {e}")
        
        return blacklist_regions
    

class IntermediateReader:
    """Reads events in internal format"""
    
    @staticmethod
    def read(filepath):
        """
        Read events with metadata
        
        Returns:
            Tuple of (events_df, genome_length)
        """
        genome_length = None
        
        # Read metadata
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('# genome_length='):
                    genome_length = int(line.strip().split('=')[1])
                    break
        
        # Read full dataframe
        events = pd.read_csv(filepath, sep='\t', comment='#')
        
        return events, genome_length

def read_intermediate(filepath):
    """Convenience function"""
    return IntermediateReader.read(filepath)

class GeneAnnotationReader:
    """Read mitochondrial gene annotations from BED file"""
    
    @staticmethod
    def load_gene_annotations(bed_file):
        """
        Load gene annotations from BED file
        
        Expected format:
        chr  start  end  name  score  strand  thickStart  thickEnd  itemRgb
        chrM 576    647  MT-TF 0      +       576         647       255,255,0
        
        Returns:
            List of dicts with gene info
        """
        annotations = []
        
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Parse RGB color
                rgb = tuple(int(x)/255.0 for x in fields[8].split(','))
                
                annotations.append({
                    'chr': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'name': fields[3],
                    'strand': fields[5] if len(fields) > 5 else '+',
                    'color': rgb
                })
        
        return annotations