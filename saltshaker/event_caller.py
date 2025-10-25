"""
Event caller - core MitoSAlt algorithm

Implements the del/dup classification from Basu et al. PLoS Genetics 2020.
Direct port of delplot.R logic - must remain functionally identical.
"""

from __future__ import annotations
from typing import Tuple, List, Optional
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import warnings
import logging
import traceback

warnings.filterwarnings('ignore')

logger = logging.getLogger(__name__)


class EventCaller:
    """
    Calls mitochondrial structural alterations as deletions or duplications
    
    Core algorithm from MitoSAlt (Basu et al. 2020) that classifies events
    based on overlap with heavy and light strand replication origins.
    """
    
    @staticmethod
    def _parse_comma_separated_median(x: str) -> float:
        """
        Calculate median from comma-separated integer string
        
        Args:
            x: Comma-separated integer string
            
        Returns:
            Median value
        """
        return float(np.median([int(i) for i in str(x).split(',')]))
    
    @staticmethod
    def _parse_comma_separated_min(x: str) -> int:
        """
        Calculate minimum from comma-separated integer string
        
        Args:
            x: Comma-separated integer string
            
        Returns:
            Minimum value
        """
        return min([int(i) for i in str(x).split(',')])
    
    @staticmethod
    def _parse_comma_separated_max(x: str) -> int:
        """
        Calculate maximum from comma-separated integer string
        
        Args:
            x: Comma-separated integer string
            
        Returns:
            Maximum value
        """
        return max([int(i) for i in str(x).split(',')])
    
    def __init__(
        self,
        genome_length: int,
        ori_h: Tuple[int, int],
        ori_l: Tuple[int, int],
        heteroplasmy_limit: float = 0.01,
        flank_size: int = 15
    ) -> None:
        """
        Initialize EventCaller
        
        Args:
            genome_length: Mitochondrial genome length (e.g., 16569 for human)
            ori_h: Tuple of (start, end) for heavy strand origin
            ori_l: Tuple of (start, end) for light strand origin
            heteroplasmy_limit: Minimum heteroplasmy to include (default: 0.01 = 1%)
            flank_size: Bases to extract around breakpoints (default: 15)
        """
        self.genome_length: int = genome_length
        self.ori_h: Tuple[int, int] = ori_h
        self.ori_l: Tuple[int, int] = ori_l
        self.heteroplasmy_limit: float = heteroplasmy_limit
        self.flank_size: int = flank_size
    
    def call_events(self, cluster_file: str, breakpoint_file: str) -> pd.DataFrame:
        """
        Load cluster/breakpoint data and call events as deletions or duplications
        
        This implements the R script logic from delplot.R
        
        Args:
            cluster_file: Path to MitoSAlt cluster file
            breakpoint_file: Path to MitoSAlt breakpoint file
            
        Returns:
            DataFrame with called events (del/dup classification)
        """
        # Load and merge data (R script lines ~90-160)
        events = self._load_and_merge_data(cluster_file, breakpoint_file)
        
        if len(events) == 0:
            return pd.DataFrame()
        
        # Classify as deletion or duplication (R script lines ~160-240)
        events = self._classify_del_or_dup(events)
        
        return events
    
    def _load_and_merge_data(
        self,
        cluster_file: str,
        breakpoint_file: str
    ) -> pd.DataFrame:
        """
        Load cluster and breakpoint files and merge them
        
        Direct port of R script logic for loading and merging data
        
        Args:
            cluster_file: Path to cluster file
            breakpoint_file: Path to breakpoint file
            
        Returns:
            Merged DataFrame with events
        """
        # Check if cluster file is empty
        if Path(cluster_file).stat().st_size == 0:
            logger.warning("Empty cluster file, no events to plot")
            return pd.DataFrame()
        
        # Load cluster data
        logger.info(f"Loading cluster file: {cluster_file}")
        clusters = pd.read_csv(cluster_file, sep='\t', header=None)
        clusters.columns = ['cluster', 'read', 'del_start', 'del_end', 'lfstart', 'lfend', 'nread', 'tread', 'perc']
        clusters = clusters[~clusters['cluster'].isna()]
        
        # Extract sample name
        sample_name: str = Path(cluster_file).name.replace('.cluster', '')
        clusters['sample'] = sample_name
        
        logger.info(f"Loaded {len(clusters)} clusters")
        
        # Calculate medians and ranges
        clusters['del_start_median'] = clusters['del_start'].apply(self._parse_comma_separated_median)
        clusters['del_end_median'] = clusters['del_end'].apply(self._parse_comma_separated_median)
        clusters['del_start_min'] = clusters['del_start'].apply(self._parse_comma_separated_min)
        clusters['del_start_max'] = clusters['del_start'].apply(self._parse_comma_separated_max)
        clusters['del_end_min'] = clusters['del_end'].apply(self._parse_comma_separated_min)
        clusters['del_end_max'] = clusters['del_end'].apply(self._parse_comma_separated_max)
        
        # Create range strings (with spaces around dash)
        clusters['del_start_range'] = clusters['del_start_min'].astype(str) + ' - ' + clusters['del_start_max'].astype(str)
        clusters['del_end_range'] = clusters['del_end_min'].astype(str) + ' - ' + clusters['del_end_max'].astype(str)
        
        # Create list.reads exactly like R script
        list_reads: List[str] = []
        length_list_reads: List[int] = []
        
        for idx, row in clusters.iterrows():
            reads = row['read'].split(',')
            reads = [r.strip() for r in reads]
            list_reads.extend(reads)
            length_list_reads.append(len(reads))
        
        # Create res.read 
        res_read_data: List[dict] = []
        for i, row in clusters.iterrows():
            sample = row['sample'] 
            cluster = row['cluster']
            reads = row['read'].split(',')
            for read in reads:
                res_read_data.append({
                    'sample': sample,
                    'cluster': cluster,
                    'read': read.strip()
                })
        
        res_read = pd.DataFrame(res_read_data)
        logger.debug(f"Expanded to {len(res_read)} individual read records")
        
        # Load breakpoint data exactly like R script
        logger.info(f"Loading breakpoint file: {breakpoint_file}")
        bp_raw = pd.read_csv(breakpoint_file, sep='\t', header=None)
        
        # Take columns 2,4,5,10 (R uses 1-based indexing, so Python is 1,3,4,9)
        bp = bp_raw.iloc[:, [1, 3, 4, 9]].copy()
        bp.columns = ['read', 'del_start', 'del_end', 'dloop']
        bp = bp[~bp['read'].isna()]
        
        logger.info(f"Loaded {len(bp)} breakpoint records")
        
        # Merge res.read with bp 
        res_read_bp = pd.merge(res_read, bp, on='read', how='inner')
        logger.debug(f"After merge with breakpoints: {len(res_read_bp)} records")
        
        if len(res_read_bp) == 0:
            logger.warning("No matching reads found between clusters and breakpoints")
            return pd.DataFrame()
        
        # Get unique combinations 
        res_read_bp1 = res_read_bp[['sample', 'cluster', 'dloop']].drop_duplicates()
        logger.debug(f"Unique cluster-dloop combinations: {len(res_read_bp1)}")
        
        # Final merge 
        final_clusters = pd.merge(clusters, res_read_bp1, on=['sample', 'cluster'], how='inner')
        logger.info(f"Final clusters after merge: {len(final_clusters)}")
        
        if len(final_clusters) == 0:
            logger.warning("No clusters remained after merging with breakpoint data")
            return pd.DataFrame()
        
        # Calculate delsize 
        final_clusters['delsize'] = final_clusters['del_end_median'] - final_clusters['del_start_median']
        
        # Handle dloop wraparound
        dloop_mask = final_clusters['dloop'] == 'yes'
        if dloop_mask.any():
            final_clusters.loc[dloop_mask, 'delsize'] = (
                self.genome_length - final_clusters.loc[dloop_mask, 'del_end_median'] + 
                final_clusters.loc[dloop_mask, 'del_start_median']
            )
        
        logger.debug(f"Calculated delsize for {len(final_clusters)} events")
        
        # Filter by heteroplasmy 
        final_clusters = final_clusters[final_clusters['perc'] >= self.heteroplasmy_limit]
        logger.info(f"Events after heteroplasmy filter (>={self.heteroplasmy_limit}): {len(final_clusters)}")
        
        if len(final_clusters) == 0:
            logger.warning("No events above heteroplasmy threshold")
            return pd.DataFrame()
        
        return final_clusters
    
    def _classify_del_or_dup(self, clusters: pd.DataFrame) -> pd.DataFrame:
        """
        Classify events as deletions or duplications based on origin overlap
        
        Direct port of R script lines ~160-240 in delplot.R
        CRITICAL: This is the published MitoSAlt algorithm
        
        Algorithm:
        1. Start with all events classified as deletions
        2. Check overlap with OriH (heavy strand origin)
        3. Check overlap with OriL (light strand origin)
        4. Events overlapping origins are reclassified as duplications
        
        Args:
            clusters: DataFrame with event data
            
        Returns:
            DataFrame with final_event column added
        """
        logger.info(f"Starting classification with {len(clusters)} events")
        logger.debug(f"OriH: {self.ori_h}, OriL: {self.ori_l}")
        
        # Classification: Start with all as deletions
        clusters['final_event'] = 'del'
        
        # First loop: OriH classification (no coordinate swapping)
        for i in clusters.index:
            Rs: int = self.ori_h[0]  # ohs
            Re: int = self.ori_h[1]  # ohe
            Ds: float = clusters.loc[i, 'del_start_median']
            De: float = clusters.loc[i, 'del_end_median']
            
            # R script OriH logic (no swapping of coordinates)
            if Re >= Rs:
                # Essential region NOT covering pos 0
                if ((Ds >= Rs) and (Ds <= Re)) or ((De >= Rs) and (De <= Re)):
                    clusters.loc[i, 'final_event'] = 'dup'
                elif (De > Ds) and (Ds <= Rs) and (De >= Re):
                    clusters.loc[i, 'final_event'] = 'dup'
                elif (De < Ds) and ((De >= Re) or (Ds <= Rs)):
                    clusters.loc[i, 'final_event'] = 'dup'
            else:
                # Essential region IS covering pos 0
                if (Ds >= Rs) or (Ds <= Re) or (De >= Rs) or (De <= Re):
                    clusters.loc[i, 'final_event'] = 'dup'
                elif (De < Ds):
                    clusters.loc[i, 'final_event'] = 'dup'
        
        after_orih: int = int((clusters['final_event'] == 'dup').sum())
        logger.debug(f"After OriH: {after_orih} events classified as dup")
        
        # Second loop: OriL classification (coordinate swapping for dloop=='yes')
        for i in clusters.index:
            dloop: str = str(clusters.loc[i, 'dloop'])
            Rs = self.ori_l[0]  # ols
            Re = self.ori_l[1]  # ole
            
            # Coordinate assignment
            if dloop == 'yes':
                Ds = clusters.loc[i, 'del_end_median']    # Swapped
                De = clusters.loc[i, 'del_start_median']  # Swapped
            else:
                Ds = clusters.loc[i, 'del_start_median']
                De = clusters.loc[i, 'del_end_median']
            
            # Same overlap logic as OriH
            if Re >= Rs:
                # Essential region NOT covering pos 0
                if ((Ds >= Rs) and (Ds <= Re)) or ((De >= Rs) and (De <= Re)):
                    clusters.loc[i, 'final_event'] = 'dup'
                elif (De > Ds) and (Ds <= Rs) and (De >= Re):
                    clusters.loc[i, 'final_event'] = 'dup'
                elif (De < Ds) and ((De >= Re) or (Ds <= Rs)):
                    clusters.loc[i, 'final_event'] = 'dup'
            else:
                # Essential region IS covering pos 0
                if (Ds >= Rs) or (Ds <= Re) or (De >= Rs) or (De <= Re):
                    clusters.loc[i, 'final_event'] = 'dup'
                elif (De < Ds):
                    clusters.loc[i, 'final_event'] = 'dup'
        
        # Final counts
        del_count: int = int((clusters['final_event'] == 'del').sum())
        dup_count: int = int((clusters['final_event'] == 'dup').sum())
        logger.info(f"Final classification - Deletions: {del_count}, Duplications: {dup_count}")
        
        return clusters
    
    def add_flanking_sequences(
        self,
        events: pd.DataFrame,
        genome_fasta: str
    ) -> pd.DataFrame:
        """
        Extract flanking sequences around breakpoints
        
        Args:
            events: DataFrame with events
            genome_fasta: Path to reference genome FASTA
            
        Returns:
            DataFrame with seq1, seq2, seq columns added
        """
        if not genome_fasta or not Path(genome_fasta).exists():
            events['seq1'] = 'NA'
            events['seq2'] = 'NA'
            events['seq'] = 'NA'
            return events
        
        try:
            genome_record = next(SeqIO.parse(genome_fasta, 'fasta'))
            genome_seq: str = str(genome_record.seq).upper()
                        
            flanking_results: pd.DataFrame = self._align_breakpoints(
                events['final_start'] - 1,
                events['final_end'],
                genome_seq,
                self.flank_size
            )
                        
            # Reset index to ensure alignment
            events = events.reset_index(drop=True)
            flanking_results = flanking_results.reset_index(drop=True)
            
            events['seq1'] = flanking_results['seq1']
            events['seq2'] = flanking_results['seq2']
            events['seq'] = flanking_results['seq']
            
            logger.debug(f"  Sample seq values in events after assignment: {events['seq'].head().tolist()}")
            
        except Exception as e:
            logger.warning(f"Flanking sequence analysis failed: {e}")
            traceback.print_exc()
            events['seq1'] = 'NA'
            events['seq2'] = 'NA'
            events['seq'] = 'NA'
        
        return events
    
    def _align_breakpoints(
        self,
        starts: pd.Series,
        ends: pd.Series,
        genome_seq: str,
        flank_size: int = 15
    ) -> pd.DataFrame:
        """
        Python implementation of R align.bp function with DEBUG OUTPUT
        
        Args:
            starts: Series of start positions (0-based)
            ends: Series of end positions (0-based)
            genome_seq: Reference genome sequence
            flank_size: Number of bases to extract on each side
            
        Returns:
            DataFrame with seq1, seq2, seq columns
        """
        results: List[dict] = []
        
        for idx, (start, end) in enumerate(zip(starts, ends)):
            # Convert to 1-based coordinates like R
            start_1based: int = int(start) + 1
            end_1based: int = int(end) + 1
                    
            # Extract sequences - R: substr(mt.fa, start-nb, start+nb)
            bp_start_genome: int = max(1, start_1based - flank_size)
            bp_end_genome: int = min(len(genome_seq), start_1based + flank_size)
            bp: str = genome_seq[bp_start_genome-1:bp_end_genome]
            
            bp1_start_genome: int = max(1, end_1based - flank_size)
            bp1_end_genome: int = min(len(genome_seq), end_1based + flank_size)
            bp1: str = genome_seq[bp1_start_genome-1:bp1_end_genome]
                        
            # Build display strings
            bp_1: str = genome_seq[max(0, start_1based-flank_size-1):start_1based-1]
            bp_2: str = genome_seq[start_1based:min(len(genome_seq), start_1based+flank_size)]
            bp_res: str = f"{bp_1}*{bp_2}"
            
            bp1_1: str = genome_seq[max(0, end_1based-flank_size-1):end_1based-1]
            bp1_2: str = genome_seq[end_1based:min(len(genome_seq), end_1based+flank_size)]
            bp1_res: str = f"{bp1_1}*{bp1_2}"
            
            # Pattern matching
            a: str = bp.replace('N', 'A')
            b: str = bp1.replace('N', 'A')
            seq_result: str = "NA"
            
            if len(a) < 31 or len(b) < 31:
                logger.debug("SKIP: Sequences too short (need 31bp)")
                results.append({'seq1': bp_res, 'seq2': bp1_res, 'seq': seq_result})
                continue
                        
            size: int = 15
            found: bool = False
            
            for i in range(1, 14):
                if found:
                    break
                size = size - 1
                
                posa: int = 1
                posb: int = 1 + size - 1
                
                tested_count: int = 0
                while posb <= 31:
                    if found:
                        break
                    
                    if posa <= 17 and posb >= 15:
                        tmp: str = a[posa-1:posb]
                        
                        # Find matches in b
                        matches: List[Tuple[int, int]] = []
                        search_pos: int = 0
                        while True:
                            match_idx: int = b.find(tmp, search_pos)
                            if match_idx == -1:
                                break
                            start1: int = match_idx + 1  # Convert to R 1-based
                            end1: int = match_idx + len(tmp)  # R 1-based inclusive end
                            matches.append((start1, end1))
                            search_pos = match_idx + 1
                        
                        if tested_count == 0:  # Show first attempt at each size
                            if matches:
                                for s1, e1 in matches:
                                    overlaps: bool = s1 <= 17 and e1 >= 15
                        
                        tested_count += 1
                        
                        if len(matches) >= 1:
                            for start1, end1 in matches:
                                if start1 <= 17 and end1 >= 15:
                                    found = True
                                    seq_result = tmp
                                    break
                    
                    posa = posa + 1
                    posb = posb + 1
            
            results.append({
                'seq1': bp_res,
                'seq2': bp1_res,
                'seq': seq_result
            })
        
        return pd.DataFrame(results)