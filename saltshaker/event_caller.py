"""
Event Caller - Core MitoSAlt Algorithm

Implements the del/dup classification from Basu et al. PLoS Genetics 2020.
Direct port of delplot.R logic - must remain functionally identical.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import warnings
import logging

warnings.filterwarnings('ignore')

logger = logging.getLogger(__name__)


class EventCaller:
    """
    Calls mitochondrial structural alterations as deletions or duplications
    
    Core algorithm from MitoSAlt (Basu et al. 2020) that classifies events
    based on overlap with heavy and light strand replication origins.
    """
    
    def __init__(self, genome_length, ori_h, ori_l, heteroplasmy_limit=0.01, flank_size=15):
        """
        Initialize EventCaller
        
        Args:
            genome_length: Mitochondrial genome length (e.g., 16569 for human)
            ori_h: Tuple of (start, end) for heavy strand origin
            ori_l: Tuple of (start, end) for light strand origin
            heteroplasmy_limit: Minimum heteroplasmy to include (default: 0.01 = 1%)
            flank_size: Bases to extract around breakpoints (default: 15)
        """
        self.genome_length = genome_length
        self.ori_h = ori_h
        self.ori_l = ori_l
        self.heteroplasmy_limit = heteroplasmy_limit
        self.flank_size = flank_size
    
    def call_events(self, cluster_file, breakpoint_file):
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
    
    def _load_and_merge_data(self, cluster_file, breakpoint_file):
        """
        Load cluster and breakpoint files and merge them
        
        Direct port of R script logic for loading and merging data
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
        sample_name = Path(cluster_file).name.replace('.cluster', '')
        clusters['sample'] = sample_name
        
        logger.info(f"Loaded {len(clusters)} clusters")
        
        # Calculate medians and ranges
        def calculate_median(x):
            return np.median([int(i) for i in str(x).split(',')])
        
        def calculate_min(x):
            return min([int(i) for i in str(x).split(',')])
            
        def calculate_max(x):
            return max([int(i) for i in str(x).split(',')])
        
        clusters['del_start_median'] = clusters['del_start'].apply(calculate_median)
        clusters['del_end_median'] = clusters['del_end'].apply(calculate_median)
        clusters['del_start_min'] = clusters['del_start'].apply(calculate_min)
        clusters['del_start_max'] = clusters['del_start'].apply(calculate_max)
        clusters['del_end_min'] = clusters['del_end'].apply(calculate_min)
        clusters['del_end_max'] = clusters['del_end'].apply(calculate_max)
        
        # Create range strings (with spaces around dash)
        clusters['del_start_range'] = clusters['del_start_min'].astype(str) + ' - ' + clusters['del_start_max'].astype(str)
        clusters['del_end_range'] = clusters['del_end_min'].astype(str) + ' - ' + clusters['del_end_max'].astype(str)
        
        # Create list.reads exactly like R script
        list_reads = []
        length_list_reads = []
        
        for idx, row in clusters.iterrows():
            reads = row['read'].split(',')
            reads = [r.strip() for r in reads]
            list_reads.extend(reads)
            length_list_reads.append(len(reads))
        
        # Create res.read 
        res_read_data = []
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
    
    def _classify_del_or_dup(self, clusters):
        """
        Classify events as deletions or duplications based on origin overlap
        
        Direct port of R script lines ~160-240 in delplot.R
        CRITICAL: This is the published MitoSAlt algorithm
        
        Algorithm:
        1. Start with all events classified as deletions
        2. Check overlap with OriH (heavy strand origin)
        3. Check overlap with OriL (light strand origin)
        4. Events overlapping origins are reclassified as duplications
        """
        logger.info(f"Starting classification with {len(clusters)} events")
        logger.debug(f"OriH: {self.ori_h}, OriL: {self.ori_l}")
        
        # Classification: Start with all as deletions
        clusters['final_event'] = 'del'
        
        # First loop: OriH classification (no coordinate swapping)
        for i in clusters.index:
            Rs = self.ori_h[0]  # ohs
            Re = self.ori_h[1]  # ohe
            Ds = clusters.loc[i, 'del_start_median']
            De = clusters.loc[i, 'del_end_median']
            
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
        
        after_orih = (clusters['final_event'] == 'dup').sum()
        logger.debug(f"After OriH: {after_orih} events classified as dup")
        
        # Second loop: OriL classification (coordinate swapping for dloop=='yes')
        for i in clusters.index:
            dloop = clusters.loc[i, 'dloop']
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
        del_count = (clusters['final_event'] == 'del').sum()
        dup_count = (clusters['final_event'] == 'dup').sum()
        logger.info(f"Final classification - Deletions: {del_count}, Duplications: {dup_count}")
        
        return clusters
    
    def add_flanking_sequences(self, events, genome_fasta):
        """Extract flanking sequences around breakpoints"""
        if not genome_fasta or not Path(genome_fasta).exists():
            events['seq1'] = 'NA'
            events['seq2'] = 'NA'
            events['seq'] = 'NA'
            return events
        
        try:
            genome_record = next(SeqIO.parse(genome_fasta, 'fasta'))
            genome_seq = str(genome_record.seq).upper()
                        
            flanking_results = self._align_breakpoints(
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
            import traceback
            traceback.print_exc()
            events['seq1'] = 'NA'
            events['seq2'] = 'NA'
            events['seq'] = 'NA'
        
        return events
    
    def _align_breakpoints(self, starts, ends, genome_seq, flank_size=15):
        """
        Python implementation of R align.bp function with DEBUG OUTPUT
        """
        results = []
        
        for idx, (start, end) in enumerate(zip(starts, ends)):
            # Convert to 1-based coordinates like R
            start_1based = int(start) + 1
            end_1based = int(end) + 1
                    
            # Extract sequences - R: substr(mt.fa, start-nb, start+nb)
            bp_start_genome = max(1, start_1based - flank_size)
            bp_end_genome = min(len(genome_seq), start_1based + flank_size)
            bp = genome_seq[bp_start_genome-1:bp_end_genome]
            
            bp1_start_genome = max(1, end_1based - flank_size)
            bp1_end_genome = min(len(genome_seq), end_1based + flank_size)
            bp1 = genome_seq[bp1_start_genome-1:bp1_end_genome]
                        
            # Build display strings
            bp_1 = genome_seq[max(0, start_1based-flank_size-1):start_1based-1]
            bp_2 = genome_seq[start_1based:min(len(genome_seq), start_1based+flank_size)]
            bp_res = f"{bp_1}*{bp_2}"
            
            bp1_1 = genome_seq[max(0, end_1based-flank_size-1):end_1based-1]
            bp1_2 = genome_seq[end_1based:min(len(genome_seq), end_1based+flank_size)]
            bp1_res = f"{bp1_1}*{bp1_2}"
            
            # Pattern matching
            a = bp.replace('N', 'A')
            b = bp1.replace('N', 'A')
            seq_result = "NA"
            
            if len(a) < 31 or len(b) < 31:
                logger.debug("SKIP: Sequences too short (need 31bp)")
                results.append({'seq1': bp_res, 'seq2': bp1_res, 'seq': seq_result})
                continue
                        
            size = 15
            found = False
            
            for i in range(1, 14):
                if found:
                    break
                size = size - 1
                
                posa = 1
                posb = 1 + size - 1
                
                tested_count = 0
                while posb <= 31:
                    if found:
                        break
                    
                    if posa <= 17 and posb >= 15:
                        tmp = a[posa-1:posb]
                        
                        # Find matches in b
                        matches = []
                        search_pos = 0
                        while True:
                            match_idx = b.find(tmp, search_pos)
                            if match_idx == -1:
                                break
                            start1 = match_idx + 1  # Convert to R 1-based
                            end1 = match_idx + len(tmp)  # R 1-based inclusive end
                            matches.append((start1, end1))
                            search_pos = match_idx + 1
                        
                        if tested_count == 0:  # Show first attempt at each size
                            if matches:
                                for s1, e1 in matches:
                                    overlaps = s1 <= 17 and e1 >= 15
                        
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