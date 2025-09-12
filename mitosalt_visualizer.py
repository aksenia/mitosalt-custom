#!/usr/bin/env python3
"""
MitoSAlt Circular Plot Generator - enhanced rewrite of 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import argparse
from pathlib import Path
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')

class MitoPlotter:
    def __init__(self, genome_length, ori_h_start, ori_h_end, ori_l_start, ori_l_end, 
                 heteroplasmy_limit=0.01, flank_size=15):
        """Initialize MitoPlotter with mitochondrial genome parameters"""
        self.genome_length = genome_length
        self.ori_h = (ori_h_start, ori_h_end)
        self.ori_l = (ori_l_start, ori_l_end)
        self.heteroplasmy_limit = heteroplasmy_limit
        self.flank_size = flank_size
        
        # Create color maps for deletions (blues) and duplications (reds)
        self.del_cmap = LinearSegmentedColormap.from_list(
            'deletions', ['#4169E1', '#1E90FF', '#0000CD', '#000080', '#191970'])
        self.dup_cmap = LinearSegmentedColormap.from_list(
            'duplications', ['#FF6347', '#DC143C', '#B22222', '#8B0000', '#800000'])
    
    def load_data(self, cluster_file, breakpoint_file, blacklist_file=None):
        """Load and process cluster and breakpoint data following R script exactly"""
        
        # Check if cluster file is empty
        if Path(cluster_file).stat().st_size == 0:
            print("Empty cluster file, no events to plot")
            return pd.DataFrame()
        
        # Load cluster data exactly like R script
        print(f"Loading cluster file: {cluster_file}")
        clusters = pd.read_csv(cluster_file, sep='\t', header=None)
        clusters.columns = ['cluster', 'read', 'del.start', 'del.end', 'lfstart', 'lfend', 'nread', 'tread', 'perc']
        clusters = clusters[~clusters['cluster'].isna()]
        
        # Extract sample name
        sample_name = Path(cluster_file).name.replace('.cluster', '')
        clusters['sample'] = sample_name
        
        print(f"Loaded {len(clusters)} clusters")
        
        # Calculate medians and ranges exactly like R script
        def calculate_median(x):
            return np.median([int(i) for i in str(x).split(',')])
        
        def calculate_min(x):
            return min([int(i) for i in str(x).split(',')])
            
        def calculate_max(x):
            return max([int(i) for i in str(x).split(',')])
        
        clusters['del.start.median'] = clusters['del.start'].apply(calculate_median)
        clusters['del.end.median'] = clusters['del.end'].apply(calculate_median)
        clusters['del.start.min'] = clusters['del.start'].apply(calculate_min)
        clusters['del.start.max'] = clusters['del.start'].apply(calculate_max)
        clusters['del.end.min'] = clusters['del.end'].apply(calculate_min)
        clusters['del.end.max'] = clusters['del.end'].apply(calculate_max)
        
        # Create range strings exactly like R script (with spaces around dash)
        clusters['del.start.range'] = clusters['del.start.min'].astype(str) + ' - ' + clusters['del.start.max'].astype(str)
        clusters['del.end.range'] = clusters['del.end.min'].astype(str) + ' - ' + clusters['del.end.max'].astype(str)
        
        # Create list.reads exactly like R script: strsplit(res$read,",",fixed=T)
        list_reads = []
        length_list_reads = []
        
        for idx, row in clusters.iterrows():
            reads = row['read'].split(',')
            reads = [r.strip() for r in reads]  # Remove any whitespace
            list_reads.extend(reads)
            length_list_reads.append(len(reads))
        
        # Create res.read exactly like R script
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
        print(f"Expanded to {len(res_read)} individual read records")
        
        # Load breakpoint data exactly like R script
        # bp <- read.delim(bpfile, header=FALSE)[,c(2,4,5,10)]
        # colnames(bp)<-c("read","del.start","del.end","dloop")
        print(f"Loading breakpoint file: {breakpoint_file}")
        bp_raw = pd.read_csv(breakpoint_file, sep='\t', header=None)
        
        # Take columns 2,4,5,10 (R uses 1-based indexing, so Python is 1,3,4,9)
        bp = bp_raw.iloc[:, [1, 3, 4, 9]].copy()
        bp.columns = ['read', 'del.start', 'del.end', 'dloop']
        bp = bp[~bp['read'].isna()]
        
        print(f"Loaded {len(bp)} breakpoint records")
        print(f"Sample breakpoint data:\n{bp.head()}")
        print(f"Unique dloop values: {bp['dloop'].unique()}")
        
        # Merge res.read with bp exactly like R script
        # res.read.bp<-merge(res.read,bp,by="read")
        res_read_bp = pd.merge(res_read, bp, on='read', how='inner')
        print(f"After merge with breakpoints: {len(res_read_bp)} records")
        
        if len(res_read_bp) == 0:
            print("No matching reads found between clusters and breakpoints")
            return pd.DataFrame()
        
        # Get unique combinations exactly like R script
        # res.read.bp1<-unique(res.read.bp[,c(2,3,6)])  # sample, cluster, dloop
        res_read_bp1 = res_read_bp[['sample', 'cluster', 'dloop']].drop_duplicates()
        print(f"Unique cluster-dloop combinations: {len(res_read_bp1)}")
        
        # Final merge exactly like R script
        # res<-merge(res,res.read.bp1,by=c("sample","cluster"))
        final_clusters = pd.merge(clusters, res_read_bp1, on=['sample', 'cluster'], how='inner')
        print(f"Final clusters after merge: {len(final_clusters)}")
        
        if len(final_clusters) == 0:
            print("No clusters remained after merging with breakpoint data")
            return pd.DataFrame()
        
        # Calculate delsize exactly like R script
        final_clusters['delsize'] = final_clusters['del.end.median'] - final_clusters['del.start.median']
        
        # Handle dloop wraparound exactly like R script
        dloop_mask = final_clusters['dloop'] == 'yes'
        if dloop_mask.any():
            final_clusters.loc[dloop_mask, 'delsize'] = (
                self.genome_length - final_clusters.loc[dloop_mask, 'del.end.median'] + 
                final_clusters.loc[dloop_mask, 'del.start.median']
            )
        
        print(f"Calculated delsize for {len(final_clusters)} events")
        print(f"dloop=='yes' events: {dloop_mask.sum()}")
        
        # Filter by heteroplasmy exactly like R script
        # if(nrow(res[res$perc>=hp.limit,])>0)
        final_clusters = final_clusters[final_clusters['perc'] >= self.heteroplasmy_limit]
        print(f"Events after heteroplasmy filter (>={self.heteroplasmy_limit}): {len(final_clusters)}")
        
        if len(final_clusters) == 0:
            print("No events above heteroplasmy threshold")
            return pd.DataFrame()
        
        # Apply blacklist filter if provided
        if blacklist_file and Path(blacklist_file).exists():
            try:
                blacklist = pd.read_csv(blacklist_file, sep='\t', header=None)
                if len(blacklist.columns) >= 3:
                    blacklist = blacklist.iloc[:, :3]
                    blacklist.columns = ['chr', 'start', 'end']
                    final_clusters = self._filter_blacklist(final_clusters, blacklist)
                    print(f"Applied blacklist filter")
            except Exception as e:
                print(f"Warning: Could not load blacklist file: {e}")
        
        # Classification: Start with all as deletions, then apply R script logic
        final_clusters['final.event'] = 'del'
        final_clusters = self._apply_origin_classification(final_clusters)
        
        return final_clusters
    
    def _filter_blacklist(self, clusters, blacklist):
        """Filter out events that overlap with blacklisted regions"""
        def overlaps_blacklist(row):
            start, end = row['del.start.median'], row['del.end.median']
            for _, region in blacklist.iterrows():
                bl_start, bl_end = region['start'], region['end']
                if not (end < bl_start or start > bl_end):
                    return True
            return False
        
        mask = ~clusters.apply(overlaps_blacklist, axis=1)
        filtered_count = len(clusters) - mask.sum()
        if filtered_count > 0:
            print(f"Filtered out {filtered_count} events overlapping blacklisted regions")
        
        return clusters[mask]
    
    def _apply_origin_classification(self, clusters):
        """
        Apply origin overlap classification exactly like R script
        """
        print(f"Starting classification with {len(clusters)} events")
        print(f"OriH: {self.ori_h}, OriL: {self.ori_l}")
        
        # Track changes for debugging
        del_count_start = (clusters['final.event'] == 'del').sum()
        
        # First loop: OriH classification (no coordinate swapping)
        for i in clusters.index:
            Rs = self.ori_h[0]  # ohs
            Re = self.ori_h[1]  # ohe
            Ds = clusters.loc[i, 'del.start.median']
            De = clusters.loc[i, 'del.end.median']
            
            # R script OriH logic (no swapping of coordinates)
            if Re >= Rs:
                # Essential region NOT covering pos 0
                if ((Ds >= Rs) and (Ds <= Re)) or ((De >= Rs) and (De <= Re)):
                    clusters.loc[i, 'final.event'] = 'dup'
                elif (De > Ds) and (Ds <= Rs) and (De >= Re):
                    clusters.loc[i, 'final.event'] = 'dup'
                elif (De < Ds) and ((De >= Re) or (Ds <= Rs)):
                    clusters.loc[i, 'final.event'] = 'dup'
            else:
                # Essential region IS covering pos 0
                if (Ds >= Rs) or (Ds <= Re) or (De >= Rs) or (De <= Re):
                    clusters.loc[i, 'final.event'] = 'dup'
                elif (De < Ds):
                    clusters.loc[i, 'final.event'] = 'dup'
        
        after_orih = (clusters['final.event'] == 'dup').sum()
        print(f"After OriH: {after_orih} events classified as dup")
        
        # Second loop: OriL classification (coordinate swapping for dloop=='yes')
        for i in clusters.index:
            dloop = clusters.loc[i, 'dloop']
            Rs = self.ori_l[0]  # ols
            Re = self.ori_l[1]  # ole
            
            # Coordinate assignment exactly like R script
            if dloop == 'yes':
                Ds = clusters.loc[i, 'del.end.median']    # Swapped
                De = clusters.loc[i, 'del.start.median']  # Swapped
            else:
                Ds = clusters.loc[i, 'del.start.median']
                De = clusters.loc[i, 'del.end.median']
            
            # Same overlap logic as OriH
            if Re >= Rs:
                # Essential region NOT covering pos 0
                if ((Ds >= Rs) and (Ds <= Re)) or ((De >= Rs) and (De <= Re)):
                    clusters.loc[i, 'final.event'] = 'dup'
                elif (De > Ds) and (Ds <= Rs) and (De >= Re):
                    clusters.loc[i, 'final.event'] = 'dup'
                elif (De < Ds) and ((De >= Re) or (Ds <= Rs)):
                    clusters.loc[i, 'final.event'] = 'dup'
            else:
                # Essential region IS covering pos 0
                if (Ds >= Rs) or (Ds <= Re) or (De >= Rs) or (De <= Re):
                    clusters.loc[i, 'final.event'] = 'dup'
                elif (De < Ds):
                    clusters.loc[i, 'final.event'] = 'dup'
        
        # Final counts
        del_count = (clusters['final.event'] == 'del').sum()
        dup_count = (clusters['final.event'] == 'dup').sum()
        print(f"Final classification - Deletions: {del_count}, Duplications: {dup_count}")
        
        return clusters

    def _classify_event_pattern(self, events):
        """
        Classify events as Single or Multiple based on MitoSAlt literature criteria:
        - Single: Dominated by one or few high-heteroplasmy events (typically >30-35%)
        - Multiple: Many low-heteroplasmy events, often representing different mechanisms
        """
        if len(events) == 0:
            return "Unknown", "No events detected", {}
        
        # Key biological thresholds based on literature
        HIGH_HETEROPLASMY_THRESHOLD = 0.30  # 30% - major pathogenic events
        LOW_HETEROPLASMY_THRESHOLD = 0.01   # 1% - potential artifacts/minor events
        MAJOR_EVENT_COUNT_THRESHOLD = 3     # Few major events = single pattern
        TOTAL_EVENT_COUNT_THRESHOLD = 20    # Many events = multiple pattern
        
        # Heteroplasmy-based classification
        high_het_events = events[events['perc'] >= HIGH_HETEROPLASMY_THRESHOLD]
        medium_het_events = events[(events['perc'] >= LOW_HETEROPLASMY_THRESHOLD) & 
                                (events['perc'] < HIGH_HETEROPLASMY_THRESHOLD)]
        low_het_events = events[events['perc'] < LOW_HETEROPLASMY_THRESHOLD]
        
        # Count event types
        del_count = (events['final.event'] == 'del').sum()
        dup_count = (events['final.event'] == 'dup').sum()
        
        # Spatial metrics
        positions = []
        for _, event in events.iterrows():
            start_pos = event['del.start.median']
            end_pos = event['del.end.median']
            positions.extend([start_pos, end_pos])
        
        positions = np.array(positions)
        position_range = np.max(positions) - np.min(positions)
        genome_coverage = position_range / self.genome_length
        
        # Size variation
        size_std = np.std(events['delsize'])
        size_cv = size_std / np.mean(events['delsize']) if np.mean(events['delsize']) > 0 else 0
        
        # Classification criteria - focus on literature-validated metrics only
        criteria = {
            # Core biological metrics from MitoSAlt literature
            'total_events': len(events),
            'high_het_count': len(high_het_events),
            'medium_het_count': len(medium_het_events), 
            'low_het_count': len(low_het_events),
            'max_heteroplasmy': events['perc'].max(),
            'median_heteroplasmy': events['perc'].median(),
            'del_count': del_count,
            'dup_count': dup_count,
            
            # Secondary metrics (may be useful for detailed analysis)
            'size_coefficient_variation': size_cv,
            'position_range': position_range,
            'size_std': size_std,
            'many_events': len(events) > TOTAL_EVENT_COUNT_THRESHOLD,
            'mixed_types': del_count > 0 and dup_count > 0 and min(del_count, dup_count) >= 3
        }
        
        # Decision logic based on MitoSAlt literature patterns
        
        # SINGLE EVENT PATTERN criteria:
        # 1. Dominated by high-heteroplasmy events (like patient samples in literature)
        # 2. Few major events with potential low-level companions
        # 3. Clear dominant event(s)
        
        single_pattern_indicators = [
            # Criterion 1: Dominant high-heteroplasmy event(s)
            len(high_het_events) >= 1 and len(high_het_events) <= MAJOR_EVENT_COUNT_THRESHOLD,
            
            # Criterion 2: High max heteroplasmy suggests pathogenic event
            criteria['max_heteroplasmy'] >= HIGH_HETEROPLASMY_THRESHOLD,
            
            # Criterion 3: Most events are artifacts if many low-het events with few high-het
            len(high_het_events) >= 1 and len(events) <= TOTAL_EVENT_COUNT_THRESHOLD,
            
            # Criterion 4: Single type dominance with high heteroplasmy
            (del_count == 0 or dup_count == 0) and criteria['max_heteroplasmy'] >= HIGH_HETEROPLASMY_THRESHOLD
        ]
        
        # MULTIPLE EVENT PATTERN criteria:
        # 1. Many events at low heteroplasmy (like mouse models in literature)
        # 2. No clear dominant high-heteroplasmy event
        # 3. Complex patterns with mixed types
        
        multiple_pattern_indicators = [
            # Criterion 1: Many total events (mouse model pattern)
            len(events) > TOTAL_EVENT_COUNT_THRESHOLD,
            
            # Criterion 2: Many events but no high-heteroplasmy dominant event
            len(events) > 10 and len(high_het_events) == 0,
            
            # Criterion 3: Mixed deletion/duplication pattern with substantial counts
            del_count >= 3 and dup_count >= 3 and len(events) > 10,
            
            # Criterion 4: High complexity - many medium heteroplasmy events
            len(medium_het_events) > 10 and criteria['median_heteroplasmy'] < HIGH_HETEROPLASMY_THRESHOLD,
            
            # Criterion 5: Very high event count suggests multiple mechanisms
            len(events) > 50
        ]
        
        # Classification decision
        single_score = sum(single_pattern_indicators)
        multiple_score = sum(multiple_pattern_indicators)
        
        if single_score > multiple_score:
            classification = "Single"
            
            # Detailed reasoning for single pattern
            reasons = []
            if len(high_het_events) >= 1:
                reasons.append(f"dominant high-heteroplasmy event(s) ({len(high_het_events)} at ≥{HIGH_HETEROPLASMY_THRESHOLD:.0%})")
            if criteria['max_heteroplasmy'] >= HIGH_HETEROPLASMY_THRESHOLD:
                reasons.append(f"max heteroplasmy {criteria['max_heteroplasmy']:.1%}")
            if len(low_het_events) > 0:
                reasons.append(f"{len(low_het_events)} low-level events (likely artifacts)")
            
            reason_str = "; ".join(reasons) if reasons else f"clustered pattern ({len(events)} events)"
            
        elif multiple_score > single_score:
            classification = "Multiple"
            
            # Detailed reasoning for multiple pattern
            reasons = []
            if len(events) > TOTAL_EVENT_COUNT_THRESHOLD:
                reasons.append(f"{len(events)} events")
            if len(high_het_events) == 0:
                reasons.append("no dominant high-heteroplasmy events")
            if del_count >= 3 and dup_count >= 3:
                reasons.append(f"mixed types ({del_count} del, {dup_count} dup)")
            if genome_coverage > 0.5:
                reasons.append(f"wide distribution ({genome_coverage:.1%} genome)")
            
            reason_str = "; ".join(reasons)
            
        else:
            # Ambiguous case - use conservative approach
            if criteria['max_heteroplasmy'] >= HIGH_HETEROPLASMY_THRESHOLD:
                classification = "Single"
                reason_str = f"ambiguous but high max heteroplasmy ({criteria['max_heteroplasmy']:.1%})"
            else:
                classification = "Multiple" 
                reason_str = f"ambiguous, multiple low-heteroplasmy events (median {criteria['median_heteroplasmy']:.2%})"
        
        # Add biological pattern classification
        criteria['classification_scores'] = {
            'single_score': single_score,
            'multiple_score': multiple_score
        }
        
        # Pattern subtype based on literature
        if classification == "Single":
            if len(high_het_events) == 1 and len(low_het_events) > 0:
                criteria['subtype'] = "Classic single deletion/duplication with artifacts"
            elif len(high_het_events) == 1:
                criteria['subtype'] = "Pure single event"
            else:
                criteria['subtype'] = "Few major events"
        else:
            if len(events) > 100:
                criteria['subtype'] = "Complex multiple (mouse model-like)"
            elif del_count > 0 and dup_count > 0:
                criteria['subtype'] = "Mixed deletion-duplication pattern"
            else:
                criteria['subtype'] = "Multiple single-type events"
        
        return classification, reason_str, criteria


    def _align_breakpoints(self, starts, ends, genome_seq, flank_size=15):
        """Python implementation of R align.bp function"""
        results = []
        
        for start, end in zip(starts, ends):
            # Convert to 1-based coordinates like R
            start_1based = int(start) + 1
            end_1based = int(end) + 1
            
            # Extract sequences around start breakpoint
            bp_start = max(0, start_1based - flank_size - 1)
            bp_end = min(len(genome_seq), start_1based + flank_size)
            bp = genome_seq[bp_start:bp_end]
            
            bp_1_start = max(0, start_1based - flank_size - 1)
            bp_1_end = start_1based - 1
            bp_1 = genome_seq[bp_1_start:bp_1_end] if bp_1_end > bp_1_start else ""
            
            bp_2_start = start_1based
            bp_2_end = min(len(genome_seq), start_1based + flank_size)
            bp_2 = genome_seq[bp_2_start:bp_2_end]
            
            bp_res = f"{bp_1}*{bp_2}"
            
            # Extract sequences around end breakpoint
            bp1_start = max(0, end_1based - flank_size - 1)
            bp1_end = min(len(genome_seq), end_1based + flank_size)
            bp1 = genome_seq[bp1_start:bp1_end]
            
            bp1_1_start = max(0, end_1based - flank_size - 1)
            bp1_1_end = end_1based - 1
            bp1_1 = genome_seq[bp1_1_start:bp1_1_end] if bp1_1_end > bp1_1_start else ""
            
            bp1_2_start = end_1based
            bp1_2_end = min(len(genome_seq), end_1based + flank_size)
            bp1_2 = genome_seq[bp1_2_start:bp1_2_end]
            
            bp1_res = f"{bp1_1}*{bp1_2}"
            
            # Pattern matching logic (simplified)
            a = bp.replace('N', 'A')
            b = bp1.replace('N', 'A')
            seq_result = "NA"  # Default value
            
            # Store results
            results.append({
                'seq1': bp_res,
                'seq2': bp1_res,
                'seq': seq_result
            })
        
        return pd.DataFrame(results)
    
    def _load_blacklist_regions(self, blacklist_file):
        """Robustly load blacklist regions from file, returns list of dicts"""
        blacklist_regions = []
        if blacklist_file and Path(blacklist_file).exists():
            try:
                blacklist = None
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
                print(f"Warning: Could not load blacklist file: {e}")
        return blacklist_regions
    
    def create_plot(self, events, output_file, figsize=(16, 10), blacklist_file=None):  # Wider figure
        """Create circular plot of mitochondrial events"""
        if len(events) == 0:
            print("No events to plot")
            return
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Load blacklist regions if provided
        try:
            blacklist_regions = self._load_blacklist_regions(blacklist_file)
        except Exception as e:
            print(f"Warning: Could not load blacklist file: {e}")

        # Create data structure for plotting
        dat = pd.DataFrame({
            'chr': 'MT',
            'start': events['del.start.median'],
            'end': events['del.end.median'],
            'value': events['perc'],
            'dloop': events['dloop'],
            'delsize': events['delsize'],
            'final.event': events['final.event']
        })
        
        # Enhanced blacklist crossing detection
        def crosses_blacklist(start_pos, end_pos):
            if not blacklist_regions:
                return False
            for region in blacklist_regions:
                bl_start, bl_end = int(region['start']), int(region['end'])
                if (bl_start <= start_pos <= bl_end) or (bl_start <= end_pos <= bl_end):
                    return True
            return False
        
        dat['blacklist_crossing'] = False
        if blacklist_regions:
            for idx, row in dat.iterrows():
                dat.loc[idx, 'blacklist_crossing'] = crosses_blacklist(row['start'], row['end'])
        
        # Count events
        del_count = (dat['final.event'] == 'del').sum()
        dup_count = (dat['final.event'] == 'dup').sum()
        bl_del_count = ((dat['final.event'] == 'del') & dat['blacklist_crossing']).sum()
        bl_dup_count = ((dat['final.event'] == 'dup') & dat['blacklist_crossing']).sum()
        
        # Add degrees
        dat['deg1'] = 358 * dat['start'] / self.genome_length
        dat['deg2'] = 358 * dat['end'] / self.genome_length
        
        dup_no_dloop_mask = (dat['final.event'] == 'dup') & (dat['dloop'] == 'no')
        if dup_no_dloop_mask.any():
            dat.loc[dup_no_dloop_mask, 'deg1'] = 360 + dat.loc[dup_no_dloop_mask, 'deg1']
        
        # Process duplications
        dat_dup = dat[dat['final.event'] == 'dup'].copy()
        if len(dat_dup) > 0:
            dat_dup['delsize'] = self.genome_length - dat_dup['delsize']
            dat.loc[dat['final.event'] == 'dup', 'delsize'] = dat_dup['delsize']
        
        # Grouped overlap elimination - dels and dups in separate subcircles
        def assign_radii_by_type(data, base_radius=380, radius_diff=5):
            if len(data) == 0:
                return data
            
            data = data.sort_values('delsize', ascending=False).reset_index(drop=True)
            data['radius'] = 0
            
            def events_overlap(event1, event2):
                start1, end1 = event1['deg1'], event1['deg2']
                start2, end2 = event2['deg1'], event2['deg2']
                
                start1, end1, start2, end2 = start1 % 360, end1 % 360, start2 % 360, end2 % 360
                
                if start1 > end1:
                    if start2 > end2:
                        return True
                    else:
                        return (start2 <= end1) or (end2 >= start1)
                elif start2 > end2:
                    return (start1 <= end2) or (end1 >= start2)
                else:
                    return not (end1 < start2 or end2 < start1)
            
            for i in range(len(data)):
                assigned = False
                test_radius = base_radius
                
                while not assigned and test_radius > 120:
                    conflict_found = False
                    
                    for j in range(len(data)):
                        if j != i and data.loc[j, 'radius'] == test_radius:
                            if events_overlap(data.loc[i], data.loc[j]):
                                conflict_found = True
                                break
                    
                    if not conflict_found:
                        data.loc[i, 'radius'] = test_radius
                        assigned = True
                    else:
                        test_radius -= radius_diff
                
                if not assigned:
                    data.loc[i, 'radius'] = 130
            
            return data
        
        # Separate processing for deletions and duplications
        dat_del = dat[dat['final.event'] == 'del'].copy()
        dat_dup = dat[dat['final.event'] == 'dup'].copy()
        
        # Assign radii in separate subcircles
        if len(dat_del) > 0:
            dat_del = assign_radii_by_type(dat_del, base_radius=380, radius_diff=5)
        
        if len(dat_dup) > 0:
            max_del_radius = dat_del['radius'].max() if len(dat_del) > 0 else 380
            dup_base_radius = max_del_radius - 15
            dat_dup = assign_radii_by_type(dat_dup, base_radius=dup_base_radius, radius_diff=5)
        
        # Combine back
        dat_processed = pd.concat([dat_del, dat_dup], ignore_index=True) if len(dat_dup) > 0 else dat_del
        
        # Automatically shrink blacklist circle
        min_event_radius = dat_processed['radius'].min() if len(dat_processed) > 0 else 200
        blacklist_radius = max(50, min_event_radius - 30)
        
        # Create figure with extended area and constrained subplot for circle
        fig = plt.figure(figsize=figsize)
        
        # Create polar subplot in CENTER portion of figure (leaving space for legends)
        ax = fig.add_subplot(111, projection='polar', position=[0.15, 0.05, 0.7, 0.9])  # [left, bottom, width, height]
        ax.set_ylim(0, 400)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        
        # Draw genome circle
        circle = patches.Circle((0, 0), 390, fill=False, linewidth=3,
                            color='gray', transform=ax.transData._b)
        ax.add_patch(circle)
        
        # Add blacklist regions
        if blacklist_regions:
            separator_circle = patches.Circle((0, 0), blacklist_radius + 15, fill=False, linewidth=2, 
                                            color='lightgray', linestyle='--', alpha=0.7,
                                            transform=ax.transData._b)
            ax.add_patch(separator_circle)
            
            for region in blacklist_regions:
                start_pos = int(region['start'])
                end_pos = int(region['end'])
                
                start_deg = np.radians(358 * start_pos / self.genome_length)
                end_deg = np.radians(358 * end_pos / self.genome_length)
                
                arc_width = end_deg - start_deg
                min_width = np.radians(2)
                if arc_width < min_width:
                    mid_deg = (start_deg + end_deg) / 2
                    start_deg = mid_deg - min_width/2
                    end_deg = mid_deg + min_width/2
                
                theta = np.linspace(start_deg, end_deg, 50)
                ax.plot(theta, [blacklist_radius]*len(theta), color='black', linewidth=4, alpha=0.9, solid_capstyle='round')
                ax.plot([start_deg, start_deg], [blacklist_radius, 390], color='gray', linewidth=1, linestyle='--', alpha=0.6)
                ax.plot([end_deg, end_deg], [blacklist_radius, 390], color='gray', linewidth=1, linestyle='--', alpha=0.6)
        
        # Add position markers
        positions = np.arange(0, self.genome_length, 1000)
        for pos in positions:
            deg = np.radians(358 * pos / self.genome_length)
            ax.plot([deg, deg], [390, 395], color='gray', linewidth=1)
            ax.text(deg, 410, f'{pos//1000}', ha='center', va='center', 
                fontsize=9, color='gray')
        
        # PURE COLOR GRADIENTS - pale but clearly blue/red
        max_het = dat_processed['value'].max()
        min_het = dat_processed['value'].min()
        
        def get_pure_blue_color(het_val):
            """Pure blue gradient - baby blue to deep blue"""
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5

            # Blue goes from baby blue (~0.8, 0.9, 1) to deep blue (0, 0, 1)
            blue   = 1.0
            red    = 0.8 * (1 - norm)          # fade red out completely
            green  = 0.9 * (1 - norm)          # fade green out completely
            blue   = 1.0 * (1 - norm) + norm   # stays maxed at 1, but norm smooths transition

            return (red, green, blue)

        
        def get_pure_red_color(het_val):
            """Pure red gradient - baby pink to deep red"""
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5

            # Red goes from baby pink (~1.0, 0.8, 0.85) to deep red (1.0, 0, 0)
            red   = 1.0
            green = 0.8 * (1 - norm)   # fades green to 0
            blue  = 0.85 * (1 - norm)  # fades blue to 0

            return (red, green, blue)

        
        def get_continuous_alpha(het_val):
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5
            return 0.5 + 0.45 * norm
        
        # Plot events with pure color gradients
        for i, (_, event) in enumerate(dat_processed.iterrows()):
            deg1_rad = np.radians(event['deg1'])
            deg2_rad = np.radians(event['deg2'])
            radius = event['radius']
            het_val = event['value']
            
            if blacklist_regions and event['blacklist_crossing']:
                color = (0.2, 0.8, 0.2)  # Bright lime green
                alpha = 0.9
            else:
                if event['final.event'] == 'del':
                    color = get_pure_blue_color(het_val)
                else:
                    color = get_pure_red_color(het_val)
                alpha = get_continuous_alpha(het_val)
            
            # Draw arc
            theta = np.linspace(deg1_rad, deg2_rad, 100)
            linewidth = 2.0 if len(dat_processed) <= 100 else (1.5 if len(dat_processed) <= 200 else 1.0)
            ax.plot(theta, [radius]*len(theta), color=color, linewidth=linewidth, alpha=alpha)
        
        # LEGENDS IN SEPARATE AREAS OF THE FIGURE
        
        # 1. EVENT COUNT SUMMARY - Top left area
        count_text = f"Del: {del_count}  Dup: {dup_count}\n"
        if blacklist_regions:
            count_text += f"BL-crossing Del: {bl_del_count}\nBL-crossing Dup: {bl_dup_count}"
        
        fig.text(0.02, 0.85, count_text, fontsize=11, weight='bold',
                bbox=dict(boxstyle="round,pad=0.5", facecolor='white', alpha=0.9, edgecolor='gray'),
                verticalalignment='top')
        
        # 2. GGPLOT2-STYLE GRADIENT LEGEND - Left middle area  
        legend_x = 0.05
        legend_y = 0.4
        legend_height = 0.25
        legend_width = 0.03
        
        # Create continuous gradient bars
        n_steps = 100
        step_height = legend_height / n_steps

        bar_gap = 0.017  # small space between Del and Dup bars
        bar_width = (legend_width - bar_gap) / 2
        
        for i in range(n_steps):
            norm = i / (n_steps - 1)
            het_val = min_het + norm * (max_het - min_het)
            y_pos = legend_y - legend_height/2 + i * step_height

            # Del gradient bar (pure blue, left)
            del_color = get_pure_blue_color(het_val)
            del_alpha = get_continuous_alpha(het_val)
            del_rect = plt.Rectangle(
                (legend_x, y_pos), bar_width, step_height,
                facecolor=del_color, alpha=del_alpha, edgecolor='none',
                transform=fig.transFigure
            )
            fig.patches.append(del_rect)

            # Dup gradient bar (pure red, right, shifted by gap)
            dup_color = get_pure_red_color(het_val)
            dup_alpha = get_continuous_alpha(het_val)
            dup_rect = plt.Rectangle(
                (legend_x + bar_width + bar_gap, y_pos), bar_width, step_height,
                facecolor=dup_color, alpha=dup_alpha, edgecolor='none',
                transform=fig.transFigure
            )
            fig.patches.append(dup_rect)
    
        # Del label centered above left bar
        del_label_x = legend_x + bar_width / 2
        fig.text(del_label_x,
                 legend_y + legend_height/2 + 0.01,
                 "Del", fontsize=8, ha='center', weight='bold', color='blue')

        # Dup label centered above right bar
        dup_label_x = legend_x + bar_width + bar_gap + bar_width / 2
        fig.text(dup_label_x,
                 legend_y + legend_height/2 + 0.01,
                 "Dup", fontsize=8, ha='center', weight='bold', color='red')
        
        # Heteroplasmy value labels
        for i, (pos, val) in enumerate([(1, max_het), (0.5, min_het + 0.5*(max_het-min_het)), (0, min_het)]):
            y_pos = legend_y - legend_height/2 + pos * legend_height
            fig.text(legend_x + legend_width + 0.01, y_pos, f"{val:.3f}", 
                    fontsize=9, va='center')
        
        fig.text(legend_x - 0.05, legend_y, "Heteroplasmy", 
                fontsize=11, weight='bold', va='center', rotation=90)
        
        # 3. Blacklist crossing legend (if needed) - Left bottom
        # Compute top of Del/Dup bars including their labels
        bars_top_y = legend_y + legend_height/2        # top of bars
        label_offset = 0.04                            # offset used for Del/Dup labels
        bl_gap = 0.02                                  # extra gap between bars+labels and rectangle

        # New position for Blacklist rectangle
        bl_x = 0.02
        bl_y = bars_top_y + label_offset + bl_gap      # just above bars/labels
        bl_width = 0.1
        bl_height = 0.05

        if blacklist_regions and (bl_del_count > 0 or bl_dup_count > 0):
            bl_rect = plt.Rectangle(
                (bl_x, bl_y), bl_width, bl_height,
                facecolor=(0.2, 0.8, 0.2), alpha=0.6,
                transform=fig.transFigure
            )
            fig.patches.append(bl_rect)

            # Label inside rectangle (centered)
            fig.text(
                bl_x + bl_width/2, bl_y + bl_height/2,
                "Blacklist\ncrossing",
                fontsize=9, ha='center', va='center', weight='bold'
            )
        
        ax.grid(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        title = 'Mitochondrial DNA Structural Alterations'
        if blacklist_regions:
            title += f' (BL: {len(blacklist_regions)} regions)'
        
        fig.suptitle(title, fontsize=15, weight='bold', y=0.95)
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved to {output_file}")
        print(f"Plotted {len(dat_processed)} events with pure color gradients in separate legend areas")

    def save_results(self, events, output_file, genome_fasta=None, blacklist_file=None):
        """Save processed results following R script logic exactly"""
        if len(events) == 0:
            print("No events to save")
            return
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        res = events.copy()
        print(f"Starting save_results with {len(res)} events")
        
        # Coordinate swapping for OUTPUT only (exactly like R script lines ~360-375)
        for i in res.index:
            if res.loc[i, 'dloop'] == 'yes':
                # Swap coordinates for display
                del_start_median = res.loc[i, 'del.start.median']
                del_end_median = res.loc[i, 'del.end.median']
                del_start_range = res.loc[i, 'del.start.range']
                del_end_range = res.loc[i, 'del.end.range']
                del_start = res.loc[i, 'del.start']
                del_end = res.loc[i, 'del.end']
                lfstart = res.loc[i, 'lfstart']
                lfend = res.loc[i, 'lfend']
                
                res.loc[i, 'del.start.median'] = del_end_median
                res.loc[i, 'del.end.median'] = del_start_median
                res.loc[i, 'del.start.range'] = del_end_range
                res.loc[i, 'del.end.range'] = del_start_range
                res.loc[i, 'del.start'] = del_end
                res.loc[i, 'del.end'] = del_start
                res.loc[i, 'lfstart'] = lfend
                res.loc[i, 'lfend'] = lfstart
        
        # Format exactly like R script
        res['perc'] = res['perc'].round(4)
        
        # Calculate final coordinates exactly like R script
        res['del.start.median'] = res['del.start.median'] + 1
        res['final.event.size'] = np.where(
            res['final.event'] == 'del',
            res['delsize'],
            self.genome_length - res['delsize']
        )
        res['final.end'] = np.where(
            res['final.event'] == 'del',
            res['del.end.median'],
            res['del.start.median'] - 1
        )
        res['final.start'] = np.where(
            res['final.event'] == 'del',
            res['del.start.median'],
            res['del.end.median'] + 1
        )
        
        # Handle wraparound
        res['del.start.median'] = np.where(
            res['del.start.median'] == self.genome_length + 1,
            1,
            res['del.start.median']
        )
        res['final.start'] = np.where(
            res['final.start'] == self.genome_length + 1,
            1,
            res['final.start']
        )

        # --- Blacklist crossing flag using final coordinates ---
        # Load blacklist regions if provided
        try:
            blacklist_regions = self._load_blacklist_regions(blacklist_file)
        except Exception as e:
            print(f"Warning: Could not load blacklist file: {e}")

        def crosses_blacklist(start_pos, end_pos):
            if not blacklist_regions:
                return False
            for region in blacklist_regions:
                bl_start, bl_end = int(region['start']), int(region['end'])
                if (bl_start <= start_pos <= bl_end) or (bl_start <= end_pos <= bl_end):
                    return True
            return False

        res['blacklist_crossing'] = [
            'yes' if crosses_blacklist(row['final.start'], row['final.end']) else 'no'
            for _, row in res.iterrows()
        ]
    
        
        # Get flanking sequences
        if genome_fasta and Path(genome_fasta).exists():
            try:
                genome_record = next(SeqIO.parse(genome_fasta, 'fasta'))
                genome_seq = str(genome_record.seq).upper()
                flanking_results = self._align_breakpoints(
                    res['final.start'] - 1,
                    res['final.end'],
                    genome_seq,
                    self.flank_size
                )
            except Exception as e:
                print(f"Warning: Flanking sequence analysis failed: {e}")
                flanking_results = pd.DataFrame({
                    'seq1': ['NA'] * len(res),
                    'seq2': ['NA'] * len(res),
                    'seq': ['NA'] * len(res)
                })
        else:
            flanking_results = pd.DataFrame({
                'seq1': ['NA'] * len(res),
                'seq2': ['NA'] * len(res),
                'seq': ['NA'] * len(res)
            })
        
        # Create final output exactly like R script
        res_final = pd.DataFrame({
            'sample': res['sample'],
            'cluster.id': res['cluster'],
            'alt.reads': res['nread'].astype(int),
            'ref.reads': res['tread'].astype(int),
            'heteroplasmy': res['perc'],
            'del.start.range': res['del.start.range'],
            'del.end.range': res['del.end.range'],
            'del.size': res['delsize'].astype(int),
            'final.event': res['final.event'],
            'final.start': res['final.start'].astype(int),
            'final.end': res['final.end'].astype(int),
            'final.size': res['final.event.size'].astype(int),
            'blacklist_crossing': res['blacklist_crossing'],
            'seq1': flanking_results['seq1'],
            'seq2': flanking_results['seq2'],
            'seq': flanking_results['seq']
        })
        
        # Save results
        res_final.to_csv(output_file, sep='\t', index=False)
        print(f"Results saved to {output_file}")
        print(f"Events: {len(res_final)}")


    def save_summary(self, events, output_file, analysis_stats):
        """Save analysis summary with biologically meaningful metrics based on MitoSAlt literature"""
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Classify event pattern using literature-based approach
        if len(events) > 0:
            classification, reason, criteria = self._classify_event_pattern(events)
        else:
            classification, reason = "No events", "No events detected"
            criteria = {}
        
        # Count events by type
        del_count = (events['final.event'] == 'del').sum() if len(events) > 0 else 0
        dup_count = (events['final.event'] == 'dup').sum() if len(events) > 0 else 0
        dloop_count = (events['dloop'] == 'yes').sum() if len(events) > 0 else 0
        
        # Calculate comprehensive statistics
        if len(events) > 0:
            size_stats = {
                'min': events['delsize'].min(),
                'max': events['delsize'].max(),
                'mean': events['delsize'].mean(),
                'median': events['delsize'].median()
            }
            het_stats = {
                'min': events['perc'].min(),
                'max': events['perc'].max(),
                'mean': events['perc'].mean(),
                'median': events['perc'].median()
            }
        else:
            size_stats = {'min': 0, 'max': 0, 'mean': 0, 'median': 0}
            het_stats = {'min': 0, 'max': 0, 'mean': 0, 'median': 0}
        
        # Write summary file with literature-based analysis
        with open(output_file, 'w') as f:
            f.write("MitoSAlt Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            
            # Analysis workflow stats
            f.write("Analysis Workflow:\n")
            f.write("-" * 20 + "\n")
            for key, value in analysis_stats.items():
                f.write(f"{key}: {value}\n")
            f.write("\n")
            
            # Event classification (PRIMARY SECTION - most important)
            f.write("Event Pattern Classification:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Type: {classification}\n")
            f.write(f"Biological basis: {reason}\n")
            if 'subtype' in criteria:
                f.write(f"Subtype: {criteria['subtype']}\n")
            f.write("\n")
            
            # Heteroplasmy-based analysis (KEY BIOLOGICAL METRIC)
            f.write("Heteroplasmy Distribution (Literature-based thresholds):\n")
            f.write("-" * 55 + "\n")
            if len(events) > 0:
                f.write(f"High heteroplasmy events (≥30%): {criteria.get('high_het_count', 0)}\n")
                f.write(f"Medium heteroplasmy events (1-30%): {criteria.get('medium_het_count', 0)}\n")
                f.write(f"Low heteroplasmy events (<1%): {criteria.get('low_het_count', 0)}\n")
                f.write(f"Maximum heteroplasmy: {criteria.get('max_heteroplasmy', 0):.3f} ({criteria.get('max_heteroplasmy', 0)*100:.1f}%)\n")
                f.write(f"Median heteroplasmy: {criteria.get('median_heteroplasmy', 0):.3f} ({criteria.get('median_heteroplasmy', 0)*100:.1f}%)\n")
            else:
                f.write("No events detected\n")
            f.write("\n")
            
            # Event counts summary
            f.write("Event Summary:\n")
            f.write("-" * 15 + "\n")
            f.write(f"Total events: {len(events)}\n")
            f.write(f"Deletions: {del_count}\n")
            f.write(f"Duplications: {dup_count}\n")
            f.write(f"Events crossing origin (dloop=yes): {dloop_count}\n\n")
            
            # Basic genomic metrics (for reference)
            f.write("Basic Metrics:\n")
            f.write("-" * 15 + "\n")
            if len(events) > 0:
                f.write(f"Position range spanned: {criteria.get('position_range', 0):.0f} bp\n")
                f.write(f"Size coefficient of variation: {criteria.get('size_coefficient_variation', 0):.2f}\n")
            else:
                f.write("No events to analyze\n")
            f.write("\n")
            
            # Size statistics
            f.write("Size Statistics (bp):\n")
            f.write("-" * 20 + "\n")
            f.write(f"Min size: {size_stats['min']:.0f}\n")
            f.write(f"Max size: {size_stats['max']:.0f}\n")
            f.write(f"Mean size: {size_stats['mean']:.1f}\n")
            f.write(f"Median size: {size_stats['median']:.1f}\n")
            if len(events) > 0:
                f.write(f"Size coefficient of variation: {criteria.get('size_coefficient_variation', 0):.2f}\n")
            f.write("\n")
            
            # Heteroplasmy statistics (expanded)
            f.write("Heteroplasmy Statistics:\n")
            f.write("-" * 25 + "\n")
            f.write(f"Min heteroplasmy: {het_stats['min']:.4f} ({het_stats['min']*100:.2f}%)\n")
            f.write(f"Max heteroplasmy: {het_stats['max']:.4f} ({het_stats['max']*100:.2f}%)\n")
            f.write(f"Mean heteroplasmy: {het_stats['mean']:.4f} ({het_stats['mean']*100:.2f}%)\n")
            f.write(f"Median heteroplasmy: {het_stats['median']:.4f} ({het_stats['median']*100:.2f}%)\n\n")
            
            # Classification details (DETAILED TECHNICAL SECTION)
            if criteria and 'classification_scores' in criteria:
                f.write("Classification Algorithm Details:\n")
                f.write("-" * 35 + "\n")
                scores = criteria['classification_scores']
                f.write(f"Single pattern score: {scores['single_score']}\n")
                f.write(f"Multiple pattern score: {scores['multiple_score']}\n")
                f.write(f"High heteroplasmy threshold: ≥30%\n")
                f.write(f"Low heteroplasmy threshold: <1%\n")
                f.write(f"Multiple event threshold: >20 events\n")
                f.write(f"Major event threshold: ≤3 high-het events\n")
                f.write("\n")
            
            # Biological interpretation
            f.write("Biological Interpretation:\n")
            f.write("-" * 25 + "\n")
            if classification == "Single":
                f.write("Pattern consistent with:\n")
                f.write("- Single pathogenic deletion/duplication\n")
                f.write("- Classical mitochondrial disease patient profile\n")
                f.write("- Possible accompanying low-level artifacts\n")
                if criteria.get('max_heteroplasmy', 0) >= 0.30:
                    f.write("- High heteroplasmy suggests functional impact\n")
            elif classification == "Multiple":
                f.write("Pattern consistent with:\n")
                f.write("- Multiple structural alterations\n") 
                f.write("- Mouse model of mtDNA maintenance defect\n")
                f.write("- Complex replication/repair dysfunction\n")
                if criteria.get('many_events', False):
                    f.write("- Extensive genomic instability\n")
            else:
                f.write("- Ambiguous or unusual pattern\n")
            f.write("\n")
            
            # Literature references
            f.write("Reference Standards:\n")
            f.write("-" * 20 + "\n")
            f.write("Classification based on:\n")
            f.write("- Basu et al. PLoS Genet 2020 (MitoSAlt methodology)\n")
            f.write("- Patient samples: single high-heteroplasmy events (>35%)\n")
            f.write("- Mouse models: multiple low-heteroplasmy events (<3%)\n")
            f.write("- Clinical thresholds: >30% for pathogenic significance\n")
        
        print(f"Enhanced analysis summary saved to {output_file}")
        print(f"Event classification: {classification} ({reason})")
        if 'subtype' in criteria:
            print(f"Pattern subtype: {criteria['subtype']}")


def main():
    parser = argparse.ArgumentParser(description='Generate circular plots of mitochondrial DNA structural alterations')
    parser.add_argument('genome_length', type=int, help='Mitochondrial genome length')
    parser.add_argument('ori_h_start', type=int, help='Heavy strand origin start')
    parser.add_argument('ori_h_end', type=int, help='Heavy strand origin end') 
    parser.add_argument('ori_l_start', type=int, help='Light strand origin start')
    parser.add_argument('ori_l_end', type=int, help='Light strand origin end')
    parser.add_argument('size_limit', type=int, help='Size limit for events')
    parser.add_argument('cluster_file', help='Cluster file path')
    parser.add_argument('breakpoint_file', help='Breakpoint file path')
    parser.add_argument('output_name', help='Output file prefix')
    parser.add_argument('heteroplasmy_limit', type=float, help='Heteroplasmy threshold')
    parser.add_argument('genome_fasta', help='Genome FASTA file')
    parser.add_argument('flank_size', type=int, help='Flanking sequence size')
    parser.add_argument('--blacklist', help='BED file with regions to exclude', default=None)
    parser.add_argument('--output-dir', help='Output directory (default: current directory)', default='.')
    
    args = parser.parse_args()
    
    # Initialize plotter
    plotter = MitoPlotter(
        genome_length=args.genome_length,
        ori_h_start=args.ori_h_start, 
        ori_h_end=args.ori_h_end,
        ori_l_start=args.ori_l_start,
        ori_l_end=args.ori_l_end,
        heteroplasmy_limit=args.heteroplasmy_limit,
        flank_size=args.flank_size
    )
    
    # Load and process data
    result = plotter.load_data(args.cluster_file, args.breakpoint_file, args.blacklist)
    if isinstance(result, tuple):
        events, stats = result
    else:
        events = result
        stats = {}
    
    if len(events) > 0:
        # Create output directories
        output_dir = Path(args.output_dir)
        plot_dir = output_dir / "plot"
        indel_dir = output_dir / "indel"
        
        plot_dir.mkdir(parents=True, exist_ok=True)
        indel_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Processing {len(events)} events")
        
        # Create plot
        plot_file = plot_dir / f"{args.output_name}.png"
        plotter.create_plot(events, str(plot_file), figsize=(10, 10), blacklist_file=args.blacklist)
        
        # Save results
        results_file = indel_dir / f"{args.output_name}.tsv"
        plotter.save_results(events, str(results_file), args.genome_fasta, args.blacklist)
        
        # Save summary
        summary_file = indel_dir / f"{args.output_name}_summary.txt"
        print(f"Attempting to save summary to: {summary_file}")
        print(f"Stats available: {len(stats) if stats else 'None'}")
        plotter.save_summary(events, str(summary_file), stats)
    else:
        print("No events to process")


if __name__ == "__main__":
    main()
