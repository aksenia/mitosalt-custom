#!/usr/bin/env python3
"""
MitoSAlt Circular Plot Generator - enhanced rewrite with spatial grouping analysis
"""

import pandas as pd
import numpy as np
import re
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
    
    def load_data(self, cluster_file, breakpoint_file):
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
        
        # Classification: Start with all as deletions, then apply R script logic
        final_clusters['final.event'] = 'del'
        final_clusters = self._apply_origin_classification(final_clusters)
        
        return final_clusters
    
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

    def _perform_spatial_grouping(self, events_df, radius=1000, high_het_threshold=0.30, sig_het_threshold=0.05, min_group_size=2):  # Increased radius
        """
        Perform spatial grouping of events within radius distance, separating by event type
        Returns grouping results with group IDs assigned to events
        """
        if len(events_df) == 0:
            return self._empty_grouping_result(events_df)

        # Separate by event type first to avoid mixed groups
        del_events = events_df[events_df['final.event'] == 'del'].copy()
        dup_events = events_df[events_df['final.event'] == 'dup'].copy()
        
        # Group each type separately
        del_groups = self._group_events_by_type(del_events, radius, 'del') if len(del_events) > 0 else []
        dup_groups = self._group_events_by_type(dup_events, radius, 'dup') if len(dup_events) > 0 else []
        
        # Combine and assign global group IDs
        all_groups = del_groups + dup_groups
        
        # Sort by median size (largest first) for outside-to-inside plotting
        all_groups.sort(key=lambda x: x['median_size'], reverse=True)

        for i, group in enumerate(all_groups):
            group['id'] = i
            group['group_id'] = f'G{i+1}'
    
        return self._build_grouping_results(all_groups, events_df, significant_groups=all_groups, 
                                        high_het_groups=[g for g in all_groups if g['high_het_count'] > 0],
                                        high_het_threshold=high_het_threshold, sig_het_threshold=sig_het_threshold)

    def _circular_distance(self, pos1, pos2):
        """Calculate minimum distance on circular genome"""
        direct = abs(pos1 - pos2)
        wraparound = self.genome_length - direct
        return min(direct, wraparound)

    def _events_are_close(self, event1, event2, radius):
        """Check if two events are within grouping radius using circular distance"""
        # Check multiple distance metrics
        start_dist = self._circular_distance(event1['start'], event2['start'])
        end_dist = self._circular_distance(event1['end'], event2['end'])
        center_dist = self._circular_distance(event1['center'], event2['center'])
        
        # Events are close if any breakpoint is within radius
        return min(start_dist, end_dist, center_dist) <= radius

    def _event_to_dict(self, event_row, idx):
        """Convert event row to dictionary format for grouping"""
        return {
            'idx': idx,
            'start': event_row['del.start.median'],
            'end': event_row['del.end.median'],
            'center': (event_row['del.start.median'] + event_row['del.end.median']) / 2,
            'heteroplasmy': event_row['perc'],
            'event_type': event_row['final.event'],
            'size': event_row['delsize']
        }

    def _group_events_by_type(self, events, radius, event_type):
        """Group events of same type using improved circular distance"""
        if len(events) == 0:
            return []
        
        groups = []
        used_indices = set()
        
        # Sort by heteroplasmy (highest first) to prioritize dominant events
        events_sorted = events.sort_values('perc', ascending=False).reset_index()
        
        for idx, (_, event) in enumerate(events_sorted.iterrows()):
            orig_idx = event['index']  # Original index before sorting
            if orig_idx in used_indices:
                continue
                
            # Start new group with this event
            group = [self._event_to_dict(event, orig_idx)]
            used_indices.add(orig_idx)
            
            # Find events close to this seed event (single-linkage clustering)
            for jdx, (_, other_event) in enumerate(events_sorted.iterrows()):
                orig_jdx = other_event['index']
                if orig_jdx in used_indices:
                    continue
                    
                other_dict = self._event_to_dict(other_event, orig_jdx)
                
                # Check if close to the seed event (not all members - allows elongated groups)
                if self._events_are_close(group[0], other_dict, radius):
                    group.append(other_dict)
                    used_indices.add(orig_jdx)
            
            # Create group info
            if len(group) >= 1:  # Keep all groups including single events
                groups.append(self._build_group_info(group, event_type))
        
        return groups

    def _build_group_info(self, group, event_type):
        """Build group information dictionary"""
        # Thresholds as same in classify events 
        HIGH_HETEROPLASMY_THRESHOLD = 30.0  # 30% - major pathogenic events
        LOW_HETEROPLASMY_THRESHOLD = 1.0    # 1% - potential artifacts/minor events
        SIGNIFICANT_HETEROPLASMY_THRESHOLD = 1.5  # 5% - only count events above this for type classification

        heteroplasmy_values = [e['heteroplasmy'] for e in group]
        sizes = [e['end'] - e['start'] for e in group]  # Use actual genomic spans
        positions = []
        for e in group:
            positions.extend([e['start'], e['end']])
        
        group_info = {
            'id': 0,  # Will be assigned later
            'group_id': 'G1',  # Will be assigned later
            'event_count': len(group),
            'event_type': event_type,
            'max_heteroplasmy': max(heteroplasmy_values),
            'mean_heteroplasmy': np.mean(heteroplasmy_values),
            'total_heteroplasmy': sum(heteroplasmy_values),
            'median_size': np.median(sizes),  
            'max_size': max(sizes),          
            'spatial_range': max(positions) - min(positions) if len(positions) > 1 else 0,
            'high_het_count': sum(1 for h in heteroplasmy_values if h >=  HIGH_HETEROPLASMY_THRESHOLD),
            'significant_count': sum(1 for h in heteroplasmy_values if h >= SIGNIFICANT_HETEROPLASMY_THRESHOLD),
            'events': group
        }
        
        # Calculate dominance score (heteroplasmy × event count)
        group_info['dominance_score'] = group_info['max_heteroplasmy'] * group_info['event_count']
        
        # Find representative event (highest heteroplasmy)
        representative = max(group, key=lambda x: x['heteroplasmy'])
        group_info['representative'] = {
            'heteroplasmy': representative['heteroplasmy'],
            'event_type': representative['event_type'],
            'start': representative['start'],
            'end': representative['end'],
            'size': representative['size']
        }
        
        return group_info

    def _build_grouping_results(self, all_groups, events_df, significant_groups, high_het_groups, high_het_threshold, sig_het_threshold):
        """Build final grouping results"""
        dominant_group = all_groups[0] if all_groups else None
        dominant_group_events = dominant_group['event_count'] if dominant_group else 0
        dominant_group_range = dominant_group['spatial_range'] if dominant_group else 0
        
        # Count outlier events (events not in significant groups)
        events_in_significant_groups = sum(g['event_count'] for g in significant_groups)
        outlier_events = len(events_df) - events_in_significant_groups
        
        # Assign group IDs to events dataframe
        events_with_groups = events_df.copy()
        events_with_groups['group'] = 'G1'  # Default for single events
        
        if len(all_groups) > 0:
            # Create mapping from event index to group ID
            idx_to_group = {}
            for group_info in all_groups:
                for event in group_info['events']:
                    idx_to_group[event['idx']] = group_info['group_id']
            
            # Assign group IDs to events
            for idx, row in events_with_groups.iterrows():
                if idx in idx_to_group:
                    events_with_groups.loc[idx, 'group'] = idx_to_group[idx]
        
        return {
            'group_analysis': all_groups,
            'significant_groups': significant_groups,
            'high_het_groups': high_het_groups,
            'dominant_group': dominant_group,
            'dominant_group_events': dominant_group_events,
            'dominant_group_range': dominant_group_range,
            'outlier_events': outlier_events,
            'events_with_groups': events_with_groups
        }

    def _empty_grouping_result(self, events_df):
        """Return empty grouping result"""
        return {
            'group_analysis': [],
            'significant_groups': [],
            'high_het_groups': [],
            'dominant_group': None,
            'dominant_group_events': 0,
            'dominant_group_range': 0,
            'outlier_events': 0,
            'events_with_groups': events_df.copy()
        }
    
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
    
    def _crosses_blacklist(self, start_pos, end_pos, blacklist_regions):
        """Check if event breakpoints cross any blacklisted regions"""
        if not blacklist_regions:
            return False
        
        for region in blacklist_regions:
            bl_start, bl_end = int(region['start']), int(region['end'])
            if (bl_start <= start_pos <= bl_end) or (bl_start <= end_pos <= bl_end):
                return True
        return False
    
    def _filter_blacklist(self, clusters, blacklist_regions):
        """Filter out events that overlap with blacklisted regions"""
        mask = ~clusters.apply(
            lambda row: self._crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions),
            axis=1
        )
        filtered_count = len(clusters) - mask.sum()
        if filtered_count > 0:
            print(f"Filtered out {filtered_count} events overlapping blacklisted regions")
        return clusters[mask]


    def _classify_event_pattern(self, events, blacklist_regions=None):
        """
        Classify events as Single or Multiple based on MitoSAlt literature criteria:
        - Single: Dominated by one or few high-heteroplasmy events (typically >30-35%)
        - Multiple: Many low-heteroplasmy events, often representing different mechanisms
        
        Note: Excludes blacklist-crossing events from classification as these are likely artifacts
        """
        if len(events) == 0:
            return "Unknown", "No events detected", {}, pd.DataFrame()
        
        # Filter out blacklist-crossing events for classification
        if blacklist_regions:
            print(f"DEBUG: Classification - checking {len(events)} events against {len(blacklist_regions)} blacklist regions")

            # Filter to non-blacklist events for classification
            clean_events = events[~events.apply(lambda row: self._crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions), axis=1)]

            if len(clean_events) == 0:
                return "Unknown", "All events cross blacklisted regions", {}, events.copy()
            
            # Use clean events for classification
            events_for_classification = clean_events
            blacklist_filtered_count = len(events) - len(clean_events)
        else:
            # No blacklist regions provided
            events_for_classification = events
            blacklist_filtered_count = 0
        
        # Key biological thresholds based on literature
        HIGH_HETEROPLASMY_THRESHOLD = 30.0  # 30% - major pathogenic events
        LOW_HETEROPLASMY_THRESHOLD = 1.0   # 1% - potential artifacts/minor events
        SIGNIFICANT_HETEROPLASMY_THRESHOLD = 1.5  # 5% - only count events above this for type classification
        MAJOR_EVENT_COUNT_THRESHOLD = 3     # Few major events = single pattern
        TOTAL_EVENT_COUNT_THRESHOLD = 20    # Many events = multiple pattern
        
        # Clustering parameters
        CLUSTER_RADIUS = 400  # bp - events within 400bp form a cluster
        MIN_CLUSTER_SIZE = 2  # minimum events to be considered a significant cluster
        
        # Heteroplasmy-based classification (using clean events only)
        high_het_events = events_for_classification[events_for_classification['perc'] >= HIGH_HETEROPLASMY_THRESHOLD]
        medium_het_events = events_for_classification[(events_for_classification['perc'] >= LOW_HETEROPLASMY_THRESHOLD) & 
                                (events_for_classification['perc'] < HIGH_HETEROPLASMY_THRESHOLD)]
        low_het_events = events_for_classification[events_for_classification['perc'] < LOW_HETEROPLASMY_THRESHOLD]
        
        # Count ALL event types (for reporting)
        del_count = (events_for_classification['final.event'] == 'del').sum()
        dup_count = (events_for_classification['final.event'] == 'dup').sum()
        
        # Count SIGNIFICANT event types only (for robust classification)
        significant_dels = (events_for_classification['final.event'] == 'del') & (events_for_classification['perc'] >= SIGNIFICANT_HETEROPLASMY_THRESHOLD)
        significant_dups = (events_for_classification['final.event'] == 'dup') & (events_for_classification['perc'] >= SIGNIFICANT_HETEROPLASMY_THRESHOLD)
        
        significant_del_count = significant_dels.sum()
        significant_dup_count = significant_dups.sum()
        
        # Also track high-heteroplasmy event types for additional robustness
        high_het_dels = (events_for_classification['final.event'] == 'del') & (events_for_classification['perc'] >= HIGH_HETEROPLASMY_THRESHOLD)
        high_het_dups = (events_for_classification['final.event'] == 'dup') & (events_for_classification['perc'] >= HIGH_HETEROPLASMY_THRESHOLD)
        
        high_het_del_count = high_het_dels.sum()
        high_het_dup_count = high_het_dups.sum()
        
        # Perform spatial grouping analysis (separate from classification)
        grouping_results = self._perform_spatial_grouping(events_for_classification, CLUSTER_RADIUS, HIGH_HETEROPLASMY_THRESHOLD, SIGNIFICANT_HETEROPLASMY_THRESHOLD, MIN_CLUSTER_SIZE)
        
        # Extract grouping metrics for classification
        group_analysis = grouping_results['group_analysis']
        significant_groups = grouping_results['significant_groups']
        high_het_groups = grouping_results['high_het_groups']
        dominant_group = grouping_results['dominant_group']
        dominant_group_events = grouping_results['dominant_group_events']
        dominant_group_range = grouping_results['dominant_group_range']
        outlier_events = grouping_results['outlier_events']
        
        # Spatial metrics - enhanced with clustering analysis
        positions = []
        breakpoint_positions = []
        for _, event in events_for_classification.iterrows():
            start_pos = event['del.start.median']
            end_pos = event['del.end.median']
            positions.extend([start_pos, end_pos])
            breakpoint_positions.extend([start_pos, end_pos])
        
        # Calculate basic range (may be useful for debugging/validation)
        if len(positions) > 0:
            positions = np.array(positions)
            position_range = np.max(positions) - np.min(positions) if len(positions) > 1 else 0
        else:
            position_range = 0
        
        # Enhanced spatial clustering analysis
        TIGHT_CLUSTER_THRESHOLD = 100   # bp - breakpoints within 100bp = tight cluster (single event)
        LOOSE_CLUSTER_THRESHOLD = 1000  # bp - breakpoints within 1kb = loose cluster
        SCATTERED_THRESHOLD = 3000      # bp - breakpoints >3kb apart = scattered (multiple events)
        
        if len(breakpoint_positions) > 0:
            bp_array = np.array(breakpoint_positions)
            bp_range = np.max(bp_array) - np.min(bp_array) if len(bp_array) > 1 else 0
            bp_std = np.std(bp_array) if len(bp_array) > 1 else 0
            
            # Clustering classifications
            tight_clustering = bp_range <= TIGHT_CLUSTER_THRESHOLD
            loose_clustering = bp_range <= LOOSE_CLUSTER_THRESHOLD
            scattered_pattern = bp_range >= SCATTERED_THRESHOLD
            
            # Calculate clustering density (events per kb)
            clustering_density = len(events_for_classification) / (bp_range / 1000) if bp_range > 0 else float('inf')
        else:
            bp_range = 0
            bp_std = 0
            tight_clustering = False
            loose_clustering = False
            scattered_pattern = False
            clustering_density = 0
        
        # Size variation
        if len(events_for_classification) > 0:
            size_std = np.std(events_for_classification['delsize'])
            mean_size = np.mean(events_for_classification['delsize'])
            size_cv = size_std / mean_size if mean_size > 0 else 0
        else:
            size_std = 0
            size_cv = 0
        
        # Calculate values needed for both criteria and reasoning (BEFORE using them anywhere)
        max_heteroplasmy = events_for_classification['perc'].max() if len(events_for_classification) > 0 else 0
        median_heteroplasmy = events_for_classification['perc'].median() if len(events_for_classification) > 0 else 0
        mixed_types_significant = significant_del_count > 0 and significant_dup_count > 0 and min(significant_del_count, significant_dup_count) >= 3
        mixed_types_all = del_count > 0 and dup_count > 0 and min(del_count, dup_count) >= 3
        many_events = len(events_for_classification) > TOTAL_EVENT_COUNT_THRESHOLD
        
        # PRE-CLASSIFICATION: Check for no significant events
        if significant_del_count + significant_dup_count == 0:
            if len(events_for_classification) == 0:
                classification = "No events"
                reason_str = "No events detected"
            else:
                classification = "No significant events" 
                reason_str = f"only low-level events (<5% heteroplasmy): {len(events_for_classification)} events below significance threshold"
                if blacklist_filtered_count > 0:
                    reason_str += f" [excluded {blacklist_filtered_count} blacklist-crossing events]"
            
            # Build criteria for no-event case
            criteria = {
                'total_events': len(events_for_classification),
                'total_raw_events': len(events),
                'blacklist_filtered_count': blacklist_filtered_count,
                'significant_del_count': significant_del_count,
                'significant_dup_count': significant_dup_count,
                'max_heteroplasmy': max_heteroplasmy,
                'subtype': "Artifacts/noise only" if len(events_for_classification) > 0 else "No events detected"
            }
            
            # For no-event cases, just return events as-is (no grouping needed)
            events_with_groups = events_for_classification.copy() if len(events_for_classification) > 0 else events.copy()
            if len(events_with_groups) > 0 and 'group' not in events_with_groups.columns:
                events_with_groups['group'] = 'G1'  # Default group for any remaining events
            
            return classification, reason_str, criteria, events_with_groups

        # Decision logic enhanced with proper spatial grouping analysis
        # SINGLE EVENT PATTERN criteria (ENHANCED)
        single_pattern_indicators = [
            # Criterion 1: Dominant high-heteroplasmy event(s)
            len(high_het_events) >= 1 and len(high_het_events) <= MAJOR_EVENT_COUNT_THRESHOLD,
            
            # Criterion 2: High max heteroplasmy suggests pathogenic event
            max_heteroplasmy >= HIGH_HETEROPLASMY_THRESHOLD,
            
            # Criterion 3: Few significant events total (ignores low-het noise)
            significant_del_count + significant_dup_count <= MAJOR_EVENT_COUNT_THRESHOLD and len(high_het_events) >= 1,
            
            # Criterion 4: Single type dominance among SIGNIFICANT events with high heteroplasmy
            (significant_del_count == 0 or significant_dup_count == 0) and max_heteroplasmy >= HIGH_HETEROPLASMY_THRESHOLD,
            
            # Criterion 5: Strong single-type dominance at high heteroplasmy
            (high_het_del_count > 0 and high_het_dup_count == 0) or (high_het_dup_count > 0 and high_het_del_count == 0),
            
            # Criterion 6: Dominant single type even with mixed low-het noise
            (significant_del_count >= 1 and significant_dup_count == 0) or (significant_dup_count >= 1 and significant_del_count == 0),
            
            # ENHANCED: Dominant group with most events (≥70% in main group)
            dominant_group_events >= 0.7 * len(events_for_classification) if dominant_group else False,
            
            # ENHANCED: Single significant group dominates
            len(significant_groups) == 1 and significant_groups[0]['high_het_count'] > 0 if significant_groups else False,
            
            # ENHANCED: Few outliers compared to dominant group
            outlier_events <= dominant_group_events if dominant_group else False,
            
            # Criterion 7: Tight spatial clustering + high heteroplasmy = single underlying event
            tight_clustering and len(high_het_events) >= 1,
            
            # Criterion 8: Loose clustering with few significant events = single
            loose_clustering and significant_del_count + significant_dup_count <= MAJOR_EVENT_COUNT_THRESHOLD,
            
            # Criterion 9: High clustering density suggests same underlying event
            clustering_density > 5.0 and len(high_het_events) >= 1
        ]
        
        # MULTIPLE EVENT PATTERN criteria (ENHANCED)
        multiple_pattern_indicators = [
            # Criterion 1: Many total events (mouse model pattern)
            len(events_for_classification) > TOTAL_EVENT_COUNT_THRESHOLD,
            
            # Criterion 2: Many events but no high-heteroplasmy dominant event
            len(events_for_classification) > 10 and len(high_het_events) == 0,
            
            # Criterion 3: Mixed SIGNIFICANT deletion/duplication pattern
            mixed_types_significant and len(events_for_classification) > 10,
            
            # Criterion 4: High complexity - many medium heteroplasmy events
            len(medium_het_events) > 10 and median_heteroplasmy < HIGH_HETEROPLASMY_THRESHOLD,
            
            # Criterion 5: Very high significant event count suggests multiple mechanisms
            significant_del_count + significant_dup_count > 20,
            
            # Criterion 6: Multiple high-heteroplasmy events of different types
            high_het_del_count >= 2 and high_het_dup_count >= 2,
            
            # ENHANCED: Multiple significant groups (≥3)
            len(significant_groups) >= 3,
            
            # ENHANCED: Multiple high-heteroplasmy groups (≥2)
            len(high_het_groups) >= 2,
            
            # ENHANCED: No dominant group (scattered pattern)
            not dominant_group or dominant_group_events < 0.5 * len(events_for_classification),
            
            # ENHANCED: Many outliers vs grouped events
            outlier_events > dominant_group_events if dominant_group else len(events_for_classification) > 10,
            
            # Criterion 7: Scattered spatial pattern suggests multiple independent events
            scattered_pattern and len(events_for_classification) > 5,
            
            # Criterion 8: Wide spatial distribution with multiple significant events
            bp_range > LOOSE_CLUSTER_THRESHOLD and significant_del_count + significant_dup_count > 5,
            
            # Criterion 9: Low clustering density with many events = multiple mechanisms
            clustering_density < 2.0 and len(events_for_classification) > 15
        ]
        
        # Classification decision
        single_score = sum(single_pattern_indicators)
        multiple_score = sum(multiple_pattern_indicators)
        
        if single_score > multiple_score:
            classification = "Single"
            
            # Detailed reasoning for single pattern (with group-based analysis)
            reasons = []
            if len(high_het_events) >= 1:
                reasons.append(f"dominant high-heteroplasmy event(s) ({len(high_het_events)} at ≥{HIGH_HETEROPLASMY_THRESHOLD:.0%})")
            if max_heteroplasmy >= HIGH_HETEROPLASMY_THRESHOLD:
                reasons.append(f"max heteroplasmy {max_heteroplasmy:.1%}")
            if significant_del_count + significant_dup_count <= MAJOR_EVENT_COUNT_THRESHOLD:
                reasons.append(f"few significant events ({significant_del_count} del, {significant_dup_count} dup ≥5%)")
            
            # Add group-based spatial analysis
            if dominant_group_events > 0:
                group_percentage = (dominant_group_events / len(events_for_classification)) * 100
                reasons.append(f"dominant group ({dominant_group_events}/{len(events_for_classification)} events, {group_percentage:.0f}%)")
                    
                # Only report spatial range for multiple events
                if dominant_group_events > 1:
                    if dominant_group_range <= TIGHT_CLUSTER_THRESHOLD:
                        reasons.append(f"tightly grouped (range: {dominant_group_range:.0f}bp)")
                    elif dominant_group_range <= LOOSE_CLUSTER_THRESHOLD:
                        reasons.append(f"loosely grouped (range: {dominant_group_range:.0f}bp)")
                        
                if outlier_events > 0:
                    reasons.append(f"{outlier_events} outlier events (likely artifacts)")
                        
            if clustering_density > 10.0:
                reasons.append(f"high clustering density ({clustering_density:.1f} events/kb)")
            elif clustering_density > 5.0:
                reasons.append(f"good clustering density ({clustering_density:.1f} events/kb)")
                    
            if len(low_het_events) > 0:
                reasons.append(f"{len(low_het_events)} low-level events (<1%)")
            
            reason_str = "; ".join(reasons) if reasons else f"grouped pattern ({len(events_for_classification)} events)"
            
        elif multiple_score > single_score:
            classification = "Multiple"
            
            # Detailed reasoning for multiple pattern (with group-based analysis)
            reasons = []
            if len(events_for_classification) > TOTAL_EVENT_COUNT_THRESHOLD:
                reasons.append(f"{len(events_for_classification)} events")
            if len(high_het_events) == 0:
                reasons.append("no dominant high-heteroplasmy events")
            if mixed_types_significant:
                reasons.append(f"mixed significant types ({significant_del_count} del, {significant_dup_count} dup ≥5%)")
            if high_het_del_count >= 2 and high_het_dup_count >= 2:
                reasons.append(f"multiple high-het types ({high_het_del_count} del, {high_het_dup_count} dup ≥30%)")
                    
            # Add group-based spatial analysis
            if len(significant_groups) > 1:
                reasons.append(f"{len(significant_groups)} significant groups")
            if len(high_het_groups) > 1:
                reasons.append(f"{len(high_het_groups)} high-heteroplasmy groups")
            if scattered_pattern and outlier_events > dominant_group_events:
                reasons.append(f"scattered pattern (outliers > dominant group)")
            elif bp_range > LOOSE_CLUSTER_THRESHOLD:
                reasons.append(f"wide distribution ({bp_range:.0f}bp span)")
                    
            if clustering_density < 2.0 and len(events_for_classification) > 15:
                reasons.append(f"low clustering density ({clustering_density:.1f} events/kb)")
            
            reason_str = "; ".join(reasons)
            
        else:
            # Ambiguous case - use conservative approach
            if max_heteroplasmy >= HIGH_HETEROPLASMY_THRESHOLD:
                classification = "Single"
                reason_str = f"ambiguous but high max heteroplasmy ({max_heteroplasmy:.1%})"
            else:
                classification = "Multiple" 
                reason_str = f"ambiguous, multiple low-heteroplasmy events (median {median_heteroplasmy:.2%})"
        
        # Add blacklist filtering information to reason if applicable
        if blacklist_filtered_count > 0:
            reason_str += f" [excluded {blacklist_filtered_count} blacklist-crossing events]"
        
        # Build criteria dictionary AFTER all calculations are done
        criteria = {
            # Core biological metrics from MitoSAlt literature (based on clean events)
            'total_events': len(events_for_classification),
            'total_raw_events': len(events),
            'blacklist_filtered_count': blacklist_filtered_count,
            'high_het_count': len(high_het_events),
            'medium_het_count': len(medium_het_events), 
            'low_het_count': len(low_het_events),
            'max_heteroplasmy': max_heteroplasmy,
            'median_heteroplasmy': median_heteroplasmy,
            
            # Event type counts
            'del_count': del_count,
            'dup_count': dup_count,
            'significant_del_count': significant_del_count,
            'significant_dup_count': significant_dup_count,
            'high_het_del_count': high_het_del_count,
            'high_het_dup_count': high_het_dup_count,
            
            # Spatial grouping analysis metrics (ENHANCED)
            'total_groups': len(group_analysis),
            'significant_groups': len(significant_groups),
            'high_het_groups': len(high_het_groups),
            'dominant_group_events': dominant_group_events,
            'dominant_group_range': dominant_group_range,
            'outlier_events': outlier_events,
            'group_analysis': group_analysis,
            
            # Spatial clustering metrics (legacy)
            'breakpoint_range': bp_range,
            'breakpoint_std': bp_std,
            'tight_clustering': tight_clustering,
            'loose_clustering': loose_clustering,
            'scattered_pattern': scattered_pattern,
            'clustering_density': clustering_density,
            
            # Secondary metrics
            'size_coefficient_variation': size_cv,
            'position_range': position_range,
            'size_std': size_std,
            'many_events': many_events,
            'mixed_types_all': mixed_types_all,
            'mixed_types_significant': mixed_types_significant
        }
        
        # Add biological pattern classification scores
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
            if len(events_for_classification) > 100:
                criteria['subtype'] = "Complex multiple (mouse model-like)"
            elif del_count > 0 and dup_count > 0:
                criteria['subtype'] = "Mixed deletion-duplication pattern"
            else:
                criteria['subtype'] = "Multiple single-type events"
        
        # Return events with group assignments
        events_with_groups = grouping_results['events_with_groups']

        # NOW ADD BLACKLIST EVENTS BACK for visualization
        if blacklist_regions:
            # Find blacklist-crossing events
            blacklist_events = events[events.apply(lambda row: self._crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions), axis=1)].copy()
            
            if len(blacklist_events) > 0:
                # Add blacklist flag and unique group info for each event
                for idx, (_, event) in enumerate(blacklist_events.iterrows()):
                    blacklist_events.loc[event.name, 'group'] = f'BL{idx+1}' 
                blacklist_events['blacklist_crossing'] = True
                
                # Add back to events_with_groups
                events_with_groups = pd.concat([events_with_groups, blacklist_events], ignore_index=True)
                print(f"DEBUG: Added {len(blacklist_events)} blacklist events back for plotting")

        return classification, reason_str, criteria, events_with_groups

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
    
    def group_sort_key(self, group_id):
        """Convert group ID to sortable tuple (priority, number)"""
        match = re.match(r'^([A-Z]+)(\d+)$', group_id)
        if match:
            prefix, number = match.groups()
            if prefix == 'G':
                return (0, int(number))  # Regular groups first
            elif prefix == 'BL':
                return (1000, int(number))  # Blacklist groups after regular groups
            else:
                # Unexpected format - log warning and put at end
                print(f"WARNING: Unexpected group ID format: {group_id}")
                return (9999, int(number))
        else:
            print(f"WARNING: Could not parse group ID: {group_id}")
            return (9999, 0)
    
    def create_plot(self, events, output_file, figsize=(16, 10), blacklist_regions=None):
        """Create circular plot of mitochondrial events"""
        if len(events) == 0:
            print("No events to plot")
            return
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Create data structure for plotting
        dat = pd.DataFrame({
            'chr': 'MT',
            'start': events['del.start.median'],
            'end': events['del.end.median'],
            'value': events['perc'],
            'dloop': events['dloop'],
            'delsize': events['delsize'],
            'final.event': events['final.event'],
            'group': events.get('group', 'G1')
        })
        
        # Order events by group for better visualization
        if 'group' in events.columns:
            dat = dat.sort_values(['group', 'value'], ascending=[True, False]).reset_index(drop=True)
        
        # Enhanced blacklist crossing detection      
        dat['blacklist_crossing'] = False
        if blacklist_regions:
            for idx, row in dat.iterrows():
                dat.loc[idx, 'blacklist_crossing'] = self._crosses_blacklist(row['start'], row['end'], blacklist_regions)

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
        
        # NESTED FUNCTIONS
        def assign_radii_by_type(data, base_radius=380, radius_diff=8):
            """Assign radii to events of a single type within their allocated radius range"""
            if len(data) == 0:
                return data
            
            data = data.sort_values(['group', 'deg1'], ascending=[True, True]).reset_index(drop=True)
            data['radius'] = 0
            
            def events_overlap(event1, event2):
                start1, end1 = event1['deg1'] % 360, event1['deg2'] % 360
                start2, end2 = event2['deg1'] % 360, event2['deg2'] % 360
                min_gap = 4
                
                def normalize_arc(start, end):
                    return [(start, 360), (0, end)] if start > end else [(start, end)]
                
                arcs1, arcs2 = normalize_arc(start1, end1), normalize_arc(start2, end2)
                for arc1_start, arc1_end in arcs1:
                    for arc2_start, arc2_end in arcs2:
                        if not (arc1_end + min_gap <= arc2_start or arc2_end + min_gap <= arc1_start):
                            return True
                return False
            
            unique_groups = data['group'].unique()
            group_counts = data['group'].value_counts().to_dict()
            
            single_event_groups = [g for g in unique_groups if group_counts[g] == 1]
            multi_event_groups = [g for g in unique_groups if group_counts[g] > 1]
            
            group_band_size = 18
            group_gap = 6
            current_radius = base_radius
            assignments = {}
            
            # Assign dedicated bands to multi-event groups
            for group_id in sorted(multi_event_groups, key=self.group_sort_key):
                band_top = current_radius
                band_bottom = band_top - group_band_size
                assignments[group_id] = {'band_top': band_top, 'band_bottom': band_bottom, 'shared': False}
                current_radius = band_bottom - group_gap
            
            # Pack single-event groups on shared radii  
            shared_levels = []
            for group_id in sorted(single_event_groups, key=self.group_sort_key):
                event_deg = data[data['group'] == group_id].iloc[0]['deg1']
                
                # Try to share with existing level
                placed = False
                for level in shared_levels:
                    if all(abs(event_deg - deg) >= 90 and abs(event_deg - deg) <= 270 
                        for deg in level['degrees']):
                        level['degrees'].append(event_deg)
                        assignments[group_id] = {'radius': level['radius'], 'shared': True}
                        placed = True
                        break
                
                if not placed:
                    new_radius = current_radius - group_gap
                    shared_levels.append({'radius': new_radius, 'degrees': [event_deg]})
                    assignments[group_id] = {'radius': new_radius, 'shared': True}
                    current_radius = new_radius - group_gap
            
            # Assign actual radii to events
            for group_id, assignment in assignments.items():
                group_indices = data[data['group'] == group_id].index.tolist()
                
                if assignment['shared']:
                    for idx in group_indices:
                        data.loc[idx, 'radius'] = assignment['radius']
                else:
                    band_top, band_bottom = assignment['band_top'], assignment['band_bottom']
                    for i, idx in enumerate(group_indices):
                        test_radius = band_top
                        while test_radius >= band_bottom:
                            if not any(data.loc[prev_idx, 'radius'] == test_radius and 
                                    events_overlap(data.loc[idx], data.loc[prev_idx])
                                    for prev_idx in group_indices[:i]):
                                data.loc[idx, 'radius'] = test_radius
                                break
                            test_radius -= 2
                        else:
                            data.loc[idx, 'radius'] = band_bottom
            
            return data

        def calculate_dynamic_radius_layout(dat_del, dat_dup, base_radius=400, separator_frac=0.15):
            """Calculate dynamic radius layout with proportional space allocation"""
            
            def calculate_space_needed(data):
                if len(data) == 0:
                    return 0
                groups = data['group'].unique()
                group_counts = data['group'].value_counts().to_dict()
                
                multi_event_groups = len([g for g in groups if group_counts[g] > 1])
                single_events = len([g for g in groups if group_counts[g] == 1])
                shared_levels_needed = max(1, (single_events + 3) // 4)
                
                return multi_event_groups + shared_levels_needed
            
            del_space_needed = calculate_space_needed(dat_del)
            dup_space_needed = calculate_space_needed(dat_dup)
            
            total_group_space = del_space_needed + dup_space_needed
            
            if total_group_space == 0:
                return dat_del, dat_dup, base_radius * separator_frac, base_radius
            
            available_frac = 1.0 - separator_frac
            del_frac = (del_space_needed / total_group_space) * available_frac
            dup_frac = (dup_space_needed / total_group_space) * available_frac
            
            # Determine outer vs inner based on group numbers
            del_min_group = int(min(dat_del['group'].str[1:]).replace('', '999')) if len(dat_del) > 0 else 999
            dup_min_group = int(min(dat_dup['group'].str[1:]).replace('', '999')) if len(dat_dup) > 0 else 999
            
            if dup_min_group < del_min_group:
                outer_frac, inner_frac = dup_frac, del_frac
                outer_data, inner_data = dat_dup, dat_del
                outer_type = 'dup'
            else:
                outer_frac, inner_frac = del_frac, dup_frac  
                outer_data, inner_data = dat_del, dat_dup
                outer_type = 'del'
            
            # Calculate radius ranges
            inner_max = base_radius * inner_frac
            outer_min = base_radius * (inner_frac + separator_frac)
            blacklist_radius = (inner_max + outer_min) / 2
            
            print(f"DEBUG: Fractions - Inner: {inner_frac:.3f}, Separator: {separator_frac:.3f}, Outer: {outer_frac:.3f}")
            print(f"DEBUG: Ranges - Inner: [0-{inner_max:.1f}], Outer: [{outer_min:.1f}-{base_radius}], BL: {blacklist_radius:.1f}")
            
            # Assign radii within ranges
            if len(inner_data) > 0:
                inner_data = assign_radii_by_type(inner_data, base_radius=inner_max, radius_diff=6)
            if len(outer_data) > 0:
                outer_data = assign_radii_by_type(outer_data, base_radius=base_radius, radius_diff=6)
            
            # Return in correct order
            if outer_type == 'dup':
                return inner_data, outer_data, blacklist_radius, base_radius
            else:
                return outer_data, inner_data, blacklist_radius, base_radius

        # MAIN PROCESSING
        # Separate data by type
        dat_del = dat[dat['final.event'] == 'del'].copy()
        dat_dup = dat[dat['final.event'] == 'dup'].copy()

        # Process duplication delsize
        if len(dat_dup) > 0:
            dat_dup['delsize'] = self.genome_length - dat_dup['delsize']

        # Calculate dynamic layout
        dat_del, dat_dup, blacklist_radius, dynamic_radius = calculate_dynamic_radius_layout(dat_del, dat_dup, base_radius=400)

        # Combine for processing
        dat_processed = pd.concat([dat_del, dat_dup], ignore_index=True) if len(dat_dup) > 0 else dat_del
        
        # Calculate color scales
        del_events = dat_processed[dat_processed['final.event'] == 'del']
        dup_events = dat_processed[dat_processed['final.event'] == 'dup']
        del_max = del_events['value'].max() if len(del_events) > 0 else 0
        del_min = del_events['value'].min() if len(del_events) > 0 else 0
        dup_max = dup_events['value'].max() if len(dup_events) > 0 else 0  
        dup_min = dup_events['value'].min() if len(dup_events) > 0 else 0
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='polar', position=[0.15, 0.05, 0.7, 0.9])
        ax.set_ylim(0, dynamic_radius + 30)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        
        # Draw genome circle (use dynamic radius instead of hardcoded 390)
        circle = patches.Circle((0, 0), dynamic_radius, fill=False, linewidth=3,
                            color='gray', transform=ax.transData._b)
        ax.add_patch(circle)
        
        # Add blacklist regions
        if blacklist_regions:
            separator_circle = patches.Circle((0, 0), blacklist_radius + 15, fill=False, linewidth=2, 
                                            color='lightgray', linestyle='--', alpha=0.7,
                                            transform=ax.transData._b)
            ax.add_patch(separator_circle)
            
            # Find actual minimum event radius for inward lines
            min_event_radius = dat_processed['radius'].min() if len(dat_processed) > 0 else 50
            
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
                
                # Radial lines from blacklist to outer circle and inner events
                ax.plot([start_deg, start_deg], [blacklist_radius, dynamic_radius], color='gray', linewidth=1, linestyle='--', alpha=0.6)
                ax.plot([end_deg, end_deg], [blacklist_radius, dynamic_radius], color='gray', linewidth=1, linestyle='--', alpha=0.6)
                
                # Lines extending inward if there are inner events
                if min_event_radius < blacklist_radius:
                    ax.plot([start_deg, start_deg], [min_event_radius, blacklist_radius], color='gray', linewidth=1, linestyle='--', alpha=0.6)
                    ax.plot([end_deg, end_deg], [min_event_radius, blacklist_radius], color='gray', linewidth=1, linestyle='--', alpha=0.6)
        
        # Add position markers
        positions = np.arange(0, self.genome_length, 1000)
        for pos in positions:
            deg = np.radians(358 * pos / self.genome_length)
            ax.plot([deg, deg], [dynamic_radius, dynamic_radius + 5], color='gray', linewidth=1)
            ax.text(deg, dynamic_radius + 15, f'{pos//1000}', ha='center', va='center', 
                fontsize=9, color='gray')
        
        # Color functions
        def get_pure_blue_color(het_val, min_het, max_het):
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5
            blue = 1.0
            red = 0.8 * (1 - norm)
            green = 0.9 * (1 - norm)
            return (red, green, blue)

        def get_pure_red_color(het_val, min_het, max_het):
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5
            red = 1.0
            green = 0.8 * (1 - norm)
            blue = 0.85 * (1 - norm)
            return (red, green, blue)

        def get_continuous_alpha(het_val, min_het, max_het):
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5
            return 0.75 + 0.2 * norm
        
        # Plot events
        for i, (_, event) in enumerate(dat_processed.iterrows()):
            deg1_rad = np.radians(event['deg1'])
            deg2_rad = np.radians(event['deg2'])
            radius = event['radius']
            het_val = event['value']
            
            if blacklist_regions and event['blacklist_crossing']:
                color = (0.2, 0.8, 0.2)
                alpha = 0.9
                print(f"DEBUG: Plotting blacklist event at {event['start']}-{event['end']}")
            else:
                if event['final.event'] == 'del':
                    color = get_pure_blue_color(het_val, del_min, del_max)
                    alpha = get_continuous_alpha(het_val, del_min, del_max)
                else:
                    color = get_pure_red_color(het_val, dup_min, dup_max) 
                    alpha = get_continuous_alpha(het_val, dup_min, dup_max)
                
            theta = np.linspace(deg1_rad, deg2_rad, 100)
            linewidth = 2.0 if len(dat_processed) <= 100 else (1.5 if len(dat_processed) <= 200 else 1.0)
            ax.plot(theta, [radius]*len(theta), color=color, linewidth=linewidth, alpha=alpha)

        # Group labeling
        if len(dat_processed) > 0:
            group_representatives = {}
            for _, event in dat_processed.iterrows():
                group_id = event['group']
                het_val = event['value']
                # pick up the leftmost one for labeling
                if group_id not in group_representatives or event['deg1'] < group_representatives[group_id]['deg']:
                    group_representatives[group_id] = {
                        'deg': event['deg1'],  # Now using leftmost position
                        'radius': event['radius'],
                        'het_val': het_val,
                        'event_type': event['final.event']
                    }
            
            for group_id, info in group_representatives.items():
                breakpoint_deg_rad = np.radians(info['deg'])
                breakpoint_radius = info['radius']
                label_radius = breakpoint_radius + 17
                
                label_color = 'blue' if info['event_type'] == 'del' else 'red'
                
                ax.plot([breakpoint_deg_rad, breakpoint_deg_rad], 
                        [breakpoint_radius + 1.5, label_radius - 4], 
                        color='grey', linewidth=1, alpha=0.7, linestyle='-')
                
                ax.plot(breakpoint_deg_rad, breakpoint_radius, 
                        marker='o', markersize=3, color='grey', alpha=0.8)
                
                ax.text(breakpoint_deg_rad, label_radius, group_id, 
                        ha='center', va='center', fontsize=6, weight='normal',
                        color=label_color,
                        bbox=dict(boxstyle="round,pad=0.15", facecolor='white', 
                                alpha=0.9, edgecolor=label_color, linewidth=0.8))

        # LEGENDS IN SEPARATE AREAS OF THE FIGURE
            
        # 1. EVENT COUNT SUMMARY - Top left area
        count_text = f"Del: {del_count}  Dup: {dup_count}\n"
        if blacklist_regions:
            count_text += f"BL-crossing Del: {bl_del_count}\nBL-crossing Dup: {bl_dup_count}"
        
        fig.text(0.02, 0.85, count_text, fontsize=11, weight='bold',
                bbox=dict(boxstyle="round,pad=0.5", facecolor='white', alpha=0.9, edgecolor='gray'),
                verticalalignment='top')
        
        # 2. SEPARATE GRADIENT LEGENDS - independently scaled
        legend_x = 0.05
        legend_y = 0.4
        legend_height = 0.25
        legend_width = 0.03

        # Create separate gradient bars
        n_steps = 100
        step_height = legend_height / n_steps
        bar_gap = 0.02
        bar_width = (legend_width - bar_gap) / 2

        # Deletion gradient (left bar)
        if len(del_events) > 0:
            for i in range(n_steps):
                norm = i / (n_steps - 1)
                het_val = del_min + norm * (del_max - del_min) if del_max > del_min else del_min
                y_pos = legend_y - legend_height/2 + i * step_height

                del_color = get_pure_blue_color(het_val, del_min, del_max)
                del_alpha = get_continuous_alpha(het_val, del_min, del_max)
                del_rect = plt.Rectangle(
                    (legend_x, y_pos), bar_width, step_height,
                    facecolor=del_color, alpha=del_alpha, edgecolor='none',
                    transform=fig.transFigure
                )
                fig.patches.append(del_rect)

        # Duplication gradient (right bar) 
        if len(dup_events) > 0:
            for i in range(n_steps):
                norm = i / (n_steps - 1)
                het_val = dup_min + norm * (dup_max - dup_min) if dup_max > dup_min else dup_min
                y_pos = legend_y - legend_height/2 + i * step_height

                dup_color = get_pure_red_color(het_val, dup_min, dup_max)
                dup_alpha = get_continuous_alpha(het_val, dup_min, dup_max)
                dup_rect = plt.Rectangle(
                    (legend_x + bar_width + bar_gap, y_pos), bar_width, step_height,
                    facecolor=dup_color, alpha=dup_alpha, edgecolor='none',
                    transform=fig.transFigure
                )
                fig.patches.append(dup_rect)

        # Labels and values
        del_label_x = legend_x + bar_width / 2
        fig.text(del_label_x, legend_y + legend_height/2 + 0.01, "Del", 
                fontsize=8, ha='center', weight='bold', color='blue')

        dup_label_x = legend_x + bar_width + bar_gap + bar_width / 2
        fig.text(dup_label_x, legend_y + legend_height/2 + 0.01, "Dup", 
                fontsize=8, ha='center', weight='bold', color='red')

        # Separate value labels for each scale
        if len(del_events) > 0:
            for pos, val in [(1, del_max), (0, del_min)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height
                fig.text(legend_x - 0.02, y_pos, f"{val:.1f}%", 
                        fontsize=8, va='center', ha='right', color='blue')

        if len(dup_events) > 0:
            for pos, val in [(1, dup_max), (0, dup_min)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height  
                fig.text(legend_x + legend_width + 0.01, y_pos, f"{val:.1f}%", 
                        fontsize=8, va='center', color='red')

        fig.text(legend_x + legend_width/2, legend_y - legend_height/2 - 0.03, 
                "Heteroplasmy (%)", fontsize=11, weight='bold', ha='center')
        
        # 3. Blacklist crossing legend (if needed)
        bars_top_y = legend_y + legend_height/2        
        label_offset = 0.04                            
        bl_gap = 0.02                                  

        bl_x = 0.02
        bl_y = bars_top_y + label_offset + bl_gap      
        bl_width = 0.1
        bl_height = 0.05

        if blacklist_regions and (bl_del_count > 0 or bl_dup_count > 0):
            bl_rect = plt.Rectangle(
                (bl_x, bl_y), bl_width, bl_height,
                facecolor=(0.2, 0.8, 0.2), alpha=0.6,
                transform=fig.transFigure
            )
            fig.patches.append(bl_rect)

            fig.text(
                bl_x + bl_width/2, bl_y + bl_height/2,
                "Blacklist\ncrossing",
                fontsize=9, ha='center', va='center', weight='bold'
            )
    
        ax.grid(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        # LEGENDS IN SEPARATE AREAS OF THE FIGURE
        
        # 1. EVENT COUNT SUMMARY - Top left area
        count_text = f"Del: {del_count}  Dup: {dup_count}\n"
        if blacklist_regions:
            count_text += f"BL-crossing Del: {bl_del_count}\nBL-crossing Dup: {bl_dup_count}"
        
        fig.text(0.02, 0.85, count_text, fontsize=11, weight='bold',
                bbox=dict(boxstyle="round,pad=0.5", facecolor='white', alpha=0.9, edgecolor='gray'),
                verticalalignment='top')
        
        # 2. SEPARATE GRADIENT LEGENDS - independently scaled
        legend_x = 0.05
        legend_y = 0.4
        legend_height = 0.25
        legend_width = 0.03

        # Create separate gradient bars
        n_steps = 100
        step_height = legend_height / n_steps
        bar_gap = 0.02
        bar_width = (legend_width - bar_gap) / 2

        # Deletion gradient (left bar)
        if len(del_events) > 0:
            for i in range(n_steps):
                norm = i / (n_steps - 1)
                het_val = del_min + norm * (del_max - del_min) if del_max > del_min else del_min
                y_pos = legend_y - legend_height/2 + i * step_height

                del_color = get_pure_blue_color(het_val, del_min, del_max)
                del_alpha = get_continuous_alpha(het_val, del_min, del_max)
                del_rect = plt.Rectangle(
                    (legend_x, y_pos), bar_width, step_height,
                    facecolor=del_color, alpha=del_alpha, edgecolor='none',
                    transform=fig.transFigure
                )
                fig.patches.append(del_rect)

        # Duplication gradient (right bar) 
        if len(dup_events) > 0:
            for i in range(n_steps):
                norm = i / (n_steps - 1)
                het_val = dup_min + norm * (dup_max - dup_min) if dup_max > dup_min else dup_min
                y_pos = legend_y - legend_height/2 + i * step_height

                dup_color = get_pure_red_color(het_val, dup_min, dup_max)
                dup_alpha = get_continuous_alpha(het_val, dup_min, dup_max)
                dup_rect = plt.Rectangle(
                    (legend_x + bar_width + bar_gap, y_pos), bar_width, step_height,
                    facecolor=dup_color, alpha=dup_alpha, edgecolor='none',
                    transform=fig.transFigure
                )
                fig.patches.append(dup_rect)

        # Labels and values
        del_label_x = legend_x + bar_width / 2
        fig.text(del_label_x, legend_y + legend_height/2 + 0.01, "Del", 
                fontsize=8, ha='center', weight='bold', color='blue')

        dup_label_x = legend_x + bar_width + bar_gap + bar_width / 2
        fig.text(dup_label_x, legend_y + legend_height/2 + 0.01, "Dup", 
                fontsize=8, ha='center', weight='bold', color='red')

        # Separate value labels for each scale
        if len(del_events) > 0:
            for pos, val in [(1, del_max), (0, del_min)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height
                fig.text(legend_x - 0.02, y_pos, f"{val:.1f}%", 
                        fontsize=8, va='center', ha='right', color='blue')

        if len(dup_events) > 0:
            for pos, val in [(1, dup_max), (0, dup_min)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height  
                fig.text(legend_x + legend_width + 0.01, y_pos, f"{val:.1f}%", 
                        fontsize=8, va='center', color='red')

        fig.text(legend_x + legend_width/2, legend_y - legend_height/2 - 0.03, 
                "Heteroplasmy (%)", fontsize=11, weight='bold', ha='center')
        
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

    def save_results(self, events, output_file, genome_fasta=None, blacklist_regions=None):
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
        res['blacklist_crossing'] = [
            'yes' if self._crosses_blacklist(row['final.start'], row['final.end'], blacklist_regions) else 'no'
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
        
        # Create final output exactly like R script PLUS group information
        res_final = pd.DataFrame({
            'sample': res['sample'],
            'cluster.id': res['cluster'],
            'group': res.get('group', 'G1'),  # Add group column with default
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

    def save_summary(self, events, output_file, analysis_stats, blacklist_regions=None):
        """Save analysis summary with biologically meaningful metrics based on MitoSAlt literature"""
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Classify event pattern using literature-based approach (excluding blacklist events)
        if len(events) > 0:
            classification, reason, criteria, events_with_groups = self._classify_event_pattern(events, blacklist_regions)
        else:
            classification, reason = "No events", "No events detected"
            criteria = {}
            events_with_groups = pd.DataFrame()
        
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
            
            # Event pattern classification (PRIMARY SECTION - most important)
            f.write("Event pattern classification (blacklist-filtered):\n")
            f.write("-" * 50 + "\n")
            f.write(f"Type: {classification}\n")
            f.write(f"Biological basis: {reason}\n")
            if 'subtype' in criteria:
                f.write(f"Subtype: {criteria['subtype']}\n")
            if criteria.get('blacklist_filtered_count', 0) > 0:
                f.write(f"Blacklist-crossing events excluded: {criteria['blacklist_filtered_count']}\n")
                f.write(f"Events used for classification: {criteria.get('total_events', 0)} (of {criteria.get('total_raw_events', 0)} total)\n")
            f.write("\n")
            
            # Spatial groups analysis (NEW SECTION)
            f.write("Spatial Groups Analysis:\n")
            f.write("-" * 25 + "\n")
            if 'group_analysis' in criteria and criteria['group_analysis']:
                # Group by event type for better reporting
                del_groups = [g for g in criteria['group_analysis'] if g.get('event_type') == 'del']
                dup_groups = [g for g in criteria['group_analysis'] if g.get('event_type') == 'dup']
                
                if del_groups:
                    f.write("Deletion groups:\n")
                    for group_info in del_groups:
                        rep = group_info['representative']
                        actual_size = rep['end'] - rep['start']
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1f} "
                            f"(at {rep['start']:.0f}-{rep['end']:.0f}bp, size {actual_size:.0f}bp)\n")
                
                if dup_groups:
                    f.write("Duplication groups:\n")
                    for group_info in dup_groups:
                        rep = group_info['representative']
                        actual_size = rep['end'] - rep['start']
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1f} "
                            f"(at {rep['start']:.0f}-{rep['end']:.0f}bp, size {actual_size:.0f}bp)\n")
                
                if not del_groups and not dup_groups:
                    f.write("Groups detected but event types not properly classified\n")
                    # Fallback - show all groups regardless of type
                    for group_info in criteria['group_analysis']:
                        rep = group_info['representative']
                        event_type = group_info.get('event_type', rep.get('event_type', 'unknown'))
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} {event_type} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1%}\n")
            else:
                f.write("No spatial groups identified\n")
            f.write("\n")
            
            # Heteroplasmy-based analysis (KEY BIOLOGICAL METRIC)
            f.write("Heteroplasmy distribution (literature-based thresholds):\n")
            f.write("-" * 55 + "\n")
            if len(events) > 0:
                f.write(f"High heteroplasmy events (≥30%): {criteria.get('high_het_count', 0)}\n")
                f.write(f"Medium heteroplasmy events (1-30%): {criteria.get('medium_het_count', 0)}\n")
                f.write(f"Low heteroplasmy events (<1%): {criteria.get('low_het_count', 0)}\n")
                f.write(f"Maximum heteroplasmy: {criteria.get('max_heteroplasmy', 0):.3f}%\n")
                f.write(f"Median heteroplasmy: {criteria.get('median_heteroplasmy', 0):.3f}%\n")
            else:
                f.write("No events detected\n")
            f.write("\n")
            
            # Event counts summary
            f.write("Event summary:\n")
            f.write("-" * 15 + "\n")
            f.write(f"Total events: {len(events)}\n")
            f.write(f"Deletions: {del_count}\n")
            f.write(f"Duplications: {dup_count}\n")
            f.write(f"Events crossing origin (dloop=yes): {dloop_count}\n\n")
            
            # Basic genomic metrics (for reference)
            f.write("Basic metrics:\n")
            f.write("-" * 15 + "\n")
            if len(events) > 0:
                f.write(f"Position range spanned: {criteria.get('position_range', 0):.0f} bp\n")
                f.write(f"Size coefficient of variation: {criteria.get('size_coefficient_variation', 0):.2f}\n")
            else:
                f.write("No events to analyze\n")
            f.write("\n")
            
            # Size statistics
            f.write("Size statistics (bp):\n")
            f.write("-" * 20 + "\n")
            f.write(f"Min size: {size_stats['min']:.0f}\n")
            f.write(f"Max size: {size_stats['max']:.0f}\n")
            f.write(f"Mean size: {size_stats['mean']:.1f}\n")
            f.write(f"Median size: {size_stats['median']:.1f}\n")
            if len(events) > 0:
                f.write(f"Size coefficient of variation: {criteria.get('size_coefficient_variation', 0):.2f}\n")
            f.write("\n")
            
            # Heteroplasmy statistics
            f.write("Heteroplasmy statistics:\n")
            f.write("-" * 25 + "\n")
            f.write(f"Min heteroplasmy: {het_stats['min']:.4f}%\n")
            f.write(f"Max heteroplasmy: {het_stats['max']:.4f}%\n")
            f.write(f"Mean heteroplasmy: {het_stats['mean']:.4f}%\n")
            f.write(f"Median heteroplasmy: {het_stats['median']:.4f}%\n\n")

            # Classification details (DETAILED TECHNICAL SECTION)
            if criteria and 'classification_scores' in criteria:
                f.write("Classification Algorithm Details:\n")
                f.write("-" * 35 + "\n")
                scores = criteria['classification_scores']
                f.write(f"Single pattern score: {scores['single_score']} / 13 possible criteria\n")
                f.write(f"Multiple pattern score: {scores['multiple_score']} / 13 possible criteria\n")
                f.write(f"Decision: {'Single' if scores['single_score'] > scores['multiple_score'] else 'Multiple'} pattern (higher score wins)\n")
                f.write("\n")
                f.write("Score calculation:\n")
                f.write("- Each pattern type has 13 biological criteria\n")
                f.write("- Single pattern criteria: dominant high-het events, few total events, tight clustering, etc.\n")
                f.write("- Multiple pattern criteria: many events, mixed types, scattered distribution, etc.\n")
                f.write("- Score = count of criteria met for each pattern type\n")
                f.write("\n")
                f.write("Biological thresholds used:\n")
                f.write(f"- High heteroplasmy threshold: ≥30% (pathogenic significance)\n")
                f.write(f"- Significance threshold: ≥5% (above noise level)\n")
                f.write(f"- Low heteroplasmy threshold: <1% (likely artifacts)\n")
                f.write(f"- Multiple event threshold: >20 events (mouse model pattern)\n")
                f.write(f"- Major event threshold: ≤3 high-het events (single pattern)\n")
                f.write(f"- Spatial clustering radius: 500bp (biologically relevant)\n")
                f.write("\n")
            
            # Biological interpretation
            f.write("Biological interpretation:\n")
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
            f.write("Reference standards:\n")
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

    # Load blacklist regions
    blacklist_regions = None
    if args.blacklist:
        try:
            blacklist_regions = plotter._load_blacklist_regions(args.blacklist)
            print(f"Loaded {len(blacklist_regions)} blacklist regions")
        except Exception as e:
            print(f"Warning: Could not load blacklist file: {e}")
            blacklist_regions = []
    
    # Load and process data
    result = plotter.load_data(args.cluster_file, args.breakpoint_file)
    if isinstance(result, tuple):
        events, stats = result
    else:
        events = result
        stats = {}
    
    if len(events) > 0:
        # Get events with group assignments from classification
        classification, reason, criteria, events_with_groups = plotter._classify_event_pattern(events, blacklist_regions=blacklist_regions)
        
        # Create output directories
        output_dir = Path(args.output_dir)
        plot_dir = output_dir / "plot"
        indel_dir = output_dir / "indel"
        
        plot_dir.mkdir(parents=True, exist_ok=True)
        indel_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Processing {len(events)} events")
        
        # Use the SAME filtered dataset for all operations to avoid index mismatch
        plot_file = plot_dir / f"{args.output_name}.png"
        plotter.create_plot(events_with_groups, str(plot_file), figsize=(10, 10), blacklist_regions=blacklist_regions)
        
        results_file = indel_dir / f"{args.output_name}.grouped.tsv"
        plotter.save_results(events_with_groups, str(results_file), args.genome_fasta, blacklist_regions=blacklist_regions)
        
        # Use filtered events for summary too to maintain consistency
        summary_file = indel_dir / f"{args.output_name}_summary.txt"
        print(f"Attempting to save summary to: {summary_file}")
        print(f"Stats available: {len(stats) if stats else 'None'}")
        plotter.save_summary(events_with_groups, str(summary_file), stats, blacklist_regions=blacklist_regions)
    else:
        print("No events to process")


if __name__ == "__main__":
    main()
