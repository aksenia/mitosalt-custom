"""
Spatial Group Analyzer

Analyzes spatial clustering and grouping of mitochondrial events.
Groups events by proximity on the circular genome.
"""

import pandas as pd
import numpy as np
import logging
from .config import ClassificationConfig

logger = logging.getLogger(__name__)


class SpatialGroupAnalyzer:
    """
    Analyzes spatial clustering of mitochondrial events
    
    Groups events that occur within close proximity on the circular genome,
    calculates clustering metrics, and identifies dominant spatial patterns.
    """
    
    def __init__(self, genome_length, config=None):
        """
        Initialize SpatialGroupAnalyzer
        
        Args:
            genome_length: Mitochondrial genome length
            config: ClassificationConfig instance (uses defaults if None)
        """
        self.genome_length = genome_length
        self.config = config or ClassificationConfig()
    
    def group_events(self, events_df, radius=None, high_het_threshold=None, 
                    sig_het_threshold=None, min_group_size=None):
        """
        Perform spatial grouping of events within radius distance
        
        Args:
            events_df: DataFrame with events
            radius: Grouping radius in bp (uses config default if None)
            high_het_threshold: High heteroplasmy threshold (uses config if None)
            sig_het_threshold: Significance threshold (uses config if None)
            min_group_size: Minimum cluster size (uses config if None)
            
        Returns:
            Dictionary with grouping results including group_analysis, 
            significant_groups, high_het_groups, etc.
        """
        # Use config defaults if not provided
        if radius is None:
            radius = self.config.CLUSTER_RADIUS
        if high_het_threshold is None:
            high_het_threshold = self.config.HIGH_HET_THRESHOLD
        if sig_het_threshold is None:
            sig_het_threshold = self.config.NOISE_THRESHOLD
        if min_group_size is None:
            min_group_size = self.config.MIN_CLUSTER_SIZE
        
        if len(events_df) == 0:
            return self._empty_grouping_result(events_df)

        logger.debug(f"Spatial Grouping: {len(events_df)} events, radius={radius}bp, min_size={min_group_size}")


        # Separate by event type first to avoid mixed groups
        del_events = events_df[events_df['final.event'] == 'del'].copy()
        dup_events = events_df[events_df['final.event'] == 'dup'].copy()

        logger.debug(f"{len(del_events)} deletions, {len(dup_events)} duplications")
        
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
    
        return self._build_grouping_results(
            all_groups, events_df, 
            significant_groups=all_groups, 
            high_het_groups=[g for g in all_groups if g['high_het_count'] > 0],
            high_het_threshold=high_het_threshold, 
            sig_het_threshold=sig_het_threshold
        )
    
    def _circular_distance(self, pos1, pos2):
        """Calculate minimum distance on circular genome"""
        direct = abs(pos1 - pos2)
        wraparound = self.genome_length - direct
        return min(direct, wraparound)

    def _events_are_close(self, event1, event2, radius):
        """Check if two events are within grouping radius"""
        start_dist = self._circular_distance(event1['start'], event2['start'])
        end_dist = self._circular_distance(event1['end'], event2['end'])
        center_dist = self._circular_distance(event1['center'], event2['center'])
        
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
        """Group events of same type using circular distance"""
        if len(events) == 0:
            return []
        
        groups = []
        used_indices = set()
        
        # Sort by heteroplasmy (highest first) - KEEP ORIGINAL INDEX
        events_sorted = events.sort_values('perc', ascending=False)
        
        for idx in events_sorted.index:  # Use original index directly
            if idx in used_indices:
                continue
            
            event = events_sorted.loc[idx]
            
            # Start new group
            group = [self._event_to_dict(event, idx)]
            used_indices.add(idx)
            
            # Find close events (single-linkage clustering)
            for other_idx in events_sorted.index:
                if other_idx in used_indices:
                    continue
                
                other_event = events_sorted.loc[other_idx]
                other_dict = self._event_to_dict(other_event, other_idx)
                
                # Check if close to ANY event in group (true single-linkage)
                if any(self._events_are_close(group_event, other_dict, radius) 
                    for group_event in group):
                    group.append(other_dict)
                    used_indices.add(other_idx)
            
            # Always create group (even single events)
            groups.append(self._build_group_info(group, event_type))
        
        return groups

    def _build_group_info(self, group, event_type):
        """Build group information dictionary"""
        cfg = self.config
        
        heteroplasmy_values = [e['heteroplasmy'] for e in group]
        sizes = [e['end'] - e['start'] for e in group]
        positions = []
        for e in group:
            positions.extend([e['start'], e['end']])
        
        group_info = {
            'id': 0,
            'group_id': 'G1',
            'event_count': len(group),
            'event_type': event_type,
            'max_heteroplasmy': max(heteroplasmy_values),
            'mean_heteroplasmy': np.mean(heteroplasmy_values),
            'total_heteroplasmy': sum(heteroplasmy_values),
            'median_size': np.median(sizes),  
            'max_size': max(sizes),          
            'spatial_range': max(positions) - min(positions) if len(positions) > 1 else 0,
            'high_het_count': sum(1 for h in heteroplasmy_values if h >= cfg.HIGH_HET_THRESHOLD),
            'significant_count': sum(1 for h in heteroplasmy_values if h >= cfg.NOISE_THRESHOLD),
            'events': group
        }
        
        group_info['dominance_score'] = group_info['max_heteroplasmy'] * group_info['event_count']
        
        # Find representative event
        representative = max(group, key=lambda x: x['heteroplasmy'])
        group_info['representative'] = {
            'heteroplasmy': representative['heteroplasmy'],
            'event_type': representative['event_type'],
            'start': representative['start'],
            'end': representative['end'],
            'size': representative['size']
        }
        
        return group_info

    def _build_grouping_results(self, all_groups, events_df, significant_groups, 
                            high_het_groups, high_het_threshold, sig_het_threshold):
        """Build final grouping results with defensive index handling"""
        dominant_group = all_groups[0] if all_groups else None
        dominant_group_events = dominant_group['event_count'] if dominant_group else 0
        dominant_group_range = dominant_group['spatial_range'] if dominant_group else 0
        
        events_in_significant_groups = sum(g['event_count'] for g in significant_groups)
        outlier_events = len(events_df) - events_in_significant_groups
        
        # Assign group IDs to events - DEFENSIVE APPROACH
        events_with_groups = events_df.copy()
        events_with_groups['group'] = 'UNGROUPED'  # Default for debugging
        
        if len(all_groups) > 0:
            # Build mapping: original_index -> group_id
            idx_to_group = {}
            for group_info in all_groups:
                for event in group_info['events']:
                    idx_to_group[event['idx']] = group_info['group_id']
            
            # Verify all events are accounted for
            assigned_count = 0
            for idx in events_with_groups.index:
                if idx in idx_to_group:
                    events_with_groups.loc[idx, 'group'] = idx_to_group[idx]
                    assigned_count += 1
                else:
                    # Defensive: Should never happen, but flag it
                    logger.warning(f"Event at index {idx} not assigned to any group!")
                    events_with_groups.loc[idx, 'group'] = 'G999'
            
            # Sanity check
            if assigned_count != len(events_df):
                logger.warning(f"Only {assigned_count}/{len(events_df)} events assigned to groups!")
        else:
            # No groups - assign all to default
            events_with_groups['group'] = 'G1'
        
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