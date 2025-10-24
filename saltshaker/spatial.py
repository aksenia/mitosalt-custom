"""
Spatial group analyzer

Analyzes spatial clustering and grouping of mitochondrial events.
Groups events by proximity on the circular genome.
"""

from __future__ import annotations
from typing import List, Dict, Any, Optional
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
    
    def __init__(
        self,
        genome_length: int,
        config: Optional[ClassificationConfig] = None
    ) -> None:
        """
        Initialize SpatialGroupAnalyzer
        
        Args:
            genome_length: Mitochondrial genome length
            config: ClassificationConfig instance (uses defaults if None)
        """
        self.genome_length: int = genome_length
        self.config: ClassificationConfig = config or ClassificationConfig()
    
    def group_events(
        self,
        events_df: pd.DataFrame,
        radius: Optional[int] = None,
        high_het_threshold: Optional[float] = None, 
        sig_het_threshold: Optional[float] = None,
        min_group_size: Optional[int] = None
    ) -> Dict[str, Any]:
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
        
        if events_df.empty:
            return self._empty_grouping_result(events_df)

        logger.debug(f"Spatial Grouping: {len(events_df)} events, radius={radius}bp, min_size={min_group_size}")

        # Separate by event type first to avoid mixed groups
        del_events = events_df[events_df['final_event'] == 'del'].copy()
        dup_events = events_df[events_df['final_event'] == 'dup'].copy()

        logger.debug(f"{len(del_events)} deletions, {len(dup_events)} duplications")
        
        # Group each type separately
        del_groups: List[Dict[str, Any]] = self._group_events_by_type(del_events, radius, 'del') if not del_events.empty else []
        dup_groups: List[Dict[str, Any]] = self._group_events_by_type(dup_events, radius, 'dup') if not dup_events.empty else []
        
        # Combine and assign global group IDs
        all_groups: List[Dict[str, Any]] = del_groups + dup_groups
        
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
    
    def _circular_distance(self, pos1: float, pos2: float) -> float:
        """
        Calculate minimum distance on circular genome
        
        Args:
            pos1: Position 1
            pos2: Position 2
            
        Returns:
            Minimum circular distance
        """
        direct: float = abs(pos1 - pos2)
        wraparound: float = self.genome_length - direct
        return min(direct, wraparound)

    def _events_are_close(
        self,
        event1: Dict[str, Any],
        event2: Dict[str, Any],
        radius: int
    ) -> bool:
        """
        Check if two events are within grouping radius
        
        Args:
            event1: First event dict
            event2: Second event dict
            radius: Maximum distance for grouping
            
        Returns:
            True if events are close enough to group
        """
        start_dist: float = self._circular_distance(event1['start'], event2['start'])
        end_dist: float = self._circular_distance(event1['end'], event2['end'])
        center_dist: float = self._circular_distance(event1['center'], event2['center'])
        
        return min(start_dist, end_dist, center_dist) <= radius

    def _event_to_dict(self, event_row: pd.Series, idx: int) -> Dict[str, Any]:
        """
        Convert event row to dictionary format for grouping
        
        Args:
            event_row: DataFrame row with event data
            idx: Index of event
            
        Returns:
            Dictionary with event information
        """
        return {
            'idx': idx,
            'start': float(event_row['del_start_median']),
            'end': float(event_row['del_end_median']),
            'center': (float(event_row['del_start_median']) + float(event_row['del_end_median'])) / 2,
            'heteroplasmy': float(event_row['perc']),
            'event_type': str(event_row['final_event']),
            'size': float(event_row['delsize'])
        }

    def _group_events_by_type(
        self,
        events: pd.DataFrame,
        radius: int,
        event_type: str
    ) -> List[Dict[str, Any]]:
        """
        Group events of same type using circular distance
        
        Args:
            events: DataFrame with events of single type
            radius: Grouping radius in bp
            event_type: Type of event ('del' or 'dup')
            
        Returns:
            List of group info dictionaries
        """
        if events.empty:
            return []
        
        groups: List[Dict[str, Any]] = []
        used_indices: set = set()
        
        # Sort by heteroplasmy (highest first) - KEEP ORIGINAL INDEX
        events_sorted = events.sort_values('perc', ascending=False)
        
        for idx in events_sorted.index:  # Use original index directly
            if idx in used_indices:
                continue
            
            event = events_sorted.loc[idx]
            
            # Start new group
            group: List[Dict[str, Any]] = [self._event_to_dict(event, idx)]
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

    def _build_group_info(
        self,
        group: List[Dict[str, Any]],
        event_type: str
    ) -> Dict[str, Any]:
        """
        Build group information dictionary
        
        Args:
            group: List of event dicts in this group
            event_type: Type of events ('del' or 'dup')
            
        Returns:
            Dictionary with group statistics and metadata
        """
        cfg = self.config
        
        heteroplasmy_values: List[float] = [e['heteroplasmy'] for e in group]
        sizes: List[float] = [e['end'] - e['start'] for e in group]
        positions: List[float] = []
        for e in group:
            positions.extend([e['start'], e['end']])
        
        group_info: Dict[str, Any] = {
            'id': 0,
            'group_id': 'G1',
            'event_count': len(group),
            'event_type': event_type,
            'max_heteroplasmy': max(heteroplasmy_values),
            'mean_heteroplasmy': float(np.mean(heteroplasmy_values)),
            'total_heteroplasmy': sum(heteroplasmy_values),
            'median_size': float(np.median(sizes)),
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

    def _build_grouping_results(
        self,
        all_groups: List[Dict[str, Any]],
        events_df: pd.DataFrame,
        significant_groups: List[Dict[str, Any]],
        high_het_groups: List[Dict[str, Any]],
        high_het_threshold: float,
        sig_het_threshold: float
    ) -> Dict[str, Any]:
        """
        Build final grouping results with defensive index handling
        
        Args:
            all_groups: All spatial groups
            events_df: Original events DataFrame
            significant_groups: Groups with significant events
            high_het_groups: Groups with high heteroplasmy
            high_het_threshold: High heteroplasmy threshold
            sig_het_threshold: Significance threshold
            
        Returns:
            Dictionary with complete grouping analysis
        """
        dominant_group: Optional[Dict[str, Any]] = all_groups[0] if all_groups else None
        dominant_group_events: int = dominant_group['event_count'] if dominant_group else 0
        dominant_group_range: float = dominant_group['spatial_range'] if dominant_group else 0.0
        
        events_in_significant_groups: int = sum(g['event_count'] for g in significant_groups)
        outlier_events: int = len(events_df) - events_in_significant_groups
        
        # Assign group IDs to events - DEFENSIVE APPROACH
        events_with_groups = events_df.copy()
        events_with_groups['group'] = 'UNGROUPED'  # Default for debugging
        
        if len(all_groups) > 0:
            # Build mapping: original_index -> group_id
            idx_to_group: Dict[int, str] = {}
            for group_info in all_groups:
                for event in group_info['events']:
                    idx_to_group[event['idx']] = group_info['group_id']
            
            # Verify all events are accounted for
            assigned_count: int = 0
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

    def _empty_grouping_result(self, events_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Return empty grouping result
        
        Args:
            events_df: Empty events DataFrame
            
        Returns:
            Empty grouping result dictionary
        """
        return {
            'group_analysis': [],
            'significant_groups': [],
            'high_het_groups': [],
            'dominant_group': None,
            'dominant_group_events': 0,
            'dominant_group_range': 0.0,
            'outlier_events': 0,
            'events_with_groups': events_df.copy()
        }