"""
Event Classifier

Naive rule-based classification of mitochondrial event patterns as Single or Multiple.
"""

import pandas as pd
import numpy as np
from .config import ClassificationConfig
from .spatial import SpatialGroupAnalyzer
from .utils import crosses_blacklist


class EventClassifier:
    """Classifies overall pattern of mitochondrial events"""
    
    def __init__(self, genome_length, config=None):
        self.genome_length = genome_length
        self.config = config or ClassificationConfig()
        self.spatial_analyzer = SpatialGroupAnalyzer(genome_length, config)
    
    def classify(self, events: pd.DataFrame, blacklist_regions=None):
        """
        Classify events as Single or Multiple pattern
        
        Single Pattern (patient-like):
        - Few high-heteroplasmy events (≥20%)
        - Spatially clustered (≥70% in one group)
        
        Multiple Pattern (mouse model-like):
        - Many events (>10)
        - Spatially dispersed
        - No dominant high-het event
        
        Args:
            events: DataFrame with mitochondrial events
            blacklist_regions: Optional list of blacklist regions to exclude
            
        Returns:
            Tuple of (classification, reason, criteria, events_with_groups)
        """
        if len(events) == 0:
            return "Unknown", "No events detected", {}, pd.DataFrame()
        
        # Load config
        cfg = self.config
        HIGH_HET = cfg.HIGH_HET_THRESHOLD
        NOISE = cfg.NOISE_THRESHOLD
        CLUSTER_RADIUS = cfg.CLUSTER_RADIUS
        MIN_CLUSTER_SIZE = cfg.MIN_CLUSTER_SIZE
        MULTIPLE_THRESHOLD = cfg.MULTIPLE_EVENT_THRESHOLD
        DOMINANT_FRACTION = cfg.DOMINANT_GROUP_FRACTION
        
        # Filter blacklist-crossing events
        if blacklist_regions:
            clean_events = events[~events.apply(
                lambda row: crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions),
                axis=1
            )]
            blacklist_filtered = len(events) - len(clean_events)
            if len(clean_events) == 0:
                return "Unknown", "All events cross blacklisted regions", {}, events.copy()
        else:
            clean_events = events
            blacklist_filtered = 0
        
        # Basic event metrics
        total_events = len(clean_events)
        high_het_events = clean_events[clean_events['perc'] >= HIGH_HET]
        significant_events = clean_events[clean_events['perc'] >= NOISE]
        
        max_het = clean_events['perc'].max()
        median_het = clean_events['perc'].median()
        
        # Count by type
        del_count = (significant_events['final.event'] == 'del').sum()
        dup_count = (significant_events['final.event'] == 'dup').sum()
        significant_count = len(significant_events)
        
        # Check for no significant events
        if significant_count == 0:
            classification = "No significant events"
            reason = f"only noise-level events (<{NOISE:.0f}%): {total_events} below threshold"
            if blacklist_filtered > 0:
                reason += f" [excluded {blacklist_filtered} blacklist-crossing]"
            
            criteria = {
                'total_events': total_events,
                'significant_count': 0,
                'max_heteroplasmy': max_het,
                'blacklist_filtered': blacklist_filtered
            }
            
            events_with_groups = clean_events.copy()
            events_with_groups['group'] = 'G1'
            return classification, reason, criteria, events_with_groups
        
        # Spatial grouping analysis
        grouping_results = self.spatial_analyzer.group_events(
            clean_events,
            CLUSTER_RADIUS,
            HIGH_HET,
            NOISE,
            MIN_CLUSTER_SIZE
        )
        
        events_with_groups = grouping_results['events_with_groups']
        group_analysis = grouping_results['group_analysis']
        dominant_group = grouping_results['dominant_group']
        dominant_group_count = grouping_results['dominant_group_events']
        
        # Calculate dominant group fraction
        dominant_fraction = dominant_group_count / total_events if total_events > 0 else 0
        
        # ============================================================
        # SIMPLIFIED CLASSIFICATION LOGIC
        # ============================================================
        
        # SINGLE PATTERN indicators:
        # 1. Few high-het events
        # 2. Dominant spatial group (≥70%)
        # 3. Low total count
        has_high_het = len(high_het_events) > 0
        few_events = total_events <= MULTIPLE_THRESHOLD
        dominant_group_pattern = dominant_fraction >= DOMINANT_FRACTION
        
        # MULTIPLE PATTERN indicators:
        # 1. Many events (>threshold)
        # 2. Dispersed (no dominant group)
        # 3. No high-het events
        many_events = total_events > MULTIPLE_THRESHOLD
        dispersed = dominant_fraction < DOMINANT_FRACTION
        no_high_het = len(high_het_events) == 0
        
        # Decision logic
        if has_high_het and (few_events or dominant_group_pattern):
            classification = "Single"
            reasons = []
            
            if len(high_het_events) <= 3:
                reasons.append(f"{len(high_het_events)} high-het event(s) (≥{HIGH_HET:.0f}%)")
            reasons.append(f"max heteroplasmy {max_het:.1f}%")
            
            if dominant_group_pattern:
                reasons.append(f"dominant spatial group ({dominant_group_count}/{total_events}, {dominant_fraction*100:.0f}%)")
            
            if significant_count <= 5:
                reasons.append(f"few significant events ({del_count} del, {dup_count} dup)")
            
            reason = "; ".join(reasons)
            
        elif many_events or (dispersed and no_high_het):
            classification = "Multiple"
            reasons = []
            
            if many_events:
                reasons.append(f"{total_events} events (>{MULTIPLE_THRESHOLD})")
            
            if no_high_het:
                reasons.append(f"no high-het events (<{HIGH_HET:.0f}%)")
            else:
                reasons.append(f"{len(high_het_events)} high-het, median {median_het:.1f}%")
            
            if dispersed:
                reasons.append(f"dispersed pattern ({len(group_analysis)} groups)")
            
            if del_count > 0 and dup_count > 0:
                reasons.append(f"mixed types ({del_count} del, {dup_count} dup)")
            
            reason = "; ".join(reasons)
            
        else:
            # Ambiguous - use max heteroplasmy as tiebreaker
            if max_het >= HIGH_HET:
                classification = "Single"
                reason = f"ambiguous, but high max heteroplasmy ({max_het:.1f}%)"
            else:
                classification = "Multiple"
                reason = f"ambiguous, multiple low-het events (median {median_het:.1f}%)"
        
        # Add blacklist info
        if blacklist_filtered > 0:
            reason += f" [excluded {blacklist_filtered} blacklist-crossing]"
        
        # Build criteria dictionary
        criteria = {
            'total_events': total_events,
            'total_raw_events': len(events),
            'blacklist_filtered': blacklist_filtered,
            'significant_count': significant_count,
            'high_het_count': len(high_het_events),
            'max_heteroplasmy': max_het,
            'median_heteroplasmy': median_het,
            'del_count': del_count,
            'dup_count': dup_count,
            'total_groups': len(group_analysis),
            'dominant_group_fraction': dominant_fraction,
            'dominant_group_count': dominant_group_count,
            'group_analysis': group_analysis
        }
        
        # Pattern subtype
        if classification == "Single":
            if len(high_het_events) == 1:
                criteria['subtype'] = "Classic single deletion/duplication"
            else:
                criteria['subtype'] = "Few major events"
        else:
            if total_events > 50:
                criteria['subtype'] = "Complex multiple (mouse model-like)"
            elif del_count > 0 and dup_count > 0:
                criteria['subtype'] = "Mixed deletion-duplication"
            else:
                criteria['subtype'] = "Multiple single-type events"
        
        # Add blacklist events back for visualization
        if blacklist_regions and blacklist_filtered > 0:
            blacklist_events = events[events.apply(
                lambda row: crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions),
                axis=1
            )].copy()
            
            # Give each blacklist event a unique group ID
            for idx, (_, event) in enumerate(blacklist_events.iterrows()):
                blacklist_events.loc[event.name, 'group'] = f'BL{idx+1}'
            
            events_with_groups = pd.concat([events_with_groups, blacklist_events], ignore_index=True)
        
        return classification, reason, criteria, events_with_groups