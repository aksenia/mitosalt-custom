"""
Event Classifier

Classifies mitochondrial event patterns as Single or Multiple.
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
    
    def classify(self, events, blacklist_regions=None):
        """
        Classify events as Single or Multiple based on literature criteria:
        - Single: Dominated by one or few high-heteroplasmy events (typically >30-35%)
        - Multiple: Many low-heteroplasmy events, often representing different mechanisms
        
        Note: Excludes blacklist-crossing events from classification as these are likely artifacts
        """
        if len(events) == 0:
            return "Unknown", "No events detected", {}, pd.DataFrame()
        
        # Filter out blacklist-crossing events for classification
        if blacklist_regions:
            print(f"DEBUG: Classification - checking {len(events)} events against {len(blacklist_regions)} blacklist regions")
            clean_events = events[~events.apply(lambda row: crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions), axis=1)]
            if len(clean_events) == 0:
                return "Unknown", "All events cross blacklisted regions", {}, events.copy()
            events_for_classification = clean_events
            blacklist_filtered_count = len(events) - len(clean_events)
        else:
            events_for_classification = events
            blacklist_filtered_count = 0
        
        # ============================================================
        # LOAD ALL CONFIG CONSTANTS
        # ============================================================
        cfg = self.config
        
        # Heteroplasmy thresholds
        HIGH_HETEROPLASMY_THRESHOLD = cfg.HIGH_HETEROPLASMY_THRESHOLD
        LOW_HETEROPLASMY_THRESHOLD = cfg.LOW_HETEROPLASMY_THRESHOLD
        SIGNIFICANT_HETEROPLASMY_THRESHOLD = cfg.SIGNIFICANT_HETEROPLASMY_THRESHOLD
        
        # Event count thresholds
        MAJOR_EVENT_COUNT_THRESHOLD = cfg.MAJOR_EVENT_COUNT_THRESHOLD
        TOTAL_EVENT_COUNT_THRESHOLD = cfg.TOTAL_EVENT_COUNT_THRESHOLD
        MIN_EVENTS_NO_DOMINANT = cfg.MIN_EVENTS_NO_DOMINANT
        MIN_MIXED_TYPE_COUNT = cfg.MIN_MIXED_TYPE_COUNT
        MIN_MEDIUM_HET_FOR_MULTIPLE = cfg.MIN_MEDIUM_HET_FOR_MULTIPLE
        MIN_HIGH_HET_BOTH_TYPES = cfg.MIN_HIGH_HET_BOTH_TYPES
        MIN_SCATTERED_EVENTS = cfg.MIN_SCATTERED_EVENTS
        MIN_WIDE_SIGNIFICANT = cfg.MIN_WIDE_SIGNIFICANT
        MIN_LOW_DENSITY_EVENTS = cfg.MIN_LOW_DENSITY_EVENTS
        
        # Spatial thresholds
        CLUSTER_RADIUS = cfg.CLUSTER_RADIUS
        MIN_CLUSTER_SIZE = cfg.MIN_CLUSTER_SIZE
        TIGHT_CLUSTER_THRESHOLD = cfg.TIGHT_CLUSTER_THRESHOLD
        LOOSE_CLUSTER_THRESHOLD = cfg.LOOSE_CLUSTER_THRESHOLD
        SCATTERED_THRESHOLD = cfg.SCATTERED_THRESHOLD
        
        # Group thresholds
        DOMINANT_GROUP_THRESHOLD = cfg.DOMINANT_GROUP_THRESHOLD
        NO_DOMINANT_THRESHOLD = cfg.NO_DOMINANT_THRESHOLD
        
        # Density thresholds
        HIGH_CLUSTERING_DENSITY = cfg.HIGH_CLUSTERING_DENSITY
        LOW_CLUSTERING_DENSITY = cfg.LOW_CLUSTERING_DENSITY
        MIN_GROUPS_FOR_MULTIPLE = cfg.MIN_GROUPS_FOR_MULTIPLE
        MIN_HIGH_HET_GROUPS = cfg.MIN_HIGH_HET_GROUPS
        
        # ============================================================
        # CLASSIFICATION LOGIC STARTS HERE
        # ============================================================
        
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
        spatial_analyzer = SpatialGroupAnalyzer(self.genome_length, self.config)
        grouping_results = spatial_analyzer.group_events(
            events_for_classification, 
            CLUSTER_RADIUS, 
            HIGH_HETEROPLASMY_THRESHOLD,
            SIGNIFICANT_HETEROPLASMY_THRESHOLD,
            MIN_CLUSTER_SIZE
)

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
        mixed_types_significant = significant_del_count > 0 and significant_dup_count > 0 and min(significant_del_count, significant_dup_count) >= MIN_MIXED_TYPE_COUNT
        mixed_types_all = del_count > 0 and dup_count > 0 and min(del_count, dup_count) >= MIN_MIXED_TYPE_COUNT
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
                events_with_groups['group'] = 'G1'
            
            return classification, reason_str, criteria, events_with_groups

        # Decision logic enhanced with proper spatial grouping analysis
        # SINGLE EVENT PATTERN criteria
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
            
            # Criterion 7: Dominant group with most events (≥70% in main group)
            dominant_group_events >= DOMINANT_GROUP_THRESHOLD * len(events_for_classification) if dominant_group else False,
            
            # Criterion 8: Single significant group dominates
            len(significant_groups) == 1 and significant_groups[0]['high_het_count'] > 0 if significant_groups else False,
            
            # Criterion 9: Few outliers compared to dominant group
            outlier_events <= dominant_group_events if dominant_group else False,
            
            # Criterion 10: Tight spatial clustering + high heteroplasmy = single underlying event
            tight_clustering and len(high_het_events) >= 1,
            
            # Criterion 11: Loose clustering with few significant events = single
            loose_clustering and significant_del_count + significant_dup_count <= MAJOR_EVENT_COUNT_THRESHOLD,
            
            # Criterion 12: High clustering density suggests same underlying event
            clustering_density > HIGH_CLUSTERING_DENSITY and len(high_het_events) >= 1
        ]
        
        # MULTIPLE EVENT PATTERN criteria
        multiple_pattern_indicators = [
            # Criterion 1: Many total events (mouse model pattern)
            len(events_for_classification) > TOTAL_EVENT_COUNT_THRESHOLD,
            
            # Criterion 2: Many events but no high-heteroplasmy dominant event
            len(events_for_classification) > MIN_EVENTS_NO_DOMINANT and len(high_het_events) == 0,
            
            # Criterion 3: Mixed SIGNIFICANT deletion/duplication pattern
            mixed_types_significant and len(events_for_classification) > MIN_EVENTS_NO_DOMINANT,
            
            # Criterion 4: High complexity - many medium heteroplasmy events
            len(medium_het_events) > MIN_MEDIUM_HET_FOR_MULTIPLE and median_heteroplasmy < HIGH_HETEROPLASMY_THRESHOLD,
            
            # Criterion 5: Very high significant event count suggests multiple mechanisms
            significant_del_count + significant_dup_count > TOTAL_EVENT_COUNT_THRESHOLD,
            
            # Criterion 6: Multiple high-heteroplasmy events of different types
            high_het_del_count >= MIN_HIGH_HET_BOTH_TYPES and high_het_dup_count >= MIN_HIGH_HET_BOTH_TYPES,
            
            # Criterion 7: Multiple significant groups (≥3)
            len(significant_groups) >= MIN_GROUPS_FOR_MULTIPLE,
            
            # Criterion 8: Multiple high-heteroplasmy groups (≥2)
            len(high_het_groups) >= MIN_HIGH_HET_GROUPS,
            
            # Criterion 9: No dominant group (scattered pattern)
            not dominant_group or dominant_group_events < NO_DOMINANT_THRESHOLD * len(events_for_classification),
            
            # Criterion 10: Many outliers vs grouped events
            outlier_events > dominant_group_events if dominant_group else len(events_for_classification) > MIN_EVENTS_NO_DOMINANT,
            
            # Criterion 11: Scattered spatial pattern suggests multiple independent events
            scattered_pattern and len(events_for_classification) > MIN_SCATTERED_EVENTS,
            
            # Criterion 12: Wide spatial distribution with multiple significant events
            bp_range > LOOSE_CLUSTER_THRESHOLD and significant_del_count + significant_dup_count > MIN_WIDE_SIGNIFICANT,
            
            # Criterion 13: Low clustering density with many events = multiple mechanisms
            clustering_density < LOW_CLUSTERING_DENSITY and len(events_for_classification) > MIN_LOW_DENSITY_EVENTS
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
                        reasons.append(f"tightly grouped (≤{TIGHT_CLUSTER_THRESHOLD}bp: {dominant_group_range:.0f}bp)")
                    elif dominant_group_range <= LOOSE_CLUSTER_THRESHOLD:
                        reasons.append(f"loosely grouped (≤{LOOSE_CLUSTER_THRESHOLD}bp: {dominant_group_range:.0f}bp)")
                        
                if outlier_events > 0:
                    reasons.append(f"{outlier_events} outlier events (likely artifacts)")
                    
            if len(low_het_events) > 0:
                reasons.append(f"{len(low_het_events)} low-level events (<{LOW_HETEROPLASMY_THRESHOLD:.0f}%)")
            
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
            if high_het_del_count >= MIN_HIGH_HET_BOTH_TYPES and high_het_dup_count >= MIN_HIGH_HET_BOTH_TYPES:
                reasons.append(f"multiple high-het types ({high_het_del_count} del, {high_het_dup_count} dup ≥30%)")
                    
            # Add group-based spatial analysis
            if len(significant_groups) > 1:
                reasons.append(f"{len(significant_groups)} significant groups")
            if len(high_het_groups) > 1:
                reasons.append(f"{len(high_het_groups)} high-heteroplasmy groups")
            if scattered_pattern and outlier_events > dominant_group_events:
                reasons.append(f"scattered pattern (>{SCATTERED_THRESHOLD}bp: outliers > dominant group)")
            elif bp_range > LOOSE_CLUSTER_THRESHOLD:
                reasons.append(f"wide distribution (>{LOOSE_CLUSTER_THRESHOLD}bp: {bp_range:.0f}bp)")
            
            reason_str = "; ".join(reasons)
            
        else:
            # Ambiguous case - use conservative approach
            if max_heteroplasmy >= HIGH_HETEROPLASMY_THRESHOLD:
                classification = "Single"
                reason_str = f"ambiguous but high max heteroplasmy ({max_heteroplasmy:.1%})"
            else:
                classification = "Multiple" 
                reason_str = f"ambiguous, multiple low-heteroplasmy events (median {median_heteroplasmy:.2%})"

        # Report clustering density for all patterns (outside the if/elif)
        if clustering_density > HIGH_CLUSTERING_DENSITY:
            reason_str += f"; high clustering density ({clustering_density:.1f} events/kb, >{HIGH_CLUSTERING_DENSITY:.0f})"
        elif clustering_density < LOW_CLUSTERING_DENSITY:
            reason_str += f"; low clustering density ({clustering_density:.1f} events/kb, <{LOW_CLUSTERING_DENSITY:.0f})"

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
            
            # Spatial grouping analysis metrics
            'total_groups': len(group_analysis),
            'significant_groups': len(significant_groups),
            'high_het_groups': len(high_het_groups),
            'dominant_group_events': dominant_group_events,
            'dominant_group_range': dominant_group_range,
            'outlier_events': outlier_events,
            'group_analysis': group_analysis,
            
            # Spatial clustering metrics
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
            blacklist_events = events[events.apply(lambda row: crosses_blacklist(row['del.start.median'], row['del.end.median'], blacklist_regions), axis=1)].copy()
            
            if len(blacklist_events) > 0:
                # Add blacklist flag and unique group info for each event
                for idx, (_, event) in enumerate(blacklist_events.iterrows()):
                    blacklist_events.loc[event.name, 'group'] = f'BL{idx+1}' 
                blacklist_events['blacklist_crossing'] = True
                
                # Add back to events_with_groups
                events_with_groups = pd.concat([events_with_groups, blacklist_events], ignore_index=True)
                print(f"DEBUG: Added {len(blacklist_events)} blacklist events back for plotting")

        return classification, reason_str, criteria, events_with_groups