"""
Layout Engine for SaltShaker - SIMPLIFIED VERSION
Pure layout logic with simple radial stacking

Key simplifications:
- No sector-based allocation (causes bunching)
- Simple radial stacking from outer to inner
- Single-pass algorithm
- Group-aware positioning with consistent spacing
- No dynamic expansion or complex budgets
"""
from __future__ import annotations
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np
import logging
from math import ceil, log

from .types import (
    LayoutResult,
    SectorBudget,  # Keep for compatibility but won't use sectors
    GroupBandLayout,
    SingleEventLayout
)

logger = logging.getLogger(__name__)


class LayoutEngine:
    """
    Simplified layout engine with radial stacking
    
    Algorithm:
    1. Determine which type (del/dup) goes outer based on median size
    2. For each type, stack groups from outer to inner radius
    3. Multi-event groups get consistent spacing within bands
    4. Single events share radius levels when possible
    """
    
    def __init__(self, config, genome_length: int = 16569):
        """
        Initialize layout engine
        
        Args:
            config: Plot configuration
            genome_length: Mitochondrial genome length (bp)
        """
        self.config = config
        self.layout_config = config.layout
        self.genome_length = genome_length
        
        logger.info(f"LayoutEngine initialized (simplified radial stacking)")
    
    def calculate_layout(
        self,
        del_events: pd.DataFrame,
        dup_events: pd.DataFrame,
        total_radius: float
    ) -> LayoutResult:
        """
        Calculate layout for all events
        
        Args:
            del_events: DataFrame with deletion events
            dup_events: DataFrame with duplication events  
            total_radius: Total radius available (px)
        
        Returns:
            LayoutResult with all positions calculated
        """
        logger.info(f"Calculating layout for {len(del_events)} dels, {len(dup_events)} dups")
        
        # Ensure events have required columns
        for df in [del_events, dup_events]:
            if not df.empty and 'radius' not in df.columns:
                df['radius'] = 0.0
        
        # Step 1: Determine ring assignment (outer vs inner)
        # Put events with SHORTER arcs on OUTER ring for better visibility
        del_median_size = del_events['delsize'].median() if not del_events.empty else 0
        dup_median_size = dup_events['delsize'].median() if not dup_events.empty else 0
        
        # Reverse logic: smaller arcs go to outer ring
        if dup_median_size < del_median_size:
            outer_type, outer_events = 'dup', dup_events
            inner_type, inner_events = 'del', del_events
        else:
            outer_type, outer_events = 'del', del_events
            inner_type, inner_events = 'dup', dup_events
        
        logger.info(f"Ring assignment: {outer_type} outer ({len(outer_events)} events), "
                   f"{inner_type} inner ({len(inner_events)} events)")
        
        # Step 2: Calculate separator position
        separator_size = total_radius * self.layout_config.separator_fraction
        available_radius = total_radius - separator_size
        
        # Allocate radius proportionally to event counts  
        total_events = len(outer_events) + len(inner_events)
        if total_events == 0:
            return self._create_empty_layout(total_radius)
        
        # Give each ring space proportional to its event count
        outer_fraction = len(outer_events) / total_events if total_events > 0 else 0.5
        inner_fraction = 1.0 - outer_fraction
        
        # Ensure minimum space for each ring if it has events
        min_ring_size = 50  # Minimum pixels per ring
        outer_radius = max(available_radius * outer_fraction, 
                          min_ring_size if len(outer_events) > 0 else 0)
        inner_radius = max(available_radius * inner_fraction,
                          min_ring_size if len(inner_events) > 0 else 0)
        
        # Normalize if total exceeds available
        if outer_radius + inner_radius > available_radius:
            scale = available_radius / (outer_radius + inner_radius)
            outer_radius *= scale
            inner_radius *= scale
        
        # Calculate ring ranges with buffer zones around separator circle
        buffer_zone = 8  # pixels of buffer around separator circle
        outer_range = (separator_size + inner_radius + buffer_zone, total_radius - buffer_zone)
        inner_range = (separator_size + buffer_zone, separator_size + inner_radius - buffer_zone)
        
        # FIX: Blacklist radius should be at the boundary between rings (separator circle)
        # This is where the blacklist regions will be drawn
        blacklist_radius = separator_size + inner_radius  # = outer_range[0] - buffer_zone
        
        logger.info(f"Outer ring ({outer_type}): {outer_range[0]:.1f}-{outer_range[1]:.1f} px (buffered)")
        logger.info(f"Inner ring ({inner_type}): {inner_range[0]:.1f}-{inner_range[1]:.1f} px (buffered)")
        logger.info(f"Blacklist separator radius: {blacklist_radius:.1f} px")
        
        # Step 3: Layout each ring
        outer_layout = self._layout_events_in_ring(outer_events, outer_range, outer_type)
        inner_layout = self._layout_events_in_ring(inner_events, inner_range, inner_type)
        
        # Step 4: Combine results
        all_events = pd.concat([outer_events, inner_events], ignore_index=False)
        
        # Combine group bands and single events
        all_group_bands = outer_layout['group_bands'] + inner_layout['group_bands']
        all_single_events = outer_layout['single_events'] + inner_layout['single_events']
        
        # Determine radius ranges
        if outer_type == 'del':
            del_radius_range = outer_range
            dup_radius_range = inner_range
        else:
            del_radius_range = inner_range
            dup_radius_range = outer_range
        
        # Calculate actual used radius
        total_radius_used = all_events['radius'].max() if not all_events.empty else 0
        
        return LayoutResult(
            events=all_events,
            sector_budgets={},  # No sectors in simplified version
            group_bands=all_group_bands,
            single_events=all_single_events,
            total_radius_used=total_radius_used,
            del_radius_range=del_radius_range,
            dup_radius_range=dup_radius_range,
            blacklist_radius=blacklist_radius,
            layout_algorithm='simplified_radial',
            layout_stats={
                'outer_type': outer_type,
                'inner_type': inner_type,
                'n_group_bands': len(all_group_bands),
                'n_single_events': len(all_single_events)
            }
        )
    
    def _layout_events_in_ring(
        self,
        events: pd.DataFrame,
        radius_range: Tuple[float, float],
        event_type: str
    ) -> Dict:
        """
        Layout all events within a ring using simple stacking
        
        Args:
            events: Events to layout (modified in place)
            radius_range: (min_radius, max_radius) for this ring
            event_type: 'del' or 'dup' for logging
        
        Returns:
            Dict with 'group_bands' and 'single_events'
        """
        if events.empty:
            return {'group_bands': [], 'single_events': []}
        
        min_radius, max_radius = radius_range
        available_space = max_radius - min_radius
        
        # Separate multi-event groups from single events
        # SPECIAL CASE: BL groups (blacklist-crossing) should ALWAYS get group bands
        # even if they have only 1 event, for visibility
        group_counts = events['group'].value_counts().to_dict()
        
        # Multi-event groups: 2+ events OR any BL group (even with 1 event)
        multi_groups = sorted(
            [g for g, c in group_counts.items() if c > 1 or g.startswith('BL')],
            key=lambda g: group_counts[g]
        )  # Smallest first
        
        # Single events: 1 event AND not a BL group
        single_groups = [g for g, c in group_counts.items() if c == 1 and not g.startswith('BL')]
        
        logger.info(f"Ring for {event_type}: {len(multi_groups)} multi-event groups "
                   f"(including BL groups), {len(single_groups)} single events")
        
        group_bands = []
        current_radius = max_radius  # Start from outer edge
        
        # Layout multi-event groups
        for group_id in multi_groups:
            group_mask = events['group'] == group_id
            group_events = events[group_mask]
            n_events = len(group_events)
            
            # Calculate band size
            band_size = self._calculate_band_size(n_events)
            
            # Check if we have space
            if current_radius - band_size < min_radius:
                # Compress band to fit
                band_size = max(current_radius - min_radius, n_events * 2.0)
                logger.warning(f"Group {group_id} compressed to fit: {band_size:.1f} px")
            
            band_top = current_radius
            band_bottom = max(band_top - band_size, min_radius)
            
            # Layout events within band
            band_layouts = self._layout_group_band(
                events, 
                group_events.index.tolist(),
                group_id,
                band_top,
                band_bottom
            )
            group_bands.extend(band_layouts)
            
            # Move to next band with gap
            current_radius = band_bottom - self.layout_config.group_gap
            
            if current_radius <= min_radius:
                logger.warning(f"Ran out of space in {event_type} ring after {group_id}")
                current_radius = min_radius
        
        # Layout single events on shared levels
        single_event_indices = [events[events['group'] == g].index[0] for g in single_groups]
        single_layouts = self._layout_single_events(
            events,
            single_event_indices,
            current_radius,
            min_radius
        )
        
        return {
            'group_bands': group_bands,
            'single_events': single_layouts
        }
    
    def _calculate_band_size(self, n_events: int) -> float:
        """
        Calculate band size for a group
        
        Simple formula: base size + scale factor * number of events
        More conservative scaling to fit large groups like G2 (196 events)
        """
        base = self.layout_config.base_band_size
        scale = self.layout_config.band_size_scale_factor
        max_size = self.layout_config.max_band_size
        
        # Apply logarithmic scaling for very large groups to conserve space
        # This prevents a single large group (like G2 with 196 events) from taking all space
        if n_events > 50:
            # Use log scaling for large groups: base + scale * (50 + log(n-50))
            size = base + (scale * (50 + (10 * log(n_events - 49, 2))))
        else:
            size = base + (n_events * scale)
        
        return min(size, max_size)
    
    def _layout_group_band(
        self,
        events: pd.DataFrame,
        group_indices: List[int],
        group_id: str,
        band_top: float,
        band_bottom: float
    ) -> List[GroupBandLayout]:
        """
        Layout events within a group band with consistent spacing
        
        Args:
            events: Full DataFrame (modified in place)
            group_indices: Indices of events in this group
            group_id: Group identifier
            band_top: Outer radius of band
            band_bottom: Inner radius of band
        
        Returns:
            List with single GroupBandLayout (no splitting in simplified version)
        """
        n_events = len(group_indices)
        band_height = band_top - band_bottom
        
        if band_height <= 0:
            logger.error(f"Invalid band for group {group_id}: top={band_top}, bottom={band_bottom}")
            band_height = n_events * 2.0
            band_bottom = band_top - band_height
        
        # Calculate spacing - UNIFORM across all groups for consistent visual appearance
        min_spacing = self.layout_config.min_event_spacing
        
        if self.layout_config.uniform_within_group_spacing:
            # Use TARGET spacing for consistent dense appearance across all groups
            # This makes small groups (G3, G4) look as dense as large groups (G2, G5)
            spacing = self.layout_config.target_event_spacing
            
            # Check if we have enough space with target spacing
            required_height = n_events * spacing
            if required_height > band_height:
                # Not enough space - scale down proportionally
                spacing = band_height / n_events
                logger.warning(f"Group {group_id}: target spacing {self.layout_config.target_event_spacing:.1f} "
                             f"too large, using {spacing:.1f} px")
            
            # Still respect minimum
            if spacing < min_spacing:
                logger.warning(f"Group {group_id}: spacing {spacing:.1f} < min {min_spacing}")
                spacing = min_spacing
        else:
            # Original proportional spacing (varies by group size)
            if n_events > 1:
                spacing = band_height / (2 * n_events - 1)
                if spacing < min_spacing:
                    logger.warning(f"Group {group_id}: calculated spacing {spacing:.1f} < min {min_spacing}")
                    spacing = min_spacing
            else:
                spacing = band_height

        group_df = events.loc[group_indices].copy()
        group_df = group_df.sort_values('deg1')
        sorted_indices = group_df.index.tolist()
        
        # Assign radii with consistent spacing
        for i, idx in enumerate(sorted_indices):
            # Center events in their slots
            radius = band_top - (i * spacing) - (spacing / 2)
            # Ensure within band bounds
            radius = max(band_bottom, min(band_top, radius))
            events.loc[idx, 'radius'] = radius

        # Create layout descriptor
        return [GroupBandLayout(
            group_id=group_id,
            sub_band_index=0,  # No splitting in simplified version
            band_top=band_top,
            band_bottom=band_bottom,
            n_events=n_events,
            event_indices=sorted_indices,
            within_spacing=spacing,
            sector_id=0  # No sectors
        )]
    
    def _layout_single_events(
        self,
        events: pd.DataFrame,
        single_indices: List[int],
        start_radius: float,
        min_radius: float
    ) -> List[SingleEventLayout]:
        """
        Layout single events on shared radius levels
        
        Args:
            events: Full DataFrame (modified in place)
            single_indices: Indices of single events
            start_radius: Starting radius for first level
            min_radius: Minimum allowed radius
        
        Returns:
            List of SingleEventLayout objects
        """
        if not single_indices:
            return []
        
        layouts = []
        shared_levels = []
        current_radius = start_radius
        gap = self.layout_config.group_gap
        
        for idx in single_indices:
            event_deg = events.loc[idx, 'deg1']
            
            # Try to place on existing level
            placed = False
            for level in shared_levels:
                # Check angular spacing with other events on this level
                can_share = all(
                    self._angular_distance(event_deg, deg) >= 45.0  # 45 degrees min spacing
                    for deg in level['degrees']
                )
                
                if can_share and len(level['events']) < 8:  # Max 8 per level
                    level['events'].append(idx)
                    level['degrees'].append(event_deg)
                    events.loc[idx, 'radius'] = level['radius']
                    placed = True
                    break
            
            if not placed:
                # Create new level
                new_radius = max(current_radius - gap, min_radius)
                shared_levels.append({
                    'radius': new_radius,
                    'events': [idx],
                    'degrees': [event_deg]
                })
                events.loc[idx, 'radius'] = new_radius
                current_radius = new_radius
        
        # Create SingleEventLayout objects
        for level_id, level in enumerate(shared_levels):
            for idx in level['events']:
                layouts.append(SingleEventLayout(
                    event_index=idx,
                    radius=level['radius'],
                    sector_id=0,  # No sectors
                    shared_level_id=level_id,
                    n_events_on_level=len(level['events'])
                ))
        
        return layouts
    
    def _angular_distance(self, deg1: float, deg2: float) -> float:
        """Calculate minimum angular distance between two angles"""
        diff = abs(deg1 - deg2)
        return min(diff, 360 - diff)
    
    def _create_empty_layout(self, total_radius: float) -> LayoutResult:
        """Create empty layout for when there are no events"""
        return LayoutResult(
            events=pd.DataFrame(),
            sector_budgets={},
            group_bands=[],
            single_events=[],
            total_radius_used=0,
            del_radius_range=(0, 0),
            dup_radius_range=(0, 0),
            blacklist_radius=total_radius / 2,
            layout_algorithm='simplified_radial',
            layout_stats={}
        )