"""
Circular Visualizer

Creates circular genome plots for mitochondrial structural alterations.
"""

from asyncio import events
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import re
import logging
from .utils import crosses_blacklist

logger = logging.getLogger(__name__)


class CircularPlotter:
    """
    Creates circular genome visualizations of mitochondrial events
    
    Generates publication-quality circular plots showing deletions and 
    duplications around the mitochondrial genome with spatial grouping.
    """
    
    def __init__(self, genome_length):
        """
        Initialize CircularPlotter
        
        Args:
            genome_length: Mitochondrial genome length
        """
        self.genome_length = genome_length

        # Color maps will be set in plot() based on parameters
        self.del_cmap = None
        self.dup_cmap = None

    def plot(self, events, output_file, blacklist_regions=None, figsize=(16, 10),
             direction='counterclockwise', del_color='red', dup_color='blue',
             gene_annotations=None, scale='dynamic'):
        """
        Create circular plot of mitochondrial events
        
        Args:
            events: DataFrame with events (must have group column)
            output_file: Path to output PNG file
            blacklist_regions: List of blacklist region dicts (optional)
            figsize: Figure size tuple (default: (16, 10))
            direction: 'clockwise' or 'counterclockwise' (default: 'counterclockwise')
            del_color: Color scheme for deletions - 'red' or 'blue' (default: 'blue')
            dup_color: Color scheme for duplications - 'red' or 'blue' (default: 'red')
            gene_annotations: List of gene annotation dicts from BED file (optional)
            scale: Heteroplasmy color scale - 'dynamic' (min-max per category) or 'fixed' (0-100%) (default: 'dynamic')

        """
        
        if len(events) == 0:
            logger.warning("No events to plot")
            return
        
        # Validate scale parameter
        if scale not in ['dynamic', 'fixed']:
            raise ValueError(f"Invalid scale: {scale}. Use 'dynamic' or 'fixed'")
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)

        # Set color maps based on parameters
        if del_color == 'red':
            self.del_cmap = LinearSegmentedColormap.from_list(
                'deletions', ['#FF6347', '#DC143C', '#B22222', '#8B0000', '#800000'])
            del_label_color = 'red'
        elif del_color == 'blue':
            self.del_cmap = LinearSegmentedColormap.from_list(
                'deletions', ['#4169E1', '#1E90FF', '#0000CD', '#000080', '#191970'])
            del_label_color = 'blue'
        else:
            raise ValueError(f"Invalid del_color: {del_color}. Use 'red' or 'blue'")

        if dup_color == 'red':
            self.dup_cmap = LinearSegmentedColormap.from_list(
                'duplications', ['#FF6347', '#DC143C', '#B22222', '#8B0000', '#800000'])
            dup_label_color = 'red'
        elif dup_color == 'blue':
            self.dup_cmap = LinearSegmentedColormap.from_list(
                'duplications', ['#4169E1', '#1E90FF', '#0000CD', '#000080', '#191970'])
            dup_label_color = 'blue'
        else:
            raise ValueError(f"Invalid dup_color: {dup_color}. Use 'red' or 'blue'")
        
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
            # Add temporary sort key column
            dat['_sort_key'] = dat['group'].apply(lambda g: self.group_sort_key(g))
            dat = dat.sort_values(['_sort_key', 'value'], ascending=[True, False]).reset_index(drop=True)
            dat = dat.drop(columns=['_sort_key'])
        
        # Enhanced blacklist crossing detection      
        dat['blacklist_crossing'] = False
        if blacklist_regions:
            for idx, row in dat.iterrows():
                dat.loc[idx, 'blacklist_crossing'] = crosses_blacklist(row['start'], row['end'], blacklist_regions)

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
                min_gap = 5
                
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
            
            group_band_size = 25
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
                        data.loc[idx, 'radius'] = max(20, assignment['radius'])  # Add minimum bound to prevent overflow to the boundary
                else:
                    band_top, band_bottom = assignment['band_top'], assignment['band_bottom']
                    band_bottom = max(20, band_bottom)  # Ensure minimum radius
                    for i, idx in enumerate(group_indices):
                        test_radius = band_top
                        while test_radius >= band_bottom:
                            if not any(data.loc[prev_idx, 'radius'] == test_radius and 
                                    events_overlap(data.loc[idx], data.loc[prev_idx])
                                    for prev_idx in group_indices[:i]):
                                data.loc[idx, 'radius'] = test_radius
                                break
                            test_radius -= 5
                        else:
                            data.loc[idx, 'radius'] = max(20, band_bottom)  # Add minimum bound
            
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
            
            # Determine outer vs inner based on MEDIAN SIZE (largest type goes outside)
            # Keep type separation: all dels in one ring, all dups in another
            del_median_size = dat_del['delsize'].median() if len(dat_del) > 0 else 0
            dup_median_size = dat_dup['delsize'].median() if len(dat_dup) > 0 else 0

            if dup_median_size > del_median_size:
                # Duplications are larger → put ALL dups outside, ALL dels inside
                outer_frac, inner_frac = dup_frac, del_frac
                outer_data, inner_data = dat_dup, dat_del
                outer_type = 'dup'
            else:
                # Deletions are larger (or equal) → put ALL dels outside, ALL dups inside
                outer_frac, inner_frac = del_frac, dup_frac  
                outer_data, inner_data = dat_del, dat_dup
                outer_type = 'del'

            logger.debug(f"Layout - {outer_type.upper()} outside (median size: {max(del_median_size, dup_median_size):.0f}bp)")
            
            # Calculate radius ranges
            inner_max = base_radius * inner_frac
            outer_min = base_radius * (inner_frac + separator_frac)
            blacklist_radius = (inner_max + outer_min) / 2
            
            logger.debug(f"Fractions - Inner: {inner_frac:.3f}, Separator: {separator_frac:.3f}, Outer: {outer_frac:.3f}")
            logger.debug(f"Ranges - Inner: [0-{inner_max:.1f}], Outer: [{outer_min:.1f}-{base_radius}], BL: {blacklist_radius:.1f}")
            
            # Assign radii within ranges
            if len(inner_data) > 0:
                inner_data = assign_radii_by_type(inner_data, base_radius=inner_max, radius_diff=6)
            if len(outer_data) > 0:
                outer_data = assign_radii_by_type(outer_data, base_radius=base_radius, radius_diff=6)

            circle_radius = base_radius + 12  # Circle drawn slightly outside events
            
            # Return in correct order
            if outer_type == 'dup':
                return inner_data, outer_data, blacklist_radius, circle_radius
            else:
                return outer_data, inner_data, blacklist_radius, circle_radius

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
        
        # Calculate color scales based on scale parameter
        del_events = dat_processed[dat_processed['final.event'] == 'del']
        dup_events = dat_processed[dat_processed['final.event'] == 'dup']
        
        if scale == 'fixed':
            # Fixed scale: 0-100% for all categories
            del_min = 0.0
            del_max = 100.0
            dup_min = 0.0
            dup_max = 100.0
            bl_min = 0.0
            bl_max = 100.0
            logger.info("Using fixed heteroplasmy scale: 0-100%")
        else:  # dynamic
            # Dynamic scale: min-max within each category
            del_max = del_events['value'].max() if len(del_events) > 0 else 0
            del_min = del_events['value'].min() if len(del_events) > 0 else 0
            dup_max = dup_events['value'].max() if len(dup_events) > 0 else 0  
            dup_min = dup_events['value'].min() if len(dup_events) > 0 else 0
            # Calculate BL min/max for gradient coloring (needed for event plotting)
            if blacklist_regions:
                bl_events_temp = dat_processed[dat_processed['blacklist_crossing'] == True]
                if len(bl_events_temp) > 0:
                    bl_min = bl_events_temp['value'].min()
                    bl_max = bl_events_temp['value'].max()
                else:
                    bl_min = 0.0
                    bl_max = 100.0
            else:
                bl_min = 0.0
                bl_max = 100.0
            logger.info(f"Using dynamic heteroplasmy scale - Del: {del_min:.1f}-{del_max:.1f}%, Dup: {dup_min:.1f}-{dup_max:.1f}%")
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='polar', position=[0.15, 0.05, 0.7, 0.88])
        # Set y-limit to accommodate gene track if present
        max_radius = dynamic_radius + 90 if gene_annotations else dynamic_radius + 30
        ax.set_ylim(0, max_radius)
        ax.set_theta_zero_location('N')

        # Set direction based on parameter
        if direction == 'counterclockwise':
            ax.set_theta_direction(1)  # Counterclockwise
        elif direction == 'clockwise':
            ax.set_theta_direction(-1)  # Clockwise
        else:
            raise ValueError(f"Invalid direction: {direction}. Use 'clockwise' or 'counterclockwise'")
        
        # Draw genome circle (use dynamic radius instead of hardcoded)
        circle = patches.Circle((0, 0), dynamic_radius, fill=False, linewidth=3,
                            color='gray', transform=ax.transData._b)
        ax.add_patch(circle)

        # Draw gene annotation track (outer ring) - NEW CODE STARTS HERE
        if gene_annotations:
            gene_track_inner = dynamic_radius + 25  # Start 25 units outside main circle (was 15)
            gene_track_outer = dynamic_radius + 40  # 15 units tall (was 20)
            gene_track_mid = (gene_track_inner + gene_track_outer) / 2
            
            # Draw track background circle
            track_bg = patches.Circle((0, 0), gene_track_outer, fill=False, linewidth=1,
                                     color='lightgray', linestyle='--', alpha=0.5,
                                     transform=ax.transData._b)
            ax.add_patch(track_bg)
            
            # Draw each gene
            for gene in gene_annotations:
                start_deg = np.radians(358 * gene['start'] / self.genome_length)
                end_deg = np.radians(358 * gene['end'] / self.genome_length)
                
                # Gene arc
                theta = np.linspace(start_deg, end_deg, 50)
                ax.fill_between(theta, gene_track_inner, gene_track_outer,
                                color=gene['color'], alpha=0.8, linewidth=0)
                
                # Gene name (for larger genes only to avoid clutter)
                gene_size = gene['end'] - gene['start']
                if gene_size > 200:  # Only label genes >200bp
                    mid_deg = (start_deg + end_deg) / 2
                    label_radius = gene_track_outer + 12  # Place just outside the track
                    
                    # Calculate rotation for text to be horizontal and readable
                    # Convert back to degrees for rotation calculation
                    rotation_deg = np.degrees(mid_deg)
                    
                    # When in clockwise mode, we need to mirror the angle
                    # because set_theta_direction(-1) mirrors the coordinate system
                    if direction == 'clockwise':
                        # Mirror the angle: 0° stays 0°, 90° becomes 270°, etc.
                        rotation_deg = 360 - rotation_deg
                    
                    # Now apply the standard logic: flip text on bottom half
                    if rotation_deg > 90 and rotation_deg < 270:
                        rotation_deg = rotation_deg + 180
                    
                    ax.text(mid_deg, label_radius, gene['name'],
                           rotation=rotation_deg, ha='center', va='center',
                           fontsize=7, weight='bold', color='black')

        
        # Draw separator circle between del and dup radii (always present)
        separator_circle = patches.Circle((0, 0), blacklist_radius + 15, fill=False, linewidth=2, 
                                        color='lightgray', linestyle='--', alpha=0.7,
                                        transform=ax.transData._b)
        ax.add_patch(separator_circle)
        
        # Add blacklist regions
        if blacklist_regions:
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
        
        # Add position markers - adjusted for gene track spacing
        positions = np.arange(0, self.genome_length, 1000)
        for pos in positions:
            deg = np.radians(358 * pos / self.genome_length)
            marker_start = gene_track_outer + 5 if gene_annotations else dynamic_radius
            text_offset = 18 if gene_annotations else 12  # More space when genes present
            ax.plot([deg, deg], [marker_start, marker_start + 5], color='gray', linewidth=1)
            ax.text(deg, marker_start + text_offset, f'{pos//1000}', ha='center', va='center', 
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

        def get_lime_green_color(het_val, min_het, max_het):
            """Lime-green gradient for BL-crossing events"""
            if max_het > min_het:
                norm = (het_val - min_het) / (max_het - min_het)
            else:
                norm = 0.5
            # Light lime-green (0.6, 1.0, 0.4) to dark green (0.0, 0.5, 0.0)
            red = 0.6 * (1 - norm)
            green = 1.0 - 0.5 * norm
            blue = 0.4 * (1 - norm)
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
                # Use gradient coloring for BL events based on heteroplasmy
                color = get_lime_green_color(het_val, bl_min, bl_max)
                alpha = get_continuous_alpha(het_val, bl_min, bl_max)
                logger.debug(f"Plotting blacklist event at {event['start']}-{event['end']}, het={het_val:.1f}%")
            else:
                if event['final.event'] == 'del':
                    # Use the color function based on del_color parameter
                    del_color_func = get_pure_red_color if del_color == 'red' else get_pure_blue_color
                    color = del_color_func(het_val, del_min, del_max)
                    alpha = get_continuous_alpha(het_val, del_min, del_max)
                else:
                    # Use the color function based on dup_color parameter
                    dup_color_func = get_pure_red_color if dup_color == 'red' else get_pure_blue_color
                    color = dup_color_func(het_val, dup_min, dup_max) 
                    alpha = get_continuous_alpha(het_val, dup_min, dup_max)
                
            theta = np.linspace(deg1_rad, deg2_rad, 100)
            linewidth = 2.5 if len(dat_processed) <= 100 else (2.0 if len(dat_processed) <= 200 else 1.5)
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
        # 1. EVENT COUNT SUMMARY - Unified box with all counts
        # Draw the box first
        from matplotlib.patches import FancyBboxPatch
        
        if blacklist_regions and (bl_del_count > 0 or bl_dup_count > 0):
            # Box for 2 lines
            box_height = 0.065
            box_y = 0.80
            text_y = 0.85
            box = FancyBboxPatch((0.055, box_y), 0.17, box_height,
                                boxstyle="round,pad=0.008", 
                                facecolor='white', edgecolor='gray',
                                alpha=0.9, transform=fig.transFigure)
            fig.patches.append(box)
            
            # First line: Del and Dup counts
            fig.text(0.06, text_y, f"Del: {del_count}  Dup: {dup_count}", 
                    fontsize=13, weight='bold', verticalalignment='top')
            
            # Second line - all in lime-green
            y_bl = 0.820
            x_start = 0.06

            # Entire BL line in green for clarity
            fig.text(x_start, y_bl, f"BL Del: {bl_del_count}  BL Dup: {bl_dup_count}", 
                    fontsize=13, weight='bold', color=(0.2, 0.8, 0.2), verticalalignment='top')
        else:
            # Single line box - adjusted to align properly
            box_height = 0.035
            box_y = 0.83
            text_y = 0.855
            box = FancyBboxPatch((0.055, box_y), 0.14, box_height,
                                boxstyle="round,pad=0.008", 
                                facecolor='white', edgecolor='gray',
                                alpha=0.9, transform=fig.transFigure)
            fig.patches.append(box)
            
            fig.text(0.06, text_y, f"Del: {del_count}  Dup: {dup_count}", 
                    fontsize=13, weight='bold', verticalalignment='top')
        
        
        def draw_gradient_legend(fig, legend_x, legend_y, legend_height, bar_width, 
                                events, min_val, max_val, label, color_func, label_color, offset_x=0, is_bl=False):
            """Draw a single gradient legend bar with labels"""
            if len(events) == 0:
                return
            
            n_steps = 100
            step_height = legend_height / n_steps
            
            # Draw gradient rectangles
            for i in range(n_steps):
                norm = i / (n_steps - 1)
                het_val = min_val + norm * (max_val - min_val) if max_val > min_val else min_val
                y_pos = legend_y - legend_height/2 + i * step_height
                
                color = color_func(het_val, min_val, max_val)
                alpha = get_continuous_alpha(het_val, min_val, max_val)
                
                rect = plt.Rectangle(
                    (legend_x + offset_x, y_pos), bar_width, step_height,
                    facecolor=color, alpha=alpha, edgecolor='none',
                    transform=fig.transFigure
                )
                fig.patches.append(rect)
            
            # Label at top
            label_x = legend_x + offset_x + bar_width / 2
            fig.text(label_x, legend_y + legend_height/2 + 0.01, label, 
                    fontsize=10, ha='center', weight='bold', color=label_color)
            
            # Min/max values - all on the right side, closer to bars
            # Use 2 decimal places for BL values, 1 for Del/Dup
            decimal_places = 2 if is_bl else 1
            for pos, val in [(1, max_val), (0, min_val)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height
                fig.text(legend_x + offset_x + bar_width + 0.005, y_pos, f"{val:.{decimal_places}f}", 
                        fontsize=9, va='center', ha='left', color=label_color)


        # 2. SEPARATE GRADIENT LEGENDS - independently scaled
        legend_x = 0.08
        legend_y = 0.48
        legend_height = 0.30
        
        # Calculate BL events for legend display (bl_min and bl_max already calculated above)
        bl_events = pd.DataFrame()
        show_bl_bar = False
        
        if blacklist_regions and (bl_del_count > 0 or bl_dup_count > 0):
            bl_events = dat_processed[dat_processed['blacklist_crossing'] == True]
            if len(bl_events) > 0:
                show_bl_bar = True
        
        # Determine number of bars and calculate widths with more spacing
        num_bars = 3 if show_bl_bar else 2
        bar_gap = 0.030  # Increased for even more space
        legend_width = 0.04 if num_bars == 2 else 0.070  # Adjusted for new spacing
        bar_width = (legend_width - (num_bars - 1) * bar_gap) / num_bars

        # Choose color functions based on parameters
        del_color_func = get_pure_red_color if del_color == 'red' else get_pure_blue_color
        dup_color_func = get_pure_red_color if dup_color == 'red' else get_pure_blue_color

        # Draw deletion legend (first bar)
        draw_gradient_legend(fig, legend_x, legend_y, legend_height, bar_width,
                            del_events, del_min, del_max, "Del", del_color_func, del_label_color, offset_x=0)

        # Draw duplication legend (second bar)
        draw_gradient_legend(fig, legend_x, legend_y, legend_height, bar_width,
                            dup_events, dup_min, dup_max, "Dup", dup_color_func, dup_label_color, 
                            offset_x=bar_width + bar_gap)

        # Draw BL legend (third bar) if applicable
        if show_bl_bar:
            bl_offset = 2 * (bar_width + bar_gap)
            draw_gradient_legend(fig, legend_x, legend_y, legend_height, bar_width,
                                bl_events, bl_min, bl_max, "BL", get_lime_green_color, (0.2, 0.8, 0.2), 
                                offset_x=bl_offset, is_bl=True)

        # Heteroplasmy label above all bars
        fig.text(legend_x + legend_width/2, legend_y + legend_height/2 + 0.05, 
                "Heteroplasmy (%)", fontsize=11, weight='bold', ha='center')
        
        # Gene annotation legend (if present) - moved to left side
        if gene_annotations:
            legend_items = [
                ('Protein-coding', (0, 0.5, 0)),      # Green
                ('tRNA', (1, 1, 0)),                  # Yellow
                ('rRNA', (0.5, 0, 0.5))               # Purple
            ]
            
            # Position on left side, below the heteroplasmy legend
            legend_y_start = legend_y - legend_height/2 - 0.1  # Below heteroplasmy label
            legend_box_size = 0.025
            
            # Header for gene legend - centered and aligned with Heteroplasmy title
            fig.text(legend_x + legend_width/2, legend_y_start + 0.04, 
                    "Genes", fontsize=11, weight='bold', ha='center')
            
            for i, (label, color) in enumerate(legend_items):
                y_pos = legend_y_start - i * 0.04
                
                # Color box
                rect = plt.Rectangle((legend_x, y_pos), legend_box_size, 0.025,
                                    facecolor=color, alpha=0.8,
                                    transform=fig.transFigure)
                fig.patches.append(rect)
                
                # Label
                fig.text(legend_x + legend_box_size + 0.01, y_pos + 0.0125,
                        label, fontsize=9, va='center')
        
        ax.grid(False)
        ax.set_xticklabels([])

        
        ax.grid(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        title = 'Mitochondrial DNA Structural Alterations'
        if blacklist_regions:
            title += f' (BL: {len(blacklist_regions)} regions)'
        
        fig.suptitle(title, fontsize=15, weight='bold', y=0.98)
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Plot saved to {output_file}")
        logger.info(f"Plotted {len(dat_processed)} events")

        
    def group_sort_key(self, group_id):
        """Convert group ID to sortable tuple (priority, number)"""
        match = re.match(r'^([A-Z]+)(\d+)$', group_id)
        if match:
            prefix, number = match.groups()
            if prefix == 'G':
                return (0, int(number))
            elif prefix == 'BL':
                return (1000, int(number))
            else:
                logger.warning(f"Unexpected group ID format: {group_id}")
                return (9999, int(number))
        else:
            logger.warning(f"Could not parse group ID: {group_id}")
            return (9999, 0)


def plot_circular(events, output_file, genome_length, blacklist_regions=None, 
                  figsize=(16, 10), direction='counterclockwise', 
                  del_color='red', dup_color='blue', gene_annotations=None, scale='dynamic'):
    """
    Convenience function to create circular plot
    
    Args:
        events: DataFrame with events
        output_file: Output PNG file path
        genome_length: Mitochondrial genome length
        blacklist_regions: List of blacklist regions
        figsize: Figure size tuple
        direction: 'clockwise' or 'counterclockwise' (default: 'counterclockwise')
        del_color: Color for deletions - 'red' or 'blue' (default: 'red')
        dup_color: Color for duplications - 'red' or 'blue' (default: 'blue')
        gene_annotations: List of gene annotation dicts (optional)
        scale: Heteroplasmy color scale - 'dynamic' (min-max per category) or 'fixed' (0-100%) (default: 'dynamic')
    """
    plotter = CircularPlotter(genome_length)
    plotter.plot(events, output_file, blacklist_regions, figsize, direction, del_color, dup_color, gene_annotations, scale)