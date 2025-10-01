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
from saltshaker.config import ClassificationConfig
from saltshaker.event_caller import EventCaller
from saltshaker.spatial import SpatialGroupAnalyzer
from saltshaker.classifier import EventClassifier

import warnings
warnings.filterwarnings('ignore')

class MitoPlotter:
    def __init__(self, genome_length, ori_h_start, ori_h_end, ori_l_start, ori_l_end, 
                 heteroplasmy_limit=0.01, flank_size=15, config=None):
        """Initialize MitoPlotter with mitochondrial genome parameters"""
        self.genome_length = genome_length
        self.ori_h = (ori_h_start, ori_h_end)
        self.ori_l = (ori_l_start, ori_l_end)
        self.heteroplasmy_limit = heteroplasmy_limit
        self.flank_size = flank_size
        self.config = config or ClassificationConfig()
        
        # Create color maps for deletions (blues) and duplications (reds)
        self.del_cmap = LinearSegmentedColormap.from_list(
            'deletions', ['#4169E1', '#1E90FF', '#0000CD', '#000080', '#191970'])
        self.dup_cmap = LinearSegmentedColormap.from_list(
            'duplications', ['#FF6347', '#DC143C', '#B22222', '#8B0000', '#800000'])
    
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
        
        # Draw genome circle (use dynamic radius instead of hardcoded)
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
        print(f"Plotted {len(dat_processed)} events")

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
            'seq1': res['seq1'],  # Changed from flanking_results
            'seq2': res['seq2'],  # Changed from flanking_results
            'seq': res['seq']     # Changed from flanking_results
        })
        
        # Save results
        res_final.to_csv(output_file, sep='\t', index=False)
        print(f"Results saved to {output_file}")
        print(f"Events: {len(res_final)}")

    def save_summary(self, events, output_file, analysis_stats, classification_result=None, blacklist_regions=None):
        """Save analysis summary with biologically meaningful metrics based on MitoSAlt literature"""
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Load config for threshold reporting
        cfg = self.config
        
        # Use passed classification results
        if classification_result:
            classification, reason, criteria, events_with_groups = classification_result
        else:
            # Fallback if not provided
            if len(events) > 0:
                classifier = EventClassifier(self.genome_length, self.config)
                classification, reason, criteria, events_with_groups = classifier.classify(events, blacklist_regions)
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
            
            # Spatial groups analysis
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
                            f"max heteroplasmy {rep['heteroplasmy']:.1f}% "
                            f"(at {rep['start']:.0f}-{rep['end']:.0f}bp, size {actual_size:.0f}bp)\n")
                
                if dup_groups:
                    f.write("Duplication groups:\n")
                    for group_info in dup_groups:
                        rep = group_info['representative']
                        actual_size = rep['end'] - rep['start']
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1f}% "
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
                f.write(f"High heteroplasmy events (≥{cfg.HIGH_HETEROPLASMY_THRESHOLD:.0f}%): {criteria.get('high_het_count', 0)}\n")
                f.write(f"Medium heteroplasmy events ({cfg.LOW_HETEROPLASMY_THRESHOLD:.0f}-{cfg.HIGH_HETEROPLASMY_THRESHOLD:.0f}%): {criteria.get('medium_het_count', 0)}\n")
                f.write(f"Low heteroplasmy events (<{cfg.LOW_HETEROPLASMY_THRESHOLD:.0f}%): {criteria.get('low_het_count', 0)}\n")
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
                f.write(f"Single pattern score: {scores['single_score']} / 12 possible criteria\n")
                f.write(f"Multiple pattern score: {scores['multiple_score']} / 13 possible criteria\n")
                f.write(f"Decision: {'Single' if scores['single_score'] > scores['multiple_score'] else 'Multiple'} pattern (higher score wins)\n")
                f.write("\n")
                f.write("Score calculation:\n")
                f.write("- Each pattern type has biological criteria (12 for Single, 13 for Multiple)\n")
                f.write("- Single pattern criteria: dominant high-het events, few total events, tight clustering, etc.\n")
                f.write("- Multiple pattern criteria: many events, mixed types, scattered distribution, etc.\n")
                f.write("- Score = count of criteria met for each pattern type\n")
                f.write("\n")
                f.write("Biological thresholds used:\n")
                f.write(f"- High heteroplasmy threshold: ≥{cfg.HIGH_HETEROPLASMY_THRESHOLD:.0f}% (pathogenic significance)\n")
                f.write(f"- Significance threshold: ≥{cfg.SIGNIFICANT_HETEROPLASMY_THRESHOLD:.0f}% (above noise level)\n")
                f.write(f"- Low heteroplasmy threshold: <{cfg.LOW_HETEROPLASMY_THRESHOLD:.0f}% (likely artifacts)\n")
                f.write(f"- Multiple event threshold: >{cfg.TOTAL_EVENT_COUNT_THRESHOLD} events (mouse model pattern)\n")
                f.write(f"- Major event threshold: ≤{cfg.MAJOR_EVENT_COUNT_THRESHOLD} high-het events (single pattern)\n")
                f.write(f"- Spatial clustering radius: {cfg.CLUSTER_RADIUS}bp (biologically relevant)\n")
                f.write("\n")
            
            # Biological interpretation
            f.write("Biological interpretation:\n")
            f.write("-" * 25 + "\n")
            if classification == "Single":
                f.write("Pattern consistent with:\n")
                f.write("- Single pathogenic deletion/duplication\n")
                f.write("- Classical mitochondrial disease patient profile\n")
                f.write("- Possible accompanying low-level artifacts\n")
                if criteria.get('max_heteroplasmy', 0) >= cfg.HIGH_HETEROPLASMY_THRESHOLD:
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
            f.write(f"- Clinical thresholds: >{cfg.HIGH_HETEROPLASMY_THRESHOLD:.0f}% for pathogenic significance\n")
        
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
    parser.add_argument('--output-vcf', action='store_true', 
                    help='Also output events in VCF format')
    
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
    
    # Load and process data using EventCaller
    event_caller = EventCaller(
        genome_length=args.genome_length,
        ori_h=(args.ori_h_start, args.ori_h_end),
        ori_l=(args.ori_l_start, args.ori_l_end),
        heteroplasmy_limit=args.heteroplasmy_limit,
        flank_size=args.flank_size
    )

    events = event_caller.call_events(args.cluster_file, args.breakpoint_file)
    # ADD THIS - calculate final coordinates first (needed for flanking sequences)
    events['perc'] = events['perc'].round(4)
    events['del.start.median'] = events['del.start.median'] + 1
    events['final.event.size'] = np.where(
        events['final.event'] == 'del',
        events['delsize'],
        args.genome_length - events['delsize']
    )
    events['final.end'] = np.where(
        events['final.event'] == 'del',
        events['del.end.median'],
        events['del.start.median'] - 1
    )
    events['final.start'] = np.where(
        events['final.event'] == 'del',
        events['del.start.median'],
        events['del.end.median'] + 1
    )

    # ADD THIS - add flanking sequences
    events = event_caller.add_flanking_sequences(events, args.genome_fasta)

    stats = {}  # Will populate later
    
    if len(events) > 0:
        # Get events with group assignments from classification
    #    classification, reason, criteria, events_with_groups = plotter._classify_event_pattern(events, blacklist_regions=blacklist_regions)
        classifier = EventClassifier(args.genome_length, plotter.config)
        classification, reason, criteria, events_with_groups = classifier.classify(events, blacklist_regions=blacklist_regions)
        
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
        plotter.save_summary(
            events_with_groups, 
            str(summary_file), 
            stats, 
            classification_result=(classification, reason, criteria, events_with_groups), 
            blacklist_regions=blacklist_regions
            )

        # Add VCF output if requested
        if args.output_vcf:
            from saltshaker.io import write_vcf
            
            vcf_file = indel_dir / f"{args.output_name}.vcf"
            write_vcf(
                events_with_groups,
                (classification, reason, criteria, events_with_groups),
                str(vcf_file),
                reference_name="chrM",
                sample_name=args.output_name,
                genome_length=args.genome_length
            )
    else:
        print("No events to process")


if __name__ == "__main__":
    main()
