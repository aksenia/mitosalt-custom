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
from .utils import crosses_blacklist


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
        
        # Create color maps for deletions (blues) and duplications (reds)
        self.del_cmap = LinearSegmentedColormap.from_list(
            'deletions', ['#4169E1', '#1E90FF', '#0000CD', '#000080', '#191970'])
        self.dup_cmap = LinearSegmentedColormap.from_list(
            'duplications', ['#FF6347', '#DC143C', '#B22222', '#8B0000', '#800000'])
    
    def plot(self, events, output_file, blacklist_regions=None, figsize=(16, 10)):
        """
        Create circular plot of mitochondrial events
        
        Args:
            events: DataFrame with events (must have group column)
            output_file: Path to output PNG file
            blacklist_regions: List of blacklist region dicts (optional)
            figsize: Figure size tuple (default: (16, 10))
        """
        
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
                            test_radius -= 2
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
            
            # Determine outer vs inner based on group numbers
            # Extract minimum group number for ordering (handles G1, BL1, etc.)
            def extract_group_num(group_id):
                """Extract numeric part from group ID"""
                match = re.match(r'^[A-Z]+(\d+)$', str(group_id))
                return int(match.group(1)) if match else 999

            del_min_group = min([extract_group_num(g) for g in dat_del['group'].unique()]) if len(dat_del) > 0 else 999
            dup_min_group = min([extract_group_num(g) for g in dat_dup['group'].unique()]) if len(dat_dup) > 0 else 999
            
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
        
        # Calculate color scales
        del_events = dat_processed[dat_processed['final.event'] == 'del']
        dup_events = dat_processed[dat_processed['final.event'] == 'dup']
        del_max = del_events['value'].max() if len(del_events) > 0 else 0
        del_min = del_events['value'].min() if len(del_events) > 0 else 0
        dup_max = dup_events['value'].max() if len(dup_events) > 0 else 0  
        dup_min = dup_events['value'].min() if len(dup_events) > 0 else 0
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='polar', position=[0.15, 0.05, 0.7, 0.88])
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
        # 1. EVENT COUNT SUMMARY - Top left area with colored BL text
        fig.text(0.06, 0.85, f"Del: {del_count}  Dup: {dup_count}", 
                fontsize=13, weight='bold',
                bbox=dict(boxstyle="round,pad=0.5", facecolor='white', alpha=0.9, edgecolor='gray'),
                verticalalignment='top')

        if blacklist_regions and (bl_del_count > 0 or bl_dup_count > 0):
            # Add BL-crossing text in green below main counts
            fig.text(0.06, 0.8, f"BL-crossing Del: {bl_del_count}\nBL-crossing Dup: {bl_dup_count}", 
                    fontsize=13, weight='bold', color=(0.2, 0.8, 0.2),
                    verticalalignment='top')
        
        # 2. SEPARATE GRADIENT LEGENDS - independently scaled
        legend_x = 0.08
        legend_y = 0.5
        legend_height = 0.35
        legend_width = 0.045

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
                fontsize=9, ha='center', weight='bold', color='blue')

        dup_label_x = legend_x + bar_width + bar_gap + bar_width / 2
        fig.text(dup_label_x, legend_y + legend_height/2 + 0.01, "Dup", 
                fontsize=9, ha='center', weight='bold', color='red')

        # Separate value labels for each scale
        if len(del_events) > 0:
            for pos, val in [(1, del_max), (0, del_min)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height
                fig.text(legend_x - 0.015, y_pos, f"{val:.1f}%", 
                        fontsize=10, va='center', ha='right', color='blue')

        if len(dup_events) > 0:
            for pos, val in [(1, dup_max), (0, dup_min)]:
                y_pos = legend_y - legend_height/2 + pos * legend_height
                fig.text(legend_x + legend_width + 0.015, y_pos, f"{val:.1f}%",
                        fontsize=10, va='center', color='red')

        fig.text(legend_x + legend_width/2, legend_y - legend_height/2 - 0.03, 
                "Heteroplasmy (%)", fontsize=12, weight='bold', ha='center')

        
        ax.grid(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        title = 'Mitochondrial DNA Structural Alterations'
        if blacklist_regions:
            title += f' (BL: {len(blacklist_regions)} regions)'
        
        fig.suptitle(title, fontsize=15, weight='bold', y=0.98)
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved to {output_file}")
        print(f"Plotted {len(dat_processed)} events")

        
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
                print(f"WARNING: Unexpected group ID format: {group_id}")
                return (9999, int(number))
        else:
            print(f"WARNING: Could not parse group ID: {group_id}")
            return (9999, 0)


def plot_circular(events, output_file, genome_length, blacklist_regions=None, figsize=(16, 10)):
    """
    Convenience function to create circular plot
    
    Args:
        events: DataFrame with events
        output_file: Output PNG file path
        genome_length: Mitochondrial genome length
        blacklist_regions: List of blacklist regions
        figsize: Figure size tuple
    """
    plotter = CircularPlotter(genome_length)
    plotter.plot(events, output_file, blacklist_regions, figsize)