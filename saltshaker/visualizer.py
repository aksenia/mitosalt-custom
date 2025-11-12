"""
Circular visualizer - SIMPLIFIED VERSION

Creates circular genome plots for mitochondrial structural alterations.
All legacy code removed - uses only simplified radial layout.
"""

from __future__ import annotations
from typing import List, Tuple, Dict, Any, Optional, Callable, Union, Literal
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.projections.polar import PolarAxes
from pathlib import Path
import re
import logging

from .utils import crosses_blacklist
from .types import BlacklistRegion, GeneAnnotation
from .config import PlotConfig
from .layout import LayoutEngine

logger = logging.getLogger(__name__)

class CircularPlotter:
    """
    Creates circular genome visualizations of mitochondrial events
    
    Simplified version with clean radial layout and no legacy code.
    """
    
    def __init__(self, genome_length: int, config: Optional[PlotConfig] = None) -> None:
        """
        Initialize CircularPlotter
        
        Args:
            genome_length: Mitochondrial genome length (typically 16569 for human)
            config: Visual configuration for plot styling. If None, uses default settings.
        
        Example:
            >>> plotter = CircularPlotter(16569)
            >>> plotter = CircularPlotter(16569, PlotConfig.publication())
        """
        self.genome_length: int = genome_length
        self.config: PlotConfig = config or PlotConfig()

        # Color maps will be set in plot() based on parameters
        self.del_cmap: Optional[LinearSegmentedColormap] = None
        self.dup_cmap: Optional[LinearSegmentedColormap] = None
        self.layout_engine = LayoutEngine(self.config, genome_length)

    def _calculate_local_density(
        self,
        event_position: float,
        all_events: pd.DataFrame,
        window_degrees: Optional[float] = None
    ) -> int:
        """
        Calculate local density of events around a position
        
        Used for adaptive line width calculation.
        
        Args:
            event_position: Angular position in degrees (0-360)
            all_events: DataFrame with all events (must have 'deg1' column)
            window_degrees: Angular window size in degrees (default: from config)
            
        Returns:
            Number of events within the window
        """
        if window_degrees is None:
            window_degrees = self.config.visualization.density_window_degrees
        
        half_window = window_degrees / 2
        
        # Handle wraparound at 0/360 boundary
        if event_position < half_window:
            # Window crosses 360 -> 0
            lower_bound = 360 + (event_position - half_window)
            upper_bound = event_position + half_window
            in_window = ((all_events['deg1'] >= lower_bound) | 
                        (all_events['deg1'] <= upper_bound))
        elif event_position > 360 - half_window:
            # Window crosses 0 -> 360
            lower_bound = event_position - half_window
            upper_bound = (event_position + half_window) - 360
            in_window = ((all_events['deg1'] >= lower_bound) | 
                        (all_events['deg1'] <= upper_bound))
        else:
            # Normal case - no wraparound
            lower_bound = event_position - half_window
            upper_bound = event_position + half_window
            in_window = ((all_events['deg1'] >= lower_bound) & 
                        (all_events['deg1'] <= upper_bound))
        
        return in_window.sum()
    
    def _get_adaptive_linewidth(
        self,
        event_position: float,
        all_events: pd.DataFrame
    ) -> float:
        """
        Calculate adaptive line width based on local density
        
        Args:
            event_position: Angular position of the event
            all_events: DataFrame with all events
        
        Returns:
            Line width in pixels
        """
        # Use FIXED line width for consistent appearance - same as outer red arcs
        # The layout algorithm will handle spacing intelligently
        return 2.5  # Fixed thickness for all arcs
    
    def _compress_gradient_fixed_scale(self, value: float, vmin: float, vmax: float) -> float:
        """
        Compress gradient for fixed scale (0-100%) to optimize visibility
        
        Maps values as follows:
        - 0-50%: Uses 0-85% of gradient (good visibility)
        - 50-100%: Uses 85-100% of gradient (compressed dark shades)
        
        This gives better visibility for the common 0-50% range while
        still allowing theoretical 50-100% values.
        
        Args:
            value: Heteroplasmy value
            vmin: Minimum value (should be 0 for fixed scale)
            vmax: Maximum value (should be 100 for fixed scale)
            
        Returns:
            Compressed normalized value (0-1) for colormap lookup
        """
        if vmax <= vmin:
            return 0.5
        
        # First normalize to 0-1 based on full range
        norm = (value - vmin) / (vmax - vmin)
        
        # For fixed scale (0-100%), compress the gradient:
        # 0-50% (0-0.5 normalized) → 0-0.85 of gradient
        # 50-100% (0.5-1.0 normalized) → 0.85-1.0 of gradient
        if norm <= 0.5:
            # Map 0-0.5 to 0-0.85 (linear)
            compressed = norm * 1.7  # 0.5 * 1.7 = 0.85
        else:
            # Map 0.5-1.0 to 0.85-1.0 (compressed)
            compressed = 0.85 + (norm - 0.5) * 0.3  # Remaining 0.15 for top half
        
        return min(1.0, max(0.0, compressed))
    
    def group_sort_key(self, group: str) -> Tuple:
        """
        Generate sort key for group names
        
        Ensures consistent ordering: G1, G2, ..., G10, G11, ..., BL1, BL2, ...
        """
        match = re.match(r'^([A-Z]+)(\d+)$', group)
        if match:
            prefix, number = match.groups()
            return (prefix, int(number))
        return (group, 0)
    
    def calculate_pattern(
        self,
        events: pd.DataFrame,
        classification_config: Optional[Any] = None
    ) -> Literal['Single', 'Multiple', 'Unknown']:
        """
        Calculate Single vs Multiple pattern classification
        
        Args:
            events: DataFrame with events (must have 'group' column)
            classification_config: Optional classification config
        
        Returns:
            'Single', 'Multiple', or 'Unknown'
        """
        if 'group' not in events.columns:
            return 'Unknown'
        
        # Get configuration thresholds
        if classification_config is None:
            from .config import ClassificationConfig
            classification_config = ClassificationConfig()
        
        # Count events by group
        group_counts = events['group'].value_counts()
        n_groups = len(group_counts)
        n_events = len(events)
        
        # Multiple pattern: >10 events OR multiple groups with no dominant group
        if n_events > classification_config.MULTIPLE_EVENT_THRESHOLD:
            return 'Multiple'
        
        # Single pattern: one dominant group with ≥50% of events
        if n_groups == 1:
            return 'Single'
        
        largest_group_fraction = group_counts.iloc[0] / n_events
        if largest_group_fraction >= classification_config.DOMINANT_GROUP_FRACTION:
            return 'Single'
        
        return 'Multiple'
    
    def plot(
        self,
        events: pd.DataFrame,
        output_file: str = 'mitochondrial_plot.png',
        title: Optional[str] = None,
        metadata_text: Optional[str] = None,
        scale: Literal['dynamic', 'fixed'] = 'dynamic',
        blacklist_regions: Optional[List[BlacklistRegion]] = None,
        gene_annotations: Optional[List[GeneAnnotation]] = None,
        del_color: Literal['red', 'blue'] = 'red',
        dup_color: Literal['blue', 'red'] = 'blue',
        figsize: Optional[Tuple[float, float]] = None,
        direction: Literal['clockwise', 'counterclockwise'] = 'counterclockwise',
        show: bool = False
    ) -> Figure:
        """
        Generate circular mitochondrial plot
        
        Args:
            events: DataFrame with columns:
                - del_start_median, del_end_median: Event boundaries (bp)
                - perc: Heteroplasmy percentage (0-100)
                - final_event: Event type ('del' or 'dup')
                - delsize: Event size (bp)
                - dloop: D-loop involvement ('yes'/'no')
                - group: Optional group identifier
            output_file: Path to save figure
            title: Plot title (auto-generated if None)
            metadata_text: Additional metadata to display
            scale: 'dynamic' (scale to data range) or 'fixed' (0-100%)
            blacklist_regions: Optional list of blacklist regions
            gene_annotations: Optional list of gene annotations
            del_color: Color for deletions ('red' or 'blue')
            dup_color: Color for duplications ('red' or 'blue')
            figsize: Figure size in inches (width, height)
            direction: Plot direction ('clockwise' or 'counterclockwise')
            show: Whether to display the plot
        
        Returns:
            matplotlib Figure object
        
        Example:
            >>> events_df = pd.read_csv('events.csv')
            >>> fig = plotter.plot(events_df, 'output.png', scale='dynamic')
        """
        if del_color not in ['red', 'blue'] or dup_color not in ['red', 'blue']:
            raise ValueError("Colors must be 'red' or 'blue'")
        
        if scale not in ['dynamic', 'fixed']:
            raise ValueError(f"Invalid scale: {scale}. Use 'dynamic' or 'fixed'")
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)

        # Set color maps based on parameters - with stronger gradients
        if del_color == 'red':
            # From very light red to dark red (no pink)
            self.del_cmap = LinearSegmentedColormap.from_list(
                'deletions', ['#FFCCCC', '#FF9999', '#FF6666', '#CC0000', '#660000'])
            del_label_color = 'red'
        else:
            # From very pale blue to dark blue
            self.del_cmap = LinearSegmentedColormap.from_list(
                'deletions', ['#E6F3FF', '#B3D9FF', '#66B2FF', '#0066CC', '#003366'])
            del_label_color = 'blue'

        if dup_color == 'red':
            # From very light red to dark red (no pink)
            self.dup_cmap = LinearSegmentedColormap.from_list(
                'duplications', ['#FFCCCC', '#FF9999', '#FF6666', '#CC0000', '#660000'])
            dup_label_color = 'red'
        else:
            # From very pale blue to dark blue
            self.dup_cmap = LinearSegmentedColormap.from_list(
                'duplications', ['#E6F3FF', '#B3D9FF', '#66B2FF', '#0066CC', '#003366'])
            dup_label_color = 'blue'
        
        # Create data structure for plotting
        dat = pd.DataFrame({
            'chr': 'MT',
            'start': events['del_start_median'],
            'end': events['del_end_median'],
            'value': events['perc'],
            'dloop': events['dloop'],
            'delsize': events['delsize'],
            'final_event': events['final_event'],
            'group': events.get('group', 'G1')
        })
        
        # Order events by group for better visualization
        if 'group' in events.columns:
            dat['_sort_key'] = dat['group'].apply(lambda g: self.group_sort_key(g))
            dat = dat.sort_values(['_sort_key', 'value'], ascending=[True, False]).reset_index(drop=True)
            dat = dat.drop(columns=['_sort_key'])
        
        # Detect blacklist crossing
        dat['blacklist_crossing'] = False
        if blacklist_regions:
            for idx, row in dat.iterrows():
                dat.loc[idx, 'blacklist_crossing'] = crosses_blacklist(row['start'], row['end'], blacklist_regions)

        # Count events
        del_count = (dat['final_event'] == 'del').sum()
        dup_count = (dat['final_event'] == 'dup').sum()
        bl_del_count = ((dat['final_event'] == 'del') & dat['blacklist_crossing']).sum()
        bl_dup_count = ((dat['final_event'] == 'dup') & dat['blacklist_crossing']).sum()
        
        # Add degrees
        dat['deg1'] = self.config.degrees_per_genome * dat['start'] / self.genome_length
        dat['deg2'] = self.config.degrees_per_genome * dat['end'] / self.genome_length
        
        # Normalize angles into [0,360)
        for col in ('deg1', 'deg2'):
            dat[col] = pd.to_numeric(dat[col], errors='coerce')
            if dat[col].isna().any():
                logger.warning("%d invalid %s values", int(dat[col].isna().sum()), col)
            mask = dat[col].notna()
            if mask.any():
                dat.loc[mask, col] = (dat.loc[mask, col] % 360.0).astype(float)
        
        # CRITICAL FIX: For duplications that DON'T cross D-loop, swap coordinates
        # because the arc represents the PRESERVED region, not the duplicated region
        # This matches the original R script behavior
        
        dup_no_dloop_mask = (dat['final_event'] == 'dup') & (dat['dloop'] == 'no')
        if dup_no_dloop_mask.any():
            # Swap deg1 and deg2 for these events
            temp = dat.loc[dup_no_dloop_mask, 'deg1'].copy()
            dat.loc[dup_no_dloop_mask, 'deg1'] = dat.loc[dup_no_dloop_mask, 'deg2']
            dat.loc[dup_no_dloop_mask, 'deg2'] = temp
            logger.info(f"Swapped coordinates for {dup_no_dloop_mask.sum()} non-dloop duplications")
        
        # Separate data by type
        dat_del = dat[dat['final_event'] == 'del'].copy()
        dat_dup = dat[dat['final_event'] == 'dup'].copy()

        # Process duplication delsize
        if not dat_dup.empty:
            dat_dup['delsize'] = self.genome_length - dat_dup['delsize']

        # Calculate layout using simplified engine 
        logger.info(f"Using simplified layout engine for {len(dat_del)} dels, {len(dat_dup)} dups")

        layout_result = self.layout_engine.calculate_layout(
            dat_del,
            dat_dup,
            self.config.layout.total_radius
        )

        # Extract events with assigned radii from layout engine
        # Note: Margin enforcement is now handled in engine.py, not here
        all_events_with_radius = layout_result.events
        
        dat_del = all_events_with_radius[all_events_with_radius['final_event'] == 'del'].copy()
        dat_dup = all_events_with_radius[all_events_with_radius['final_event'] == 'dup'].copy()

        logger.info(f"Layout complete: {len(dat_del)} dels, {len(dat_dup)} dups")
        logger.info(f"  Del radius range: {layout_result.del_radius_range}")
        logger.info(f"  Dup radius range: {layout_result.dup_radius_range}")
        logger.info(f"  Blacklist radius: {layout_result.blacklist_radius:.1f}")

        dat_processed = pd.concat([dat_del, dat_dup], ignore_index=False)
        
        # Calculate color scales
        del_events = dat_processed[dat_processed['final_event'] == 'del']
        dup_events = dat_processed[dat_processed['final_event'] == 'dup']
        
        if scale == 'fixed':
            # Fixed scale: 0-100% range but with optimized gradient
            # 0-50% uses full color range for better visibility
            # 50-100% uses compressed dark shades (rare but theoretically possible)
            del_min, del_max = 0.0, 100.0
            dup_min, dup_max = 0.0, 100.0
            bl_min, bl_max = 0.0, 100.0
            logger.info("Using fixed heteroplasmy scale: 0-100% (gradient optimized for 0-50%)")
            
            # Store actual max values for legend indicator
            del_actual_max = del_events['value'].max() if not del_events.empty else 0
            dup_actual_max = dup_events['value'].max() if not dup_events.empty else 0
            if blacklist_regions:
                bl_events_temp = dat_processed[dat_processed['blacklist_crossing'] == True]
                bl_actual_max = bl_events_temp['value'].max() if not bl_events_temp.empty else 0
            else:
                bl_actual_max = 0
        else:  # dynamic
            # Dynamic scale: min-max within each category
            del_max = del_events['value'].max() if not del_events.empty else 0
            del_min = del_events['value'].min() if not del_events.empty else 0
            dup_max = dup_events['value'].max() if not dup_events.empty else 0  
            dup_min = dup_events['value'].min() if not dup_events.empty else 0
            
            # Calculate BL min/max for gradient coloring
            if blacklist_regions:
                bl_events = dat_processed[dat_processed['blacklist_crossing'] == True]
                if not bl_events.empty:
                    bl_min = bl_events['value'].min()
                    bl_max = bl_events['value'].max()
                else:
                    bl_min, bl_max = 0, 0
            else:
                bl_min, bl_max = 0, 0
            
            # For dynamic scale, actual max = scale max
            del_actual_max = del_max
            dup_actual_max = dup_max
            bl_actual_max = bl_max
            
            logger.info(f"Dynamic scale - Del: {del_min:.1f}-{del_max:.1f}%, "
                       f"Dup: {dup_min:.1f}-{dup_max:.1f}%")

        # Create figure
        if figsize is None:
            figsize = (self.config.figure_size, self.config.figure_size)
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='polar')
        ax.set_theta_zero_location('N')
        
        # Set direction based on parameter
        if direction == 'clockwise':
            ax.set_theta_direction(-1)  # Clockwise
        else:
            ax.set_theta_direction(1)   # Counterclockwise (default)
        
        # Calculate radius info
        blacklist_radius = layout_result.blacklist_radius
        genome_radius = self.config.layout.total_radius  # Genome circle at full boundary
        
        # Add offset for gene track and labels
        gene_track_offset = 35  # Space for gene track (15px track + margins)
        label_offset = 20       # Space for coordinate labels
        circle_offset = self.config.circle_offset + gene_track_offset + label_offset
        max_plot_radius = genome_radius + circle_offset
        
        # Draw genome circle at full radius boundary
        theta_genome = np.linspace(0, self.config.degrees_per_genome * np.pi / 180, 
                                  self.config.arc_resolution)
        radius_genome = [genome_radius] * len(theta_genome)
        ax.plot(theta_genome, radius_genome, 'k-', 
               linewidth=self.config.genome_circle_linewidth)
        
        # Draw separator circle ONLY if both deletions AND duplications exist
        # When only one event type exists, no separator is needed
        has_dels = not dat_del.empty and len(dat_del[dat_del['blacklist_crossing'] == False]) > 0
        has_dups = not dat_dup.empty and len(dat_dup[dat_dup['blacklist_crossing'] == False]) > 0
        
        if has_dels and has_dups:
            theta_sep = np.linspace(0, 2 * np.pi, 100)
            radius_sep = [blacklist_radius] * len(theta_sep)
            ax.plot(theta_sep, radius_sep, 'k--', 
                   linewidth=self.config.separator_circle_linewidth,
                   alpha=self.config.separator_circle_alpha)
        
        # Plot deletions (including blacklist events in green at their assigned radii)
        for idx, row in dat_del.iterrows():
            theta = np.radians(row['deg1'])
            theta2 = np.radians(row['deg2'])
            
            if theta2 < theta:
                theta2 += 2 * np.pi
            
            theta_arc = np.linspace(theta, theta2, self.config.arc_resolution)
            radius_arc = [row['radius']] * len(theta_arc)  # Use ASSIGNED radius only
            
            # Check if blacklist-crossing event - color GREEN, otherwise use normal color
            is_blacklist = row.get('blacklist_crossing', False)
            
            if is_blacklist:
                # Green gradient for blacklist events
                if bl_max > bl_min:
                    raw_norm = (row['value'] - bl_min) / (bl_max - bl_min)
                    # Apply gradient compression for fixed scale
                    if scale == 'fixed':
                        norm_value = self._compress_gradient_fixed_scale(row['value'], bl_min, bl_max)
                    else:
                        norm_value = raw_norm
                else:
                    norm_value = 0.5
                green_cmap = LinearSegmentedColormap.from_list(
                    'bl_green', ['#E6FFE6', '#B3FFB3', '#66FF66', '#32CD32', '#228B22'])
                color = green_cmap(norm_value)
            else:
                # Normal deletion color
                if del_max > del_min:
                    raw_norm = (row['value'] - del_min) / (del_max - del_min)
                    # Apply gradient compression for fixed scale
                    if scale == 'fixed':
                        norm_value = self._compress_gradient_fixed_scale(row['value'], del_min, del_max)
                    else:
                        norm_value = raw_norm
                else:
                    norm_value = 0.5
                color = self.del_cmap(norm_value)
            
            # Get alpha based on normalized value
            alpha = self.config.visualization.alpha_min + (norm_value * self.config.visualization.alpha_range)
            
            # Get line width
            linewidth = self._get_adaptive_linewidth(row['deg1'], dat_processed)
            
            # Draw arc at assigned radius (no overlay!)
            ax.plot(theta_arc, radius_arc, color=color, 
                   linewidth=linewidth, alpha=alpha, solid_capstyle='round', zorder=5)
        
        # Plot duplications (including blacklist events in green at their assigned radii)
        for idx, row in dat_dup.iterrows():
            theta = np.radians(row['deg1'])
            theta2 = np.radians(row['deg2'])
            
            if theta2 < theta:
                theta2 += 2 * np.pi
            
            theta_arc = np.linspace(theta, theta2, self.config.arc_resolution)
            radius_arc = [row['radius']] * len(theta_arc)  # Use ASSIGNED radius only
            
            # Check if blacklist-crossing event - color GREEN, otherwise use normal color
            is_blacklist = row.get('blacklist_crossing', False)
            
            if is_blacklist:
                # Green gradient for blacklist events
                if bl_max > bl_min:
                    raw_norm = (row['value'] - bl_min) / (bl_max - bl_min)
                    # Apply gradient compression for fixed scale
                    if scale == 'fixed':
                        norm_value = self._compress_gradient_fixed_scale(row['value'], bl_min, bl_max)
                    else:
                        norm_value = raw_norm
                else:
                    norm_value = 0.5
                green_cmap = LinearSegmentedColormap.from_list(
                    'bl_green', ['#E6FFE6', '#B3FFB3', '#66FF66', '#32CD32', '#228B22'])
                color = green_cmap(norm_value)
            else:
                # Normal duplication color
                if dup_max > dup_min:
                    raw_norm = (row['value'] - dup_min) / (dup_max - dup_min)
                    # Apply gradient compression for fixed scale
                    if scale == 'fixed':
                        norm_value = self._compress_gradient_fixed_scale(row['value'], dup_min, dup_max)
                    else:
                        norm_value = raw_norm
                else:
                    norm_value = 0.5
                color = self.dup_cmap(norm_value)
            
            # Get alpha based on normalized value
            alpha = self.config.visualization.alpha_min + (norm_value * self.config.visualization.alpha_range)
            
            # Get line width
            linewidth = self._get_adaptive_linewidth(row['deg1'], dat_processed)
            
            # Draw arc at assigned radius (no overlay!)
            ax.plot(theta_arc, radius_arc, color=color, 
                   linewidth=linewidth, alpha=alpha, solid_capstyle='round', zorder=5)
        
        # Initialize gene track variables for later use
        gene_track_outer = genome_radius  # Default if no genes
        
        # Add gene annotations if provided - with margin from genome circle
        if gene_annotations:
            gene_margin = 8  # Margin between genome circle and gene track
            gene_track_width = 15  # Track width
            gene_track_inner = genome_radius + gene_margin
            gene_track_outer = gene_track_inner + gene_track_width  # Outer edge of gene track
            gene_radius = gene_track_inner + (gene_track_width / 2)  # Center of track for drawing
            
            for gene in gene_annotations:
                # Handle both dict and GeneAnnotation object formats
                if hasattr(gene, 'start'):
                    # It's a GeneAnnotation object
                    start = gene.start
                    end = gene.end
                    color = gene.color
                    name = getattr(gene, 'name', None)
                elif isinstance(gene, dict):
                    # It's a dict
                    start = gene['start']
                    end = gene['end']
                    color = gene.get('color', (0, 0.5, 0))  # Default green
                    name = gene.get('name')
                else:
                    logger.warning(f"Unknown gene annotation format: {type(gene)}")
                    continue
                
                start_deg = self.config.degrees_per_genome * start / self.genome_length
                end_deg = self.config.degrees_per_genome * end / self.genome_length
                
                theta_start = np.radians(start_deg)
                theta_end = np.radians(end_deg)
                
                if theta_end < theta_start:
                    theta_end += 2 * np.pi
                
                # Draw gene bar - thicker for wider track
                theta_gene = np.linspace(theta_start, theta_end, 50)
                radius_gene = [gene_radius] * len(theta_gene)
                
                ax.plot(theta_gene, radius_gene, color=color, 
                       linewidth=gene_track_width, alpha=0.8, solid_capstyle='butt')
                
                # Add gene name - for larger genes only to avoid clutter
                # Skip short yellow tRNA genes for cleaner appearance
                gene_size_bp = end - start
                if name and gene_size_bp > 200 and color != (1.0, 1.0, 0.0):  # Only label genes >200bp, skip yellow tRNA
                    mid_theta = (theta_start + theta_end) / 2
                    # Position text just outside the gene track for better visibility
                    label_radius = float(gene_track_outer) + 12.0
                    
                    # Calculate rotation for text to be horizontal and readable
                    rotation_deg = np.degrees(mid_theta)
                    
                    # When in clockwise mode, mirror the angle
                    if direction == 'clockwise':
                        rotation_deg = 360 - rotation_deg
                    
                    # Flip text on bottom half to keep readable
                    if rotation_deg > 90 and rotation_deg < 270:
                        rotation_deg = rotation_deg + 180
                    
                    ax.text(mid_theta, label_radius, name,
                           rotation=rotation_deg, ha='center', va='center',
                           fontsize=7, weight='bold', color='black')
        
        # Add blacklist regions visualization - CLEAR ARC with EDGE LINES
        if blacklist_regions:
            logger.info(f"Drawing {len(blacklist_regions)} blacklist regions at radius {blacklist_radius:.1f}")
            # Find the actual min/max radius of events for proper edge line rendering
            min_event_radius = dat_processed['radius'].min() if not dat_processed.empty else 20
            max_event_radius = dat_processed['radius'].max() if not dat_processed.empty else genome_radius
            
            for bl in blacklist_regions:
                # Handle both dict and BlacklistRegion object formats
                if hasattr(bl, 'name'):
                    name = bl.name
                    start = bl.start
                    end = bl.end
                elif isinstance(bl, dict):
                    name = bl.get('name')
                    start = bl['start']
                    end = bl['end']
                else:
                    continue
                
                # Calculate angles
                start_deg = self.config.degrees_per_genome * start / self.genome_length
                end_deg = self.config.degrees_per_genome * end / self.genome_length
                theta_start = np.radians(start_deg)
                theta_end = np.radians(end_deg)
                
                if theta_end < theta_start:
                    theta_end += 2 * np.pi
                
                # REMOVED: Thick black arc on separator circle - too cluttering
                # Blacklist regions are only indicated by the dashed edge lines
                
                # Draw THIN EDGE LINES at both boundaries from config
                for theta in [theta_start, theta_end]:
                    # From minimum event radius to genome circle
                    ax.plot([theta, theta], [min_event_radius, genome_radius],
                           color='gray', linestyle='--', 
                           linewidth=self.config.layout.blacklist_boundary_linewidth, 
                           alpha=self.config.layout.blacklist_boundary_alpha, 
                           zorder=11)
                
                # Add label if name exists
                if name:
                    mid_pos = (start + end) / 2
                    mid_deg = self.config.degrees_per_genome * mid_pos / self.genome_length
                    theta = np.radians(mid_deg)
                    
                    # Draw connecting line from BL region to label
                    label_radius = blacklist_radius + 25
                    ax.plot([theta, theta], [blacklist_radius + 3, label_radius - 3],
                           'k-', linewidth=0.8, alpha=0.6, zorder=10)
                    
                    ax.text(theta, label_radius, name, 
                           ha='center', va='center', fontsize=9, weight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', 
                                   facecolor='white', edgecolor='black', 
                                   linewidth=1.5, alpha=1.0),
                           zorder=11)
        
        logger.info("=== LABEL CREATION ===")
        logger.info(f"All groups in processed data: {sorted(dat_processed['group'].unique())}")
        logger.info(f"Groups in bands: {[b.group_id for b in layout_result.group_bands]}")
        logger.info(f"Groups in singles: {[layout_result.events.iloc[s.event_index]['group'] for s in layout_result.single_events]}")

        # Add group labels with SYSTEMATIC collision avoidance
        # Include multi-event bands, single events, AND fallback for any unlabeled groups
        if layout_result.group_bands or layout_result.single_events or not dat_processed.empty:
            seen_groups = set()  # Track which groups we've already labeled
            label_positions = []
            
            # Add labels for multi-event group bands (ONE per group)
            for group_band in layout_result.group_bands:
                if group_band.n_events > 0:
                    # Skip if we already labeled this base group
                    base_group_id = group_band.group_id  # e.g., "G2" even if band is "G2a"
                    if base_group_id in seen_groups:
                        continue
                    seen_groups.add(base_group_id)
                    
                    # Smart label positioning with collision avoidance
                    # Try leftmost first, then next events if collision detected
                    candidate_indices = group_band.event_indices[:min(5, len(group_band.event_indices))]
                    
                    best_idx = None
                    min_conflict = float('inf')
                    
                    min_separation = self.config.layout.label_min_angular_separation
                    
                    for idx in candidate_indices:
                        event = layout_result.events.loc[idx]
                        label_deg = event['deg1']
                        
                        # Calculate total conflict score with existing labels
                        total_conflict = 0
                        for existing in label_positions:
                            # Calculate angular distance (handling wraparound)
                            ang_dist = abs(label_deg - existing['deg'])
                            if ang_dist > 180:
                                ang_dist = 360 - ang_dist
                            
                            # Score conflict based on proximity
                            if ang_dist < min_separation:
                                # Too close - high penalty
                                total_conflict += (min_separation - ang_dist) * 10
                        
                        # Choose event with minimum conflict
                        if total_conflict < min_conflict:
                            min_conflict = total_conflict
                            best_idx = idx
                    
                    # Use the best position found
                    leftmost_event = layout_result.events.loc[best_idx]
                    
                    # Get position and radius from engine-assigned data
                    label_deg = leftmost_event['deg1']
                    event_radius = leftmost_event['radius']
                    event_type = leftmost_event['final_event']
                    
                    # LIGHTWEIGHT COLLISION RESOLUTION: If still conflicts after trying all candidates,
                    # apply small radial offset (nudge inward) to create vertical separation
                    radial_offset = 0
                    if min_conflict > 0:  # Still has conflict
                        # Check if there's a nearby label (within min_separation)
                        for existing in label_positions:
                            ang_dist = abs(label_deg - existing['deg'])
                            if ang_dist > 180:
                                ang_dist = 360 - ang_dist
                            
                            if ang_dist < min_separation:
                                # Apply radial nudge: alternate between inward/outward
                                # Use group index to alternate direction
                                nudge_direction = 1 if len(label_positions) % 2 == 0 else -1
                                radial_offset = nudge_direction * self.config.layout.label_radial_nudge
                                logger.info(f"  Applying radial nudge {radial_offset}px to {base_group_id} "
                                          f"(conflict with nearby label at {ang_dist:.1f}°)")
                                break
                    
                    # Label color based on event TYPE (red=del, blue=dup)
                    label_color = del_label_color if event_type == 'del' else dup_label_color
                    
                    label_positions.append({
                        'deg': label_deg,
                        'radius': event_radius,
                        'radial_offset': radial_offset,  # Store offset for rendering
                        'label': base_group_id,  # Use base group ID only
                        'color': label_color
                    })
                    logger.info(f"  Band label {base_group_id}: deg={label_deg:.1f}, radius={event_radius:.1f}, type={event_type}, conflict={min_conflict:.1f}")
            
            # Add labels for single-event groups (if not already labeled)
            for single_event in layout_result.single_events:
                event = layout_result.events.iloc[single_event.event_index]
                group_id = event['group']
                
                # Skip if already labeled
                if group_id in seen_groups:
                    continue
                seen_groups.add(group_id)
                
                # Label color based on event TYPE (red=del, blue=dup)
                label_color = del_label_color if event['final_event'] == 'del' else dup_label_color
                
                label_positions.append({
                    'deg': event['deg1'],
                    'radius': single_event.radius,
                    'radial_offset': 0,  # No offset for single events (no conflict expected)
                    'label': group_id,
                    'color': label_color
                })
                logger.info(f"  Single label {group_id}: deg={event['deg1']:.1f}, radius={single_event.radius:.1f}")

            # FALLBACK: Label any remaining groups that weren't in bands or single_events
            all_groups = layout_result.events['group'].unique()
            for group_id in all_groups:
                if group_id in seen_groups:
                    continue
                
                # Find representative event for this group from layout_result.events
                group_events = layout_result.events[layout_result.events['group'] == group_id]
                if not group_events.empty:
                    # Use leftmost event
                    leftmost_idx = group_events['deg1'].idxmin()
                    leftmost_event = group_events.loc[leftmost_idx]
                    
                    # Label color based on event TYPE (red=del, blue=dup)
                    label_color = del_label_color if leftmost_event['final_event'] == 'del' else dup_label_color
                    
                    seen_groups.add(group_id)
                    label_positions.append({
                        'deg': leftmost_event['deg1'],
                        'radius': leftmost_event['radius'],
                        'radial_offset': 0,  # No offset for fallback
                        'label': group_id,
                        'color': label_color
                    })
                    logger.info(f"Added fallback label for group {group_id} (not in layout bands)")
            
            # Log all labels being drawn
            logger.info(f"Drawing {len(label_positions)} group labels: {[p['label'] for p in label_positions]}")
            
            # Third pass: Draw labels at adjusted positions with config styling
            # HIGH Z-ORDER ensures labels draw on top of all arcs and blacklist lines
            logger.info(f"TOTAL LABELS TO DRAW: {len(label_positions)}")
            for pos in label_positions:
                logger.info(f"  {pos['label']}: deg={pos['deg']:.2f}, radius={pos['radius']:.2f}, color={pos['color']}")
                
            for pos in label_positions:
                theta = np.radians(pos['deg'])
                event_radius = pos['radius']
                radial_offset = pos.get('radial_offset', 0)  # Get offset if present
                
                # Keep labels inside plot area with margin
                label_offset = 12
                max_label_radius = genome_radius - 20  # 20px margin from edge
                
                # Apply radial nudge if there was a conflict
                label_radius = min(event_radius + label_offset + radial_offset, max_label_radius)
                
                # Draw thin pointer line from event to label using config
                ax.plot([theta, theta], 
                       [event_radius + 1, label_radius - 3],
                       color='gray', 
                       linewidth=self.config.layout.label_connector_linewidth, 
                       alpha=0.6, zorder=100)  # High z-order
                
                # Draw small marker at event position
                ax.plot(theta, event_radius, 'o', 
                       markersize=1.5, color='gray', alpha=0.8, zorder=101)
                
                # Add label with clean colored box - using config parameters
                ax.text(theta, label_radius, pos['label'],
                       ha='center', va='center', 
                       fontsize=self.config.layout.label_fontsize,
                       weight='normal',
                       color=pos['color'],
                       bbox=dict(boxstyle=f'round,pad={self.config.layout.label_box_padding}',
                               facecolor='white', edgecolor=pos['color'], 
                               linewidth=self.config.layout.label_box_linewidth, 
                               alpha=1.0),
                       zorder=102)  # Highest z-order - on top of everything
        
        # Configure axes
        ax.set_ylim(0, max_plot_radius)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)
        
        # Add position markers - adjusted for gene track spacing
        positions = np.arange(0, self.genome_length, 1000)
        for pos in positions:
            deg = np.radians(self.config.degrees_per_genome * pos / self.genome_length)
            # Position markers outside genes if present
            marker_start = gene_track_outer + 5 if gene_annotations else genome_radius + 5
            text_offset = 18 if gene_annotations else 12  # More space when genes present
            
            # Small tick mark
            ax.plot([deg, deg], [marker_start, marker_start + 5], 
                   color='gray', linewidth=1, alpha=0.7)
            
            # Position label (kb)
            ax.text(deg, marker_start + text_offset, f'{pos//1000}', 
                   ha='center', va='center', fontsize=9, color='gray')
        
        # Add title - position from config for easy adjustment
        if title is None:
            bl_text = f" (BL: {len(blacklist_regions)} regions)" if blacklist_regions else ""
            title = f"Mitochondrial DNA Structural Alterations{bl_text}"
        
        # Title at top - y position configurable
        fig.text(0.5, self.config.title_y_position, title, 
                ha='center', va='top', fontsize=self.config.title_fontsize, 
                weight='bold', transform=fig.transFigure)
        
        # Add count box - centered at same X as other legends for perfect alignment
        count_text_lines = [f"Del: {del_count}  Dup: {dup_count}"]
        if blacklist_regions and (bl_del_count > 0 or bl_dup_count > 0):
            count_text_lines.append(f"BL Del: {bl_del_count}  BL Dup: {bl_dup_count}")
        
        # Tight, compact box sizing
        if len(count_text_lines) == 2:
            # Two lines - compact vertical spacing
            box_width = 0.14  # Fixed width for consistency
            box_height = 0.052  # Compact height for 2 lines
            box_y_pos = 0.78  # Position
            line_spacing = 0.022  # Tight line spacing
        else:
            # Single line - even more compact
            box_width = 0.12
            box_height = 0.03
            box_y_pos = 0.80
            line_spacing = 0
        
        # Center box at same x as other legends (0.12)
        count_box_x = 0.12 - box_width / 2  # Center box at 0.12
        
        # Draw count box with minimal padding
        count_box = FancyBboxPatch(
            (count_box_x, box_y_pos),
            box_width, box_height,
            boxstyle="round,pad=0.004",  # Minimal padding
            facecolor='white',
            edgecolor='gray',
            linewidth=1.2,
            transform=fig.transFigure,
            zorder=100
        )
        fig.patches.append(count_box)
        
        # Add count text with tight spacing - centered at 0.12
        if len(count_text_lines) == 2:
            # Two lines - centered vertically in box
            y_start = box_y_pos + box_height * 0.68
            for i, line in enumerate(count_text_lines):
                text_color = 'black' if i == 0 else 'limegreen'
                fig.text(0.12, y_start - (i * line_spacing), line,  # Centered at 0.12
                        fontsize=12, weight='bold', color=text_color,
                        ha='center', va='center', transform=fig.transFigure, zorder=101)
        else:
            # Single line - centered at 0.12
            fig.text(0.12, box_y_pos + box_height/2, count_text_lines[0],  # Centered at 0.12
                    fontsize=12, weight='bold', color='black',
                    ha='center', va='center', transform=fig.transFigure, zorder=101)
        
        # Add heteroplasmy scale legend - PROPERLY CENTERED and LARGER
        # Calculate legend positioning dynamically for better centering
        legend_x = 0.12  # Base X position for legend group
        legend_y = 0.30  # Base Y position (raised for better spacing)
        legend_height = 0.22  # Taller bars
        bar_width = 0.020  # Slightly wider bars for visibility
        bar_gap = 0.024  # Better spacing between bars
        
        # Title for heteroplasmy legend - LARGER and bold
        fig.text(legend_x, legend_y + legend_height + 0.045, "Heteroplasmy (%)",
                fontsize=14, weight='bold', ha='center', va='bottom', 
                transform=fig.transFigure)
        
        # Determine which scale bars to show
        scale_types = []
        if not del_events.empty:
            # Filter out blacklist-crossing deletions for the Del bar
            del_non_bl = del_events[del_events.get('blacklist_crossing', False) == False]
            if not del_non_bl.empty:
                scale_types.append(('Del', del_label_color, self.del_cmap, del_min, del_max))
        
        if not dup_events.empty:
            # Filter out blacklist-crossing duplications for the Dup bar  
            dup_non_bl = dup_events[dup_events.get('blacklist_crossing', False) == False]
            if not dup_non_bl.empty:
                scale_types.append(('Dup', dup_label_color, self.dup_cmap, dup_min, dup_max))
        
        if blacklist_regions and bl_max > bl_min:
            # Use proper lime green gradient for BL
            scale_types.append(('BL', 'limegreen', LinearSegmentedColormap.from_list(
                'blacklist', ['#E6FFE6', '#B3FFB3', '#66FF66', '#32CD32', '#228B22']), bl_min, bl_max))
        
        # Calculate total width and starting position for better centering
        total_bars = len(scale_types)
        if total_bars > 0:
            total_width = total_bars * bar_width + (total_bars - 1) * bar_gap
            start_x = legend_x - total_width / 2
            
            # Draw gradient bars
            for i, (label, color, cmap, vmin, vmax) in enumerate(scale_types):
                bar_x = start_x + i * (bar_width + bar_gap)
                
                # Get actual max for this type
                if label == 'Del':
                    actual_max = del_actual_max
                elif label == 'Dup':
                    actual_max = dup_actual_max
                elif label == 'BL':
                    actual_max = bl_actual_max
                else:
                    actual_max = vmax
                
                # Create gradient using rectangles
                n_colors = 50  # Smooth gradient
                for j in range(n_colors):
                    y_frac = j / n_colors
                    y = legend_y + y_frac * legend_height
                    height = legend_height / n_colors
                    
                    # Color based on position in gradient
                    # For fixed scale, apply the same compression used in plotting
                    if scale == 'fixed':
                        # Calculate the actual value this y position represents
                        value_at_y = vmin + y_frac * (vmax - vmin)
                        # Apply compression function to get colormap lookup
                        norm_val = self._compress_gradient_fixed_scale(value_at_y, vmin, vmax)
                    else:
                        # Dynamic scale uses linear gradient
                        norm_val = j / (n_colors - 1)
                    
                    bar_color = cmap(norm_val)
                    
                    rect = plt.Rectangle((bar_x, y), bar_width, height,
                                        facecolor=bar_color, transform=fig.transFigure,
                                        linewidth=0)
                    fig.patches.append(rect)
                
                # Determine precision based on event type (ONCE for all labels)
                # BL events are typically low heteroplasmy (artifacts) - need 2 decimals
                # Del/Dup events use 1 decimal
                decimal_places = 2 if label == 'BL' else 1
                
                # Add dashed line at actual max value (if different from scale max)
                if scale == 'fixed' and actual_max > 0 and actual_max < vmax:
                    # Calculate y position for actual max
                    actual_frac = actual_max / vmax
                    actual_y = legend_y + actual_frac * legend_height
                    
                    # Draw dashed line across the bar
                    line_x = [bar_x - 0.003, bar_x + bar_width + 0.003]
                    line_y = [actual_y, actual_y]
                    
                    from matplotlib.lines import Line2D
                    line = Line2D(line_x, line_y, color='black', linewidth=1.5, 
                                 linestyle='--', alpha=0.8, transform=fig.transFigure, zorder=10)
                    fig.lines.append(line)
                    
                    # Add label for actual max (small, to the right) with appropriate precision
                    fig.text(bar_x + bar_width + 0.006, actual_y, f"{actual_max:.{decimal_places}f}",
                            fontsize=8, ha='left', va='center', transform=fig.transFigure,
                            color='black', weight='bold')
                
                # NO border around gradient - cleaner appearance
                
                # Add min/max labels with appropriate precision
                fig.text(bar_x + bar_width/2, legend_y - 0.01, f"{vmin:.{decimal_places}f}",
                        fontsize=11, ha='center', va='top', transform=fig.transFigure)
                fig.text(bar_x + bar_width/2, legend_y + legend_height + 0.01, f"{vmax:.{decimal_places}f}",
                        fontsize=11, ha='center', va='bottom', transform=fig.transFigure)
                
                # Add type label with LARGER font (12pt)
                fig.text(bar_x + bar_width/2, legend_y - 0.03, label,
                        fontsize=12, color=color, weight='bold', ha='center', va='top',
                        transform=fig.transFigure)
        
        # Add metadata if provided
        if metadata_text:
            fig.text(0.5, 0.02, metadata_text,
                    fontsize=8, ha='center', va='bottom', transform=fig.transFigure,
                    style='italic', wrap=True)
        
        # Add genes legend - LEFT side, properly CENTERED with LARGER fonts
        if gene_annotations:
            genes_x = 0.12  # Same x as other legends for perfect alignment
            genes_y = 0.16  # Below heteroplasmy scale with proper spacing
            
            # Title - LARGER and bold
            fig.text(genes_x, genes_y + 0.015, "Genes",
                    fontsize=14, weight='bold', ha='center', va='bottom',
                    transform=fig.transFigure)
            
            # Gene type legend with better spacing and larger fonts
            gene_types = [
                ('Protein-coding', (0.0, 0.5, 0.0)),  # Green
                ('tRNA', (1.0, 1.0, 0.0)),  # Yellow
                ('rRNA', (0.5, 0.0, 0.5))   # Purple
            ]
            
            # Center the legend items properly
            legend_box_width = 0.15  # Slightly wider for better appearance
            start_x = genes_x - legend_box_width/2 + 0.01
            
            for i, (gene_type, color) in enumerate(gene_types):
                y_pos = genes_y - 0.020 - i * 0.030  # Better vertical spacing
                
                # Color box - NO border for cleaner appearance
                rect = plt.Rectangle((start_x, y_pos - 0.010), 0.032, 0.020,
                                    facecolor=color, edgecolor='none', linewidth=0,
                                    transform=fig.transFigure)
                fig.patches.append(rect)
                
                # Label - LARGER font with proper spacing
                fig.text(start_x + 0.038, y_pos, gene_type,
                        fontsize=11, va='center', ha='left', 
                        transform=fig.transFigure)
        
        plt.tight_layout()
        
        # Save figure
        fig.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        logger.info(f"Plot saved to {output_file}")
        
        if show:
            plt.show()
        
        return fig