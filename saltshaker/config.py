"""
SaltShaker Configuration - SIMPLIFIED VERSION
Clean configuration with only necessary parameters
NO LEGACY CODE - NO HARDCODED VALUES
"""
from dataclasses import dataclass, field
from typing import Tuple, List, Literal, Optional


@dataclass
class ClassificationConfig:
    """
    Classification thresholds for Single vs Multiple patterns
    Based on Basu et al. PLoS Genetics 2020
    """
    
    HIGH_HET_THRESHOLD: float = 10.0
    """High heteroplasmy threshold - pathogenic significance (≥10%)"""
    
    NOISE_THRESHOLD: float = 0.3
    """Noise threshold - below this likely artifacts (<0.3%)"""
    
    CLUSTER_RADIUS: int = 600
    """Spatial grouping radius in base pairs"""
    
    MIN_CLUSTER_SIZE: int = 2
    """Minimum events to form a spatial cluster"""
    
    MULTIPLE_EVENT_THRESHOLD: int = 10
    """Event count threshold for Multiple pattern (>10 = Multiple)"""
    
    DOMINANT_GROUP_FRACTION: float = 0.5
    """Fraction of events in dominant group for Single pattern (≥50%)"""


@dataclass
class LayoutConfig:
    """
    Simplified layout configuration for radial positioning
    
    Uses simple radial stacking without sectors or complex budgets.
    """
    
    # ============================================================
    # SPACE ALLOCATION
    # ============================================================
    total_radius: int = 400
    """Total radius available for event placement (px)"""
    
    separator_fraction: float = 0.15
    """Fraction of total radius for del/dup separator band"""
    
    # ============================================================
    # GROUP BAND SIZING
    # ============================================================
    base_band_size: int = 25
    """Base band size for multi-event groups (px)"""
    
    band_size_scale_factor: float = 0.8
    """Linear scaling factor: pixels per event"""
    
    max_band_size: int = 200
    """Maximum band size for very large groups (px)"""
    
    # ============================================================
    # SPACING
    # ============================================================
    min_event_spacing: float = 3.0
    """Minimum spacing between event centers within a group (px)"""
    
    uniform_within_group_spacing: bool = True
    """Use uniform dense spacing within all groups (True) or proportional spacing (False)"""
    
    target_event_spacing: float = 3.5
    """Target spacing between events within a group when uniform spacing enabled (px)"""
    
    group_gap: int = 6
    """Vertical gap between adjacent group bands (px)"""
    
    # ============================================================
    # BLACKLIST VISUALIZATION
    # ============================================================
    blacklist_boundary_linewidth: float = 0.8
    """Line width for dashed boundary lines at blacklist edges (px)"""
    
    blacklist_boundary_alpha: float = 0.4
    """Transparency for blacklist boundary lines"""
    
    # ============================================================
    # LABEL STYLING
    # ============================================================
    label_fontsize: int = 7
    """Font size for group labels"""
    
    label_box_padding: float = 0.15
    """Padding inside label boxes"""
    
    label_box_linewidth: float = 0.8
    """Border width for label boxes (px)"""
    
    label_connector_linewidth: float = 0.4
    """Line width for label connector lines (px)"""
    
    label_min_angular_separation: float = 12.0
    """Minimum angular separation between labels (degrees) to avoid overlap"""
    
    label_radial_nudge: float = 8.0
    """Radial offset applied to conflicting labels (px) for vertical separation"""


@dataclass
class VisualizationConfig:
    """
    Adaptive visualization parameters
    Controls line width and transparency based on density
    """
    
    # ============================================================
    # ADAPTIVE LINE WIDTH
    # ============================================================
    adaptive_linewidth_enabled: bool = True
    """Enable adaptive line width based on local density"""
    
    density_window_degrees: float = 30.0
    """Angular window for calculating local density (degrees)"""
    
    density_linewidth_min: float = 1.5
    """Minimum line width in very dense areas (px)"""
    
    density_linewidth_max: float = 3.5
    """Maximum line width in sparse areas (px)"""
    
    density_threshold_high: int = 20
    """Event count threshold for 'high density' in window"""
    
    density_threshold_medium: int = 10
    """Event count threshold for 'medium density' in window"""
    
    # ============================================================
    # FIXED LINE WIDTH (Fallback)
    # ============================================================
    event_linewidth_small: float = 3.0
    """Line width when ≤150 events (used when adaptive disabled)"""
    
    event_linewidth_medium: float = 2.5
    """Line width when ≤300 events (used when adaptive disabled)"""
    
    event_linewidth_large: float = 2.0
    """Line width when >300 events (used when adaptive disabled)"""
    
    event_linewidth_threshold_medium: int = 150
    """Event count threshold for medium line width"""
    
    event_linewidth_threshold_large: int = 300
    """Event count threshold for large line width"""
    
    # ============================================================
    # TRANSPARENCY
    # ============================================================
    alpha_min: float = 0.70
    """Minimum alpha (transparency) for low heteroplasmy events"""
    
    alpha_range: float = 0.25
    """Alpha range from min to max (max = alpha_min + alpha_range)"""


@dataclass
class PlotConfig:
    """
    Complete plot configuration - SIMPLIFIED
    
    Clean configuration without legacy parameters.
    """
    
    # ============================================================
    # SUB-CONFIGURATIONS
    # ============================================================
    layout: LayoutConfig = field(default_factory=LayoutConfig)
    """Layout configuration"""
    
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    """Visualization configuration"""
    
    # ============================================================
    # RENDERING PARAMETERS
    # ============================================================
    arc_resolution: int = 100
    """Number of points for drawing event arcs (smoothness)"""
    
    # ============================================================
    # GENOME CIRCLE
    # ============================================================
    degrees_per_genome: int = 358
    """Total degrees around circle (360° leaves 2° gap at top)"""
    
    genome_circle_linewidth: int = 3
    """Line width for main genome circle (px)"""
    
    separator_circle_linewidth: int = 2
    """Line width for del/dup separator circle (px)"""
    
    separator_circle_alpha: float = 0.7
    """Transparency for separator circle"""
    
    # ============================================================
    # FIGURE SETTINGS
    # ============================================================
    figure_size: float = 12.0
    """Figure size in inches (square plot)"""
    
    dpi: int = 300
    """DPI for saved figures"""
    
    title_fontsize: int = 14
    """Font size for main title"""
    
    title_y_position: float = 0.98
    """Vertical position of title in figure coordinates (0-1, where 1 is top)"""
    
    circle_offset: float = 20.0
    """Offset from outermost radius to figure edge (px)"""
    
    # ============================================================
    # PRESET CONFIGURATIONS
    # ============================================================
    
    @classmethod
    def publication(cls) -> 'PlotConfig':
        """
        High-quality settings for publication figures
        
        - 600 DPI
        - Larger figure size (14 inches)
        - Thicker lines for better visibility in print
        
        Example:
            >>> config = PlotConfig.publication()
            >>> plotter = CircularPlotter(16569, config)
        """
        config = cls()
        config.dpi = 600
        config.figure_size = 14.0
        config.genome_circle_linewidth = 4
        config.separator_circle_linewidth = 3
        config.visualization.density_linewidth_min = 2.0
        config.visualization.density_linewidth_max = 4.0
        return config
    
    @classmethod
    def presentation(cls) -> 'PlotConfig':
        """
        Settings optimized for presentations
        
        - Lower DPI (150) for smaller file size
        - Larger fonts and thicker lines for screen viewing
        - Higher contrast
        
        Example:
            >>> config = PlotConfig.presentation()
            >>> plotter = CircularPlotter(16569, config)
        """
        config = cls()
        config.dpi = 150
        config.figure_size = 10.0
        config.title_fontsize = 16
        config.genome_circle_linewidth = 5
        config.separator_circle_linewidth = 4
        config.visualization.density_linewidth_min = 2.5
        config.visualization.density_linewidth_max = 5.0
        config.visualization.alpha_min = 0.85
        return config
    
    @classmethod
    def compact(cls) -> 'PlotConfig':
        """
        Compact settings for many events
        
        - Smaller bands and gaps
        - Thinner lines
        - More events per band
        
        Example:
            >>> config = PlotConfig.compact()
            >>> plotter = CircularPlotter(16569, config)
        """
        config = cls()
        config.layout.base_band_size = 15
        config.layout.band_size_scale_factor = 0.5
        config.layout.group_gap = 3
        config.layout.min_event_spacing = 2.0
        config.visualization.density_linewidth_min = 1.0
        config.visualization.density_linewidth_max = 2.5
        return config
    
    @classmethod
    def debug(cls) -> 'PlotConfig':
        """
        Settings for debugging layout issues
        
        - Larger spacing for clarity
        - Maximum transparency
        - Thicker lines
        
        Example:
            >>> config = PlotConfig.debug()
            >>> plotter = CircularPlotter(16569, config)
        """
        config = cls()
        config.layout.min_event_spacing = 5.0
        config.layout.group_gap = 10
        config.visualization.alpha_min = 1.0
        config.visualization.alpha_range = 0.0
        config.visualization.density_linewidth_min = 3.0
        config.visualization.density_linewidth_max = 5.0
        return config