"""
Layout types for SaltShaker
Data structures for layout engine results

All types are immutable (frozen) for safety and testability.
"""
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Literal, Any
import pandas as pd


@dataclass(frozen=True)
class SectorBudget:
    """
    Radial budget allocation for a single sector
    
    Attributes:
        sector_id: Sector index (0-11 for 12 sectors)
        start_deg: Starting angle in degrees
        end_deg: Ending angle in degrees
        event_count: Number of events in this sector
        initial_budget: Initially allocated radial space (px)
        final_budget: Final radial space after dynamic expansion (px)
        expansion_iterations: Number of times budget was expanded
    """
    sector_id: int
    start_deg: float
    end_deg: float
    event_count: int
    initial_budget: float
    final_budget: float
    expansion_iterations: int = 0
    
    @property
    def angular_width(self) -> float:
        """Angular width of sector in degrees"""
        return self.end_deg - self.start_deg
    
    @property
    def was_expanded(self) -> bool:
        """Whether sector budget was dynamically expanded"""
        return self.final_budget > self.initial_budget


@dataclass(frozen=True)
class GroupBandLayout:
    """
    Layout result for a single group's band
    
    Represents a group of events positioned in a radial band.
    All events in the band have consistent spacing.
    
    Attributes:
        group_id: Group identifier (e.g., 'G2', 'G5', 'BL1')
        sub_band_index: Sub-band index if group was split (0 for unsplit)
        band_top: Outer radius of band (px)
        band_bottom: Inner radius of band (px)
        n_events: Number of events in this band
        event_indices: DataFrame indices of events in this band
        within_spacing: Consistent spacing between events (px)
        sector_id: Sector this band belongs to
    """
    group_id: str
    sub_band_index: int
    band_top: float
    band_bottom: float
    n_events: int
    event_indices: List[int]
    within_spacing: float
    sector_id: int
    
    @property
    def band_height(self) -> float:
        """Height of band in pixels"""
        return self.band_top - self.band_bottom
    
    @property
    def is_sub_band(self) -> bool:
        """Whether this is a sub-band (group was split)"""
        return self.sub_band_index > 0
    
    @property
    def display_id(self) -> str:
        """Display ID for labeling (e.g., 'G2', 'G2a', 'G2b')"""
        if self.is_sub_band:
            # Convert 0,1,2 -> a,b,c
            suffix = chr(ord('a') + self.sub_band_index)
            return f"{self.group_id}{suffix}"
        return self.group_id


@dataclass(frozen=True)
class SingleEventLayout:
    """
    Layout result for a single event on shared radius
    
    Attributes:
        event_index: DataFrame index of event
        radius: Assigned radius (px)
        sector_id: Sector this event belongs to
        shared_level_id: ID of shared radius level
        n_events_on_level: Total events sharing this radius
    """
    event_index: int
    radius: float
    sector_id: int
    shared_level_id: int
    n_events_on_level: int


@dataclass
class LayoutResult:
    """
    Complete layout solution for all events
    
    Contains all information needed for rendering:
    - Event positions (radii)
    - Sector budgets
    - Group band layouts
    - Single event layouts
    
    This is the output of LayoutEngine and input to PlotRenderer.
    
    Attributes:
        events: DataFrame with all events (includes 'radius' column)
        sector_budgets: Budget allocation for each sector
        group_bands: Layout info for multi-event group bands
        single_events: Layout info for single events
        total_radius_used: Maximum radius used
        del_radius_range: (min, max) radius for deletions
        dup_radius_range: (min, max) radius for duplications
        blacklist_radius: Radius of del/dup separator
        layout_algorithm: Algorithm used ('legacy' or 'smart_dynamic')
        layout_stats: Statistics about the layout
    """
    events: pd.DataFrame
    sector_budgets: Dict[int, SectorBudget]
    group_bands: List[GroupBandLayout]
    single_events: List[SingleEventLayout]
    total_radius_used: float
    del_radius_range: Tuple[float, float]
    dup_radius_range: Tuple[float, float]
    blacklist_radius: float
    layout_algorithm: str
    layout_stats: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def n_sectors(self) -> int:
        """Number of sectors"""
        return len(self.sector_budgets)
    
    @property
    def n_group_bands(self) -> int:
        """Number of group bands (including sub-bands)"""
        return len(self.group_bands)
    
    @property
    def n_single_events(self) -> int:
        """Number of single events"""
        return len(self.single_events)
    
    @property
    def total_events(self) -> int:
        """Total number of events"""
        return len(self.events)
    
    def get_events_in_band(self, band: GroupBandLayout) -> pd.DataFrame:
        """Get events DataFrame for a specific band"""
        return self.events.iloc[band.event_indices]
    
    def get_band_for_group(self, group_id: str) -> List[GroupBandLayout]:
        """Get all bands for a group (may be split)"""
        return [b for b in self.group_bands if b.group_id == group_id]
    
    def get_sector_budget(self, sector_id: int) -> SectorBudget:
        """Get budget for a sector"""
        return self.sector_budgets[sector_id]


@dataclass(frozen=True)
class BlacklistRegion:
    """
    Blacklist region definition
    
    Attributes:
        start: Start position (bp)
        end: End position (bp)
        name: Optional region name
    """
    start: int
    end: int
    name: Optional[str] = None
    
    @property
    def size(self) -> int:
        """Size of region in bp"""
        return self.end - self.start


@dataclass(frozen=True)
class GeneAnnotation:
    """
    Gene annotation for display
    
    Attributes:
        name: Gene name (e.g., 'MT-CYB', 'MT-ND1')
        start: Start position (bp)
        end: End position (bp)
        gene_type: Type of gene
        color: RGB color tuple for display
    """
    name: str
    start: int
    end: int
    gene_type: Literal['protein-coding', 'tRNA', 'rRNA']
    color: Tuple[float, float, float]
    
    @property
    def size(self) -> int:
        """Size of gene in bp"""
        return self.end - self.start
    
    @property
    def midpoint(self) -> float:
        """Midpoint position (bp)"""
        return (self.start + self.end) / 2