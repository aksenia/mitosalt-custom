"""
Type definitions for SaltShaker

Common types used throughout the package for type checking and documentation.
"""

from __future__ import annotations
from typing import TypedDict, Literal, List, Dict, Union, Tuple, TYPE_CHECKING
from pathlib import Path

if TYPE_CHECKING:
    import pandas as pd

# Type aliases
PathLike = Union[str, Path]
"""File path as string or Path object"""

EventType = Literal['del', 'dup']
"""Type of structural variant event"""

ClassificationType = Literal['Single', 'Multiple', 'Blacklist-only', 'No significant events', 'Unknown']
"""Classification result for event pattern"""

DLoopStatus = Literal['yes', 'no']
"""Whether event crosses D-loop region"""

GeneType = Literal['protein_coding', 'tRNA', 'rRNA']
"""Type of mitochondrial gene"""

# Complex return types
OriginTuple = Tuple[int, int]
"""Origin region as (start, end) tuple"""

ClassificationResultTuple = Tuple[ClassificationType, str, Dict[str, Union[int, float, bool, str, List]], 'pd.DataFrame']
"""Full classification result: (classification, reason, criteria, events_with_groups)"""

IntermediateReadResult = Tuple['pd.DataFrame', int]
"""Result from reading intermediate format: (events_df, genome_length)"""


# Structured data types

class BlacklistRegion(TypedDict):
    """Genomic region to exclude from analysis"""
    chr: str
    start: int
    end: int


class GeneAnnotation(TypedDict):
    """Gene annotation from BED file"""
    chr: str
    start: int
    end: int
    name: str
    strand: str  # '+' or '-'
    color: Tuple[float, float, float]  # RGB tuple (normalized 0-1)


class EventRecord(TypedDict, total=False):
    """
    Single mitochondrial structural variant event
    
    Required fields marked with total=False to allow partial dicts.
    In practice, events from EventCaller will have all core fields.
    """
    # Core identification
    cluster: str
    sample: str
    
    # Position and size
    del_start_median: float
    del_end_median: float
    delsize: int
    
    # Heteroplasmy
    perc: float  # Heteroplasmy percentage (0-100)
    nread: int  # Alternative reads
    tread: int  # Total reads
    
    # Event classification
    final_event: EventType
    dloop: DLoopStatus
    
    # Spatial grouping (added by classifier)
    group: str  # e.g., 'G1', 'G2', 'BL1'
    blacklist_crossing: bool
    
    # Calculated fields (added by call.py)
    final_event_size: int
    final_start: float
    final_end: float
    
    # Flanking sequences (added by EventCaller)
    left_flank: str
    right_flank: str
    
    # Degrees for plotting (added by visualizer)
    deg1: float
    deg2: float
    radius: float


class ClassificationResult(TypedDict):
    """Result from event pattern classification"""
    classification: ClassificationType
    reason: str
    criteria: Dict[str, Union[int, float, bool, str]]
    # Note: events_with_groups would be pd.DataFrame, typed separately


class OriginRegion(TypedDict):
    """Replication origin region coordinates"""
    start: int
    end: int


class GroupInfo(TypedDict):
    """Spatial group information from SpatialGroupAnalyzer"""
    id: int
    group_id: str  # e.g., 'G1', 'G2'
    event_count: int
    event_type: EventType
    max_heteroplasmy: float
    mean_heteroplasmy: float
    total_heteroplasmy: float
    median_size: float
    max_size: float
    spatial_range: float
    high_het_count: int
    significant_count: int
    dominance_score: float
    representative: Dict[str, Union[float, str]]  # Event with highest heteroplasmy
    events: List[Dict[str, Union[int, float, str]]]  # Individual events in group


# Configuration types

class ClassificationConfigDict(TypedDict, total=False):
    """Configuration for event classification"""
    HIGH_HET_THRESHOLD: float
    NOISE_THRESHOLD: float
    CLUSTER_RADIUS: int
    MIN_CLUSTER_SIZE: int
    MULTIPLE_EVENT_THRESHOLD: int
    DOMINANT_GROUP_FRACTION: float