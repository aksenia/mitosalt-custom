"""
SaltShaker Configuration
Constants and thresholds for mitochondrial event classification
"""
from dataclasses import dataclass


@dataclass
class ClassificationConfig:
    """
    Classification thresholds based on Basu et al. PLoS Genetics 2020
    
    These values distinguish Single (patient-like) from Multiple (mouse model-like)
    patterns of mitochondrial DNA structural alterations.
    """
    
    # ============================================================
    # HETEROPLASMY THRESHOLDS
    # ============================================================
    HIGH_HETEROPLASMY_THRESHOLD: float = 20.0
    """Threshold for pathogenic significance (≥30%)"""
    
    SIGNIFICANT_HETEROPLASMY_THRESHOLD: float = 1.0
    """Above technical noise, reliable detection (≥5%)"""
    
    LOW_HETEROPLASMY_THRESHOLD: float = 1.0
    """Detection limit, likely artifacts (<1%)"""
    
    # ============================================================
    # EVENT COUNT THRESHOLDS
    # ============================================================
    MAJOR_EVENT_COUNT_THRESHOLD: int = 3
    """Maximum high-het events for single pattern (≤3)"""
    
    TOTAL_EVENT_COUNT_THRESHOLD: int = 10
    """Minimum total events for multiple pattern (>20)"""
    
    MIN_EVENTS_NO_DOMINANT: int = 5
    """Events without dominant high-het for multiple pattern"""
    
    MIN_MIXED_TYPE_COUNT: int = 3
    """Minimum of each type to count as mixed (≥3)"""
    
    MIN_MEDIUM_HET_FOR_MULTIPLE: int = 10
    """Medium-het events suggesting multiple pattern"""
    
    MIN_HIGH_HET_BOTH_TYPES: int = 2
    """High-het events of both types for multiple"""
    
    MIN_SCATTERED_EVENTS: int = 4
    """Minimum events for scattered pattern"""
    
    MIN_WIDE_SIGNIFICANT: int = 5
    """Significant events with wide distribution"""
    
    MIN_LOW_DENSITY_EVENTS: int = 15
    """Events for low clustering density threshold"""
    
    # ============================================================
    # SPATIAL CLUSTERING THRESHOLDS
    # ============================================================
    CLUSTER_RADIUS: int = 600
    """Spatial grouping radius (bp)"""
    
    MIN_CLUSTER_SIZE: int = 2
    """Minimum events per cluster"""
    
    TIGHT_CLUSTER_THRESHOLD: int = 100
    """Tight clustering - same molecular event (≤100 bp)"""
    
    LOOSE_CLUSTER_THRESHOLD: int = 1000
    """Loose clustering threshold (≤1000 bp)"""
    
    SCATTERED_THRESHOLD: int = 3000
    """Scattered pattern - independent events (≥3000 bp)"""
    
    # ============================================================
    # GROUP DOMINANCE THRESHOLDS
    # ============================================================
    DOMINANT_GROUP_THRESHOLD: float = 0.7
    """Fraction for dominant group (≥70%)"""
    
    NO_DOMINANT_THRESHOLD: float = 0.5
    """Below this = no dominant group (<50%)"""
    
    # ============================================================
    # CLUSTERING DENSITY THRESHOLDS
    # ============================================================
    HIGH_CLUSTERING_DENSITY: float = 3.0
    """High clustering density (>5 events/kb)"""
    
    LOW_CLUSTERING_DENSITY: float = 1.0
    """Low clustering density (<2 events/kb)"""
    
    MIN_GROUPS_FOR_MULTIPLE: int = 3
    """Multiple spatial groups threshold (≥3)"""
    
    MIN_HIGH_HET_GROUPS: int = 2
    """Multiple high-het groups threshold (≥2)"""