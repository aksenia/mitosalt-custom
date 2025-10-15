"""
SaltShaker Configuration
Naive classification parameters for mitochondrial event patterns
"""
from dataclasses import dataclass


@dataclass
class ClassificationConfig:
    """
    Classification thresholds for Single vs Multiple patterns
    Based on Basu et al. PLoS Genetics 2020
    """
    
    # ============================================================
    # HETEROPLASMY THRESHOLDS
    # ============================================================
    HIGH_HET_THRESHOLD: float = 20.0
    """High heteroplasmy threshold - pathogenic significance (≥20%)"""
    
    NOISE_THRESHOLD: float = 0.7
    """Noise threshold - below this likely artifacts (<1%)"""
    
    # ============================================================
    # SPATIAL CLUSTERING
    # ============================================================
    CLUSTER_RADIUS: int = 600
    """Spatial grouping radius in base pairs"""
    
    MIN_CLUSTER_SIZE: int = 2
    """Minimum events to form a spatial cluster"""
    
    # ============================================================
    # PATTERN CLASSIFICATION
    # ============================================================
    MULTIPLE_EVENT_THRESHOLD: int = 5
    """Event count threshold for Multiple pattern (>10 = Multiple)"""
    
    DOMINANT_GROUP_FRACTION: float = 0.70
    """Fraction of events in dominant group for Single pattern (≥70%)"""