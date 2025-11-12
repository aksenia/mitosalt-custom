"""
Layout Module for SaltShaker
Smart dynamic layout engine for circular mitochondrial plots

Public API:
    - LayoutEngine: Main layout calculation engine
    - LayoutResult: Complete layout solution
    - SectorBudget: Sector budget allocation
    - GroupBandLayout: Multi-event group layout
    - SingleEventLayout: Single event layout
"""

from .engine import LayoutEngine
from .types import (
    LayoutResult,
    SectorBudget,
    GroupBandLayout,
    SingleEventLayout,
    BlacklistRegion,
    GeneAnnotation,
)

__all__ = [
    'LayoutEngine',
    'LayoutResult',
    'SectorBudget',
    'GroupBandLayout',
    'SingleEventLayout',
    'BlacklistRegion',
    'GeneAnnotation',
]
