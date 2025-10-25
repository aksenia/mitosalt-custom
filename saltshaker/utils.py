"""
Utility functions

General-purpose utilities used across SaltShaker modules.
"""

from __future__ import annotations
from typing import List, Optional, Union
import pandas as pd

from .types import BlacklistRegion


def crosses_blacklist(
    start_pos: Union[int, float],
    end_pos: Union[int, float],
    blacklist_regions: Optional[List[BlacklistRegion]]
) -> bool:
    """
    Check if genomic coordinates overlap with any blacklisted region
    
    Args:
        start_pos: Event start position (bp)
        end_pos: Event end position (bp)
        blacklist_regions: List of blacklist regions with 'start' and 'end' keys,
                          or None if no blacklist is active
        
    Returns:
        True if coordinates overlap any blacklist region, False otherwise
    """
    if not blacklist_regions:
        return False
    
    for region in blacklist_regions:
        bl_start: int = int(region['start'])
        bl_end: int = int(region['end'])
        if (bl_start <= start_pos <= bl_end) or (bl_start <= end_pos <= bl_end):
            return True
    return False