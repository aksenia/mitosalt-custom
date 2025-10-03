import pandas as pd

"""
Utility Functions

General-purpose utilities used across SaltShaker modules.
"""


def crosses_blacklist(start_pos, end_pos, blacklist_regions):
    """
    Check if genomic coordinates overlap with any blacklisted region
    
    Args:
        start_pos: Event start position
        end_pos: Event end position
        blacklist_regions: List of dicts with 'start' and 'end' keys
        
    Returns:
        Boolean indicating if coordinates cross any blacklist region
    """
    if not blacklist_regions:
        return False
    
    for region in blacklist_regions:
        bl_start, bl_end = int(region['start']), int(region['end'])
        if (bl_start <= start_pos <= bl_end) or (bl_start <= end_pos <= bl_end):
            return True
    return False