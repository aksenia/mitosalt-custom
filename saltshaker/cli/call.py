"""Call subcommand - event calling (del/dup classification)"""

from __future__ import annotations
from typing import Optional, List, TYPE_CHECKING
import numpy as np
from pathlib import Path
import logging
from argparse import ArgumentParser, Namespace, _SubParsersAction

# Runtime imports
from ..event_caller import EventCaller
from ..io import BlacklistReader, write_tsv, write_intermediate
from ..data import DEFAULT_MT_BLACKLIST

# Type checking imports (not used at runtime)
if TYPE_CHECKING:
    from ..types import BlacklistRegion

logger = logging.getLogger(__name__)


def add_parser(subparsers: _SubParsersAction) -> ArgumentParser:
    """
    Add call subcommand parser
    
    Args:
        subparsers: Subparser action from main argument parser
        
    Returns:
        Configured ArgumentParser for call subcommand
    """
    parser = subparsers.add_parser(
        'call',
        help='Call events as deletions or duplications'
    )
    
    # Sample identification
    parser.add_argument('--prefix', required=True,
                       help='Sample prefix for output files')
    parser.add_argument('--output-dir', default='.',
                       help='Output directory (default: current directory)')
    
    # Input files
    parser.add_argument('-c', '--cluster', required=True,
                       help='Cluster file path')
    parser.add_argument('-p', '--breakpoint', required=True,
                       help='Breakpoint file path')
    parser.add_argument('-r', '--reference', required=True,
                       help='Reference genome FASTA file')
    
    # Genome parameters
    parser.add_argument('-g', '--genome-length', type=int, required=True,
                       help='Mitochondrial genome length')
    parser.add_argument('--ori-h-start', type=int, required=True,
                       help='Heavy strand origin start')
    parser.add_argument('--ori-h-end', type=int, required=True,
                       help='Heavy strand origin end')
    parser.add_argument('--ori-l-start', type=int, required=True,
                       help='Light strand origin start')
    parser.add_argument('--ori-l-end', type=int, required=True,
                       help='Light strand origin end')
    
    # Analysis parameters
    parser.add_argument('-H', '--het-limit', type=float, default=0.01,
                       help='Heteroplasmy threshold (default: 0.01)')
    parser.add_argument('-f', '--flank-size', type=int, default=15,
                       help='Flanking sequence size in bp (default: 15)')
    
    # Optional blacklist
    parser.add_argument('--blacklist', nargs='?', const='default', metavar='BED_FILE',
                       help='Enable blacklist regions. Use built-in default if no file specified, or provide custom BED file path')
    
    return parser  # type: ignore[no-any-return]


def run(args: Namespace) -> None:
    """
    Execute call subcommand
    
    Args:
        args: Parsed command-line arguments from argparse
    """
    logger.info("=== SaltShaker: Event Calling ===")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate output filenames
    display_tsv = output_dir / f"{args.prefix}.saltshaker_call.tsv"
    intermediate_tsv = output_dir / f"{args.prefix}.saltshaker_call_metadata.tsv"
    
    logger.info(f"Sample prefix: {args.prefix}")
    logger.info(f"Output directory: {output_dir}")
    
    # Initialize EventCaller
    event_caller = EventCaller(
        genome_length=args.genome_length,
        ori_h=(args.ori_h_start, args.ori_h_end),
        ori_l=(args.ori_l_start, args.ori_l_end),
        heteroplasmy_limit=args.het_limit,
        flank_size=args.flank_size
    )
    
    # Load blacklist regions
    blacklist_regions: Optional[List[BlacklistRegion]] = None
    if args.blacklist is not None:
        if args.blacklist == 'default':
            # Use built-in default blacklist
            blacklist_file: str = DEFAULT_MT_BLACKLIST
            logger.info("Using default MT blacklist regions")
        else:
            # Use user-provided file
            blacklist_file = args.blacklist
            logger.info(f"Using custom blacklist: {blacklist_file}")
        
        # Validate file exists
        if not Path(blacklist_file).exists():
            raise FileNotFoundError(f"Blacklist file not found: {blacklist_file}")
        
        try:
            blacklist_regions = BlacklistReader.load_blacklist_regions(blacklist_file)
            logger.info(f"Loaded {len(blacklist_regions)} blacklist regions")
        except Exception as e:
            logger.warning(f"Could not load blacklist file: {e}")
            blacklist_regions = []
    
    # Call events
    events = event_caller.call_events(args.cluster, args.breakpoint)
    
    if events.empty:
        logger.warning("No events called")
        return
    
    # Calculate final coordinates
    events['perc'] = events['perc'].round(4)
    events['del_start_median'] = events['del_start_median'] + 1
    events['final_event_size'] = np.where(
        events['final_event'] == 'del',
        events['delsize'],
        args.genome_length - events['delsize']
    )
    events['final_end'] = np.where(
        events['final_event'] == 'del',
        events['del_end_median'],
        events['del_start_median'] - 1
    )
    events['final_start'] = np.where(
        events['final_event'] == 'del',
        events['del_start_median'],
        events['del_end_median'] + 1
    )
    
    # Add flanking sequences
    events = event_caller.add_flanking_sequences(events, args.reference)
    
    # Write intermediate format (all columns) for downstream processing
    write_intermediate(events, str(intermediate_tsv), args.genome_length)
    
    # Write display TSV (formatted, human-readable)
    write_tsv(events, str(display_tsv), args.genome_length, blacklist_regions)
    
    logger.info(f"Called {len(events)} events")
    logger.info(f"Display TSV: {display_tsv}")
    logger.info(f"Intermediate TSV: {intermediate_tsv} (use for classify)")