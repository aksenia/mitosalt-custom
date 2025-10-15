"""Call subcommand - event calling (del/dup classification)"""

import numpy as np
from pathlib import Path

from ..event_caller import EventCaller
from ..io import BlacklistReader, write_tsv, write_intermediate


def add_parser(subparsers):
    """Add call subcommand parser"""
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
    
    # Optional filters
    parser.add_argument('-b', '--blacklist',
                       help='BED file with regions to exclude')
    
    return parser


def run(args):
    """Execute call subcommand"""
    print("=== SaltShaker: Event Calling ===\n")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate output filenames
    display_tsv = output_dir / f"{args.prefix}.saltshaker_call.tsv"
    intermediate_tsv = output_dir / f"{args.prefix}.saltshaker_call_metadata.tsv"
    
    print(f"Sample prefix: {args.prefix}")
    print(f"Output directory: {output_dir}")
    
    # Initialize EventCaller
    event_caller = EventCaller(
        genome_length=args.genome_length,
        ori_h=(args.ori_h_start, args.ori_h_end),
        ori_l=(args.ori_l_start, args.ori_l_end),
        heteroplasmy_limit=args.het_limit,
        flank_size=args.flank_size
    )
    
    # Load blacklist if provided
    blacklist_regions = None
    if args.blacklist:
        try:
            blacklist_regions = BlacklistReader.load_blacklist_regions(args.blacklist)
            print(f"Loaded {len(blacklist_regions)} blacklist regions")
        except Exception as e:
            print(f"Warning: Could not load blacklist file: {e}")
            blacklist_regions = []
    
    # Call events
    events = event_caller.call_events(args.cluster, args.breakpoint)
    
    if len(events) == 0:
        print("No events called")
        return
    
    # Calculate final coordinates
    events['perc'] = events['perc'].round(4)
    events['del.start.median'] = events['del.start.median'] + 1
    events['final.event.size'] = np.where(
        events['final.event'] == 'del',
        events['delsize'],
        args.genome_length - events['delsize']
    )
    events['final.end'] = np.where(
        events['final.event'] == 'del',
        events['del.end.median'],
        events['del.start.median'] - 1
    )
    events['final.start'] = np.where(
        events['final.event'] == 'del',
        events['del.start.median'],
        events['del.end.median'] + 1
    )
    
    # Add flanking sequences
    events = event_caller.add_flanking_sequences(events, args.reference)
    
    # Write intermediate format (all columns) for downstream processing
    write_intermediate(events, str(intermediate_tsv), args.genome_length)
    
    # Write display TSV (formatted, human-readable)
    write_tsv(events, str(display_tsv), args.genome_length, blacklist_regions)
    
    print(f"\nCalled {len(events)} events")
    print(f"Display TSV: {display_tsv}")
    print(f"Intermediate TSV: {intermediate_tsv} (use for classify)")