"""Plot subcommand - visualization"""

import pandas as pd
from pathlib import Path

from ..visualizer import plot_circular
from ..io import BlacklistReader, read_intermediate

def add_parser(subparsers):
    """Add plot subcommand parser"""
    parser = subparsers.add_parser(
        'plot',
        help='Create circular genome visualization'
    )
    
    # Input
    parser.add_argument('-i', '--input', required=True,
                       help='Input classified intermediate TSV file from classify step (*.classified.tsv)')
    
    # Output
    parser.add_argument('-o', '--output', required=True,
                       help='Output PNG file')
    
    # Optional
    parser.add_argument('-b', '--blacklist',
                       help='BED file with regions to exclude')
    parser.add_argument('--figsize', nargs=2, type=int, default=[16, 10],
                       help='Figure size (width height), default: 16 10')
    
    return parser


def run(args):
    """Execute plot subcommand"""
    print("=== SaltShaker: Circular Plot ===\n")
    
    # Load events
    events, genome_length = read_intermediate(args.input)
    print(f"Loaded {len(events)} events from {args.input}")
    
    # Load blacklist if provided
    blacklist_regions = None
    if args.blacklist:
        blacklist_regions = BlacklistReader.load_blacklist_regions(args.blacklist)
        print(f"Loaded {len(blacklist_regions)} blacklist regions")
    
    # Create plot
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    plot_circular(
        events,
        args.output,
        genome_length,
        blacklist_regions,
        figsize=tuple(args.figsize)
    )
    
    print(f"Plot saved to: {args.output}")