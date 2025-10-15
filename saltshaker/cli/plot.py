"""Plot subcommand - visualization"""

from pathlib import Path

from ..visualizer import plot_circular
from ..io import BlacklistReader, read_intermediate


def add_parser(subparsers):
    """Add plot subcommand parser"""
    parser = subparsers.add_parser(
        'plot',
        help='Create circular genome visualization'
    )
    
    # Sample identification
    parser.add_argument('--prefix', required=True,
                       help='Sample prefix (matches classify output)')
    parser.add_argument('--input-dir', required=True,
                       help='Input directory containing .saltshaker_classify_metadata.tsv from classify')
    parser.add_argument('--output-dir',
                       help='Output directory (default: same as input-dir)')
    
    # Optional
    parser.add_argument('-b', '--blacklist',
                       help='BED file with regions to exclude')
    parser.add_argument('--figsize', nargs=2, type=int, default=[16, 10],
                       help='Figure size (width height), default: 16 10')
    
    return parser


def run(args):
    """Execute plot subcommand"""
    print("=== SaltShaker: Circular Plot ===\n")
    
    # Setup directories
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filenames
    input_file = input_dir / f"{args.prefix}.saltshaker_classify_metadata.tsv"
    plot_file = output_dir / f"{args.prefix}.saltshaker.png"
    
    print(f"Sample prefix: {args.prefix}")
    print(f"Input: {input_file}")
    print(f"Output: {plot_file}")
    
    # Check input exists
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}\n"
                              f"Did you run 'saltshaker classify --prefix {args.prefix}' first?")
    
    # Load events
    events, genome_length = read_intermediate(str(input_file))
    print(f"Loaded {len(events)} events")
    
    # Load blacklist
    blacklist_regions = None
    if args.blacklist:
        blacklist_regions = BlacklistReader.load_blacklist_regions(args.blacklist)
        print(f"Loaded {len(blacklist_regions)} blacklist regions")
    
    # Create plot
    plot_circular(
        events,
        str(plot_file),
        genome_length,
        blacklist_regions,
        figsize=tuple(args.figsize)
    )
    
    print(f"\nPlot saved: {plot_file}")