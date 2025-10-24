"""Plot subcommand - visualization"""

from __future__ import annotations
from typing import Optional, List, Tuple
from pathlib import Path
import logging
from argparse import ArgumentParser, Namespace, _SubParsersAction
import pandas as pd

from ..types import BlacklistRegion, GeneAnnotation
from ..visualizer import plot_circular
from ..io import BlacklistReader, read_intermediate, GeneAnnotationReader

logger = logging.getLogger(__name__)


def add_parser(subparsers: _SubParsersAction) -> ArgumentParser:
    """
    Add plot subcommand parser
    
    Args:
        subparsers: Subparser action from main argument parser
        
    Returns:
        Configured ArgumentParser for plot subcommand
    """
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
    parser.add_argument('--figsize', nargs=2, type=int, default=[16, 10],
                       help='Figure size (width height), default: 16 10')    
    parser.add_argument('--direction', choices=['clockwise', 'counterclockwise'],
                       default='counterclockwise',
                       help='Plot direction (default: counterclockwise, field standard)')
    parser.add_argument('--del-color', choices=['red', 'blue'], default='red',
                       help='Color scheme for deletions (default: red)')
    parser.add_argument('--dup-color', choices=['red', 'blue'], default='blue',
                       help='Color scheme for duplications (default: blue)')
    parser.add_argument('--scale', choices=['dynamic', 'fixed'], default='dynamic',
                       help='Heteroplasmy color scale: dynamic (min-max per category) or fixed (0-100%%) (default: dynamic)')
    parser.add_argument('--genes', nargs='?', const='default', metavar='BED_FILE',
                        help='Enable gene annotations. Use built-in default if no file specified, or provide custom BED file path')
    parser.add_argument('--blacklist', nargs='?', const='default', metavar='BED_FILE',
                        help='Enable blacklist regions. Use built-in default if no file specified, or provide custom BED file path')

    
    return parser


def run(args: Namespace) -> None:
    """
    Execute plot subcommand
    
    Args:
        args: Parsed command-line arguments from argparse
    """
    logger.info("=== SaltShaker: Circular Plot ===")
    
    # Setup directories
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filenames
    input_file = input_dir / f"{args.prefix}.saltshaker_classify_metadata.tsv"
    plot_file = output_dir / f"{args.prefix}.saltshaker.png"
    
    logger.info(f"Sample prefix: {args.prefix}")
    logger.info(f"Input: {input_file}")
    logger.info(f"Output: {plot_file}")
    logger.info(f"Direction: {args.direction}")
    logger.info(f"Del color: {args.del_color}, Dup color: {args.dup_color}")
    logger.info(f"Scale: {args.scale}")
    
    # Check input exists
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}\n"
                              f"Did you run 'saltshaker classify --prefix {args.prefix}' first?")
    
    # Load events
    events: pd.DataFrame
    genome_length: int
    events, genome_length = read_intermediate(str(input_file))
    logger.info(f"Loaded {len(events)} events")
    
    # Load blacklist regions
    blacklist_regions: Optional[List[BlacklistRegion]] = None
    if args.blacklist is not None:
        if args.blacklist == 'default':
            # Use built-in default blacklist
            from ..data import DEFAULT_MT_BLACKLIST
            blacklist_file: str = DEFAULT_MT_BLACKLIST
            logger.info("Using default MT blacklist regions")
        else:
            # Use user-provided file
            blacklist_file = args.blacklist
            logger.info(f"Using custom blacklist: {blacklist_file}")
        
        # Validate file exists
        if not Path(blacklist_file).exists():
            raise FileNotFoundError(f"Blacklist file not found: {blacklist_file}")
        
        blacklist_regions = BlacklistReader.load_blacklist_regions(blacklist_file)
        logger.info(f"Loaded {len(blacklist_regions)} blacklist regions")

    # Load gene annotations
    gene_annotations: Optional[List[GeneAnnotation]] = None
    if args.genes is not None:
        if args.genes == 'default':
            # Use built-in default annotations
            from ..data import DEFAULT_MT_GENES
            gene_file: str = DEFAULT_MT_GENES
            logger.info("Using default hg38 MT gene annotations")
        else:
            # Use user-provided file
            gene_file = args.genes
            logger.info(f"Using custom gene annotations: {gene_file}")
        
        # Validate file exists
        if not Path(gene_file).exists():
            raise FileNotFoundError(f"Gene annotation file not found: {gene_file}")
        
        gene_annotations = GeneAnnotationReader.load_gene_annotations(gene_file)
        logger.info(f"Loaded {len(gene_annotations)} gene annotations")

    # Create plot
    logger.info("Generating plot...")
    
    figsize: Tuple[int, int] = (args.figsize[0], args.figsize[1])
    plot_circular(
        events,
        str(plot_file),
        genome_length,
        blacklist_regions,
        figsize=figsize,
        direction=args.direction,
        del_color=args.del_color,
        dup_color=args.dup_color,
        gene_annotations=gene_annotations,
        scale=args.scale
    )
    
    
    logger.info(f"✓ Plot saved: {plot_file}")