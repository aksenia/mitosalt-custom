"""Classify subcommand - pattern classification"""

from __future__ import annotations
from typing import Optional, List, Tuple, TYPE_CHECKING
from pathlib import Path
import logging
from argparse import ArgumentParser, Namespace, _SubParsersAction
import pandas as pd

# Type checking imports
if TYPE_CHECKING:
    from ..types import BlacklistRegion, ClassificationType

# These would be actual imports in the real module
# from ..config import ClassificationConfig
# from ..classifier import EventClassifier
# from ..io import BlacklistReader, write_summary, write_vcf, write_intermediate, read_intermediate

logger = logging.getLogger(__name__)


def add_parser(subparsers: _SubParsersAction) -> ArgumentParser:
    """
    Add classify subcommand parser
    
    Args:
        subparsers: Subparser action from main argument parser
        
    Returns:
        Configured ArgumentParser for classify subcommand
    """
    parser = subparsers.add_parser(
        'classify',
        help='Classify event pattern as Single or Multiple'
    )
    
    # Sample identification
    parser.add_argument('--prefix', required=True,
                       help='Sample prefix (matches call output)')
    parser.add_argument('--input-dir', required=True,
                       help='Input directory containing .saltshaker_call_metadata.tsv from call')
    parser.add_argument('--output-dir',
                       help='Output directory (default: same as input-dir)')
    
    # Output options
    parser.add_argument('--vcf', action='store_true',
                       help='Also output classified events in VCF format')
    
    # Optional blacklist
    parser.add_argument('--blacklist', nargs='?', const='default', metavar='BED_FILE',
                       help='Enable blacklist regions. Use built-in default if no file specified, or provide custom BED file path')
    
    # Classification parameters (optional, use config defaults)
    parser.add_argument('--high-het', type=float,
                       help='High heteroplasmy threshold %% (default: 20)')
    parser.add_argument('--noise', type=float,
                       help='Noise threshold %% (default: 1)')
    parser.add_argument('--radius', type=int,
                       help='Spatial clustering radius bp (default: 600)')
    parser.add_argument('--multiple-threshold', type=int,
                       help='Event count for Multiple pattern (default: 10)')
    parser.add_argument('--dominant-fraction', type=float,
                       help='Fraction for dominant group (default: 0.70)')
    
    return parser  # type: ignore[no-any-return]


def run(args: Namespace) -> None:
    """
    Execute classify subcommand
    
    Args:
        args: Parsed command-line arguments from argparse
    """
    logger.info("=== SaltShaker: Pattern Classification ===")
    
    # Setup directories
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filenames
    input_file = input_dir / f"{args.prefix}.saltshaker_call_metadata.tsv"
    summary_file = output_dir / f"{args.prefix}.saltshaker_classify.txt"
    classified_file = output_dir / f"{args.prefix}.saltshaker_classify_metadata.tsv"
    vcf_file = output_dir / f"{args.prefix}.saltshaker.vcf"
    
    logger.info(f"Sample prefix: {args.prefix}")
    logger.info(f"Input: {input_file}")
    logger.info(f"Output directory: {output_dir}")
    
    # Check input exists
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}\n"
                              f"Did you run 'saltshaker call --prefix {args.prefix}' first?")
    
    # Load events
    from ..io import read_intermediate  # type: ignore
    events: pd.DataFrame
    genome_length: int
    events, genome_length = read_intermediate(str(input_file))
    logger.info(f"Loaded {len(events)} events")
    logger.info(f"Genome length: {genome_length}")
    
    # Load blacklist regions
    blacklist_regions: Optional[List[BlacklistRegion]] = None
    if args.blacklist is not None:
        if args.blacklist == 'default':
            # Use built-in default blacklist
            from ..data import DEFAULT_MT_BLACKLIST  # type: ignore
            blacklist_file: str = DEFAULT_MT_BLACKLIST
            logger.info("Using default MT blacklist regions")
        else:
            # Use user-provided file
            blacklist_file = args.blacklist
            logger.info(f"Using custom blacklist: {blacklist_file}")
        
        # Validate file exists
        if not Path(blacklist_file).exists():
            raise FileNotFoundError(f"Blacklist file not found: {blacklist_file}")
        
        from ..io import BlacklistReader  # type: ignore
        blacklist_regions = BlacklistReader.load_blacklist_regions(blacklist_file)
        logger.info(f"Loaded {len(blacklist_regions)} blacklist regions")
    
    # Create config with CLI overrides
    from ..config import ClassificationConfig  # type: ignore
    config = ClassificationConfig()
    if args.high_het is not None:
        config.HIGH_HET_THRESHOLD = args.high_het
    if args.noise is not None:
        config.NOISE_THRESHOLD = args.noise
    if args.radius is not None:
        config.CLUSTER_RADIUS = args.radius
    if args.multiple_threshold is not None:
        config.MULTIPLE_EVENT_THRESHOLD = args.multiple_threshold
    if args.dominant_fraction is not None:
        config.DOMINANT_GROUP_FRACTION = args.dominant_fraction
    
    # Print parameters
    logger.info("Classification Parameters:")
    logger.info(f"  High heteroplasmy threshold: {config.HIGH_HET_THRESHOLD:.1f}%")
    logger.info(f"  Noise threshold: {config.NOISE_THRESHOLD:.1f}%")
    logger.info(f"  Spatial clustering radius: {config.CLUSTER_RADIUS} bp")
    logger.info(f"  Min cluster size: {config.MIN_CLUSTER_SIZE} events")
    logger.info(f"  Multiple event threshold: {config.MULTIPLE_EVENT_THRESHOLD} events")
    logger.info(f"  Dominant group fraction: {config.DOMINANT_GROUP_FRACTION:.0%}")
    
    # Classify
    from ..classifier import EventClassifier  # type: ignore
    classifier = EventClassifier(genome_length, config)
    
    classification: ClassificationType
    reason: str
    criteria: dict
    events_with_groups: pd.DataFrame
    classification, reason, criteria, events_with_groups = classifier.classify(
        events, blacklist_regions=blacklist_regions
    )
    
    logger.info(f"Classification: {classification}")
    logger.info(f"Reason: {reason}")
    
    # Write summary
    from ..io import write_summary, write_vcf, write_intermediate  # type: ignore
    write_summary(
        events_with_groups,
        str(summary_file),
        {},
        (classification, reason, criteria, events_with_groups),
        config=config,
        blacklist_regions=blacklist_regions
    )
    logger.info(f"Summary: {summary_file}")
    
    # Always save classified intermediate for plotting
    write_intermediate(events_with_groups, str(classified_file), genome_length)
    logger.info(f"Classified events: {classified_file} (use for plot)")
    
    # Optional VCF
    if args.vcf:
        write_vcf(
            events_with_groups,
            str(vcf_file),
            genome_length,
            reference_name="chrM",
            sample_name=args.prefix
        )
        logger.info(f"VCF: {vcf_file}")