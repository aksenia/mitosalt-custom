"""Classify subcommand - pattern classification"""

import pandas as pd
from pathlib import Path

from ..config import ClassificationConfig
from ..classifier import EventClassifier
from ..io import BlacklistReader, write_summary, write_vcf, read_intermediate, write_intermediate


def add_parser(subparsers):
    """Add classify subcommand parser"""
    parser = subparsers.add_parser(
        'classify',
        help='Classify event pattern as Single or Multiple'
    )
    
    # Input/Output
    parser.add_argument('-i', '--input', required=True,
                       help='Input intermediate TSV file from call step')
    parser.add_argument('-o', '--output', required=True,
                       help='Output summary text file')
    parser.add_argument('--vcf', action='store_true',
                       help='Also output classified events in VCF format')
    parser.add_argument('-b', '--blacklist',
                       help='BED file with regions to exclude')
    
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
    
    return parser


def run(args):
    """Execute classify subcommand"""
    print("=== SaltShaker: Pattern Classification ===\n")
    
    # Load events
    events, genome_length = read_intermediate(args.input)
    print(f"Loaded {len(events)} events from {args.input}")
    
    # Load blacklist
    blacklist_regions = None
    if args.blacklist:
        blacklist_regions = BlacklistReader.load_blacklist_regions(args.blacklist)
        print(f"Loaded {len(blacklist_regions)} blacklist regions")
    
    # Create config with CLI overrides
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

    # Print classification parameters
    print("\nClassification Parameters:")
    print(f"  High heteroplasmy threshold: {config.HIGH_HET_THRESHOLD:.1f}%")
    print(f"  Noise threshold: {config.NOISE_THRESHOLD:.1f}%")
    print(f"  Spatial clustering radius: {config.CLUSTER_RADIUS} bp")
    print(f"  Min cluster size: {config.MIN_CLUSTER_SIZE} events")
    print(f"  Multiple event threshold: {config.MULTIPLE_EVENT_THRESHOLD} events")
    print(f"  Dominant group fraction: {config.DOMINANT_GROUP_FRACTION:.0%}")
    print()
    
    # Classify
    classifier = EventClassifier(genome_length, config)
    classification, reason, criteria, events_with_groups = classifier.classify(
        events, blacklist_regions=blacklist_regions
    )
    
    print(f"\nClassification: {classification}")
    print(f"Reason: {reason}")
    
    # Write outputs (rest unchanged)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    write_summary(
        events_with_groups,
        args.output,
        {},
        (classification, reason, criteria, events_with_groups),
        config=config,
        blacklist_regions=blacklist_regions
    )
    
    print(f"Summary saved to: {args.output}")
    
    # Auto-save classified intermediate
    output_path = Path(args.output)
    classified_intermediate = output_path.parent / f"{output_path.stem}.classified.intermediate.tsv"
    write_intermediate(events_with_groups, str(classified_intermediate), genome_length)
    print(f"Classified events saved to: {classified_intermediate}")
    
    # Optional VCF
    if args.vcf:
        vcf_output = str(Path(args.output).with_suffix('.vcf'))
        write_vcf(
            events_with_groups,
            vcf_output,
            genome_length,
            reference_name="chrM",
            sample_name=Path(args.output).stem
        )
        print(f"VCF saved to: {vcf_output}")