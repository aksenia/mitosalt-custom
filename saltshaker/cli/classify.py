"""Classify subcommand - pattern classification"""

from pathlib import Path

from ..config import ClassificationConfig
from ..classifier import EventClassifier
from ..io import BlacklistReader, write_summary, write_vcf, write_intermediate, read_intermediate


def add_parser(subparsers):
    """Add classify subcommand parser"""
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
    
    # Setup directories
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filenames
    input_file = input_dir / f"{args.prefix}.saltshaker_call.tsv"
    summary_file = output_dir / f"{args.prefix}.saltshaker_classify.txt"
    classified_file = output_dir / f"{args.prefix}.saltshaker_classify_metadata.tsv"
    vcf_file = output_dir / f"{args.prefix}.vcf"
    
    print(f"Sample prefix: {args.prefix}")
    print(f"Input: {input_file}")
    print(f"Output directory: {output_dir}")
    
    # Check input exists
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}\n"
                              f"Did you run 'saltshaker call --prefix {args.prefix}' first?")
    
    # Load events
    events, genome_length = read_intermediate(str(input_file))
    print(f"Loaded {len(events)} events")
    print(f"Genome length: {genome_length}")
    
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
    
    # Print parameters
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
    
    print(f"Classification: {classification}")
    print(f"Reason: {reason}\n")
    
    # Write summary
    write_summary(
        events_with_groups,
        str(summary_file),
        {},
        (classification, reason, criteria, events_with_groups),
        config=config,
        blacklist_regions=blacklist_regions
    )
    print(f"Summary: {summary_file}")
    
    # Always save classified intermediate for plotting
    write_intermediate(events_with_groups, str(classified_file), genome_length)
    print(f"Classified events: {classified_file} (use for plot)")
    
    # Optional VCF
    if args.vcf:
        write_vcf(
            events_with_groups,
            str(vcf_file),
            genome_length,
            reference_name="chrM",
            sample_name=args.prefix
        )
        print(f"VCF: {vcf_file}")