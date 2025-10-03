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
    
    # Input
    parser.add_argument('-i', '--input', required=True,
                       help='Input TSV file from call step')
    
    # Output
    parser.add_argument('-o', '--output', required=True,
                       help='Output summary text file')
    parser.add_argument('--vcf', action='store_true',
                       help='Also output classified events in VCF format')
    
    # Optional filters
    parser.add_argument('-b', '--blacklist',
                       help='BED file with regions to exclude')
    
    return parser


def run(args):
    """Execute classify subcommand"""
    print("=== SaltShaker: Pattern Classification ===\n")
    
    # Load events with metadata
    events, genome_length = read_intermediate(args.input)
    print(f"Loaded {len(events)} events from {args.input}")
    print(f"Genome length: {genome_length}")

    
    # Load blacklist if provided
    blacklist_regions = None
    if args.blacklist:
        blacklist_regions = BlacklistReader.load_blacklist_regions(args.blacklist)
        print(f"Loaded {len(blacklist_regions)} blacklist regions")
    
    # Classify pattern
    config = ClassificationConfig()
    classifier = EventClassifier(genome_length, config)
    classification, reason, criteria, events_with_groups = classifier.classify(
        events, blacklist_regions=blacklist_regions
    )
    
    print(f"\nClassification: {classification}")
    print(f"Reason: {reason}")
    
    # Write summary
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_summary(
        events_with_groups,
        args.output,
        {},  # stats
        (classification, reason, criteria, events_with_groups),
        config=config,
        blacklist_regions=blacklist_regions
    )
    
    print(f"Summary saved to: {args.output}")
    
    # Always save classified intermediate TSV for downstream processing (plot needs this)
    output_path = Path(args.output)
    classified_intermediate = output_path.parent / f"{output_path.stem}.classified.tsv"

    write_intermediate(events_with_groups, str(classified_intermediate), genome_length)
    print(f"Classified events saved to: {classified_intermediate} (use for plot)")
        
    # Optionally save VCF with groups and classification
    if args.vcf:
        vcf_output = str(Path(args.output).with_suffix('.vcf'))
        
        write_vcf(
            events_with_groups,
            vcf_output,
            genome_length,
            reference_name="chrM",
            sample_name=Path(args.output).stem
        )
        print(f"VCF with classification saved to: {vcf_output}")