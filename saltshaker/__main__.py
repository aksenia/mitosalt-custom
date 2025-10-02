"""
SaltShaker CLI

Command-line interface for SaltShaker event calling, classification, and visualization.
"""

import argparse
import numpy as np
from pathlib import Path

from .config import ClassificationConfig
from .event_caller import EventCaller
from .classifier import EventClassifier
from .io import BlacklistReader, write_tsv, write_vcf, write_summary
from .visualizer import plot_circular


def main():
    parser = argparse.ArgumentParser(
        description='SaltShaker: Pattern classification and visualization for MitoSAlt',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('genome_length', type=int, help='Mitochondrial genome length')
    parser.add_argument('ori_h_start', type=int, help='Heavy strand origin start')
    parser.add_argument('ori_h_end', type=int, help='Heavy strand origin end') 
    parser.add_argument('ori_l_start', type=int, help='Light strand origin start')
    parser.add_argument('ori_l_end', type=int, help='Light strand origin end')
    parser.add_argument('size_limit', type=int, help='Size limit for events')
    parser.add_argument('cluster_file', help='Cluster file path')
    parser.add_argument('breakpoint_file', help='Breakpoint file path')
    parser.add_argument('output_name', help='Output file prefix')
    parser.add_argument('heteroplasmy_limit', type=float, help='Heteroplasmy threshold')
    parser.add_argument('genome_fasta', help='Genome FASTA file')
    parser.add_argument('flank_size', type=int, help='Flanking sequence size')
    
    # Optional arguments
    parser.add_argument('--blacklist', help='BED file with regions to exclude', default=None)
    parser.add_argument('--output-dir', help='Output directory (default: current directory)', default='.')
    parser.add_argument('--output-vcf', action='store_true', help='Also output events in VCF format')
    
    args = parser.parse_args()
    
    # Create config
    config = ClassificationConfig()
    
    # Initialize EventCaller
    event_caller = EventCaller(
        genome_length=args.genome_length,
        ori_h=(args.ori_h_start, args.ori_h_end),
        ori_l=(args.ori_l_start, args.ori_l_end),
        heteroplasmy_limit=args.heteroplasmy_limit,
        flank_size=args.flank_size
    )
    
    # Load blacklist regions
    blacklist_regions = None
    if args.blacklist:
        try:
            blacklist_regions = BlacklistReader.load_blacklist_regions(args.blacklist)
            print(f"Loaded {len(blacklist_regions)} blacklist regions")
        except Exception as e:
            print(f"Warning: Could not load blacklist file: {e}")
            blacklist_regions = []
    
    # Call events (del/dup classification)
    events = event_caller.call_events(args.cluster_file, args.breakpoint_file)
    stats = {}
    
    if len(events) > 0:
        # Calculate final coordinates (needed for flanking sequences)
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
        events = event_caller.add_flanking_sequences(events, args.genome_fasta)
        
        # Classify pattern (Single/Multiple)
        classifier = EventClassifier(args.genome_length, config)
        classification, reason, criteria, events_with_groups = classifier.classify(
            events, blacklist_regions=blacklist_regions
        )
        
        # Create output directories
        output_dir = Path(args.output_dir)
        plot_dir = output_dir / "plot"
        indel_dir = output_dir / "indel"
        
        plot_dir.mkdir(parents=True, exist_ok=True)
        indel_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Processing {len(events)} events")
        
        # Generate circular plot
        plot_file = plot_dir / f"{args.output_name}.png"
        plot_circular(
            events_with_groups, 
            str(plot_file), 
            args.genome_length, 
            blacklist_regions, 
            figsize=(10, 10)
        )
        
        # Write TSV results
        results_file = indel_dir / f"{args.output_name}.grouped.tsv"
        write_tsv(events_with_groups, str(results_file), args.genome_length, blacklist_regions)
        
        # Write summary
        summary_file = indel_dir / f"{args.output_name}_summary.txt"
        write_summary(
            events_with_groups, 
            str(summary_file), 
            stats, 
            (classification, reason, criteria, events_with_groups),
            config=config,
            blacklist_regions=blacklist_regions
        )
        
        # Write VCF if requested
        if args.output_vcf:
            vcf_file = indel_dir / f"{args.output_name}.vcf"
            write_vcf(
                events_with_groups,
                (classification, reason, criteria, events_with_groups),
                str(vcf_file),
                reference_name="MT",
                sample_name=args.output_name,
                genome_length=args.genome_length
            )
    else:
        print("No events to process")


if __name__ == "__main__":
    main()