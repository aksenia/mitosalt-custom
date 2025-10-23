"""
I/O Writers

Handles writing of analysis results in various formats.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from ..utils import crosses_blacklist

logger = logging.getLogger(__name__)


class TSVWriter:
    """Writes event results in TSV format"""
    
    def __init__(self, genome_length):
        """
        Initialize TSV writer
        
        Args:
            genome_length: Mitochondrial genome length
        """
        self.genome_length = genome_length
    
    def write(self, events, output_file, blacklist_regions=None):
        """
        Write events to TSV file with coordinate swapping for output
        
        This handles the R script's coordinate swapping logic for display.
        NOTE: Does NOT include 'group' column - that's added by classify step.
        
        Args:
            events: DataFrame with events (from call step, no groups)
            output_file: Path to output TSV file
            blacklist_regions: List of blacklist regions for flagging
        """
        if len(events) == 0:
            logger.warning("No events to save")
            return
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
#        events.to_pickle(output_file + '.pkl')  # Save intermediate pickle for debugging
        
        res = events.copy()
        logger.debug(f"Starting save_results with {len(res)} events")
        
        # Coordinate swapping for OUTPUT only (R script logic)
        for i in res.index:
            if res.loc[i, 'dloop'] == 'yes':
                # Swap coordinates for display
                del_start_median = res.loc[i, 'del.start.median']
                del_end_median = res.loc[i, 'del.end.median']
                del_start_range = res.loc[i, 'del.start.range']
                del_end_range = res.loc[i, 'del.end.range']
                del_start = res.loc[i, 'del.start']
                del_end = res.loc[i, 'del.end']
                lfstart = res.loc[i, 'lfstart']
                lfend = res.loc[i, 'lfend']
                
                res.loc[i, 'del.start.median'] = del_end_median
                res.loc[i, 'del.end.median'] = del_start_median
                res.loc[i, 'del.start.range'] = del_end_range
                res.loc[i, 'del.end.range'] = del_start_range
                res.loc[i, 'del.start'] = del_end
                res.loc[i, 'del.end'] = del_start
                res.loc[i, 'lfstart'] = lfend
                res.loc[i, 'lfend'] = lfstart
        
        # Format exactly like R script
        res['perc'] = res['perc'].round(4)
        
        # Calculate final coordinates exactly like R script
        res['del.start.median'] = res['del.start.median'] + 1
        res['final.event.size'] = np.where(
            res['final.event'] == 'del',
            res['delsize'],
            self.genome_length - res['delsize']
        )
        res['final.end'] = np.where(
            res['final.event'] == 'del',
            res['del.end.median'],
            res['del.start.median'] - 1
        )
        res['final.start'] = np.where(
            res['final.event'] == 'del',
            res['del.start.median'],
            res['del.end.median'] + 1
        )
        
        # Handle wraparound
        res['del.start.median'] = np.where(
            res['del.start.median'] == self.genome_length + 1,
            1,
            res['del.start.median']
        )
        res['final.start'] = np.where(
            res['final.start'] == self.genome_length + 1,
            1,
            res['final.start']
        )

        # Blacklist crossing flag using final coordinates
        res['blacklist_crossing'] = [
            'yes' if crosses_blacklist(row['final.start'], row['final.end'], blacklist_regions) else 'no'
            for _, row in res.iterrows()
        ]
        
        # Create final output - NO GROUP COLUMN (added by classify)
        res_final = pd.DataFrame({
            'sample': res['sample'],
            'cluster.id': res['cluster'],
            'alt.reads': res['nread'].astype(int),
            'ref.reads': res['tread'].astype(int),
            'heteroplasmy': res['perc'],
            'del.start.range': res['del.start.range'],
            'del.end.range': res['del.end.range'],
            'del.size': res['delsize'].astype(int),
            'final.event': res['final.event'],
            'final.start': res['final.start'].astype(int),
            'final.end': res['final.end'].astype(int),
            'final.size': res['final.event.size'].astype(int),
            'blacklist_crossing': res['blacklist_crossing'],
            'seq1': res['seq1'],
            'seq2': res['seq2'],
            'seq': res['seq']
        })
        
        # Save results
        res_final.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Results saved to {output_file}")
        logger.info(f"Events: {len(res_final)}")


def write_tsv(events, output_file, genome_length, blacklist_regions=None):
    """
    Convenience function to write TSV
    
    Args:
        events: DataFrame with events
        output_file: Output TSV file path
        genome_length: Mitochondrial genome length
        blacklist_regions: List of blacklist regions
    """
    writer = TSVWriter(genome_length)
    writer.write(events, output_file, blacklist_regions)


class IntermediateWriter:
    """Writes events in full internal format for downstream processing"""
    
    @staticmethod
    def write(events, output_file, genome_length):
        """
        Write events with all columns preserved for downstream processing
        
        Args:
            events: DataFrame with all internal columns
            output_file: Path to output file
            genome_length: Genome length for metadata
        """
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Write metadata + full dataframe
        with open(output_file, 'w') as f:
            f.write(f"# genome_length={genome_length}\n")
            events.to_csv(f, sep='\t', index=False)
        
        logger.info(f"Intermediate results saved to {output_file}")
    
def write_intermediate(events, output_file, genome_length):
    """Convenience function"""
    IntermediateWriter.write(events, output_file, genome_length)


class SummaryWriter:
    """Writes analysis summary in human-readable text format"""
    
    def __init__(self, config=None):
        """
        Initialize summary writer
        
        Args:
            config: ClassificationConfig instance for threshold reporting
        """
        from ..config import ClassificationConfig
        self.config = config or ClassificationConfig()
    
    def write(self, events, output_file, analysis_stats, classification_result, blacklist_regions=None):
        """
        Write analysis summary with classification metrics
        
        Args:
            events: DataFrame with events
            output_file: Path to output summary file
            analysis_stats: Dict with workflow statistics
            classification_result: Tuple of (classification, reason, criteria, events_with_groups)
            blacklist_regions: List of blacklist regions
        """
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Load config for threshold reporting
        cfg = self.config
        
        # Use passed classification results
        if classification_result:
            classification, reason, criteria, events_with_groups = classification_result
        else:
            # Fallback if not provided
            if len(events) > 0:
                classifier = EventClassifier(self.genome_length, self.config)
                classification, reason, criteria, events_with_groups = classifier.classify(events, blacklist_regions)
            else:
                classification, reason = "No events", "No events detected"
                criteria = {}
                events_with_groups = pd.DataFrame()
        
        # Count events by type
        del_count = (events['final.event'] == 'del').sum() if len(events) > 0 else 0
        dup_count = (events['final.event'] == 'dup').sum() if len(events) > 0 else 0
        dloop_count = (events['dloop'] == 'yes').sum() if len(events) > 0 else 0
        
        # Calculate comprehensive statistics
        if len(events) > 0:
            size_stats = {
                'min': events['delsize'].min(),
                'max': events['delsize'].max(),
                'mean': events['delsize'].mean(),
                'median': events['delsize'].median()
            }
            het_stats = {
                'min': events['perc'].min(),
                'max': events['perc'].max(),
                'mean': events['perc'].mean(),
                'median': events['perc'].median()
            }
        else:
            size_stats = {'min': 0, 'max': 0, 'mean': 0, 'median': 0}
            het_stats = {'min': 0, 'max': 0, 'mean': 0, 'median': 0}
        
        # Write summary file with literature-based analysis
        with open(output_file, 'w') as f:
            f.write("MitoSAlt Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            
            # Analysis workflow stats
            f.write("Analysis Workflow:\n")
            f.write("-" * 20 + "\n")
            for key, value in analysis_stats.items():
                f.write(f"{key}: {value}\n")
            f.write("\n")
            
            # Event pattern classification (PRIMARY SECTION - most important)
            f.write("Event pattern classification (blacklist-filtered):\n")
            f.write("-" * 50 + "\n")
            f.write(f"Type: {classification}\n")
            f.write(f"Biological basis: {reason}\n")
            if 'subtype' in criteria:
                f.write(f"Subtype: {criteria['subtype']}\n")
            if criteria.get('blacklist_filtered_count', 0) > 0:
                f.write(f"Blacklist-crossing events excluded: {criteria['blacklist_filtered_count']}\n")
                f.write(f"Events used for classification: {criteria.get('total_events', 0)} (of {criteria.get('total_raw_events', 0)} total)\n")
            f.write("\n")
            
            # Spatial groups analysis
            f.write("Spatial Groups Analysis:\n")
            f.write("-" * 25 + "\n")
            if 'group_analysis' in criteria and criteria['group_analysis']:
                # Group by event type for better reporting
                del_groups = [g for g in criteria['group_analysis'] if g.get('event_type') == 'del']
                dup_groups = [g for g in criteria['group_analysis'] if g.get('event_type') == 'dup']
                
                if del_groups:
                    f.write("Deletion groups:\n")
                    for group_info in del_groups:
                        rep = group_info['representative']
                        actual_size = rep['end'] - rep['start']
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1f}% "
                            f"(at {rep['start']:.0f}-{rep['end']:.0f}bp, size {actual_size:.0f}bp)\n")
                
                if dup_groups:
                    f.write("Duplication groups:\n")
                    for group_info in dup_groups:
                        rep = group_info['representative']
                        actual_size = rep['end'] - rep['start']
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1f}% "
                            f"(at {rep['start']:.0f}-{rep['end']:.0f}bp, size {actual_size:.0f}bp)\n")
                
                if not del_groups and not dup_groups:
                    f.write("Groups detected but event types not properly classified\n")
                    # Fallback - show all groups regardless of type
                    for group_info in criteria['group_analysis']:
                        rep = group_info['representative']
                        event_type = group_info.get('event_type', rep.get('event_type', 'unknown'))
                        f.write(f"  {group_info['group_id']}: {group_info['event_count']} {event_type} events, "
                            f"max heteroplasmy {rep['heteroplasmy']:.1%}\n")
            else:
                f.write("No spatial groups identified\n")
            f.write("\n")
            
            # Heteroplasmy-based analysis (KEY BIOLOGICAL METRIC)
            f.write("Heteroplasmy distribution (literature-based thresholds):\n")
            f.write("-" * 55 + "\n")
            if len(events) > 0:
                f.write(f"High heteroplasmy events (≥{cfg.HIGH_HET_THRESHOLD:.0f}%): {criteria.get('high_het_count', 0)}\n")
                f.write(f"Low heteroplasmy events (<{cfg.NOISE_THRESHOLD:.0f}%): {criteria.get('low_het_count', 0)}\n")
                f.write(f"Maximum heteroplasmy: {criteria.get('max_heteroplasmy', 0):.3f}%\n")
                f.write(f"Median heteroplasmy: {criteria.get('median_heteroplasmy', 0):.3f}%\n")
            else:
                f.write("No events detected\n")
            f.write("\n")
            
            # Event counts summary
            f.write("Event summary:\n")
            f.write("-" * 15 + "\n")
            f.write(f"Total events: {len(events)}\n")
            f.write(f"Deletions: {del_count}\n")
            f.write(f"Duplications: {dup_count}\n")
            f.write(f"Events crossing origin (dloop=yes): {dloop_count}\n\n")
            
            # Basic genomic metrics (for reference)
            f.write("Basic metrics:\n")
            f.write("-" * 15 + "\n")
            if len(events) > 0:
                f.write(f"Position range spanned: {criteria.get('position_range', 0):.0f} bp\n")
                f.write(f"Size coefficient of variation: {criteria.get('size_coefficient_variation', 0):.2f}\n")
            else:
                f.write("No events to analyze\n")
            f.write("\n")
            
            # Size statistics
            f.write("Size statistics (bp):\n")
            f.write("-" * 20 + "\n")
            f.write(f"Min size: {size_stats['min']:.0f}\n")
            f.write(f"Max size: {size_stats['max']:.0f}\n")
            f.write(f"Mean size: {size_stats['mean']:.1f}\n")
            f.write(f"Median size: {size_stats['median']:.1f}\n")
            if len(events) > 0:
                f.write(f"Size coefficient of variation: {criteria.get('size_coefficient_variation', 0):.2f}\n")
            f.write("\n")
            
            # Heteroplasmy statistics
            f.write("Heteroplasmy statistics:\n")
            f.write("-" * 25 + "\n")
            f.write(f"Min heteroplasmy: {het_stats['min']:.4f}%\n")
            f.write(f"Max heteroplasmy: {het_stats['max']:.4f}%\n")
            f.write(f"Mean heteroplasmy: {het_stats['mean']:.4f}%\n")
            f.write(f"Median heteroplasmy: {het_stats['median']:.4f}%\n\n")

            # Classification details (DETAILED TECHNICAL SECTION)
            if criteria and 'classification_scores' in criteria:
                f.write("Classification Algorithm Details:\n")
                f.write("-" * 35 + "\n")
                scores = criteria['classification_scores']
                f.write(f"Single pattern score: {scores['single_score']} / 12 possible criteria\n")
                f.write(f"Multiple pattern score: {scores['multiple_score']} / 13 possible criteria\n")
                f.write(f"Decision: {'Single' if scores['single_score'] > scores['multiple_score'] else 'Multiple'} pattern (higher score wins)\n")
                f.write("\n")
                f.write("Score calculation:\n")
                f.write("- Each pattern type has biological criteria (12 for Single, 13 for Multiple)\n")
                f.write("- Single pattern criteria: dominant high-het events, few total events, tight clustering, etc.\n")
                f.write("- Multiple pattern criteria: many events, mixed types, scattered distribution, etc.\n")
                f.write("- Score = count of criteria met for each pattern type\n")
                f.write("\n")
                f.write("Biological thresholds used:\n")
                f.write(f"- High heteroplasmy threshold: ≥{cfg.HIGH_HET_THRESHOLD:.0f}% (pathogenic significance)\n")
                f.write(f"- Significance threshold: ≥{cfg.NOISE_THRESHOLD:.0f}% (above noise level)\n")
                f.write(f"- Low heteroplasmy threshold: <{cfg.NOISE_THRESHOLD:.0f}% (likely artifacts)\n")
                f.write(f"- Multiple event threshold: >{cfg.MULTIPLE_EVENT_THRESHOLD} events (mouse model pattern)\n")
                f.write(f"- Dominant group threshold: ≥{cfg.DOMINANT_GROUP_FRACTION*100:.0f}% of events in main group\n")
                f.write(f"- Spatial clustering radius: {cfg.CLUSTER_RADIUS}bp (biologically relevant)\n")
                f.write("\n")
            
            # Biological interpretation
            f.write("Biological interpretation:\n")
            f.write("-" * 25 + "\n")
            if classification == "Single":
                f.write("Pattern consistent with:\n")
                f.write("- Single pathogenic deletion/duplication\n")
                f.write("- Classical mitochondrial disease patient profile\n")
                f.write("- Possible accompanying low-level artifacts\n")
                if criteria.get('max_heteroplasmy', 0) >= cfg.HIGH_HET_THRESHOLD:
                    f.write("- High heteroplasmy suggests functional impact\n")
            elif classification == "Multiple":
                f.write("Pattern consistent with:\n")
                f.write("- Multiple structural alterations\n") 
                f.write("- Mouse model of mtDNA maintenance defect\n")
                f.write("- Complex replication/repair dysfunction\n")
                if criteria.get('many_events', False):
                    f.write("- Extensive genomic instability\n")
            else:
                f.write("- Ambiguous or unusual pattern\n")
            f.write("\n")
            
            # Literature references
            f.write("Reference standards:\n")
            f.write("-" * 20 + "\n")
            f.write("Classification based on:\n")
            f.write("- Basu et al. PLoS Genet 2020 (MitoSAlt methodology)\n")
            f.write("- Patient samples: single high-heteroplasmy events (>35%)\n")
            f.write("- Mouse models: multiple low-heteroplasmy events (<3%)\n")
            f.write(f"- Clinical thresholds: >{cfg.HIGH_HET_THRESHOLD:.0f}% for pathogenic significance\n")
        
        logger.info(f"Enhanced analysis summary saved to {output_file}")
        logger.info(f"Event classification: {classification} ({reason})")
        if 'subtype' in criteria:
            logger.info(f"Pattern subtype: {criteria['subtype']}")

def write_summary(events, output_file, analysis_stats, classification_result, config=None, blacklist_regions=None):
    """
    Convenience function to write summary
    
    Args:
        events: DataFrame with events
        output_file: Output summary file path
        analysis_stats: Workflow statistics dict
        classification_result: Classification tuple
        config: ClassificationConfig instance
        blacklist_regions: List of blacklist regions
    """
    writer = SummaryWriter(config)
    writer.write(events, output_file, analysis_stats, classification_result, blacklist_regions)