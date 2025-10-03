# SAltShaker

A Python package for classifying and visualizing mitochondrial structural variants from MitoSAlt pipeline output.

## Overview

SaltShaker is a Python port and extension of the original MitoSAlt `delplot.R` visualization script. The package provides three modular commands for a flexible analysis workflow:

- **Event calling**: Direct Python port of the original R script's deletion/duplication classification logic
- **Pattern classification**: Enhanced analysis distinguishing patient-like single events from mouse model-like multiple events
- **Visualization**: Circular genome plotting based on the original R script visualization with spatial grouping enhancements

The core deletion/duplication classification algorithm (`saltshaker call`) faithfully replicates the original R implementation. The additional pattern classification (`saltshaker classify`) extends this with literature-based criteria to distinguish between pathogenic single high-heteroplasmy events and mtDNA maintenance defect patterns following [Basu et al. PLoS Genetics 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009242).

## Installation

```bash
# From source
git clone https://github.com/aksenia/mitosalt_custom.git
cd mitosalt_custom
pip install -e .
```

## Commands

### 1. `saltshaker call` - Event calling (R script port)

**Python port of the original MitoSAlt R script event classification logic.** Classifies detected breakpoint clusters as deletions or duplications based on replication origin overlap, following the exact algorithm from `delplot.R`.

**Usage:**
```bash
saltshaker call \
    -g 16569 \                           # genome length
    --ori-h-start 16081 --ori-h-end 407 \  # heavy strand origin
    --ori-l-start 5730 --ori-l-end 5763 \  # light strand origin
    -c sample.cluster \                  # cluster file from MitoSAlt
    -p sample.breakpoint \               # breakpoint file from MitoSAlt
    -r reference.fasta \                 # mitochondrial reference genome
    -o sample.tsv \                      # output display TSV
    -H 0.01 \                            # heteroplasmy threshold
    -f 15 \                              # flanking sequence size (bp)
    -b blacklist.bed                     # optional: regions to exclude
```

**Outputs:**
- `sample.tsv` - Human-readable results matching original R script output format
- `sample.intermediate.tsv` - Internal format for downstream processing (classify/plot)

**Algorithm (from original R script):**

The classification logic is a direct port of the original MitoSAlt R implementation:

1. **Data loading**: Parses cluster and breakpoint files, merges data to link clusters with D-loop crossing information
2. **Deletion size calculation**: Handles circular genome wraparound for events crossing position 1
3. **Origin-based classification**: Events overlapping replication origins (OriH or OriL) are classified as duplications; non-overlapping events are deletions
4. **Coordinate handling**: Implements the R script's coordinate swapping logic for D-loop crossing events
5. **Flanking sequence analysis**: Uses Biostrings-equivalent pattern matching to find microhomology sequences near breakpoints

This approach identifies the arc complementary to the actual structural change, consistent with the biological interpretation that origin-overlapping events represent duplications of the non-deleted arc.

**Original R script functionality preserved:**
- Exact deletion/duplication classification algorithm
- D-loop crossing detection and coordinate handling
- Flanking sequence extraction and microhomology analysis
- Output TSV format and column names
- Heteroplasmy calculation and filtering

**Python enhancements in call step:**
- Blacklist region detection and flagging

### 2. `saltshaker classify` - Pattern classification

**Extended analysis beyond the original R script.** Performs spatial grouping and classifies the overall pattern as Single or Multiple based on heteroplasmy distribution and event clustering.

**Usage:**
```bash
saltshaker classify \
    -i sample.intermediate.tsv \     # input from call step
    -o sample_summary.txt \          # analysis summary
    -b blacklist.bed \               # optional: blacklist regions
    --vcf                            # optional: output VCF format
```

**Outputs:**
- `sample_summary.txt` - Detailed analysis report with classification reasoning
- `sample_summary.classified.tsv` - Events with spatial group assignments (for plotting)
- `sample_summary.vcf` - VCF format with groups and heteroplasmy (if `--vcf` specified)

**Classification criteria:**

**Single Pattern** (patient-like):
- One or few high-heteroplasmy events (≥30%)
- Dominant cluster with high total heteroplasmy
- Spatially clustered events
- Consistent with pathogenic single deletion/duplication

**Multiple Pattern** (mouse model-like):
- Many low-heteroplasmy events (<30%)
- High event count (typically >20)
- Dispersed spatial distribution
- Consistent with mtDNA maintenance defects

**Spatial grouping:**
Events within 400bp (configurable) are grouped together. 

### 3. `saltshaker plot` - Visualization

Generates circular genome plots based on the original R script visualization with enhanced spatial grouping. Creates publication-quality figures showing spatial distribution of events with color-coded heteroplasmy intensity.

**Usage:**

```bash
saltshaker plot \
    -i sample_summary.classified.intermediate.tsv \  # input from classify
    -o sample_plot.png \                            # output plot
    -b blacklist.bed \                              # optional: blacklist regions
    --figsize 16 10                                 # figure size (width height)
```

**Output:**
- Circular genome visualization with:
  - Blue arcs: deletions
  - Red arcs: duplications  
  - Color intensity: heteroplasmy level
  - Spatial grouping: events grouped by proximity
  - Blacklisted regions: marked if provided

**Visualization features from original R script:**
- Circular genome representation
- Color coding by event type (blue/red)
- Heteroplasmy intensity mapping
- Arc-based event display
- Coordinate labeling

**Enhanced features:**
- Spatial group organization for non-overlapping display
- Size-based layering (largest to smallest radius)
- Blacklist region visualization
- Configurable figure dimensions

## Input files

### Required (from MitoSAlt pipeline)

**Cluster file** (`.cluster`)
- Tab-separated clustered breakpoint data
- Generated by MitoSAlt's clustering step
- Columns: cluster ID, read counts, positions, heteroplasmy levels

**Breakpoint file** (`.breakpoint`)
- Tab-separated raw breakpoint data
- Generated by MitoSAlt's breakpoint detection step
- Columns: read names, positions, D-loop crossing flags

**Reference genome** (`.fasta`)
- Mitochondrial reference sequence
- Used for flanking sequence extraction and coordinate validation

### Optional

**Blacklist file** (`.bed`)
- BED format regions to exclude (e.g., artifacts, repetitive sequences)
- Format: chromosome, start position, end position, name

## Output formats

### TSV files

**Display TSV** (`sample.tsv`):
Human-readable format matching original R script output with columns:
- `sample`: Sample identifier
- `cluster.id`: Cluster identifier
- `alt.reads`, `ref.reads`: Read counts
- `heteroplasmy`: Heteroplasmy percentage
- `del.start.range`, `del.end.range`: Coordinate ranges
- `del.size`: Event size in base pairs
- `final.event`: Event type (del/dup)
- `final.start`, `final.end`: Final coordinates
- `blacklist_crossing`: Flag for blacklist overlap
- `seq1`, `seq2`, `seq`: Flanking sequences and microhomology

**Intermediate TSV** (`.intermediate.tsv`):
Internal format preserving all columns for downstream processing. Contains metadata header with genome length.

**Classified TSV** (`.classified.tsv`):
Intermediate format with additional `group` column for spatial group assignments.

### VCF format

Standard VCF 4.3 format with structural variant fields:
- `SVTYPE`: DEL or DUP
- `END`: Variant end position
- `SVLEN`: Variant length
- `HF`: Heteroplasmy fraction (0-1)
- `GROUP`: Spatial group identifier
- `CLUSTER`: Original cluster ID
- `DLOOP`: Flag for D-loop crossing
- `BLCROSS`: Flag for blacklist crossing

### Summary text

Analysis report including:
- Pattern classification (Single/Multiple) with reasoning
- Event statistics and heteroplasmy distribution
- Spatial clustering metrics
- Classification criteria scores

## Complete workflow example

```bash
# Step 1: Call events from MitoSAlt output (R script port)
saltshaker call \
    -g 16569 --ori-h-start 16081 --ori-h-end 407 \
    --ori-l-start 5730 --ori-l-end 5763 \
    -c sample.cluster -p sample.breakpoint -r reference.fasta \
    -o results/sample.tsv \
    -H 0.01 -f 15 -b blacklist.bed

# Step 2: Classify pattern and perform spatial grouping (extended analysis)
saltshaker classify \
    -i results/sample.intermediate.tsv \
    -o results/sample_summary.txt \
    -b blacklist.bed --vcf

# Step 3: Generate visualization (enhanced R script plotting)
saltshaker plot \
    -i results/sample_summary.classified.tsv \
    -o results/sample_plot.png \
    -b blacklist.bed --figsize 16 10
```

## Configuration

Default thresholds are defined in `saltshaker/config.py`:

```python
# Classification thresholds
HIGH_HETEROPLASMY_THRESHOLD = 0.30      # 30%
SIGNIFICANT_HETEROPLASMY_THRESHOLD = 0.05  # 5%
CLUSTER_RADIUS = 500                     # bp
MIN_CLUSTER_SIZE = 1

# Scoring weights for pattern classification
SINGLE_HIGH_HET_WEIGHT = 0.35
SINGLE_DOMINANT_GROUP_WEIGHT = 0.30
SINGLE_LOW_COUNT_WEIGHT = 0.20
SINGLE_CLUSTERING_WEIGHT = 0.15

MULTIPLE_COUNT_WEIGHT = 0.35
MULTIPLE_LOW_HET_WEIGHT = 0.30
MULTIPLE_DISPERSION_WEIGHT = 0.20
MULTIPLE_NO_DOMINANT_WEIGHT = 0.15
```

These can be customized by modifying the configuration file.

## Command reference

### Global options

All commands support:
- `-h, --help`: Show help message
- `-b, --blacklist FILE`: BED file with regions to exclude

### `call` command

**Required:**
- `-g, --genome-length INT`: Mitochondrial genome length
- `--ori-h-start INT`: Heavy strand origin start
- `--ori-h-end INT`: Heavy strand origin end
- `--ori-l-start INT`: Light strand origin start
- `--ori-l-end INT`: Light strand origin end
- `-c, --cluster FILE`: Cluster file from MitoSAlt
- `-p, --breakpoint FILE`: Breakpoint file from MitoSAlt
- `-r, --reference FILE`: Reference genome FASTA
- `-o, --output FILE`: Output TSV file

**Optional:**
- `-H, --het-limit FLOAT`: Heteroplasmy threshold (default: 0.01)
- `-f, --flank-size INT`: Flanking sequence size in bp (default: 15)

### `classify` command

**Required:**
- `-i, --input FILE`: Intermediate TSV from call step
- `-o, --output FILE`: Output summary text file

**Optional:**
- `--vcf`: Also output VCF format

### `plot` command

**Required:**
- `-i, --input FILE`: Classified intermediate TSV from classify step
- `-o, --output FILE`: Output PNG file

**Optional:**
- `--figsize WIDTH HEIGHT`: Figure dimensions (default: 16 10)

## Dependencies

```
pandas>=1.3.0
numpy>=1.20.0
matplotlib>=3.3.0
biopython>=1.78
```

## Package Structure

```bash
saltshaker/
├── __init__.py
├── __main__.py          # CLI entry point with subcommands
├── config.py            # Configuration and thresholds
├── event_caller.py      # Event calling (R script port)
├── classifier.py        # Pattern classification (extended)
├── spatial.py           # Spatial grouping (extended)
├── visualizer.py        # Circular plotting (R script based)
├── utils.py             # Utility functions
├── cli/
│   ├── call.py          # Call subcommand
│   ├── classify.py      # Classify subcommand
│   └── plot.py          # Plot subcommand
└── io/
    ├── readers.py       # File input (blacklist BED files)
    ├── writers.py       # TSV and summary output
    └── vcf_writer.py    # VCF format output
```