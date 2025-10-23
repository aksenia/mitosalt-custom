# SAltShaker

A Python package for classifying and visualizing mitochondrial structural variants from MitoSAlt pipeline output.

<p align="center">
  <img src="assets/logo.png" alt="SAltShaker Logo" width="250"/>
</p>

## Overview

SAltShaker is a Python port and extension of the original MitoSAlt `delplot.R` visualization script. The package provides three modular commands for a flexible analysis workflow:

- **Event calling**: Direct Python port of the original R script's deletion/duplication classification logic
- **Pattern classification**: Rule-based decision tree to distinguish single from mouse multiple type of events
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
    --prefix sample \                       # Sample identifier
    --output-dir results/ \                 # Output directory
    -c sample.cluster \                     # cluster file from MitoSAlt
    -p sample.breakpoint \                  # breakpoint file from MitoSAlt
    -r reference.fasta \                    # mitochondrial reference genome
    -g 16569 \                              # genome length
    --ori-h-start 16081 --ori-h-end 407 \   # heavy strand origin
    --ori-l-start 5730 --ori-l-end 5763 \   # light strand origin
    -H 0.01 \                               # heteroplasmy threshold
    -f 15 \                                 # flanking sequence size (bp)
    --blacklist                             # Optional: enable with default MT blacklist, OR
    --blacklist custom_bl.bed               # Optional: enable with custom BED file
Out
```

**Outputs:**

- `results/sample.saltshaker_call.tsv` - Human-readable results matching original R script output format
- `results/sample.saltshaker_call_metadata.tsv` - Metadata for downstream processing (classify/plot)


**Algorithm (from original R script):**

The  logic is a direct port of the original MitoSAlt R implementation:

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
    --prefix sample1 \              # Sample identifier (matches call output)
    --input-dir results/ \          # Input directory with .saltshaker_call.tsv
    --output-dir results/ \         # Output directory (default: same as input-dir)
    --blacklist                     # Optional: enable with default MT blacklist, OR
    --blacklist custom_bl.bed       # Optional: enable with custom BED file
   --vcf \                         # Optional: also output VCF format
    --high-het 20 \                 # Optional: high heteroplasmy threshold % (default: 20)
    --noise 1.0 \                   # Optional: noise threshold % (default: 1.0)
    --radius 600 \                  # Optional: spatial clustering radius bp (default: 600)
    --multiple-threshold 10 \       # Optional: event count for Multiple pattern (default: 10)
    --dominant-fraction 0.70        # Optional: fraction for dominant group (default: 0.70)
```

**Outputs:**

- `results/sample.saltshaker_classify.txt` - Detailed analysis report with classification reasoning
- `results/sample.saltshaker_classify_metadata.tsv` - Events with spatial group assignments (for plotting)
- `results/sample.vcf` - VCF format with groups and heteroplasmy (if `--vcf` specified)


**Classification criteria:**

**Single pattern** (patient-like):

- One or few high-heteroplasmy events (≥20%)
- Dominant spatial group (≥70% of events)
- Few total events (≤10)
- Consistent with pathogenic single deletion/duplication

**Multiple pattern** (mouse model-like):

- Many events (>10)
- Dispersed spatial distribution (no dominant group)
- No high-heteroplasmy events
- Consistent with mtDNA maintenance defects

**Spatial grouping:**
Events within 600bp (configurable) are grouped together. All classification thresholds can be customized via CLI (see `saltshaker classify --help`).

### 3. `saltshaker plot` - Visualization

Generates circular genome plots based on the original R script visualization with enhanced spatial grouping. 

**Usage:**

```bash
saltshaker plot \
    --prefix sample1 \
    --input-dir results/ \
    --output-dir results/plots/ \  # Optional: default is input-dir
    --genes \                       # Optional: enable with default MT genes, OR
    --genes custom_genes.bed \      # Optional: enable with custom BED file
    --blacklist \                   # Optional: enable with default MT blacklist, OR
    --blacklist custom_bl.bed \     # Optional: enable with custom BED file
    --figsize 16 10 \               # Optional: width height (default: 16 10)
    --direction clockwise \         # Optional: clockwise or counterclockwise (default: counterclockwise)
    --del-color red \               # Optional: red or blue (default: blues)
    --dup-color blue                # Optional: red or blue (default: red)

```

**Output:**

- `results/plot/sample.saltshaker.png` - Circular genome visualization

**Visualization features:**

- Circular genome with arc-based event display
- Heteroplasmy gradient coloring
- Spatial grouping for non-overlapping events
- Optional gene annotations with track and labels
- Optional blacklist region marking (BL-crossing events in lime-green)
- Configurable colors and polar direction


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

### Display TSV (`{prefix}.saltshaker_call.tsv`)

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

### Metadata files

**Call metadata** (`{prefix}.saltshaker_call_metadata.tsv`):
Internal format preserving all columns for downstream processing. Contains metadata header with genome length.

**Classify metadata** (`{prefix}.saltshaker_classify_metadata.tsv`):
Internal format with additional `group` column for spatial group assignments. Used by plot command.

### Analysis summary (`{prefix}.saltshaker_classify.txt`)

Human-readable analysis report including:

- Pattern classification (Single/Multiple) with reasoning
- Event statistics and heteroplasmy distribution
- Spatial clustering metrics
- Classification criteria scores

### VCF format (`{prefix}.vcf`)

Standard VCF 4.3 format with structural variant fields:

- `SVTYPE`: DEL or DUP
- `END`: Variant end position
- `SVLEN`: Variant length
- `HF`: Heteroplasmy fraction (0-1)
- `GROUP`: Spatial group identifier
- `CLUSTER`: Original cluster ID
- `DLOOP`: Flag for D-loop crossing
- `BLCROSS`: Flag for blacklist crossing

### Circular Plot (`{prefix}.saltshaker.png`)

## Complete workflow example

**Single-sample pipeline:**

```bash
# Step 1: Call events from MitoSAlt output (R script port)
saltshaker call \
    --prefix sample1 \
    --output-dir results/ \
    -c sample1.cluster -p sample1.breakpoint \
    -r reference.fasta \
    -g 16569 --ori-h-start 16081 --ori-h-end 407 \
    --ori-l-start 5730 --ori-l-end 5763 \
    --blacklist

# Step 2: Classify pattern and perform spatial grouping (extended analysis)
saltshaker classify \
    --prefix sample1 \
    --input-dir results/ \
    --blacklist \
    --vcf

# Step 3: Generate visualization (enhanced R script plotting)
saltshaker plot \
    --prefix sample1 \
    --input-dir results/ \
    --output-dir results/plot/ \
    --blacklist \
    --genes
```

**Batch processing multiple samples:**

```bash
for sample in sample1 sample2 sample3; do
    mkdir -p results/${sample}
    # Call events
    saltshaker call --prefix ${sample }--output-dir results/${sample} \
        -c ${sample}.cluster -p ${sample}.breakpoint \
        -r reference.fasta
        -g 16569 --ori-h-start 16081 --ori-h-end 407 \
        --ori-l-start 5730 --ori-l-end 5763 \
        --blacklist
    # Classify events
    saltshaker classify --prefix ${sample} \
        --input-dir results/${sample} \
        --blacklist
    # Plot events
    saltshaker plot --prefix ${sample} \
        --input-dir results/${sample} \
        --blacklist \
        --genes
done
```

## Configuration

Default classification thresholds are defined in `saltshaker/config.py`:

```python
# Heteroplasmy thresholds
HIGH_HET_THRESHOLD = 20.0        # High heteroplasmy threshold (%), --high-het
NOISE_THRESHOLD = 1.0            # Noise threshold (%), --noise

# Spatial clustering
CLUSTER_RADIUS = 600             # Spatial grouping radius (bp), --radius
MIN_CLUSTER_SIZE = 2             # Minimum events per cluster (not configurable via CLI)

# Pattern classification
MULTIPLE_EVENT_THRESHOLD = 10    # Event count for Multiple pattern, --multiple-threshold
DOMINANT_GROUP_FRACTION = 0.70   # Fraction for dominant group (70%), --dominant-fraction
```

These can be customized by modifying the configuration file or via CLI arguments (see `saltshaker classify --help`).

## Command reference

### Global options

All commands support:

- `-h, --help`: Show help message
- `-b, --blacklist FILE`: BED file with regions to exclude

### `call` command

**Required:**

- `--prefix STR`: Sample prefix for output files
- `-c, --cluster FILE`: Cluster file from MitoSAlt
- `-p, --breakpoint FILE`: Breakpoint file from MitoSAlt
- `-r, --reference FILE`: Reference genome FASTA
- `-g, --genome-length INT`: Mitochondrial genome length
- `--ori-h-start INT`: Heavy strand origin start
- `--ori-h-end INT`: Heavy strand origin end
- `--ori-l-start INT`: Light strand origin start
- `--ori-l-end INT`: Light strand origin end

**Optional:**

- `--output-dir DIR`: Output directory (default: .)
- `-H, --het-limit FLOAT`: Heteroplasmy threshold (default: 0.01)
- `-f, --flank-size INT`: Flanking sequence size in bp (default: 15)
- `--blacklist [FILE]`: Enable blacklist regions (default: built-in MT blacklist; optional: custom BED file)

### `classify` command

**Required:**

- `--prefix STR`: Sample prefix (matches call output)
- `--input-dir DIR`: Input directory containing saltshaker_call_metadata.tsv from call

**Optional:**

- `--output-dir DIR`: Output directory (default: input-dir)
- `--blacklist [FILE]`: Enable blacklist regions (default: built-in MT blacklist; optional: custom BED file)
- `--vcf`: Also output VCF format
- `--high-het FLOAT`: High heteroplasmy threshold % (default: 20)
- `--noise FLOAT`: Noise threshold % (default: 1)
- `--radius INT`: Spatial clustering radius bp (default: 600)
- `--multiple-threshold INT`: Event count for Multiple pattern (default: 10)
- `--dominant-fraction FLOAT`: Fraction for dominant group (default: 0.70)

### `plot` command

**Required:**

- `--prefix STR`: Sample prefix (matches classify output)
- `--input-dir DIR`: Input directory containing saltshaker_classify_metadata.tsv from classify

**Optional:**

- `--output-dir DIR`: Output directory (default: input-dir)
- `--genes [FILE]`: Enable gene annotations (default: built-in MT genes; optional: custom BED file)
- `--blacklist [FILE]`: Enable blacklist regions (default: built-in MT blacklist; optional: custom BED file)
- `--figsize WIDTH HEIGHT`: Figure dimensions (default: 16 10)
- `--direction STR`: Plot direction - 'clockwise' or 'counterclockwise' (default: counterclockwise, MitoSAlt original)
- `--del-color STR`: Deletion color - 'red' or 'blue' (default: blue, MitoSAlt original)
- `--dup-color STR`: Duplication color - 'red' or 'blue' (default: red, MitoSAlt original)

## Dependencies

```bash
pandas>=1.3.0
numpy>=1.20.0
matplotlib>=3.3.0
biopython>=1.78
```

## Package structure

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
└── data/
    ├── __init__.py      # Default file paths
    ├── gencode.v49.annotation.MT_genes.bed      # Default MT gene annotations
    └── mt_blacklist_regions.bed                 # Default MT blacklist regions
docs/
└── classification_algorithm.md  # Detailed classification algorithm documentation
```

## Documentation

- [Classification algorithm](docs/classification_algorithm.md) - Detailed explanation of the Single vs Multiple pattern classification logic
