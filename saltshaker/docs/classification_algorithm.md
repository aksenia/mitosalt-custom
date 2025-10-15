# SaltShaker Classification algorithm

## Rationale

Based on Basu et al. PLoS Genetics 2020, the algorithm distinguishes two patterns observed in mitochondrial DNA studies:

**SINGLE PATTERN (patient-like):**

- Few events with high heteroplasmy (≥20%)
- One dominant deletion/duplication
- Spatially clustered breakpoints (≥70% in one group)
- Reflects: clonal expansion of single pathogenic event

**MULTIPLE PATTERN (mouse model-like):**

- Many low-heteroplasmy events
- Scattered across genome (multiple spatial groups)
- No dominant high-het event
- Reflects: ongoing mtDNA maintenance defect

---

## Classification parameters

All parameters are configurable via CLI with sensible defaults:

```python
HIGH_HET_THRESHOLD = 20.0          # High heteroplasmy threshold (%)
NOISE_THRESHOLD = 1.0              # Below this = artifacts (%)
CLUSTER_RADIUS = 600               # Spatial grouping radius (bp)
MIN_CLUSTER_SIZE = 2               # Minimum events per cluster
MULTIPLE_EVENT_THRESHOLD = 10      # Event count for Multiple pattern
DOMINANT_GROUP_FRACTION = 0.70     # Fraction in dominant group for Single (70%)
```

**CLI usage:**

```bash
saltshaker classify \
    -i events.tsv \
    -o summary.txt \
    --high-het 20 \
    --noise 1.0 \
    --radius 600 \
    --multiple-threshold 10 \
    --dominant-fraction 0.70
```

---

## Classification logic

### Step 1: Filter and count

```python
# Filter blacklist-crossing events
clean_events = [events not crossing blacklist regions]

# Count significant events (≥ NOISE_THRESHOLD)
significant_events = [events ≥ NOISE_THRESHOLD]
high_het_events = [events ≥ HIGH_HET_THRESHOLD]

# Check for no significant events
IF significant_count == 0:
    RETURN "No significant events"
```

### Step 2: Spatial grouping

```python
# Group events by proximity (CLUSTER_RADIUS)
groups = spatial_grouping(events, CLUSTER_RADIUS)

# Find dominant group
dominant_group = largest_group
dominant_fraction = dominant_group_size / total_events
```

### Step 3: Pattern classification

**SINGLE PATTERN if:**

```python
(has_high_het_events AND few_events) 
    OR 
(has_high_het_events AND dominant_group_pattern)

WHERE:
    has_high_het_events = len(high_het_events) > 0
    few_events = total_events ≤ MULTIPLE_EVENT_THRESHOLD
    dominant_group_pattern = dominant_fraction ≥ DOMINANT_GROUP_FRACTION
```

**MULTIPLE PATTERN if:**

```python
many_events 
    OR 
(dispersed AND no_high_het)

WHERE:
    many_events = total_events > MULTIPLE_EVENT_THRESHOLD
    dispersed = dominant_fraction < DOMINANT_GROUP_FRACTION
    no_high_het = len(high_het_events) == 0
```

**Ambiguous cases (tie-breaker):**

```python
IF max_heteroplasmy ≥ HIGH_HET_THRESHOLD:
    RETURN "Single"
ELSE:
    RETURN "Multiple"
```

---

## Output: Criteria dictionary

The classifier returns a criteria dictionary with classification metrics:

### Core event metrics

```python
'total_events'              # Events used for classification (after blacklist filtering)
'total_raw_events'          # All events before filtering
'blacklist_filtered'        # Number of events excluded (blacklist crossing)
'significant_count'         # Events ≥ NOISE_THRESHOLD
'high_het_count'            # Events ≥ HIGH_HET_THRESHOLD
'max_heteroplasmy'          # Maximum heteroplasmy value (%)
'median_heteroplasmy'       # Median heteroplasmy value (%)
```

### Event type counts

```python
'del_count'                 # Deletions (≥ NOISE_THRESHOLD)
'dup_count'                 # Duplications (≥ NOISE_THRESHOLD)
```

### Spatial organization

```python
'total_groups'              # Number of spatial groups formed
'dominant_group_fraction'   # Fraction of events in largest group (0-1)
'dominant_group_count'      # Number of events in largest group
'group_analysis'            # Detailed list of all groups with coordinates
```

### Classification result

```python
'subtype'                   # Pattern subtype (see below)
```

---

## Pattern subtypes

**For Single pattern:**

- `"Classic single deletion/duplication"` - One high-het event
- `"Few major events"` - 2-3 high-het events

**For Multiple pattern:**

- `"Complex multiple (mouse model-like)"` - >50 events
- `"Mixed deletion-duplication"` - Both event types present
- `"Multiple single-type events"` - Many events of one type

---