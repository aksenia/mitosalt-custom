# SaltShaker Classification algorithm

## Rationale

Based on Basu et al. PLoS Genetics 2020, the algorithm distinguishes two patterns observed in mitochondrial DNA studies:

**SINGLE PATTERN:**

- Few events with high heteroplasmy (≥20%)
- One dominant deletion/duplication
- Spatially clustered breakpoints (≥70% in one group)
- Reflects: clonal expansion of single pathogenic event

**MULTIPLE PATTERN:**

- Many low-heteroplasmy events
- Scattered across genome (multiple spatial groups)
- No dominant high-het event
- Reflects: ongoing mtDNA maintenance defect

---

## Classification parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `HIGH_HET_THRESHOLD` | 20.0% | High heteroplasmy threshold |
| `NOISE_THRESHOLD` | 1.0% | Noise level threshold |
| `MIN_CLUSTER_SIZE` | 2 | Minimum events for spatial grouping |
| `MULTIPLE_EVENT_THRESHOLD` | 10 | Threshold for "many events" |
| `DOMINANT_GROUP_FRACTION` | 0.70 | Threshold for dominant group (70%) |
| `CLUSTER_RADIUS` | 600 bp | Spatial grouping radius |

**CLI usage:**

```bash
saltshaker classify \
    --prefix sample \.      # prefix used in call command
    --input-dir results/ \  # outputs folder from call command
    --high-het 20 \
    --noise 1.0 \
    --radius 600 \
    --multiple-threshold 10 \
    --dominant-fraction 0.70 \
    -b blacklist.bed \    # will exclude blacklist-crossing events from classification
    --vcf
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

```bash
START
  │
  ├─ All events < NOISE (0.7%)?
  │   └─ YES → "No significant events" (noise only)
  │
  ├─ Events < MIN_CLUSTER_SIZE (2) AND max_het < HIGH_HET (20%)?
  │   └─ YES → "No significant events" (below clinical threshold)
  │
  ├─ Perform spatial grouping
  │   │
  │   ├─ has_high_het (≥20%) AND (few_events ≤10 OR dominant_group ≥70%)?
  │   │   └─ YES → "Single"
  │   │
  │   ├─ many_events (>10) OR (dispersed <70% AND no_high_het)?
  │   │   └─ YES → "Multiple"
  │   │
  │   └─ AMBIGUOUS → Tiebreaker by max_het
  │       ├─ max_het ≥ 20% → "Single"
  │       └─ max_het < 20% → "Multiple"
  │
END
```

## Classification outcomes

### No significant events

- All events < 0.7% (noise threshold)
- **OR** <2 events AND max heteroplasmy <20%
- **Result:** Not clinically actionable

### Single pattern

- High heteroplasmy (≥20%) **AND** (≤10 events **OR** ≥70% in dominant group)
- **Result:** Pathogenic single deletion/duplication

### Multiple pattern

- Many events (>10) **OR** (Dispersed <70% **AND** no high-het)
- **Result:** Complex mtDNA maintenance defect

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