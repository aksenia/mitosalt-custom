"""
SaltShaker Integration Tests

Validates SaltShaker output against R script ground truth.

These tests automatically run SaltShaker and identify:
- Coordinate calculation bugs (off-by-one errors)
- Event classification errors
- Heteroplasmy calculation issues
- Sequence extraction problems

Run: pytest tests/integration/ -v
"""
import pytest


# ============================================================================
# CLUSTER PRESENCE & BASIC VALIDATION
# ============================================================================

@pytest.mark.integration
def test_all_clusters_present(saltshaker_output, ground_truth):
    """
    Verify all 8 test clusters are processed
    """
    expected_clusters = set(ground_truth.keys())
    actual_clusters = set(saltshaker_output['cluster_id'].values)
    
    missing = expected_clusters - actual_clusters
    extra = actual_clusters - expected_clusters
    
    assert not missing, f"Missing clusters: {missing}"
    assert not extra, f"Unexpected clusters: {extra}"
    assert len(saltshaker_output) == 8, f"Expected 8 clusters, got {len(saltshaker_output)}"


# ============================================================================
# EVENT CLASSIFICATION
# ============================================================================

@pytest.mark.integration
def test_event_classification(saltshaker_output, ground_truth):
    """
    Test: del/dup classification matches R script
    
    Validates EventCaller's classification logic.
    All 8 clusters should have correct del/dup labels.
    """
    # DEBUG: Print what we're actually seeing
    print("\n" + "="*70)
    print("DEBUG: SaltShaker Output Columns:")
    print(saltshaker_output.columns.tolist())
    print("\nDEBUG: Cluster IDs and Classifications:")
    print(saltshaker_output[['cluster_id', 'final_event']].to_string())
    print("="*70)
    
    mismatches = []
    
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        print(f"\nChecking {cluster_id}:")
        print(f"  Expected: {expected['final_event']}")
        print(f"  Actual: {row['final_event']}")
        
        if row['final_event'] != expected['final_event']:
            mismatch = {
                'cluster': cluster_id,
                'expected': expected['final_event'],
                'actual': row['final_event']
            }
            mismatches.append(mismatch)
    
    if mismatches:
        msg = "\n".join([
            f"  {m['cluster']}: expected '{m['expected']}', got '{m['actual']}'"
            for m in mismatches
        ])
        pytest.fail(f"Event classification mismatches:\n{msg}")


# ============================================================================
# HETEROPLASMY CALCULATION
# ============================================================================

@pytest.mark.integration
def test_heteroplasmy_calculation(saltshaker_output, ground_truth):
    """
    Test: heteroplasmy percentages match R script (within ±0.01)
    
    Validates: alt_reads / (alt_reads + ref_reads) calculation
    """
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        assert row['heteroplasmy'] == pytest.approx(expected['heteroplasmy'], abs=0.01), \
            f"{cluster_id}: heteroplasmy expected {expected['heteroplasmy']:.4f}, " \
            f"got {row['heteroplasmy']:.4f}"


# ============================================================================
# READ COUNTS
# ============================================================================

@pytest.mark.integration
def test_read_counts(saltshaker_output, ground_truth):
    """
    Test: alt_reads and ref_reads match R script exactly
    
    Validates read counting logic.
    """
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        assert row['alt_reads'] == expected['alt_reads'], \
            f"{cluster_id}: alt_reads expected {expected['alt_reads']}, " \
            f"got {row['alt_reads']}"
        
        assert row['ref_reads'] == expected['ref_reads'], \
            f"{cluster_id}: ref_reads expected {expected['ref_reads']}, " \
            f"got {row['ref_reads']}"


# ============================================================================
# COORDINATE TESTS - MAIN BUG DETECTION
# ============================================================================

@pytest.mark.integration
@pytest.mark.coordinates
@pytest.mark.parametrize("cluster_id", [
    "cluster_001",  # Standard deletion (9499-13735)
    "cluster_002",  # Small deletion with microhomology (10349-10690)
    "cluster_003",  # Deletion with CCCTG microhomology (1465-4476)
    "cluster_006",  # Large deletion (5787.5-16072)
])
def test_deletion_final_start(saltshaker_output, ground_truth, cluster_id):
    """
    Test: final_start coordinate for deletions
    
    **THIS TEST CATCHES THE MAIN BUG**
    Expected to FAIL before fix (off by +1)
    Expected to PASS after fix
    """
    expected = ground_truth[cluster_id]
    assert expected['final_event'] == 'del', f"{cluster_id} should be deletion"
    
    row = saltshaker_output[saltshaker_output['cluster_id'] == cluster_id].iloc[0]
    
    difference = row['final_start'] - expected['final_start']
    
    assert row['final_start'] == pytest.approx(expected['final_start'], abs=0.1), \
        f"{cluster_id}: final_start expected {expected['final_start']}, " \
        f"got {row['final_start']} (OFF BY {difference:+.1f})"


@pytest.mark.integration
@pytest.mark.coordinates
@pytest.mark.parametrize("cluster_id", [
    "cluster_004",  # Large duplication (11039-11046)
    "cluster_005",  # D-loop crossing duplication (12299-3263)
    "cluster_007",  # D-loop crossing duplication (14815-499)
    "cluster_008",  # Small duplication (303-495)
])
def test_duplication_final_end(saltshaker_output, ground_truth, cluster_id):
    """
    Test: final_end coordinate for duplications
    
    **THIS TEST CATCHES THE MAIN BUG**
    Expected to FAIL before fix (off by +1)
    Expected to PASS after fix
    """
    expected = ground_truth[cluster_id]
    
    row = saltshaker_output[saltshaker_output['cluster_id'] == cluster_id].iloc[0]
    
    difference = row['final_end'] - expected['final_end']
    
    assert row['final_end'] == pytest.approx(expected['final_end'], abs=0.1), \
        f"{cluster_id}: final_end expected {expected['final_end']}, " \
        f"got {row['final_end']} (OFF BY {difference:+.1f})"


# ============================================================================
# EVENT SIZE CALCULATIONS
# ============================================================================

@pytest.mark.integration
def test_final_size(saltshaker_output, ground_truth):
    """
    Test: final_size calculations match R script
    
    Validates size calculation for final reported events.
    """
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        assert row['final_size'] == expected['final_size'], \
            f"{cluster_id}: final_size expected {expected['final_size']}, " \
            f"got {row['final_size']}"


@pytest.mark.integration
def test_deletion_size(saltshaker_output, ground_truth):
    """
    Test: del_size matches R script
    
    Validates initial deletion size calculation.
    """
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        assert row['del_size'] == expected['del_size'], \
            f"{cluster_id}: del_size expected {expected['del_size']}, " \
            f"got {row['del_size']}"


# ============================================================================
# COORDINATE RANGES
# ============================================================================

@pytest.mark.integration
def test_coordinate_ranges(saltshaker_output, ground_truth):
    """
    Test: del_start_range and del_end_range match R script
    
    Validates range reporting (e.g., "9495 - 9498")
    """
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        assert row['del_start_range'] == expected['del_start_range'], \
            f"{cluster_id}: del_start_range expected '{expected['del_start_range']}', " \
            f"got '{row['del_start_range']}'"
        
        assert row['del_end_range'] == expected['del_end_range'], \
            f"{cluster_id}: del_end_range expected '{expected['del_end_range']}', " \
            f"got '{row['del_end_range']}'"


# ============================================================================
# SEQUENCE-BASED TESTS (Require reference FASTA)
# ============================================================================

@pytest.mark.integration
@pytest.mark.sequences
@pytest.mark.parametrize("cluster_id", [
    "cluster_001",  # TTCGCAGGATTT
    "cluster_002",  # CCTAGCCCTA
    "cluster_003",  # CCCTG
    "cluster_004",  # AAAAAAACTCTACC
    "cluster_005",  # CTTA
    "cluster_007",  # CCCATCC
    "cluster_008",  # AACCCCC
])
def test_microhomology_detection(saltshaker_output, ground_truth, cluster_id):
    """
    Test: microhomology sequence detection
    
    **REQUIRES REFERENCE FASTA**
    Will be skipped if reference.fasta not present.
    
    Expected to FAIL before coordinate bugs are fixed.
    """
    if not saltshaker_output._has_sequences:
        pytest.skip("Requires reference FASTA for sequence extraction")
    
    expected = ground_truth[cluster_id]
    
    # cluster_006 has no microhomology
    if expected['seq'] in ['nan', 'NA']:
        pytest.skip(f"{cluster_id} has no microhomology")
    
    row = saltshaker_output[saltshaker_output['cluster_id'] == cluster_id].iloc[0]
    
    assert row['seq'] == expected['seq'], \
        f"{cluster_id}: microhomology expected '{expected['seq']}', " \
        f"got '{row['seq']}'"


@pytest.mark.integration
@pytest.mark.sequences
def test_flanking_sequence_seq1(saltshaker_output, ground_truth):
    """
    Test: seq1 (upstream flanking sequence) matches R script
    
    **REQUIRES REFERENCE FASTA**
    
    This catches coordinate bugs affecting upstream extraction.
    Expected to FAIL before coordinate bugs are fixed.
    """
    if not saltshaker_output._has_sequences:
        pytest.skip("Requires reference FASTA for sequence extraction")
    
    mismatches = []
    
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        if row['seq1'] != expected['seq1']:
            mismatch = {
                'cluster': cluster_id,
                'expected': expected['seq1'],
                'actual': row['seq1']
            }
            mismatches.append(mismatch)
    
    if mismatches:
        msg = "\n".join([
            f"  {m['cluster']}:\n"
            f"    Expected: {m['expected']}\n"
            f"    Got:      {m['actual']}"
            for m in mismatches
        ])
        pytest.fail(
            f"seq1 (upstream flanking) mismatches:\n{msg}\n"
            f"This indicates coordinate bug affecting sequence extraction"
        )


@pytest.mark.integration
@pytest.mark.sequences
def test_flanking_sequence_seq2(saltshaker_output, ground_truth):
    """
    Test: seq2 (downstream flanking sequence) matches R script
    
    **REQUIRES REFERENCE FASTA**
    
    This catches coordinate bugs affecting downstream extraction.
    Expected to FAIL before coordinate bugs are fixed.
    """
    if not saltshaker_output._has_sequences:
        pytest.skip("Requires reference FASTA for sequence extraction")
    
    mismatches = []
    
    for _, row in saltshaker_output.iterrows():
        cluster_id = row['cluster_id']
        expected = ground_truth[cluster_id]
        
        if row['seq2'] != expected['seq2']:
            mismatch = {
                'cluster': cluster_id,
                'expected': expected['seq2'],
                'actual': row['seq2']
            }
            mismatches.append(mismatch)
    
    if mismatches:
        msg = "\n".join([
            f"  {m['cluster']}:\n"
            f"    Expected: {m['expected']}\n"
            f"    Got:      {m['actual']}"
            for m in mismatches
        ])
        pytest.fail(
            f"seq2 (downstream flanking) mismatches:\n{msg}\n"
            f"This indicates coordinate bug affecting sequence extraction"
        )
