"""
Shared pytest fixtures for SaltShaker tests

Supports both development mode (python -m saltshaker) and installed mode (pip install -e .)
"""
import pytest
import pandas as pd
import json
from pathlib import Path
from argparse import Namespace
import sys


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    """Return path to fixtures directory"""
    return Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def ground_truth(fixtures_dir):
    """
    Load ground truth expected values from R script outputs
    """
    ground_truth_file = fixtures_dir / "expected" / "ground_truth.json"
    with open(ground_truth_file) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def reference_fasta(fixtures_dir):
    """
    Path to reference FASTA file (optional)
    """
    fasta_file = fixtures_dir / "inputs" / "human_mt_rCRS.fasta"
    return str(fasta_file) if fasta_file.exists() else None


@pytest.fixture(scope="session", autouse=True)
def setup_saltshaker_path():
    """
    Add repository root to Python path for development mode
    
    Structure:
      mitosalt-custom/              <- repo root (need to add this to sys.path)
      └── saltshaker/               <- package
          ├── __init__.py
          ├── cli/
          │   └── call.py
          └── tests/
              └── conftest.py       <- we are here
    """
    # From conftest.py, go up 2 levels to reach repo root
    # Path(__file__).parent = .../saltshaker/tests/
    # Path(__file__).parent.parent = .../saltshaker/ (package)
    # Path(__file__).parent.parent.parent = .../mitosalt-custom/ (repo root)
    
    repo_root = Path(__file__).parent.parent.parent
    
    # Add repo root to sys.path if not already there
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
        print(f"✓ Added {repo_root} to Python path")


@pytest.fixture(scope="session")
def saltshaker_output(fixtures_dir, reference_fasta, tmp_path_factory):
    """
    Run SaltShaker programmatically by calling the run() function from cli/call.py
    """
    # Import the run function from cli/call.py (not commands/call.py!)
    try:
        from saltshaker.cli.call import run
    except ImportError as e:
        pytest.fail(
            f"Cannot import saltshaker.cli.call.run: {e}\n"
            f"\nCurrent working directory: {Path.cwd()}\n"
            f"Python path: {sys.path}\n"
            f"\nMake sure you run from repository root:\n"
            f"  cd mitosalt-custom/\n"
            f"  pytest saltshaker/tests/ -v\n"
        )
    
    # Input files
    cluster_file = str(fixtures_dir / "inputs" / "test_clusters.cluster")
    breakpoint_file = str(fixtures_dir / "inputs" / "test_breakpoints.breakpoint")
    
    # Skip test if no reference FASTA
    if reference_fasta is None:
        pytest.skip("Reference FASTA required for SaltShaker call command")
    
    # Output directory
    output_dir = tmp_path_factory.mktemp("saltshaker_output")
    
    # Create args namespace matching what argparse would create
    args = Namespace(
        prefix='test_sample',
        output_dir=str(output_dir),
        cluster=cluster_file,
        breakpoint=breakpoint_file,
        reference=reference_fasta,
        genome_length=16569,
        ori_h_start=16081,
        ori_h_end=407,
        ori_l_start=5730,
        ori_l_end=5763,
        het_limit=0.01,
        flank_size=15,
        blacklist=None
    )
    
    # Run SaltShaker
    try:
        run(args)
    except Exception as e:
        import traceback
        pytest.fail(
            f"SaltShaker execution failed:\n"
            f"{traceback.format_exc()}"
        )
    
    # Find the display TSV output file
    display_tsv = output_dir / "test_sample.saltshaker_call.tsv"
    
    if not display_tsv.exists():
        pytest.fail(
            f"SaltShaker did not create output TSV: {display_tsv}\n"
            f"Files in output directory: {list(output_dir.iterdir())}"
        )
    
    # Load the output
    df = pd.read_csv(display_tsv, sep='\t')
    
    # Normalize column names
    df.columns = [col.replace('.', '_') for col in df.columns]
    
    # Mark that we have sequences
    df._has_sequences = True
    
    return df


# Pytest configuration
def pytest_configure(config):
    """Register custom markers"""
    config.addinivalue_line(
        "markers", "unit: Unit tests for individual functions"
    )
    config.addinivalue_line(
        "markers", "integration: Integration tests running full SaltShaker pipeline"
    )
    config.addinivalue_line(
        "markers", "coordinates: Tests validating coordinate calculations (main bug detection)"
    )
    config.addinivalue_line(
        "markers", "sequences: Tests requiring reference FASTA file"
    )
