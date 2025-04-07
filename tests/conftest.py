from unittest.mock import Mock
from os.path import join as pjoin
import gzip
import tempfile
import pytest
import os
import pathlib

@pytest.fixture
def project_test_data_directory(scope="session", autouse=True):
    return pjoin(pathlib.Path(__file__).parent, "data")

@pytest.fixture
def setup_temp_dir(tmp_path, scope="session", autouse=True):
    """
    Fixture to set up a temporary dir for testing.
    """
    temp_dir = pjoin(tmp_path, "test_dir")
    os.makedirs(temp_dir)
    return temp_dir


@pytest.fixture
def path_to_hg38_genome(scope = 'session', autouse = True):
	"""
	"""
	return '/data/GenomicTrackRepository/data/processed/hg38/hg38.fasta.gz'


@pytest.fixture
def sample_fasta_gz(scope = 'session', autouse = True):
    """Creates a temporary gzipped FASTA file for testing."""
    fasta_content = '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA'
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fa.gz") as tmp_fasta:
        with gzip.open(tmp_fasta.name, "wt") as f:
            f.write(fasta_content)
        return tmp_fasta.name  # Return the file path


@pytest.fixture
def example_hg38_genome_size(scope = 'session', autouse = True):
	"""
	"""
	