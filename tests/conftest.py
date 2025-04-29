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
def path_to_hg38_refseq(scope = 'session', autouse = True):
	"""
	"""
	return '/data/GenomicTrackRepository/data/processed/hg38/ncbiRefSeqCurated_expanded.csv'

@pytest.fixture
def path_to_hg38_encode(scope = 'session', autouse = True):
	"""
	"""
	return '/data/GenomicTrackRepository/data/processed/hg38/encode_ccre_all.csv'

@pytest.fixture
def sample_fasta_gz(scope = 'session', autouse = True):
    """Creates a temporary gzipped FASTA file for testing."""
    fasta_content = '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA'
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fa.gz") as tmp_fasta:
        with gzip.open(tmp_fasta.name, "wt") as f:
            f.write(fasta_content)
        return tmp_fasta.name  # Return the file path

@pytest.fixture
def path_to_hg38_gene_names(scope = 'session', autouse = True):
	return '/data/GenomicTrackRepository/data/external/hg38/gene_names.csv'

@pytest.fixture
def example_hg38_genome_size(scope = 'session', autouse = True):
	"""
	"""
	return {
		'chr1': 248956422,
		'chr10': 133797422,
		'chr11': 135086622,
		'chr12': 133275309,
		'chr13': 114364328,
		'chr14': 107043718,
		'chr15': 101991189,
		'chr16': 90338345,
		'chr17': 83257441,
		'chr18': 80373285,
		'chr19': 58617616,
		'chr2': 242193529,
		'chr20': 64444167,
		'chr21': 46709983,
		'chr22': 50818468,
		'chr3': 198295559,
		'chr4': 190214555,
		'chr5': 181538259,
		'chr6': 170805979,
		'chr7': 159345973,
		'chr8': 145138636,
		'chr9': 138394717,
		'chrM': 16569,
		'chrX': 156040895,
		'chrY': 57227415
		}



