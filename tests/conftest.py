from unittest.mock import Mock
from os.path import join as pjoin
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
	return '/data/GenomicTrackRepository/data/processed/hg38'