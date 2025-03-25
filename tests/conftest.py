from unittest.mock import Mock
from os.path import join as pjoin
import pytest
import os
import pathlib

@pytest.fixture
def project_test_data_directory(scope="session", autouse=True):
    return pathlib.Path(__file__).parent / "data"

@pytest.fixture
def setup_temp_dir(tmp_path, scope="session", autouse=True):
    """
    Fixture to set up a temporary dir for testing.
    """
    temp_dir = tmp_path / "test_dir"
    temp_dir.mkdir()
    return temp_dir