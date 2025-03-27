from CRISPRito.Utils import greedy_clustering_first, greedy_clustering_incremental
from os.path import join as pjoin
import pytest
import os
import pandas as pd


def test_greedy_clustering_first():
	"""
	pytest -sv tests/unit/test_greedy_clustering.py::test_greedy_clustering_first
	"""	
	numbers = [1, 3, 2, 8, 10, 15, 12, 18, 20, 1000]
	range_threshold = 2
	result = greedy_clustering_first(numbers, range_threshold)
	assert result == {0: [1, 2, 3], 1: [8, 10], 2: [12], 3: [15], 4: [18, 20], 5: [1000]}

def test_greedy_clustering_incremental():
	"""
	pytest -sv tests/unit/test_greedy_clustering.py::test_greedy_clustering_incremental
	"""	
	numbers = [1, 3, 2, 8, 10, 15, 12, 18, 20, 1000]
	range_threshold = 2
	result = greedy_clustering_incremental(numbers, range_threshold)
	assert result == {0: [1, 2, 3], 1: [8, 10, 12], 2: [15], 3: [18, 20], 4: [1000]}