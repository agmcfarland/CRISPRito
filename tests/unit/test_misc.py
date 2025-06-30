import pytest
import pandas as pd
from CRISPRito.Utils import df_long_to_wide, sliding_windows

def test_df_long_to_wide():
	"""
	"""
	df = pd.DataFrame({
		"sample": ["s1", "s1", "s2", "s2", "s2"],
		"method": ["m1", "m1", "m1", "m1", "m1"],
		"cut_cluster": [100, 101, 100, 101, 102],
	})

	df_test = df[["sample", "cut_cluster"]]

	result = df_long_to_wide(df_test, to_rows="sample", to_columns="cut_cluster", column_prefix="cc")

	assert result.shape == (2, 4)  # 2 samples, 3 clusters + 1 index
	assert "cc_100" in result.columns
	assert result.loc[result["sample"] == "s2", "cc_102"].values[0] == 1
	assert result.loc[result["sample"] == "s1", "cc_101"].values[0] == 1

def test_sliding_windows_basic():
	# Sequence of length 20, window size 5, step size 5
	result = sliding_windows(seq_length=20, window_size=5, step_size=5)
	expected = [(0, 5), (5, 10), (10, 15), (15, 20)]
	assert result == expected

def test_sliding_windows_overlap():
	# Overlapping windows
	result = sliding_windows(seq_length=10, window_size=4, step_size=2)
	expected = [(0, 4), (2, 6), (4, 8), (6, 10)]
	assert result == expected

def test_sliding_windows_no_windows():
	# Window size larger than sequence
	result = sliding_windows(seq_length=5, window_size=10, step_size=1)
	expected = []
	assert result == expected

def test_sliding_windows_exact_fit():
	# One window exactly fits the sequence
	result = sliding_windows(seq_length=6, window_size=6, step_size=2)
	expected = [(0, 6)]
	assert result == expected


	
