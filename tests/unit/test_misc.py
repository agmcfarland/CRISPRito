import pytest
import pandas as pd
from CRISPRito.Utils import df_long_to_wide

def test_df_long_to_wide():
	"""
	written by chatgpt
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