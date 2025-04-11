import pytest
import pandas as pd
import numpy as np
from CRISPRito.Utils import extract_annotations, get_closest_annotation

@pytest.fixture
def sample_df():
	return pd.DataFrame({
		'chrom': ['chr1', 'chr1', 'chr1', 'chr2'],
		'start': [100, 200, 300, 100],
		'end':   [150, 250, 400, 200],
		'label': ['a', 'b', 'c', 'd']
	})


def test_extract_annotations_match(sample_df):
	result = extract_annotations(sample_df, 'chr1', 225)
	expected_labels = ['b']
	assert result['label'].tolist() == expected_labels


def test_extract_annotations_no_match(sample_df):
	result = extract_annotations(sample_df, 'chr1', 500)
	assert result.empty


def test_extract_annotations_multiple_matches(sample_df):
	result = extract_annotations(sample_df, 'chr1', 125)
	expected_labels = ['a']
	assert result['label'].tolist() == expected_labels


def test_get_closest_annotation_inside_region(sample_df):
	label, distance = get_closest_annotation(sample_df, 'chr1', 225, 'label')
	assert label == 'b'
	assert distance == 0


def test_get_closest_annotation_outside_region(sample_df):
	label, distance = get_closest_annotation(sample_df, 'chr1', 260, 'label')
	assert label == 'b'
	assert distance == 10  # 300 - 260


def test_get_closest_annotation_ties(sample_df):
	"""
	pytest -sv tests/unit/test_various_annotations.py::test_get_closest_annotation_ties
	"""
	# Add another region at the same distance as existing one
	df = pd.concat([sample_df, pd.DataFrame({'chrom': ['chr1'], 'start': [420], 'end': [430], 'label': ['e']})])
	label, distance = get_closest_annotation(df, 'chr1', 410, 'label')

	assert distance == 10
	assert label in ['c', 'e']  # either is acceptable due to tie


def test_get_closest_annotation_chromosome_not_found(sample_df):
	with pytest.raises(IndexError):
		get_closest_annotation(sample_df, 'chrX', 150, 'label')



# get_closest_annotation(df = dfx, chromosome = 'chr1', position=260, column_name='label')


# df = dfx
# df = pd.concat([dfx, pd.DataFrame({'chrom': ['chr1'], 'start': [420], 'end': [430], 'label': ['e']})])
# chromosome = 'chr1'
# position=410
# column_name='label'

# df = df[df['chrom'] == chromosome].copy()

# df["distance"] = np.where(
# 	(df["start"] <= position) & (df["end"] >= position),
# 	0,
# 	np.minimum(np.abs(df["start"] - position), np.abs(df["end"] - position))
# )

# df = df[df["distance"] == df["distance"].min()]

# return list(df[column_name].unique())[0], df['distance'].tolist()[0]


