import pytest
import pandas as pd
import numpy as np
import pyranges as pr
from CRISPRito.Utils import extract_annotations, get_closest_annotation, slice_annotation

@pytest.fixture
def sample_df():
	return pr.PyRanges(pd.DataFrame({
		'Chromosome': ['chr1', 'chr1', 'chr1', 'chr2', 'chr1', 'chr1'],
		'Start': [100, 200, 300, 100, 1226, 100000],
		'End':   [150, 250, 400, 200, 9800, 100010],
		'label': ['a', 'b', 'c', 'd', 'y', 'z']
	}))


def test_slice_annotation(sample_df):
	"""
	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_slice_annotation
	"""
	result = slice_annotation(sample_df, 'chr1', 225, tolerance=1000)
	assert len(result) == 3


def test_extract_annotations_match(sample_df):
	"""
	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_extract_annotations_match
	"""
	sample_df = sample_df[sample_df.Chromosome == 'chr1']
	result = extract_annotations(sample_df, 225)
	expected_labels = ['b']
	assert result['label'].tolist() == expected_labels


def test_extract_annotations_no_match(sample_df):
	sample_df = sample_df[sample_df['chrom'] == 'chr1']
	result = extract_annotations(sample_df, 500)
	assert result.empty


def test_extract_annotations_multiple_matches(sample_df):
	sample_df = sample_df[sample_df['chrom'] == 'chr1']
	result = extract_annotations(sample_df, 125)
	expected_labels = ['a']
	assert result['label'].tolist() == expected_labels


def test_get_closest_annotation_inside_region(sample_df):
	sample_df = sample_df[sample_df['chrom'] == 'chr1']
	label, distance = get_closest_annotation(sample_df, 225, 'label')
	assert label == 'b'
	assert distance == 0


def test_get_closest_annotation_outside_region(sample_df):
	sample_df = sample_df[sample_df['chrom'] == 'chr1']
	label, distance = get_closest_annotation(sample_df, 260, 'label')
	assert label == 'b'
	assert distance == 10  # 300 - 260


def test_get_closest_annotation_ties(sample_df):
	"""
	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_get_closest_annotation_ties
	"""
	sample_df = sample_df[sample_df['chrom'] == 'chr1']
	# Add another region at the same distance as existing one
	df = pd.concat([sample_df, pd.DataFrame({'chrom': ['chr1'], 'start': [420], 'end': [430], 'label': ['e']})])
	label, distance = get_closest_annotation(df, 410, 'label')

	assert distance == 10
	assert label in ['c', 'e']  # either is acceptable due to tie


# def test_get_closest_annotation_chromosome_not_found(sample_df):
# 	with pytest.raises(IndexError):
# 		get_closest_annotation(sample_df, 'chrX', 150, 'label')

