import pytest
import pandas as pd
import numpy as np
import pyranges as pr
from CRISPRito.Utils import (
	batch_overlaps,
	batch_nearest_feature
	)

@pytest.fixture
def sample_df():
	return pr.PyRanges(pd.DataFrame({
		'Chromosome': ['chr1', 'chr1', 'chr1', 'chr2', 'chr1', 'chr1'],
		'Start': [100, 200, 300, 100, 1226, 100000],
		'End':   [150, 250, 400, 200, 9800, 100010],
		'label': ['a', 'b', 'c', 'd', 'y', 'z']
	}))

def test_batch_overlaps():
	"""
	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_batch_overlaps
	"""
	#  name    start      end chrom strand name2 feature
	# 1242836  YP_003024037.1_87376  14148.0  14673.0  chrM      -   ND6    exon
	# 1242837  YP_003024037.1_87376  14148.0  14148.0  chrM      -   ND6    3UTR
	# 1242838  YP_003024037.1_87376  14673.0  14673.0  chrM      -   ND6    5UTR
	expanded_df = pr.from_dict({'Chromosome': ['chr1', 'chr1', 'chr1','chr1'],
	'Start': [100, 100, 100, 200],
	'End': [300, 301, 110, 300],
	'label': ['exon', 'exon2', '3UTR', '5UTR']})

	query_sites = pr.from_dict({
		"Chromosome": ['chr1', 'chr1', 'chr1'],
		"Start": 	[105, 115, 305],
		"End": 		[105, 115, 305],
		'label': ['A1', 'A2', 'A3']
	})

	
	assert batch_overlaps(expanded_df, query_sites).shape == (5,7)


def test_batch_nearest_feature():
	"""
	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_batch_nearest_feature
	"""
	#  name    start      end chrom strand name2 feature
	# 1242836  YP_003024037.1_87376  14148.0  14673.0  chrM      -   ND6    exon
	# 1242837  YP_003024037.1_87376  14148.0  14148.0  chrM      -   ND6    3UTR
	# 1242838  YP_003024037.1_87376  14673.0  14673.0  chrM      -   ND6    5UTR
	expanded_df = pr.from_dict({'Chromosome': ['chr1', 'chr1', 'chr1','chr1'],
	'Start': [100, 100, 100, 200],
	'End': [300, 301, 110, 300],
	'label': ['exon', 'exon2', '3UTR', '5UTR']})

	query_sites = pr.from_dict({
		"Chromosome": ['chr1', 'chr1', 'chr1'],
		"Start": 	[105, 115, 305],
		"End": 		[105, 115, 305],
		'label': ['A1', 'A2', 'A3']
	})

	# print(batch_nearest_feature(expanded_df, query_sites))
	print(batch_nearest_feature(expanded_df, query_sites).shape)
	assert batch_nearest_feature(expanded_df, query_sites).shape == (3,8)



# def test_slice_annotation(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_slice_annotation
# 	"""
# 	result = slice_annotation(sample_df, 'chr1', 225, tolerance=1000)
# 	assert len(result) == 3


# def test_extract_annotations_match(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_extract_annotations_match
# 	"""
# 	sample_df = sample_df[sample_df.Chromosome == 'chr1']
# 	result = extract_annotations(sample_df, 'chr1', 225)
# 	print(result)
# 	# expected_labels = ['b']
# 	# assert result['label'].tolist() == expected_labels


# def test_extract_annotations_no_match(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_extract_annotations_no_match
# 	"""
# 	sample_df = sample_df[sample_df.Chromosome == 'chr1']
# 	result = extract_annotations(sample_df, 'chr1', 500)
# 	assert result.empty


# def test_extract_annotations_multiple_matches(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_extract_annotations_multiple_matches
# 	"""
# 	sample_df = sample_df[sample_df.Chromosome == 'chr1']
# 	result = extract_annotations(sample_df, 'chr1', 125)
# 	expected_labels = ['a']
# 	assert result['label'].tolist() == expected_labels


# def test_get_closest_annotation_inside_region(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_get_closest_annotation_inside_region
# 	"""
# 	sample_df = sample_df[sample_df.Chromosome == 'chr1']
# 	label, distance = get_closest_annotation(sample_df, 'chr1', 225, 'label')
# 	assert label == 'b'
# 	assert distance == 0


# def test_get_closest_annotation_outside_region(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_get_closest_annotation_outside_region
# 	"""
# 	sample_df = sample_df[sample_df.Chromosome == 'chr1']
# 	label, distance = get_closest_annotation(sample_df, 'chr1', 260, 'label')
# 	assert label == 'b'
# 	assert distance == 11  # 300 - 260


# def test_get_closest_annotation_ties(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_get_closest_annotation_ties
# 	"""
# 	expanded_df = pr.PyRanges(pd.DataFrame({
# 		'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr1', 'chr1'],
# 		'Start': [100, 200, 300, 100, 1226, 100000, 420,],
# 		'End':   [150, 250, 400, 200, 9800, 100010, 430,],
# 		'label': ['a', 'b', 'c', 'd', 'y', 'z', 'e']
# 	}))
# 	expanded_df = expanded_df[expanded_df.Chromosome == 'chr1']
	
# 	label, distance = get_closest_annotation(expanded_df, 'chr1', 410, 'label')

# 	assert distance == 11
# 	assert label in ['c', 'e']  # either is acceptable due to tie

# def test_batch_annotation(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_batch_annotation
# 	"""
# 	expanded_df = pr.PyRanges(pd.DataFrame({
# 		'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr1', 'chr1'],
# 		'Start': [100, 100, 300, 100, 1226, 100000, 420,],
# 		'End':   [150, 250, 400, 200, 9800, 100010, 430,],
# 		'label': ['a', 'a', 'c', 'd', 'y', 'z', 'e']
# 	}))

# 	query_sites = pr.from_dict({
# 		"Chromosome": ['chr1', 'chr1', 'chr1', 'chr1'],
# 		"Start": [99, 105, 260, 295],
# 		"End": [99, 105, 260, 295]
# 	})

# 	print(batch_annotation(expanded_df, query_sites))

# def test_batch_nearest(sample_df):
# 	"""
# 	pytest -sv tests/unit/test_various_annotations_pyranges.py::test_batch_nearest
# 	"""
# 	expanded_df = pr.PyRanges(pd.DataFrame({
# 		'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr1', 'chr1'],
# 		'Start': [100, 100, 101, 300, 100, 1226, 100000, 420,],
# 		'End':   [150, 250, 110, 400, 200, 9800, 100010, 430,],
# 		'label': ['a', 'a', 'b', 'c', 'd', 'y', 'z', 'e']
# 	}))

# 	query_sites = pr.from_dict({
# 		"Chromosome": ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
# 		"Start": 	[99, 105, 260, 295, 105, 2_000_000],
# 		"End": 		[99, 105, 260, 295, 105, 2_000_000],
# 		'label': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']
# 	})

# 	print('\n')
# 	print(batch_nearest(expanded_df, query_sites))
# 	# print('\n')
# 	# print(batch_annotation(expanded_df, query_sites))





