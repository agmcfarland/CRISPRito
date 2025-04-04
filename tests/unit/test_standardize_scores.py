from Utils import scale_min_max, scale_zscore
import pytest
import numpy as np



@pytest.fixture
def example_scores():
	data_groups = [
		np.array([0, 10, 20, 30, 34, 35, 40]),  # Max=40, Min=0
		np.array([0, 0.5, .6, .7, .8, 1]),          # Max=1, Min=0
		np.array([0, 0.5, .6, .7, .8, 0.89]),          # Max=0.89, Min=0
		np.array([40, 20, 0]),          # Max=40, Min=0 (but reversed order)
		np.array([1, 500, 600, 700, 800, 1000]),        # Max=1000, Min=1
		np.array([1, 3, 3, 3, 4, 500, 1000]),        # Max=1000, Min=1
		np.array([0,0,0,0]),        # Max=1000, Min=1
		np.array([1,1,1,1])        # Max=1000, Min=1
	]
	return data_groups

@pytest.fixture
def expected_scores_min_max():
	return [[0.0, 0.25, 0.5, 0.75, 0.8500000000000001, 0.875, 1.0],
	[0.0, 0.5, 0.6, 0.7, 0.8, 1.0],
	[0.0, 0.5617977528089888, 0.6741573033707865, 0.7865168539325843, 0.8988764044943821, 1.0],
	[1.0, 0.5, 0.0],
	[0.0, 0.49949949949949957, 0.5995995995995996, 0.6996996996996997, 0.7997997997997998, 1.0],
	[0.0, 0.002002002002002002, 0.002002002002002002, 0.002002002002002002, 0.003003003003003003, 0.49949949949949957, 1.0],
	[0.0, 0.0, 0.0, 0.0],
	[0.0, 0.0, 0.0, 0.0]]


@pytest.fixture
def expected_scores_zscore():
	return [[-1.7710176183055835, -1.0374600249245725, -0.3039024315435616, 0.4296551618374493, 0.7230781991898536, 0.7964339585279547, 1.1632127552184601],
	[-1.9298025627080304, -0.32163376045133835, 0.0, 0.32163376045133835, 0.643267520902677, 1.2865350418053538],
	[-2.0107676946857143, -0.28231408893868193, 0.06337663221072445, 0.40906735336013084, 0.7547580745095376, 1.0658797235440034],
	[1.224744871391589, 0.0, -1.224744871391589],
	[-1.9291175662330593, -0.3225033817263055, -0.00053661128407027, 0.321430159158165, 0.6433969296004002, 1.2873304704848707],
	[-0.5930078822799147, -0.5874988508586017, -0.5874988508586017, -0.5874988508586017, -0.5847443351479451, 0.781495457337698, 2.1587533126659673],
	[np.nan, np.nan, np.nan, np.nan],
	[np.nan, np.nan, np.nan, np.nan]]


def test_scale_to_min_max(example_scores, expected_scores_min_max):
	"""
	pytest -sv tests/unit/test_standardize_scores.py::test_scale_to_min_max
	"""
	print('\n')
	for e, i in enumerate(example_scores):
		# print('unscaled:', i)
		# print('scaled  :', scale_min_max(i))
		result = scale_min_max(i)
		print(i)
		assert result == expected_scores_min_max[e]
		# print(expected_scores_min_max[e])


def test_scale_to_min_max_list(example_scores, expected_scores_min_max):
	"""
	pytest -sv tests/unit/test_standardize_scores.py::test_scale_to_min_max_list
	"""
	print('\n')
	for e, i in enumerate(example_scores):
		i = list(i)
		# print('unscaled:', i)
		# print('scaled  :', scale_min_max(i))
		result = scale_min_max(i)
		assert result == expected_scores_min_max[e]


def test_z_standardize(example_scores, expected_scores_zscore):
	"""
	pytest -sv tests/unit/test_standardize_scores.py::test_z_standardize
	"""
	for e, i in enumerate(example_scores[:6]):
		# print('unscaled:', i)
		# print('scaled  :', scale_zscore(i))
		result = scale_zscore(i)
		assert result == expected_scores_zscore[e]


	for e, i in enumerate(example_scores[-2:]):
		result = scale_zscore(i)
		for e1, z in enumerate(result):
			assert np.isnan(z) == np.isnan(expected_scores_zscore[-2:][e][e1])

def test_z_standardize_list(example_scores, expected_scores_zscore):
	"""
	pytest -sv tests/unit/test_standardize_scores.py::test_z_standardize_list
	"""
	for e, i in enumerate(example_scores[:6]):
		i = list(i)
		# print('unscaled:', i)
		# print('scaled  :', scale_zscore(i))
		result = scale_zscore(i)
		assert result == expected_scores_zscore[e]


	for e, i in enumerate(example_scores[-2:]):
		i = list(i)
		result = scale_zscore(i)
		for e1, z in enumerate(result):
			assert np.isnan(z) == np.isnan(expected_scores_zscore[-2:][e][e1])





