import pytest
from unittest import mock
from os.path import join as pjoin
import pandas as pd
from CRISPRito.StandardCutRank import StandardCutRank, RankOperator
from CRISPRito.FeatureManager import FeatureManager

@pytest.fixture
def load_1_group_samplesheet_ptprc(project_test_data_directory):

	df_gr1 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '1_group_samplesheet_ptprc.csv'))

	df_gr1['standard_format_file_path'] = df_gr1['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr1

@pytest.fixture
def mock_registry():
	return {
		'featureA': {
			'type': 'annotation',
			'file_path': '/mock/path/featureA.csv'
		},
		'featureB': {
			'type': 'presence_absence',
			'file_path': '/mock/path/featureB.csv'
		}
	}

def generate_feature_table(mock_registry):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_feature_manager_with_mock_registry
	"""
	mock_df = pd.DataFrame()  # Not used since we're patching

	with mock.patch.object(FeatureManager, '_validate_and_build'):
		fm = FeatureManager(mock_df)
		fm.registry = mock_registry
	return fm


def test_load_default_rank_list():
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_default_rank_list
	"""
	sc = StandardCutRank()
	sc.load_default_rank_list()
	assert len(sc.ranking_criteria) == 17

def test_load_user_sample_rank_list(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_sample_rank_list
	"""
	sc = StandardCutRank()
	sc.load_user_sample_rank_list(load_1_group_samplesheet_ptprc)
	print(sc.ranking_criteria)
	assert len(sc.ranking_criteria) == 3


def test_load_user_feature_rank_list(mock_registry):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_feature_rank_list
	"""
	feature_manager = generate_feature_table(mock_registry)
	
	sc = StandardCutRank()
	sc.load_user_feature_rank_list(feature_manager)
	assert len(sc.ranking_criteria) == 3


def test_load_user_sample_rank_list(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_sample_rank_list
	"""
	sc = StandardCutRank()
	sc.load_user_sample_rank_list(load_1_group_samplesheet_ptprc)
	print(sc.ranking_criteria)
	assert len(sc.ranking_criteria) == 3


def test_load_user_method_rank_list(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_method_rank_list
	"""
	sc = StandardCutRank()
	sc.load_user_method_rank_list(load_1_group_samplesheet_ptprc)
	assert len(sc.ranking_criteria) == 3


def test_generate_ranking_skeleton(load_1_group_samplesheet_ptprc, mock_registry):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_generate_ranking_skeleton
	"""
	feature_manager = generate_feature_table(mock_registry)
	sc = StandardCutRank()
	sc.generate_ranking_skeleton(sample_sheet = load_1_group_samplesheet_ptprc, feature_manager = feature_manager)
	assert len(sc.ranking_criteria) == 26
	assert sc.df_skeleton.shape == (26, 6)


