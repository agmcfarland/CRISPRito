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
	assert sc.df_skeleton.shape == (26, 7)


def test_input_dataframes_to_standardcutrank(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_input_dataframes_to_standardcutrank
	"""
	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)
	print(sc.id_counts)
	print('\nnext')
	print(sc.method_counts)
	# print(sc.samplesheet)
	# print(sc.rank_table_weights)

def test_load_rank_from_table_row(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_rank_from_table_row
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(df)
	for _, row in sc.rank_table_weights.iterrows():
		rank_criteria = RankOperator.load_from_rank_table_row(row)
		print(rank_criteria.__dict__)
		assert True
		break
	pass



def test_temp_1(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_1
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 1:
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.__dict__)
			print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass


def test_temp_2(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_2
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 0:
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass

def test_temp_3(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_3
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 17:
			# print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass


def test_temp_4(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_4
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 25:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass

def test_temp_5(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_5
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 23:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass


def test_temp_6(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_6
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 24:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass

def test_temp_7(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_7
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 30:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass

def test_temp_go_through_all(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_go_through_all
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	for e, row in sc.rank_table_weights.iterrows():
		rank_criteria = RankOperator.load_from_rank_table_row(row)
		rank_criteria.score_criteria(datasets = sc.datasets)

def test_get_score_from_rank_table(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_get_score_from_rank_table
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	sc.get_score_from_rank_table()


def test_tally_cut_cluster_scores(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::tally_cut_cluster_scores
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	sc.get_score_from_rank_table()

	sc.tally_cut_cluster_scores()

	assert sc.cut_cluster_scores.shape == (10, 33)
























