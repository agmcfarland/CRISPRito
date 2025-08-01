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
	Helper function
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
	for i in sc.ranking_criteria:
		assert len(i) == 8 # check number of values used for ranking

def test_load_user_sample_rank_list(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_sample_rank_list
	"""
	sc = StandardCutRank()
	sc.load_user_sample_rank_list(load_1_group_samplesheet_ptprc)
	print(sc.ranking_criteria)
	assert len(sc.ranking_criteria) == 3
	for i in sc.ranking_criteria:
		assert len(i) == 8

def test_load_user_feature_rank_list(mock_registry):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_feature_rank_list
	"""
	feature_manager = generate_feature_table(mock_registry)
	
	sc = StandardCutRank()
	sc.load_user_feature_rank_list(feature_manager)
	assert len(sc.ranking_criteria) == 3
	for i in sc.ranking_criteria:
		assert len(i) == 8

def test_load_user_sample_rank_list(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_sample_rank_list
	"""
	sc = StandardCutRank()
	sc.load_user_sample_rank_list(load_1_group_samplesheet_ptprc)
	print(sc.ranking_criteria)
	assert len(sc.ranking_criteria) == 3
	for i in sc.ranking_criteria:
		assert len(i) == 8

def test_load_user_method_rank_list(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_load_user_method_rank_list
	This one is a bit weird because only iGUIDE is a true method and the rest are just
	the sampel names instead of something like computational which would be better.
	"""
	sc = StandardCutRank()
	sc.load_user_method_rank_list(load_1_group_samplesheet_ptprc)
	assert len(sc.ranking_criteria) == 6
	for i in sc.ranking_criteria:
		assert len(i) == 8


def test_generate_ranking_skeleton(load_1_group_samplesheet_ptprc, mock_registry):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_generate_ranking_skeleton
	"""
	feature_manager = generate_feature_table(mock_registry)
	sc = StandardCutRank()
	sc.generate_ranking_skeleton(sample_sheet = load_1_group_samplesheet_ptprc, feature_manager = feature_manager)
	print(sc.ranking_criteria)
	print('\n')
	print(sc.df_skeleton)
	assert len(sc.ranking_criteria) == 29
	assert sc.df_skeleton.shape == (29, 8)


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

	Testing that the score_criteria() method for feature works
	Testing that the upper and lower limit works
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	sc.datasets['cut_profiles']['zscore_mean'] = [100, 90, 80, 10, 1, 2, 3, 4, -0.5, -1]

	# print(sc.datasets['cut_profiles']['zscore_mean'])
	# return
	row_id = 2

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			break
	print(row)
	# print(row)
	rank_criteria = RankOperator.load_from_rank_table_row(row)
	rank_criteria.score_criteria(datasets = sc.datasets)
	print(rank_criteria.score)
	# return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [2, 2, 2, 2, 0, 2, 2, 2, 0, 0]
	# return
	print(row)
	row.upper_threshold = 10
	rank_criteria = RankOperator.load_from_rank_table_row(row)
	rank_criteria.score_criteria(datasets = sc.datasets)
	print(rank_criteria.score)
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [0, 0, 0, 2, 0, 2, 2, 2, 0, 0]



def test_temp_2(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_2
	test overlap on the df_cut_sites df works
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)

	# print(sc.rank_table_weights)
	# return
	row_id = 0

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			break
	rank_criteria = RankOperator.load_from_rank_table_row(row)
	rank_criteria.score_criteria(datasets = sc.datasets)
	print(rank_criteria.score)
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]

	print(row)
	row2 = row
	row2.upper_threshold = 90.0
	print(row2)

	rank_criteria2 = RankOperator.load_from_rank_table_row(row2)
	# print(rank_criteria.__dict__)
	rank_criteria2.score_criteria(datasets = sc.datasets)
	print(rank_criteria2.score)
	assert rank_criteria2.score[f'rank_criteria_{row_id}'].tolist() == [0, 0, 0, 1, 1, 1, 0, 0, 0, 0]


def test_temp_3(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_3
	alignment_length >= 17 and 22
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 16

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			break
	print(row)
	# return
	rank_criteria = RankOperator.load_from_rank_table_row(row)
	rank_criteria.score_criteria(datasets = sc.datasets)
	print(rank_criteria.score)
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [2, 2, 2, 2, 2, 2, 2, 2, 0, 2]



def test_temp_4(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_4
	nearest_gene_distance 1-1000
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 28

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 28:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
			# return
		pass
	# return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [0, 1, 1, 0, 0 ,0 ,0 ,0, 1 ,1]

def test_temp_5(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_5
	Overlap of method (biochemical) with percentage
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 24

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			# return
			break
		pass
	# return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [1, 1, 1, 1, 1, 1, 0, 0, 1, 0]

def test_temp_5_presence(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_5_presence
	Overlap of method (biochemical) with raw number
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 26

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			# return
			break
		pass
	return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [0,0,0,0,0,-1,0,0,-1, 0]


def test_temp_6(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_6
	Testing biosample4
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 20

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			# print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

def test_temp_7(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_7
	nearest oncogene distance == 0,0
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 31

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == row_id:
			print(row)
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass
	# return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [2, 2, 0, 0, 0, 2, 0, 0, 0, 0]

def test_temp_8(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_8
	overlap method biochemical percentage
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 23

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 23:
			print(row)
			# return
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass
	# return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [2, 2, 2, 2, 2, 2, 2, 0, 0, 0]

def test_temp_9(ranking_inputs):
	"""
	pytest -sv tests/unit/test_StandardCutRank.py::test_temp_9
	overlap presence
	"""

	sc = StandardCutRank(
		rank_table_weights=pd.read_csv(ranking_inputs['weight_skeleton']),
		cut_profiles=pd.read_csv(ranking_inputs['cut_profiles']),
		id_counts=pd.read_csv(ranking_inputs['id_counts']),
		method_counts=pd.read_csv(ranking_inputs['method_counts']),
		samplesheet=pd.read_csv(ranking_inputs['samplesheet'])
		)

	# print(sc.datasets)
	row_id = 1

	# print(df)
	for e, row in sc.rank_table_weights.iterrows():
		if e == 1:
			print(row)
			# return
			# break
			rank_criteria = RankOperator.load_from_rank_table_row(row)
			print(rank_criteria.score)
			rank_criteria.score_criteria(datasets = sc.datasets)
			print(rank_criteria.score)
			break
		pass
	# return
	assert rank_criteria.score[f'rank_criteria_{row_id}'].tolist() == [0,0,0,0,0,0,0,1,1,1]


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
	pytest -sv tests/unit/test_StandardCutRank.py::test_tally_cut_cluster_scores
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

	print(sc.cut_cluster_scores['total_score'])
	print(sc.cut_cluster_scores.shape)
	# return

	assert sc.cut_cluster_scores.shape == (10, 35)

	# print(sc.cut_cluster_scores[['cut_cluster','rank_criteria_0']])

	# return

	assert sc.cut_cluster_scores['total_score'].tolist() == [20.5,17.5,16.5,14.0,14.0,17.5,12.5,13.5, 3.0, 4.0]























