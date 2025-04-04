from unittest.mock import patch
import pytest
import pandas as pd
from os.path import join as pjoin
from CRISPRito.StandardCuts import StandardCuts

@pytest.fixture
def load_1_group_samplesheet_ptprc(project_test_data_directory):

	df_gr1 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '1_group_samplesheet_ptprc.csv'))

	df_gr1['standard_format_file_path'] = df_gr1['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr1

@pytest.fixture
def load_2_group_samplesheet_ptprc(project_test_data_directory):

	df_gr2 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '2_group_samplesheet_ptprc.csv'))

	df_gr2['standard_format_file_path'] = df_gr2['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr2


def test_load_1_group_samplesheet_ptprc(load_1_group_samplesheet_ptprc):

	df_gr1 = load_1_group_samplesheet_ptprc

	assert True

def test_load_cut_sites_group1(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_load_cut_sites_group1
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	# print('\n')

	# print(standard_group.df_cut_sites)

	standard_group.update_cut_cluster_id()

	# print(standard_group.df_cut_sites)

	assert len(standard_group.df_cut_sites) == 19


def test_load_cut_sites_group2(load_2_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_load_cut_sites_group2
	"""

	standard_group = StandardCuts(sample_sheet = load_2_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	# print('\n')

	# print(standard_group.df_cut_sites)

	standard_group.update_cut_cluster_id()

	# print(standard_group.df_cut_sites)

	assert len(standard_group.df_cut_sites) == 5


def test_standarize_scores(load_1_group_samplesheet_ptprc):
	"""

	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	# print('\n')

	# print(standard_group.df_cut_sites)

	standard_group.update_cut_cluster_id()

	

def skip_test_extract_cut_region_group1(load_1_group_samplesheet_ptprc, path_to_hg38_genome):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_extract_cut_region_group1
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.load_genome(genome_path = path_to_hg38_genome)

	standard_group.extract_cut_region()

	# print(standard_group.df_reference_cut_sites)

	assert len(standard_group.df_reference_cut_sites) == 15

	# standard_group.df_reference_cut_sites.to_csv('/data/CRISPRito/1_group_cluster_regions_ptprc.csv', index = None)

def skip_test_extract_cut_region_group2(load_2_group_samplesheet_ptprc, path_to_hg38_genome):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_extract_cut_region_group2
	"""

	standard_group = StandardCuts(sample_sheet = load_2_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.load_genome(genome_path = path_to_hg38_genome)

	standard_group.extract_cut_region()

	# print(standard_group.df_reference_cut_sites)

	# standard_group.df_reference_cut_sites.to_csv('/data/CRISPRito/2_group_cluster_regions_ptprc.csv', index = None)

	assert len(standard_group.df_reference_cut_sites) == 5



def test_build_cut_sites(load_1_group_samplesheet_ptprc, project_test_data_directory):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_build_cut_sites
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	# standard_group.load_genome(genome_path = path_to_hg38_genome)

	# standard_group.extract_cut_region()
	with patch.object(standard_group, 'extract_cut_region', return_value = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites', '1_group_cluster_regions_ptprc.csv'))):
		standard_group.df_reference_cut_sites = standard_group.extract_cut_region()	

	print(standard_group.df_reference_cut_sites)

	print(standard_group.df_cut_sites)

	# standard_group.df_reference_cut_sites = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites', '1_group_cluster_regions_ptprc.csv'))

	# print(standard_group.df_reference_cut_sites)

	# print(standard_group.df_reference_cut_sites)






















