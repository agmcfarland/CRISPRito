import pytest
import pandas as pd
from unittest import mock
from os.path import join as pjoin
from CRISPRito.ProcessGroup import process_group

@pytest.fixture
def load_1_group_samplesheet_ptprc(project_test_data_directory):

	df_gr1 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '1_group_samplesheet_ptprc.csv'))

	df_gr1['standard_format_file_path'] = df_gr1['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr1

def side_effect_read_csv_factory(path_to_df_map):
	
	real_read_csv = __import__('pandas').read_csv # Unmocked read_csv straight from pandas module

	def _side_effect_read_csv(path, *args, **kwargs):
		if path in path_to_df_map:
			return path_to_df_map[path]
		return real_read_csv(path, *args, **kwargs)

	return _side_effect_read_csv

@pytest.mark.usefixtures("load_1_group_samplesheet_ptprc")
def test_process_group_1(
	project_test_data_directory,
	load_1_group_samplesheet_ptprc,
	retrieve_feature_input,
	path_to_hg38_genome
):
	"""
	Testing strategy is to mock some of the read_csv inputs so that the proper file_paths within them can be constructed.
	"""
	output_dir = '/fake/path'
	input_samplesheet_path = "group1_fake_samplesheet.csv"
	input_featuresheet_path = 'fake_featuresheet.csv'

	preloaded_csvs = {
		input_samplesheet_path: load_1_group_samplesheet_ptprc,
		input_featuresheet_path: retrieve_feature_input
	}

	side_effect = side_effect_read_csv_factory(preloaded_csvs)

	with mock.patch('CRISPRito.ProcessGroup.pd.read_csv', side_effect=side_effect), \
	     mock.patch('CRISPRito.FeatureManager.pd.read_csv', side_effect=side_effect), \
	     mock.patch('CRISPRito.ProcessGroup.pd.core.frame.DataFrame.to_csv') as mock_to_csv:

		process_group(
			group_samplesheet_path=input_samplesheet_path,
			output_path=output_dir,
			genome_path=path_to_hg38_genome,
			flank_size=30,
			sgRNA='AAAATATGCAAACATCACTG',
			PAM_alignment='-GG',
			feature_table_path=input_featuresheet_path
		)

		assert mock_to_csv.call_count == 5
		filepaths = [call.args[0] for call in mock_to_csv.call_args_list]
		assert any('1_group_cut_profiles' in path for path in filepaths)
		assert any('1_group_method_counts' in path for path in filepaths)
		assert any('1_group_id_counts' in path for path in filepaths)
		assert any('1_group_id_rank_weight_skeleton' in path for path in filepaths)

