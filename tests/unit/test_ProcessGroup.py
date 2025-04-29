import pytest
import pandas as pd
from unittest import mock
from os.path import join as pjoin
from CRISPRito.ProcessGroup import process_samples

@pytest.fixture
def load_1_group_samplesheet_ptprc(project_test_data_directory):

	df_gr1 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '1_group_samplesheet_ptprc.csv'))

	df_gr1['standard_format_file_path'] = df_gr1['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr1

real_read_csv = pd.read_csv

def side_effect_read_csv_factory(preloaded_df):
	def _side_effect_read_csv(path, *args, **kwargs):
		if path == "group1_fake_samplesheet.csv":
			return preloaded_df
		else:
			return real_read_csv(path, *args, **kwargs)
			
	return _side_effect_read_csv

@mock.patch('CRISPRito.ProcessGroup.pd.core.frame.DataFrame.to_csv')
@mock.patch('CRISPRito.ProcessGroup.pd.read_csv')
def test_process_samples_1(
	mock_read_csv,
	mock_to_csv,
	project_test_data_directory,
	load_1_group_samplesheet_ptprc,
	path_to_hg38_genome,
	path_to_hg38_refseq,
	path_to_hg38_gene_names):

	"""
	pytest -sv tests/unit/test_ProcessGroup.py::test_process_samples_1
	"""
	output_dir = '/fake/path'
	input_samplesheet_path = "group1_fake_samplesheet.csv"
	mock_read_csv.side_effect = side_effect_read_csv_factory(load_1_group_samplesheet_ptprc)

	process_samples(
		group_samplesheet_path=input_samplesheet_path,
		output_path=output_dir,
		genome_path= path_to_hg38_genome,
		feature_path= path_to_hg38_refseq,
		gene_names_path= path_to_hg38_gene_names,
		flank_size=30,
		sgRNA='AAAATATGCAAACATCACTG',
		PAM_alignment='-GG'
	)

	mock_to_csv.call_count == 3

	filepaths = [call.args[0] for call in mock_to_csv.call_args_list]

	assert any('1_group_cut_profiles' in path for path in filepaths)
	assert any('1_group_method_counts' in path for path in filepaths)
	assert any('1_group_id_counts' in path for path in filepaths)
