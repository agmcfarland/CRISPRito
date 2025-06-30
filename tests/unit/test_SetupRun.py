import pytest
import os
from os.path import join as pjoin
from unittest import mock
import pathlib
from CRISPRito.SetupRun import setup_run
import pandas as pd

real_read_csv = pd.read_csv

@mock.patch('CRISPRito.SampleManager.pd.core.frame.DataFrame.to_csv')
@mock.patch('CRISPRito.SetupRun.RunParameters.manage_output_dir')
@mock.patch('CRISPRito.SampleManager.pd.read_csv')
def test_setup_run(
	mock_read_csv,
	mock_manage_dir,
	mock_to_csv,
	project_test_data_directory,
	path_to_hg38_genome,
	path_to_feature_table
	):
	"""
	pytest -sv tests/unit/test_SetupRun.py::test_setup_run
	"""

	mock_read_csv.return_value = real_read_csv(pjoin(project_test_data_directory, 'input_samplesheet_ptprc_reduced.csv'))

	setup_run(
		sample_sheet_path = pjoin(project_test_data_directory, 'input_samplesheet_ptprc_reduced.csv'),
		output_dir = 'dummy_output',
		genome_path = path_to_hg38_genome,
		feature_table_path = path_to_feature_table,
		overwrite_output_dir = True)

	mock_manage_dir.assert_called_once()

	# mock_to_csv.

	print(mock_to_csv.call_args_list)

	assert mock_to_csv.call_count == 2

