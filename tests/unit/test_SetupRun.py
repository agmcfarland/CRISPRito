import pytest
import os
from os.path import join as pjoin
from unittest import mock
from CRISPRito.SetupRun import setup_run
import pandas as pd

real_read_csv = pd.read_csv

@mock.patch('CRISPRito.SetupRun.RunParameters.manage_output_dir')
@mock.patch('CRISPRito.SampleManager.pd.read_csv')
def test_setup_run(
	mock_read_csv,
	mock_manage_dir,
	project_test_data_directory,
	path_to_hg38_genome,
	path_to_hg38_refseq,
	path_to_hg38_gene_names
	):
	"""
	pytest -sv tests/unit/test_SetupRun.py::test_setup_run
	"""

	mock_read_csv.return_value = real_read_csv(pjoin(project_test_data_directory, 'input_samplesheet_ptprc_reduced.csv'))

	# Act
	setup_run(
		input_file = 'dummy_input.csv',
		output_dir = 'dummy_output',
		
		overwrite_output_dir = True)

	# Assert
	mock_manage_dir.assert_called_once()


