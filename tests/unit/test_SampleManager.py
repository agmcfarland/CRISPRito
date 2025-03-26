from CRISPRito.SampleManager import SampleManager
from os.path import join as pjoin
import pytest
import os
import pandas as pd

# samples names need to be unique
# at least one cluster group assigned
# measurement must be one of the accepted types
# standard format file path must exist

@pytest.fixture
def load_default_sample_manager_test_input(project_test_data_directory):
	samples = SampleManager(
		input_file = pjoin(project_test_data_directory, 'input_samplesheet_ptprc_reduced.csv')
		)
	samples.load_samplesheet()

	samples.table['standard_format_file_path'] = samples.table['standard_format_file_path'].apply(
		lambda x: pjoin(project_test_data_directory, x)
		)

	return samples


def test_initialize_samplemanager(load_default_sample_manager_test_input):

	samples = load_default_sample_manager_test_input

	for i in samples.table['standard_format_file_path']:
		assert os.path.exists(i)


def test_uuid_generation(load_default_sample_manager_test_input):
	"""
	pytest -sv tests/unit/test_SampleManager.py::test_uuid_generation
	"""

	samples = load_default_sample_manager_test_input

	samples.assign_unique_id()

	print(samples.table)


def test_enforce_measurement_type(load_default_sample_manager_test_input):
	"""
	pytest -sv tests/unit/test_SampleManager.py::test_enforce_measurement_type
	"""

	samples = load_default_sample_manager_test_input

	samples.enforce_measurement_type()

	with pytest.raises(ValueError, match = 'not a valid measurment type.'):
		samples.table['measurement_type'] = 'wrong'
		samples.enforce_measurement_type()






# input_crista_ptprc_reduced.csv
# input_crisprme_ptprc_reduced.csv
# input_iguide_gtsp6619_ptprc_reduced.csv
# input_iguide_gtsp6616_ptprc_reduced.csv

