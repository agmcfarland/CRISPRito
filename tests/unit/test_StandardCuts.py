import pytest
import pandas as pd
from os.path import join as pjoin

@pytest.fixture
def load_1_group_samplesheet_ptprc(project_test_data_directory):

	df_gr1 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '1_group_samplesheet_ptprc.csv'))

	df_gr1['standard_format_file_path'] = df_gr1['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr1

def test_it_worked(load_1_group_samplesheet_ptprc):

	df_gr1 = load_1_group_samplesheet_ptprc

	print(df_gr1)