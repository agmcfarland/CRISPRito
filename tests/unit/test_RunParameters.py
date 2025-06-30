import pytest
from CRISPRito.RunParameters import RunParameters
from unittest import mock


@mock.patch('CRISPRito.RunParameters.os.path.exists')
def test_check_inputs_exist(
	mock_os_path_exists
	):
	"""
	pytest -sv tests/unit/test_RunParameters.py::test_check_inputs_exist
	"""
	mock_os_path_exists.return_value = True

	rp = RunParameters(output_dir = '')

	rp.check_inputs_exist(
		sample_sheet_path = 'path1',
		genome_path = 'path2',
		feature_table_path = 'path3',
		)
	exist_calls = [call.args[0] for call in mock_os_path_exists.call_args_list]

	assert set(exist_calls) == {'path1', 'path2', 'path3'}


@mock.patch('CRISPRito.RunParameters.os.path.exists')
def test_check_inputs_donot_exist(mock_os_path_exists):
	"""
	pytest -sv tests/unit/test_RunParameters.py::test_check_inputs_donot_exist
	"""
	mock_os_path_exists.return_value = False

	rp = RunParameters(output_dir='')

	with pytest.raises(ValueError, match="does not exist"):
		rp.check_inputs_exist(
			sample_sheet_path = 'path1',
			genome_path = 'path2',
			feature_table_path = 'path3',
			)
	exist_calls = [call.args[0] for call in mock_os_path_exists.call_args_list]
	
	assert set(exist_calls) == {'path1'}#, 'path2', 'path3', 'path4'}

