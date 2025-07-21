import os
import tempfile
import pytest
import pandas as pd
from unittest import mock
from os.path import join as pjoin
from CRISPRito.RankSites import rank_sites


def test_rank_sites_1(ranking_inputs):
	"""
	pytest -sv tests/unit/test_RankSites.py::test_rank_sites_1
	"""

	with tempfile.TemporaryDirectory() as tmpdir:
		with mock.patch('CRISPRito.RankSites.pd.core.frame.DataFrame.to_csv') as mock_to_csv:

			rank_sites(
				group_samplesheet_path = ranking_inputs['samplesheet'],
				rank_table_weights_path = ranking_inputs['weight_skeleton'],
				cut_profiles_path = ranking_inputs['cut_profiles'],
				id_counts_path = ranking_inputs['id_counts'],
				method_counts_path = ranking_inputs['method_counts'],
				output_dir = tmpdir,
				output_name = 'tempout'
				)

			assert mock_to_csv.call_count == 2

