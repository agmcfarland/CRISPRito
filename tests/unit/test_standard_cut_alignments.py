from unittest.mock import patch
import pytest
import time
import pandas as pd
import pyranges as pr
from os.path import join as pjoin
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from CRISPRito.StandardCuts import *
from CRISPRito.Utils import convert_df_to_granges


@pytest.fixture
def sgRNA_info():
		sgRNA = 'AAAATATGCAAACATCACTG'
		PAM_alignment = '-GG'

		return [{
			'fwd_NGG' : DNA(sgRNA + PAM_alignment),
			'fwd' : DNA(sgRNA)
		},

		{
			'fwd' : range(18,21),
			'fwd_NGG' : range(20,25)
				},

		{
					'fwd' : 0,
					'fwd_NGG' : 3
				}]


def test_alignment1(sgRNA_info):
	"""
	pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment1
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [1375]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(17697666, 17697726)):
		mock_cutsite = CutSite(
			chromosome='chr16',
			strand='+',
			ref_position=17697696,
			cut_region='ATCATCCCTCGCTGATAACCACTAGCCCAAATGCAAACATCACTGAGGCACGAGGGCCTT',
			sgRNA=sgRNA_info[0],
			sgRNA_alignment_tolerance=sgRNA_info[1],
			sgRNA_alignment_start_offset=sgRNA_info[2],
			cut_distance=3,
			detail=mock_detail,
			flank_size=30
		)

	mock_cutsite.find_best_sgRNA_alignment()

	print(mock_cutsite.alignment)

	mock_cutsite.calculate_global_positions()

	print(mock_cutsite.global_position)

	assert mock_cutsite.global_position['protospacer_stop'] == 17697694-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == 17697711-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 17697708-1 #results taken from blat -1 to account for 1-index



def test_alignment2(sgRNA_info):
	"""
	pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment2
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [3607]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(137762432, 137762492)):
		mock_cutsite = CutSite(
			chromosome='chr5',
			strand='+',
			ref_position=137762462,
			cut_region='AATTTGAACCCCATACTTGTCAGAATGGCAAAAATGTCAACATTACTGGGGTTCCTTCCA',
			sgRNA=sgRNA_info[0],
			sgRNA_alignment_tolerance=sgRNA_info[1],
			sgRNA_alignment_start_offset=sgRNA_info[2],
			cut_distance=3,
			detail=mock_detail,
			flank_size=30
		)

	mock_cutsite.find_best_sgRNA_alignment()

	print(mock_cutsite.alignment)

	mock_cutsite.calculate_global_positions()

	print(mock_cutsite.global_position)

	assert mock_cutsite.global_position['protospacer_stop'] == 137_762_462-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == 137_762_480-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 137_762_477-1 #results taken from blat -1 to account for 1-index



def test_alignment3(sgRNA_info):
	"""
	pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment3
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [3607]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(137762432, 137762492)):
		mock_cutsite = CutSite(
			chromosome='chr5',
			strand='+',
			ref_position=137762462,
			cut_region='AATTTGAACCCCATACTTGTCAGAATGGCAAAAATGTCAACATTACTGGGGTTCCTTCCA',
			sgRNA=sgRNA_info[0],
			sgRNA_alignment_tolerance=sgRNA_info[1],
			sgRNA_alignment_start_offset=sgRNA_info[2],
			cut_distance=3,
			detail=mock_detail,
			flank_size=30
		)

	mock_cutsite.find_best_sgRNA_alignment()

	print(mock_cutsite.alignment)

	mock_cutsite.calculate_global_positions()

	print(mock_cutsite.global_position)

	assert mock_cutsite.global_position['protospacer_stop'] == 137_762_462-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == 137_762_480-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 137_762_477-1 #results taken from blat -1 to account for 1-index







