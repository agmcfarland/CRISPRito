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
	{'chromosome': 'chr5',
	 'strand': '+',
	 'ref_position': 137762462,
	 'cut_region': {'start': 137762432,
	  'stop': 137762492,
	  'sequence': 'AATTTGAACCCCATACTTGTCAGAATGGCAAAAATGTCAACATTACTGGGGTTCCTTCCA'},
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [3120]})

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
	CutSite(chrom=chr5, strand=+, ref_pos=32948416, cut=32948433 diversity=2)
	3033
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [3033]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(32948386, 32948446)):
		mock_cutsite = CutSite(
			chromosome='chr5',
			strand='+',
			ref_position=137762462,
			cut_region='AATAGATAAGCAGCTATTAGCGCTAGTTTAAAATACTTGCAACATCACTGAGGGATTCCT',
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

	assert mock_cutsite.global_position['protospacer_stop'] == 32_948_416-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == 32_948_436-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 32_948_433-1 #results taken from blat -1 to account for 1-index



def test_alignment4(sgRNA_info):
	"""
	pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment4
	1380, # aligned sequence normal, gRNA one gap, plus strand, aln length 21
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [1380]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(30498437, 30498497)):
		mock_cutsite = CutSite(
			chromosome='chr16',
			strand='+',
			ref_position=30498467,
			cut_region='AAGCTGCAACTTTATTATGAAGCATAATAAAATTCATGCAAGCATCACTGAAGGCTTAAA',
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

	assert mock_cutsite.global_position['protospacer_stop'] == 30_498_467-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == 30_498_488-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 30_498_485-1 #results taken from blat -1 to account for 1-index




def test_alignment5(sgRNA_info):
	"""
	pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment5
	1380, # aligned sequence normal, gRNA one gap, plus strand, aln length 21
	923	chr13	+	25001778	0	0	1	1	CENPJ	78918	-0.674640176	-0.674640176	-0.674640176	-0.674640176	0.117857297	0.117857297	0.117857297	0.117857297	2195	2195	2195	2195	25001781	25001762	29	48	ATAAAATGCAATCTTCACTG	AAAATATGCAAACATCACTG	GGA	20	CTCCACTGCTGATCCCCAGCTTAGGCAAAATAAAATGCAATCTTCACTGGGAAATATGCG
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [923]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(25001733, 25001793)):
		mock_cutsite = CutSite(
			chromosome='chr13',
			strand='+',
			ref_position=25001763,
			cut_region='CTCCACTGCTGATCCCCAGCTTAGGCAAAATAAAATGCAATCTTCACTGGGAAATATGCG',
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

	assert mock_cutsite.global_position['protospacer_stop'] == 25_001_763-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == 25_001_781-1 #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 25_001_778-1 #results taken from blat -1 to account for 1-index



def test_alignment6(sgRNA_info):
	"""
	pytest --disable-warnings -sv tests/unit/test_standard_cut_alignments.py::test_alignment6
	1980, # aligned sequence normal, gRNA one gap, plus strand, aln length 21
	str(Seq(standard_group.genome['chr2'][21561760:21561777+1].tobytes().decode()).reverse_complement())
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [1980]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(21561727, 21561787)):
		mock_cutsite = CutSite(
			chromosome='chr2',
			strand='-',
			ref_position=21561757,
			cut_region= Seq('GGAACTCAAAAAATATAAGCATCACTGAGGATAATAATGGCAATGGATAGGATTTTAAAT').reverse_complement(),
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

	expected_protospacer_stop = 21_561_778-1 #results taken from blat -1 to account for 1-index

	expected_protospacer_start = 21_561_761-1 #results taken from blat -1 to account for 1-index

	print(f"{mock_cutsite.global_position['protospacer_stop']} vs {expected_protospacer_stop}")

	print(f"{mock_cutsite.global_position['protospacer_start']} vs {expected_protospacer_start}")

	assert mock_cutsite.global_position['protospacer_stop'] == expected_protospacer_stop #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == expected_protospacer_start #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 21_561_764-1 #results taken from blat -1 to account for 1-index


def test_alignment7(sgRNA_info):
	"""
	pytest --disable-warnings -sv tests/unit/test_standard_cut_alignments.py::test_alignment7
	1054	chr13	-	72141571	0	0	1	3	DACH1	274368	1.207177456	-0.044419777	3.625503467	0.040448677	0.23420457	0.001574803	0.509038908	0.192	110.6666667	10	302	20	72141568	72141587	9	28	AAAATGAGCAAACATTACTG	AAAATATGCAAACATCACTG	AGG	20	ACTTTTTAGAAAATGAGCAAACATTACTGAGGAACTAGTAATTTTATTTAATGCCATTTC
	str(Seq(standard_group.genome['chr13'][72141567:72141586+1].tobytes().decode()).reverse_complement())
	"""
	print(sgRNA_info)

	mock_detail = pd.DataFrame({'cut_cluster': [1054]})

	with patch('CRISPRito.StandardCuts.sequence_slice_locations', return_value=(72141536, 72141596)):
		mock_cutsite = CutSite(
			chromosome='chr13',
			strand='-',
			ref_position=72141566,
			cut_region= Seq('ACTTTTTAGAAAATGAGCAAACATTACTGAGGAACTAGTAATTTTATTTAATGCCATTTC').reverse_complement(),
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

	expected_protospacer_stop = 72_141_587-1 #results taken from blat -1 to account for 1-index

	expected_protospacer_start = 72_141_568-1 #results taken from blat -1 to account for 1-index

	print(f"{mock_cutsite.global_position['protospacer_stop']} vs {expected_protospacer_stop}")

	print(f"{mock_cutsite.global_position['protospacer_start']} vs {expected_protospacer_start}")

	assert mock_cutsite.global_position['protospacer_stop'] == expected_protospacer_stop #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['protospacer_start'] == expected_protospacer_start #results taken from blat -1 to account for 1-index

	assert mock_cutsite.global_position['cut'] == 72_141_571-1 #results taken from blat -1 to account for 1-index





# AAAATGAGCAAACATTACTG

# 'GTCATTACAAACGAGTAAAA'



