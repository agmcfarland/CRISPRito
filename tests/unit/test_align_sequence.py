import pytest
import gzip
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
from CRISPRito.Utils import parse_global_alignment
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import DNA

import time

@pytest.fixture
def target_list():
	targets = [
	['chr1:+:198706754', 0, 'tttcttttagatgaaaaatatgcaaacatcactgtggattacttatataacaaggaaacta'.upper()],
	['chr2:-:143961596', 2,  'ttatgaaaatgacagagcattttccctcagtgatgttttcatatttataatttttgagaac'.upper()],
	['chr3:+:138494328', 4,  'ctcattatttgtcagtaatatacaagcatcactgaggacacttatgtttggaaattcttta'.upper()],
	['chr8:+:74745562', 5,  'attggagatgatagtctttatgtaaacatcactgtgggtttttttttcactgtaaataggc'.upper()],
	['chr7:+:50816217', 6,  'cccaaagtataggcaacataggcataaataaatgggtttagatcacactataaagcttcta'.upper()]
	]
	return targets


def test_aligner_works(target_list):
	"""
	pytest -sv tests/unit/test_align_sequence.py::test_aligner_works
	"""
	query_sequence = DNA('AAAATATGCAAACATCACTG')

	print('\n')

	result_list = [(14, 33), (-1, -1), (14, 33), (14, 33), (14, 33)]

	for e, t in enumerate(target_list):
		print(t)
		result = global_pairwise_align_nucleotide(query_sequence, DNA(t[2]))
		# print(result)
		parsed_result = parse_global_alignment(result)
		# print(parsed_result)
		assert parsed_result == result_list[e]

def test_aligner_works_reverse_complement(target_list):
	"""
	pytest -sv tests/unit/test_align_sequence.py::test_aligner_works_reverse_complement
	"""
	query_sequence = str(Seq('AAAATATGCAAACATCACTG').reverse_complement())

	query_sequence = DNA(query_sequence)

	print('\n')

	result_list = [(17, 19), (27, 46), (38, 57), (60, 60), (-1, -1)]

	for e, t in enumerate(target_list):
		print(t)
		result = global_pairwise_align_nucleotide(query_sequence, DNA(t[2]))
		print(result)
		parsed_result = parse_global_alignment(result)
		print(parsed_result)

		assert parsed_result == result_list[e]


def test_time_complexity(target_list):
	"""
	pytest -sv tests/unit/test_align_sequence.py::test_speed
	
	Testing how processing time increases within increasing number of sequences. Single processor.

	In general, 10 sequences = 0.5 seconds and this holds as sequence increases.
	"""
	query_sequence = DNA('AAAATATGCAAACATCACTG')
	
	for t in target_list:
		print('\n',t)
		start_time = time.time()
		for z in [1, 10]:#, 100, 1000, 10000]:
			for i in range(z):
				global_pairwise_align_nucleotide(query_sequence, DNA(t[2]))

				assert time.time()-start_time < 1


