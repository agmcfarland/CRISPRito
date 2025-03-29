import pytest
import gzip
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
from CRISPRito.Utils import align_sequence
from skbio.alignment import StripedSmithWaterman

@pytest.fixture
def target_list():
	targets = [['chr1:+:198706754', 0, 'tttcttttagatgaaaaatatgcaaacatcactgtggattacttatataacaaggaaacta'],
	['chr2:-:143961596', 2,  'ttatgaaaatgacagagcattttccctcagtgatgttttcatatttataatttttgagaac'],
	['chr3:+:138494328', 4,  'ctcattatttgtcagtaatatacaagcatcactgaggacacttatgtttggaaattcttta'],
	['chr8:+:74745562', 5,  'attggagatgatagtctttatgtaaacatcactgtgggtttttttttcactgtaaataggc'],
	['chr7:+:50816217', 6,  'cccaaagtataggcaacataggcataaataaatgggtttagatcacactataaagcttcta']]
	return targets


def test_aligner_works(target_list):
	"""
	pytest -sv tests/unit/test_align_sequence.py::test_aligner_works
	"""
	targets = target_list

	aligner = StripedSmithWaterman('AAAATATGCAAACATCACTG', gap_open_penalty=2, gap_extend_penalty=1, score_size = 0)
	# result = aligner(targets[0][2])
	# for k,v in result.items():
	# 	print(k, v)


	print('\n')

	for t in targets:
		query = 'GAGTAGCGCGAGCACAGCTA'
		result = align_sequence(aligner, t[2])
		print(t)

		print('optimal_alignment_score: ', result.optimal_alignment_score)
		print('suboptimal_alignment_score: ', result.suboptimal_alignment_score)
		print('query_begin: ', result.query_begin)
		print('query_end: ', result.query_end)
		print('target_begin: ', result.target_begin)
		print('target_end_optimal: ', result.target_end_optimal)
		print('target_end_suboptimal: ', result.target_end_suboptimal)
		print('cigar: ', result.cigar)
		print('query_sequence: ', result.query_sequence)
		print('target_sequence: ', result.target_sequence)
		print('\n')





	# result = align_sequence(aligner, target)

	# for i in 

	# print(result)


# for i in ['chr1:+:198706754 ', 'chr2:-:143961596', 'chr3:+:138494328', 'chr8:+:74745562', 'chr7:+:50816217']:
# 	chromosome = i.split(':')[0]
# 	strand = i.split(':')[1]
# 	position = int(i.split(':')[2])
# 	upstream = position - 30
# 	downstream = position + 30
# 	print(f'https://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment={chromosome}:{upstream},{downstream}')

# [['chr1:+:198706754', 0, 'tttcttttagatgaaaaatatgcaaacatcactgtggattacttatataacaaggaaacta'],
# ['chr2:-:143961596', 2,  'ttatgaaaatgacagagcattttccctcagtgatgttttcatatttataatttttgagaac'],
# ['chr3:+:138494328', 4,  'ctcattatttgtcagtaatatacaagcatcactgaggacacttatgtttggaaattcttta'],
# ['chr8:+:74745562', 5,  'attggagatgatagtctttatgtaaacatcactgtgggtttttttttcactgtaaataggc'],
# ['chr7:+:50816217', 6,  'cccaaagtataggcaacataggcataaataaatgggtttagatcacactataaagcttcta']]

# https://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr1:198706724,198706784