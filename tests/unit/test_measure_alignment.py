

# aln_units = [
# 	['AAAATATGCAAACATCACTG',	'AAAATATGCAAACATCACTG'],
# 	['AAAATAAGCAAACATTCACTG',	'AAAATATGCAAACA-TCACTG'],
# 	['AAAATGTGAAAACATCACTG',	'AAAATATGCAAACATCACTG'],
# 	['AAAACTGTGAAAACATTCACTG',	'AAAA-TATGCAAACAT-CACTG'],
# 	]


# for i in aln_units:
# 	aligned_sequence = i[0]
# 	aligned_gRNA = i[1]
# 	mismatches = sum(a != b for a, b in zip(aligned_sequence, aligned_gRNA) if a != '-' and b != '-')
# 	gaps = aligned_gRNA.count('-')
# 	print(mismatches)
# 	print(gaps)

# import pytest
# import gzip
# import tempfile
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from io import StringIO
# from CRISPRito.Utils import parse_global_alignment
# from skbio.alignment import global_pairwise_align_nucleotide
# from skbio import DNA



# sgrna=DNA('AAAATATGCAAACATCACTG-GG')

# region=DNA('CCCTAATTGAGGTAAAATTAGCAAACAAAAATTGTGTATATTCAAGGTGTACAAAATGAT')

# global_pairwise_align_nucleotide(sgrna, region)

# for standard_cut in standard_group.cut_sites:
# 	if standard_cut.profile['nearest_gene'] == 'XYLT1':
# 		break

# """
# In [9]: standard_cut.alignment
# Out[9]:
# {'sgRNA': 'fwd',
#  'alignment': (TabularMSA[DNA]
#   --------------------------------------------------------------
#   Stats:
#       sequence count: 2
#       position count: 62
#   --------------------------------------------------------------
#   ---------------------------AAAATATGCAAACATCACTG---------------
#   ATCATCCCTCGCTGATAACCACTAGCCCAAAT--GCAAACATCACTGAGGCACGAGGGCCTT,
#   np.float64(8.0),
#   [(0, 19), (0, 59)]),
#  'local_start': 27,
#  'local_stop': 46,
#  'aligned_sequence': 'CAAATGCAAACATCACTGAG', # Im extracting from genomic and not from local
#  'aligned_gRNA': 'AAAATATGCAAACATCACTG',
#  'alignment_length': 20,
#  'PAM': 'GCA',
#  'local_cut': 43}

# In [10]: 'ATCATCCCTCGCTGATAACCACTAGCCCAAAT--GCAAACATCACTGAGGCACGAGGGCCTT'[27:46+1]
# Out[10]: 'CAAAT--GCAAACATCACTG'

# """