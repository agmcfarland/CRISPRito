

aln_units = [['AAAATATGCAAACATCACTG',	'AAAATATGCAAACATCACTG'],
	['AAAATAAGCAAACATTCACTG',	'AAAATATGCAAACA-TCACTG'],
	['AAAATGTGAAAACATCACTG',	'AAAATATGCAAACATCACTG'],
	['AAAACTGTGAAAACATTCACTG',	'AAAA-TATGCAAACAT-CACTG'],
	]


for i in aln_units:
	aligned_sequence = i[0]
	aligned_gRNA = i[1]
	mismatches = sum(a != b for a, b in zip(aligned_sequence, aligned_gRNA) if a != '-' and b != '-')
	print(mismatches)