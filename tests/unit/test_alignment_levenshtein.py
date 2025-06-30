from CRISPRito.Utils import alignment_levenshtein

def test_one():
	"""
	pytest -sv tests/unit/test_alignment_levenshtein.py::test_one
	"""

	test_pairs = [
		['CT-GACCTAAATCCTGGACAAG', 'CTTGACCAATAGCCTTGACA-G'],
		['CTGTAGGCAGAGCCTTGACAAG', 'CTTGACCAATAGCCTTGACA-G'],
		['CT--ACCAATAGCCTTGACA' , 'CTTGACCAATAGCCTTGACA']
	]

	test_pairs_score = [7, 6, 2]

	for e, pair_ in enumerate(test_pairs):

		result = alignment_levenshtein(
			observed = pair_[0],
			reference = pair_[1]
			)

		print(result)

		assert result == test_pairs_score[e]


