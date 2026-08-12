import pytest
import numpy as np
import pandas as pd
from CRISPRito.AutoRankCuts import AutoRankCuts


# ---------------------------------------------------------------------------
# calculate_rra_scores
# ---------------------------------------------------------------------------

def test_calculate_rra_scores_top_consensus_cluster_ranks_first():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_calculate_rra_scores_top_consensus_cluster_ranks_first

	cluster 1 is top-scoring in both methods -> should get the smallest p_rra
	"""
	cut_id_detail = pd.DataFrame({
		'cut_cluster': [1, 2, 3, 1, 2, 3],
		'method':      ['A', 'A', 'A', 'B', 'B', 'B'],
		'score':       [10, 5, 1, 10, 1, 5]
	})

	arc = AutoRankCuts(cut_id_detail=cut_id_detail)
	rra = arc.calculate_rra_scores()

	assert set(rra.columns) >= {'cut_cluster', 'p_rra', 'log_p_rra'}
	assert len(rra) == 3
	# rra_scores is sorted ascending by p_rra, so first row is the top hit
	assert rra.iloc[0]['cut_cluster'] == 1
	assert (rra['p_rra'] > 0).all() and (rra['p_rra'] <= 1).all()


def test_calculate_rra_scores_undetected_cluster_treated_as_worst_rank():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_calculate_rra_scores_undetected_cluster_treated_as_worst_rank

	cluster 2 is only ever detected by method A, never by method B
	-> should be penalized relative to a cluster detected strongly by both
	"""
	cut_id_detail = pd.DataFrame({
		'cut_cluster': [1, 2, 1],
		'method':      ['A', 'A', 'B'],
		'score':       [5, 10, 5]
	})

	arc = AutoRankCuts(cut_id_detail=cut_id_detail)
	rra = arc.calculate_rra_scores()

	# cluster 1 is detected by both methods; cluster 2 only by one
	cluster1_p = rra.loc[rra['cut_cluster'] == 1, 'p_rra'].item()
	cluster2_p = rra.loc[rra['cut_cluster'] == 2, 'p_rra'].item()
	assert cluster1_p <= cluster2_p


# ---------------------------------------------------------------------------
# calculate_magnitude_score
# ---------------------------------------------------------------------------

def test_calculate_magnitude_score_requires_rra_scores_first():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_calculate_magnitude_score_requires_rra_scores_first
	"""
	cut_id_detail = pd.DataFrame({
		'cut_cluster': [1, 1, 2],
		'method':      ['A', 'B', 'A'],
		'score':       [10, 8, 3]
	})

	arc = AutoRankCuts(cut_id_detail=cut_id_detail)
	arc.calculate_rra_scores()
	mag = arc.calculate_magnitude_score(transform='percentile', aggregation='max')

	assert set(mag.columns) == {'cut_cluster', 'magnitude_score'}
	assert len(mag) == 2
	assert (mag['magnitude_score'] > 0).all()


def test_calculate_magnitude_score_missing_cluster_gets_floor_value():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_calculate_magnitude_score_missing_cluster_gets_floor_value

	rra_scores comes from clusters 1 and 2, but cut_id_detail only has scores
	for cluster 1 by the time magnitude is calculated separately
	"""
	cut_id_detail_full = pd.DataFrame({
		'cut_cluster': [1, 2],
		'method':      ['A', 'A'],
		'score':       [10, 5]
	})

	arc = AutoRankCuts(cut_id_detail=cut_id_detail_full)
	arc.calculate_rra_scores()

	# now simulate a magnitude table missing cluster 2 entirely
	arc.cut_id_detail = cut_id_detail_full[cut_id_detail_full['cut_cluster'] == 1]
	mag = arc.calculate_magnitude_score()

	missing_row = mag.loc[mag['cut_cluster'] == 2, 'magnitude_score']
	assert missing_row.item() == 0.1


# ---------------------------------------------------------------------------
# compute_method_support_summary
# ---------------------------------------------------------------------------

def test_compute_method_support_summary_counts_methods_and_replicates():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_compute_method_support_summary_counts_methods_and_replicates
	"""
	cut_id_detail = pd.DataFrame({
		'cut_cluster': [1, 1, 1, 2],
		'method':      ['A', 'A', 'B', 'A'],
		'id':          ['rep1', 'rep2', 'rep1', 'rep1'],
		'score':       [1, 2, 3, 4]
	})

	arc = AutoRankCuts(cut_id_detail=cut_id_detail)
	summary = arc.compute_method_support_summary()

	row1 = summary.loc[summary['cut_cluster'] == 1].iloc[0]
	row2 = summary.loc[summary['cut_cluster'] == 2].iloc[0]

	assert row1['total_methods'] == 2       # A and B both detected cluster 1
	assert row1['A'] == 2                    # 2 unique replicate ids for method A
	assert row1['B'] == 1
	assert row2['total_methods'] == 1
	assert row2['A'] == 1


# ---------------------------------------------------------------------------
# calculate_power_law_distance_decay
# ---------------------------------------------------------------------------

def test_power_law_decay_zero_distance_gives_weight_one():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_power_law_decay_zero_distance_gives_weight_one
	"""
	cut_profile = pd.DataFrame({
		'cut_cluster': [1, 2, 3],
		'nearest_gene_distance': [0, 5000, 10000]
	})

	arc = AutoRankCuts(cut_profile=cut_profile)
	out = arc.calculate_power_law_distance_decay('nearest_gene_distance', tau=5000.0, gamma=0.5, max_distance_bp=10000.0)

	weights = out.set_index('cut_cluster')['nearest_gene_distance_decay_weight']
	assert weights[1] == pytest.approx(1.0)
	# tau is defined as the distance where weight = 0.5
	assert weights[2] == pytest.approx(0.5, abs=1e-6)
	# weight should be strictly decreasing with distance
	assert weights[1] > weights[2] > weights[3]


def test_power_law_decay_unannotated_distance_treated_as_max():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_power_law_decay_unannotated_distance_treated_as_max

	-1 is the "not detected" sentinel and should decay like max_distance_bp
	"""
	cut_profile = pd.DataFrame({
		'cut_cluster': [1, 2],
		'nearest_gene_distance': [-1, 10000]
	})

	arc = AutoRankCuts(cut_profile=cut_profile)
	out = arc.calculate_power_law_distance_decay('nearest_gene_distance', max_distance_bp=10000.0)

	weights = out.set_index('cut_cluster')['nearest_gene_distance_decay_weight']
	assert weights[1] == pytest.approx(weights[2])


def test_power_law_decay_populates_feature_weight_across_multiple_calls():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_power_law_decay_populates_feature_weight_across_multiple_calls
	"""
	cut_profile = pd.DataFrame({
		'cut_cluster': [1, 2],
		'gene_distance': [0, 1000],
		'oncogene_distance': [2000, 0]
	})

	arc = AutoRankCuts(cut_profile=cut_profile)
	arc.calculate_power_law_distance_decay('gene_distance')
	arc.calculate_power_law_distance_decay('oncogene_distance')

	assert 'gene_distance_decay_weight' in arc.feature_weight.columns
	assert 'oncogene_distance_decay_weight' in arc.feature_weight.columns
	assert len(arc.feature_weight) == 2


# ---------------------------------------------------------------------------
# construct_crisprito_rank
# ---------------------------------------------------------------------------

def _build_full_autorankcuts():
	"""Helper: builds an AutoRankCuts instance with all upstream steps run."""
	cut_id_detail = pd.DataFrame({
		'cut_cluster': [1, 1, 2, 2],
		'method':      ['A', 'B', 'A', 'B'],
		'id':          ['rep1', 'rep1', 'rep1', 'rep1'],
		'score':       [10, 10, 2, 1]
	})
	cut_profile = pd.DataFrame({
		'cut_cluster': [1, 2],
		'gene_distance': [0, 9000]
	})

	arc = AutoRankCuts(cut_profile=cut_profile, cut_id_detail=cut_id_detail)
	arc.calculate_rra_scores()
	arc.calculate_magnitude_score()
	arc.compute_method_support_summary()
	arc.calculate_power_law_distance_decay('gene_distance')
	return arc


def test_construct_crisprito_rank_rejects_invalid_score_col():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_construct_crisprito_rank_rejects_invalid_score_col
	"""
	arc = _build_full_autorankcuts()
	with pytest.raises(ValueError):
		arc.construct_crisprito_rank(score_col='not_a_real_column')


def test_construct_crisprito_rank_rejects_mismatched_feature_weights():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_construct_crisprito_rank_rejects_mismatched_feature_weights
	"""
	arc = _build_full_autorankcuts()
	with pytest.raises(ValueError):
		arc.construct_crisprito_rank(feature_weights=[1.0, 2.0])  # only 1 decay col exists


def test_construct_crisprito_rank_top_cluster_ranked_first():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_construct_crisprito_rank_top_cluster_ranked_first

	cluster 1: strong consensus (both methods), max score, closest to gene
	cluster 2: weak on all three fronts
	-> cluster 1 should be crisprito_rank 1
	"""
	arc = _build_full_autorankcuts()
	result = arc.construct_crisprito_rank(score_col='rra_weighted_magnitude_score')

	assert set(result.columns) == {'cut_cluster', 'crisprito_rank', 'combined_priority_index', 'composite_decay_weight'}
	top_row = result.loc[result['crisprito_rank'] == 1].iloc[0]
	assert top_row['cut_cluster'] == 1


# ---------------------------------------------------------------------------
# prepare_combined_rank_output
# ---------------------------------------------------------------------------

def test_prepare_combined_rank_output_merges_all_components():
	"""
	pytest -sv tests/unit/test_AutoRankCuts.py::test_prepare_combined_rank_output_merges_all_components
	"""
	arc = _build_full_autorankcuts()
	arc.cut_profile['chromosome'] = ['chr1', 'chr1']
	arc.cut_profile['strand'] = ['+', '+']
	arc.cut_profile['cut'] = [100, 200]
	arc.cut_profile['nearest_gene'] = ['GENE1', 'GENE2']

	arc.construct_crisprito_rank(score_col='rra_weighted_magnitude_score')
	out = arc.prepare_combined_rank_output()

	expected_cols = {
		'cut_cluster', 'chromosome', 'strand', 'cut', 'nearest_gene',
		'crisprito_rank', 'combined_priority_index', 'composite_decay_weight',
		'p_rra', 'log_p_rra', 'magnitude_score',
		'gene_distance', 'gene_distance_decay_weight', 'total_methods', 'A', 'B'
	}
	assert expected_cols.issubset(set(out.columns))
	assert len(out) == 2
	assert out.iloc[0]['crisprito_rank'] == 1