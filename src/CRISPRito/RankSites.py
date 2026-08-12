import traceback
import logging
import sys
import argparse
import os
from os.path import join as pjoin
import pandas as pd
from CRISPRito.StandardCutRank import StandardCutRank
from CRISPRito.AutoRankCuts import AutoRankCuts
from CRISPRito.Utils import setup_logging


def rank_sites_standard(
	group_samplesheet_path,
	rank_table_weights_path,
	cut_profiles_path,
	id_counts_path,
	method_counts_path,
	output_dir,
	output_name
	):
	"""
	StandardCutRank workflow: rule-based additive scoring from a user-supplied weight table.
	"""

	df_group_samplesheet = pd.read_csv(group_samplesheet_path)
	df_rank_table_weights = pd.read_csv(rank_table_weights_path)
	df_cut_profiles = pd.read_csv(cut_profiles_path)

	df_id_counts = None
	df_method_counts = None
	if id_counts_path is not None:
		df_id_counts = pd.read_csv(id_counts_path)
	if method_counts_path is not None:
		df_method_counts = pd.read_csv(method_counts_path)

	sc_ranks = StandardCutRank(
		rank_table_weights = df_rank_table_weights,
		cut_profiles = df_cut_profiles,
		id_counts = df_id_counts,
		method_counts = df_method_counts,
		samplesheet = df_group_samplesheet
		)

	sc_ranks.get_score_from_rank_table()
	sc_ranks.tally_cut_cluster_scores()

	sc_ranks.cut_cluster_scores.to_csv(pjoin(output_dir, output_name + '.csv'), index = None)
	sc_ranks.rank_table_weights.to_csv(pjoin(output_dir, output_name + '_weights_metadata.csv'), index = None)


def rank_sites_auto(
	cut_profiles_path,
	cut_id_detail_path,
	distance_cols,
	magnitude_transform,
	magnitude_aggregation,
	tau,
	gamma,
	max_distance_bp,
	feature_weights,
	spatial_aggregation,
	score_col,
	output_dir,
	output_name
	):
	"""
	AutoRankCuts workflow: RRA-based consensus ranking with optional magnitude and
	spatial-proximity modulation. score_col='log_p_rra' gives a fully parameter-free
	ranking; the other options layer on optional, transparent enrichment.
	"""

	df_cut_profiles = pd.read_csv(cut_profiles_path)
	df_cut_id_detail = pd.read_csv(cut_id_detail_path)

	autoranker = AutoRankCuts(
		cut_id_detail = df_cut_id_detail,
		cut_profile = df_cut_profiles
		)

	autoranker.calculate_rra_scores()

	autoranker.calculate_magnitude_score(
		transform = magnitude_transform,
		aggregation = magnitude_aggregation
		)

	if distance_cols:
		for distance_col_ in distance_cols:
			autoranker.calculate_power_law_distance_decay(
				distance_col = distance_col_,
				tau = tau,
				gamma = gamma,
				max_distance_bp = max_distance_bp
				)

	autoranker.construct_crisprito_rank(
		feature_weights = feature_weights,
		aggregation = spatial_aggregation,
		score_col = score_col
		)

	df_output = autoranker.prepare_combined_rank_output()

	df_output.to_csv(pjoin(output_dir, output_name + '.csv'), index = None)


def main():
	parser = argparse.ArgumentParser(description="Rank sites generated from a CRISPRito output")

	parser.add_argument(
		"--workflow",
		choices = ["standard", "auto"],
		required = True,
		help = "Which ranking workflow to run: 'standard' (StandardCutRank, user-weighted rule-based scoring) "
		       "or 'auto' (AutoRankCuts, RRA-based consensus ranking with optional magnitude/spatial enrichment)."
		)

	# Shared args
	parser.add_argument("--cut_profiles_path", required = True)
	parser.add_argument("--output_dir", default = "CRISPRito_output", help = "Directory to save output CSV.")
	parser.add_argument("--output_name", default = "ranked_cut_sites", help = "Name of ranked sites table.")

	# StandardCutRank-specific args
	parser.add_argument("--group_samplesheet_path", help = "[standard] Path to the sample sheet CSV.")
	parser.add_argument("--rank_table_weights_path", help = "[standard] Path to the weights table.")
	parser.add_argument("--id_counts_path", help = "[standard] Path to id_counts CSV.")
	parser.add_argument("--method_counts_path", help = "[standard] Path to method_counts CSV.")

	# AutoRankCuts-specific args
	parser.add_argument("--cut_id_detail_path", help = "[auto] Path to the per-detection cut_id_detail CSV.")
	parser.add_argument(
		"--distance_cols",
		nargs = "*",
		default = [],
		help = "[auto] One or more distance columns in cut_profiles to compute power-law decay weights for "
		       "(e.g. nearest_gene_distance nearest_oncogene_distance)."
		)
	parser.add_argument(
		"--magnitude_transform",
		choices = ["zscore", "minmax", "percentile", "raw"],
		default = "percentile",
		help = "[auto] Per-method score normalization used in calculate_magnitude_score."
		)
	parser.add_argument(
		"--magnitude_aggregation",
		choices = ["mean", "max", "sum"],
		default = "max",
		help = "[auto] How per-method magnitude scores are aggregated per cut_cluster."
		)
	parser.add_argument("--tau", type = float, default = 5000.0, help = "[auto] Power-law decay midpoint distance (bp).")
	parser.add_argument("--gamma", type = float, default = 0.5, help = "[auto] Power-law decay tail-shape parameter.")
	parser.add_argument("--max_distance_bp", type = float, default = 10000.0, help = "[auto] Distance cap (bp) applied before decay.")
	parser.add_argument(
		"--feature_weights",
		nargs = "*",
		type = float,
		default = None,
		help = "[auto] Weights for combining decay columns, in the same order as --distance_cols. "
		       "Defaults to equal weighting if omitted."
		)
	parser.add_argument(
		"--spatial_aggregation",
		choices = ["harmonic", "weighted_mean", "max"],
		default = "weighted_mean",
		help = "[auto] How per-feature decay weights are combined into a single composite_decay_weight."
		)
	parser.add_argument(
		"--score_col",
		choices = ["log_p_rra", "rra_weighted_magnitude_score"],
		default = "rra_weighted_magnitude_score",
		help = "[auto] Base score multiplied by composite_decay_weight to form the final rank. "
		       "'log_p_rra' gives a fully parameter-free ranking using RRA consensus alone."
		)

	args = parser.parse_args()

	os.makedirs(args.output_dir, exist_ok = True)
	log_path = pjoin(args.output_dir, args.output_name + ".log")
	setup_logging(log_path)

	try:
		print(f"Starting ranking process ({args.workflow})...")

		if args.workflow == "standard":
			required = {
				"group_samplesheet_path": args.group_samplesheet_path,
				"rank_table_weights_path": args.rank_table_weights_path,
				}
			missing = [k for k, v in required.items() if v is None]
			if missing:
				parser.error(f"--workflow standard requires: {', '.join('--' + m for m in missing)}")

			rank_sites_standard(
				group_samplesheet_path = args.group_samplesheet_path,
				rank_table_weights_path = args.rank_table_weights_path,
				cut_profiles_path = args.cut_profiles_path,
				id_counts_path = args.id_counts_path,
				method_counts_path = args.method_counts_path,
				output_dir = args.output_dir,
				output_name = args.output_name
				)

		elif args.workflow == "auto":
			if args.cut_id_detail_path is None:
				parser.error("--workflow auto requires --cut_id_detail_path")

			if args.feature_weights and len(args.feature_weights) != len(args.distance_cols):
				parser.error("--feature_weights length must match --distance_cols length")

			rank_sites_auto(
				cut_profiles_path = args.cut_profiles_path,
				cut_id_detail_path = args.cut_id_detail_path,
				distance_cols = args.distance_cols,
				magnitude_transform = args.magnitude_transform,
				magnitude_aggregation = args.magnitude_aggregation,
				tau = args.tau,
				gamma = args.gamma,
				max_distance_bp = args.max_distance_bp,
				feature_weights = args.feature_weights,
				spatial_aggregation = args.spatial_aggregation,
				score_col = args.score_col,
				output_dir = args.output_dir,
				output_name = args.output_name
				)

		print("Finished ranking process.")

	except Exception as e:
		logging.error("Unhandled exception occurred:")
		logging.error(traceback.format_exc())
		sys.exit(1)


if __name__ == '__main__':
	main()