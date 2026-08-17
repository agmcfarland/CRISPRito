import traceback
import logging
import sys
import argparse
import os
from os.path import join as pjoin
import pandas as pd
from CRISPRito.StandardCutRank import StandardCutRank
from CRISPRito.AutoRankCuts import AutoRankCuts
from CRISPRito.Utils import setup_logging, report_time


def resolve_group_paths(cluster_group, data_dir, workflow):
	"""
	Resolve the standard step-1 pipeline output filenames for a given cluster_group,
	so downstream steps only need --cluster_group and --data_dir rather than five
	separate explicit path arguments.

	Naming convention (must match step 1 output):
		{cluster_group}_group_cut_profiles.csv
		{cluster_group}_group_method_counts.csv
		{cluster_group}_group_id_counts.csv
		{cluster_group}_group_id_cut_detail.csv
		{cluster_group}_group_samplesheet.csv

	Only the files actually used by `workflow` are required to exist: 'auto' doesn't
	use method_counts/id_counts, and 'standard' doesn't use cut_id_detail.
	"""

	paths = {
		'cut_profiles_path': pjoin(data_dir, f'{cluster_group}_group_cut_profiles.csv'),
		'method_counts_path': pjoin(data_dir, f'{cluster_group}_group_method_counts.csv'),
		'id_counts_path': pjoin(data_dir, f'{cluster_group}_group_id_counts.csv'),
		'cut_id_detail_path': pjoin(data_dir, f'{cluster_group}_group_id_cut_detail.csv'),
		'group_samplesheet_path': pjoin(data_dir, f'{cluster_group}_group_samplesheet.csv'),
		}

	required_by_workflow = {
		'standard': ['cut_profiles_path', 'method_counts_path', 'id_counts_path', 'group_samplesheet_path'],
		'auto': ['cut_profiles_path', 'cut_id_detail_path', 'group_samplesheet_path'],
		}
	required = required_by_workflow[workflow]

	missing = [name for name in required if not os.path.isfile(paths[name])]
	if missing:
		raise FileNotFoundError(
			f"Could not find expected file(s) for cluster_group '{cluster_group}' in '{data_dir}': "
			+ ", ".join(f"{name} -> {paths[name]}" for name in missing)
			)

	return paths

@report_time
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

@report_time
def rank_sites_auto(
	cut_profiles_path,
	cut_id_detail_path,
	group_samplesheet_path,
	feature_driver_path,
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
	df_group_samplesheet = pd.read_csv(group_samplesheet_path)
	if feature_driver_path is not None:
		df_feature_driver = pd.read_csv(feature_driver_path)

	df_cut_id_detail = df_cut_id_detail.merge(
		df_group_samplesheet[['id', 'method']],
		how = 'left',
		on = 'id'
		)

	autoranker = AutoRankCuts(
		cut_id_detail = df_cut_id_detail,
		cut_profile = df_cut_profiles
		)

	autoranker.calculate_rra_scores()

	autoranker.calculate_magnitude_score(
		transform = magnitude_transform,
		aggregation = magnitude_aggregation
		)

	# automatically extract feature by annotation type
	if feature_driver_path is not None:
		distance_cols = df_feature_driver[df_feature_driver['type'] == 'annotation']['feature'].tolist()
		distance_cols = [f'nearest_{i}_distance' for i in distance_cols]

		for distance_col_ in distance_cols:
			autoranker.calculate_power_law_distance_decay(
				distance_col = distance_col_,
				tau = tau,
				gamma = gamma,
				max_distance_bp = max_distance_bp
				)

	autoranker.compute_method_support_summary()

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
		default = 'auto',
		help = "Which ranking workflow to run: 'standard' (StandardCutRank, user-weighted rule-based scoring) "
		       "or 'auto' (AutoRankCuts, RRA-based consensus ranking with optional magnitude/spatial enrichment). [Default auto]."
		)

	# Shared args
	parser.add_argument(
		"--cluster_group",
		required = True,
		help = "cluster_group identifier from step 1 of the pipeline. Used with --data_dir to auto-resolve "
		       "cut_profiles, method_counts, cut_id_detail, and samplesheet paths."
		)
	parser.add_argument(
		"--data_dir",
		required = True,
		help = "Directory containing step 1 pipeline outputs for this cluster_group "
		       "(e.g. '{cluster_group}_group_cut_profiles.csv', etc.)."
		)
	parser.add_argument("--output_dir", default = "CRISPRito_output", help = "Directory to save output CSV.")
	parser.add_argument("--output_name", default = "ranked_cut_sites", help = "Name of ranked sites table.")

	# StandardCutRank-specific args
	parser.add_argument("--rank_table_weights_path", help = "[standard] Path to the weights table.")

	# AutoRankCuts-specific args
	parser.add_argument("--feature_driver_path", default = None, help = "['auto] Path to the feature annotation CSV.")
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

		group_paths = resolve_group_paths(args.cluster_group, args.data_dir, args.workflow)

		if args.workflow == "standard":

			rank_sites_standard(
				group_samplesheet_path = group_paths['group_samplesheet_path'],
				rank_table_weights_path = args.rank_table_weights_path,
				cut_profiles_path = group_paths['cut_profiles_path'],
				id_counts_path = group_paths['id_counts_path'],
				method_counts_path = group_paths['method_counts_path'],
				output_dir = args.output_dir,
				output_name = args.output_name
				)

		elif args.workflow == "auto":

		
			rank_sites_auto(
				cut_profiles_path = group_paths['cut_profiles_path'],
				cut_id_detail_path = group_paths['cut_id_detail_path'],
				group_samplesheet_path = group_paths['group_samplesheet_path'],
				feature_driver_path = args.feature_driver_path,
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