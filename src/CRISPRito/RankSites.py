import argparse
import os
from os.path import join as pjoin
import pandas as pd
from CRISPRito.StandardCutRank import StandardCutRank

def rank_sites(
	group_samplesheet_path,
	rank_table_weights_path,
	cut_profiles_path,
	id_counts_path,
	method_counts_path,
	output_dir,
	output_name
	):
	
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

	sc_ranks.cut_cluster_scores.to_csv(pjoin(output_dir, output_name +'.csv'), index = None)

	sc_ranks.rank_table_weights.to_csv(pjoin(output_dir, output_name +'_weights_metadata.csv'), index = None)


def main():

	parser = argparse.ArgumentParser(description="Rank sites generated from a CRISPRito output")

	parser.add_argument("--group_samplesheet_path", help="Path to the sample sheet CSV.")
	parser.add_argument("--rank_table_weights_path", help="Path to the weights.")
	parser.add_argument("--cut_profiles_path")
	parser.add_argument("--id_counts_path")
	parser.add_argument("--method_counts_path")
	parser.add_argument("--output_dir", default="CRISPRito_output", help="Directory to save output CSV.")
	parser.add_argument("--output_name", default="ranked_cut_sites", help="Name of ranked sites table.")

	args = parser.parse_args()

	rank_sites(
		group_samplesheet_path = args.group_samplesheet_path,
		rank_table_weights_path = args.rank_table_weights_path,
		cut_profiles_path = args.cut_profiles_path,
		id_counts_path = args.id_counts_path,
		method_counts_path = args.method_counts_path,
		output_dir = args.output_dir,
		output_name = args.output_name
	)


if __name__ == '__main__':
	main()

