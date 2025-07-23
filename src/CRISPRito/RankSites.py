import traceback
import argparse
import os
from os.path import join as pjoin
import pandas as pd
from CRISPRito.StandardCutRank import StandardCutRank
from CRISPRito.Utils import setup_logging

# from CRISPRito.StandardCutRank import RankOperator

def rank_sites(
	group_samplesheet_path,
	rank_table_weights_path,
	cut_profiles_path,
	id_counts_path,
	method_counts_path,
	output_dir,
	output_name
	):
	"""
	main_dir='/data/friederike_herbst_nowrouzi_project/projects/base_editor_ptprc_project/data/processed/1-crisprito_analysis/0-make_standard_inputs/score_methods_CRISPRito_0_2_0'
	group_samplesheet_path = pjoin(main_dir,'1_group_samplesheet.csv')
	rank_table_weights_path = pjoin(main_dir,'1_group_id_rank_weight_skeleton_modified.csv')
	cut_profiles_path = pjoin(main_dir,'1_group_cut_profiles.csv')
	id_counts_path = pjoin(main_dir,'1_group_id_counts.csv')
	method_counts_path = pjoin(main_dir,'1_group_method_counts.csv')
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

	log_path = pjoin(args.output_dir, args.output_name + ".log")
	
	setup_logging(log_path)

	try:
		print("Starting ranking process...")

		rank_sites(
			group_samplesheet_path=args.group_samplesheet_path,
			rank_table_weights_path=args.rank_table_weights_path,
			cut_profiles_path=args.cut_profiles_path,
			id_counts_path=args.id_counts_path,
			method_counts_path=args.method_counts_path,
			output_dir=args.output_dir,
			output_name=args.output_name
		)

		print("Finished ranking process.")

	except Exception as e:
		logging.error("Unhandled exception occurred:")
		logging.error(traceback.format_exc())
		sys.exit(1)


if __name__ == '__main__':
	main()

