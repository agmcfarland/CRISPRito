import argparse
import os
from os.path import join as pjoin
import pandas as pd
from CRISPRito.StandardCuts import StandardCuts
from CRISPRito.FeatureManager import FeatureManager

def process_group(
	group_samplesheet_path,
	output_path,
	genome_path, 
	# feature_path,
	# gene_names_path,
	flank_size:int = 30, 
	sgRNA:str = '', 
	PAM_alignment:str = 'NGG',
	range_threshold = 20,
	feature_table_path = None
	):
	# troubleshooting inputs start
	# group_samplesheet_path = '/data/friederike_herbst_nowrouzi_project/projects/base_editor_ptprc_project/data/processed/1-crisprito_analysis/0-make_standard_inputs/iguide_bfx_rhamp/1_group_samplesheet.csv'
	# output_path = None
	# genome_path = '/data/GenomicTrackRepository/data/processed/hg38/hg38.fasta.gz'
	# feature_path = '/data/GenomicTrackRepository/data/processed/hg38/ncbiRefSeqCurated_expanded.csv'
	# gene_names_path = '/data/GenomicTrackRepository/data/external/hg38/gene_names.csv'
	# flank_size = 30
	# sgRNA = 'AAAATATGCAAACATCACTG'
	# PAM_alignment:str = 'NGG'

	# standard_group.df_cut_sites[standard_group.df_cut_sites['position'] == 28359447]
	# 3050	chr21	-	5525780	0	0	1	1	LOC102724159

	# troubleshooting inputs end

	# Update this
	PAM_alignment = PAM_alignment.replace('N', '-')

	# Update thisExtracting cut regions:
	allowed_chromosome = [f'chr{i}' for i in range(1,22)]
	allowed_chromosome.append('chrY')
	allowed_chromosome.append('chrX')

	# Update this variable name
	df_group_samplesheet = pd.read_csv(group_samplesheet_path) #'/data/CRISPRito/tests/temp/new/1_group_samplesheet.csv')

	standard_group = StandardCuts(sample_sheet = df_group_samplesheet, flank_size = flank_size, sgRNA = sgRNA, PAM_alignment = PAM_alignment)	

	feature_manager = FeatureManager(feature_table = pd.read_csv(feature_table_path))

	standard_group.load_cut_sites()
	# Update this
	standard_group.df_cut_sites = standard_group.df_cut_sites[standard_group.df_cut_sites['chromosome'].isin(allowed_chromosome)]

	standard_group.cluster_cut_sites(range_threshold = range_threshold)

	standard_group.update_cut_cluster_id()

	standard_group.remove_cluster_duplicates()

	standard_group.standardize_scores()

	standard_group.cut_cluster_representation()

	standard_group.load_genome(genome_path = genome_path)

	standard_group.extract_cut_region()

	standard_group.build_cut_sites()

	standard_group.parallel_build_cut_site_alignment()

	# Update this
	standard_group.assign_features(feature_manager = feature_manager)

	standard_group.build_cut_profile()

	standard_group.cut_profiles_to_df()

	standard_group.cut_detail_to_df()

	standard_group.df_cut_profiles = standard_group.df_cut_profiles.sort_values(['local_rank_median', 'zscore_mean'], ascending=[True, False])

	standard_group.df_cut_profiles.to_csv(pjoin(output_path, f'{standard_group.cluster_group}_group_cut_profiles.csv'), index = None)

	standard_group.method_counts.to_csv(pjoin(output_path, f'{standard_group.cluster_group}_group_method_counts.csv'), index = None)

	standard_group.id_counts.to_csv(pjoin(output_path, f'{standard_group.cluster_group}_group_id_counts.csv'), index = None)

	standard_group.df_cut_detail.to_csv(pjoin(output_path, f'{standard_group.cluster_group}_group_id_cut_detail.csv'), index = None)

	# Update this

	print(standard_group.method_counts)

	print(standard_group.id_counts)

	print(standard_group.df_cut_profiles)


def main():

	parser = argparse.ArgumentParser(description="Process CRISPRito sample group.")

	parser.add_argument("--group_samplesheet_path", help="Path to the sample sheet CSV.")
	parser.add_argument("--output_dir", default="CRISPRito_output", help="Directory to save output CSV.")
	parser.add_argument("--genome_path", required=True, help="Path to the gzipped genome FASTA")
	# parser.add_argument("--feature_path", required=True, help="Path to the features CSV file.")
	# parser.add_argument("--gene_names_path", required=True, help="Path to the gene names CSV file.")
	parser.add_argument("--flank_size", type=int, default=30, help="Flank size around cut sites.")
	parser.add_argument("--sgRNA", type=str, default="", help="sgRNA sequence (optional).")
	parser.add_argument("--PAM_alignment", type=str, default="-GG", help="PAM sequence for alignment.")
	parser.add_argument("--range_threshold", type=int, default=20, help="Distance between clusters")
	parser.add_argument("--feature_table_path", type=str, help="Path to feature table with features data.")

	args = parser.parse_args()

	process_group(
		group_samplesheet_path=args.group_samplesheet_path,
		output_dir=args.output_dir,
		genome_path=args.genome_path,
		# feature_path=args.feature_path,
		# gene_names_path=args.gene_names_path,
		flank_size=args.flank_size,
		sgRNA=args.sgRNA,
		PAM_alignment=args.PAM_alignment,
		range_threshold=args.range_threshold,
		feature_table_path = args.feature_table_path
	)


if __name__ == '__main__':
	main()

