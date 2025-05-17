from itertools import repeat
import pandas as pd
import os
from os.path import join as pjoin
from Bio.Seq import Seq
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from tqdm import tqdm
import pyranges as pr
from CRISPRito.Utils import (
	greedy_clustering_incremental,
	retrieve_genome_slices_memoryview,
	genome_to_dict_memoryview,
	scale_zscore,
	scale_min_max,
	parse_global_alignment,
	sequence_slice_locations,
	central_tendency,
	report_time,
	convert_df_to_granges,
	batch_overlaps,
	batch_nearest_feature,
	find_revised_pam,
	df_long_to_wide
	)
import statistics
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from CRISPRito.StandardCuts import CutSite, StandardCuts
from CRISPRito.Utils import (
	genome_to_dict_memoryview
	)


def range_to_site(
	group_samplesheet_path,
	output_filename,
	output_path,
	genome_path, 
	flank_size:int = 30, 
	sgRNA:str = '', 
	PAM_alignment:str = 'NGG',
	range_threshold = 1
	):

	# Update this
	PAM_alignment = PAM_alignment.replace('N', '-')

	# Update thisExtracting cut regions:
	allowed_chromosome = [f'chr{i}' for i in range(1,22)]
	allowed_chromosome.append('chrY')
	allowed_chromosome.append('chrX')

	strand_search = {
		'plus' : '+',
		'minus' : '-'
	}

	for strand_ in strand_search.keys():

		print(f'searching {strand_search[strand_]}')

		df_group_samplesheet = pd.read_csv(group_samplesheet_path) 

		df_group_samplesheet['strand'] = strand_search[strand_]

		standard_group = StandardCuts(sample_sheet = df_group_samplesheet, flank_size = flank_size, sgRNA = sgRNA, PAM_alignment = PAM_alignment)	

		standard_group.load_cut_sites()
		# Update this
		standard_group.df_cut_sites = standard_group.df_cut_sites[standard_group.df_cut_sites['chromosome'].isin(allowed_chromosome)]

		standard_group.cluster_cut_sites(range_threshold = range_threshold)

		standard_group.update_cut_cluster_id()

		standard_group.load_genome(genome_path = genome_path)

		standard_group.extract_cut_region()

		standard_group.build_cut_sites()

		standard_group.parallel_build_cut_site_alignment()

		standard_group.cut_detail_to_df(include_cut_location = True)
		
		standard_group.df_cut_detail.to_csv(pjoin(output_path, f'{standard_group.cluster_group}_group_id_cut_detail_{strand_}.csv'), index = None)


def main():

	parser = argparse.ArgumentParser(description="Process CRISPRito sample group.")

	parser.add_argument("--group_samplesheet_path", help="Path to the sample sheet CSV.")
	parser.add_argument("--output_dir", default="CRISPRito_output", help="Directory to save output CSV.")
	parser.add_argument("--genome_path", required=True, help="Path to the gzipped genome FASTA")
	parser.add_argument("--flank_size", type=int, default=30, help="Flank size around cut sites.")
	parser.add_argument("--sgRNA", type=str, default="", help="sgRNA sequence (optional).")
	parser.add_argument("--PAM_alignment", type=str, default="-GG", help="PAM sequence for alignment.")
	parser.add_argument("--range_threshold", type=int, default=20, help="Distance between clusters")

	args = parser.parse_args()

	range_to_site(
		group_samplesheet_path=args.group_samplesheet_path,
		output_dir=args.output_dir,
		genome_path=args.genome_path,
		flank_size=args.flank_size,
		sgRNA=args.sgRNA,
		PAM_alignment=args.PAM_alignment,
		range_threshold=args.range_threshold
	)


if __name__ == '__main__':
	main()







