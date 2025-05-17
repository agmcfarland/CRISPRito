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

		df_group_samplesheet = pd.read_csv(group_samplesheet_path) 

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

		standard_group.cut_detail_to_df()
		
		standard_group.df_cut_detail.to_csv(pjoin(output_path, f'{standard_group.cluster_group}_group_id_cut_detail_{strand_}.csv'), index = None)




# df_ranges = pd.read_csv(pjoin("/data/friederike_herbst_nowrouzi_project/projects/base_editor_ptprc_project/data/external/rhAMP_seq/nominated_rhamp_seq_coordinates_ptprc.csv"))

# df_ranges = df_ranges.rename(columns = {'Chrom': 'chromosome', 'ChromStart': 'start', 'ChromEnd': 'end', 'Name': 'id'})

# df_ranges = df_ranges[['id', 'chromosome', 'start', 'end']]

# df_ranges['ref_position'] = df_ranges.apply(lambda x: statistics.median((x['start'], x['end'])), axis = 1)

# df_ranges['size'] = df_ranges.apply(lambda x: x['end'] - x['start'], axis = 1)

# genome = genome_to_dict_memoryview(genome_path = '/data/GenomicTrackRepository/data/processed/hg38/hg38.fasta.gz')


# # alignments processed in order of self.sgRNA 
# sgRNA = 'AAAATATGCAAACATCACTG'
# PAM_alignment = '-GG'
# sgRNA = {
# 	'fwd_NGG' : DNA(sgRNA + PAM_alignment),
# 	'fwd' : DNA(sgRNA)
# }

# sgRNA_alignment_tolerance = {
# 	'fwd' : range(18,21),
# 	'fwd_NGG' : range(18,25)
# }

# sgRNA_alignment_start_offset = {
# 	'fwd' : 0,
# 	'fwd_NGG' : len(PAM_alignment)
# }


# for _, row in df_ranges.iterrows():break

# 	for strand in ('+', '-'):break

# 		sequence_region = 

# 		cut_site = CutSite(
# 			chromosome = row.chromosome,
# 			strand = strand,
# 			ref_position = row.ref_position,
# 			cut_region = {},
# 			sgRNA = sgRNA,
# 			sgRNA_alignment_tolerance = sgRNA_alignment_tolerance,
# 			sgRNA_alignment_start_offset = sgRNA_alignment_start_offset,
# 			cut_distance = 3,
# 			detail = pd.DataFrame({'cut_cluster': [1]}),
# 			flank_size = 30,
# 			)

# 		cut_site.find_best_sgRNA_alignment()
	
# 	break

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
		range_threshold=args.range_threshold
	)


if __name__ == '__main__':
	main()







