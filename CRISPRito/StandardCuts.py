import pandas as pd
import os
from os.path import join as pjoin
from Bio.Seq import Seq
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from CRISPRito.Utils import (
	greedy_clustering_incremental,
	retrieve_genome_slices_memoryview,
	genome_to_dict_memoryview,
	zscore,
	scale_min_max,
	parse_global_alignment
	)


class StandardCuts:

	def __init__(self, sample_sheet:pd.DataFrame, flank_size:int = 30, sgRNA:str = ''):
		self.sample_sheet = sample_sheet
		self.flank_size = flank_size
		self.sgRNA = {
			'fwd' : DNA(sgRNA),
			'fwd_NGG' : DNA(sgRNA + '-GG')
		}
		self.sgRNA_alignment_tolerance = {
			'fwd' : range(18,21),
			'fwd_NGG' : range(19,24)
		}
		

	def load_cut_sites(self):

		self.df_cut_sites = pd.DataFrame()
		
		for _, row in self.sample_sheet.iterrows():
			
			df_sample = pd.read_csv(pjoin(row.standard_format_file_path))
			
			df_sample['id'] = row['id']
			
			self.df_cut_sites = pd.concat([self.df_cut_sites, df_sample])


	def update_cut_cluster_id(self):
		"""
		"""
		self.df_cut_sites = self.df_cut_sites.sort_values(by=["chromosome", "strand", "position"])

		self.df_cut_sites['cut_cluster'] = self.df_cut_sites.groupby(['chromosome', 'strand', 'cut_cluster']).ngroup()

	
	def cluster_cut_sites(self, range_threshold = 20):

		self.df_cut_sites = self.df_cut_sites.sort_values(by=["chromosome", "strand", "position"])

		grouped_cuts = self.df_cut_sites.groupby(['chromosome', 'strand'])

		df_clustered_cuts = pd.DataFrame()

		cluster_id_start_ = 0
		
		for _, df_group in grouped_cuts:

			cluster_groups = greedy_clustering_incremental(df_group['position'].tolist(), range_threshold = range_threshold)

			position_to_group = {pos: group for group, positions in cluster_groups.items() for pos in positions}

			df_group['cut_cluster'] = df_group['position'].map(position_to_group)

			df_clustered_cuts = pd.concat([df_clustered_cuts, df_group])


		if len(self.df_cut_sites) != len(df_clustered_cuts):

			raise ValueError('Not all cut sites had a cut group assigned')

		self.df_cut_sites = df_clustered_cuts 

	def standardize_scores(self):

		self.df_cut_sites['zscore'] = self.df_cut_sites.groupby('id')['score'].transform(lambda x: zscore(x.tolist(), ddof = 0))

		self.df_cut_sites['min_max'] = self.df_cut_sites.groupby('id')['score'].transform(lambda x: scale_min_max(x.tolist()))


	def load_genome(self, genome_path):
		if not os.path.exists(genome_path):
			raise FileNotFoundError(f"Genome file not found: {genome_path}")
		self.genome = genome_to_dict_memoryview(genome_path)

	def get_genome_size(self):
		self.genome_size = {}
		for k, v in self.genome.items():
			self.genome_size[k] = len(v)

	def extract_cut_region(self):
		"""
		Iterate through list of positions, grouped by chromosome, and extract self.flank size bp from the mean predicted position 
		"""

		df_reference = self.df_cut_sites.copy()

		df_reference['reference_position'] = df_reference.groupby('cut_cluster')['position'].transform('mean').astype(int)

		df_reference = df_reference[['chromosome', 'strand', 'cut_cluster', 'reference_position']].drop_duplicates()

		chromosome_groups = df_reference.groupby(['chromosome'])

		df_sequence = pd.DataFrame()

		for _, df_group in chromosome_groups:

			chromosome_ = df_group['chromosome'].unique().item()

			slices = retrieve_genome_slices_memoryview(sequence = self.genome[chromosome_], positions = df_group['reference_position'], flank_size = self.flank_size)

			df_group['cut_region'] = df_group['reference_position'].map(slices)

			df_sequence = pd.concat([df_sequence, df_group])

		self.df_reference_cut_sites = df_sequence


	def build_cut_sites(self):
		"""
		"""

		self.cut_sites = []

		for _, row in self.df_reference_cut_sites.iterrows():
			
			df_cluster_subset = self.df_cut_sites[self.df_cut_sites['cut_cluster'] == row.cut_cluster]
			
			df_cluster_subset = df_cluster_subset[['position', 'score', 'id']]

			standard_cut = CutSite(
				chromosome = row.chromosome,
				strand = row.strand,
				ref_position = row.reference_position,
				cut_region = row.cut_region,
				sgRNA = self.sgRNA,
				sgRNA_alignment_tolerance = self.sgRNA_alignment_tolerance,
				detail = df_cluster_subset
				)

			self.cut_sites.append(standard_cut)




	# def extract_cut_region():
		# pass

		


class CutSite:

	def __init__(self, chromosome, strand, ref_position, cut_region, sgRNA, sgRNA_alignment_tolerance, detail):
		self.chromosome = chromosome
		self.strand = strand
		self.ref_position = ref_position
		self.cut_region = cut_region
		self.sgRNA = sgRNA
		self.sgRNA_alignment_tolerance = sgRNA_alignment_tolerance
		self.detail = detail
		self.alignment = {}

	def __len__(self):
		return len(self.detail) 

	def __repr__(self):
		return f"CutSite(chrom={self.chromosome}, strand={self.strand}, ref_pos={self.ref_position}, cuts={len(self)})"

	def find_best_sgRNA_alignment(self):

		if self.strand == '-':
			sequence = DNA(str(Seq(self.cut_region).reverse_complement()))
		else:
			sequence = DNA(self.cut_region)

		for seqname, query_seq in self.sgRNA.items():
			self.alignment['sgRNA'] = seqname

			self.alignment['alignment'] = global_pairwise_align_nucleotide(query_seq, sequence)

			self.alignment['alignment_start'], self.alignment['alignment_stop'] = parse_global_alignment(self.alignment['alignment'])

			self.alignment['alignment_length'] = self.alignment['alignment_stop'] - self.alignment['alignment_start']

			if self.alignment['alignment_length'] in self.sgRNA_alignment_tolerance:
				break

	def print_alignment_stats(self):
		for k,v in self.alignment.items():
			print(f'{k}: {v}')


	def calculate():
		pass

	# def 






