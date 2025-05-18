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


class StandardCuts:

	def __init__(self, sample_sheet:pd.DataFrame, flank_size:int = 30, sgRNA:str = '', PAM_alignment:str = '-GG', cut_distance = 3):
		self.sample_sheet = sample_sheet
		self.flank_size = flank_size

		self.cluster_group = self.sample_sheet.cluster_group.unique().item()

		# alignments processed in order of self.sgRNA 
		self.sgRNA = {
			'fwd_NGG' : DNA(sgRNA + PAM_alignment),
			'fwd' : DNA(sgRNA)
		}

		self.sgRNA_alignment_tolerance = {
			'fwd' : range(18,21),
			'fwd_NGG' : range(18,25)
		}

		self.sgRNA_alignment_start_offset = {
			'fwd' : 0,
			'fwd_NGG' : len(PAM_alignment)
		}
		self.cut_distance = cut_distance

		self.df_cut_sites = pd.DataFrame()

		self.df_cut_profiles = None

		self.cut_sites = []
		
		self.method_counts = pd.DataFrame()
		
		self.id_counts = pd.DataFrame()

		self.df_cut_detail = pd.DataFrame()

	def __repr__(self):
		return f'StandardCuts\nSample Sheet:\n{self.sample_sheet}\nCut Sites\n{self.df_cut_sites}'

	def __len__(self):
		return len(self.df_cut_sites)

	def load_cut_sites(self):
		
		for _, row in self.sample_sheet.iterrows():
			
			df_sample = pd.read_csv(pjoin(row.standard_format_file_path))
			
			df_sample['id'] = row['id']
			
			self.df_cut_sites = pd.concat([self.df_cut_sites, df_sample])

		self.df_cut_sites = self.df_cut_sites.merge(self.sample_sheet[['id', 'measurement_type']], on = 'id')

	def remove_cluster_duplicates(self):
		"""
		Keep highest scoring 
		"""
		self.df_cut_sites = (
			self.df_cut_sites.sort_values('score', ascending=False)
			.drop_duplicates(subset=['cut_cluster', 'id'])
			.reset_index(drop=True)
		)

	def update_cut_cluster_id(self):
		"""
		"""
		self.df_cut_sites = self.df_cut_sites.sort_values(by=["chromosome", "strand", "position"])

		self.df_cut_sites['cut_cluster'] = self.df_cut_sites.groupby(['chromosome', 'strand', 'cut_cluster']).ngroup()


	def cut_cluster_representation(self):
		# List of (attribute name, grouping key) pairs
		groupings = [
			('method_counts', 'method'),
			('id_counts', 'id'),
		]

		for attr_name, group_key in groupings:
			df = self.sample_sheet[['sample', 'method', 'id']].copy()
			df = df.merge(self.df_cut_sites[['cut_cluster', 'id']].drop_duplicates(), on='id')

			df_wide = df_long_to_wide(df, to_rows=group_key, to_columns='cut_cluster', column_prefix = 'cut')
			setattr(self, attr_name, df_wide)
			
	
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

		self.df_cut_sites['zscore'] = self.df_cut_sites.groupby('id')['score'].transform(lambda x: scale_zscore(x.tolist(), degrees_of_freedom = 0))

		self.df_cut_sites['min_max'] = self.df_cut_sites.groupby('id')['score'].transform(lambda x: scale_min_max(x.tolist()))

		self.df_cut_sites['min_max'] = self.df_cut_sites.apply(lambda x: x['score'] if x['measurement_type'] == 'one_scaled' else x['min_max'], axis = 1)

		self.df_cut_sites = (self.df_cut_sites.sort_values(['id', 'score'], ascending=[True, False]).assign(local_rank=lambda x: x.groupby('id').cumcount() + 1))

	@report_time
	def load_genome(self, genome_path):
		if not os.path.exists(genome_path):
			raise FileNotFoundError(f"Genome file not found: {genome_path}")
		self.genome = genome_to_dict_memoryview(genome_path)

	def get_genome_size(self):
		self.genome_size = {}
		for k, v in self.genome.items():
			self.genome_size[k] = len(v)


	@report_time
	def extract_cut_region(self):
		"""
		Iterate through list of positions, grouped by chromosome, and extract self.flank size bp from the mean predicted position 
		"""

		df_reference = self.df_cut_sites.copy()

		df_reference['reference_position'] = df_reference.groupby('cut_cluster')['position'].transform('mean').astype(int)

		df_reference = df_reference[['chromosome', 'strand', 'cut_cluster', 'reference_position']].drop_duplicates()

		chromosome_groups = df_reference.groupby(['chromosome'])

		df_sequence = pd.DataFrame()

		for _, df_group in tqdm(chromosome_groups, desc="Extracting cut regions"):

			chromosome_ = df_group['chromosome'].unique().item()

			slices = retrieve_genome_slices_memoryview(sequence = self.genome[chromosome_], positions = df_group['reference_position'], flank_size = self.flank_size)

			df_group['cut_region'] = df_group['reference_position'].map(slices)

			df_sequence = pd.concat([df_sequence, df_group])

		self.df_reference_cut_sites = df_sequence

	def assign_features(self, all_features, gene_names):
		"""
		all_features must have "feature" column
		gene_names must have "gene_id" as column for gene name
		"""
		all_standard_cuts = []
		for standard_cut in self.cut_sites:
			all_standard_cuts.append({'Chromosome': standard_cut.chromosome , 'Start': standard_cut.global_position['cut'], 'End': standard_cut.global_position['cut'], 'cut_cluster': standard_cut.cut_cluster})

		all_standard_cuts = convert_df_to_granges(pd.DataFrame(all_standard_cuts))

		all_features = convert_df_to_granges(all_features)

		gene_names = convert_df_to_granges(gene_names)		

		feature_cut_overlaps = batch_overlaps(gr = all_features, sites_gr = all_standard_cuts)

		feature_cut_overlaps = feature_cut_overlaps.drop(columns = ['Start', 'End', 'name2']).rename(columns = {'Start_b': 'Start', 'End_b': 'End'})

		feature_cut_nearest = batch_nearest_feature(gr = gene_names, sites_gr = all_standard_cuts)

		feature_cut_nearest = feature_cut_nearest.drop(columns = ['Start', 'End']).rename(columns = {'Start_b': 'Start', 'End_b': 'End', 'name2': 'name'})

		for standard_cut in self.cut_sites:
			standard_cut.extract_features(df = feature_cut_overlaps)
			standard_cut.extract_nearest_gene(df = feature_cut_nearest)

	@report_time
	def build_cut_profile(self):
		for standard_cut in self.cut_sites:
			standard_cut.build_cut_profile()


	@report_time
	def build_simple_cut_profile(self):
		for standard_cut in self.cut_sites:
			standard_cut.build_simple_cut_profile()

	@report_time
	def cut_profiles_to_df(self):
		cut_profiles = []
		for standard_cut in self.cut_sites:
			cut_profiles.append(standard_cut.profile)

		self.df_cut_profiles = pd.DataFrame(cut_profiles)


	def build_cut_sites(self):
		"""
		"""
		for _, row in self.df_reference_cut_sites.iterrows():
			
			df_cluster_subset = self.df_cut_sites[self.df_cut_sites['cut_cluster'] == row.cut_cluster]

			standard_cut = CutSite(
				chromosome = row.chromosome,
				strand = row.strand,
				ref_position = row.reference_position,
				cut_region = row.cut_region,
				sgRNA = self.sgRNA,
				sgRNA_alignment_tolerance = self.sgRNA_alignment_tolerance,
				sgRNA_alignment_start_offset = self.sgRNA_alignment_start_offset,
				flank_size = self.flank_size,
				cut_distance = self.cut_distance,
				detail = df_cluster_subset
				)

			self.cut_sites.append(standard_cut)

	@report_time
	def cut_detail_to_df(self, include_cut_location = False):
		"""
		include_cut_location adds self.cut_site from the CutSite object to self.detail
		"""
		cut_detail = []
		for standard_cut in self.cut_sites:
			if not include_cut_location:
				standard_cut.detail['cut'] = standard_cut.cut_site

			cut_detail.append(standard_cut.detail)

		self.df_cut_detail = pd.concat(cut_detail)

	@report_time
	def parallel_build_cut_site_alignment(self, max_workers=None):
		"""
		Parallelize sgRNA alignment across all CutSite objects
		"""
		if max_workers is None:
			max_workers = multiprocessing.cpu_count()

		with ProcessPoolExecutor(max_workers=max_workers) as executor:
			futures = {executor.submit(build_cut_site_alignment_worker, cut_site): cut_site for cut_site in self.cut_sites}

			results = []
			for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning"):
				try:
					result = future.result()
					results.append(result)
				except Exception as e:
					print(f"Error during alignment: {e}")

		self.cut_sites = results

	@report_time
	def single_build_cut_site_alignment(self):
		self.cut_sites = [build_cut_site_alignment_worker(i) for i in self.cut_sites]


def build_cut_site_alignment_worker(cut_site):
	cut_site.find_best_sgRNA_alignment()
	cut_site.calculate_global_positions()
	return cut_site


class CutSite:

	def __init__(self, chromosome, strand, ref_position, sgRNA, sgRNA_alignment_tolerance, sgRNA_alignment_start_offset, cut_distance, detail, flank_size, cut_region = {}):
		self.chromosome = chromosome
		self.strand = strand
		self.ref_position = ref_position

		if cut_region != {}:
			
			self.cut_region = {}
			self.cut_region['start'], self.cut_region['stop'] = sequence_slice_locations(pos = ref_position, flank_size = flank_size)
			
			if self.strand == '-':
				self.cut_region['sequence'] = str(Seq(cut_region).reverse_complement())
			
			else:
				self.cut_region['sequence'] = cut_region

		self.sgRNA = sgRNA
		self.sgRNA_alignment_tolerance = sgRNA_alignment_tolerance
		self.sgRNA_alignment_start_offset = sgRNA_alignment_start_offset

		self.cut_distance = cut_distance

		self.detail = detail

		self.cut_cluster = self.detail.cut_cluster.unique().item()

		self.flank_size = flank_size

		self.alignment = {}
		self.features = {}
		self.global_position = {
			'cut': None
		}

		self.cut_site = None

	def __len__(self):
		return len(self.detail) 

	def __repr__(self):
		return f"CutSite(chrom={self.chromosome}, strand={self.strand}, ref_pos={self.ref_position}, cut={self.cut_site} diversity={len(self)})"

	def find_best_sgRNA_alignment(self):
		"""
		Loop through different alignment parameters starting with the preferred alignment method.
		Start Stop are in reference to the sgRNA alignment 5'->3' on the plus strand 
		"""
		sequence = DNA(self.cut_region['sequence'])

		pam_length = 3

		for seqname, query_seq in self.sgRNA.items():

			self.alignment['sgRNA'] = seqname

			self.alignment['alignment'] = global_pairwise_align_nucleotide(query_seq, sequence)

			self.alignment['local_start'], self.alignment['local_stop'] = parse_global_alignment(self.alignment['alignment'])

			self.alignment['local_stop'] = self.alignment['local_stop'] - self.sgRNA_alignment_start_offset[self.alignment['sgRNA']]

			self.alignment['local_cut'] = self.alignment['local_stop'] - self.cut_distance

			# PAM search refinement

			self.alignment['PAM'] = str(self.alignment['alignment'][0][1])[self.alignment['local_stop'] + 1 : self.alignment['local_stop'] + 1 + pam_length]

			self.alignment['PAM_gaps'] = 0

			if self.alignment['PAM'][1:] != 'GG': # hardcoded for NGG pams

				pam_search_results = find_revised_pam(
					sequence = str(self.alignment['alignment'][0][1]),
					protospacer_start = self.alignment['local_stop'],
					search_direction = 'forward',
					max_pam_gaps_allowed = 2,
					pam_type = 'NGG')

				if not pam_search_results['pam_found']:
					pam_search_results = find_revised_pam(
						sequence = str(self.alignment['alignment'][0][1]),
						protospacer_start = self.alignment['local_stop'],
						search_direction = 'backward',
						max_pam_gaps_allowed = 2,
						pam_type = 'NGG')

				if pam_search_results['pam_found']:

					self.alignment['PAM'] = pam_search_results['revised_pam']

					self.alignment['local_stop'] = pam_search_results['revised_protospacer_start']

					self.alignment['PAM_gaps'] = pam_search_results['pam_n_gap']

			self.alignment['aligned_sequence'] = str(self.alignment['alignment'][0][1])[self.alignment['local_start'] : self.alignment['local_stop'] + 1]

			self.alignment['aligned_gRNA'] = str(self.alignment['alignment'][0][0])[self.alignment['local_start'] : self.alignment['local_stop'] + 1]

			self.alignment['alignment_length'] = len(self.alignment['aligned_sequence'])

			self.alignment['aligned_sequence_gaps'] = self.alignment['aligned_sequence'].count('-')

			self.alignment['aligned_gRNA_gaps'] = self.alignment['aligned_gRNA'].count('-')

			if self.alignment['alignment_length'] in self.sgRNA_alignment_tolerance[seqname]:
				break

	def print_alignment_stats(self):
		for k,v in self.alignment.items():
			print(f'{k}: {v}')


	def calculate_global_positions(self):
		"""
		mixed
		"""
		if self.strand == '-':
			# print('negative')
			self.global_position = {
				# In negative strand, 'start' is at the higher genomic position (cut_region['stop']) minus local_stop
				'protospacer_start': self.cut_region['stop'] - self.alignment['local_stop'] + self.alignment['aligned_sequence_gaps'] - 1,
				'protospacer_stop': self.cut_region['stop'] - self.alignment['local_start'] - 1
			}
			self.global_position['cut'] = self.global_position['protospacer_start'] + self.cut_distance
		else:
			# print('positive')
			self.global_position = {
				'protospacer_stop': self.cut_region['start'] + self.alignment['local_start'],
				'protospacer_start': self.cut_region['start'] + self.alignment['local_stop'] - self.alignment['aligned_sequence_gaps']
				}
			self.global_position['cut'] = self.global_position['protospacer_start'] - self.cut_distance	

		self.cut_site = self.global_position['cut']

	def extract_features(self, df):
		self.features['feature_full'] = df[df['cut_cluster'] == self.cut_cluster]


	def extract_nearest_gene(self, df):
		df = df[df['cut_cluster'] == self.cut_cluster]

		self.features['nearest_gene'] = df.gene.item()

		self.features['nearest_gene_distance'] = df.Distance.item()


	def build_simple_cut_profile(self):
		self.profile = {
			'cut_cluster': self.cut_cluster,
			'chromosome': self.chromosome,
			'strand': self.strand,
			'cut': 0,
			'overlap': len(self.detail)
			}

		# global position
		for feature in ['cut', 'protospacer_start', 'protospacer_stop']:
			self.profile[feature] = self.global_position[feature]

		# local position
		for feature in ['local_start', 'local_stop']:
			self.profile[feature] = self.alignment[feature]

		# alignment
		for feature in ['aligned_sequence', 'aligned_gRNA', 'alignment_length', 'PAM', 'PAM_gaps']:
			self.profile[feature] = self.alignment[feature]

		self.profile['cut_region_sequence'] = self.cut_region['sequence']


	def build_cut_profile(self):

		self.profile = {
			'cut_cluster': self.cut_cluster,
			'chromosome': self.chromosome,
			'strand': self.strand,
			'cut': 0,
			'exon': 0,
			'intron': 0,
			'intergenic': 0,
			'overlap': len(self.detail)
			}

		# genomic feature
		cut_site_features = self.features.get('feature_full')

		if cut_site_features is not None and not cut_site_features.empty:
			present_features = set(cut_site_features['feature'].unique())
			for feature in ('exon', 'intron'):
				self.profile[feature] = int(feature in present_features)

		if sum([self.profile['exon'], self.profile['intron']]) == 0:
			self.profile['intergenic'] = 1

		self.profile['nearest_gene'] = self.features['nearest_gene']

		self.profile['nearest_gene_distance'] = self.features['nearest_gene_distance']

		# measurements
		for measurement in ['zscore', 'min_max', 'local_rank']:
			result = central_tendency(self.detail[measurement])
			self.profile[f'{measurement}_mean'] = result['mean']
			self.profile[f'{measurement}_min'] = result['min']
			self.profile[f'{measurement}_max'] = result['max']
			self.profile[f'{measurement}_median'] = result['50%']

		# global position
		for feature in ['cut', 'protospacer_start', 'protospacer_stop']:
			self.profile[feature] = self.global_position[feature]

		# local position
		for feature in ['local_start', 'local_stop']:
			self.profile[feature] = self.alignment[feature]

		# alignment
		for feature in ['aligned_sequence', 'aligned_gRNA', 'alignment_length', 'PAM', 'PAM_gaps']:
			self.profile[feature] = self.alignment[feature]

		self.profile['cut_region_sequence'] = self.cut_region['sequence']

		

	# def score_cut(self):
	# 	"""
	# 	"""
	# 	{'intron': [1 if self.features['genomic_full']['feature']]
	# 	}
	# 	if 'intron':






