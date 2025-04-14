from itertools import repeat
import pandas as pd
import os
from os.path import join as pjoin
from Bio.Seq import Seq
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
import multiprocessing
from tqdm import tqdm
from CRISPRito.Utils import (
	greedy_clustering_incremental,
	retrieve_genome_slices_memoryview,
	genome_to_dict_memoryview,
	scale_zscore,
	scale_min_max,
	parse_global_alignment,
	sequence_slice_locations,
	extract_annotations,
	get_closest_annotation,
	slice_annotation
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
			'fwd_NGG' : range(20,25)
		}

		self.sgRNA_alignment_start_offset = {
			'fwd' : 0,
			'fwd_NGG' : 3
		}
		self.cut_distance = 3

		self.df_cut_sites = pd.DataFrame()

		self.cut_sites = []

	def __repr__(self):
		return f'StandardCuts\n{self.sample_sheet}\n{self.df_cut_sites}'

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


	def standardize_cuts(self, genome_path, df_genomic_features):

		self.load_cut_sites()

		self.cluster_cut_sites()

		self.update_cut_cluster_id()

		self.remove_cluster_duplicates()

		self.standardize_scores()

		self.load_genome(genome_path = genome_path)

		self.extract_cut_region()

		self.build_cut_sites()

		self.parallel_build_cut_site_alignment()

		self.multithread_build_cut_site_annotation(df_genomic_features = df_genomic_features)

		cut_profile = []
		for standard_cut in self.cut_sites:
			cut_profile.append(standard_cut)

		self.df_cut_profile = pd.DataFrame(cut_profile)


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


	def multithread_build_cut_site_annotation(self, df_genomic_features, max_workers = None):
		if max_workers is None:
			max_workers = multiprocessing.cpu_count()

		with ThreadPoolExecutor(max_workers= max_workers) as executor:
			# self.cut_sites = list(tqdm(executor.map(build_cut_site_profile_worker, self.cut_sites), total=len(self.cut_sites)))
			self.cut_sites = list(executor.map(
				build_cut_site_annotation_worker,
				self.cut_sites,
				repeat(df_genomic_features)
				)
			)

	def parallel_build_cut_site_alignment(self, max_workers=None):
		"""
		Parallelize sgRNA alignment across all CutSite objects
		"""
		if max_workers is None:
			max_workers = multiprocessing.cpu_count()

		with ProcessPoolExecutor(max_workers=max_workers) as executor:
			futures = {executor.submit(build_cut_site_alignment_worker, cut_site): cut_site for cut_site in self.cut_sites}

			results = []
			for future in as_completed(futures):
				try:
					result = future.result()
					results.append(result)
				except Exception as e:
					print(f"Error during alignment: {e}")

		self.cut_sites = results

	def single_build_cut_site_alignment(self):
		self.cut_sites = [build_cut_site_alignment_worker(i) for i in self.cut_sites]

	# def collect_cut_sit


def build_cut_site_alignment_worker(cut_site):
	cut_site.find_best_sgRNA_alignment()
	cut_site.calculate_global_positions()
	return cut_site


def build_cut_site_annotation_worker(cut_site, df_genomic_features):
	cut_site.identify_genomic_features(df = df_genomic_features)
	cut_site.build_local_score()
	return cut_site


class CutSite:

	def __init__(self, chromosome, strand, ref_position, cut_region, sgRNA, sgRNA_alignment_tolerance, sgRNA_alignment_start_offset, cut_distance, detail, flank_size):
		self.chromosome = chromosome
		self.strand = strand
		self.ref_position = ref_position
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

		self.flank_size = flank_size

		self.alignment = {}
		self.features = {}
		self.global_position = {
			'cut': None
		}

	def __len__(self):
		return len(self.detail) 

	def __repr__(self):
		return f"CutSite(chrom={self.chromosome}, strand={self.strand}, ref_pos={self.ref_position}, cut={self.global_position['cut']} diversity={len(self)})"

	def find_best_sgRNA_alignment(self):
		"""
		Loop through different alignment parameters starting with the preferred alignment method.
		Start Stop are in reference to the sgRNA alignment 5'->3' on the plus strand 
		"""
		sequence = DNA(self.cut_region['sequence'])

		for seqname, query_seq in self.sgRNA.items():
			self.alignment['sgRNA'] = seqname

			self.alignment['alignment'] = global_pairwise_align_nucleotide(query_seq, sequence)

			self.alignment['local_start'], self.alignment['local_stop'] = parse_global_alignment(self.alignment['alignment'])

			self.alignment['local_stop'] = self.alignment['local_stop'] - self.sgRNA_alignment_start_offset[self.alignment['sgRNA']]

			self.alignment['aligned_sequence'] = self.cut_region['sequence'][self.alignment['local_start'] : self.alignment['local_stop'] + 1]

			self.alignment['aligned_gRNA'] = str(self.alignment['alignment'][0][0])[self.alignment['local_start'] : self.alignment['local_stop'] + 1]

			self.alignment['alignment_length'] = len(self.alignment['aligned_sequence'])

			self.alignment['PAM'] = self.cut_region['sequence'][self.alignment['local_stop'] + 1 : self.alignment['local_stop'] + 1 + 3]

			self.alignment['local_cut'] = self.alignment['local_stop'] - self.cut_distance

			if self.alignment['alignment_length'] in self.sgRNA_alignment_tolerance[seqname]:
				break

	def print_alignment_stats(self):
		for k,v in self.alignment.items():
			print(f'{k}: {v}')

	def calculate_global_positions(self):
		if self.strand == '-':
			self.global_position = {
				'protospacer_stop': self.cut_region['stop'] - self.alignment['local_start'],
				'protospacer_start': self.cut_region['stop'] - self.alignment['local_stop']
				}
			self.global_position['cut'] = self.global_position['protospacer_start'] + self.cut_distance	
		else:
			self.global_position = {
				'protospacer_stop': self.cut_region['start'] + self.alignment['local_start'],
				'protospacer_start': self.cut_region['start'] + self.alignment['local_stop']
				}
			self.global_position['cut'] = self.global_position['protospacer_start'] - self.cut_distance	

	def identify_genomic_features(self, df, tolerance = 3_000_000):

		df = slice_annotation(df, chromosome = self.chromosome, position = self.global_position['cut'])

		self.features['genomic_full'] = extract_annotations(df = df, position = self.global_position['cut'])

		self.features['genomic_summary'] = self.features['genomic_full'].drop_duplicates(subset=["name2", "feature"])

		self.features['nearest_gene'], self.features['nearest_gene_distance'] = get_closest_annotation(df = df, position = self.global_position['cut'], column_name = 'name2')


	def build_local_score(self):
		cut_site_features = self.features.get('genomic_full')

		self.local_score = {
			'exon': 0,
			'intron': 0,
			'intergenic': 0,
			'overlap': len(self.detail),
			'mean_z_score': self.detail['zscore'].mean().astype(float),
			'mean_min_max_score': self.detail['min_max'].mean().astype(float),
			'mean_ordered_rank': self.detail['local_rank'].mean().astype(float)
		}

		if cut_site_features is not None and not cut_site_features.empty:
			present_features = set(cut_site_features['feature'].unique())
			for feature in ('exon', 'intron'):
				self.local_score[feature] = int(feature in present_features)

		if sum([self.local_score['exon'], self.local_score['intron']]) == 0:
			self.local_score['intergenic'] = 1

		


	# def score_cut(self):
	# 	"""
	# 	"""
	# 	{'intron': [1 if self.features['genomic_full']['feature']]
	# 	}
	# 	if 'intron':






