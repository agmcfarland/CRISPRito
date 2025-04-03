import pandas as pd
from os.path import join as pjoin
from CRISPRito.Utils import greedy_clustering_incremental, retrieve_genome_slices_memoryview, genome_to_dict_memoryview

class StandardCuts:

	def __init__(self, sample_sheet, flank_size = 30):
		self.sample_sheet = sample_sheet
		self.flank_size = flank_size
		pass

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
		
		for group_name, df_group in grouped_cuts:

			cluster_groups = greedy_clustering_incremental(df_group['position'].tolist(), range_threshold = range_threshold)

			position_to_group = {pos: group for group, positions in cluster_groups.items() for pos in positions}

			df_group['cut_cluster'] = df_group['position'].map(position_to_group)

			df_clustered_cuts = pd.concat([df_clustered_cuts, df_group])


		if len(self.df_cut_sites) != len(df_clustered_cuts):

			raise ValueError('Not all cut sites had a cut group assigned')

		self.df_cut_sites = df_clustered_cuts 

	def load_genome(self, genome_path):
		self.genome = genome_to_dict_memoryview(genome_path)

	def extract_cut_region(self):
		"""
		Iterate through list of positions, grouped by chromosome, and extract self.flank size bp from the mean predicted position 
		"""

		df_reference = self.df_cut_sites.copy()

		df_reference['reference_position'] = df_reference.groupby('cut_cluster')['position'].transform('mean').astype(int)

		df_reference = df_reference[['chromosome', 'strand', 'cut_cluster', 'reference_position']].drop_duplicates()

		chromosome_groups = df_reference.groupby(['chromosome'])

		df_sequence = pd.DataFrame()

		for group_name, df_group in chromosome_groups:

			chromosome_ = df_group['chromosome'].unique().item()

			slices = retrieve_genome_slices_memoryview(sequence = self.genome[chromosome_], positions = df_group['reference_position'], flank_size = self.flank_size)

			df_group['cut_region'] = df_group['reference_position'].map(slices)

			df_sequence = pd.concat([df_sequence, df_group])

		self.df_reference_cut_sites = df_sequence


	def build_cut_sites(self):
		"""
		"""





	# def extract_cut_region():
		# pass

		


class CutSite:

	def __init__(self, chromosome, strand, sequence, sgRNA, ):
		self.chromosome = chromosome
		self.strand = strand
		self.sequence = sequence
		self.sgRNA = sgRNA


	# def 






