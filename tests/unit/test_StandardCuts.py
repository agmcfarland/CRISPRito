from unittest.mock import patch
import pytest
import time
import pandas as pd
from os.path import join as pjoin
from CRISPRito.StandardCuts import StandardCuts

@pytest.fixture
def load_1_group_samplesheet_ptprc(project_test_data_directory):

	df_gr1 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '1_group_samplesheet_ptprc.csv'))

	df_gr1['standard_format_file_path'] = df_gr1['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr1

@pytest.fixture
def load_2_group_samplesheet_ptprc(project_test_data_directory):

	df_gr2 = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites',  '2_group_samplesheet_ptprc.csv'))

	df_gr2['standard_format_file_path'] = df_gr2['standard_format_file_path'].apply(lambda x: pjoin(project_test_data_directory, x))

	return df_gr2


def test_load_1_group_samplesheet_ptprc(load_1_group_samplesheet_ptprc):

	df_gr1 = load_1_group_samplesheet_ptprc

	assert True

def test_load_cut_sites_group1(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_load_cut_sites_group1
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	# print('\n')

	# print(standard_group.df_cut_sites)

	standard_group.update_cut_cluster_id()

	# print(standard_group.df_cut_sites)

	assert len(standard_group.df_cut_sites) == 19


def test_load_cut_sites_group2(load_2_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_load_cut_sites_group2
	"""

	standard_group = StandardCuts(sample_sheet = load_2_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	# print('\n')

	# print(standard_group.df_cut_sites)

	standard_group.update_cut_cluster_id()

	# print(standard_group.df_cut_sites)

	assert len(standard_group.df_cut_sites) == 5


def test_standarize_scores(load_1_group_samplesheet_ptprc):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_standarize_scores
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.standardize_scores()

	assert 'zscore' in standard_group.df_cut_sites.columns

	assert 'min_max' in standard_group.df_cut_sites.columns


def test_load_genome_and_get_size(load_1_group_samplesheet_ptprc, sample_fasta_gz):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_load_genome_and_get_size
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_genome(sample_fasta_gz)
	standard_group.get_genome_size()

	print(standard_group.genome_size)



def skip_test_extract_cut_region_group1(load_1_group_samplesheet_ptprc, path_to_hg38_genome):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_extract_cut_region_group1
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.load_genome(genome_path = path_to_hg38_genome)

	standard_group.extract_cut_region()

	# print(standard_group.df_reference_cut_sites)

	assert len(standard_group.df_reference_cut_sites) == 15

	# standard_group.df_reference_cut_sites.to_csv('/data/CRISPRito/1_group_cluster_regions_ptprc.csv', index = None)

def skip_test_extract_cut_region_group2(load_2_group_samplesheet_ptprc, path_to_hg38_genome):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_extract_cut_region_group2
	"""

	standard_group = StandardCuts(sample_sheet = load_2_group_samplesheet_ptprc)	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.load_genome(genome_path = path_to_hg38_genome)

	standard_group.extract_cut_region()

	# print(standard_group.df_reference_cut_sites)

	# standard_group.df_reference_cut_sites.to_csv('/data/CRISPRito/2_group_cluster_regions_ptprc.csv', index = None)

	assert len(standard_group.df_reference_cut_sites) == 5


def test_get_genome_size(load_1_group_samplesheet_ptprc, sample_fasta_gz):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_get_genome_size
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc, sgRNA = 'AAAATATGCAAACATCACTG')	
	standard_group.load_genome(sample_fasta_gz)
	standard_group.get_genome_size()
	assert standard_group.genome_size == {'chr1': 64, 'chr2': 64}

def test_build_cut_sites(load_1_group_samplesheet_ptprc, project_test_data_directory, example_hg38_genome_size):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_build_cut_sites
	"""

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc, sgRNA = 'AAAATATGCAAACATCACTG')	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.standardize_scores()

	# standard_group.load_genome(genome_path = path_to_hg38_genome)

	# standard_group.extract_cut_region()
	with patch.object(standard_group, 'extract_cut_region', return_value = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites', '1_group_cluster_regions_ptprc.csv'))):
		standard_group.df_reference_cut_sites = standard_group.extract_cut_region()	

	with patch.object(standard_group, 'get_genome_size', return_value = example_hg38_genome_size):
		standard_group.genome_size = standard_group.get_genome_size()	

	print(standard_group.genome_size)

	# print(standard_group.df_reference_cut_sites)

	# print(standard_group.df_cut_sites)

	standard_group.build_cut_sites()

	for standard_cut in standard_group.cut_sites:
		print(standard_cut)

	assert len(standard_group.cut_sites) == 15

@pytest.fixture
def standard_group_1_ptprc_cut_sites(load_1_group_samplesheet_ptprc, project_test_data_directory, example_hg38_genome_size):

	standard_group = StandardCuts(sample_sheet = load_1_group_samplesheet_ptprc, sgRNA = 'AAAATATGCAAACATCACTG')	

	standard_group.load_cut_sites()

	standard_group.cluster_cut_sites()

	standard_group.update_cut_cluster_id()

	standard_group.standardize_scores()

	# standard_group.load_genome(genome_path = path_to_hg38_genome)

	# standard_group.extract_cut_region()
	with patch.object(standard_group, 'extract_cut_region', return_value = pd.read_csv(pjoin(project_test_data_directory, 'cluster_sites', '1_group_cluster_regions_ptprc.csv'))):
		standard_group.df_reference_cut_sites = standard_group.extract_cut_region()	

	with patch.object(standard_group, 'get_genome_size', return_value = example_hg38_genome_size):
		standard_group.genome_size = standard_group.get_genome_size()	

	standard_group.build_cut_sites()

	return standard_group


def test_cut_site_alignment(standard_group_1_ptprc_cut_sites):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_cut_site_alignment
	Using PAM sequence to identify correct since this relies on local_start and local_stop being correct
	"""

	expected =[
		'AGG',
		'TGG',
		'TGG',
		'CGG',
		'TGG',
		'AGG',
		'TGG',
		'GTG',
		'AGG',
		'AGT',
		'AGC',
		'AGA',
		'AGG',
		'GGG',
		'TGG']

	standard_group = standard_group_1_ptprc_cut_sites

	for e, standard_cut in enumerate(standard_group.cut_sites):
		# print('\n')

		# print(standard_cut)

		standard_cut.find_best_sgRNA_alignment()

		assert standard_cut.alignment['PAM'] == expected[e]

		# standard_cut.print_alignment_stats()
		# break


def check_cut_sites_for_pam(standard_group, expected):

	for e, standard_cut in enumerate(standard_group.cut_sites):

		assert standard_cut.alignment['PAM'] == expected[str(standard_cut)]


def test_multiprocessing_cut_site_alignment(standard_group_1_ptprc_cut_sites):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_multiprocessing_cut_site_alignment
	Using PAM sequence to identify correct since this relies on local_start and local_stop being correct
	"""

	expected = {
	'CutSite(chrom=chr1, strand=+, ref_pos=183741771, cut=183741786 diversity=1)': 'AGG',
	'CutSite(chrom=chr1, strand=+, ref_pos=198706743, cut=198706753 diversity=3)': 'TGG',
	'CutSite(chrom=chr14, strand=+, ref_pos=40300864, cut=40300879 diversity=1)': 'TGG',
	'CutSite(chrom=chr18, strand=+, ref_pos=57630524, cut=57630539 diversity=1)': 'CGG',
	'CutSite(chrom=chr18, strand=-, ref_pos=41108133, cut=41108140 diversity=1)': 'TGG',
	'CutSite(chrom=chr2, strand=-, ref_pos=143961591, cut=143961596 diversity=3)': 'AGG',
	'CutSite(chrom=chr2, strand=-, ref_pos=184581387, cut=184581394 diversity=1)': 'TGG',
	'CutSite(chrom=chr3, strand=+, ref_pos=65866797, cut=65866814 diversity=1)': 'GTG',
	'CutSite(chrom=chr3, strand=+, ref_pos=138494328, cut=138494327 diversity=1)': 'AGG',
	'CutSite(chrom=chr6, strand=-, ref_pos=100224884, cut=100224884 diversity=1)': 'AGT',
	'CutSite(chrom=chr6, strand=-, ref_pos=134850663, cut=134850668 diversity=1)': 'AGC',
	'CutSite(chrom=chr7, strand=-, ref_pos=28169223, cut=28169229 diversity=1)': 'AGA',
	'CutSite(chrom=chr7, strand=-, ref_pos=115239484, cut=115239490 diversity=1)': 'AGG',
	'CutSite(chrom=chr8, strand=+, ref_pos=6301238, cut=6301255 diversity=1)': 'GGG',
	'CutSite(chrom=chr8, strand=+, ref_pos=28930554, cut=28930553 diversity=1)': 'TGG'}

	print('\n')
	print('Multithreaded')
	standard_group = standard_group_1_ptprc_cut_sites
	start = time.time()
	standard_group.multithread_sgRNA_alignment()
	print(time.time()-start)
	check_cut_sites_for_pam(standard_group, expected)

	print('Multiprocessor')
	standard_group = standard_group_1_ptprc_cut_sites
	start = time.time()
	standard_group.parallel_sgRNA_alignment()
	print(time.time()-start)
	check_cut_sites_for_pam(standard_group, expected)

	print('Single threaded')
	standard_group = standard_group_1_ptprc_cut_sites
	start = time.time()
	standard_group.single_sgRNA_alignment()
	print(time.time()-start)
	check_cut_sites_for_pam(standard_group, expected)



def test_calculate_global_positions(standard_group_1_ptprc_cut_sites):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_calculate_global_positions
	"""

	standard_group = standard_group_1_ptprc_cut_sites

	print('\n')
	for e, standard_cut in enumerate(standard_group.cut_sites):
		
		# print('\n')

		# print(standard_cut)

		standard_cut.find_best_sgRNA_alignment()

		standard_cut.calculate_global_positions()

		if standard_cut.strand == '-':
			assert standard_cut.global_position['protospacer_stop'] > standard_cut.global_position['protospacer_start']
		else:
			assert standard_cut.global_position['protospacer_start'] > standard_cut.global_position['protospacer_stop']
		
		print(standard_cut.global_position)


def test_extract_in_refseq_feature(standard_group_1_ptprc_cut_sites, path_to_hg38_refseq):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_extract_in_refseq_feature
	"""
	
	standard_group = standard_group_1_ptprc_cut_sites

	df_genomic_positions = pd.read_csv(path_to_hg38_refseq)

	expected_results = {
	'CutSite(chrom=chr1, strand=+, ref_pos=183741771, cut=183741786 diversity=1)': {'genomic_full': 4,
	'genomic_summary':2,
	'nearest_gene':'RGL1',
	'nearest_gene_distance':0.0},


	'CutSite(chrom=chr1, strand=+, ref_pos=198706743, cut=198706753 diversity=3)': {'genomic_full': 2,
	'genomic_summary':1,
	'nearest_gene':'PTPRC',
	'nearest_gene_distance':0.0},


	'CutSite(chrom=chr14, strand=+, ref_pos=40300864, cut=40300879 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'FBXO33',
	'nearest_gene_distance':868446.0},


	'CutSite(chrom=chr18, strand=+, ref_pos=57630524, cut=57630539 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'NARS1',
	'nearest_gene_distance':8703.0},


	'CutSite(chrom=chr18, strand=-, ref_pos=41108133, cut=41108140 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'PIK3C3',
	'nearest_gene_distance':847093.0},


	'CutSite(chrom=chr2, strand=-, ref_pos=143961591, cut=143961596 diversity=3)': {'genomic_full': 27,
	'genomic_summary':1,
	'nearest_gene':'GTDC1',
	'nearest_gene_distance':0.0},


	'CutSite(chrom=chr2, strand=-, ref_pos=184581387, cut=184581394 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'ZNF804A',
	'nearest_gene_distance':17134.0},


	'CutSite(chrom=chr3, strand=+, ref_pos=65866797, cut=65866814 diversity=1)': {'genomic_full': 3,
	'genomic_summary':1,
	'nearest_gene':'MAGI1',
	'nearest_gene_distance':0.0},


	'CutSite(chrom=chr3, strand=+, ref_pos=138494328, cut=138494327 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'CEP70',
	'nearest_gene_distance':16.0},


	'CutSite(chrom=chr6, strand=-, ref_pos=100224884, cut=100224884 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'SIM1',
	'nearest_gene_distance':160124.0},


	'CutSite(chrom=chr6, strand=-, ref_pos=134850663, cut=134850668 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'ALDH8A1',
	'nearest_gene_distance':66724.0},


	'CutSite(chrom=chr7, strand=-, ref_pos=28169223, cut=28169229 diversity=1)': {'genomic_full': 1,
	'genomic_summary':1,
	'nearest_gene':'JAZF1',
	'nearest_gene_distance':0.0},


	'CutSite(chrom=chr7, strand=-, ref_pos=115239484, cut=115239490 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'MDFIC',
	'nearest_gene_distance':219573.0},


	'CutSite(chrom=chr8, strand=+, ref_pos=6301238, cut=6301255 diversity=1)': {'genomic_full': 0,
	'genomic_summary':0,
	'nearest_gene':'MCPH1',
	'nearest_gene_distance':105371.0},


	'CutSite(chrom=chr8, strand=+, ref_pos=28930554, cut=28930553 diversity=1)': {'genomic_full': 16,
	'genomic_summary':2,
	'nearest_gene':'HMBOX1',
	'nearest_gene_distance':0.0}
	}

	def check_cut_site_correct(standard_group, expected_results):
		
		for standard_cut in standard_group.cut_sites:

			for k, v in standard_cut.features.items():
				if type(v) == pd.DataFrame:
					assert len(v) == expected_results[str(standard_cut)][k]
				else:
					assert v == expected_results[str(standard_cut)][k]

				assert standard_cut.alignment['PAM'] == expected[str(standard_cut)]


	print('\n')
	print('Multithreaded')
	standard_group = standard_group_1_ptprc_cut_sites
	start = time.time()
	standard_group.multithread_sgRNA_alignment()
	print(time.time()-start)
	check_cut_site_correct(standard_group, expected_results)

	print('Multiprocessor')
	standard_group = standard_group_1_ptprc_cut_sites
	start = time.time()
	standard_group.parallel_sgRNA_alignment()
	print(time.time()-start)
	check_cut_site_correct(standard_group, expected_results)

	print('Single threaded')
	standard_group = standard_group_1_ptprc_cut_sites
	start = time.time()
	standard_group.single_sgRNA_alignment()
	print(time.time()-start)
	check_cut_site_correct(standard_group, expected_results)


def test_multiprocessing_extract_in_refseq_feature(standard_group_1_ptprc_cut_sites, path_to_hg38_refseq):
	"""
	pytest -sv tests/unit/test_StandardCuts.py::test_extract_in_refseq_feature
	"""


	pass












