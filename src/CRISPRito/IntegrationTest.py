import sys
import pandas as pd
from os.path import join as pjoin
import argparse
import os
import subprocess


def pull_data(test_dir, data_dir, server_path):
	"""
	Performs all operations needed for the base inputs to cluster_cuts and rank_sites
	"""
	if os.path.exists(test_dir):
		sys.exit(f'Error.\nTest dir: "{test_dir}" already exists')
	
	else:
		print(f'Making test dir at: {test_dir}')
		os.makedirs(test_dir)

	# # Download data
	# os.system(f'wget https://microb191.med.upenn.edu/crisprito_integration_test.tar.gz -P {test_dir}')
	os.system(f'wget {server_path} -P {test_dir}')

	# stand in for download data
	# os.system(
	# 	'ln -s /data/test_crisprito/crisprito_integration_test.tar.gz /data/test_crisprito/setup_test/crisprito_integration_test.tar.gz'
	# 	)
	# # Untar
	os.system(f'tar -xzf {pjoin(test_dir, "crisprito_integration_test.tar.gz")} -C {test_dir}')

	# Set paths in drivers
	# feature
	df_feature = pd.read_csv(pjoin(data_dir, 'crisprito_feature_driver.csv'))
	df_feature['file_path'] = df_feature['file_path'].apply(lambda x: pjoin(data_dir, os.path.basename(x)))
	df_feature.to_csv(pjoin(test_dir, 'feature_driver.csv'), index = None)

	# samplesheet
	df_samplesheet = pd.read_csv(pjoin(data_dir, 'example_samplesheet_driver.csv'))
	df_samplesheet['standard_format_file_path'] = df_samplesheet['standard_format_file_path'].apply(lambda x: pjoin(data_dir, os.path.basename(x)))
	df_samplesheet = df_samplesheet[['sample', 'method', 'cluster_group', 'standard_format_file_path']]
	df_samplesheet.to_csv(pjoin(test_dir, 'samplesheet_driver.csv'), index = None)

	assert df_samplesheet.shape == (27,4)

	assert df_feature.shape == (2,3)

	assert os.path.exists(pjoin(data_dir, 'hg38.fasta.gz'))

	print('-----------------------------\nCRISPRito test data loaded successfully!\n-----------------------------')	
	



def main():
	parser = argparse.ArgumentParser(description="CRISPRito integration test")
	parser.add_argument("--test_dir", default = None, required = True, help = "Directory for testing inputs and outputs")
	parser.add_argument(
		"--test_workflow", 
		default = None, 
		choices = ['pull_data', 'cluster_cuts', 'auto_rank_sites'], 
		required = True,
		help = "Workflow to test"
		)
	parser.add_argument(
		'--server_path',
		default = 'https://microb191.med.upenn.edu/crisprito_integration_test.tar.gz',
		help = 'Path to test data host.'
		)

	args = parser.parse_args()

	test_dir = os.path.abspath(args.test_dir)

	# test_dir = '/data/test_crisprito/setup_test' # temp hardcode

	data_dir = pjoin(test_dir, 'crisprito_integration_test_data')

	crisprito_output_dir = pjoin(test_dir, 'crisprito_output')

	if args.test_workflow == 'pull_data':
		pull_data(test_dir = test_dir, data_dir = data_dir, server_path = args.server_path)

	elif args.test_workflow == 'cluster_cuts':

		## Run first step
		subprocess.run(
			[
			'crisprito-cluster-cuts',
			'--sample_sheet_path', pjoin(test_dir, 'samplesheet_driver.csv'),
			'--output_dir', crisprito_output_dir,
			'--genome_path', pjoin(data_dir, 'hg38.fasta.gz'), 
			'--feature_table_path', pjoin(test_dir, 'feature_driver.csv'),
			'--overwrite_output_dir',
			'--flank_size', '40',
			'--sgRNA', 'AAAATATGCAAACATCACTG',
			'--PAM_alignment', 'NGG',
			'--range_threshold', '60',
			'--workflow', 'cluster_cuts'
			]
			)

		expected_outputs = [
			'1_group_samplesheet.csv',
			'CRISPRito.log',
			'1_group_cut_profiles.csv',
			'1_group_method_counts.csv',
			'1_group_id_counts.csv',
			'1_group_id_cut_detail.csv',
			'1_group_id_rank_weight_skeleton.csv'
			]

		for file_ in expected_outputs:
			assert os.path.exists(pjoin(crisprito_output_dir, file_))

		df_cut_profiles = pd.read_csv(pjoin(crisprito_output_dir, '1_group_cut_profiles.csv'))

		assert df_cut_profiles.shape == (14352, 32)

		df_cut_detail = pd.read_csv(pjoin(crisprito_output_dir, '1_group_id_cut_detail.csv'))

		assert df_cut_detail.shape == (16451, 9)

		print('-----------------------------\ncrisprito-cluster-cuts completed successfully!\n-----------------------------')


	## Run second step
	elif args.test_workflow == 'auto_rank_sites':
		subprocess.run(
			[
			'crisprito-rank-sites',
			'--workflow', 'auto', 
			'--data_dir', crisprito_output_dir,
			'--cluster_group', '1',
			'--feature_driver_path', pjoin(test_dir, 'feature_driver.csv'),
			'--magnitude_transform', 'percentile',
			'--magnitude_aggregation', 'max',
			'--tau', '5000',
			'--gamma', '0.5',
			'--max_distance_bp', '10000',
			'--spatial_aggregation', 'weighted_mean',
			'--score_col', 'rra_weighted_magnitude_score',
			'--output_dir', crisprito_output_dir,
			'--output_name', 'ranked_cut_sites'
			]
			)

		expected_outputs = [
			'ranked_cut_sites.csv',
			'ranked_cut_sites.log'
			]

		for file_ in expected_outputs:
			assert os.path.exists(pjoin(crisprito_output_dir, file_))

		df_ranked = pd.read_csv(pjoin(crisprito_output_dir, 'ranked_cut_sites.csv'))

		assert df_ranked.shape == (14352, 24)

		print('-----------------------------\ncrisprito-rank-sites completed successfully!\n-----------------------------')

	else:
		sys.exit('No useable workflow specified')


if __name__ == '__main__':
	main()