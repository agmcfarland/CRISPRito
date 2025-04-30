from CRISPRito.SampleManager import SampleManager
from CRISPRito.RunParameters import RunParameters
import argparse
import os

def setup_run(
	sample_sheet_path,
	output_dir,
	genome_path,
	feature_path,
	gene_names_path,
	overwrite_output_dir):
	"""
	"""
	run_parameters = RunParameters(
		output_dir = output_dir
		)

	run_parameters.manage_output_dir(
		overwrite = overwrite_output_dir
		)

	run_parameters.check_inputs_exist(
		sample_sheet_path = sample_sheet_path,
		genome_path = genome_path,
		feature_path = feature_path,
		gene_names_path = gene_names_path
		)

	manager = SampleManager(
		sample_sheet_path = sample_sheet_path,
		output_dir = output_dir
		)

	manager.set_up()

	print(os.listdir(output_dir))

def main():
	parser = argparse.ArgumentParser(description="Set up a CRISPRito run.")

	parser.add_argument("--sample_sheet_path", help="Path to the group sample sheet CSV.")
	parser.add_argument("--output_dir", help="Directory to store run outputs.")
	parser.add_argument("--genome_path", help="Path to the reference genome file.")
	parser.add_argument("--feature_path", help="Path to the genomic features CSV file.")
	parser.add_argument("--gene_names_path", help="Path to the gene names CSV file.")
	parser.add_argument("--overwrite_output_dir", action="store_true", help="Overwrite output directory if it exists.")

	args = parser.parse_args()

	print(args)

	setup_run(
		sample_sheet_path=args.sample_sheet_path,
		output_dir=args.output_dir,
		genome_path=args.genome_path,
		feature_path=args.feature_path,
		gene_names_path=args.gene_names_path,
		overwrite_output_dir=args.overwrite_output_dir
	)



if __name__ == '__main__':
	main()
