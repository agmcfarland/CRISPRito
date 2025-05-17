import argparse
import os
from os.path import join as pjoin
import pandas as pd
from glob import glob
from CRISPRito.ProcessGroup import process_group
from CRISPRito.SetupRun import setup_run


def main():
	parser = argparse.ArgumentParser(description="Run full CRISPRito cut \nstandardization and annotation\n pipline")
	parser.add_argument("--sample_sheet_path", required=True)
	parser.add_argument("--output_dir", required=True)
	parser.add_argument("--genome_path", required=True, default = None)
	parser.add_argument("--feature_path", required=False, default = None)
	parser.add_argument("--gene_names_path", required=False, default = None)
	parser.add_argument("--overwrite_output_dir", action="store_true")
	parser.add_argument("--flank_size", type=int, default=30)
	parser.add_argument("--sgRNA", required=True)
	parser.add_argument("--PAM_alignment", default='-GG', required=True)
	parser.add_argument("--range_threshold", type=int, default=20, help="Distance between clusters")
	parser.add_argument("--workflow", type = str, help = 'One of [cluster_cuts, detect_cut]')

	args = parser.parse_args()

	# Step 1: Setup
	setup_run(
		sample_sheet_path=args.sample_sheet_path,
		output_dir=args.output_dir,
		genome_path=args.genome_path,
		feature_path=args.feature_path,
		gene_names_path=args.gene_names_path,
		overwrite_output_dir=args.overwrite_output_dir
	)

	# # Step 2: Process each *_group_samplesheet.csv
	group_files = glob(pjoin(args.output_dir, "*_group_samplesheet.csv"))

	print(group_files)

	# if not group_files:
	# 	raise RuntimeError("No *_group_samplesheet.csv files found after setup.")

	for group_file in group_files:

		if args.workflow == 'detect_cut':
			range_to_site(
				group_samplesheet_path=group_file,
				output_path=args.output_dir,
				genome_path=args.genome_path,
				feature_path=args.feature_path,
				gene_names_path=args.gene_names_path,
				flank_size=args.flank_size,
				sgRNA=args.sgRNA,
				PAM_alignment=args.PAM_alignment,
				range_threshold = args.range_threshold
				)

		if args.workflow == 'cluster_cuts':
			process_group(
				group_samplesheet_path=group_file,
				output_path=args.output_dir,
				genome_path=args.genome_path,
				feature_path=args.feature_path,
				gene_names_path=args.gene_names_path,
				flank_size=args.flank_size,
				sgRNA=args.sgRNA,
				PAM_alignment=args.PAM_alignment,
				range_threshold = args.range_threshold
			)

	print(f'Output\n{args.output_dir}')


if __name__ == "__main__":
	main()