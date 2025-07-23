import argparse
import os
from os.path import join as pjoin
import pandas as pd
import traceback
from glob import glob
from CRISPRito.Utils import report_time, setup_logging
from CRISPRito.ProcessGroup import process_group
from CRISPRito.RangeToSite import range_to_site
from CRISPRito.SetupRun import setup_run

@report_time
def main():
	parser = argparse.ArgumentParser(description="Run CRISPRito worfklows")
	parser.add_argument("--sample_sheet_path", required=True)
	parser.add_argument("--output_dir", required=True)
	parser.add_argument("--genome_path", required=True, default = None)
	parser.add_argument("--feature_table_path", required=False)
	parser.add_argument("--overwrite_output_dir", action="store_true")
	parser.add_argument("--flank_size", type=int, default=30)
	parser.add_argument("--sgRNA", required=True)
	parser.add_argument("--PAM_alignment", default='-GG', required=True)
	parser.add_argument("--range_threshold", type=int, default=20, help="Distance between clusters")
	parser.add_argument("--workflow", type = str, help = 'One of [cluster_cuts, detect_cut]')

	args = parser.parse_args()

	log_path = pjoin(args.output_dir, "CRISPRito.log")
	
	setup_logging(log_path)

	try:
		print("Starting CRISPRito...")

		# Step 1: Setup
		setup_run(
			sample_sheet_path=args.sample_sheet_path,
			output_dir=args.output_dir,
			genome_path=args.genome_path,
			feature_table_path=args.feature_table_path,
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
					feature_table_path=args.feature_table_path,
					flank_size=args.flank_size,
					sgRNA=args.sgRNA,
					PAM_alignment=args.PAM_alignment,
					range_threshold = args.range_threshold
				)

		print(f'Output\n{args.output_dir}')
		print('Finished CRISPRito.')
	except Exception as e:
		logging.error("Unhandled exception occurred:")
		logging.error(traceback.format_exc())
		sys.exit(1)


if __name__ == "__main__":
	main()