import os
import shutil
import sys

class RunParameters:
	def __init__(self, output_dir):
		"""
		Initializes RunParameters and manages the output directory.

		:param output_dir: Path to the output directory.
		"""
		self.output_dir = output_dir

	def manage_output_dir(self, overwrite):
		"""Creates or resets the output directory based on overwrite setting."""
		if os.path.exists(self.output_dir):
			if overwrite:
				shutil.rmtree(self.output_dir)  # Delete existing directory
				print(f"Output directory '{self.output_dir}' exists. Overwriting...")
				os.makedirs(self.output_dir)
			else:
				print(f"Output directory '{self.output_dir}' already exists. Stopping run.")
				sys.exit()
		else:
			os.makedirs(self.output_dir, exist_ok = True)
			print(f"Created output directory: {self.output_dir}")

	def check_inputs_exist(
		self,
		sample_sheet_path,
		feature_table_path,
		genome_path,
		# feature_path,
		# gene_names_path
		):
		for i in [sample_sheet_path, genome_path, feature_table_path]:
			if i !=  None:
				if not os.path.exists(i):
					raise ValueError(f'{i} does not exist')
