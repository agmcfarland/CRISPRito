import pandas as pd
import uuid
import os
from os.path import join as pjoin

class SampleManager:
	allowed_measurements = [
		'one_scaled',
		'abundance'
	]
	def __init__(self, input_file, output_dir="CRISPRito_output"):
		"""
		Initializes the SampleManager with an input file path.
		
		:param input_file: Path to the sample sheet containing sample data.
		"""
		self.input_file = input_file
		self.output_dir = output_dir
		self.table = None

	def load_samplesheet(self):
		self.table = pd.read_csv(self.input_file)

	def assign_unique_id(self):
		"""Reads the input file into a pandas DataFrame and assigns unique IDs."""
		self.table["id"] = [str(uuid.uuid4()) for _ in range(len(self.table))]

	def enforce_measurement_type(self):
		for i in self.table['measurement_type'].unique():
			if i not in self.allowed_measurements:
				raise ValueError(f'{i} not a valid measurment type.')

	def write_samples_by_group(self):
		unique_groups = list(self.table['cluster_group'].unique())
		for ug in unique_groups:
			self.table[self.table['cluster_group'] == ug].to_csv(pjoin(self.output_dir, f'{ug}_group_samplesheet.csv'), index = None)

	def set_up(self):
		self.load_samplesheet()
		self.assign_unique_id()
		self.enforce_measurement_type()
		self.test_write_samples_by_group()


	# def write_samples(self):
	#     """Writes each row to a separate text file using its unique ID as the filename."""
	#     if self.samples is None:
	#         raise ValueError("Samples must be loaded first using load_samples().")

	#     for _, row in self.samples.iterrows():
	#         file_path = os.path.join(self.output_dir, f"{row['id']}.txt")
			
	#         with open(file_path, "w") as f:
	#             f.write(row.to_json(indent=4))  # Write row data as JSON for readability

	# def process(self):
	# 	"""Runs the full pipeline: loading samples and writing output files."""
	# 	self.load_samples()
	# 	self.write_samples()

