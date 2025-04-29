from CRISPRito.SampleManager import SampleManager
from CRISPRito.RunParameters import RunParameters

def setup_run(
	input_file,
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

	run_parameters.manage_output_dir(overwrite = overwrite_output_dir)

	run_parameters.check_inputs_exist()

	manager = SampleManager(
		input_file = input_file,
		output_dir = output_dir
		)

	manager.load_samplesheet(
		
		)

	print(manager.table)

	# manager = manager.setup()


if __name__ == '__main__':

	setup_run()