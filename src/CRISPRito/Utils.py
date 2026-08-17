import logging
import sys
import time
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from scipy.stats import zscore
from sklearn.preprocessing import MinMaxScaler
import pyranges as pr


def setup_logging(log_path):
	log_formatter = logging.Formatter('%(asctime)s — %(levelname)s — %(message)s')

	root_logger = logging.getLogger()
	root_logger.setLevel(logging.INFO)

	file_handler = logging.FileHandler(log_path)
	file_handler.setFormatter(log_formatter)
	root_logger.addHandler(file_handler)

	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setFormatter(log_formatter)
	root_logger.addHandler(console_handler)

	# Redirect stdout and stderr to logging
	class StreamToLogger:
		def __init__(self, logger, level):
			self.logger = logger
			self.level = level

		def write(self, message):
			if message.strip():  # avoid empty messages
				self.logger.log(self.level, message.strip())

		def flush(self):
			pass

	sys.stdout = StreamToLogger(root_logger, logging.INFO)
	sys.stderr = StreamToLogger(root_logger, logging.INFO)

def alignment_levenshtein(observed, reference):
    """
    Compute Levenshtein distance from aligned sequences.

    Parameters:
    - observed (str): Aligned observed sequence (may contain '-')
    - reference (str): Aligned reference sequence (may contain '-')

    Returns:
    - int: Total edit distance (insertions + deletions + substitutions)
    """
    if len(observed) != len(reference):
        raise ValueError("Aligned sequences must be the same length.")

    insertions = deletions = substitutions = 0

    for o, r in zip(observed, reference):
        if o == '-' and r != '-':
            deletions += 1
        elif o != '-' and r == '-':
            insertions += 1
        elif o != r:
            substitutions += 1

    return insertions + deletions + substitutions


def greedy_clustering_first(numbers, range_threshold):
	"""
	The current number is compared to the first item in the current cluster when 
	comparing ranges.
	"""
	numbers = sorted(numbers)
	clusters = {}
	cluster_id = 0
	current_cluster = []
	for number in numbers:
		if not current_cluster:
			current_cluster.append(number)
			continue
		# The current number is compared against the first number of the group
		if number - current_cluster[0] <= range_threshold:
			current_cluster.append(number)
		else:
			clusters[cluster_id] = current_cluster
			cluster_id += 1
			current_cluster = [number]
	if current_cluster:
		clusters[cluster_id] = current_cluster
	return clusters

def greedy_clustering_incremental(numbers, range_threshold):
	"""
	The current number is compared to the last item in the current cluster when 
	comparing ranges.
	"""
	numbers = sorted(numbers)
	clusters = {}
	cluster_id = 0
	current_cluster = []
	for number in numbers:
		if not current_cluster:
			current_cluster.append(number)
			continue
		# The current number is compared against the last added number to the group
		if number - current_cluster[-1] <= range_threshold:
			current_cluster.append(number)
		else:
			clusters[cluster_id] = current_cluster
			cluster_id += 1
			current_cluster = [number]
	if current_cluster:
		clusters[cluster_id] = current_cluster
	return clusters


def genome_to_dict_memoryview(genome_path):
	with gzip.open(genome_path, 'rt') as infile:
		record_dict = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))
		for k, v in record_dict.items():
			record_dict[k] = memoryview(str(v.seq).encode())
		return record_dict


def sequence_slice_locations(pos, flank_size):
	start, end = max(0, pos - flank_size), pos + flank_size
	return start, end

def retrieve_genome_slices_memoryview(sequence, positions, flank_size):
	slices = {}
	for pos in positions:
		start, end = sequence_slice_locations(pos, flank_size)
		slices[pos] = sequence[start:end].tobytes().decode()
	return slices


def parse_global_alignment(aln):
	"""
	aln is an is the raw output of global_pairwise_align_nucleotide(seq1, seq2)
	"""
	sequence_alignment = aln[0].conservation()

	alignment = np.where(~np.isnan(sequence_alignment))[0]

	# this allows for counting even if gaps exist
	record = [e for e, i in enumerate(sequence_alignment) if not np.isnan(i)]

	if len(record) != 0:
		return min(record), max(record)

	return -1, -1


def numpy_list_converter(func):
	"""
	Decorator to ensure list is converted to numpy array and returned as a flattened list
	"""
	def wrapper(score, *args, **kwargs):
		if isinstance(score, list):
			score = np.array(score)
		result = func(score, *args, **kwargs)
		return result.flatten().tolist()
	return wrapper

@numpy_list_converter
def scale_min_max(score:list):

	scaler = MinMaxScaler(feature_range=(0, 1))

	scaled_score = scaler.fit_transform(score.reshape(-1,1))

	return scaled_score

@numpy_list_converter
def scale_zscore(score:list, degrees_of_freedom = 0):

	scaled_score = zscore(score, ddof = degrees_of_freedom)

	return scaled_score

def batch_overlaps(gr, sites_gr):
	return sites_gr.join(gr).df

def batch_nearest_feature(gr, sites_gr):
	return sites_gr.nearest(gr, overlap = True).df


def central_tendency(measurement: pd.core.series.Series):
	return measurement.describe()[['mean', 'min', 'max', '50%']].to_dict()


def report_time(func):
	def wrapper(*args, **kwargs):
		time_start = time.time()
		print(f'Starting: {func}')
		result = func(*args, **kwargs)
		end_time = time.time()-time_start
		print(f'Finished: {func} in {round(end_time/60, 2)} minutes')
		return result
	return wrapper

def convert_df_to_granges(df):
	"""
	Pyranges requires these columns so renaming columns to the correct format falls under converting the df to granges falls int
	"""
	df = df.rename(columns = {'chrom': 'Chromosome', 'start': 'Start', 'end': 'End', 'chromosome': 'Chromosome'})
	return pr.PyRanges(df)


def find_revised_pam(sequence, protospacer_start, search_direction = 'forward', max_pam_gaps_allowed = 2, pam_type = 'NGG'):
	if search_direction not in ['forward', 'backward']:
		raise ValueError('search_direction must be "forward" or "backward"')

	pam_length = len(pam_type)

	starting_pam = sequence[protospacer_start + 1: protospacer_start + 1 + pam_length]

	revised_pam = starting_pam

	pam_n_gap = 0 

	direction_modifier = 1

	if search_direction == 'backward':
		direction_modifier = -1

	max_pam_gaps_allowed = max_pam_gaps_allowed * direction_modifier

	while revised_pam[1:] != pam_type[1:] and pam_n_gap != max_pam_gaps_allowed:

		revised_pam = sequence[protospacer_start + 1 + pam_n_gap + direction_modifier : protospacer_start + 1 + pam_length + pam_n_gap + direction_modifier]

		pam_n_gap += direction_modifier

	return {
	'starting_pam': starting_pam,
	'revised_pam': revised_pam,
	'pam_n_gap': pam_n_gap,
	'search_direction': search_direction,
	'revised_protospacer_start': protospacer_start + pam_n_gap,
	'pam_found': revised_pam[1:] == pam_type[1:]
	}

def df_long_to_wide(df, to_rows, to_columns, column_prefix = None):
	"""
	"""
	df_wide = (
		df
		.groupby([to_rows, to_columns])
		.size()  
		.unstack(fill_value=0)  
		.reset_index()
		)
	df_wide.index.name = 'index'

	if column_prefix:

		df_wide = df_wide.rename(
			columns={col: f"{column_prefix}_{col}" for e, col in enumerate(df_wide.columns) if e != 0})

	return df_wide

def sliding_windows(seq_length, window_size, step_size):
	"""
	Generate (start, end) index pairs for subsetting a sequence.
	
	Parameters:
	- seq_length (int): Length of the sequence.
	- window_size (int): Length of each window.
	- step_size (int): Number of bases to move the window each step.

	Returns:
	- List of (start, end) tuples.
	"""
	windows = []
	for start in range(0, seq_length - window_size + 1, step_size):
		end = start + window_size
		windows.append((start, end))
	return windows
