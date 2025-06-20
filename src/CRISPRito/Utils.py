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
		print('Start:', func)
		result = func(*args, **kwargs)
		print('End:', func, time.time()-time_start)
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

def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

# # def extract_annotations(df, position):
# # 	return df[(df['start'] <= position) & (df['end'] >= position)]


# def slice_annotation(gr, chromosome, position, tolerance=3_000_000):
# 	"""
# 	Slice a PyRanges object to a nearby window.
# 	"""
# 	gr_chr = gr[gr.Chromosome == chromosome] 
# 	gr_chr = gr_chr[position - tolerance: position + tolerance]
# 	return gr_chr


# def extract_annotations(gr, chromosome, position):
# 	"""Return all rows where the cut site overlaps a genomic feature."""
# 	site = pr.from_dict({
# 		"Chromosome": [chromosome],
# 		"Start": [position],
# 		"End": [position]
# 	})

# 	return gr.overlap(site).df

# def get_closest_annotation(gr, chromosome, position, column_name="name2"):
# 	"""Find nearest annotation and its distance."""
# 	site = pr.from_dict({
# 		"Chromosome": [chromosome],
# 		"Start": [position],
# 		"End": [position]
# 	})

# 	nearest = site.k_nearest(gr, k = 1)

# 	return nearest.df[column_name].values[0], abs(nearest.df["Distance"].values[0])




# def batch_overlaps(gr, sites_gr):
# 	return sites_gr.join(gr).df

# def batch_nearest_feature(gr, sites_gr):
# 	return sites_gr.nearest(gr, overlap = True).df
# 	# return gr.k_nearest(sites_gr, overlap = True, k = 1).df


# 	if len(overlaps) > 0:
# 		result = overlaps
# 	else:
# 		# Step 2: Get closest if no overlaps
# 		# result = sites_gr.k_nearest(gr, k=1)
# 		overlaps = sites_gr.nearest(gr, overlap = True).df
# 		print('this happened')

# 	print('\n')
# 	print(result)
# 	# print('\n')
# 	# print(overlaps)
# 	# Display

# def extract_annotations(df, position):
# 	start_le = df['start'].values <= position
# 	end_ge = df['end'].values >= position
# 	return df[start_le & end_ge]

# def slice_annotation(df, chromosome, position, tolerance=3_000_000):
# 	"""
# 	Slice the genomic feature dataframe to a nearby window to reduce overhead.
# 	"""
# 	return df[
# 		(df['chrom'] == chromosome) &
# 		(df['end'] >= position - tolerance) &
# 		(df['start'] <= position + tolerance)
# 	].copy()

# def get_closest_annotation(df, position, column_name):

# 	df["distance"] = np.where(
# 		(df["start"] <= position) & (df["end"] >= position),
# 		0,
# 		np.minimum(np.abs(df["start"] - position), np.abs(df["end"] - position))
# 	)

# 	df = df[df["distance"] == df["distance"].min()]
	
# 	return list(df[column_name].unique())[0], df['distance'].tolist()[0]







