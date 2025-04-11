import time
import gzip
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
from scipy.stats import zscore
from sklearn.preprocessing import MinMaxScaler

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


def extract_annotations(df, position):
	# return df[(df['chrom'] == chromosome) & (df['start'] <= position) & (df['end'] >= position)]
	return df[(df['start'] <= position) & (df['end'] >= position)]


def slice_annotation(df, chromosome, position, tolerance=3_000_000):
	"""
	Slice the genomic feature dataframe to a nearby window to reduce overhead.
	"""
	return df[
		(df['chrom'] == chromosome) &
		(df['end'] >= position - tolerance) &
		(df['start'] <= position + tolerance)
	].copy()

def get_closest_annotation(df, position, column_name):

	# df = df[df['chrom'] == chromosome].copy()

	df["distance"] = np.where(
		(df["start"] <= position) & (df["end"] >= position),
		0,
		np.minimum(np.abs(df["start"] - position), np.abs(df["end"] - position))
	)

	df = df[df["distance"] == df["distance"].min()]
	
	return list(df[column_name].unique())[0], df['distance'].tolist()[0]








