# import timeit

# Function 1: Dictionary-based Clustering
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
		# The current number is compared against the last added number to the group
		# if number - current_cluster[-1] <= range_threshold:
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
		# if number - current_cluster[0] <= range_threshold:
			current_cluster.append(number)
		else:
			clusters[cluster_id] = current_cluster
			cluster_id += 1
			current_cluster = [number]
	if current_cluster:
		clusters[cluster_id] = current_cluster
	return clusters


def convert_sequence_to_memory_view():pass


def retrieve_genome_slices(sequence, positions, flank_size, sequence_memoryview):
	sequence_memoryview #seq_view = memoryview(sequence.encode())  # Convert string to memoryview
	slices = {}
	for pos in positions:
		start, end = max(0, pos - flank_size), pos + flank_size
		slices[pos] = sequence_memoryview[start:end].tobytes().decode()
	return slices






def fast_fasta_dict(genome_path):
	"""Uses zcat to decompress and stream into SeqIO.parse()"""


for 

import subprocess
genome_path = '/data/GenomicTrackRepository/data/processed/hg38/hg38.fasta.gz'
with subprocess.Popen(["zcat", genome_path], stdout=subprocess.PIPE, text = True) as proc:
	print(proc.stdout)
	
	break
	records = SeqIO.to_dict(SeqIO.parse(proc.stdout, "fasta"))


from Bio import SeqIO
import gzip
# records = list(SeqIO.parse("example.fasta", "fasta"))

genome_path = '/data/GenomicTrackRepository/data/processed/hg38/hg38.fasta.gz'

with gzip.open(genome_path, 'rt') as infile:



	record_dict = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))

	record_dict_2 = 

for 

for chromosome, sequence in record_dict.items():
	print(chromosome)
	sequence = str(sequence.seq)
	break



def get_flanking_sequences(positions, reference_sequence, flank_size=30):
	"""
	Retrieves the upstream and downstream sequences (default: 30bp) for given positions.

	Parameters:
	- positions (list): List of integer positions (0-based index).
	- reference_sequence (str): The full reference sequence as a string.
	- flank_size (int): Number of bases to include upstream and downstream.

	Returns:
	- dict: A dictionary where keys are positions, and values are tuples (upstream, downstream).
	"""
	flanking_sequences = {}

	for pos in positions:
		upstream_start = max(0, pos - flank_size)  # Ensure start doesn't go negative
		downstream_end = min(len(reference_sequence), pos + flank_size)  # Ensure end doesn't exceed sequence length
		
		upstream = reference_sequence[upstream_start:pos]
		downstream = reference_sequence[pos:downstream_end]
		
		flanking_sequences[pos] = (upstream, downstream)
	
	return flanking_sequences

# Example usage:
reference_seq = "ACTG" * 100  # Example reference sequence of 400 bases
positions = [50, 100, 150, 200]  # Example positions

result = get_flanking_sequences(positions, reference_seq)
for pos, (up, down) in result.items():
	print(f"Position {pos}:\nUpstream: {up}\nDownstream: {down}\n")

import time

short_list = [10**3, 10**4, 10**5, 10**6, 10**7]
start_time = time.time()
for f in [1, 10, 100, 1000]:
	long_list = short_list*f
	long_list = sorted(long_list)
	iterator_start_time = time.time()
	for i in long_list:
		position = 10**i
		upstream_position = 30-position
		downstream_position = 30+position
		reduced_sequence = sequence[upstream_position:downstream_position]
	print(f'iterator {f}:', time.time() - iterator_start_time)
print('final:', print(time.time() - start_time))

10**7

sequence_memoryview = memoryview(str(sequence).encode())

long_list = short_list*10

def memoryview_slicing(sequence, positions, flank_size, sequence_memoryview):
	sequence_memoryview #seq_view = memoryview(sequence.encode())  # Convert string to memoryview
	slices = {}
	for pos in positions:
		start, end = max(0, pos - flank_size), pos + flank_size
		slices[pos] = sequence_memoryview[start:end].tobytes().decode()  # Extract efficiently
		print(sequence_memoryview[start:end].tobytes().decode())
	return slices


memoryview_slicing(sequence=sequence, positions=long_list, flank_size=30, sequence_memoryview=sequence_memoryview)

sequence_memoryview[0:10].tobytes().decode()

import numpy as np

sequence_np = np.array(list(sequence))  # Convert to NumPy array
positions = [10**3, 10**4, 10**5, 10**6, 7]

def numpy_slicing(sequence, positions, flank_size=30):
	return {pos: "".join(sequence[max(0, pos - flank_size): pos + flank_size]) for pos in positions}





for e, i in enumerate()

# # Test data
# numbers = [1, 1, 1, 3, 2, 8, 10, 15, 12, 18, 20, 25, 30, 35, 40, 45, 50, 1000]
# range_threshold = 20

# # Timeit comparison for 1000 executions of each function
# dict_time = timeit.timeit('greedy_clustering_dict(numbers, range_threshold)', 
# 						  globals=globals(), 
# 						  number=1000)

# list_time = timeit.timeit('greedy_clustering_dict_last(numbers, range_threshold)', 
# 						  globals=globals(), 
# 						  number=1000)

# print(f"Dictionary-based clustering time: {dict_time} seconds")
# print(f"greedy_clustering_dict_last clustering time: {list_time} seconds")

# print(greedy_clustering_dict(numbers, range_threshold))
# print(greedy_clustering_dict_last(numbers, range_threshold))