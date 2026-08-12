import numpy as np
import pandas as pd
from scipy.stats import beta, norm

class AutoRankCuts:

	def __init__(self, cut_profile=None, cut_id_detail=None):
		self.cut_profile = cut_profile
		self.cut_id_detail = cut_id_detail
		self.rra_scores = None
		self.magnitude_score = None
		self.feature_weight = None
		self.method_presence = None

	def calculate_rra_scores(self):
		"""
		Performs Standard Robust Rank Aggregation (RRA) on cut data.
		Collects the max score per method per cut site
		Ranks each cut site per method
		Calculates p-value under beta CDF distribution
		Returns cut_cluster, p_rra, and log_p_rra 

		Notes: 
			Cuts not detected by a method are assigned 1 (worst score)
			Uses all cut sites detected over all methods as denominator. 
		"""
		df = self.cut_id_detail.copy()

		df['score'] = pd.to_numeric(df['score'], errors='coerce')

		# Max score per cut_cluster-method pair
		method_collapsed = (
			df.groupby(['cut_cluster', 'method'], as_index=False)
			.agg(best_score=('score', 'max'))
		)

		# Per-method ranking of cut_clusters
		method_collapsed['method_rank'] = (
			method_collapsed.groupby('method')['best_score']
			.rank(ascending=False, method='min')
		)

		all_clusters = np.sort(df['cut_cluster'].unique())
		all_methods = np.sort(df['method'].unique())
		N_total = float(len(all_clusters))

		# Normalized rank positions u_im in (0, 1]
		method_collapsed['u_im'] = method_collapsed['method_rank'] / N_total

		# Create matrix of all cut sites and methods, fill NA (missed cut clusters by methods)
		# with 1 (max proportion possible)
		full_grid = pd.MultiIndex.from_product(
			[all_clusters, all_methods], 
			names=['cut_cluster', 'method']
		).to_frame().reset_index(drop=True)

		eval_matrix = full_grid.merge(
			method_collapsed[['cut_cluster', 'method', 'u_im']], 
			on=['cut_cluster', 'method'], 
			how='left'
		)

		# Assign non-detected entries upper-bound quantile 1.0
		eval_matrix['u_im'] = eval_matrix['u_im'].fillna(1.0)

		def _calc_rra_pvalue(u_series):
			"""
			Takes as input u_series, which is the normalized rank position from each method (m) for a cut site.
			Notes:
				u_vals: normalized rank quantiles (0-1)
				m: total number of methods evaluated
				k_vec: ordered 1-D array of methods of length 1 to m
				a: count of successs/top ranks up to position K
				b: count of remaiing lower-ranked methods
				p_beta: order statistic p value. Probability of observing a quantile as extreme as u(k) at rank position k by chance.
				rho: minimal probility value in p_beta
			"""
			u_vals = np.sort(u_series.values)
			m = len(u_vals)
			k_vec = np.arange(1, m + 1)
			
			p_beta = beta.cdf(u_vals, a=k_vec, b=m - k_vec + 1)
			rho = np.min(p_beta)
			rra_p = min(rho * m, 1.0)
			return rra_p

		# calculate per cut_cluster significance under beta cdf
		rra_summary = (
			eval_matrix.groupby('cut_cluster')['u_im']
			.apply(_calc_rra_pvalue)
			.reset_index(name='p_rra')
		)

		# Transform to log10 p value for interpretability and sort on this column
		rra_summary['log_p_rra'] = -np.log10(np.maximum(rra_summary['p_rra'], 1e-300))
		rra_summary = rra_summary.sort_values('p_rra', ascending=True).reset_index(drop=True)

		self.rra_scores = rra_summary
		return self.rra_scores

	def calculate_magnitude_score(self, transform='percentile', aggregation='max'):
		"""
		Computes a method-level magnitude score per cut_cluster and scales standard 
		Returns cut_cluster and magnitude_score
		
		Notes: 
			The highest raw score from each cut_cluster-method combo is kept to represent the method
			Normalization happens on a per method scale over all cut_clusters
			An aggregation statistic is computed per cut_cluster accross methods


		RRA scores by multiplying log_p_rra by the aggregated magnitude weight.

		"""
		df = self.cut_id_detail.copy()
		# df['score'] = pd.to_numeric(df['score'], errors='coerce').fillna(0.0)

		# Highest score per cut_cluster-method combination
		method_collapsed = (
			df.groupby(['cut_cluster', 'method'], as_index=False)
			.agg(best_score=('score', 'max'))
		)

		# Generate norm_score, which can be calculated from a variety of transform calls
		def _scale_series(series):
			scores = series.values
			if transform == 'zscore':
				std_val = np.std(scores) if np.std(scores) > 0 else 1.0
				z = (scores - np.mean(scores)) / std_val
				return np.maximum(0.1, 1.0 + z)
			elif transform == 'minmax':
				rng = np.ptp(scores) if np.ptp(scores) > 0 else 1.0
				return (scores - np.min(scores)) / rng
			elif transform == 'percentile':
				return series.rank(pct=True, method='average')
			elif transform == 'raw':
				return scores
			else:
				raise ValueError("transform must be 'zscore', 'minmax', 'percentile', or 'raw'")

		method_collapsed['norm_score'] = (
			method_collapsed.groupby('method')['best_score']
			.transform(_scale_series)
		)

		# Per cut_cluster score summarization
		if aggregation == 'mean':
			mag_summary = method_collapsed.groupby('cut_cluster')['norm_score'].mean()
		elif aggregation == 'max':
			mag_summary = method_collapsed.groupby('cut_cluster')['norm_score'].max()
		elif aggregation == 'sum':
			mag_summary = method_collapsed.groupby('cut_cluster')['norm_score'].sum()
		else:
			raise ValueError("aggregation must be 'mean', 'max', or 'sum'")

		mag_summary = mag_summary.reset_index(name='magnitude_score')

		# self.rra_scores table already contains a full list of all cut_clusters detected
		rra_weighted = self.rra_scores.merge(mag_summary, on='cut_cluster', how='left')
		
		# Replace NA with 0.1 as a neutral missing value
		rra_weighted['magnitude_score'] = rra_weighted['magnitude_score'].fillna(0.1)
		
		rra_weighted = rra_weighted.reset_index(drop=True)

		rra_weighted = rra_weighted[['cut_cluster', 'magnitude_score']]

		self.magnitude_score = rra_weighted
		return self.magnitude_score


	def compute_method_support_summary(self):
		"""
		Computes a summary DataFrame per cut_cluster detailing the total number of 
		detecting methods and the replicate count for each method.

		Notes:
			Get replicate counts per method per cut_cluster
			Get total unique method observations per cut_cluster
		"""

		df = self.cut_id_detail.copy()

		# Count total replicates per cut_cluster and method group
		rep_counts = (
			df.groupby(['cut_cluster', 'method'])['id']
			.nunique()
			.reset_index(name='replicate_count')
			)

		# Each method is turned into a column, with its name, with replicate counts
		df_method = rep_counts.pivot(
			index= 'cut_cluster', 
			columns= 'method', 
			values= 'replicate_count'
		).fillna(0).astype(int)

		# Calculate total unique methods detecting each cut_cluster
		df_method.insert(0, 'total_methods', (df_method > 0).sum(axis=1))

		self.method_presence = df_method.reset_index().sort_values('total_methods', ascending=False).reset_index(drop=True)

		return self.method_presence

	def calculate_power_law_distance_decay(
		self, 
		distance_col, 
		tau=5000.0, 
		gamma=0.5, 
		max_distance_bp=10000.0
	):
		"""
		Scales distance metrics into a [0, 1] proximity weight using a power-law 
		decay function. Distance 0 bp yields weight 1.0; large distances decay toward 0.0.

		Power law is used as opposed to another function due to the control over decay. Don't want 
		values to go to zero quickly.

		returns cut_cluster, distance col, weighted distance col. 

		Merges onto an exisiting self.feature_weight col
		 
		Notes:
			Tau: distance where weight is 0.5
				Tau at 5,000 bp
			Gamma: tail shape and decay 
				Gamma at 0.5 is a faster initial drop and a slower decay
			max_distance_bp: Truncates scores greater than that this distance.
		"""
		df = self.cut_profile[['cut_cluster', distance_col]].copy()

		df = df.drop_duplicates()

		# Parse numeric series and replace unannotated -1 with max_distance_bp
		distances = pd.to_numeric(df[distance_col].replace(-1, max_distance_bp), errors='coerce').fillna(max_distance_bp)

		# Enforce physical distance cap
		if max_distance_bp is not None:
			distances = np.clip(distances, a_min=0.0, a_max=float(max_distance_bp))

		# Power-Law Proximity Weight: 1.0 at d=0, dropping toward 0.0 as distance increases
		normalized_dist = distances / float(tau)
		decay_weights = 1.0 / (1.0 + np.power(normalized_dist, float(gamma)))

		# Assign weight column
		output_col = f'{distance_col}_decay_weight'
		df[output_col] = decay_weights

		df = df[['cut_cluster', distance_col, output_col]]

		if self.feature_weight is None:
			self.feature_weight = df
		else:
			self.feature_weight = self.feature_weight.merge(
				df,
				how = 'left',
				on = 'cut_cluster'
				)

		return df

	def construct_crisprito_rank(
		self, 
		feature_weights=None, 
		aggregation='weighted_mean', 
		score_col='rra_weighted_magnitude_score'
	):
		"""
		Computes the CRISPRito Priority Index (CPI). 
		Expects decay_cols where 1.0 = highly proximate/high priority and 0.0 = distant.

		score_col can be 'log_p_rra' or 'rra_weighted_magnitude_score'
		"""

		if score_col not in ['log_p_rra', 'rra_weighted_magnitude_score']:
			raise ValueError('Inputted score col not accepted.')
		
		# Combine all ranking data in order of generation
		self.priority_index = self.rra_scores.copy()

		self.priority_index = self.priority_index.merge(
			self.magnitude_score,
			how = 'left',
			on = 'cut_cluster'
			)

		self.priority_index = self.priority_index.merge(
			self.feature_weight,
			how = 'left',
			on = 'cut_cluster'
			)

		self.priority_index['rra_weighted_magnitude_score'] = self.priority_index['log_p_rra'] * self.priority_index['magnitude_score'] 

		decay_cols = [c for c in list(self.feature_weight.columns) if c.find('decay_weight')>-1]

		if feature_weights is None:
			feature_weights = [1.0] * len(decay_cols)

		if len(feature_weights) != len(decay_cols):
			raise ValueError('feature_weights input length mismatch with decay_cols length')
		
		weights_arr = np.array(feature_weights, dtype=float)
		weights_arr /= weights_arr.sum()

		# Bound decay weights between 0.0001 and 1.0
		decay_matrix = self.priority_index[decay_cols].to_numpy(dtype=float)
		# decay_matrix = np.clip(self.priority_index[decay_cols].to_numpy(dtype=float), a_min=1e-4, a_max=1.0)

		# Aggregation across spatial features
		if aggregation == 'harmonic':
			weighted_reciprocals = weights_arr[None, :] / decay_matrix
			denom = np.sum(weighted_reciprocals, axis=1)
			composite_decay = 1.0 / denom
		elif aggregation == 'weighted_mean':
			composite_decay = np.sum(decay_matrix * weights_arr[None, :], axis=1)
		elif aggregation == 'max':
			# Takes the strongest proximity signal (e.g., if close to gene OR oncogene)
			composite_decay = np.max(decay_matrix, axis=1)
		else:
			raise ValueError("aggregation must be 'harmonic', 'weighted_mean', or 'max'")

		# Create column with the composite decay weight generated above
		composite_decay = np.clip(composite_decay, a_min=0.0, a_max=1.0)
		self.priority_index['composite_decay_weight'] = composite_decay

		# Priority index and rank
		self.priority_index['combined_priority_index'] = self.priority_index[score_col] * self.priority_index['composite_decay_weight']
		self.priority_index['crisprito_rank'] = (
			self.priority_index['combined_priority_index']
			.rank(ascending=False, method='min')
			.astype(int)
		)

		self.priority_index = self.priority_index.sort_values('crisprito_rank').reset_index(drop=True)

		self.priority_index = self.priority_index[['cut_cluster', 'crisprito_rank', 'combined_priority_index', 'composite_decay_weight']]

		return self.priority_index

	def prepare_combined_rank_output(self):
		"""
		"""
		# Prepare output
		df = self.cut_profile.copy()

		df = df[['cut_cluster', 'chromosome', 'strand', 'cut', 'nearest_gene']]

		df = df.merge(
			self.priority_index,
			how = 'left',
			on = 'cut_cluster'
			)

		df = df.merge(
			self.rra_scores,
			how = 'left',
			on = 'cut_cluster'
			)

		df = df.merge(
			self.magnitude_score,
			how = 'left',
			on = 'cut_cluster'
			)

		df = df.merge(
			self.feature_weight,
			how = 'left',
			on = 'cut_cluster'
			)

		df = df.merge(
			self.method_presence,
			how = 'left',
			on = 'cut_cluster'
			)

		df = df.sort_values('crisprito_rank', ascending = True)

		return df
