import pandas as pd
import operator

class RankOperator:


	def __init__(self, variable, variable_type, condition, outcome, weight, source, type, rank_criteria_id):
		self.variable = variable
		self.variable_type = variable_type
		self.condition = condition
		self.outcome = outcome
		self.weight = weight
		self.source = source
		self.type = type
		self.rank_criteria_id = rank_criteria_id
		self.condition_map = {
			'==': operator.eq,
			'!=': operator.ne,
			'>': operator.gt,
			'>=': operator.ge,
			'<': operator.lt,
			'<=': operator.le
		}
		self.score = pd.DataFrame()
		

	@classmethod
	def load_from_rank_table_row(cls, row):
		return cls(
			variable=row['variable'],
			variable_type=row['variable_type'],
			condition=row['condition'],  # or row['modifier'] if using older naming
			outcome=row['outcome'],
			weight=row['weight'],
			source=row['source'],
			type=row['type'],
			rank_criteria_id = row['rank_criteria_id']
		)

	def score_criteria(self, datasets):

		if self.type == 'feature':
			
			if self.variable_type == 'presence' or self.variable_type == 'distance':

				df = datasets['cut_profiles'].copy()

				df[self.rank_criteria_id] = df[self.variable].apply(lambda x: self.weight if self.condition_map[self.condition](x, self.outcome) else 0)

				self.score = df[['cut_cluster', self.rank_criteria_id, self.variable]]

		# this is a special case for feature and overlap
		if self.variable_type == 'overlap':

			if self.variable == 'overlap':

				if self.source == 'default':

					max_samples = datasets['samplesheet'].shape[0]

					df = datasets['cut_profiles']

					df[self.variable] = 100 * (df[self.variable]/max_samples)

					df[self.rank_criteria_id] = df[self.variable].apply(lambda x: self.weight if self.condition_map[self.condition](x, self.outcome) else 0)

					self.score = df[['cut_cluster', self.rank_criteria_id, self.variable]]

		if self.type == 'sample':

			if self.variable_type == 'presence':

				df = datasets['id_counts'].copy()

				df = df[df[self.type] == self.variable]

				df[self.variable] = df['detected']

				df[self.rank_criteria_id] = df[self.variable].apply(lambda x: self.weight if self.condition_map[self.condition](x, self.outcome) else 0)

				self.score = df[['cut_cluster', self.rank_criteria_id, self.variable]]

		if self.type == 'method':

			if self.variable_type == 'overlap':

				max_samples = datasets['samplesheet'].shape[0]

				df = datasets['method_counts'].copy()

				df = df[df[self.type] == self.variable]

				df[self.variable] = 100 * (df['detected']/max_samples)

				df[self.rank_criteria_id] = df[self.variable].apply(lambda x: self.weight if self.condition_map[self.condition](x, self.outcome) else 0)

				self.score = df[['cut_cluster', self.rank_criteria_id, self.variable]]


class StandardCutRank:

	def __init__(self, rank_table_weights = None, cut_profiles=None, id_counts=None, method_counts=None, samplesheet=None):
		self.ranking_criteria = []
		self.rank_table_weights = rank_table_weights
		self.cut_profiles = cut_profiles
		self.id_counts = id_counts
		self.method_counts = method_counts
		self.samplesheet = samplesheet

		if self.samplesheet is not None:
			if self.id_counts is not None:
				self.id_counts = self.id_counts.merge(self.samplesheet[['id', 'sample']], on = 'id')
				self.id_counts = self.id_counts.drop(columns = 'id')
				self.id_counts = pd.melt(self.id_counts, id_vars = ['sample'], var_name = 'cut_cluster', value_name = 'detected')
				self.id_counts['cut_cluster'] = self.id_counts['cut_cluster'].apply(lambda x: int(x.replace('cut_', '')))

		if self.method_counts is not None:
			self.method_counts = pd.melt(self.method_counts, id_vars = ['method'], var_name = 'cut_cluster', value_name = 'detected')
			self.method_counts['cut_cluster'] = self.method_counts['cut_cluster'].apply(lambda x: int(x.replace('cut_', '')))

		if self.rank_table_weights is not None:
			self.rank_table_weights['rank_criteria_id'] = [f'rank_criteria_{i}' for i in range(0, len(self.rank_table_weights))]

		self.datasets = {
			'rank_table_weights' : self.rank_table_weights,
			'cut_profiles' : self.cut_profiles,
			'id_counts' : self.id_counts,
			'method_counts' : self.method_counts,
			'samplesheet' : self.samplesheet
			}

	def get_score_from_rank_table(self):
		
		self.scored_criteria = []

		for _, row in self.rank_table_weights.iterrows():

			rank_criteria = RankOperator.load_from_rank_table_row(row)

			rank_criteria.score_criteria(datasets = self.datasets)

			self.scored_criteria.append(rank_criteria)


	def tally_cut_cluster_scores(self):

		self.cut_cluster_scores = pd.DataFrame({'cut_cluster' : self.datasets['cut_profiles']['cut_cluster'].tolist()})

		for rank_criteria in self.scored_criteria:

			self.cut_cluster_scores = self.cut_cluster_scores.merge(rank_criteria.score[['cut_cluster', rank_criteria.rank_criteria_id]], on = 'cut_cluster')
			
		self.cut_cluster_scores['total_score'] = self.cut_cluster_scores.drop(columns='cut_cluster').sum(axis=1)
		

	def generate_ranking_skeleton(self, sample_sheet, feature_manager=None):
		self.load_default_rank_list()
		self.load_user_sample_rank_list(sample_sheet = sample_sheet)
		self.load_user_method_rank_list(sample_sheet = sample_sheet)
		if feature_manager != None:
			self.load_user_feature_rank_list(feature_manager = feature_manager)

		self.df_skeleton = pd.DataFrame(self.ranking_criteria)
		self.df_skeleton.columns = ['variable', 'variable_type', 'condition', 'outcome', 'weight', 'source', 'type']
		self.df_skeleton = self.df_skeleton[['variable', 'variable_type', 'source', 'type', 'condition', 'outcome', 'weight']]

	def load_default_rank_list(self):
		self.ranking_criteria =	[
				['overlap',	'overlap',	'>=',	50,	1,	'default', 'feature'],
				['zscore_mean',	'presence',	'>=',	1.5,	2,	'default', 'feature'],
				['zscore_min',	'presence',	'>=',	1.5,	0,	'default', 'feature'],
				['zscore_max',	'presence',	'>=',	1.5,	0,	'default', 'feature'],
				['zscore_median',	'presence',	'>=',	1.5,	0,	'default', 'feature'],
				['min_max_mean',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['min_max_min',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['min_max_max',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['min_max_median',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['local_rank_mean',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['local_rank_min',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['local_rank_max',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['local_rank_median',	'presence',	'>=',	1,	0,	'default', 'feature'],
				['lev_distance',	'distance',	'<=',	3,	5,	'default', 'feature'],
				['lev_distance',	'distance',	'<=',	7,	3,	'default', 'feature'],
				['alignment_length',	'distance',	'>=',	17,	2,	'default', 'feature'],
				['alignment_length',	'distance',	'<=',	22,	2,	'default', 'feature']
				]

	def load_user_feature_rank_list(self, feature_manager):
		"""
		"""
		for feature_, info_ in feature_manager.registry.items():
			if info_['type'] == 'presence_absence':
				self.ranking_criteria.append([f'in_{feature_}', 'presence', '==', 1, 2, 'user', 'feature'])

			if info_['type'] == 'annotation':
				self.ranking_criteria.append([f'nearest_{feature_}_distance', 'distance', '==', 0, 2, 'user', 'feature'])
				self.ranking_criteria.append([f'nearest_{feature_}_distance', 'distance', '<=', 1000, 1, 'user', 'feature'])

	def load_user_sample_rank_list(self, sample_sheet):
		for _, row in sample_sheet.iterrows():
			self.ranking_criteria.append([row['sample'], 'presence', '==', 1, 0, 'user', 'sample'])


	def load_user_method_rank_list(self, sample_sheet):
		for method_ in sample_sheet['method'].drop_duplicates().tolist():
			self.ranking_criteria.append([method_, 'overlap', '>=', 50, 4, 'user', 'method'])





