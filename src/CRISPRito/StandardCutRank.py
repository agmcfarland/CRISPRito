import pandas as pd

class RankOperator:

	def __init__(self, variable, variable_type, condition, outcome, weight, source):
		self.variable = variable
		self.variable_type = variable_type
		self.condition = condition
		self.outcome = outcome
		self.weight = weight
		self.source = source


class StandardCutRank:

	def __init__(self, ranking_criteria = None, rank_table_weights = None, cut_sites, cut_sites_detail):
		if ranking_criteria is None:
			ranking_criteria = []
		self.ranking_criteria = ranking_criteria

		self.rank_table_weights = rank_table_weights

	def score_presence_variables(self):
		pass


	def generate_ranking_skeleton(self, sample_sheet, feature_manager=None):
		self.load_default_rank_list()
		self.load_user_sample_rank_list(sample_sheet = sample_sheet)
		self.load_user_method_rank_list(sample_sheet = sample_sheet)
		if feature_manager != None:
			self.load_user_feature_rank_list(feature_manager = feature_manager)

		self.df_skeleton = pd.DataFrame(self.ranking_criteria)
		self.df_skeleton.columns = ['variable', 'variable_type', 'condition', 'outcome', 'weight', 'source']


	def load_default_rank_list(self):
		self.ranking_criteria =	[
				['overlap',	'overlap',	'>=',	50,	1,	'default'],
				['zscore_mean',	'presence',	'>=',	1.5,	2,	'default'],
				['zscore_min',	'presence',	'>=',	1.5,	0,	'default'],
				['zscore_max',	'presence',	'>=',	1.5,	0,	'default'],
				['zscore_median',	'presence',	'>=',	1.5,	0,	'default'],
				['min_max_mean',	'presence',	'>=',	1,	0,	'default'],
				['min_max_min',	'presence',	'>=',	1,	0,	'default'],
				['min_max_max',	'presence',	'>=',	1,	0,	'default'],
				['min_max_median',	'presence',	'>=',	1,	0,	'default'],
				['local_rank_mean',	'presence',	'>=',	1,	0,	'default'],
				['local_rank_min',	'presence',	'>=',	1,	0,	'default'],
				['local_rank_max',	'presence',	'>=',	1,	0,	'default'],
				['local_rank_median',	'presence',	'>=',	1,	0,	'default'],
				['lev_distance',	'distance',	'<=',	3,	5,	'default'],
				['lev_distance',	'distance',	'<=',	7,	3,	'default'],
				['alignment_length',	'distance',	'>=',	17,	2,	'default'],
				['alignment_length',	'distance',	'<=',	22,	2,	'default']
				]

	def load_user_feature_rank_list(self, feature_manager):
		"""
		"""
		for feature_, info_ in feature_manager.registry.items():
			if info_['type'] == 'presence_absence':
				self.ranking_criteria.append([f'in_{feature_}', 'presence', '==', 1, 2, 'user'])

			if info_['type'] == 'annotation':
				self.ranking_criteria.append([f'nearest_{feature_}_distance', 'distance', '==', 0, 2, 'user'])
				self.ranking_criteria.append([f'nearest_{feature_}_distance', 'distance', '<=', 1000, 1, 'user'])

	def load_user_sample_rank_list(self, sample_sheet):
		for _, row in sample_sheet.iterrows():
			self.ranking_criteria.append([row['sample'], 'presence', '==', 1, 0, 'user'])


	def load_user_method_rank_list(self, sample_sheet):
		for method_ in sample_sheet['method'].drop_duplicates().tolist():
			self.ranking_criteria.append([method_, 'overlap', '>=', 50, 4, 'user'])





