import pandas as pd
import os

class FeatureManager:
	def __init__(self, feature_table: pd.DataFrame):
		self.registry = {}
		self._validate_and_build(feature_table)

	def _validate_and_build(self, df: pd.DataFrame):
		seen = set()
		for _, row in df.iterrows():
			feature, ftype, path = row['feature'], row['type'], row['file_path']

			if feature in seen:
				raise ValueError(f"Duplicate feature: {feature}")
			seen.add(feature)

			if ftype not in {"annotation", "presence_absence"}:
				raise ValueError(f"Invalid type: {ftype}")

			if not os.path.exists(path):
				raise FileNotFoundError(f"{path} does not exist")

			sheet = pd.read_csv(path, nrows = 5)

			for col in ['chromosome', 'start', 'end']:
				if col not in sheet.columns:
					raise ValueError(f"Missing required column '{col}' in {feature}")

			if ftype == 'annotation' and 'annotation' not in sheet.columns:
				raise ValueError(f"Missing 'annotation' column in annotation-type feature: {feature}")

			self.registry[feature] = {
				'type': ftype,
				'file_path': path
			}





# def assign_features(self, cut_sites):
#         self.cut_sites = cut_sites
#         all_standard_cuts = [{
#             'Chromosome': cut.chromosome,
#             'Start': cut.global_position['cut'],
#             'End': cut.global_position['cut'],
#             'cut_cluster': cut.cut_cluster
#         } for cut in cut_sites]

#         cuts_gr = convert_df_to_granges(pd.DataFrame(all_standard_cuts))

#         for feature_name, info in self.registry.items():
#             feature_type = info['type']
#             feature_df = info['data'].copy()
#             feature_df['Start'] = feature_df['position']
#             feature_df['End'] = feature_df['position']

#             feature_gr = convert_df_to_granges(feature_df)

#             # Overlap with cut sites
#             overlaps = batch_overlaps(gr=feature_gr, sites_gr=cuts_gr)
#             overlaps = overlaps.drop(columns=['Start', 'End', 'name2']).rename(columns={'Start_b': 'Start', 'End_b': 'End'})
            
#             # Assign overlaps to each cut site
#             for cut in cut_sites:
#                 cut.extract_features(df=overlaps)

#             # If annotation, also calculate nearest
#             if feature_type == 'annotation':
#                 nearest = batch_nearest_feature(gr=feature_gr, sites_gr=cuts_gr)
#                 nearest = nearest.drop(columns=['Start', 'End']).rename(columns={'Start_b': 'Start', 'End_b': 'End', 'name2': 'name'})
                
#                 for cut in cut_sites:
#                     cut.extract_nearest_gene(df=nearest)









