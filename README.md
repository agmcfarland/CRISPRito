<img src="docs/logo.png" alt="Logo" width=1000>


`CRISPRito` is a pipeline that takes CRISPR cut coordinates from many methods and/or replicates and does two things:

1. Clusters all cut coordinates to calculate overlapping observations and scoring statistics

2. Produces a ranked table of cut sites, roughly ordered by number of overlapping observations and aggregated strength of scoring

# Installation

Download the latest distribution from anaconda.

```sh
conda create -n crisprito_env -c agmcfarland -c conda-forge -c bioconda crisprito
```

# Verify installation

Run all three commands. Please ensure `/path/to/testdir` does not exist prior to running.

* All input files will be in `/path/to/testdir/`

All CRISPRito outputs will be in `/path/to/testdir/crisprito_output`

```sh
conda activate crisprito_env

test_workdir=/path/to/testdir

crisprito-test --test_dir $test_workdir --test_workflow pull_data

crisprito-test --test_dir $test_workdir --test_workflow cluster_cuts

crisprito-test --test_dir $test_workdir --test_workflow auto_rank_sites
```


# Usage

## Cluster cuts

```
usage: crisprito-cluster-cuts [-h]

Run CRISPRito worfklows

options:
  -h, --help            show this help message and exit
  --sample_sheet_path SAMPLE_SHEET_PATH
  --output_dir OUTPUT_DIR
  --genome_path GENOME_PATH
  --feature_table_path FEATURE_TABLE_PATH
  --overwrite_output_dir
  --flank_size FLANK_SIZE
  --sgRNA SGRNA
  --PAM_alignment PAM_ALIGNMENT
  --range_threshold RANGE_THRESHOLD
                        Distance between clusters
  --workflow WORKFLOW   One of [cluster_cuts, detect_cut]
```

## Rank cuts

```
usage: crisprito-rank-sites [-h]

Rank sites generated from a CRISPRito output

options:
  -h, --help            show this help message and exit
  --workflow {standard,auto}
                        Which ranking workflow to run: 'standard' (StandardCutRank, user-weighted rule-based scoring) or 'auto' (AutoRankCuts, RRA-based consensus ranking with optional magnitude/spatial enrichment). [Default auto].
  --cluster_group CLUSTER_GROUP
                        cluster_group identifier from step 1 of the pipeline. Used with --data_dir to auto-resolve cut_profiles, method_counts, cut_id_detail, and samplesheet paths.
  --data_dir DATA_DIR   Directory containing step 1 pipeline outputs for this cluster_group (e.g. '{cluster_group}_group_cut_profiles.csv', etc.).
  --output_dir OUTPUT_DIR
                        Directory to save output CSV.
  --output_name OUTPUT_NAME
                        Name of ranked sites table.
  --rank_table_weights_path RANK_TABLE_WEIGHTS_PATH
                        [standard] Path to the weights table.
  --feature_driver_path FEATURE_DRIVER_PATH
                        ['auto] Path to the feature annotation CSV.
  --magnitude_transform {zscore,minmax,percentile,raw}
                        [auto] Per-method score normalization used in calculate_magnitude_score.
  --magnitude_aggregation {mean,max,sum}
                        [auto] How per-method magnitude scores are aggregated per cut_cluster.
  --tau TAU             [auto] Power-law decay midpoint distance (bp).
  --gamma GAMMA         [auto] Power-law decay tail-shape parameter.
  --max_distance_bp MAX_DISTANCE_BP
                        [auto] Distance cap (bp) applied before decay.
  --feature_weights [FEATURE_WEIGHTS ...]
                        [auto] Weights for combining decay columns, in the same order as --distance_cols. Defaults to equal weighting if omitted.
  --spatial_aggregation {harmonic,weighted_mean,max}
                        [auto] How per-feature decay weights are combined into a single composite_decay_weight.
  --score_col {log_p_rra,rra_weighted_magnitude_score}
                        [auto] Base score multiplied by composite_decay_weight to form the final rank. 'log_p_rra' gives a fully parameter-free ranking using RRA consensus alone.
```