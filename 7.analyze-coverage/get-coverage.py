"""
2018 Gregory Way
7.analyze-coverage/get-coverage.py

For a given dataset and gene set compendia extract coverage metrics. The coverage
metrics include the propotion of genesets falling as the top ranked feature in
independent and ensemble models. We also extract the z score for these top features.

Usage:

    python get-coverage.py

Output:
Three dataframes storing gene set coverage and z score values for the given input
dataset and geneset. The three dataframes include:

1. Individual models
2. Ensemble models - aggregate seeds within algorithms and z dimensions
3. All model - aggregate all seeds across algorithms but within z dimension
"""

import os
import sys
import argparse
import pandas as pd

sys.path.append("../scripts")
from latent import parse_gmt

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dataset', choices=['TCGA', 'GTEX', 'TARGET'],
                    help='the dataset used to explore results')
parser.add_argument('-g', '--geneset', help='the geneset to evaluate coverage')
args = parser.parse_args()

# Load command arguments
dataset = args.dataset.lower()
geneset = args.geneset.lower()


def get_minmax_geneset(row, minmax='min'):
    """
    Determine the highest z score and corresponding geneset. Used in pandas apply().

    Arguments:
    row - a pandas dataframe row - used in apply
    minmax - a string (either "min" or "max") to find highest or lowest feature

    Output:
    A pandas series of the top z score and geneset combination
    """
    if minmax == 'min':
        idx = row.z_score.idxmin()
    else:
        idx = row.z_score.idxmax()

    z_score = row.z_score[idx]
    geneset = row.variable[idx]

    return pd.Series([z_score, geneset], index=['top_z_score', 'geneset'])


def assign_attributes(df, dataset, z_dim, geneset):
    """
    Create new columns for the given input data frame - will be used for plotting

    Arguments:
    df - a pandas dataframe storing individual, ensemble, or all results
    dataset - a string indicating the dataset of interest
    z_dim - a string of the internal dimensionality of the given model
    geneset - a string indicating the geneset of interest

    Output:
    The same dataset with additional columns to be used for plotting downstream
    """
    df_copy = df.copy()
    return df_copy.assign(dataset=dataset,
                          z_dim=z_dim,
                          geneset_name=geneset)


# Load gmt file
base_gmt_dir = os.path.join('..', '3.build-hetnets', 'data')
gmt_file_dict = {
    'gpc1': 'c1.all.v6.1.entrez.gmt',
    'gpc2cpg': 'c2.cgp.v6.1.entrez.gmt',
    'gpc2creactome': 'c2.cp.reactome.v6.1.entrez.gmt',
    'gpc3mir': 'c3.mir.v6.1.entrez.gmt',
    'gpc3tft': 'c3.tft.v6.1.entrez.gmt',
    'gpc4cgn': 'c4.cgn.v6.1.entrez.gmt',
    'gpc4cm': 'c4.cm.v6.1.entrez.gmt',
    'gpc5bp': 'c5.bp.v6.1.entrez.gmt',
    'gpc5cc': 'c5.cc.v6.1.entrez.gmt',
    'gpc5mf': 'c5.mf.v6.1.entrez.gmt',
    'gpc6': 'c6.all.v6.1.entrez.gmt',
    'gpc7': 'c7.all.v6.1.entrez.gmt',
    'gph': 'h.all.v6.1.entrez.gmt',
    'gpxcell': 'xcell_all_entrez.gmt'
}

gmt_file = os.path.join(base_gmt_dir, gmt_file_dict[geneset])
gmt_data = parse_gmt(gene_sets=[gmt_file])
n_gmt = len(gmt_data)

# Build common groupby lists
individual_features = ['algorithm', 'feature', 'seed', 'model_type']
individual_models = ['algorithm', 'seed', 'model_type']
ensemble_models = ['algorithm', 'model_type']
all_models = ['model_type']

# What is the directory to search for results in
base_dir = os.path.join('..', '6.analyze-weights', 'results', dataset, geneset, 'signal')

model_results_list = []
ensemble_results_list = []
all_results_list = []
top_list = []

for file in os.listdir(base_dir):
    # Extract z dimension
    z_dim = file.split('_')[2]

    # Build full file path
    file = os.path.join(base_dir, file)
    print("Now processing... {}".format(file))

    # Read data
    df = pd.read_table(file)

    # Select the ceiling and floor of each feature
    # within each model (model_type, algorithm, and seed)
    min_df = (
        df.groupby(individual_features)
        .apply(lambda x: get_minmax_geneset(x, minmax='min'))
        .reset_index()
    )

    max_df = (
        df.groupby(individual_features)
        .apply(lambda x: get_minmax_geneset(x, minmax='max'))
        .reset_index()
    )

    top_df = pd.concat([min_df, max_df])
    top_df = top_df.assign(abs_top_z_score=top_df.top_z_score.abs())

    # TEST 1 - Coverage
    # Get coverage per individual model (algorithm, seed, model_type)
    model_coverage_df = (
        (
            top_df.groupby(individual_models)
            .geneset
            .nunique() / n_gmt
        )
        .reset_index()
        .rename({'geneset': 'geneset_coverage'}, axis='columns')
    )

    # Also get coverage per ensemble model (algorithm, model_type)
    ensemble_coverage_df = (
        (
            top_df.groupby(ensemble_models)
            .geneset
            .nunique() / n_gmt
        )
        .reset_index()
        .rename({'geneset': 'geneset_coverage'}, axis='columns')
    )

    # And, lastly, get coverage for ALL models
    all_coverage_df = (
        (
            top_df.groupby(all_models)
            .geneset
            .nunique() / n_gmt
        )
        .reset_index()
        .rename({'geneset': 'geneset_coverage'}, axis='columns')
    )

    # TEST 2 - Strength of signal (by z-score)
    model_zscore_summary_df = (
        top_df.groupby(individual_models)
        .apply(lambda x: x.abs_top_z_score
               .sort_values(ascending=False).head(1))
    ).reset_index()

    # Also get coverage per ensemble model (algorithm, model_type)
    ensemble_zscore_summary_df = (
        top_df.groupby(ensemble_models)
        .apply(lambda x: x.abs_top_z_score
               .sort_values(ascending=False).head(1))
    ).reset_index()

    # And, lastly, get coverage for ALL models
    all_zscore_summary_df = (
        top_df.groupby(all_models)
        .apply(lambda x: x.abs_top_z_score
               .sort_values(ascending=False).head(1))
    )
    all_zscore_summary_df.columns = ['abs_top_z_score']

    # Compile all results together
    model_results_df = (
        model_coverage_df
        .merge(model_zscore_summary_df, on=individual_models)
    )

    ensemble_results_df = (
        ensemble_coverage_df
        .merge(ensemble_zscore_summary_df, on=ensemble_models)
    )

    all_results_df = (
        all_coverage_df
        .merge(all_zscore_summary_df, on=all_models)
    )

    # Track additional metrics
    model_results_df = (
        assign_attributes(df=model_results_df,
                          dataset=dataset,
                          z_dim=z_dim,
                          geneset=geneset)
    )
    ensemble_results_df = (
        assign_attributes(df=ensemble_results_df,
                          dataset=dataset,
                          z_dim=z_dim,
                          geneset=geneset)
    )
    all_results_df = (
        assign_attributes(df=all_results_df,
                          dataset=dataset,
                          z_dim=z_dim,
                          geneset=geneset)
    )

    top_df = (
        assign_attributes(df=top_df,
                          dataset=dataset,
                          z_dim=z_dim,
                          geneset=geneset)
    )

    # Save results
    model_results_list.append(model_results_df)
    ensemble_results_list.append(ensemble_results_df)
    all_results_list.append(all_results_df)
    top_list.append(top_df)

# Compile results
complete_model_results_df = pd.concat(model_results_list).reset_index(drop=True)
complete_ensemble_results_df = pd.concat(ensemble_results_list).reset_index(drop=True)
complete_all_results_df = pd.concat(all_results_list).reset_index(drop=True)
complete_top_results_df = pd.concat(top_list).reset_index(drop=True)

# Save results to file
file = 'model_results_{}_{}.tsv'.format(dataset, geneset)
file = os.path.join('results', file)
complete_model_results_df.to_csv(file, sep='\t', index=False)

file = 'ensemble_results_{}_{}.tsv'.format(dataset, geneset)
file = os.path.join('results', file)
complete_ensemble_results_df.to_csv(file, sep='\t', index=False)

file = 'all_results_{}_{}.tsv'.format(dataset, geneset)
file = os.path.join('results', file)
complete_all_results_df.to_csv(file, sep='\t', index=False)

file = 'top_results_{}_{}.tsv'.format(dataset, geneset)
file = os.path.join('results', file)
complete_top_results_df.to_csv(file, sep='\t', index=False)
