"""
2018 Gregory Way
5.analyze-components/geneset_tracking.py

Apply interpret-compression approach to weight matrices

Usage:

    python geneset_tracking.py

    with the following required flags:

        --dataset       string of dataset of interest (TCGA, TARGET, or GTEX)
        --metaedge      the metaedge to build the adjacency matrix using

    and the following optional flags:

        --subset_num_genes  the number of genes to subset (default: 8000)
        --shuffled          decision to use real data or shuffled data

Output:
A single long dataframe storing the values and z_scores of all input gene set
analysis results for the given input dataset across all weight matrices for
all bottleneck dimensions and all compression algorithms.
"""

import os
import sys
import glob
import argparse
import pandas as pd

sys.path.append('../scripts')
from latent import latentModel, load_hetnets

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dataset', choices=['TCGA', 'GTEX', 'TARGET'],
                    help='the dataset used to explore results')
parser.add_argument('-m', '--metaedge',
                    help='the metaedge to to construct adjacency matrix')
parser.add_argument('-n', '--subset_num_genes', default=8000,
                    help='subset num genes based on mean absolute deviation')
parser.add_argument('-s', '--shuffled', action='store_true',
                    help='compile results with shuffled data or not')
args = parser.parse_args()

# Load command arguments
dataset = args.dataset
metaedge = args.metaedge
subset_mad_genes = int(args.subset_num_genes)
shuffled = args.shuffled

output_dir = os.path.join('results', dataset.lower(), metaedge.lower())

if shuffled:
    output_dir = os.path.join(output_dir, 'shuffled')
else:
    output_dir = os.path.join(output_dir, 'signal')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load data
file = '{}_mad_genes.tsv'.format(dataset.lower())
file = os.path.join('..', '0.expression-download', 'data', file)
genes_df = pd.read_table(file)
genes = genes_df.iloc[0:subset_mad_genes, ].gene_id.astype(str)

# Load hetnets for the given metaedge for rapid latent feature interpretation
hetnets = load_hetnets(
    hetnet_file='../3.build-hetnets/hetnets/interpret_hetnet.json.bz2',
    permuted_directory='../3.build-hetnets/hetnets/permuted/',
    subset_genes=genes,
    metaedge_abbrev=metaedge
   )

# Load weight matrices
if shuffled:
    results = '{}_shuffled_results'.format(dataset)
else:
    results = '{}_results'.format(dataset)

base_dir = os.path.join('..', '2.ensemble-z-analysis', 'results',
                        results, 'ensemble_z_matrices')
weight_matrices = {}
for path in os.listdir(base_dir):
    z_dim = path.split('_')[-1]
    weight_matrices[z_dim] = {}
    z_dir = os.path.join(base_dir, path)
    for file_name in glob.glob('{}/*_weight_matrix*'.format(z_dir)):
        seed = file_name.split('/')[-1].split('_')[1]
        weight_matrices[z_dim][seed] = file_name

# Extract results from the weight weight matrices
for z_dim in weight_matrices.keys():
    print('processing dimension {}...'.format(z_dim))
    weight_matrix_seed_dict = weight_matrices[z_dim]

    # Build output file
    if shuffled:
        out_base = '_geneset_scores_shuffled.tsv.gz'
    else:
        out_base = '_geneset_scores.tsv.gz'

    file = '{}_z_{}_{}_{}'.format(dataset.lower(), z_dim, metaedge, out_base)
    file = os.path.join(output_dir, file)

    weight_matrix_seed_results = []
    for seed in weight_matrix_seed_dict.keys():
        print('     for seed {}'.format(seed))
        weight_matrix = weight_matrix_seed_dict[seed]

        # Load the latent model
        lm = latentModel(
            filename=weight_matrix,
            z_dim=z_dim,
            dataset_name=dataset.lower(),
            algorithm_name=None,
            weight_seed=seed,
            shuffled_true=False
            )

        # Perform the matrix multiplication
        mult_results = {}
        all_list = []
        for model in hetnets.keys():
            print('        and model {}'.format(model))
            hetnet = hetnets[model]
            mult_results[model] = lm.w_df.T @ hetnet
            long_result = (
                mult_results[model]
                .reset_index()
                .melt(id_vars=['index'])
                )
            long_result = long_result.assign(model=model)

            if model != 'real':
                long_result = long_result.assign(model_type='permuted')
            else:
                long_result = long_result.assign(model_type='real')

            all_list.append(long_result)

        all_df = pd.concat(all_list)
        all_df.value = all_df.value.astype(float)

        # Obtain the z-scores for the real data against permutations
        permuted_mean = (
            all_df
            .groupby(['model_type', 'index', 'variable'])
            .mean()
            .reset_index()
            .query("model_type == 'permuted'")
        )

        permuted_std = (
            all_df
            .groupby(['model_type', 'index', 'variable'])
            .std()
            .reset_index()
            .query("model_type == 'permuted'")
        )

        real_df = (
            all_df
            .groupby(['model_type', 'index', 'variable'])
            .mean()
            .reset_index()
            .query("model_type == 'real'")
        )

        z_score = (
            real_df.reset_index(drop=True).value - permuted_mean.value
            ) / permuted_std.value
        real_df = real_df.reset_index(drop=True).assign(z_score=z_score)

        # Process matrix for output
        algorithm_feature_df = real_df['index'].str.split('_', expand=True)
        algorithm_feature_df.columns = ['algorithm', 'feature']

        real_df = pd.concat([real_df, algorithm_feature_df], axis='columns')
        real_df = real_df.drop(['index'], axis='columns')
        real_df = real_df.assign(
            z=z_dim,
            seed=seed
            )
        weight_matrix_seed_results.append(real_df)

    weight_seed_df = pd.concat(weight_matrix_seed_results)
    weight_seed_df.value = weight_seed_df.value.astype(float)
    (
        weight_seed_df
        .sort_values(by='z_score')
        .to_csv(file, sep='\t', index=False, compression='gzip',
                float_format='%.6g')
    )
