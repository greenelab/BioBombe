"""
2018 Gregory Way
2.ensemble-z-analysis/scripts/train_models_given_z.py

This script will train various compression models given a specific z dimension.
Each model will train several times with different initializations.

The script pulls hyperparameters from a parameter file that was determined
after initial hyperparameter sweeps testing latent dimensionality.

Usage:

    python train_models_given_z.py

    With required command line arguments:

        --dataset           The focus dataset to apply compression models to
        --num_components    The z dimensionality we're testing
        --param_config      A tsv file (param by z dimension) indicating the
                            specific parameter combination for the z dimension
        --out_dir           The directory to store the output files

    And optional command line arguments

        --num_seeds         The number of specific models to generate
                                default: 5
        --shuffle           If provided, shuffle the input expression matrix

Output:

The script will save associated weights and z matrices for each permutation as
well as reconstruction costs for all algorithms, sample specific correlations,
and training histories for Tybalt and ADAGE models. The population of models
(ensemble) are saved, across dimensions z, for downstream evaluation.
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from tybalt.data_models import DataModel

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dataset', choices=['TCGA', 'GTEX', 'TARGET'],
                    help='the dataset used for compression')
parser.add_argument('-n', '--num_components', help='dimensionality of z')
parser.add_argument('-p', '--param_config',
                    help='text file optimal hyperparameter assignment for z')
parser.add_argument('-o', '--out_dir', help='where to save the output files')
parser.add_argument('-s', '--num_seeds', default=5,
                    help='number of different seeds to run on current data')
parser.add_argument('-r', '--shuffle', action='store_true',
                    help='randomize gene expression data for negative control')
parser.add_argument('-m', '--subset_mad_genes', default=8000,
                    help='subset num genes based on mean absolute deviation')
args = parser.parse_args()

# Load command arguments
dataset = args.dataset.lower()
num_components = int(args.num_components)
param_config = args.param_config
out_dir = args.out_dir
num_seeds = int(args.num_seeds)
shuffle = args.shuffle
subset_mad_genes = int(args.subset_mad_genes)


def get_recon_correlation(df, recon_mat_dict, algorithm, cor_type,
                          genes=False):
    """
    Get gene or sample correlations between input and reconstructed input

    Arguments:
    df - the input dataframe
    recon_mat_dict - dictionary of different algorithms reconstructions
    algorithm - string representing the compression algorithm
    cor_type - string representing Pearson or Spearman correlations
    genes - boolean if to calculate correaltion over genes (sample by default)
    """
    recon_mat = recon_mat_dict[algorithm]
    if genes:
        df = df.T
        recon_mat = recon_mat.T
    if cor_type == 'pearson':
        r = [pearsonr(recon_mat.iloc[x, :],
                      df.iloc[x, :])[0] for x in range(df.shape[0])]
    elif cor_type == 'spearman':
        r = [spearmanr(recon_mat.iloc[x, :],
                       df.iloc[x, :])[0] for x in range(df.shape[0])]
    return r


def compile_corr_df(pearson_list, spearman_list, algorithm_list, column_names,
                    seed, data_type):
    """
    Compile together correlations across algorithms

    Arguments:
    pearson_list - a list of pearson correlations across algorithms
    spearman_list - a list of spearman correlations across algorithms
    algorithm_list - list of algorithm names
    column_names - list of names supplied to the compiled dataframe
    seed - the current random seed
    data_type - training or testing set
    """
    pearson_df = pd.DataFrame(pearson_list,
                              index=algorithm_list,
                              columns=column_names)
    pearson_df.index.name = 'algorithm'
    spearman_df = pd.DataFrame(spearman_list,
                               index=algorithm_list,
                               columns=column_names)
    spearman_df.index.name = 'algorithm'
    pearson_df = pearson_df.reset_index().melt(id_vars=['algorithm'],
                                               var_name='id',
                                               value_name='correlation')
    spearman_df = spearman_df.reset_index().melt(id_vars=['algorithm'],
                                                 var_name='id',
                                                 value_name='correlation')
    corr_df = pd.concat([pearson_df.assign(cor_type='pearson'),
                         spearman_df.assign(cor_type='spearman')])
    corr_df = corr_df.assign(seed=seed)
    corr_df = corr_df.assign(data=data_type)
    return corr_df

# Extract parameters from parameter configuation file
param_df = pd.read_table(param_config, index_col=0)

component_error = ' '.join(str(x) for x in
                           [num_components, 'is not found in', param_config,
                            '- either add it to the file or choose a',
                            'different number of components'])
assert str(num_components) in param_df.columns, component_error

vae_epochs = param_df.loc['vae_epochs', str(num_components)]
dae_epochs = param_df.loc['dae_epochs', str(num_components)]
vae_lr = param_df.loc['vae_lr', str(num_components)]
dae_lr = param_df.loc['dae_lr', str(num_components)]
vae_batch_size = param_df.loc['vae_batch_size', str(num_components)]
dae_batch_size = param_df.loc['dae_batch_size', str(num_components)]
dae_noise = param_df.loc['dae_noise', str(num_components)]
dae_sparsity = param_df.loc['dae_sparsity', str(num_components)]
vae_kappa = param_df.loc['vae_kappa', str(num_components)]

# Set output directory and file names
train_dir = os.path.join(out_dir, 'ensemble_z_results',
                         '{}_components'.format(num_components))
if not os.path.exists(train_dir):
    os.makedirs(train_dir)

if shuffle:
    file_pre = '{}_{}_components_shuffled_'.format(dataset, num_components)
else:
    file_pre = '{}_{}_components_'.format(dataset, num_components)

recon_file = os.path.join(train_dir, '{}reconstruction.tsv'.format(file_pre))
co_file = os.path.join(train_dir, '{}sample_corr.tsv.gz'.format(file_pre))
co_g_file = os.path.join(train_dir, '{}gene_corr.tsv.gz'.format(file_pre))
tybalt_hist_file = os.path.join(train_dir,
                                '{}tybalt_training_hist.tsv'.format(file_pre))
adage_hist_file = os.path.join(train_dir,
                               '{}adage_training_hist.tsv'.format(file_pre))

# Load Data
base_dir = os.path.join('..', '0.expression-download', 'data')
rnaseq_train = (
    os.path.join(base_dir,
    'train_{}_expression_matrix_processed.tsv.gz'.format(dataset))
    )
rnaseq_test = (
    os.path.join(base_dir,
    'test_{}_expression_matrix_processed.tsv.gz'.format(dataset))
    )

rnaseq_train_df = pd.read_table(rnaseq_train, index_col=0)
rnaseq_test_df = pd.read_table(rnaseq_test, index_col=0)

# Determine most variably expressed genes and subset
if subset_mad_genes is not None:
    mad_genes = rnaseq_train_df.mad(axis=0).sort_values(ascending=False)
    top_mad_genes = mad_genes.iloc[0:subset_mad_genes, ].index
    rnaseq_train_df = rnaseq_train_df.reindex(top_mad_genes, axis='columns')
    rnaseq_test_df = rnaseq_test_df.reindex(top_mad_genes, axis='columns')

# Initialize DataModel class with pancancer data
dm = DataModel(df=rnaseq_train_df, test_df=rnaseq_test_df)
dm.transform(how='zeroone')

# Set seed and list of algorithms for compression
random_seeds = np.random.randint(0, high=1000000, size=num_seeds)
algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']

# Save population of models in specific folder
comp_out_dir = os.path.join(out_dir, 'ensemble_z_matrices',
                            '{}_components_{}'.format(dataset, num_components))
if not os.path.exists(comp_out_dir):
    os.makedirs(comp_out_dir)

reconstruction_results = []
test_reconstruction_results = []
sample_correlation_results = []
tybalt_training_histories = []
adage_training_histories = []
for seed in random_seeds:
    seed_file = os.path.join(comp_out_dir, 'model_{}'.format(seed))

    if shuffle:
        seed_file = '{}_shuffled'.format(seed_file)

        # randomly permute genes of each sample in the rnaseq matrix
        shuf_df = rnaseq_train_df.apply(lambda x:
                                        np.random.permutation(x.tolist()),
                                        axis=1)

        # Initiailze a new DataModel, with different shuffling each permutation
        dm = DataModel(df=shuf_df, test_df=rnaseq_test_df)
        dm.transform(how='zeroone')

    # Fit models
    dm.pca(n_components=num_components, transform_test_df=True)
    dm.ica(n_components=num_components, transform_test_df=True)
    dm.nmf(n_components=num_components, transform_test_df=True)

    dm.nn(n_components=num_components,
          model='tybalt',
          loss='binary_crossentropy',
          epochs=int(vae_epochs),
          batch_size=int(vae_batch_size),
          learning_rate=float(vae_lr),
          separate_loss=True,
          verbose=True,
          transform_test_df=True)

    dm.nn(n_components=num_components,
          model='adage',
          loss='binary_crossentropy',
          epochs=int(dae_epochs),
          batch_size=int(dae_batch_size),
          learning_rate=float(dae_lr),
          noise=float(dae_noise),
          sparsity=float(dae_sparsity),
          verbose=True,
          transform_test_df=True)

    # Obtain z matrix (sample scores per latent space feature) for all models
    full_z_file = '{}_z_matrix.tsv.gz'.format(seed_file)
    dm.combine_models().to_csv(full_z_file, sep='\t', compression='gzip')

    full_test_z_file = '{}_z_test_matrix.tsv.gz'.format(seed_file)
    dm.combine_models(test_set=True).to_csv(full_test_z_file, sep='\t',
                                            compression='gzip')

    # Obtain weight matrices (gene by latent space feature) for all models
    full_weight_file = '{}_weight_matrix.tsv.gz'.format(seed_file)
    dm.combine_weight_matrix().to_csv(full_weight_file, sep='\t',
                                      compression='gzip')

    # Store reconstruction costs and reconstructed input at training end
    full_reconstruction, reconstructed_matrices = dm.compile_reconstruction()

    # Store reconstruction evaluation and data for test set
    full_test_recon, test_recon_mat = dm.compile_reconstruction(test_set=True)

    # Get correlations across samples and genes between input and output data
    pearson_corr = []
    spearman_corr = []
    pearson_corr_test = []
    spearman_corr_test = []

    for algorithm in algorithms:
        # Training Sample Correlations
        pearson_corr.append(
            get_recon_correlation(df=dm.df,
                                  recon_mat_dict=reconstructed_matrices,
                                  algorithm=algorithm,
                                  cor_type='pearson',
                                  genes=False)
                            )
        spearman_corr.append(
            get_recon_correlation(df=dm.df,
                                  recon_mat_dict=reconstructed_matrices,
                                  algorithm=algorithm,
                                  cor_type='spearman',
                                  genes=False)
                            )

        # Testing Sample Correlations
        pearson_corr_test.append(
            get_recon_correlation(df=dm.test_df,
                                  recon_mat_dict=test_recon_mat,
                                  algorithm=algorithm,
                                  cor_type='pearson',
                                  genes=False)
                                )
        spearman_corr_test.append(
            get_recon_correlation(df=dm.test_df,
                                  recon_mat_dict=test_recon_mat,
                                  algorithm=algorithm,
                                  cor_type='spearman',
                                  genes=False)
                                 )

    # Training - Sample correlations between input and reconstruction
    sample_correlation_results.append(
        compile_corr_df(pearson_list=pearson_corr,
                        spearman_list=spearman_corr,
                        algorithm_list=algorithms,
                        column_names=dm.df.index,
                        seed=seed,
                        data_type='training')
                        )

    # Testing - Sample correlations between input and reconstruction
    sample_correlation_results.append(
        compile_corr_df(pearson_list=pearson_corr_test,
                        spearman_list=spearman_corr_test,
                        algorithm_list=algorithms,
                        column_names=dm.test_df.index,
                        seed=seed,
                        data_type='testing')
                        )

    # Store training histories and intermediate results for neural networks
    reconstruction_results.append(
        full_reconstruction.assign(seed=seed, shuffled=shuffle)
        )
    test_reconstruction_results.append(
        full_test_recon.assign(seed=seed, shuffled=shuffle)
        )
    tybalt_training_histories.append(
        dm.tybalt_fit.history_df.assign(seed=seed, shuffle=shuffle)
        )
    adage_training_histories.append(
        dm.adage_fit.history_df.assign(seed=seed, shuffle=shuffle)
        )

# Save reconstruction and neural network training results
pd.concat([
    pd.concat(reconstruction_results).assign(data_type='training'),
    pd.concat(test_reconstruction_results).assign(data_type='testing')
    ]).reset_index(drop=True).to_csv(recon_file, sep='\t')
pd.concat(sample_correlation_results).to_csv(co_file,
                                             sep='\t',
                                             index=False,
                                             float_format='%.3f',
                                             compression='gzip')
pd.concat(tybalt_training_histories).to_csv(tybalt_hist_file, sep='\t')
pd.concat(adage_training_histories).to_csv(adage_hist_file, sep='\t')
