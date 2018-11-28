"""
Interpretation of Compression Models
Gregory Way 2018

Storing helper functions to enable SVCCA analysis and loading data

Usage:
Import Only
"""

import os
import glob
import numpy as np
import pandas as pd

import scripts.cca_core as cca


def read_in_matrices(dataset, z_dim, algorithm, shuffled_data=False,
                     load_weights=False):
    """
    Read in dataset, algorithm, and dimensionality specific z or weight matrix

    Arguments:
    dataset - a string indicating either "TCGA", "TARGET", or "GTEX"
    z_dim - a string indicating the bottleneck dimensionality
    algorithm - a string indicating which algorithm to focus on (can be 'all')
    shuffled_data - a boolean if the data was first shuffled before training
    read_weights - a boolean if the weight matrices should be loaded instead

    Output:
    a list of mean SVCCA similarity scores for each cross comparison
    """

    # What string to glob for and what output to expect
    if load_weights:
        lookfor = '_weight_'
        out_dict = {}
    else:
        lookfor = '_z_'
        out_dict = {'train': {}, 'test': {}}

    # Build the directory where the results are stored
    base_dir = os.path.join('..', '2.ensemble-z-analysis', 'results')

    if shuffled_data:
        results_dir = "shuffled_results"
    else:
        results_dir = "results"

    dataset_dir = os.path.join('{}_{}'.format(dataset, results_dir),
                               'ensemble_z_matrices')
    z_dim_dir = '{}_components_{}'.format(dataset.lower(), z_dim)
    full_dir = os.path.join(base_dir, dataset_dir, z_dim_dir)

    for file_name in glob.glob('{}/*{}*'.format(full_dir, lookfor)):
        file = os.path.basename(file_name)
        seed = os.path.basename(file_name).split('_')[1]
        df = pd.read_table(file_name, index_col=0)

        if algorithm != 'all':
            df = df.loc[:, df.columns.str.contains(algorithm)]

        if load_weights:
            out_dict[seed] = df
        else:
            if 'test' in file:
                out_key = 'test'
            else:
                out_key = 'train'

            out_dict[out_key][seed] = df.assign(train_or_test=out_key)

    return out_dict


def get_svcca_across_algorithm_stability(z_dict,
                                         algorithms,
                                         train_or_test=None,
                                         threshold=0.98):
    """
    Compile SVCCA results for all combinations of within algorithms for a given
    dataset and z

    Arguments:
    z_dict - a dictionary storing the specific z dataframes
    algorithms - a list storing the algorithms to compare
    train_or_test - a string that indicates if the comparison should be
                    training or testing
    threshold - remove trailing zeros if applicable (for cca_similarity)

    Output:
    a list of mean SVCCA similarity scores for each cross comparison
    """

    if train_or_test:
        z_dict = z_dict[train_or_test]

    output_list = []
    for model_a in z_dict.keys():
        model_a_df = z_dict[model_a]
        for model_b in z_dict.keys():
            if model_a != model_b:
                model_b_df = z_dict[model_b]
                for algorithm_a in algorithms:
                    for algorithm_b in algorithms:
                        compile_list = [model_a, model_b, algorithm_a,
                                        algorithm_b]

                        a_subset = model_a_df.columns.str.contains(algorithm_a)
                        b_subset = model_b_df.columns.str.contains(algorithm_b)

                        z_a = model_a_df.loc[:, a_subset]
                        z_b = model_b_df.loc[:, b_subset]

                        result = cca.robust_cca_similarity(z_a.T,
                                                           z_b.T,
                                                           verbose=False,
                                                           threshold=threshold)

                        compile_list += [np.mean(result['mean'])]

                        output_list.append(compile_list)

    # Convert output to pandas dataframe
    out_cols = ['seed_1', 'seed_2', 'algorithm_1', 'algorithm_2',
                'svcca_mean_similarity']
    return pd.DataFrame(output_list, columns=out_cols)


def get_svcca_across_z_stability(z_dict_a,
                                 z_dict_b,
                                 algorithm,
                                 train_or_test=None,
                                 threshold=0.98):
    """
    Given two dictionaries, assess model stability between entries and output.
    The a dictionary should always contain fewer dimensions than b

    Arguments:
    z_dict_a - a dict of the output of `read_in_z` for specific z and dataset
    z_dict_b - the same structure but a different z and dataset than z_dict_a
    train_or_test - a string that indicates if the comparison should be
                    training or testing (or None if weight matrix)
    threshold - remove trailing zeros if applicable (for cca_similarity)

    Output:
    A dataframe of SVCCA similarity scores comparing all models
    """

    if train_or_test:
        z_dict_a = z_dict_a[train_or_test]
        z_dict_b = z_dict_b[train_or_test]

    output_list = []
    for model_a in z_dict_a.keys():
        model_a_df = z_dict_a[model_a]
        for model_b in z_dict_b.keys():
            if model_a != model_b:
                model_b_df = z_dict_b[model_b]

                a_subset = model_a_df.columns.str.contains(algorithm)
                b_subset = model_b_df.columns.str.contains(algorithm)

                z_a = model_a_df.loc[:, a_subset]
                z_b = model_b_df.loc[:, b_subset]

                result = cca.robust_cca_similarity(z_a.T,
                                                   z_b.T,
                                                   verbose=False,
                                                   threshold=threshold)

                output_list.append(np.mean(result['cca_coef1']))

    # Convert output to pandas dataframe
    return pd.DataFrame(output_list, columns=['svcca_mean_similarity'])
