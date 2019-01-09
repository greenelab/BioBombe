"""
2018 Gregory Way
8.gtex-interpret/utils.py

Various methods for the application of specific weight vectors to alternative
gene expression datasets

Usage:
    Import only
"""

import os
import pandas as pd
from urllib.request import urlretrieve


def download_geo(base_url, name, directory):
    """
    Download a dataset from GEO and save to file
    """

    os.makedirs(directory, exist_ok=True)

    path = os.path.join(directory, name)
    urlretrieve(base_url + name, path)


def load_weight_matrix(dataset, z_dim, seed, shuffled=False):
    """
    Load a given weight matrix from file
    """

    if shuffled:
        results_dir = '{}_shuffled_results'.format(dataset)
    else:
        results_dir = '{}_results'.format(dataset)

    base_dir = os.path.join('..', '2.ensemble-z-analysis', 'results',
                            results_dir, 'ensemble_z_matrices')
    z_dir = '{}_components_{}'.format(dataset.lower(), str(z_dim))
    weight_file = 'model_{}_weight_matrix.tsv.gz'.format(str(seed))

    full_file = os.path.join(base_dir, z_dir, weight_file)
    weight_df = pd.read_table(full_file, index_col=0)
    weight_df.index = weight_df.index.map(str)

    return weight_df


def load_enrichment_results(dataset, z_dim, metaedge, algorithm=None,
                            feature=None, seed=None, shuffled=False):
    """
    Load enrichment results for the given dimension, algorithm, and seed
    """

    dataset = dataset.lower()

    if shuffled:
        signal_dir = 'shuffled'
    else:
        signal_dir = 'signal'

    base_file = '{}_z_{}_{}__geneset_scores.tsv.gz'.format(dataset,
                                                           z_dim,
                                                           metaedge)
    base_dir = os.path.join('..', '6.analyze-weights', 'results', dataset,
                            metaedge.lower(), signal_dir)
    enr_file = os.path.join(base_dir, base_file)
    enr_df = pd.read_table(enr_file)

    if seed:
        enr_df = enr_df.query("seed == @seed")

    if algorithm:
        enr_df = enr_df.query("algorithm == @algorithm")

    if feature:
        enr_df = enr_df.query("feature == @feature")

    return enr_df


def apply_signature(weight_df, other_df, feature, align=False):
    """
    Apply a signature to alternative datasets
    """

    if align:
        other_df = other_df.reindex(weight_df.index, axis='columns').fillna(0)

    signature_df = pd.DataFrame(weight_df.loc[:, feature])
    result_df = other_df @ signature_df

    return result_df
