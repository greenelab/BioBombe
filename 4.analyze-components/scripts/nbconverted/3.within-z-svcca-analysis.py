
# coding: utf-8

# # Determining Within Z Dimension Stability with SVCCA
# 
# **Gregory Way 2018**
# 
# Here, we apply Singular Vector Canonical Correlation Analysis ([Raghu et al. 2017](https://arxiv.org/abs/1706.05806 "SVCCA: Singular Vector Canonical Correlation Analysis for Deep Learning Dynamics and Interpretability")) ([github](https://github.com/google/svcca)) to the neuron matrices Z to quantify model stability both within and across algorithms over each bottleneck dimensionality.
# 
# Briefly, SVCCA uses Singular Value Decomposition (SVC) to extract the components explaining 99% of the variation.
# This is done to remove potential dimensions described by noise.
# Next, SVCCA performs a Canonical Correlation Analysis (CCA) on the SVD matrices to identify maximum correlations of linear combinations of both input matrices. The algorithm will identify the canonical correlations of highest magnitude across and within algorithms of the same dimensionality. 
# 
# The output of the SVCCA analysis is the SVCCA mean similarity score. This single number can be interpreted as a measure of similarity, or stability, in the solutions identified from each compression algorithm.

# In[1]:


import os
import glob
import numpy as np
import pandas as pd

import scripts.cca_core as cca


# In[2]:


def get_svcca_model_stability(dataset, z_dim, algorithms, shuffled_data=False):
    """
    Compile SVCCA results for all combinations of within algorithm for given dataset and z
    
    Arguments:
    dataset - a string indicating either "TCGA", "TARGET", or "GTEX"
    z_dim - a string indicating the bottleneck dimensionality
    algorithm - a list of string indicating which algorithm to focus on
    shuffled_data - a boolean indicating if the data was first shuffled before training

    Output:
    a list of mean SVCCA similarity scores for each cross comparison
    """
    
    # Build the directory where the results are stored
    base_dir = os.path.join('..', '2.ensemble-z-analysis', 'results')
    
    if shuffled_data:
        results_dir = "shuffled_results"
        shuffled_assign = "shuffled"
    else:
        results_dir = "results"
        shuffled_assign = "signal"

    dataset_dir = os.path.join('{}_{}'.format(dataset, results_dir), 'ensemble_z_matrices')
    z_dim_dir = '{}_components_{}'.format(dataset.lower(), z_dim)
    full_dir = os.path.join(base_dir, dataset_dir, z_dim_dir)

    z_dict = {}
    for file_name in glob.glob('{}/*_z_matrix*'.format(full_dir)):
        seed = os.path.basename(file_name).split('_')[1]
        z_dict[seed] = pd.read_table(file_name, index_col=0)
        
    output_list = []
    for model_a in z_dict.keys():
        model_a_df = z_dict[model_a]
        for model_b in z_dict.keys():
            if model_a != model_b:
                model_b_df = z_dict[model_b]
                for algorithm_a in algorithms:
                    for algorithm_b in algorithms:
                        compile_list = [model_a, model_b, algorithm_a, algorithm_b]

                        z_a = model_a_df.loc[:, model_a_df.columns.str.contains(algorithm_a)]
                        z_b = model_b_df.loc[:, model_b_df.columns.str.contains(algorithm_b)]
                        result = cca.get_cca_similarity(z_a.T, z_b.T, verbose=False)

                        compile_list += [np.mean(result['mean'])]

                        output_list.append(compile_list)
    
    output_df = pd.DataFrame(output_list)
    output_df = output_df.assign(dataset=dataset,
                                 z_dim=z_dim,
                                 shuffled=shuffled_assign)

    return output_df


# In[3]:


datasets = ['TARGET', 'TCGA', 'GTEX']
algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']
z_dims = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30,
          35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200]
signals = (True, False)


# In[4]:


large_svcca_results_list = []
for dataset in datasets:
    for z in z_dims:
        for signal in signals:
            print("Calculating... dataset {} for {} dimension {}".format(dataset, signal, z))
            svcca_out = get_svcca_model_stability(dataset=dataset,
                                                  z_dim=z,
                                                  algorithms=algorithms,
                                                  shuffled_data=signal)
            large_svcca_results_list.append(svcca_out)


# In[5]:


svcca_results_df = pd.concat(large_svcca_results_list)
svcca_results_df.columns = ['seed_1', 'seed_2', 'algorithm_1', 'algorithm_2',
                            'svcca_mean_similarity', 'dataset', 'z_dim', 'shuffled']
svcca_results_df.head()


# In[6]:


out_file = os.path.join('results', 'svcca_mean_correlation_within_z.tsv.gz')
svcca_results_df.to_csv(out_file, sep='\t', index=False, compression='gzip')

