
# coding: utf-8

# # Determining Across Z Dimension Stability with SVCCA
# 
# **Gregory Way 2018**
# 
# Here, we apply Singular Vector Canonical Correlation Analysis ([Raghu et al. 2017](https://arxiv.org/abs/1706.05806 "SVCCA: Singular Vector Canonical Correlation Analysis for Deep Learning Dynamics and Interpretability")) ([github](https://github.com/google/svcca)) to the neuron matrices Z to quantify model stability both across algorithms and datasets.
# We test the ability of algorithms to identify stable compressed gene expression features across z dimensions.
# 
# Briefly, SVCCA uses Singular Value Decomposition (SVC) to extract the components explaining 99% of the variation.
# This is done to remove potential dimensions described by noise.
# Next, SVCCA performs a Canonical Correlation Analysis (CCA) on the SVD matrices to identify maximum correlations of linear combinations of both input matrices. The algorithm will identify the canonical correlations of highest magnitude across and within algorithms of the same dimensionality. 
# 
# The output of the SVCCA analysis is the SVCCA mean similarity score. This single number can be interpreted as a measure of similarity, or stability, in the solutions identified from each compression algorithm.
# 
# We perform SVCCA across dimensions, and save the mean estimated SVCCA mean similarity score for each comparison.

# In[1]:


import os
import glob
import numpy as np
import pandas as pd

from scripts.util import read_in_matrices, get_svcca_across_z_stability


# In[2]:


datasets = ['TARGET', 'TCGA', 'GTEX']
algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']
z_dims = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30,
          35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200]


# In[3]:


for dataset in datasets:
    
    dataset_svcca_results_list = []
    for algorithm in algorithms:

        for z_a in z_dims:

            # Read in the first dictionary
            z_dict_a = read_in_matrices(
                dataset=dataset,
                z_dim=z_a,
                algorithm=algorithm,
                shuffled_data=False,
                load_weights=True
            )
            
            # Only compare to higher z dimensions
            z_bs = [x for x in z_dims if x > z_a]

            for z_b in z_bs:

                # Read in the second dictionary
                z_dict_b = read_in_matrices(
                    dataset=dataset,
                    z_dim=z_b,
                    algorithm=algorithm,
                    shuffled_data=False,
                    load_weights=True
                )

                print("Calculating... dataset {}, algorithm {}, and dimension {} vs. {}"
                      .format(dataset, algorithm, z_a, z_b))

                # Perform across z dimension SVCCA
                svcca_out = get_svcca_across_z_stability(
                    z_dict_a=z_dict_a,
                    z_dict_b=z_dict_b,
                    algorithm=algorithm,
                )

                # Append info to the output dataframe
                svcca_out = svcca_out.assign(
                    dataset=dataset,
                    algorithm=algorithm,
                    z_dim_a=z_a,
                    z_dim_b=z_b
                )

                # Append to final list
                dataset_svcca_results_list.append(svcca_out)
                
    # Save dataset specific results
    svcca_results_df = pd.concat(dataset_svcca_results_list)
    svcca_results_df.columns = ['svcca_mean_similarity', 'dataset', 'algorithm',
                                'z_dim_a', 'z_dim_b']

    out_file = os.path.join('results',
                            'svcca_across_z_{}_mean_correlation.tsv.gz'.format(dataset))
    svcca_results_df.to_csv(out_file, sep='\t', index=False, compression='gzip')

