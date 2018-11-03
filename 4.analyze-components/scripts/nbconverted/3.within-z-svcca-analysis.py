
# coding: utf-8

# # Determining Within Z Dimension Stability with SVCCA
# 
# **Gregory Way 2018**
# 
# Here, we apply Singular Vector Canonical Correlation Analysis ([Raghu et al. 2017](https://arxiv.org/abs/1706.05806 "SVCCA: Singular Vector Canonical Correlation Analysis for Deep Learning Dynamics and Interpretability")) ([github](https://github.com/google/svcca)) to the neuron matrices Z to quantify model stability both within and across algorithms over each bottleneck dimensionality.
# 
# Briefly, SVCCA uses Singular Value Decomposition (SVD) to extract the components explaining 99% of the variation.
# This is done to remove potential dimensions described by noise.
# Next, SVCCA performs a Canonical Correlation Analysis (CCA) on the SVD matrices to identify maximum correlations of linear combinations of both input matrices. The algorithm will identify the canonical correlations of highest magnitude across and within algorithms of the same dimensionality. 
# 
# The output of the SVCCA analysis is the SVCCA mean similarity score. This single number can be interpreted as a measure of similarity, or stability, in the solutions identified from each compression algorithm.

# In[1]:


import os
import glob
import numpy as np
import pandas as pd

from scripts.util import read_in_z, get_svcca_across_algorithm_stability


# In[2]:


datasets = ['TARGET', 'TCGA', 'GTEX']
algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']
z_dims = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30,
          35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200]
shuffled_data = (True, False)


# In[3]:


large_svcca_results_list = []
for dataset in datasets:

    for z in z_dims:

        for signal in shuffled_data:

            if signal:
                shuffled_status = 'shuffled'
            else:
                shuffled_status = 'signal'

            print("Calculating... dataset {} for {} dimension {}"
                  .format(dataset, shuffled_status, z))

            # Read in the z matrix of interest
            z_dict = read_in_z(
                dataset=dataset,
                z_dim=z,
                algorithm='all',
                shuffled_data=signal
            )

            # Perform across algorithm SVCCA
            # For training
            svcca_out_train_df = (
                get_svcca_across_algorithm_stability(z_dict=z_dict,
                                                     algorithms=algorithms,
                                                     train_or_test='train')
            )
            svcca_out_train_df = svcca_out_train_df.assign(train_or_test='train')

            # SVCCA cannot be calculated for TARGET test data because there are
            # not enough datapoints
            if dataset != 'TARGET':
                # And for testing
                svcca_out_test_df = (
                    get_svcca_across_algorithm_stability(z_dict=z_dict,
                                                         algorithms=algorithms,
                                                         train_or_test='test')
                )
                svcca_out_test_df = svcca_out_test_df.assign(train_or_test='test')

                # Concatenate training and testing
                svcca_out_df = pd.concat([svcca_out_train_df, svcca_out_test_df])
            else:
                svcca_out_df = svcca_out_train_df

            # Append info to the output dataframe
            svcca_out_df = svcca_out_df.assign(
                dataset=dataset,
                z_dim=z,
                shuffled=shuffled_status
            )

            large_svcca_results_list.append(svcca_out_df)


# In[4]:


svcca_results_df = pd.concat(large_svcca_results_list)
svcca_results_df.columns = ['seed_1', 'seed_2', 'algorithm_1', 'algorithm_2',
                            'svcca_mean_similarity', 'train_or_test',
                            'dataset', 'z_dim', 'shuffled']
svcca_results_df.head()


# In[5]:


out_file = os.path.join('results', 'svcca_mean_correlation_within_z.tsv.gz')
svcca_results_df.to_csv(out_file, sep='\t', index=False, compression='gzip')

