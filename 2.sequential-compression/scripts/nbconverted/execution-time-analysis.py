#!/usr/bin/env python
# coding: utf-8

# ## Testing the Execution Time of Compression Algorithms
# 
# **Gregory Way, 2019**
# 
# Observing the execution time of algorithms and latent dimensionalities across datasets.

# In[1]:


import os
import pandas as pd
import tensorflow as tf
import time

from tybalt.data_models import DataModel


# In[2]:


tf.logging.set_verbosity(tf.logging.ERROR)


# In[3]:


# Setup constants
datasets = ["TARGET", "TCGA", "GTEX"]
ks = [2, 4, 10, 16, 25, 50, 80, 200]
subset_mad_genes = 8000
data_dir = os.path.join('..', '0.expression-download', 'data')


# In[ ]:


time_results = []
for dataset in datasets:
    
    # Setup constants per dataset
    dataset_id = dataset.lower()
    param_config = os.path.join("config", "z_parameter_sweep_{}.tsv".format(dataset))

    # For extracting parameters from optimized parameter configuation file
    param_df = pd.read_table(param_config, index_col=0)

    # Load MAD genes per dataset
    mad_file = os.path.join(data_dir, '{}_mad_genes.tsv'.format(dataset_id))
    mad_genes_df = pd.read_table(mad_file)
    mad_genes = mad_genes_df.iloc[0:subset_mad_genes, ].gene_id.astype(str)
    
    # Load input data per dataset and reindex to mad genes
    train_file = os.path.join(data_dir,
                              'train_{}_expression_matrix_processed.tsv.gz'.format(dataset_id))
    test_file = os.path.join(data_dir,
                             'test_{}_expression_matrix_processed.tsv.gz'.format(dataset_id))
    rnaseq_train_df = pd.read_table(train_file, index_col=0).reindex(mad_genes, axis='columns')
    rnaseq_test_df = pd.read_table(test_file, index_col=0).reindex(mad_genes, axis='columns')
    
    # Initialize DataModel class with the input data
    dm = DataModel(df=rnaseq_train_df, test_df=rnaseq_test_df)
    dm.transform(how='zeroone')
    
    # Loop over the latent dimensionalities
    for k in ks:
        print("Timing {}: k = {}".format(dataset, k))

        # Retrieve optimized parameters for neural network models
        vae_epochs = param_df.loc['vae_epochs', str(k)]
        dae_epochs = param_df.loc['dae_epochs', str(k)]
        vae_lr = param_df.loc['vae_lr', str(k)]
        dae_lr = param_df.loc['dae_lr', str(k)]
        vae_batch_size = param_df.loc['vae_batch_size', str(k)]
        dae_batch_size = param_df.loc['dae_batch_size', str(k)]
        dae_noise = param_df.loc['dae_noise', str(k)]
        dae_sparsity = param_df.loc['dae_sparsity', str(k)]
        vae_kappa = param_df.loc['vae_kappa', str(k)]
        
        # Fit models
        # 1) PCA
        start = time.time()
        dm.pca(n_components=k, transform_test_df=False)
        end = time.time()
        total_time = end - start
        
        result = [dataset, k, "PCA", total_time]
        time_results.append(result)
        
        # 2) ICA
        start = time.time()
        dm.ica(n_components=k, transform_test_df=False)
        end = time.time()
        total_time = end - start
        
        result = [dataset, k, "ICA", total_time]
        time_results.append(result)
        
        # 3) NMF
        start = time.time()
        dm.nmf(n_components=k, transform_test_df=False)
        end = time.time()
        total_time = end - start
        
        result = [dataset, k, "NMF", total_time]
        time_results.append(result)
        
        # 4) DAE
        start = time.time()
        dm.nn(n_components=k,
              model='adage',
              loss='binary_crossentropy',
              epochs=int(dae_epochs),
              batch_size=int(dae_batch_size),
              learning_rate=float(dae_lr),
              noise=float(dae_noise),
              sparsity=float(dae_sparsity),
              verbose=False,
              transform_test_df=False)
        end = time.time()
        total_time = end - start
        
        result = [dataset, k, "DAE", total_time]
        time_results.append(result)

        # 4) VAE
        start = time.time()
        dm.nn(n_components=k,
              model='tybalt',
              loss='binary_crossentropy',
              epochs=int(vae_epochs),
              batch_size=int(vae_batch_size),
              learning_rate=float(vae_lr),
              separate_loss=True,
              verbose=False,
              transform_test_df=False)
        end = time.time()
        total_time = end - start
        
        result = [dataset, k, "VAE", total_time]
        time_results.append(result)


# In[ ]:


# Combine and output the time analysis results
time_results_df = (
    pd.DataFrame(time_results,
                 columns=["dataset", "k", "algorithm", "seconds"])
    .sort_values("seconds", ascending=False)
)

time_file = os.path.join("results", "time_analysis_results.tsv")
time_results_df.to_csv(time_file, sep='\t', index=False)

print(time_results_df.shape)
time_results_df.head()

