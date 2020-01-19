#!/usr/bin/env python
# coding: utf-8

# # Apply Compression Models to Simulated Data
# 
# Here, we apply the suite of compression models over a range of latent dimensionalities (k) to the simulated data.
# 
# We apply PCA, ICA, NMF, DAE, and VAE models over a range of k (k = 1, 2, 3, 4, 5, 6).
# 
# We extract the weight matrices for every iteration.
# We will show which compressed feature captures the two groups of simulated signals.

# In[1]:


import os
import random
import pandas as pd
from sklearn import decomposition

from tybalt.data_models import DataModel


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


random.seed(123)


# In[4]:


# Setup constants
ks = list(range(1, 7))
data_file = os.path.join("data", "simulated_signal_n1000_p10.tsv")


# In[5]:


data_df = pd.read_csv(data_file, sep='\t')
data_df.index = ["sample_{}".format(x) for x in data_df.index]

print(data_df.shape)
data_df.head()


# In[6]:


# Split into training and testing sets
# (For compatibility with tybalt.DataModel)
split_prop = 0.05
test_samples = random.sample(range(0, data_df.shape[0]), int(data_df.shape[0] * split_prop))

test_df = data_df.iloc[test_samples, :]
train_df = data_df.drop(test_df.index, axis="index")


# In[7]:


# Initialize DataModel class with the input data
dm = DataModel(df=train_df, test_df=test_df)
dm.transform(how='zeroone')


# In[8]:


# Parameters selected to be similar to real data parameter sweep
epochs = 25
batch_size = 50
vae_learning_rate = 0.0015
dae_learning_rate = 0.0005
dae_noise = 0.01
dae_sparsity = 0


# In[9]:


# Loop over the latent dimensionalities
sim_results = list()
for k in ks:
    # Fit models
    # 1) PCA
    dm.pca(n_components=k, transform_test_df=False)
    result = dm.pca_weights.assign(k=k, algorithm="PCA")
    sim_results.append(result)

    # 2) ICA
    dm.ica(n_components=k, transform_test_df=False)
    result = dm.ica_weights.assign(k=k, algorithm="ICA")
    sim_results.append(result)

    # 3) NMF
    dm.nmf(n_components=k, transform_test_df=False)
    result = dm.nmf_weights.assign(k=k, algorithm="NMF")
    sim_results.append(result)

    # 4) DAE
    dm.nn(n_components=k,
          model='adage',
          loss='binary_crossentropy',
          epochs=epochs,
          batch_size=batch_size,
          learning_rate=dae_learning_rate,
          noise=dae_noise,
          sparsity=dae_sparsity,
          verbose=False,
          transform_test_df=False)
    result = dm.adage_weights.assign(k=k, algorithm="DAE")
    sim_results.append(result)

    # 4) VAE
    dm.nn(n_components=k,
          model='tybalt',
          loss='binary_crossentropy',
          epochs=epochs,
          batch_size=batch_size,
          learning_rate=vae_learning_rate,
          separate_loss=False,
          verbose=False,
          transform_test_df=False)
    result = dm.tybalt_weights.assign(k=k, algorithm="VAE")
    sim_results.append(result)


# In[10]:


# Compile and output results
full_sim_results = (
    pd.concat(sim_results)
    .reset_index()
    .rename({"index": "compressed_feature"}, axis="columns")
)

out_file = os.path.join("results", "compression_simulation_results.tsv")
full_sim_results.to_csv(out_file, sep='\t', index=False)

full_sim_results.tail(10)

