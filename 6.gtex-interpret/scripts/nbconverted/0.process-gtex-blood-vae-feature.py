
# coding: utf-8

# # Capturing Blood Tissue in GTEx Data
# 
# **Gregory Way, 2018**
# 
# When viewing sample correlation differences across z dimensions stratified by tissue-type, we often observed a rapid increase in correlation after increasing model capacity by one.
# 
# For example, the ability of variational autoencoders to capture blood tissue correlation jumps by nearly 50% between bottleneck dimensions 2 and 3 (see below).
# 
# ![sample-correlation_Blood_GTEX_signal_pearson.png](attachment:sample-correlation_Blood_GTEX_signal_pearson.png)
# 
# ## Procedure
# 
# In the following notebook, we extract two representative weight matrices for VAE latent space dimensions 2 and 3.
# We apply our matrix interpretation approach to the weight vectors for each latent space feature.
# We also apply GSEA.
# 
# In both approaches we use genesets derived in the XCELL paper that represent cell-types ([Aran et al. 2017](https://doi.org/10.1186/s13059-017-1349-1))
# 
# We output the results of both approaches and analyze the results in subsequent notebooks.

# In[1]:


import os
import glob
import requests
import pandas as pd
import numpy as np
import gseapy as gp

from scripts.latent import latentModel, load_hetnets, parse_gmt, run_gsea_prerank


# In[2]:


np.random.seed(123)


# In[3]:


metaedge = 'GpXCELL'


# In[4]:


base_dir = '../2.ensemble-z-analysis/results/GTEX_results/ensemble_z_matrices/'


# In[5]:


sample_weight_2 = os.path.join(base_dir, 'gtex_components_2', 'model_109635_weight_matrix.tsv.gz')
sample_weight_3 = os.path.join(base_dir, 'gtex_components_3', 'model_174930_weight_matrix.tsv.gz')


# In[6]:


sample_weight_2_df = pd.read_table(sample_weight_2, index_col=0)
sample_weight_3_df = pd.read_table(sample_weight_3, index_col=0)

sample_weight_2_df = sample_weight_2_df.loc[:, sample_weight_2_df.columns.str.contains('vae')]
sample_weight_3_df = sample_weight_3_df.loc[:, sample_weight_3_df.columns.str.contains('vae')]

sample_weight_2_df.head()


# In[7]:


combined_model = sample_weight_2_df.merge(sample_weight_3_df,
                                          left_index=True,
                                          right_index=True,
                                          suffixes=('_two', '_three'))
combined_model.index = combined_model.index.astype('str')
combined_model.head(2)


# In[8]:


# Load hetnets for the given metaedge for rapid latent feature interpretation
hetnets = load_hetnets(
    hetnet_file='../3.build-hetnets/hetnets/interpret_hetnet.json.bz2',
    permuted_directory='../3.build-hetnets/hetnets/permuted/',
    subset_genes=combined_model.index,
    metaedge_abbrev=metaedge
   )


# ## Apply Interpret Compression Approach

# In[9]:


get_ipython().run_cell_magic('time', '', "\nmult_results = {}\nall_list = []\nfor model in hetnets.keys():\n    hetnet = hetnets[model]\n    mult_results[model] = combined_model.T @ hetnet\n    long_result = mult_results[model].reset_index().melt(id_vars=['index'])\n    long_result = long_result.assign(model=model)\n    if model != 'real':\n        long_result = long_result.assign(model_type='permuted')\n    else:\n        long_result = long_result.assign(model_type='real')\n    all_list.append(long_result)")


# In[10]:


all_df = pd.concat(all_list)
all_df.value = all_df.value.astype(float)
all_df.head()


# In[11]:


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


# In[12]:


z_score = (real_df.reset_index(drop=True).value - permuted_mean.value) / permuted_std.value
real_df = real_df.reset_index(drop=True).assign(z_score=z_score)

real_df.head()


# In[13]:


file = os.path.join('results', 'interpret_compression_gtex_vae_example.tsv')
real_df.sort_values(by='z_score').to_csv(file, sep='\t')


# ## Apply GSEA to the same features

# In[14]:


vae_2_lm = latentModel(filename=sample_weight_2,
                       z_dim=2,
                       dataset_name='GTEX',
                       algorithm_name='VAE',
                       weight_seed='109635',
                       shuffled_true=False)

vae_3_lm = latentModel(filename=sample_weight_3,
                       z_dim=2,
                       dataset_name='GTEX',
                       algorithm_name='VAE',
                       weight_seed='174930',
                       shuffled_true=False)


# In[15]:


xcell_gmt = parse_gmt(gene_sets=['../3.build-hetnets/data/xcell_all_entrez.gmt'])
len(xcell_gmt)


# In[16]:


get_ipython().run_cell_magic('time', '', 'vae_feature_0_zdim_2 = run_gsea_prerank(gene_score_df=vae_2_lm.w_df.vae_0,\n                                        gene_sets=xcell_gmt,\n                                        permutation_num=1000)\n\nvae_feature_1_zdim_2 = run_gsea_prerank(gene_score_df=vae_2_lm.w_df.vae_1,\n                                        gene_sets=xcell_gmt,\n                                        permutation_num=1000)\n\nvae_feature_0_zdim_3 = run_gsea_prerank(gene_score_df=vae_3_lm.w_df.vae_0,\n                                        gene_sets=xcell_gmt,\n                                        permutation_num=1000)\n\nvae_feature_1_zdim_3 = run_gsea_prerank(gene_score_df=vae_3_lm.w_df.vae_1,\n                                        gene_sets=xcell_gmt,\n                                        permutation_num=1000)\n\nvae_feature_2_zdim_3 = run_gsea_prerank(gene_score_df=vae_3_lm.w_df.vae_2,\n                                        gene_sets=xcell_gmt,\n                                        permutation_num=1000)')


# In[17]:


vae_feature_0_zdim_2 = vae_feature_0_zdim_2.assign(z_dim=2, feature=0)
vae_feature_1_zdim_2 = vae_feature_1_zdim_2.assign(z_dim=2, feature=1)

vae_feature_0_zdim_3 = vae_feature_0_zdim_3.assign(z_dim=3, feature=0)
vae_feature_1_zdim_3 = vae_feature_1_zdim_3.assign(z_dim=3, feature=1)
vae_feature_2_zdim_3 = vae_feature_2_zdim_3.assign(z_dim=3, feature=2)


# In[18]:


gsea_results_df = pd.concat([
    vae_feature_0_zdim_2,
    vae_feature_1_zdim_2,
    vae_feature_0_zdim_3,
    vae_feature_1_zdim_3,
    vae_feature_2_zdim_3
])

gsea_results_df.head()


# In[19]:


file = os.path.join('results', 'gsea_gtex_vae_example.tsv')
gsea_results_df.sort_values(by='fdr').to_csv(file, sep='\t')

