
# coding: utf-8

# # Capturing Blood Tissue in GTEx Data
# 
# **Gregory Way, 2018**
# 
# When viewing sample correlation differences across z dimensions stratified by tissue-type, we often observed a rapid increase in correlation after increasing model capacity by one.
# 
# For example, the ability of variational autoencoders to capture blood tissue correlation jumps by nearly 50% between bottleneck dimensions 2 and 3 (see below).
# 
# ![sample-correlation_Blood_GTEX_signal_pearson.png](https://raw.githubusercontent.com/greenelab/interpret-compression/master/4.analyze-components/figures/GTEX/sample-correlation/sample-type/sample-correlation_Blood_GTEX_signal_pearson.png)
# 
# ## Procedure
# 
# In the following notebook, we extract two representative weight matrices for VAE latent space dimensions 2 and 3.
# We apply two compression feature interpretaion approaches to the weight vectors for each latent space feature. 
# 
# 1. Our matrix interpretation approach
# 2. Overrepresentation Tests using high weight genes
# 
# In both approaches we use genesets derived in the XCELL paper that represent cell-types ([Aran et al. 2017](https://doi.org/10.1186/s13059-017-1349-1))
# 
# We output the results of both approaches and analyze the results in subsequent notebooks.

# In[1]:


import os
import sys
import glob
import requests
import pandas as pd
import numpy as np
import gseapy as gp

sys.path.append('../scripts')
from latent import latentModel, load_hetnets, parse_gmt, run_overrepresentation


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


file = os.path.join('results', 'gtex_vae_example_interpret_compression.tsv')
real_df.sort_values(by='z_score').to_csv(file, sep='\t')


# ## Perform overrepresentation tests on the same features
# 
# ### Split into positive and negative tails and extract high weight genes

# In[14]:


geneset_file = os.path.join('..', '3.build-hetnets', 'data', 'xcell_all_entrez.gmt')
xcell_genesets_gmt = parse_gmt(gene_sets=[geneset_file])
len(xcell_genesets_gmt)


# In[15]:


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


# In[16]:


vae_2_lm.get_high_weight_genes()
vae_3_lm.get_high_weight_genes()


# In[17]:


background_genes = []
for xcell_name, xcell_gene_set in xcell_genesets_gmt.items():
    background_genes += xcell_gene_set


# In[18]:


get_ipython().run_cell_magic('time', '', '\ngene_list = vae_2_lm.w_df[vae_2_lm.pos_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_2_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_2_lm.w_df[vae_2_lm.neg_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_2_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_2_lm.w_df[vae_2_lm.pos_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_2_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_2_lm.w_df[vae_2_lm.neg_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_2_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.pos_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_3_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.neg_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_3_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.pos_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_3_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.neg_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_3_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.pos_high_w_df.vae_2].index.tolist()\nvae_feature_2_zdim_3_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.neg_high_w_df.vae_2].index.tolist()\nvae_feature_2_zdim_3_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)')


# In[19]:


vae_feature_0_zdim_2_pos = vae_feature_0_zdim_2_pos.assign(feature='vae_0_two', tailed='pos')
vae_feature_0_zdim_2_neg = vae_feature_0_zdim_2_neg.assign(feature='vae_0_two', tailed='neg')
vae_feature_1_zdim_2_pos = vae_feature_1_zdim_2_pos.assign(feature='vae_1_two', tailed='pos')
vae_feature_1_zdim_2_neg = vae_feature_1_zdim_2_neg.assign(feature='vae_1_two', tailed='neg')

vae_feature_0_zdim_3_pos = vae_feature_0_zdim_3_pos.assign(feature='vae_0_three', tailed='pos')
vae_feature_0_zdim_3_neg = vae_feature_0_zdim_3_neg.assign(feature='vae_0_three', tailed='neg')
vae_feature_1_zdim_3_pos = vae_feature_1_zdim_3_pos.assign(feature='vae_1_three', tailed='pos')
vae_feature_1_zdim_3_neg = vae_feature_1_zdim_3_neg.assign(feature='vae_1_three', tailed='neg')
vae_feature_2_zdim_3_pos = vae_feature_2_zdim_3_pos.assign(feature='vae_2', tailed='pos')
vae_feature_2_zdim_3_neg = vae_feature_2_zdim_3_neg.assign(feature='vae_2', tailed='neg')


# In[20]:


overrepresented_results_df = pd.concat([
    vae_feature_0_zdim_2_pos,
    vae_feature_0_zdim_2_neg,
    vae_feature_1_zdim_2_pos,
    vae_feature_1_zdim_2_neg,
    vae_feature_0_zdim_3_pos,
    vae_feature_0_zdim_3_neg,
    vae_feature_1_zdim_3_pos,
    vae_feature_1_zdim_3_neg,
    vae_feature_2_zdim_3_pos,
    vae_feature_2_zdim_3_neg
])

overrepresented_results_df.head()


# In[21]:


file = os.path.join('results', 'gtex_vae_example_overrepresentation.tsv')
overrepresented_results_df.sort_values(by='pval').to_csv(file, sep='\t')

