
# coding: utf-8

# # Capturing Blood Tissue in GTEx Data
# 
# **Gregory Way, 2018**
# 
# When viewing sample correlation differences across z dimensions stratified by tissue-type, we often observed a rapid increase in correlation after increasing model capacity by one.
# 
# For example, the ability of variational autoencoders to capture blood tissue correlation jumps by nearly 50% between bottleneck dimensions 2 and 3 (see below).
# 
# ![sample-correlation_Blood_GTEX_signal_pearson.png](https://raw.githubusercontent.com/greenelab/BioBombe/master/4.analyze-components/figures/GTEX/sample-correlation/sample-type/sample-correlation_Blood_GTEX_signal_pearson.png)
# 
# ## Procedure
# 
# In the following notebook, we extract two representative weight matrices for VAE latent space dimensions 2 and 3.
# We apply two compression feature interpretaion approaches to the weight vectors for each latent space feature. 
# 
# 1. Our matrix interpretation approach
# 2. Overrepresentation Tests using high weight genes
# 
# In both approaches we use genesets derived in the xCell paper that represent cell-types ([Aran et al. 2017](https://doi.org/10.1186/s13059-017-1349-1))
# 
# We output the results of both approaches and analyze the results in subsequent notebooks.

# In[1]:


import os
import sys
import glob
import pandas as pd
import numpy as np

sys.path.append('../scripts')
from latent import latentModel, load_hetnets, parse_gmt, run_overrepresentation


# In[2]:


np.random.seed(123)


# In[3]:


metaedge = 'GpXCELL'


# In[4]:


base_dir = os.path.join('..', '2.ensemble-z-analysis', 'results',
                        'GTEX_results', 'ensemble_z_matrices')


# In[5]:


sample_weight_2 = os.path.join(base_dir, 'gtex_components_2', 'model_908341_weight_matrix.tsv.gz')
sample_weight_3 = os.path.join(base_dir, 'gtex_components_3', 'model_908341_weight_matrix.tsv.gz')


# In[6]:


sample_weight_2_df = pd.read_table(sample_weight_2, index_col=0)
sample_weight_3_df = pd.read_table(sample_weight_3, index_col=0)

sample_weight_2_df = sample_weight_2_df.loc[:, sample_weight_2_df.columns.str.contains('vae')]
sample_weight_3_df = sample_weight_3_df.loc[:, sample_weight_3_df.columns.str.contains('vae')]

# Recode column names
sample_weight_2_df.columns = sample_weight_2_df.columns + '_two'
sample_weight_3_df.columns = sample_weight_3_df.columns + '_three'

sample_weight_2_df.head()


# In[7]:


combined_model = sample_weight_2_df.merge(sample_weight_3_df,
                                          left_index=True,
                                          right_index=True)
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
all_df = all_df.rename(columns={'index': 'full_feature', 'value': 'raw_score'})
all_df.head()


# In[11]:


all_group = all_df.groupby(['model_type', 'full_feature', 'variable'])
all_group_mean = all_group.mean().reset_index()

real_df = all_group_mean.query("model_type == 'real'")
permuted_mean = all_group_mean.query("model_type == 'permuted'")

permuted_std = (
    all_group
    .std()
    .reset_index()
    .query("model_type == 'permuted'")
)


# In[12]:


z_score = (real_df.reset_index(drop=True).raw_score - permuted_mean.raw_score) / permuted_std.raw_score
real_df = real_df.reset_index(drop=True).assign(z_score=z_score)

real_df.head()


# In[13]:


# Output z scores per xCell genesets for all VAE features (5 total)
file = os.path.join('results', 'gtex_vae_example_interpret_compression.tsv')
real_df.sort_values(by='z_score').to_csv(file, sep='\t', index=False)


# In[14]:


# Determine the most distinguishing features between the two models. This will test
# which genesets are the most differently enriched _in sum_ between z = 2 and z = 3
feature_info_df = (
    pd.DataFrame(real_df['full_feature'].str.split('_').values.tolist(),
                 columns=['algorithm', 'feature', 'model_z'])

)


# In[15]:


feature_info_df = pd.concat([real_df, feature_info_df], axis=1)
feature_info_df = feature_info_df.assign(abs_z_score = feature_info_df.z_score.abs())


# In[16]:


feature_info_df = (
    feature_info_df.groupby(['variable', 'model_z'])
    .mean()
    .reset_index()
    .pivot(index='variable', columns='model_z', values = 'abs_z_score')
)

feature_info_df = (
    feature_info_df.assign(abs_diff = (feature_info_df.three - feature_info_df.two).abs())
    .sort_values(by='abs_diff', ascending=False)
)

feature_info_df.head()


# In[17]:


# Output the most distinguishing features considering all features within each z dimension
file = os.path.join('results', 'gtex_vae_example_differentiating_features.tsv')
feature_info_df.to_csv(file, sep='\t')


# ## Perform overrepresentation tests on the same features
# 
# ### Split into positive and negative tails and extract high weight genes

# In[18]:


geneset_file = os.path.join('..', '3.build-hetnets', 'data', 'xcell_all_entrez.gmt')
xcell_genesets_gmt = parse_gmt(gene_sets=[geneset_file])
len(xcell_genesets_gmt)


# In[19]:


vae_2_lm = latentModel(filename=sample_weight_2,
                       z_dim=2,
                       dataset_name='GTEX',
                       algorithm_name='VAE',
                       weight_seed='908341',
                       shuffled_true=False)

vae_3_lm = latentModel(filename=sample_weight_3,
                       z_dim=3,
                       dataset_name='GTEX',
                       algorithm_name='VAE',
                       weight_seed='908341',
                       shuffled_true=False)


# In[20]:


vae_2_lm.get_high_weight_genes()
vae_3_lm.get_high_weight_genes()


# In[21]:


background_genes = []
for xcell_name, xcell_gene_set in xcell_genesets_gmt.items():
    background_genes += xcell_gene_set


# In[22]:


get_ipython().run_cell_magic('time', '', '\ngene_list = vae_2_lm.w_df[vae_2_lm.pos_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_2_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_2_lm.w_df[vae_2_lm.neg_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_2_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_2_lm.w_df[vae_2_lm.pos_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_2_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_2_lm.w_df[vae_2_lm.neg_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_2_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.pos_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_3_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.neg_high_w_df.vae_0].index.tolist()\nvae_feature_0_zdim_3_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.pos_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_3_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.neg_high_w_df.vae_1].index.tolist()\nvae_feature_1_zdim_3_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.pos_high_w_df.vae_2].index.tolist()\nvae_feature_2_zdim_3_pos = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)\n\ngene_list = vae_3_lm.w_df[vae_3_lm.neg_high_w_df.vae_2].index.tolist()\nvae_feature_2_zdim_3_neg = run_overrepresentation(gene_list=gene_list,\n                                                  gene_set_dict=xcell_genesets_gmt,\n                                                  background_genes=background_genes)')


# In[23]:


vae_feature_0_zdim_2_pos = vae_feature_0_zdim_2_pos.assign(feature='vae_0_two', tailed='pos')
vae_feature_0_zdim_2_neg = vae_feature_0_zdim_2_neg.assign(feature='vae_0_two', tailed='neg')
vae_feature_1_zdim_2_pos = vae_feature_1_zdim_2_pos.assign(feature='vae_1_two', tailed='pos')
vae_feature_1_zdim_2_neg = vae_feature_1_zdim_2_neg.assign(feature='vae_1_two', tailed='neg')

vae_feature_0_zdim_3_pos = vae_feature_0_zdim_3_pos.assign(feature='vae_0_three', tailed='pos')
vae_feature_0_zdim_3_neg = vae_feature_0_zdim_3_neg.assign(feature='vae_0_three', tailed='neg')
vae_feature_1_zdim_3_pos = vae_feature_1_zdim_3_pos.assign(feature='vae_1_three', tailed='pos')
vae_feature_1_zdim_3_neg = vae_feature_1_zdim_3_neg.assign(feature='vae_1_three', tailed='neg')
vae_feature_2_zdim_3_pos = vae_feature_2_zdim_3_pos.assign(feature='vae_2_three', tailed='pos')
vae_feature_2_zdim_3_neg = vae_feature_2_zdim_3_neg.assign(feature='vae_2_three', tailed='neg')


# In[24]:


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


# In[25]:


file = os.path.join('results', 'gtex_vae_example_overrepresentation.tsv')
(
    overrepresented_results_df
    .reset_index()
    .rename(columns={'index': 'variable',
                     'feature': 'full_feature'})
    .sort_values(by='pval')
    .to_csv(file, sep='\t', index=False)
)

