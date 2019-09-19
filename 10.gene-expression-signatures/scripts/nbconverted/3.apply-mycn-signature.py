#!/usr/bin/env python
# coding: utf-8

# # Applying MYCN Signature Learned in TARGET data to Cell Lines
# 
# **Gregory Way, 2019**

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("../8.gtex-interpret")
from scripts.utils import load_weight_matrix, apply_signature


# In[2]:


np.random.seed(123)


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')


# ## 1. Load NBL Processed Dataset

# In[4]:


file = os.path.join('data', 'nbl_celllines_processed_matrix.tsv.gz')
nbl_df = pd.read_table(file, index_col=0)

print(nbl_df.shape)
nbl_df.head(2)


# ## 2. Apply Signature Learned by TARGET data

# In[5]:


# Identify the top performing feature
file = os.path.join("results", "nbl_mycn_separation_target_t_test.tsv")
target_full_results_df = pd.read_table(file).head(5)
target_full_results_df


# In[6]:


top_seed = str(target_full_results_df.seed.values[0])
top_k = int(target_full_results_df.z_dim.values[0])
top_feature = "{}_{}".format(target_full_results_df.algorithm.values[0],
                             target_full_results_df.feature_num.values[0])

print(top_seed, top_k, top_feature)


# In[7]:


weight_df = load_weight_matrix(dataset='TARGET',
                               z_dim=top_k,
                               seed=top_seed)
weight_df.head()


# In[8]:


result_mycn_df, nbl_missing_genes = (
    apply_signature(weight_df=weight_df,
                    other_df=nbl_df,
                    feature=top_feature,
                    align=True)
)

use_genes = weight_df.shape[0] - len(nbl_missing_genes)
print('{} ({}%) of genes are used'.format(use_genes, use_genes / weight_df.shape[0] * 100 ))


# In[9]:


result_mycn_df.head()


# ## 3. Align with Phenotype Data

# In[10]:


file = os.path.join("download", "nbl_cellline_phenotype.txt")
pheno_df = pd.read_table(file)

pheno_df.head()


# In[11]:


merged_df = (
    result_mycn_df
    .merge(pheno_df,
           left_index=True,
           right_on="Cell Line")
    .reset_index(drop=True)
)

file = os.path.join("results", "mycn_nbl_scores.tsv")
merged_df.to_csv(file, sep='\t')

merged_df.head()


# In[12]:


merged_df['MYCN status'].value_counts()


# In[13]:


sns.boxplot(x="MYCN status", y=top_feature, data=merged_df);


# In[14]:


# Perform t-test on the result
amplified_scores = merged_df.loc[merged_df['MYCN status'] == "Amplified", top_feature]
notamplified_scores = merged_df.loc[merged_df['MYCN status'] != "Amplified", top_feature]

ttest_ind(amplified_scores, notamplified_scores, equal_var=False)

