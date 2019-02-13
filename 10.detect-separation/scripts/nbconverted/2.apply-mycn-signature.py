
# coding: utf-8

# # Applying MYCN Signature Learned in TARGET data to Cell Lines
# 
# **Gregory Way, 2019**

# In[1]:


import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

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


vae_seed = '451283'
vae_k = 200
vae_feature = "vae_111"


# In[6]:


weight_df = load_weight_matrix(dataset='TARGET',
                               z_dim=vae_k,
                               seed=vae_seed)
weight_df.head()


# In[7]:


result_mycn_df, nbl_missing_genes = (
    apply_signature(weight_df=weight_df,
                    other_df=nbl_df,
                    feature=vae_feature,
                    align=True)
)

use_genes = weight_df.shape[0] - len(nbl_missing_genes)
print('{} ({}%) of genes are used'.format(use_genes, use_genes / weight_df.shape[0] * 100 ))


# In[8]:


result_mycn_df.head()


# ## 3. Align with Phenotype Data

# In[9]:


file = os.path.join("download", "nbl_cellline_phenotype.txt")
pheno_df = pd.read_table(file)

pheno_df.head()


# In[10]:


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


# In[11]:


merged_df['MYCN status'].value_counts()


# In[12]:


sns.boxplot(x="MYCN status", y='vae_111', data=merged_df);

