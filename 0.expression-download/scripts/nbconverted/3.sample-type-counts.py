
# coding: utf-8

# # Describe tissue types of input data
# 
# **Gregory Way, 2019**
# 
# Load in previously identified tissue type counts and output a supplementary table.

# In[1]:


import os
import pandas as pd


# In[2]:


# Load TCGA data
file = os.path.join('data', 'tcga_sample_counts.tsv')
tcga_count_df = pd.read_table(file, sep='\t').rename({'cancertype': 'tissue'}, axis='columns')
tcga_count_df = tcga_count_df.assign(dataset="TCGA")
tcga_count_df.head()


# In[3]:


# Load GTEX data
file = os.path.join('data', 'gtex_sample_counts.tsv')
gtex_count_df = pd.read_table(file, sep='\t').rename({'tissuetype': 'tissue'}, axis='columns')
gtex_count_df = gtex_count_df.assign(dataset="GTEX")
gtex_count_df.head()


# In[4]:


# Load TARGET data
file = os.path.join('data', 'target_sample_counts.tsv')
target_count_df = pd.read_table(file, sep='\t').rename({'cancertype': 'tissue'}, axis='columns')
target_count_df = target_count_df.assign(dataset="TARGET")
target_count_df.head()


# In[5]:


# Combine all data to generate supplementary table
full_count_df = (
    pd.concat([tcga_count_df, gtex_count_df, target_count_df], axis='rows')
    .sort_values(by='tissue', ascending=True)
    .reset_index(drop=True)
)


file = os.path.join('results', 'full_sample_counts.tsv')
full_count_df.to_csv(file, sep='\t', index=False)

full_count_df

