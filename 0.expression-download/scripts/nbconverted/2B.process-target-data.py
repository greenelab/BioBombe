#!/usr/bin/env python
# coding: utf-8

# # Process TARGET PanCancer Data
# 
# Retrieve the downloaded expression data, update gene identifiers to entrez, and curate sample IDs. The script will also identify a balanced hold-out test set to compare projection performance into learned latent spaces across algorithms.

# In[1]:


import os
import random
import pandas as pd
from sklearn.model_selection import train_test_split


# In[2]:


random.seed(1234)


# ## Read Phenotype Information

# In[3]:


path = os.path.join('download', 'TARGET_phenotype.gz')
pheno_df = pd.read_table(path)

print(pheno_df.shape)
pheno_df.head(3)


# ## Read Entrez ID Curation Information
# 
# Load curated gene names from versioned resource. See https://github.com/cognoma/genes for more details

# In[4]:


# Commit from https://github.com/cognoma/genes
genes_commit = 'ad9631bb4e77e2cdc5413b0d77cb8f7e93fc5bee'


# In[5]:


url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(genes_commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
)

print(gene_df.shape)
gene_df.head(2)


# In[6]:


# Load gene updater - old to new Entrez gene identifiers
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(genes_commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# ## Read Probe Mapping Info

# In[7]:


path = os.path.join('download', 'gencode.v23.annotation.gene.probemap')
probe_map_df = pd.read_table(path)

# Inner merge gene df to get ensembl to entrez mapping
probe_map_df = probe_map_df.merge(gene_df, how='inner', left_on='gene', right_on='symbol')

# Mapping to rename gene expression index
ensembl_to_entrez = dict(zip(probe_map_df.id, probe_map_df.entrez_gene_id))

print(probe_map_df.shape)
probe_map_df.head(3)


# ## Read Gene Expression Data

# In[8]:


file = os.path.join('download', 'target_RSEM_gene_fpkm.gz')
expr_df = pd.read_table(file, index_col=0)

print(expr_df.shape)
expr_df.head(2)


# ## Process gene expression matrix
# 
# This involves updating Entrez gene ids, sorting and subsetting

# In[9]:


expr_df = (expr_df
    .dropna(axis='rows')
    .reindex(probe_map_df.id)
    .rename(index=ensembl_to_entrez)
    .rename(index=old_to_new_entrez)
    .groupby(level=0).mean()
    .transpose()
    .sort_index(axis='rows')
    .sort_index(axis='columns')
)

expr_df.index.rename('sample_id', inplace=True)

print(expr_df.shape)
expr_df.head(2)


# ## Stratify Balanced Training and Testing Sets in TARGET Gene Expression
# 
# Output training and testing gene expression datasets

# In[10]:


strat = pheno_df.set_index('sample_id').reindex(expr_df.index).primary_disease_code


# In[11]:


cancertype_count_df = (
    pd.DataFrame(strat.value_counts())
    .reset_index()
    .rename({'index': 'cancertype', 'primary_disease_code': 'n ='}, axis='columns')
)

file = os.path.join('data', 'target_sample_counts.tsv')
cancertype_count_df.to_csv(file, sep='\t', index=False)

cancertype_count_df


# In[12]:


train_df, test_df = train_test_split(expr_df,
                                     test_size=0.1,
                                     random_state=123,
                                     stratify=strat)


# In[13]:


print(train_df.shape)
test_df.shape


# In[14]:


train_file = os.path.join('data', 'train_target_expression_matrix_processed.tsv.gz')
train_df.to_csv(train_file, sep='\t', compression='gzip', float_format='%.3g')


# In[15]:


test_file = os.path.join('data', 'test_target_expression_matrix_processed.tsv.gz')
test_df.to_csv(test_file, sep='\t', compression='gzip', float_format='%.3g')


# ## Sort genes based on median absolute deviation and output to file

# In[16]:


# Determine most variably expressed genes and subset
mad_genes_df = pd.DataFrame(train_df.mad(axis=0).sort_values(ascending=False)).reset_index()
mad_genes_df.columns = ['gene_id', 'median_absolute_deviation']

file = os.path.join('data', 'target_mad_genes.tsv')
mad_genes_df.to_csv(file, sep='\t', index=False)

