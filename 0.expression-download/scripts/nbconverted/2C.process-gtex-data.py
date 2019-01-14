
# coding: utf-8

# # Process GTEx Gene Expression Data
# 
# Retrieve the downloaded expression data, update gene identifiers to entrez, and curate sample IDs. The script will also identify a balanced hold-out test set to compare projection performance into learned latent spaces across algorithms.
# 
# **Note:** GTEx version 7 was downloaded from https://www.gtexportal.org/home/datasets

# In[1]:


import os
import random
import pandas as pd
from sklearn.model_selection import train_test_split


# In[2]:


random.seed(1234)


# ## Read Phenotype Information

# In[3]:


path = os.path.join('download', 'GTEx_v7_Annotations_SampleAttributesDS.txt')
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


# ## Read Gene Expression Data

# In[7]:


file = os.path.join('download', 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz')
expr_df = pd.read_table(file, sep='\t', skiprows=2, index_col=1)

print(expr_df.shape)

expr_df.head(2)


# In[8]:


# Get GTEx gene mapping
expr_gene_ids = (
    expr_df
    .loc[:, ['Name']]
    .reset_index()
    .drop_duplicates(subset='Description')
)

# Inner merge gene df to get ensembl to entrez mapping
map_df = expr_gene_ids.merge(gene_df, how='inner', left_on='Description', right_on='symbol')

symbol_to_entrez = dict(zip(map_df.symbol, map_df.entrez_gene_id))


# ## Process gene expression matrix
# 
# This involves updating Entrez gene ids, sorting and subsetting

# In[9]:


expr_df = (expr_df
 .drop(['Name'], axis='columns')
 .dropna(axis='rows')
 .groupby(level=0).mean()
 .reindex(map_df.symbol)
 .rename(index=symbol_to_entrez)
 .rename(index=old_to_new_entrez)
 .transpose()
 .sort_index(axis='rows')
 .sort_index(axis='columns')
)

expr_df.index.rename('sample_id', inplace=True)

print(expr_df.shape)
expr_df.head(2)


# ## Stratify Balanced Training and Testing Sets in GTEx Gene Expression
# 
# Output training and testing gene expression datasets.

# In[10]:


strat = pheno_df.set_index('SAMPID').reindex(expr_df.index).SMTSD


# In[11]:


tissuetype_count_df = (
    pd.DataFrame(strat.value_counts())
    .reset_index()
    .rename({'index': 'tissuetype', 'SMTSD': 'n ='}, axis='columns')
)

file = os.path.join('data', 'gtex_sample_counts.tsv')
tissuetype_count_df.to_csv(file, sep='\t', index=False)

tissuetype_count_df


# In[12]:


train_df, test_df = train_test_split(expr_df,
                                     test_size=0.1,
                                     random_state=123,
                                     stratify=strat)


# In[13]:


print(train_df.shape)
test_df.shape


# In[14]:


train_file = os.path.join('data', 'train_gtex_expression_matrix_processed.tsv.gz')
train_df.to_csv(train_file, sep='\t', compression='gzip', float_format='%.3g')


# In[15]:


test_file = os.path.join('data', 'test_gtex_expression_matrix_processed.tsv.gz')
test_df.to_csv(test_file, sep='\t', compression='gzip', float_format='%.3g')


# ## Sort genes based on median absolute deviation and output to file

# In[16]:


# Determine most variably expressed genes and subset
mad_genes_df = pd.DataFrame(train_df.mad(axis=0).sort_values(ascending=False)).reset_index()
mad_genes_df.columns = ['gene_id', 'median_absolute_deviation']

file = os.path.join('data', 'gtex_mad_genes.tsv')
mad_genes_df.to_csv(file, sep='\t', index=False)

