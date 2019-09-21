#!/usr/bin/env python
# coding: utf-8

# # Investigating Sex Signature Features
# 
# **Gregory Way, 2019**

# In[1]:


import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from sklearn import preprocessing

sys.path.append("../8.gtex-interpret")
from scripts.utils import load_weight_matrix, apply_signature


# In[2]:


np.random.seed(123)


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')


# ## Load and Process Gene Dictionary

# In[4]:


# Load curated gene names from versioned resource 
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(commit)

gene_df = pd.read_table(url)

symbol_to_entrez = dict(zip(gene_df.symbol,
                            gene_df.entrez_gene_id))

gene_df = gene_df.dropna(axis='rows', subset=['synonyms'])
gene_df.synonyms = gene_df.synonyms.str.split('|')

all_syn = (
    gene_df.apply(lambda x: pd.Series(x.synonyms), axis=1)
    .stack()
    .reset_index(level=1, drop=True)
)

# Name the synonym series and join with rest of genes
all_syn.name = 'all_synonyms'
gene_with_syn_df = gene_df.join(all_syn)

# Remove rows that have redundant symbols in all_synonyms
gene_with_syn_df = (
    gene_with_syn_df
    
    # Drop synonyms that are duplicated - can't be sure of mapping
    .drop_duplicates(['all_synonyms'], keep=False)

    # Drop rows in which the symbol appears in the list of synonyms
    .query('symbol not in all_synonyms')
)

# Create a synonym to entrez mapping and add to dictionary
synonym_to_entrez = dict(zip(gene_with_syn_df.all_synonyms,
                             gene_with_syn_df.entrez_gene_id))

symbol_to_entrez.update(synonym_to_entrez)

# Load gene updater
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))

gene_df.entrez_gene_id = gene_df.entrez_gene_id.astype(str)

print(gene_df.shape)
gene_df.head()


# ## Load and Proces Sex Genes
# 
# Using the resource `Sex-Associated Gene Database` (SAGD) ([Shi et al. 2018](https://doi.org/10.1093/nar/gky1040))
# 
# Downloading from http://bioinfo.life.hust.edu.cn/SAGD#!/browse_gene
# 
# Selecting human species, all tissues, all stages. The downloaded file is included in this repo.

# In[5]:


sex_gene_file = os.path.join("download", "browse_gene_9606.csv")
sex_gene_df = pd.read_csv(sex_gene_file)

# Translate the symbol column to entrez_gene_id
sex_gene_map = sex_gene_df.Symbol.replace(symbol_to_entrez)
sex_gene_map = sex_gene_map.replace(old_to_new_entrez)
sex_gene_df = sex_gene_df.assign(entrez_gene_id=sex_gene_map)

# Median collapse duplicate gene IDs across SAGD groups 
sex_gene_df = (
    sex_gene_df
    .groupby(["Species", "Symbol", "entrez_gene_id"])
    .median()
    .sort_values(by="Padj")
    .reset_index()
)

sex_gene_df.entrez_gene_id = sex_gene_df.entrez_gene_id.astype(str)
sex_gene_df = sex_gene_df.assign(neg_log_p=-1 * np.log10(sex_gene_df.Padj + 1e-300))

print(sex_gene_df.shape)
sex_gene_df.head()


# In[6]:


sex_gene_df.neg_log_p.hist(bins=100)


# ## Load Sex Signatures

# In[7]:


gtex_seed = '451283'
gtex_k = 200
gtex_feature = "nmf_111"


# In[8]:


# Load the gtex weight matrix containing the best sex feature
gtex_weight_df = (
    load_weight_matrix(dataset='GTEX',
                       z_dim=gtex_k,
                       seed=gtex_seed)
    .reset_index()
)
gtex_weight_df.gene_id = gtex_weight_df.gene_id.astype(str)

gtex_weight_df.head()


# In[9]:


# Align the weight matrix to the Sex Gene Database
gtex_sex_feature = (
    gtex_weight_df
    .merge(gene_df,
           left_on="gene_id",
           right_on="entrez_gene_id",
           how="left")
    .loc[:, gene_df.columns.tolist() + [gtex_feature]]
    .assign(abs_value_feature = gtex_weight_df.loc[:, gtex_feature].abs().tolist())
)

gtex_sex_feature.entrez_gene_id = gtex_sex_feature.entrez_gene_id.astype(str)

gtex_sex_feature = (
    gtex_sex_feature
    .merge(sex_gene_df,
           left_on="entrez_gene_id",
           right_on="entrez_gene_id",
           how="left")
    .sort_values(by="abs_value_feature", ascending=False)
    .dropna(subset=["entrez_gene_id", "Padj"])
    .reset_index(drop=True)
)

print(gtex_sex_feature.shape)

# Show the top 10 genes
gtex_sex_feature.head(10)


# In[10]:


gtex_sex_feature.plot(kind="scatter", x=gtex_feature, y="neg_log_p")


# ## TCGA Sex Signature

# In[11]:


tcga_seed = '165158'
tcga_k = 200
tcga_feature = "ica_151"


# In[12]:


# Load the TCGA weight matrix containing the best sex feature
tcga_weight_df = (
    load_weight_matrix(dataset='TCGA',
                       z_dim=tcga_k,
                       seed=tcga_seed)
    .reset_index()
)

tcga_weight_df.gene_id = tcga_weight_df.gene_id.astype(str)
tcga_weight_df.head()


# In[13]:


# Align the weight matrix to the Sex Gene Database
tcga_sex_feature = (
    tcga_weight_df
    .merge(gene_df,
           left_on="gene_id",
           right_on="entrez_gene_id",
           how="left")
    .loc[:, gene_df.columns.tolist() + [tcga_feature]]
    .assign(abs_value_feature = tcga_weight_df.loc[:, tcga_feature].abs().tolist())
)

tcga_sex_feature.entrez_gene_id = tcga_sex_feature.entrez_gene_id.astype(str)

tcga_sex_feature = (
    tcga_sex_feature
    .merge(sex_gene_df,
           left_on="entrez_gene_id",
           right_on="entrez_gene_id",
           how="left")
    .sort_values(by="abs_value_feature", ascending=False)
    .dropna(subset=["entrez_gene_id", "Padj"])
    .reset_index(drop=True)
)

print(tcga_sex_feature.shape)
tcga_sex_feature.head(20)


# In[14]:


tcga_sex_feature.plot(kind="scatter", x=tcga_feature, y="neg_log_p")

