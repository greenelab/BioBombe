
# coding: utf-8

# # Determining the Top 50 mutated genes in TCGA PanCanAtlas
# 
# **Gregory Way, 2018**
# 
# Using the MC3 mutation calling effort, determine the top 50 most mutated genes in the TCGA PanCanAtlas.
# Details describing the mutation calling approach can be reviewed in [Ellrott et al. 2018](https://doi.org/10.1016/j.cels.2018.03.002).
# 
# These genes will be input into an elastic net logistic regression model to predict mutation status.

# In[1]:


import os
import pandas as pd


# In[2]:


# Load data to build y matrices
base_url = "https://github.com/greenelab/pancancer/raw"
commit = "2a0683b68017fb226f4053e63415e4356191734f"

# Load data
file = "{}/{}/data/sample_freeze.tsv".format(base_url, commit)
sample_freeze_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/pancan_mutation_freeze.tsv.gz".format(base_url, commit)
mutation_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/vogelstein_cancergenes.tsv".format(base_url, commit)
gene_type_df = pd.read_table(file, index_col=0)


# In[3]:


# Process gene classification (oncogene or tumor suppressor)
gene_type_df = gene_type_df.loc[:, 'Classification*']
gene_type_df = pd.DataFrame(gene_type_df).reset_index()
gene_type_df.columns = ['symbol', 'classification']
gene_type_df.head()


# In[4]:


# Reindex the mutation data to the frozen samples used
mutation_df = mutation_df.reindex(sample_freeze_df.SAMPLE_BARCODE)


# In[5]:


# Identify the top 50 mutated genes
top_50_mutated_genes_df = mutation_df.sum().sort_values(ascending=False).head(50).reset_index()
top_50_mutated_genes_df.columns = ['gene', 'num_mutations']

top_50_mutated_genes_df = (
    top_50_mutated_genes_df
    .merge(gene_type_df, left_on='gene', right_on='symbol', how='left')
    .drop('symbol', axis='columns')
    .fillna('neither')
)

top_50_mutated_genes_df


# In[6]:


# Write to file
file = os.path.join('data', 'top50_mutated_genes.tsv')
top_50_mutated_genes_df.to_csv(file, sep='\t', index=False)

