
# coding: utf-8

# # Processing MSigDB Gene Sets into Binary Matrix
# 
# This notebook loads the full MSigDB gene set `.gmt` file (version 6.1) and outputs a binary, gene by gene set matrix indicating gene membership in the given gene set.
# 
# **Note that we exclude gene sets with restrictive licences (KEGG, Biocarta, and The AAAS/STKE Cell Signaling Database)**

# In[1]:


import os
import csv
import numpy as np
import pandas as pd


# In[2]:


def make_template_matrix(msigdb_file, blacklist, checkblacklist=True):
    """
    Retrieve all genes and pathways from given msigdb .gmt file
    
    Output:
    sorted gene by pathways pandas dataframe. Entries indicate membership
    """
    all_db_pathways = []
    all_db_genes = []

    # Get a set of all genes and all pathways in MSigDB (not blacklisted)
    with open(msigdb_file, 'r') as msigdb_fh:
        msigdb_reader = csv.reader(msigdb_fh, delimiter='\t')

        for row in msigdb_reader:
            signature_name = row[0]
            signature_genes = row[2:]
            
            if checkblacklist:
                if signature_name.startswith(blacklist):
                    continue

            all_db_pathways.append(signature_name)
            all_db_genes += signature_genes
        
    big_msigdb_df = pd.DataFrame(0, index=set(all_db_genes), columns=all_db_pathways)
    big_msigdb_df = big_msigdb_df.sort_index()
    big_msigdb_df = big_msigdb_df.T.sort_index().T
    
    # Loop through file again to populate dataframe. This is a fast implementation
    with open(msigdb_file, 'r') as msigdb_fh:
        msigdb_reader = csv.reader(msigdb_fh, delimiter='\t')
        for row in msigdb_reader:
            signature_name = row[0]
            signature_genes = row[2:]
            if checkblacklist:
                if signature_name.startswith(blacklist):
                    continue

            for gene in signature_genes:
                big_msigdb_df.at[gene, signature_name] = 1

    return big_msigdb_df


# In[3]:


# Store .gmt files
full_msigdb_file = os.path.join('data', 'msigdb.v6.1.symbols.gmt')

# Resources with restrictive licenses
blacklist = ('KEGG', 'BIOCARTA', 'ST_')


# ## Process MSigDB gmt files into large matrix

# In[4]:


get_ipython().run_cell_magic('time', '', 'full_msigdb_df = make_template_matrix(full_msigdb_file, blacklist, checkblacklist=True)\nprint(full_msigdb_df.shape)')


# In[5]:


get_ipython().run_cell_magic('time', '', "full_msigdb_file = os.path.join('data', 'full_msigdb_binary_matrix.tsv.bz2')\nfull_msigdb_df.to_csv(full_msigdb_file, sep='\\t', compression='bz2')")

