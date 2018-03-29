
# coding: utf-8

# # Download MSigDB Gene Sets
# 
# **Gregory Way 2018**
# 
# **Modified scripts originally written by Daniel Himmelstein (@dhimmel)**
# 
# _Most_ MSigDB gene sets (_version 6.1_) are now `CC BY 4.0` (except KEGG, BioCarta, and AAAS/STKE). Download and process.

# In[1]:


import csv
import pandas as pd


# In[2]:


# MSigDB version
version = '6.1'


# In[3]:


# Download full MSigDB matrix
# NOTE - This fill is not added to the repository because it contains
# gene sets with restrictive licenses
url_prefix = 'https://www.broadinstitute.org/gsea/resources/msigdb/'
url = '{}{}/msigdb.v{}.symbols.gmt'.format(url_prefix, version, version)
get_ipython().system(" wget --timestamping --no-verbose --directory-prefix 'data' $url")


# In[4]:


# Many of the genesets have sub gene sets - process these as well
msigdb_dict = {
    'c1.all': 'positional gene sets',
    'c2.cgp': 'chemical and genetic perturbations',
    'c2.cp.reactome': 'Reactome gene sets',
    'c3.mir': 'microRNA targets',
    'c3.tft': 'transcription factor targets',
    'c4.cgn': 'cancer gene neighborhoods',
    'c4.cm': 'cancer modules',
    'c5.bp': 'GO biological processes',
    'c5.cc': 'GO cellular components',
    'c5.mf': 'GO molecular functions',
    'c6.all': 'oncogenic signatures',
    'c7.all': 'immunologic signatrues'
}

for gene_set in msigdb_dict:
    url = '{}{}/{}.v{}.symbols.gmt'.format(url_prefix, version, gene_set, version)
    get_ipython().system(" wget --timestamping --no-verbose --directory-prefix 'data' $url")

