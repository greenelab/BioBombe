#!/usr/bin/env python
# coding: utf-8

# # Download TCGA PanCanAtlas Data for Compression
# 
# The notebook downloads gene expression and clinical data from The Cancer Genome Atlas PanCanAtlas project. The data is accessed directly from the [Genome Data Commons](https://gdc.cancer.gov/about-data/publications/pancanatlas).

# In[1]:


import os
from urllib.request import urlretrieve


# In[2]:


url = 'http://api.gdc.cancer.gov/data/9a4679c3-855d-4055-8be9-3577ce10f66e'
name = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv'
path = os.path.join('download', name)


# In[3]:


urlretrieve(url, path)


# In[4]:


get_ipython().system(" sha256sum 'download/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv'")

