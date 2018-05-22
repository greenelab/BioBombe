
# coding: utf-8

# # Download TARGET Pan Cancer Data for Compression
# 
# The notebook downloads gene expression and clinical data from The TARGET project. The data is downloaded from [UCSC Xena](https://xenabrowser.net/datapages/?dataset=target_RSEM_gene_fpkm&host=https%3A%2F%2Ftoil.xenahubs.net).
# 
# The data is in `log2(FPKM)` RSEM transformed.

# In[1]:


import os
from urllib.request import urlretrieve


# In[2]:


# Get Gene Expression Data
url = 'https://toil.xenahubs.net/download/'
name = 'target_RSEM_gene_fpkm.gz'

path = os.path.join('download', name)


# In[3]:


urlretrieve('{}{}'.format(url, name), path)


# In[4]:


get_ipython().system(" sha256sum 'download/target_RSEM_gene_fpkm.gz'")


# In[5]:


# Get Probe Mappings
name = 'gencode.v23.annotation.gene.probeMap.gz'

path = os.path.join('download', name)


# In[6]:


urlretrieve('{}{}'.format(url, name), path)


# In[7]:


get_ipython().system(" sha256sum 'download/gencode.v23.annotation.gene.probeMap.gz'")


# In[8]:


# Get Sample Identifiers
name = 'TARGET_phenotype.gz'

path = os.path.join('download', name)


# In[9]:


urlretrieve('{}{}'.format(url, name), path)


# In[10]:


get_ipython().system(" sha256sum 'download/TARGET_phenotype.gz'")

