#!/usr/bin/env python
# coding: utf-8

# # Download GTEx Phenotype Data
# 
# **Gregory Way, 2019**

# In[1]:


import os
from urllib.request import urlretrieve


# In[2]:


name = "GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
url = "https://storage.googleapis.com/gtex_analysis_v7/annotations/{}".format(name)
path = os.path.join('download', name)


# In[3]:


urlretrieve(url, path)


# In[4]:


get_ipython().system(' md5sum ../0.expression-download/download/GTEx_v7_Annotations_SubjectPhenotypesDS.txt')


# In[5]:


name = "GTEx_v7_Annotations_SampleAttributesDS.txt"
url = "https://storage.googleapis.com/gtex_analysis_v7/annotations/{}".format(name)
path = os.path.join('download', name)


# In[6]:


urlretrieve(url, path)


# In[7]:


get_ipython().system(' md5sum ../0.expression-download/download/GTEx_v7_Annotations_SampleAttributesDS.txt')

