
# coding: utf-8

# # Download Publicly Available Hematopoietic Dataset
# 
# **Gregory Way, 2018**
# 
# Here, I download [GSE24759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24759) which is associated with [Novershtern et al. 2011](https://doi.org/10.1016/j.cell.2011.01.004).
# 
# This dataset includes 211 samples consisting of 38 distinct hematopoietic states in various stages of differentiation.
# 
# We hypothesized that our constructed feature identified through our interpret compression approach would have higher activation patterns in Monocytes.

# In[1]:


import os
import csv
import pandas as pd
from sklearn import preprocessing

from scripts.utils import download_geo


# In[2]:


base_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24759/suppl/'
name = 'GSE24759_data.sort.txt.gz'
directory = 'download'


# In[3]:


download_geo(base_url, name, directory)


# In[4]:


path = 'download/GSE24759_data.sort.txt.gz'
get_ipython().system(' sha256sum $path')


# ## Process the Data

# In[5]:


# Load Additional File 3
geo_df = pd.read_table(path)

print(geo_df.shape)
geo_df.head(2)


# ## Update Gene Names

# In[6]:


# Load gene updater
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# In[7]:


# Update the entrez gene IDs in the index
#entrez_ids = geo_
geo_df.index = geo_df.A_Name.replace(old_to_new_entrez)
geo_df.index.name = 'entrez_gene_id'
geo_df = geo_df.drop(['A_Name', 'A_Desc'], axis='columns')
geo_df.head(2)


# ## Scale Data and Output to File

# In[8]:


# Scale RNAseq data using zero-one normalization
geo_scaled_zeroone_df = preprocessing.MinMaxScaler().fit_transform(geo_df.transpose())
geo_scaled_zeroone_df = pd.DataFrame(geo_scaled_zeroone_df,
                                     columns=geo_df.index,
                                     index=geo_df.columns)

os.makedirs('data', exist_ok=True)

file = os.path.join('data', 'GSE24759_processed_matrix.tsv.gz')
geo_scaled_zeroone_df.to_csv(file, sep='\t', compression='gzip')

geo_scaled_zeroone_df.head()


# ## Process Cell-Type Classification
# 
# Data acquired from Supplementary Table 1 of [Novershtern et al. 2011](https://doi.org/10.1016/j.cell.2011.01.004)

# In[9]:


cell_class = {
    # Hematopoietic Stem Cells
    'HSC1': 'HSC',
    'HSC2': 'HSC',
    'HSC3': 'HSC',
    
    # Myeloid Progenitors
    'CMP': 'Myeloid',
    'MEP': 'Myeloid',
    'GMP': 'Myeloid',
    
    # Erythroid Populations
    'ERY1': 'Erythroid',
    'ERY2': 'Erythroid',
    'ERY3': 'Erythroid',
    'ERY4': 'Erythroid',
    'ERY5': 'Erythroid',
    
    # Megakaryocytic Populations
    'MEGA1': 'Megakaryocytic',
    'MEGA2': 'Megakaryocytic',
    
    # Granulocytic Populations
    'GRAN1': 'Granulocytic',
    'GRAN2': 'Granulocytic',
    'GRAN3': 'Granulocytic',
    
    # Monocyte Population
    'MONO1': 'Monocyte',
    'MONO2': 'Monocyte',
    
    # Basophil Population
    'BASO1': 'Basophil',
    
    # Eosinophil Population
    'EOS2': 'Eosinophil',
    
    # B Lymphoid Progenitors
    'PRE_BCELL2': 'B Lymphoid Progenitor',
    'PRE_BCELL3': 'B Lymphoid Progenitor',
    
    # Naive Lymphoid Progenitors
    'BCELLA1': 'Naive Lymphoid',
    'TCELLA6': 'Naive Lymphoid',
    'TCELLA2': 'Naive Lymphoid',
    
    # Differentiated B Cells
    'BCELLA2': 'Differentiated B Cell',
    'BCELLA3': 'Differentiated B Cell',
    'BCELLA4': 'Differentiated B Cell',
    
    # Differentiated T Cells
    'TCELLA7': 'Differentiated T Cell',
    'TCELLA8': 'Differentiated T Cell',
    'TCELLA1': 'Differentiated T Cell',
    'TCELLA3': 'Differentiated T Cell',
    'TCELLA4': 'Differentiated T Cell',
    
    # Natural Killer Population
    'NKA1': 'NK Cell',
    'NKA2': 'NK Cell',
    'NKA3': 'NK Cell',
    'NKA4': 'NK Cell',
    
    # Dendritic Cell
    'DENDA1': 'Dendritic',
    'DENDA2': 'Dendritic',
}


# In[10]:


cell_class_df = (
    pd.DataFrame(cell_class, index=[0])
    .transpose()
    .reset_index()
    .rename(columns={'index': 'label', 0: 'classification'})
)

cell_class_df.head()


# In[11]:


file = os.path.join('results', 'cell-type-classification.tsv')
cell_class_df.to_csv(file, sep='\t', index=False)

