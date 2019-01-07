
# coding: utf-8

# # Download Publicly Available Neutrophil Dataset
# 
# **Gregory Way, 2018**
# 
# Here, I download [GSE103706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103706) which is associated with [Rincon et al. 2018](https://doi.org/10.1186/s12864-018-4957-6).
# 
# This dataset includes two acute myeloid leukemia (AML) cell lines; PLB-985 and HL-60.
# There are 14 samples total in this dataset.
# The cell lines are exposed to two treatments - DMSO and DMSO+Nutridoma.
# The treatments are demonstrated to induce neutrophil differentiation in these cell lines.
# 
# We hypothesized that our constructed feature identified through our interpret compression approach would have higher activation patterns in the cell lines with induced neutrophil differentiation.

# In[1]:


import os
import csv
import pandas as pd
from sklearn import preprocessing

from scripts.utils import download_geo


# In[2]:


base_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103706/suppl/'
name = 'GSE103706_Merged_data_FPKM_and_normalized.xlsx'
directory = 'download'


# In[3]:


download_geo(base_url, name, directory)


# In[4]:


path = 'download/GSE103706_Merged_data_FPKM_and_normalized.xlsx'
get_ipython().system(' sha256sum $path')


# ## Process the Data

# In[5]:


# Load Data
geo_df = pd.read_excel(path, index_col=0, skiprows=1)

print(geo_df.shape)
geo_df.head(2)


# ## Update Gene Names

# In[6]:


# Load curated gene names from versioned resource 
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
)

symbol_to_entrez = dict(zip(gene_df.symbol,
                            gene_df.entrez_gene_id))


# In[7]:


# Add alternative symbols to entrez mapping dictionary
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


# In[8]:


# Create a synonym to entrez mapping and add to dictionary
synonym_to_entrez = dict(zip(gene_with_syn_df.all_synonyms,
                             gene_with_syn_df.entrez_gene_id))

symbol_to_entrez.update(synonym_to_entrez)


# In[9]:


# Load gene updater
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# In[10]:


# Update the symbol column to entrez_gene_id
geo_df.symbol = geo_df.symbol.replace(symbol_to_entrez)
geo_df = geo_df.loc[geo_df.symbol.isin(symbol_to_entrez.values()), :]
geo_df.symbol = geo_df.symbol.replace(old_to_new_entrez)
geo_df.index = geo_df.symbol
geo_df.index.name = 'entrez_gene_id'
geo_df = geo_df.drop(['ens_gene_id', 'ncbi_gene_id', 'gene_short', 'symbol'], axis='columns')


# ## Scale Data and Output to File

# In[11]:


# Scale RNAseq data using zero-one normalization
geo_scaled_zeroone_df = preprocessing.MinMaxScaler().fit_transform(geo_df.transpose())
geo_scaled_zeroone_df = pd.DataFrame(geo_scaled_zeroone_df,
                                     columns=geo_df.index,
                                     index=geo_df.columns)

os.makedirs('data', exist_ok=True)

file = os.path.join('data', 'GSE103706_processed_matrix.tsv.gz')
geo_scaled_zeroone_df.to_csv(file, sep='\t', compression='gzip')

geo_scaled_zeroone_df.head()

