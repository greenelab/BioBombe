
# coding: utf-8

# # Download and Process Neuroblastoma RNAseq Data
# 
# **Gregory Way 2019**
# 
# We are downloading the dataset associated with [Harenza et al. 2017](https://doi.org/10.1038/sdata.2017.33). The data profiles RNAseq data from 39 commonly used neuroblastoma (NBL) cell lines.
# 
# We are interested in the MYCN amplification status of these cell lines. We will test if the MYCN amplification score learned through the BioBombe signature approach applied to TARGET data generalizes to this cell line dataset.
# 
# MYCN Amplification refers to the number of copies of the _MYCN_ gene. MYCN amplification is used as a biomarker for poor prognosis in neuroblastoma patients ([Huang and Weiss 2013](https://doi.org/10.1101/cshperspect.a014415)).

# In[1]:


import os
import requests
import pandas as pd
from urllib.request import urlretrieve

from sklearn import preprocessing


# In[2]:


url = "https://ndownloader.figshare.com/files/14138792"
name = "2019-01-22-CellLineSTAR-fpkm-2pass_matrix.txt"
path = os.path.join("download", name)


# In[3]:


os.makedirs("download", exist_ok=True)


# In[4]:


urlretrieve(url, path)


# In[5]:


get_ipython().system(' md5sum "download/2019-01-22-CellLineSTAR-fpkm-2pass_matrix.txt"')


# ## Download Phenotype Data

# In[6]:


url = "https://www.nature.com/articles/sdata201733/tables/3"
name = "nbl_cellline_phenotype.txt"
path = os.path.join("download", name)


# In[7]:


html = requests.get(url).content

pheno_df = pd.read_html(html)[0]
pheno_df['Cell Line'] = pheno_df['Cell Line'].str.replace("-", "")

pheno_df.to_csv(path, sep='\t', index=False)

pheno_df.head()


# In[8]:


get_ipython().system(' md5sum "download/nbl_cellline_phenotype.txt"')


# ## Process RNAseq Data

# In[9]:


raw_file = os.path.join("download", "2019-01-22-CellLineSTAR-fpkm-2pass_matrix.txt")

raw_df = pd.read_table(raw_file, sep='\t')
raw_df.head()


# ### Update Gene Names

# In[10]:


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


# In[11]:


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


# In[12]:


# Create a synonym to entrez mapping and add to dictionary
synonym_to_entrez = dict(zip(gene_with_syn_df.all_synonyms,
                             gene_with_syn_df.entrez_gene_id))

symbol_to_entrez.update(synonym_to_entrez)


# In[13]:


# Load gene updater
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# In[14]:


gene_map = raw_df.GeneID.replace(symbol_to_entrez)
gene_map = gene_map.replace(old_to_new_entrez)


# In[15]:


raw_df.index = gene_map
raw_df.index.name = 'entrez_gene_id'
raw_df = raw_df.drop(['GeneID'], axis='columns')
raw_df = raw_df.loc[raw_df.index.isin(symbol_to_entrez.values()), :]

print(raw_df.shape)
raw_df.head()


# ## Scale Data and Output

# In[16]:


raw_scaled_df = preprocessing.MinMaxScaler().fit_transform(raw_df.transpose())
raw_scaled_df = (
    pd.DataFrame(raw_scaled_df,
                 columns=raw_df.index,
                 index=raw_df.columns)
    .sort_index(axis='columns')
    .sort_index(axis='rows')
)
raw_scaled_df.columns = raw_scaled_df.columns.astype(str)
raw_scaled_df = raw_scaled_df.loc[:, ~raw_scaled_df.columns.duplicated(keep='first')]

raw_scaled_df.head()


# In[17]:


os.makedirs('data', exist_ok=True)

file = os.path.join('data', 'nbl_celllines_processed_matrix.tsv.gz')
raw_scaled_df.to_csv(file, sep='\t', compression='gzip')

