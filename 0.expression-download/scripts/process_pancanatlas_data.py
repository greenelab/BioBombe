
# coding: utf-8

# # Process TCGA PanCanAtlas Data
# 
# Retrieve the downloaded expression data, update gene identifiers to entrez, and curate sample IDs. The script will also identify a balanced hold-out test set to compare projection performance into learned latent spaces across algorithms.

# In[1]:


import os
import random
import pandas as pd
from sklearn.model_selection import train_test_split


# In[2]:


random.seed(1234)


# ## Read TCGA Barcode Curation Information
# 
# Extract information from TCGA barcodes - `cancer-type` and `sample-type`. See https://github.com/cognoma/cancer-data for more details

# In[3]:


# Commit from https://github.com/cognoma/cancer-data/
sample_commit = 'da832c5edc1ca4d3f665b038d15b19fced724f4c'


# In[4]:


url = 'https://raw.githubusercontent.com/cognoma/cancer-data/{}/mapping/tcga_cancertype_codes.csv'.format(sample_commit)
cancer_types_df = pd.read_csv(url, dtype='str', keep_default_na=False)
cancertype_codes_dict = dict(zip(cancer_types_df['TSS Code'], cancer_types_df.acronym))
cancer_types_df.head(2)


# In[5]:


url = 'https://raw.githubusercontent.com/cognoma/cancer-data/{}/mapping/tcga_sampletype_codes.csv'.format(sample_commit)
sample_types_df = pd.read_csv(url, dtype='str')
sampletype_codes_dict = dict(zip(sample_types_df.Code, sample_types_df.Definition))
sample_types_df.head(2)


# ## Read Entrez ID Curation Information
# 
# Load curated gene names from versioned resource. See https://github.com/cognoma/genes for more details

# In[6]:


# Commit from https://github.com/cognoma/genes
genes_commit = 'ad9631bb4e77e2cdc5413b0d77cb8f7e93fc5bee'


# In[7]:


url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(genes_commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
)

print(gene_df.shape)
gene_df.head(2)


# In[8]:


# Load gene updater - old to new Entrez gene identifiers
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(genes_commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# ## Read Gene Expression Data

# In[9]:


file = os.path.join('download', 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv')
tcga_expr_df = pd.read_table(file, index_col=0)
tcga_expr_df.head()


# ## Process gene expression matrix
# 
# This involves updating Entrez gene ids, sorting and subsetting

# In[10]:


# Set index as entrez_gene_id
tcga_expr_df.index = tcga_expr_df.index.map(lambda x: x.split('|')[1])


# In[11]:


tcga_expr_df = (tcga_expr_df
    .dropna(axis='rows')
    .rename(index=old_to_new_entrez)
    .groupby(level=0).mean()
    .transpose()
    .sort_index(axis='rows')
    .sort_index(axis='columns')
)

tcga_expr_df.index.rename('sample_id', inplace=True)


# In[12]:


# Update sample IDs
tcga_expr_df.index = tcga_expr_df.index.str.slice(start=0, stop=15)
tcga_expr_df = tcga_expr_df.loc[~tcga_expr_df.index.duplicated(), :]


# In[13]:


# Filter for valid Entrez gene identifiers
tcga_expr_df = tcga_expr_df.loc[:, tcga_expr_df.columns.isin(gene_df.entrez_gene_id.astype(str))]


# In[14]:


print(tcga_expr_df.shape)
tcga_expr_df.head()


# ## Process TCGA cancer-type and sample-type info from barcodes
# 
# Cancer-type includes `OV`, `BRCA`, `LUSC`, `LUAD`, etc. while sample-type includes `Primary`, `Metastatic`, `Solid Tissue Normal`, etc.
# 
# See https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes for more details.
# 
# The goal is to use this info to stratify training (90%) and testing (10%) balanced by cancer-type and sample-type. 

# In[15]:


# Extract sample time in the order of the gene expression matrix
tcga_id = pd.DataFrame(tcga_expr_df.index)

# Extract the last two digits of the barcode and recode sample-type
tcga_id = tcga_id.assign(sample_type = tcga_id.sample_id.str[-2:])
tcga_id.sample_type = tcga_id.sample_type.replace(sampletype_codes_dict)

# Extract the first two ID numbers after `TCGA-` and recode cancer-type
tcga_id = tcga_id.assign(cancer_type = tcga_id.sample_id.str[5:7])
tcga_id.cancer_type = tcga_id.cancer_type.replace(cancertype_codes_dict)

# Append cancer-type with sample-type to generate stratification variable
tcga_id = tcga_id.assign(stratify_samples = tcga_id.cancer_type.str.cat(tcga_id.sample_type))

# Get stratification counts - function cannot work with singleton strats
stratify_counts = tcga_id.stratify_samples.value_counts().to_dict()

# Recode stratification variables if they are singletons
tcga_id = tcga_id.assign(stratify_samples_count = tcga_id.stratify_samples)
tcga_id.stratify_samples_count = tcga_id.stratify_samples_count.replace(stratify_counts)
tcga_id.loc[tcga_id.stratify_samples_count == 1, "stratify_samples"] = "other"


# In[16]:


tcga_id.head()


# ## Stratify Balanced Training and Testing Sets in TCGA Gene Expression
# 
# Output training and testing gene expression datasets

# In[17]:


train_df, test_df = train_test_split(tcga_expr_df,
                                     test_size=0.1,
                                     random_state=123,
                                     stratify=tcga_id.stratify_samples_count)


# In[18]:


print(train_df.shape)
test_df.shape


# In[19]:


train_file = os.path.join('data', 'train_tcga_expression_matrix_processed.tsv.gz')
train_df.to_csv(train_file, sep='\t', compression='gzip', float_format='%.3g')


# In[20]:


train_file = os.path.join('data', 'test_tcga_expression_matrix_processed.tsv.gz')
test_df.to_csv(train_file, sep='\t', compression='gzip', float_format='%.3g')

