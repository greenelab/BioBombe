
# coding: utf-8

# # Process xCell cell-type genesets
# 
# The genesets were extracted from additional file 3 in [Aran et al. 2017](https://doi.org/10.1186/s13059-017-1349-1 "xCell: digitally portraying the tissue cellular heterogeneity landscape"). The data consist of 489 gene signatures. Many of these signatures represent redundant cell-types with different evidence.

# In[1]:


import os
from datetime import datetime
import csv
import pandas as pd


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


# Load Additional File 3
xcell_file = os.path.join('data', '13059_2017_1349_MOESM3_ESM.xlsx')
xcell_df = pd.read_excel(xcell_file, index_col=0)
print(xcell_df.shape)
xcell_df.head()


# In[4]:


# Distribution of number of cell types
cell_types = pd.Series([x.split('_')[0] for x in xcell_df.index])
cell_types.value_counts().head(10)


# In[5]:


# Distribution of number of cell types
cell_types.value_counts().hist(bins=10)


# In[6]:


geneset_url = 'https://doi.org/10.1186/s13059-017-1349-1'


# In[7]:


# Load curated gene names from versioned resource 
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
    .drop_duplicates('entrez_gene_id')
)

symbol_to_entrez = dict(zip(gene_df.symbol,
                            gene_df.entrez_gene_id))

# Add alternative symbols to entrez mapping dictionary
gene_df = gene_df.dropna(axis='rows', subset=['synonyms'])
gene_df.synonyms = gene_df.synonyms.str.split('|')

all_syn = (
    gene_df.apply(lambda x: pd.Series(x.synonyms), axis=1)
    .stack()
    .reset_index(level=1, drop=True)
)

all_syn.name = 'all_synonyms'
gene_df = gene_df.join(all_syn)

symbol_to_entrez_alt = dict(zip(gene_df.all_synonyms,
                                gene_df.entrez_gene_id))

symbol_to_entrez.update(symbol_to_entrez_alt)


# In[8]:


len(symbol_to_entrez)


# In[9]:


# Load gene updater
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# In[10]:


xcell_gmt_file = os.path.join('data', 'xcell_all_entrez.gmt')
filtered_genes = []
with open(xcell_gmt_file, 'w') as csvfile:
    xcell_writer = csv.writer(csvfile, delimiter='\t')
    for row in xcell_df.iterrows():
        geneset_name = row[0].replace(' ', '-')
        genes = row[1][2:].dropna().tolist()

        # There are data issues of input datetime types
        # Ignore these genes for now
        # See https://doi.org/10.1186/s13059-016-1044-7
        genes = [x for x in genes if not type(x) == datetime.date ]

        # Convert genes to entrez IDs (consider alternate mappings)
        genes = [symbol_to_entrez[x] for x in genes if x in symbol_to_entrez ]

        # Update old to new entrez IDs
        genes = [old_to_new_entrez[x] if x in old_to_new_entrez else x for x in genes]

        # Remove empty cells in the cell file
        #genes = [x for x in genes if x != '']
        row_to_write = [geneset_name, geneset_url] + list(set(genes))
        xcell_writer.writerow(row_to_write)

