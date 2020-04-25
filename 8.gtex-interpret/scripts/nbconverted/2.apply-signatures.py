#!/usr/bin/env python
# coding: utf-8

# # Application of Compressed Gene Expression Signatures
# 
# **Gregory Way, 2018**
# 
# We previously identified a specific VAE feature (z = 3) that captured blood signatures.
# These signatures were not captured in VAE z = 2, and, when present, contributed to a rapid increase in the ability to capture the signal in GTEX blood tissues.
# 
# The primary differences between the two VAE models appeared to be related to neutrophil and monocyte signatures. Here, we test the ability of these signatures to generalize to external datasets.
# 
# 
# ## Part 1:
# 
# ### Enrichment of Neutrophil Signatures
# 
# Here, we apply the VAE feature enriched for neutrophil genes to a publicly available dataset capturing neutrophil differentiation in two leukemia cell lines.
# 
# ![cell_type_Neutrophils_HPCA_2.png](https://github.com/greenelab/BioBombe/raw/master/6.biobombe-projection/figures/GTEX/signal/GpXCELL/gene_set_Neutrophils_HPCA_2.png)
# 
# ## Part 2:
# 
# ### Enrichment of Monocyte Signatures
# 
# Here, we apply the VAE features enriched for monocyte genes to a different publicly available dataset that captures various cell-types undergoing hematopoiesis.
# 
# ![cell_type_Monocytes_FANTOM_2.png](https://github.com/greenelab/BioBombe/raw/master/6.biobombe-projection/figures/GTEX/signal/GpXCELL/gene_set_Monocytes_FANTOM_2.png)
# 
# ## Output
# 
# In both cases, various scores are output that will be visualized in a separate notebook.

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scripts.utils import load_weight_matrix, apply_signature, load_enrichment_results

sys.path.append('../scripts')
from latent import parse_gmt


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


np.random.seed(123)


# In[4]:


# First, load the XCELL dataset and extract genes
geneset_file = os.path.join('..', '3.build-hetnets', 'data', 'xcell_all_entrez.gmt')

xcell_genesets_gmt = parse_gmt(gene_sets=[geneset_file])
len(xcell_genesets_gmt)


# ## 1.0. Load and Process External Neutrophil Dataset

# In[5]:


file = os.path.join('data', 'GSE103706_processed_matrix.tsv.gz')
geo_scaled_zeroone_df = pd.read_table(file, index_col=0)

print(geo_scaled_zeroone_df.shape)
geo_scaled_zeroone_df.head(2)


# ## 1.1. Apply Signature from VAE z = 3 (feature 0)
# 
# We are using feature 0 from VAE z = 3 because we previously observed an enrichment of a Neutrophil signature in this specific feature.

# In[6]:


vae_z3_seed = 908341
vae_z3_feature = 'vae_0'


# In[7]:


weight_z3_df = load_weight_matrix(dataset='GTEX',
                                  z_dim=3,
                                  seed=vae_z3_seed)

result_vae_3_feat0, neutrophil_missing_genes = (
    apply_signature(weight_df=weight_z3_df,
                    other_df=geo_scaled_zeroone_df,
                    feature=vae_z3_feature,
                    align=True)
)

use_genes = weight_z3_df.shape[0] - len(neutrophil_missing_genes)
print('{} ({}%) of genes are used'.format(use_genes, use_genes / weight_z3_df.shape[0] * 100 ))


# In[8]:


neutrophil_hpca_genes = xcell_genesets_gmt['Neutrophils_HPCA_2']
neutrophil_hpca_genes


# In[9]:


# But how many of the missing genes belong to the specific Neutrophil signature
neutrophil_hpca_genes_missing = [x for x in neutrophil_hpca_genes if int(x) in neutrophil_missing_genes]
neutrophil_hpca_genes_missing


# ## 1.2. Apply Signature from VAE z = 14 (feature 10)
# 
# We are using VAE z = 14 because, for some reason, this dimension was best able to capture the `Neutrophil_HPCA_2` xCell geneset.
# This geneset was the same geneset that was enriched in feature 0 for VAE z = 3.
# 
# We identified feature 10 as the feature with the greatest enrichment by visually scanning the biobombe results in `results/gtex/gpxcell/signal/gtex_z_14_GpXCELL__geneset_scores.tsv`

# In[10]:


z_14_df = load_enrichment_results(dataset="GTEX",
                                  metaedge="GpXCELL",
                                  z_dim=14)

z_14_df = (
    z_14_df
    .query("variable == 'Neutrophils_HPCA_2'")
    .assign(abs_z_score=z_14_df.z_score.abs())
    .sort_values(by='abs_z_score', ascending=False)
)

z_14_df.head(1)


# In[11]:


vae_z14_seed = z_14_df.head(1).seed.values[0]
vae_z14_feature = '{}_{}'.format(z_14_df.head(1).algorithm.values[0],
                                 z_14_df.head(1).feature.values[0])

print(vae_z14_seed)
vae_z14_feature


# In[12]:


weight_z14_df = load_weight_matrix(dataset='GTEX',
                                   z_dim=14,
                                   seed=vae_z14_seed)

result_vae_14_feat10, _ = apply_signature(weight_df=weight_z14_df,
                                          other_df=geo_scaled_zeroone_df,
                                          feature=vae_z14_feature,
                                          align=True)


# ## 1.3. Combine Data and Add Phenotype Information

# In[13]:


full_neutrophil_results_df = result_vae_14_feat10.merge(result_vae_3_feat0,
                                                        left_index=True,
                                                        right_index=True)


# In[14]:


# Process phenotype data
cell_line = [x[0] for x in result_vae_3_feat0.index.str.split(',')]
treatment = [x[1] for x in result_vae_3_feat0.index.str.split(',')]
day = [x[2].strip(' ') if 'replicate' not in x[2] else 'day 0'
       for x in result_vae_3_feat0.index.str.split(',')]


# In[15]:


full_neutrophil_results_df = (
    full_neutrophil_results_df
    .assign(cell_line=cell_line,
            treatment=treatment,
            day=day)
    .reset_index()
    .rename(columns={'index': 'full_id'})
)

recode_labels = {' not differentiated': 'Not Differentiated',
                 ' DMSO': 'DMSO',
                 ' DMSO+Nutridoma': 'DMSO+Nutridoma'}

full_neutrophil_results_df.treatment = full_neutrophil_results_df.treatment.replace(recode_labels)

file = os.path.join('results', 'neutrophil_data_biobombe_results.tsv')
full_neutrophil_results_df.to_csv(file, index=False, sep='\t')

full_neutrophil_results_df


# ## 1.4. Quickly Visualize Signature Applied to External Dataset
# 
# Note the final figure is compiled in an alternative notebook

# In[16]:


# Quickly visualize results
plt.rcParams['figure.figsize'] = 5, 3
ax = sns.stripplot(x='cell_line',
                   y='vae_0',
                   data=full_neutrophil_results_df,
                   hue='treatment',
                   jitter=0.1)

handles, labels = ax.get_legend_handles_labels()

l = plt.legend(handles,
               labels,
               bbox_to_anchor=(1.02, 0.8),
               loc=2,
               borderaxespad=0.)
l.set_title("Treatment")

ax.set_ylabel('VAE Feature 0 (z = 3)')
ax.set_xlabel('Cell Lines')

plt.tight_layout()


# In[17]:


# Quickly visualize results
plt.rcParams['figure.figsize'] = 5, 3
ax = sns.stripplot(x='cell_line',
                   y='vae_10',
                   data=full_neutrophil_results_df,
                   hue='treatment',
                   jitter=0.1)

handles, labels = ax.get_legend_handles_labels()

l = plt.legend(handles,
               labels,
               bbox_to_anchor=(1.02, 0.8),
               loc=2,
               borderaxespad=0.)
l.set_title("Treatment")

ax.set_ylabel('VAE Feature 10 (z = 14)')
ax.set_xlabel('Cell Lines')

plt.tight_layout()


# ## 1.5. Compare the two signatures derived from both models
# 
# Use the function `load_enrichment_results` to retrieve and subset previously compiled BioBombe results.
# The files are located in `6.biobombe-projection/results/`.

# In[18]:


# What other genesets are enriched in VAE z = 3 feature?
vae_z3_p0 = load_enrichment_results(dataset='GTEX',
                                    z_dim=3,
                                    metaedge='GpXCELL',
                                    algorithm='vae',
                                    feature=0,
                                    seed=vae_z3_seed,
                                    shuffled=False)

# What other genesets are enriched in VAE z = 14 feature?
vae_z14_p10 = load_enrichment_results(dataset='GTEX',
                                      z_dim=14,
                                      metaedge='GpXCELL',
                                      algorithm='vae',
                                      feature=10,
                                      seed=vae_z14_seed,
                                      shuffled=False)

full_test_df = vae_z3_p0.merge(vae_z14_p10, on='variable', suffixes=('_z3', '_z14'))

file = os.path.join('results', 'latent_feature_enrichment_comparison_neutrophil_genesets.tsv')
full_test_df.to_csv(file, index=False, sep='\t')

full_test_df.head(3)


# In[19]:


# Save dataframe for better plotting in R, but visualize quickly here
sns.scatterplot(data=full_test_df, x='z_score_z3', y='z_score_z14');


# ### 1.5.1 Determine Gene Weights across the Two Features
# 
# Also assign labels to which genesets the genes contribute to

# In[20]:


unique_genes = []
for geneset_name, geneset in xcell_genesets_gmt.items():
    for gene in geneset:
        unique_genes.append(gene)

unique_genes = set(unique_genes)

classification_genes = []
for geneset_name, geneset in xcell_genesets_gmt.items():
    if 'neutrophil' in geneset_name.lower():
        classification = 'Neutrophils'
    elif 'keratinocytes' in geneset_name.lower():
        classification = 'Keratinocytes'
    elif 'neurons' in geneset_name.lower():
        classification = 'Neurons'
    elif 'skeletal' in geneset_name.lower():
        classification = 'Skeletal Muscle'
    elif 'monocytes' in geneset_name.lower():
        classification = 'Monocytes'
    else:
        classification = 'Other Geneset'
    for gene in geneset:
        if gene in weight_z14_df.index:
            classification_genes.append([classification, gene, geneset_name])


# In[21]:


result_df = (
    pd.DataFrame(classification_genes, columns=['classification', 'gene', 'gene_set'])
    .sort_values(by='classification')
    .reset_index(drop=True)
    .drop_duplicates(subset='gene', keep='first')
)

result_df.gene = result_df.gene.astype(str)
result_df.index = result_df.gene
result_df.head()


# In[22]:


result_df.classification.value_counts()


# In[23]:


both_weight_df = (
    weight_z3_df.merge(weight_z14_df,
                       left_index=True,
                       right_index=True,
                       suffixes=('_3', '_14'))
    .merge(result_df, left_index=True,
           right_index=True, how='left')
    .fillna('No Geneset')
)

file = os.path.join('results', 'latent_feature_enrichment_comparison_neutrophil_genes.tsv')
both_weight_df.to_csv(file, index=False, sep='\t')


# ## 2.0. Load and Process External Hematopoietic Dataset

# In[24]:


file = os.path.join('data', 'GSE24759_processed_matrix.tsv.gz')
heme_zeroone_df = pd.read_table(file, index_col=0)

print(heme_zeroone_df.shape)
heme_zeroone_df.head(2)


# In[25]:


heme_z3_seed = 908341
heme_z3_feature = 'vae_2'


# In[26]:


# Transform the external dataset with this learned feature
weight_heme_z3_df = load_weight_matrix(dataset='GTEX',
                                       z_dim=3,
                                       seed=heme_z3_seed)

result_heme_vae_3_feat2, monocyte_missing_genes = (
    apply_signature(weight_df=weight_heme_z3_df,
                    other_df=heme_zeroone_df,
                    feature=heme_z3_feature,
                    align=True)
)

use_genes = weight_heme_z3_df.shape[0] - len(monocyte_missing_genes)
print('{} ({:.2f}%) of genes are used'.format(use_genes, use_genes / weight_heme_z3_df.shape[0] * 100 ))


# In[27]:


monocyte_fantom_genes = xcell_genesets_gmt['Monocytes_FANTOM_2']
print(len(monocyte_fantom_genes))
monocyte_fantom_genes


# In[28]:


# But how many of the missing genes belong to the specific Monocyte signature
monocyte_fantom_genes_missing = [x for x in monocyte_fantom_genes if int(x) in monocyte_missing_genes]
print(len(monocyte_fantom_genes_missing))
monocyte_fantom_genes_missing


# In[29]:


print('{:.2f}% of monocyte genes are missing'.format(len(monocyte_fantom_genes_missing) / len(monocyte_fantom_genes) * 100 ))


# In[30]:


# Additionall, the top scoring feature for Monocytes_FANTOM_2 is in the nmf model with 200 features
file = os.path.join('..', '6.biobombe-projection', 'results', 'gtex',
                    'gpxcell', 'signal',
                    'gtex_z_200_GpXCELL__geneset_scores.tsv.gz')

gtex_z200_scores_df = (
    pd.read_table(file)
    .query('variable == "Monocytes_FANTOM_2"'))

gtex_z200_scores_df = (
    gtex_z200_scores_df
    .assign(abs_z_score = gtex_z200_scores_df.z_score.abs())
    .sort_values(by='abs_z_score', ascending=False)
    .head(1)
)

gtex_z200_scores_df


# In[31]:


heme_z200_feature = '{}_{}'.format(gtex_z200_scores_df.algorithm.values[0],
                                   gtex_z200_scores_df.feature.values[0])

heme_z200_feature


# In[32]:


# Obtain this transformation too
weight_heme_z200_df = load_weight_matrix(dataset='GTEX',
                                         z_dim=200,
                                         seed=gtex_z200_scores_df.seed.values[0])


result_heme_nmf_200_feat6, _ = apply_signature(weight_df=weight_heme_z200_df,
                                            other_df=heme_zeroone_df,
                                            feature=heme_z200_feature,
                                            align=True)


# In[33]:


# Combine the full scores and output for downstream visualization
full_heme_result_df = (
    result_heme_vae_3_feat2
    .merge(result_heme_nmf_200_feat6, left_index=True, right_index=True)
    .reset_index().rename(columns={'index': 'cell'})
)


# In[34]:


heme_cell_type_recode_df = (
    pd.DataFrame(full_heme_result_df.cell.str.split('_').values.tolist(),
                 columns = ['cell_type', 'replicate', 'additional'])
)

heme_cell_type_recode_df.loc[~heme_cell_type_recode_df.additional.isna(), 'cell_type'] = "PRE_BCELL2"


# In[35]:


full_heme_result_df = (
    pd.concat([heme_cell_type_recode_df.drop(['additional'], axis='columns'),
               full_heme_result_df], axis='columns')
)


# In[36]:


# Recode cell-type into larger classification
file = os.path.join('results', 'cell-type-classification.tsv')
cell_class_df = pd.read_table(file)

cell_updater = dict(zip(cell_class_df.label, cell_class_df.classification))

cell_class_df.head()


# In[37]:


full_heme_result_df = (
    full_heme_result_df
    .assign(cell_class = full_heme_result_df.cell_type.replace(cell_updater))
)

file = os.path.join('results', 'hematopoietic_data_biobombe_results.tsv')
full_heme_result_df.to_csv(file, index=False, sep='\t')

full_heme_result_df.head()


# In[38]:


# Quickly plot results for both features
plt.rcParams['figure.figsize'] = 5, 3
ax = sns.stripplot(y='vae_2', x = 'cell_class', data = full_heme_result_df, hue = 'cell_class', jitter=0.1)

handles, labels = ax.get_legend_handles_labels()

l = plt.legend(handles,
               labels,
               bbox_to_anchor=(1.02, 0.8), loc=2, borderaxespad=0.)


# In[39]:


plt.rcParams['figure.figsize'] = 5, 3
ax = sns.stripplot(y='nmf_6', x = 'cell_class', data = full_heme_result_df, hue = 'cell_class', jitter=0.1)

handles, labels = ax.get_legend_handles_labels()

l = plt.legend(handles,
               labels,
               bbox_to_anchor=(1.02, 0.8), loc=2, borderaxespad=0.)


# ## Generate Supplementary Table of Neutrophil and Monocyte Signatures

# In[40]:


# Load curated gene names from versioned resource 
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
)

entrez_to_symbol = dict(zip(gene_df.entrez_gene_id,
                            gene_df.symbol))


# In[41]:


neutrophil_df = pd.DataFrame(neutrophil_hpca_genes, columns=['entrez_gene_id'])
neutrophil_df = (
    neutrophil_df.assign(
        gene_symbol=neutrophil_df.entrez_gene_id.astype(int).replace(entrez_to_symbol),
        signature='Neutrophil_HPCA_2',
        in_external_dataset='Yes')
)
neutrophil_df.head()


# In[42]:


monocyte_df = pd.DataFrame(monocyte_fantom_genes, columns=['entrez_gene_id'])
monocyte_df = (
    monocyte_df.assign(
        gene_symbol=monocyte_df.entrez_gene_id.astype(int).replace(entrez_to_symbol),
        signature='Monocyte_FANTOM_2',
        in_external_dataset='Yes')
)

monocyte_df.loc[monocyte_df.entrez_gene_id.isin(monocyte_fantom_genes_missing), 'in_external_dataset'] = "No"
monocyte_df.head()


# In[43]:


sup_table_df = pd.concat([neutrophil_df, monocyte_df], axis='rows').reset_index(drop=True)

file = os.path.join("results", "neutrophil_and_monocyte_signature_genes.tsv")
sup_table_df.to_csv(file, sep='\t', index=False)

