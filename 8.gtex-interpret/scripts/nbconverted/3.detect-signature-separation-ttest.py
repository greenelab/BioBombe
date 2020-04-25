#!/usr/bin/env python
# coding: utf-8

# # Applying Neutrophil and Monocyte Signatures using Top Features
# 
# **Gregory Way, 2018**
# 
# Instead of focusing only on two compressed features, like we did in [notebook 3](https://github.com/greenelab/BioBombe/blob/master/8.gtex-interpret/3.apply-signatures.ipynb), we apply all of the signatures in the following plots to their external validation data.
# 
# We previously demonstrated that the enrichment of specific genesets (like Neutrophils and Monocytes) were different across algorithms and k dimensions.
# This implies that certain biological features are best captured by different dimensions and algorithms.
# Here, we test if our network projection scores were associated with strength of signature separation across different groups of samples.
# 
# ## Part 1:
# 
# ### Enrichment of Neutrophil Signatures
# 
# Publicly available dataset capturing neutrophil differentiation in two leukemia cell lines.
# 
# ![cell_type_Neutrophils_HPCA_2.png](https://github.com/greenelab/BioBombe/raw/master/6.biobombe-projection/figures/GTEX/signal/GpXCELL/gene_set_Neutrophils_HPCA_2.png)
# 
# ## Part 2:
# 
# ### Enrichment of Monocyte Signatures
# 
# Publicly available dataset that captures various cell-types, including monocytes, undergoing hematopoiesis.
# 
# ![cell_type_Monocytes_FANTOM_2.png](https://github.com/greenelab/BioBombe/raw/master/6.biobombe-projection/figures/GTEX/signal/GpXCELL/gene_set_Monocytes_FANTOM_2.png)
# 
# ## Output
# 
# We output scores for all compression feature scores for both validation data sets.
# There are a total of 5 * 28 = 140 scores per dataset

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, ttest_rel
import matplotlib.pyplot as plt
import seaborn as sns

from scripts.utils import (
    load_weight_matrix,
    apply_signature,
    load_enrichment_results,
    extract_feature,
)

sys.path.append('../scripts')
from latent import parse_gmt


# In[2]:


genesets = ['Neutrophils_HPCA_2', 'Monocytes_FANTOM_2']


# In[3]:


# Get top scores across algorithms and z dimensions
score_dir = os.path.join("..", "6.biobombe-projection", "results", "gtex",
                         "gpxcell", "signal")

full_list = []
for score_file in os.listdir(score_dir):
    score_file = os.path.join(score_dir, score_file)
    
    top_df = pd.read_table(score_file)
    top_df = (
        top_df
        .query("variable in @genesets")
        .assign(abs_z_score=top_df.z_score.abs())
        .sort_values(by='abs_z_score', ascending=False)
        .groupby(['z', 'algorithm', 'variable'])
        .first()
        .reset_index()
    )
    
    full_list.append(top_df)


# In[4]:


# Confirm extraction worked (note different colors)
full_top_df = pd.concat(full_list)

sns.relplot(x="z", y="abs_z_score", hue="algorithm", data=full_top_df,
            col='variable', kind='line')


# In[5]:


# Compile feature matrix for each geneset
neutrophil_list = []
monocyte_list = []
for idx, feature in full_top_df.iterrows():
    z_dim = feature.z
    seed = feature.seed
    algorithm = feature.algorithm
    feature_num = feature.feature
    geneset = feature.variable
    weight_df = load_weight_matrix(dataset='GTEX',
                                   z_dim=z_dim,
                                   seed=seed)
    feature_df = extract_feature(weight_df=weight_df,
                                 algorithm=algorithm,
                                 feature=feature_num)

    rename_feature = '{}_zdim_{}_seed_{}'.format(feature_df.name, z_dim, seed)
    feature_df = feature_df.rename(rename_feature)
    
    if geneset == "Neutrophils_HPCA_2":
        neutrophil_list.append(feature_df)
    else:
        monocyte_list.append(feature_df)


# In[6]:


neutrophil_df = pd.concat(neutrophil_list, axis='columns')
monocyte_df = pd.concat(monocyte_list, axis='columns')


# ## 1.0. Load External Neutrophil Dataset

# In[7]:


file = os.path.join('data', 'GSE103706_processed_matrix.tsv.gz')
geo_scaled_zeroone_df = pd.read_table(file, index_col=0)

print(geo_scaled_zeroone_df.shape)
geo_scaled_zeroone_df.head(2)


# ## 1.1. Apply Signature from All Top Features

# In[8]:


neutrophil_result_df, missing_genes = (
    apply_signature(weight_df=neutrophil_df,
                    other_df=geo_scaled_zeroone_df,
                    align=True)
)

top_compressed_features = neutrophil_result_df.columns.tolist()
len(missing_genes)


# ## 1.2. Combine Data and Add Phenotype Information

# In[9]:


# Process phenotype data
cell_line = [x[0] for x in neutrophil_result_df.index.str.split(',')]
treatment = [x[1] for x in neutrophil_result_df.index.str.split(',')]
day = [x[2].strip(' ') if 'replicate' not in x[2] else 'day 0'
       for x in neutrophil_result_df.index.str.split(',')]


# In[10]:


neutrophil_result_df = (
    neutrophil_result_df
    .assign(cell_line=cell_line,
            treatment=treatment,
            day=day)
    .reset_index()
    .rename(columns={'index': 'full_id'})
)

recode_labels = {' not differentiated': 'control',
                 ' DMSO': 'treated',
                 ' DMSO+Nutridoma': 'treated'}

neutrophil_result_df.treatment = neutrophil_result_df.treatment.replace(recode_labels)

neutrophil_result_df.head(2)


# ## 1.3. Perform t-test on Treatment vs. Control

# In[11]:


ttest_results = []
for compressed_feature in top_compressed_features:
    signature_df = neutrophil_result_df.loc[:, [compressed_feature, 'treatment']]
    treatment_values = signature_df.query("treatment == 'treated'").iloc[:, 0].values
    control_values = signature_df.query("treatment == 'control'").iloc[:, 0].values
    
    t_stat, t_p = ttest_ind(treatment_values, control_values)
    ttest_results.append(pd.Series([compressed_feature, t_stat, t_p]))


# In[12]:


t_results_df = pd.concat(ttest_results, axis='columns').transpose()
t_results_df.columns = ['feature', 't_stat', 't_p']
t_results_df = t_results_df.assign(neg_log_p = -np.log10(t_results_df.t_p.astype(np.float64)))
t_results_df.head()


# In[13]:


neutrophils_top_df = full_top_df.query("variable == 'Neutrophils_HPCA_2'")
neutrophils_top_df.head(2)


# ## 1.4. Compile and Output Plotting Data

# In[14]:


final_neutrophil_results_df = (
    pd.DataFrame(t_results_df
                 .feature
                 .str
                 .split('_')
                 .values
                 .tolist(),
                 columns=['algorithm',
                          'feature_num',
                          "drop_z",
                          "z_dim",
                          "drop_seed",
                          "seed"])
    .merge(t_results_df,
           left_index=True,
           right_index=True)
    .drop(['drop_z',
           'drop_seed',
           'feature'],
          axis='columns')
)

final_neutrophil_results_df.loc[:, 'z_dim'] = final_neutrophil_results_df.z_dim.astype(np.int64)
final_neutrophil_results_df.loc[:, 'feature_num'] = final_neutrophil_results_df.feature_num.astype(np.int64)
final_neutrophil_results_df.loc[:, 'seed'] = final_neutrophil_results_df.seed.astype(np.int64)

final_neutrophil_results_df = (
    final_neutrophil_results_df
    .merge(neutrophils_top_df,
           left_on=['algorithm', 'feature_num', 'z_dim', 'seed'],
           right_on=['algorithm', 'feature', 'z', 'seed'])
)

final_neutrophil_results_df.head()


# In[15]:


sns.scatterplot(x="z_score", y="t_stat", data=final_neutrophil_results_df)


# In[16]:


file = os.path.join("results", "all_neutrophil_top_scores_and_separation.tsv")
final_neutrophil_results_df.to_csv(file, sep='\t', index=False)


# ## 2.0. Load External Monocyte Dataset
# 
# We perform a similar procedure, but apply top monocyte signatures to an alternative publicly available dataset.

# In[17]:


file = os.path.join('data', 'GSE24759_processed_matrix.tsv.gz')
heme_zeroone_df = pd.read_table(file, index_col=0)

print(heme_zeroone_df.shape)
heme_zeroone_df.head(2)


# ## 2.1. Apply Signature from All Top Features and Process Data

# In[18]:


result_df, missing_genes = apply_signature(weight_df=monocyte_df,
                                           other_df=heme_zeroone_df,
                                           align=True)
full_heme_result_df = result_df.reset_index().rename(columns={'index': 'cell'})

heme_cell_type_recode_df = (
    pd.DataFrame(full_heme_result_df.cell.str.split('_').values.tolist(),
                 columns = ['cell_type', 'replicate', 'additional'])
)

heme_cell_type_recode_df.loc[~heme_cell_type_recode_df.additional.isna(), 'cell_type'] = "PRE_BCELL2"

full_heme_result_df = (
    pd.concat([heme_cell_type_recode_df.drop(['additional'], axis='columns'),
               full_heme_result_df], axis='columns')
)

top_compressed_features = result_df.columns.tolist()
len(missing_genes)


# ## 2.2. Add Phenotype Information

# In[19]:


# Recode cell-type into larger classification
file = os.path.join('results', 'cell-type-classification.tsv')
cell_class_df = pd.read_table(file)

cell_updater = dict(zip(cell_class_df.label, cell_class_df.classification))
monocyte_updater = dict(zip(cell_class_df.label, cell_class_df.monocyte))

cell_class_df.head()


# In[20]:


full_heme_result_df = (
    full_heme_result_df
    .assign(cell_class = full_heme_result_df.cell_type.replace(cell_updater),
            monocyte_status = full_heme_result_df.cell_type.replace(monocyte_updater))
)

full_heme_result_df.head(2)


# ## 2.3. Perform t-test on Monocyte vs. Non-Monocyte

# In[21]:


ttest_results = []
for compressed_feature in top_compressed_features:
    signature_df = full_heme_result_df.loc[:, [compressed_feature, 'monocyte_status']]
    treatment_values = signature_df.query("monocyte_status == 'Monocyte'").iloc[:, 0].values
    control_values = signature_df.query("monocyte_status == 'Non Monocyte'").iloc[:, 0].values

    t_stat, t_p = ttest_ind(treatment_values, control_values)
    ttest_results.append(pd.Series([compressed_feature, t_stat, t_p]))


# In[22]:


monocyte_top_df = full_top_df.query("variable == 'Monocytes_FANTOM_2'")
monocyte_top_df.head(2)


# In[23]:


t_results_df = pd.concat(ttest_results, axis='columns').transpose()
t_results_df.columns = ['feature', 't_stat', 't_p']
t_results_df = t_results_df.assign(neg_log_p = -np.log10(t_results_df.t_p.astype(np.float64)))
t_results_df.head()


# ## 2.4. Compile and Output Plotting Data

# In[24]:


full_heme_result_df = (
    pd.DataFrame(t_results_df
                 .feature
                 .str
                 .split('_')
                 .values
                 .tolist(),
                 columns=['algorithm',
                          'feature_num',
                          "drop_z",
                          "z_dim",
                          "drop_seed",
                          "seed"])
    .merge(t_results_df,
           left_index=True,
           right_index=True)
    .drop(['drop_z',
           'drop_seed',
           'feature'],
          axis='columns')
)

full_heme_result_df.loc[:, 'z_dim'] = full_heme_result_df.z_dim.astype(np.int64)
full_heme_result_df.loc[:, 'feature_num'] = full_heme_result_df.feature_num.astype(np.int64)
full_heme_result_df.loc[:, 'seed'] = full_heme_result_df.seed.astype(np.int64)

full_heme_result_df = (
    full_heme_result_df
    .merge(monocyte_top_df,
           left_on=['algorithm', 'feature_num', 'z_dim', 'seed'],
           right_on=['algorithm', 'feature', 'z', 'seed'])
)

full_heme_result_df.head()


# In[25]:


file = os.path.join("results", "all_monocyte_top_scores_and_separation.tsv")
full_heme_result_df.to_csv(file, sep='\t', index=False)

