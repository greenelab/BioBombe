
# coding: utf-8

# ## Generating a Table of the BioBombe Interpreted Features in a High Scoring DAE model of TP53 inactivation
# 
# **Gregory Way, 2019**
# 
# I use the model previously identified that was used to predict TP53 inactivation.
# I observe the BioBombe gene set enrichment scores for the features with high coefficients in this model.

# In[1]:


import os
import sys
import pandas as pd


# ## Load the Top Model Identified Previously

# In[2]:


model_file = os.path.join("results", "top_dae_tp53_feature_for_followup.tsv")
top_model_df = pd.read_table(model_file)
top_model_df


# In[3]:


seed = top_model_df.seed.values[0]
z_dim = top_model_df.z_dim.values[0]
algorithm = top_model_df.algorithm.values[0].lower()


# ## Load the BioBombe network projection results for Cancer Hallmarks

# In[4]:


file = os.path.join("..", "6.analyze-weights", "results", "tcga",
                    "gph", "signal", "tcga_z_200_GpH__geneset_scores.tsv.gz")

scores_df = (
    pd
    .read_table(file)
    .query("seed == @seed")
    .query("z == @z_dim")
    .query("algorithm == @algorithm")
)

scores_df = (
    scores_df
    .assign(full_feature=scores_df.algorithm.astype(str) + "_" + scores_df.feature.astype(str),
            abs_z_score=scores_df.z_score.abs())
)
scores_df.head()


# ## Load Model Coefficients

# In[5]:


file = os.path.join("results", "mutation", "TP53", "TP53_coefficients.tsv.gz")

top_n_features = 8

coef_df = (
    pd.read_table(file)
    .query("seed == @seed")
    .query("z_dim == @z_dim")
    .query("algorithm == @algorithm")
    .query("signal == 'signal'")
    .sort_values(by='abs', ascending=False)
    .head(top_n_features)
)

use_features = coef_df.feature.tolist()
coef_df


# In[6]:


# Explore the biobombe scores for specific DAE features
top_n_features = 10

biobombe_df = (
    scores_df
    .query("full_feature in @use_features")
    .merge(coef_df,
           how='left',
           left_on=['full_feature', 'algorithm', 'seed'],
           right_on=['feature', 'algorithm', 'seed'])
    .drop(['model_type', 'feature_x', 'feature_y', 'signal'], axis='columns')
    .sort_values(by=['abs', 'abs_z_score'], ascending=False)
    .reset_index(drop=True)
)

top_biobombe_df = (
    biobombe_df
    .groupby('full_feature')
    .apply(func=lambda x: x.abs_z_score.nlargest(top_n_features))
    .reset_index()
    .merge(biobombe_df
           .reset_index(),
           right_on=['index', 'abs_z_score', 'full_feature'],
           left_on=['level_1', 'abs_z_score', 'full_feature'])
    .drop(['level_1', 'index'], axis='columns')
    .sort_values(by=['weight', 'z_score'], ascending=False)
)
    
    
print(top_biobombe_df.shape)
top_biobombe_df


# In[7]:


# Output biobombe scores applied to high scoring DAE features
file = os.path.join('results', 'tcga_tp53_classify_top_biobombe_scores_dae_table.tsv')
top_biobombe_df.to_csv(file, sep='\t', index=False)

