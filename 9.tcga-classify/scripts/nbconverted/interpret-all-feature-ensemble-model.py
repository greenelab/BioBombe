
# coding: utf-8

# ## Interpreting Ensemble Compressed Features
# 
# **Gregory Way, 2019**
# 
# The following notebook will assign biological knowledge to the compressed features using the network projection approach. I use the model previously identified that was used to predict TP53 inactivation.
# I observe the BioBombe gene set enrichment scores for the features with high coefficients in this model.

# In[1]:


import os
import sys
import pandas as pd


# ## Load the `All Feature` Ensemble Model

# In[2]:


model_file = os.path.join("results", "top_model_ensemble_all_features_tp53_feature_for_followup.tsv")
top_model_df = pd.read_table(model_file)
top_model_df


# In[3]:


coef_file = os.path.join("results",
                         "mutation_ensemble_all",
                         "TP53",
                         "TP53_ensemble_all_features_coefficients.tsv.gz")
coef_df = pd.read_table(coef_file).drop(['signal', 'z_dim', 'seed', 'algorithm'], axis='columns')
coef_df.head()


# In[4]:


full_coef_id_df = (
    pd.DataFrame(coef_df.feature.str.split("_").values.tolist(),
                 columns=['algorithm', 'individual_feature', 'seed', 'k', 'signal'])
)

full_coef_id_df = pd.concat([full_coef_id_df, coef_df], axis='columns')
full_coef_id_df = full_coef_id_df.query("abs > 0").query("signal == 'signal'")

print(full_coef_id_df.shape)
full_coef_id_df.head()


# ## Load Network Projection Results

# In[5]:


gph_dir = os.path.join("..",
                       "6.biobombe-projection",
                       "results",
                       "tcga",
                       "gph",
                       "signal")
gph_files = os.listdir(gph_dir)


# In[6]:


all_scores_list = []
for file in gph_files:
    file = os.path.join(gph_dir, file)
    scores_df = pd.read_table(file)
    all_scores_list.append(scores_df)


# In[7]:


all_scores_df = pd.concat(all_scores_list, axis='rows')

print(all_scores_df.shape)
all_scores_df.head()


# In[8]:


all_scores_df = all_scores_df.assign(big_feature_id=all_scores_df.algorithm + "_" +
                                     all_scores_df.feature.astype(str) + "_" +
                                     all_scores_df.seed.astype(str) + "_" +
                                     all_scores_df.z.astype(str) + "_signal")
all_scores_df = all_scores_df.assign(abs_z_score=all_scores_df.z_score.abs())


# In[9]:


all_coef_scores_df = (
    full_coef_id_df
    .merge(all_scores_df,
           how='left',
           left_on="feature",
           right_on="big_feature_id")
    .sort_values(by=['abs', 'abs_z_score'], ascending=False)
    .reset_index(drop=True)
)

all_coef_scores_df.head()


# In[10]:


# Explore the biobombe scores for specific DAE features
top_n_features = 5

biobombe_df = (
    all_coef_scores_df
    .groupby('big_feature_id')
    .apply(func=lambda x: x.abs_z_score.nlargest(top_n_features))
    .reset_index()
    .merge(all_coef_scores_df
           .reset_index(),
           right_on=['index', 'abs_z_score', 'big_feature_id'],
           left_on=['level_1', 'abs_z_score', 'big_feature_id'])
    .drop(['level_1', 'index', 'feature_x',
           'algorithm_x', 'seed_x',
           'model_type', 'algorithm_y',
           'feature_y', 'seed_y', 'z'], axis='columns')
    .sort_values(by=['abs', 'abs_z_score'], ascending=False)
    .reset_index(drop=True)
)

print(biobombe_df.shape)
biobombe_df.head(20)


# In[11]:


# Output biobombe scores applied to the all feature ensemble model
file = os.path.join('results', 'tcga_tp53_classify_top_biobombe_scores_all_feature_ensemble_model_table.tsv')
biobombe_df.to_csv(file, sep='\t', index=False)


# ## Detect the highest contributing variables

# In[12]:


neg_biobombe_df = biobombe_df.query("weight < 0")
pos_biobombe_df = biobombe_df.query("weight > 0")

top_neg_variables_df = neg_biobombe_df.groupby("variable")['weight'].sum().sort_values(ascending=True)
top_pos_variables_df = pos_biobombe_df.groupby("variable")['weight'].sum().sort_values(ascending=False)


# In[13]:


full_result_df = pd.DataFrame(pd.concat([top_pos_variables_df, top_neg_variables_df]))
full_result_df = (
    full_result_df
    .assign(abs_weight=full_result_df.weight.abs())
    .sort_values(by='abs_weight', ascending=False)
)

full_result_df.head()


# In[14]:


# Output biobombe scores applied to the all feature ensemble model
file = os.path.join('results', 'tcga_tp53_classify_aggregate_biobombe_scores_all_feature_ensemble.tsv')
full_result_df.to_csv(file, sep='\t', index=False)

