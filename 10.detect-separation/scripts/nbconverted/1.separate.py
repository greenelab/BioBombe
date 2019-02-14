
# coding: utf-8

# # Detect Separation Between Two Phenotypes in TCGA, GTEx, and TARGET Data
# 
# **Gregory Way, 2019**
# 
# Perform a t-test between two distinct phenotypes. In TCGA and GTEx data, we perform a t-test on males and females while in TARGET data, we test MYCN amplified vs. MYCN not amplified neuroblastoma tumors. We track the t-test p values across k dimensions and algorithms to isolate the features that best distinguishes the two groups.

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, ttest_rel
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("../8.gtex-interpret")
from scripts.utils import (
    load_weight_matrix,
    apply_signature,
    load_enrichment_results,
    extract_feature,
)

sys.path.append("../9.tcga-classify/")
from scripts.tcga_util import build_feature_dictionary


# In[2]:


def ttest_difference(feature_series, group_a_ids, group_b_ids):
    """
    To be applied to a pandas dataframe by column
    """
    feature_name = feature_series.name
    feature_algorithm, feature_num = feature_name.split('_')
    
    a_activation = feature_series[feature_series.index.isin(group_a_ids)]
    b_activation = feature_series[feature_series.index.isin(group_b_ids)]
    
    # Perform t-test on two groups
    t_stat, t_p = ttest_ind(a_activation, b_activation)
    
    return([t_stat, t_p, feature_algorithm, feature_num])


def get_ttest_results(z_matrix_dict, group_a_ids, group_b_ids, train_or_test='test'):
    """
    Loop through z matrix performing t-test using the compressed feature scores.
    Output full t-test results
    """
    
    # Perform t-test for all compressed features
    full_results = []
    for signal in z_matrix_dict.keys():
        for z_dim in z_matrix_dict[signal].keys():
            for seed in z_matrix_dict[signal][z_dim].keys():
                z_df = z_matrix_dict[signal][z_dim][seed][train_or_test]

                result_df = pd.DataFrame(z_df.apply(lambda x:
                                                    ttest_difference(feature_series=x,
                                                                     group_a_ids=group_a_ids,
                                                                     group_b_ids=group_b_ids)),
                                         columns = ['result'])

                result_df = (
                    pd.DataFrame(result_df.result.values.tolist(),
                                 columns=['t_stat', 't_p', 'algorithm', 'feature_num'])
                ).fillna(1)

                result_df = result_df.assign(
                    z_dim=z_dim,
                    signal=signal,
                    seed=seed
                )
                full_results.append(result_df)
    
    full_results_df = pd.concat(full_results)
    full_results_df = full_results_df.assign(neg_log_p=-np.log10(full_results_df.t_p))
    full_results_df = full_results_df.sort_values(by='neg_log_p', ascending=False)
    return full_results_df


# ## 1. GTEx Sex Analysis

# In[3]:


# Load GTEx phenotype data
file = os.path.join("..", "0.expression-download", "download", "GTEx_v7_Annotations_SubjectPhenotypesDS.txt")
gtex_pheno_df = pd.read_table(file)
gtex_pheno_df.head()


# In[4]:


gtex_z_matrix_dict = build_feature_dictionary(dataset="GTEX",
                                              load_data=True,
                                              store_train_test='test')


# In[5]:


# Extract male and female ids from the dataset using one matrix as an example
# (All matrices are aligned with the same IDs)
example_matrix_df = gtex_z_matrix_dict['signal']['8']['451283']['test']

patient_id_df = pd.concat(
    [
    pd.DataFrame(["{}-{}".format(x[0], x[1]) for x in example_matrix_df.index.str.split('-')],
                 columns=['patient_id'])
        .merge(gtex_pheno_df,
               how='left',
               left_on='patient_id',
               right_on='SUBJID'),
    pd.DataFrame(example_matrix_df.index)
    ],
    axis='columns'
)

gtex_males = patient_id_df.query("SEX == 1").sample_id.tolist()
gtex_females = patient_id_df.query("SEX == 2").sample_id.tolist()

print(patient_id_df.shape)
patient_id_df.head()


# In[6]:


# Perform t-test for all compressed features
gtex_full_results_df = get_ttest_results(z_matrix_dict=gtex_z_matrix_dict,
                                         group_a_ids=gtex_males,
                                         group_b_ids=gtex_females)


# In[7]:


# Output results
file = os.path.join("results", "sex_separation_gtex_t_test.tsv")
gtex_full_results_df.to_csv(file, sep='\t', index=False)

print(gtex_full_results_df.shape)
gtex_full_results_df.head()


# ## 2. TCGA Sex Analysis

# In[8]:


# Load TCGA phenotype data
file = os.path.join("..", "0.expression-download", "download", "TCGA-CDR-SupplementalTableS1.xlsx")
tcga_pheno_df = pd.read_excel(file)

tcga_pheno_df.head()


# In[9]:


tcga_z_matrix_dict = build_feature_dictionary(dataset="TCGA",
                                              load_data=True,
                                              store_train_test='test')


# In[10]:


# Extract male and female ids from the dataset
example_matrix_df = tcga_z_matrix_dict['signal']['2']['451283']['test']

patient_id_df = pd.concat(
    [
    pd.DataFrame(["{}-{}-{}".format(x[0], x[1], x[2]) for x in example_matrix_df.index.str.split('-')],
                 columns=['patient_id'])
        .merge(tcga_pheno_df,
               how='left',
               left_on='patient_id',
               right_on='bcr_patient_barcode'),
    pd.DataFrame(example_matrix_df.index)
    ],
    axis='columns'
)

tcga_males = patient_id_df.query("gender == 'MALE'").sample_id.tolist()
tcga_females = patient_id_df.query("gender == 'FEMALE'").sample_id.tolist()

print(patient_id_df.shape)
patient_id_df.head()


# In[11]:


# Perform t-test for all compressed features
tcga_full_results_df = get_ttest_results(z_matrix_dict=tcga_z_matrix_dict,
                                         group_a_ids=tcga_males,
                                         group_b_ids=tcga_females)


# In[12]:


# Output results
file = os.path.join("results", "sex_separation_tcga_t_test.tsv")
tcga_full_results_df.to_csv(file, sep='\t', index=False)

print(tcga_full_results_df.shape)
tcga_full_results_df.head()


# ## 3. TARGET NBL MYCN Status Analysis

# In[13]:


# Load TARGET phenotype data
file = os.path.join("..", "0.expression-download", "data", "2017-09-30-TARGET update harmonized.txt")
nbl_pheno_df = pd.read_table(file)
nbl_pheno_df.head()


# In[14]:


# Load TARGET matrices
target_z_matrix_dict = build_feature_dictionary(dataset="TARGET",
                                                load_data=True,
                                                store_train_test='train')


# In[15]:


# Extract male and female ids from the dataset
example_matrix_df = target_z_matrix_dict['signal']['8']['451283']['train']

patient_id_df = pd.concat(
    [
    pd.DataFrame([x[2] for x in example_matrix_df.index.str.split('-')],
                 columns=['patient_id'])
        .merge(nbl_pheno_df,
               how='left',
               left_on='patient_id',
               right_on='usi'),
    pd.DataFrame(example_matrix_df.index)
    ],
    axis='columns'
).dropna(subset=['usi'])

mycn_amp = patient_id_df.loc[patient_id_df["MYCN status"] == "Amplified", "sample_id"].tolist()
mycn_nonamp = patient_id_df.loc[patient_id_df["MYCN status"] == "Not Amplified", "sample_id"].tolist()

print(patient_id_df.shape)
patient_id_df.head()


# In[16]:


# Perform t-test for all compressed features
target_full_results_df = get_ttest_results(z_matrix_dict=target_z_matrix_dict,
                                           group_a_ids=mycn_amp,
                                           group_b_ids=mycn_nonamp,
                                           train_or_test='train')


# In[17]:


file = os.path.join("results", "nbl_mycn_separation_target_t_test.tsv")
target_full_results_df.to_csv(file, sep='\t', index=False)

print(target_full_results_df.shape)
target_full_results_df.head()

