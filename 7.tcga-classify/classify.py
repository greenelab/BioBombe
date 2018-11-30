"""
2018 Gregory Way
7.tcga-classify/classify.py

Predict Ras and TP53 status in TCGA based on compressed expression features

Usage:

    python classify.py

Output:
Three DataFrames storing ROC, precision-recall, and classifier coefficients
for every compression model trained in their ability to predict Ras pathway
activation and TP53 pathway inactivation
"""

import os
import glob
import numpy as np
import pandas as pd

from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import cross_val_predict
from dask_searchcv import GridSearchCV
from sklearn.pipeline import Pipeline

from scripts.tcga_util import get_threshold_metrics, summarize_results

np.random.seed(123)

# Load constants
genes = ['TP53', 'RAS']
filter_prop = 0.05
filter_count = 15
folds = 5
algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']
signals = ['signal', 'shuffled']
alphas = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3]
l1_ratios = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4]

# Load data to build y matrices
base_url = 'https://github.com/greenelab/pancancer/raw'
commit = '9fd9afbecdbb4f855ecc85bb282fc59e75c7744f'

# Load data
file = '../../pancancer/data/sample_freeze.tsv'
sample_freeze_df = pd.read_table(file, index_col=0)

file = '{}/{}/data/pancan_mutation_freeze.tsv.gz'.format(base_url, commit)
mutation_df = pd.read_table(file, index_col=0)

file = '{}/{}/data/copy_number_loss_status.tsv.gz'.format(base_url, commit)
copy_loss_df = pd.read_table(file, index_col=0)

file = '{}/{}/data/copy_number_gain_status.tsv.gz'.format(base_url, commit)
copy_gain_df = pd.read_table(file, index_col=0)

file = '{}/{}/data/mutation_burden_freeze.tsv'.format(base_url, commit)
mut_burden_df = pd.read_table(file, index_col=0)

# What results to track
full_metrics_list = []
full_auc_list = []
full_aupr_list = []
full_coef_list = []
all_failed_params = []

# Obtain a list of locations for each feature matrix (X)
z_matrix_dict = {}
for signal in signals:
    z_matrix_dict[signal] = {}

    if signal == 'signal':
        results_dir = 'TCGA_results'
    else:
        results_dir = 'TCGA_shuffled_results'

    matrix_dir = os.path.join('..', '2.ensemble-z-analysis', 'results',
                              results_dir, 'ensemble_z_matrices')

    for comp_dir in os.listdir(matrix_dir):
        matrix_comp_dir = os.path.join(matrix_dir, comp_dir)
        z_dim = comp_dir.split('_')[2]
        z_matrix_dict[signal][z_dim] = {}

        for z_file in glob.glob('{}/*_z_*'.format(matrix_comp_dir)):
            seed = os.path.basename(z_file).split('_')[1]

            if seed not in z_matrix_dict[signal][z_dim].keys():
                z_matrix_dict[signal][z_dim][seed] = {}

            if '_test_' in z_file:
                z_matrix_dict[signal][z_dim][seed]['test'] = z_file
            else:
                z_matrix_dict[signal][z_dim][seed]['train'] = z_file

for gene in genes:
    # Process the y matrix for the given gene or pathway
    if gene == 'RAS':
        subset_gene = ['KRAS', 'HRAS', 'NRAS']
        y_copy_number_df = copy_gain_df.loc[:, subset_gene].max(axis='columns')
        y_mutation_df = mutation_df.loc[:, subset_gene].max(axis='columns')
    else:
        y_copy_number_df = copy_loss_df.loc[:, gene]
        y_mutation_df = mutation_df.loc[:, gene]

    y_df = y_copy_number_df + y_mutation_df
    y_df.loc[y_df > 1] = 1
    y_df = pd.DataFrame(y_df)
    y_df.columns = ['status']

    y_df = (
        y_df
        .merge(sample_freeze_df, how='left', left_index=True,
               right_on='SAMPLE_BARCODE')
        .set_index('SAMPLE_BARCODE')
        .merge(mut_burden_df, left_index=True, right_index=True)
    )

    # Get statistics per gene and disease
    disease_counts_df = pd.DataFrame(y_df.groupby('DISEASE').sum()['status'])

    disease_proportion_df = (
        disease_counts_df
        .divide(y_df['DISEASE']
                .value_counts(sort=False)
                .sort_index(), axis=0)
    )

    # Filter diseases with low counts or proportions for classification balance
    filter_disease_df = ((disease_counts_df > filter_count) &
                         (disease_proportion_df > filter_prop))
    filter_disease_df.columns = ['disease_included']

    disease_stats_df = (
        disease_counts_df
        .merge(disease_proportion_df,
               left_index=True,
               right_index=True,
               suffixes=('_count', '_proportion'))
        .merge(filter_disease_df,
               left_index=True,
               right_index=True)
    )

    filter_file = '{}_filtered_cancertypes.tsv'.format(gene)
    filter_file = os.path.join('results', filter_file)
    disease_stats_df.to_csv(filter_file, sep='\t')

    # Filter
    use_diseases = disease_stats_df.query("disease_included").index.tolist()
    burden_filter = y_df['log10_mut'] < 5 * y_df['log10_mut'].std()
    y_df = y_df.loc[burden_filter, :].query("DISEASE in @use_diseases")

    # Now, perform all the analyses for each X matrix
    for signal in z_matrix_dict.keys():
        z_dim_dict = z_matrix_dict[signal]
        for z_dim in z_dim_dict.keys():
            seed_z_dim_dict = z_dim_dict[z_dim]
            for seed in seed_z_dim_dict.keys():
                z_train_file = z_matrix_dict[signal][z_dim][seed]['train']
                z_test_file = z_matrix_dict[signal][z_dim][seed]['test']

                for alg in algorithms:
                    # Load Data
                    x_train_df = pd.read_table(z_train_file, index_col=0)
                    use_col = x_train_df.columns.str.contains(alg)
                    x_train_df = x_train_df.loc[:, use_col]

                    x_test_df = pd.read_table(z_test_file, index_col=0)
                    use_col = x_test_df.columns.str.contains(alg)
                    x_test_df = x_test_df.loc[:, use_col]

                    # Subset samples
                    train_samples = (
                        set(y_df.index).intersection(set(x_train_df.index))
                    )
                    test_samples = (
                        set(y_df.index).intersection(set(x_test_df.index))
                    )

                    x_train_df = x_train_df.reindex(train_samples)
                    y_train_df = y_df.reindex(train_samples)

                    x_test_df = x_test_df.reindex(test_samples)
                    y_test_df = y_df.reindex(test_samples)

                    # Add in covariate info
                    covar_train_df = pd.get_dummies(y_train_df.DISEASE)
                    covar_test_df = pd.get_dummies(y_test_df.DISEASE)

                    mut_covar_train_df = (
                        pd.DataFrame(y_train_df.loc[:, 'log10_mut'],
                                     index=y_train_df.index)
                        )
                    mut_covar_test_df = (
                        pd.DataFrame(y_test_df.loc[:, 'log10_mut'],
                                     index=y_test_df.index)
                        )

                    x_train_df = (
                        x_train_df
                        .merge(covar_train_df, left_index=True,
                               right_index=True)
                        .merge(mut_covar_train_df, left_index=True,
                               right_index=True)
                    )
                    x_test_df = (
                        x_test_df
                        .merge(covar_test_df, left_index=True,
                               right_index=True)
                        .merge(mut_covar_test_df, left_index=True,
                               right_index=True)
                    )

                    # Setup the classifier parameters
                    clf_parameters = {'classify__loss': ['log'],
                                      'classify__penalty': ['elasticnet'],
                                      'classify__alpha': alphas,
                                      'classify__l1_ratio': l1_ratios}

                    estimator = (
                        Pipeline(
                            steps=[('classify',
                                    SGDClassifier(random_state=0,
                                                  class_weight='balanced',
                                                  loss='log',
                                                  max_iter=100,
                                                  tol=1e-3))]
                                )
                    )

                    cv_pipeline = GridSearchCV(estimator=estimator,
                                               param_grid=clf_parameters,
                                               n_jobs=-1,
                                               cv=folds,
                                               scoring='roc_auc',
                                               return_train_score=True)

                    print('Training model... gene: {}, '
                          'algorithm: {}, signal: {}, z_dim: {}, '
                          'seed: {}'.format(gene, alg, signal, z_dim, seed))

                    # Fit the model
                    try:
                        output_cv = cv_pipeline.fit(X=x_train_df,
                                                    y=y_train_df.status)
                    except ValueError:
                        failed_params = [gene, alg, signal, z_dim, seed]
                        all_failed_params.append(failed_params)
                        continue

                    # Obtain cross validation results
                    y_cv_df = (
                        cross_val_predict(cv_pipeline.best_estimator_,
                                          X=x_train_df,
                                          y=y_train_df.status,
                                          cv=folds,
                                          method='decision_function')
                    )

                    # Get all performance results
                    y_predict_train_df = (
                        cv_pipeline.decision_function(x_train_df)
                    )
                    y_predict_test_df = (
                        cv_pipeline.decision_function(x_test_df)
                    )

                    # Get metric  predictions
                    y_train_results = (
                        get_threshold_metrics(y_train_df.status,
                                              y_predict_train_df,
                                              drop=False)
                    )
                    y_test_results = (
                        get_threshold_metrics(y_test_df.status,
                                              y_predict_test_df,
                                              drop=False)
                    )
                    y_cv_results = (
                        get_threshold_metrics(y_train_df.status,
                                              y_cv_df,
                                              drop=False)
                    )

                    # Get coefficients
                    final_pipeline = cv_pipeline.best_estimator_
                    final_classifier = final_pipeline.named_steps['classify']

                    coef_df = pd.DataFrame.from_dict(
                        {'feature': x_train_df.columns,
                         'weight': final_classifier.coef_[0]})

                    coef_df = (
                        coef_df
                        .assign(abs=coef_df['weight'].abs())
                        .sort_values('abs', ascending=False)
                        .reset_index(drop=True)
                        .assign(gene=gene,
                                signal=signal,
                                z_dim=z_dim,
                                seed=seed,
                                algorithm=alg)
                    )

                    # Store all results
                    train_metrics_, train_roc_df, train_pr_df = (
                        summarize_results(y_train_results, gene, signal, z_dim,
                                          seed, alg, 'train')
                    )
                    test_metrics_, test_roc_df, test_pr_df = (
                        summarize_results(y_test_results, gene, signal, z_dim,
                                          seed, alg, 'test')
                    )
                    cv_metrics_, cv_roc_df, cv_pr_df = (
                        summarize_results(y_cv_results, gene, signal, z_dim,
                                          seed, alg, 'cv')
                    )

                    # Compile summary metrics
                    cols = ['auroc', 'aupr', 'gene', 'signal', 'z_dim', 'seed',
                            'algorithm', 'data_type']
                    metrics_ = [train_metrics_, test_metrics_, cv_metrics_]
                    metric_df_ = pd.DataFrame(metrics_, columns=cols)
                    full_metrics_list.append(metric_df_)

                    full_auc_df = pd.concat(
                        [train_roc_df, test_roc_df, cv_roc_df]
                    )
                    full_auc_list.append(full_auc_df)

                    full_aupr_df = pd.concat(
                        [train_pr_df, test_pr_df, cv_pr_df]
                    )
                    full_aupr_list.append(full_aupr_df)
                    full_coef_list.append(coef_df)

# Now, compile all results and write to file
final_metrics_df = pd.concat(full_metrics_list)
final_auc_df = pd.concat(full_auc_list)
final_aupr_df = pd.concat(full_aupr_list)
final_coef_df = pd.concat(full_coef_list)

file = os.path.join('results', 'classify_metrics.tsv')
final_metrics_df.to_csv(file, sep='\t', index=False)

file = os.path.join('results', 'auc_threshold_metrics.tsv.gz')
final_auc_df.to_csv(file, sep='\t', index=False, compression='gzip',
                    float_format='%.5g')

file = os.path.join('results', 'aupr_threshold_metrics.tsv.gz')
final_aupr_df.to_csv(file, sep='\t', index=False, compression='gzip',
                     float_format='%.5g')

file = os.path.join('results', 'coefficients.tsv.gz')
final_coef_df.to_csv(file, sep='\t', index=False, compression='gzip',
                     float_format='%.5g')

failed_params_df = pd.DataFrame(failed_params)
file = os.path.join('results', 'failed_params.tsv')
failed_params_df.to_csv(file, sep='\t')
