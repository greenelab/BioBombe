"""
2018 Gregory Way
8.tcga-classify/classify-cancer-types.py

Predicting 33 different cancer-types using elastic net logistic regression and
compressed gene expression features in TCGA PanCanAtlas

Usage:

    python classify-cancer-types.py

Output:
Cancer-type specific DataFrames storing ROC, precision-recall, and classifier
coefficients for every compression model trained in their ability to predict
cancer-type.
"""

import os
import numpy as np
import pandas as pd

from scripts.tcga_util import (
    get_threshold_metrics,
    summarize_results,
    extract_coefficients,
    align_matrices,
    train_model,
    build_feature_dictionary,
    process_y_matrix_cancertype,
)

np.random.seed(123)

# Load constants
folds = 5
max_iter = 100
algorithms = ["pca", "ica", "nmf", "dae", "vae"]
signals = ["signal", "shuffled"]
alphas = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3]
l1_ratios = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4]

# Load data to build y matrices
base_url = "https://github.com/greenelab/pancancer/raw"
commit = "2a0683b68017fb226f4053e63415e4356191734f"

# Load data
file = "{}/{}/data/sample_freeze.tsv".format(base_url, commit)
sample_freeze_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/mutation_burden_freeze.tsv".format(base_url, commit)
mut_burden_df = pd.read_table(file, index_col=0)

# Track total metrics for each gene in one file
full_metrics_list = []
count_list = []

# Obtain a dictionary of file directories for loading each feature matrix (X)
z_matrix_dict = build_feature_dictionary()

# Provide one vs all classifications for all 33 different cancertypes in TCGA
for acronym in sample_freeze_df.DISEASE.unique():

    # Create list to store cancer-type specific results
    cancertype_auc_list = []
    cancertype_aupr_list = []
    cancertype_coef_list = []

    # Create directory for the cancer-type
    cancertype_dir = os.path.join("results", "cancer-type", acronym)
    os.makedirs(cancertype_dir, exist_ok=True)

    y_df, count_df = process_y_matrix_cancertype(
        acronym=acronym,
        sample_freeze=sample_freeze_df,
        mutation_burden=mut_burden_df,
        hyper_filter=5,
    )

    # Track the status counts of all classifiers
    count_list.append(count_df)

    # Now, perform all the analyses for each X matrix
    for signal in z_matrix_dict.keys():
        z_dim_dict = z_matrix_dict[signal]
        for z_dim in z_dim_dict.keys():
            seed_z_dim_dict = z_dim_dict[z_dim]
            for seed in seed_z_dim_dict.keys():
                z_train_file = z_matrix_dict[signal][z_dim][seed]["train"]
                z_test_file = z_matrix_dict[signal][z_dim][seed]["test"]

                for alg in algorithms:
                    # Load and process data
                    train_samples, x_train_df, y_train_df = align_matrices(
                        x_file=z_train_file,
                        y=y_df,
                        add_cancertype_covariate=False,
                        algorithm=alg,
                    )

                    test_samples, x_test_df, y_test_df = align_matrices(
                        x_file=z_test_file,
                        y=y_df,
                        add_cancertype_covariate=False,
                        algorithm=alg,
                    )

                    # Train the model
                    print(
                        "Training model... cancer-type: {}, "
                        "algorithm: {}, signal: {}, z_dim: {}, "
                        "seed: {}".format(acronym, alg, signal, z_dim, seed)
                    )

                    # Fit the model
                    cv_pipeline, y_pred_train_df, y_pred_test_df, y_cv_df = train_model(
                        x_train=x_train_df,
                        x_test=x_test_df,
                        y_train=y_train_df,
                        alphas=alphas,
                        l1_ratios=l1_ratios,
                        n_folds=folds,
                        max_iter=max_iter,
                    )
                    # Get metric  predictions
                    y_train_results = get_threshold_metrics(
                        y_train_df.status, y_pred_train_df, drop=False
                    )
                    y_test_results = get_threshold_metrics(
                        y_test_df.status, y_pred_test_df, drop=False
                    )
                    y_cv_results = get_threshold_metrics(
                        y_train_df.status, y_cv_df, drop=False
                    )

                    # Get coefficients
                    coef_df = extract_coefficients(
                        cv_pipeline=cv_pipeline,
                        feature_names=x_train_df.columns,
                        signal=signal,
                        z_dim=z_dim,
                        seed=seed,
                        algorithm=alg,
                    )

                    coef_df = coef_df.assign(acronym=acronym)

                    # Store all results
                    train_metrics_, train_roc_df, train_pr_df = summarize_results(
                        results=y_train_results,
                        gene_or_cancertype=acronym,
                        signal=signal,
                        z_dim=z_dim,
                        seed=seed,
                        algorithm=alg,
                        data_type="train",
                    )
                    test_metrics_, test_roc_df, test_pr_df = summarize_results(
                        results=y_test_results,
                        gene_or_cancertype=acronym,
                        signal=signal,
                        z_dim=z_dim,
                        seed=seed,
                        algorithm=alg,
                        data_type="test",
                    )
                    cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
                        results=y_cv_results,
                        gene_or_cancertype=acronym,
                        signal=signal,
                        z_dim=z_dim,
                        seed=seed,
                        algorithm=alg,
                        data_type="cv",
                    )

                    # Compile summary metrics
                    cols = [
                        "auroc",
                        "aupr",
                        "gene",
                        "signal",
                        "z_dim",
                        "seed",
                        "algorithm",
                        "data_type",
                    ]
                    metrics_ = [train_metrics_, test_metrics_, cv_metrics_]
                    metric_df_ = pd.DataFrame(metrics_, columns=cols)
                    full_metrics_list.append(metric_df_)

                    auc_df = pd.concat([train_roc_df, test_roc_df, cv_roc_df])
                    cancertype_auc_list.append(auc_df)

                    aupr_df = pd.concat([train_pr_df, test_pr_df, cv_pr_df])
                    cancertype_aupr_list.append(aupr_df)

                    cancertype_coef_list.append(coef_df)

    cancertype_auc_df = pd.concat(cancertype_auc_list)
    cancertype_aupr_df = pd.concat(cancertype_aupr_list)
    cancertype_coef_df = pd.concat(cancertype_coef_list)

    file = os.path.join(
        cancertype_dir, "{}_auc_threshold_metrics.tsv.gz".format(acronym)
    )
    cancertype_auc_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(
        cancertype_dir, "{}_aupr_threshold_metrics.tsv.gz".format(acronym)
    )
    cancertype_aupr_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(cancertype_dir, "{}_coefficients.tsv.gz".format(acronym))
    cancertype_coef_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

# Now, compile all results and write to file
final_metrics_df = pd.concat(full_metrics_list)

file = os.path.join("results", "all_cancertype_classify_metrics.tsv.gz")
final_metrics_df.to_csv(
    file, sep="\t", index=False, compression="gzip", float_format="%.5g"
)

# Write out a combined status matrix as well
final_count_list_df = pd.concat(count_list)

file = os.path.join("results", "all_cancertype_status_counts.tsv")
final_count_list_df.to_csv(file, sep="\t", index=False)
