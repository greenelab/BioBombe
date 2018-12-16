"""
2018 Gregory Way
9.tcga-classify/classify-with-raw-expression.py

Predict if specific genes are mutated across TCGA tumors based on raw RNAseq gene
expression features. Also make predictions on cancer types using raw gene expression.

Usage:

    python classify-with-raw-expression.py

Output:
Gene specific DataFrames storing ROC, precision-recall, and classifier coefficients
for raw gene expression models trained in their ability to predict mutations. The genes
used are the top 50 most mutated genes in TCGA PanCanAtlas. A gene was considered
mutated if a non-silent mutation was observed by the MC3 mutation calling effort. An
additional metrics file that stores all gene AUROC and AUPR is also saved. We also save
predictions and results for cancer-types.
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

from scripts.tcga_util import (
    get_threshold_metrics,
    summarize_results,
    extract_coefficients,
    align_matrices,
    process_y_matrix,
    train_model,
    process_y_matrix_cancertype,
    check_status,
)

np.random.seed(123)

# Load constants
filter_prop = 0.05
filter_count = 15
folds = 5
num_features = 8000
max_iter = 100
seed = "123"
algorithm = "raw"
alphas = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3]
l1_ratios = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4]

# Load genes
file = os.path.join("data", "top50_mutated_genes.tsv")
genes_df = pd.read_table(file)

# Load data to build x and y matrices
base_url = "https://github.com/greenelab/pancancer/raw"
commit = "2a0683b68017fb226f4053e63415e4356191734f"

file = "{}/{}/data/sample_freeze.tsv".format(base_url, commit)
sample_freeze_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/pancan_mutation_freeze.tsv.gz".format(base_url, commit)
mutation_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/copy_number_loss_status.tsv.gz".format(base_url, commit)
copy_loss_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/copy_number_gain_status.tsv.gz".format(base_url, commit)
copy_gain_df = pd.read_table(file, index_col=0)

file = "{}/{}/data/mutation_burden_freeze.tsv".format(base_url, commit)
mut_burden_df = pd.read_table(file, index_col=0)

# Load and process X matrix
base_dir = os.path.join("..", "0.expression-download", "data")
x_train_file = os.path.join(base_dir, "train_tcga_expression_matrix_processed.tsv.gz")
x_test_file = os.path.join(base_dir, "test_tcga_expression_matrix_processed.tsv.gz")

rnaseq_train_df = pd.read_table(x_train_file, index_col=0)
rnaseq_test_df = pd.read_table(x_test_file, index_col=0)

# Determine most variably expressed genes and subset X matrix
mad_file = os.path.join(base_dir, "tcga_mad_genes.tsv")
mad_genes_df = pd.read_table(mad_file)
mad_genes = mad_genes_df.iloc[0:num_features, ].gene_id.astype(str)

rnaseq_train_df = rnaseq_train_df.reindex(mad_genes, axis="columns")
rnaseq_test_df = rnaseq_test_df.reindex(mad_genes, axis="columns")

# Scale RNAseq matrix the same way RNAseq was scaled for compression algorithms
train_fitted_scaler = MinMaxScaler().fit(rnaseq_train_df)
rnaseq_train_df = pd.DataFrame(
    train_fitted_scaler.transform(rnaseq_train_df),
    columns=rnaseq_train_df.columns,
    index=rnaseq_train_df.index,
)

test_fitted_scaler = MinMaxScaler().fit(rnaseq_test_df)
rnaseq_test_df = pd.DataFrame(
    test_fitted_scaler.transform(rnaseq_test_df),
    columns=rnaseq_test_df.columns,
    index=rnaseq_test_df.index,
)

# Track total metrics for each gene in one file
metric_cols = [
    "auroc",
    "aupr",
    "gene_or_cancertype",
    "signal",
    "z_dim",
    "seed",
    "algorithm",
    "data_type",
]

for gene_idx, gene_series in genes_df.iterrows():

    gene_name = gene_series.gene
    classification = gene_series.classification

    # Create list to store gene specific results
    gene_auc_list = []
    gene_aupr_list = []
    gene_coef_list = []
    gene_metrics_list = []

    # Create directory for the gene
    gene_dir = os.path.join("results", "mutation", gene_name)
    os.makedirs(gene_dir, exist_ok=True)

    # Check if gene has been processed already
    check_file = os.path.join(gene_dir, "{}_raw_coefficients.tsv.gz".format(gene_name))
    if check_status(check_file):
        continue
    # Process the y matrix for the given gene or pathway
    y_mutation_df = mutation_df.loc[:, gene_name]

    # Include copy number gains for oncogenes
    # and copy number loss for tumor suppressor genes (TSG)
    include_copy = True
    if classification == "Oncogene":
        y_copy_number_df = copy_gain_df.loc[:, gene_name]
    elif classification == "TSG":
        y_copy_number_df = copy_loss_df.loc[:, gene_name]
    else:
        y_copy_number_df = pd.DataFrame()
        include_copy = False

    y_df = process_y_matrix(
        y_mutation=y_mutation_df,
        y_copy=y_copy_number_df,
        include_copy=include_copy,
        gene=gene_name,
        sample_freeze=sample_freeze_df,
        mutation_burden=mut_burden_df,
        filter_count=filter_count,
        filter_prop=filter_prop,
        output_directory=gene_dir,
        hyper_filter=5,
    )

    for signal in ["signal", "shuffled"]:
        if signal == "shuffled":
            # Shuffle training data
            x_train_raw_df = rnaseq_train_df.apply(
                lambda x: np.random.permutation(x.tolist()),
                axis=1,
                result_type="expand",
            )

            x_train_raw_df.columns = rnaseq_train_df.columns
            x_train_raw_df.index = rnaseq_train_df.index

            # Shuffle testing data
            x_test_raw_df = rnaseq_test_df.apply(
                lambda x: np.random.permutation(x.tolist()),
                axis=1,
                result_type="expand",
            )

            x_test_raw_df.columns = rnaseq_test_df.columns
            x_test_raw_df.index = rnaseq_test_df.index

        else:
            x_train_raw_df = rnaseq_train_df
            x_test_raw_df = rnaseq_test_df

        # Now, perform all the analyses for each X matrix
        train_samples, x_train_df, y_train_df = align_matrices(
            x_file_or_df=x_train_raw_df, y=y_df
        )

        test_samples, x_test_df, y_test_df = align_matrices(
            x_file_or_df=x_test_raw_df, y=y_df
        )

        # Train the model
        print(
            "Training model... gene: {}, for raw {} features".format(gene_name, signal)
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
        y_cv_results = get_threshold_metrics(y_train_df.status, y_cv_df, drop=False)

        # Get coefficients
        coef_df = extract_coefficients(
            cv_pipeline=cv_pipeline,
            feature_names=x_train_df.columns,
            signal=signal,
            z_dim=num_features,
            seed=seed,
            algorithm=algorithm,
        )

        coef_df = coef_df.assign(gene=gene_name)

        # Store all results
        train_metrics_, train_roc_df, train_pr_df = summarize_results(
            y_train_results, gene_name, signal, num_features, seed, algorithm, "train"
        )
        test_metrics_, test_roc_df, test_pr_df = summarize_results(
            y_test_results, gene_name, signal, num_features, seed, algorithm, "test"
        )
        cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
            y_cv_results, gene_name, signal, num_features, seed, algorithm, "cv"
        )

        # Compile summary metrics
        metrics_ = [train_metrics_, test_metrics_, cv_metrics_]
        metric_df_ = pd.DataFrame(metrics_, columns=metric_cols)
        gene_metrics_list.append(metric_df_)

        gene_auc_df = pd.concat([train_roc_df, test_roc_df, cv_roc_df])
        gene_auc_list.append(gene_auc_df)

        gene_aupr_df = pd.concat([train_pr_df, test_pr_df, cv_pr_df])
        gene_aupr_list.append(gene_aupr_df)

        gene_coef_list.append(coef_df)

    gene_auc_df = pd.concat(gene_auc_list)
    gene_aupr_df = pd.concat(gene_aupr_list)
    gene_coef_df = pd.concat(gene_coef_list)
    gene_metrics_df = pd.concat(gene_metrics_list)

    file = os.path.join(
        gene_dir, "{}_raw_auc_threshold_metrics.tsv.gz".format(gene_name)
    )
    gene_auc_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(
        gene_dir, "{}_raw_aupr_threshold_metrics.tsv.gz".format(gene_name)
    )
    gene_aupr_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    gene_coef_df.to_csv(
        check_file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join("results", "{}_raw_classify_metrics.tsv.gz".format(gene_name))
    gene_metrics_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

# Provide one vs all classifications for all 33 different cancertypes in TCGA
# Track total metrics for each cancer-type in one file
count_list = []

for acronym in sample_freeze_df.DISEASE.unique():

    # Create list to store cancer-type specific results
    cancertype_auc_list = []
    cancertype_aupr_list = []
    cancertype_coef_list = []
    cancertype_metrics_list = []

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

    # Check if cancer-type has been processed already
    check_file = (
        os.path.join(cancertype_dir, "{}_raw_coefficients.tsv.gz".format(acronym))
    )
    if check_status(check_file):
        continue

    # Now, perform all the analyses for each X matrix
    for signal in ["signal", "shuffled"]:
        if signal == "shuffled":
            # Shuffle training data
            x_train_raw_df = rnaseq_train_df.apply(
                lambda x: np.random.permutation(x.tolist()),
                axis=1,
                result_type="expand",
            )

            x_train_raw_df.columns = rnaseq_train_df.columns
            x_train_raw_df.index = rnaseq_train_df.index

            # Shuffle testing data
            x_test_raw_df = rnaseq_test_df.apply(
                lambda x: np.random.permutation(x.tolist()),
                axis=1,
                result_type="expand",
            )

            x_test_raw_df.columns = rnaseq_test_df.columns
            x_test_raw_df.index = rnaseq_test_df.index

        else:
            x_train_raw_df = rnaseq_train_df
            x_test_raw_df = rnaseq_test_df

        # Now, perform all the analyses for each X matrix
        train_samples, x_train_df, y_train_df = align_matrices(
            x_file_or_df=x_train_raw_df,
            y=y_df,
            add_cancertype_covariate=False
        )

        test_samples, x_test_df, y_test_df = align_matrices(
            x_file_or_df=x_test_raw_df,
            y=y_df,
            add_cancertype_covariate=False
        )

        # Train the model
        print(
            "Training model... cancertype: {}, for raw {} features".format(acronym,
                                                                           signal)
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
        y_cv_results = get_threshold_metrics(y_train_df.status, y_cv_df, drop=False)

        # Get coefficients
        coef_df = extract_coefficients(
            cv_pipeline=cv_pipeline,
            feature_names=x_train_df.columns,
            signal=signal,
            z_dim=num_features,
            seed=seed,
            algorithm=algorithm,
        )

        coef_df = coef_df.assign(acronym=acronym)

        # Store all results
        train_metrics_, train_roc_df, train_pr_df = summarize_results(
            results=y_train_results,
            gene_or_cancertype=acronym,
            signal=signal,
            z_dim=num_features,
            seed=seed,
            algorithm=algorithm,
            data_type="train",
        )
        test_metrics_, test_roc_df, test_pr_df = summarize_results(
            results=y_test_results,
            gene_or_cancertype=acronym,
            signal=signal,
            z_dim=num_features,
            seed=seed,
            algorithm=algorithm,
            data_type="test",
        )
        cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
            results=y_cv_results,
            gene_or_cancertype=acronym,
            signal=signal,
            z_dim=num_features,
            seed=seed,
            algorithm=algorithm,
            data_type="cv"
        )

        # Compile summary metrics
        metrics_ = [train_metrics_, test_metrics_, cv_metrics_]
        metric_df_ = pd.DataFrame(metrics_, columns=metric_cols)
        cancertype_metrics_list.append(metric_df_)

        auc_df = pd.concat([train_roc_df, test_roc_df, cv_roc_df])
        cancertype_auc_list.append(auc_df)

        aupr_df = pd.concat([train_pr_df, test_pr_df, cv_pr_df])
        cancertype_aupr_list.append(aupr_df)

        cancertype_coef_list.append(coef_df)

    cancertype_auc_df = pd.concat(cancertype_auc_list)
    cancertype_aupr_df = pd.concat(cancertype_aupr_list)
    cancertype_coef_df = pd.concat(cancertype_coef_list)
    cancertype_metrics_df = pd.concat(cancertype_metrics_list)

    file = os.path.join(
        cancertype_dir, "{}_raw_auc_threshold_metrics.tsv.gz".format(acronym)
    )
    cancertype_auc_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(
        cancertype_dir, "{}_raw_aupr_threshold_metrics.tsv.gz".format(acronym)
    )
    cancertype_aupr_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    cancertype_coef_df.to_csv(
        check_file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(
        cancertype_dir, "{}_raw_classify_metrics.tsv.gz".format(acronym)
    )
    cancertype_metrics_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

# Write out a combined status matrix as well
final_count_list_df = pd.concat(count_list)

file = os.path.join("results", "all_cancertype_status_counts.tsv")
final_count_list_df.to_csv(file, sep="\t", index=False)
