"""
2018 Gregory Way
9.tcga-classify/classify-with-ensemble-models.py

Predicting 33 different cancer-types using elastic net logistic regression and
compressed gene expression features in TCGA PanCanAtlas using different ensemble models

Usage:

    python classify-with-ensemble-models.py

and optional command arugments:

    --algorithm <the algorithm to build ensemble matrices>
        if None, then use only
    --use_all_features <boolean if all features used across algorithms and k>

Output:
Cancer-type specific DataFrames storing ROC, precision-recall, and classifier
coefficients for every compression model trained in their ability to predict
cancer-type.
"""
import os
import argparse
import numpy as np
import pandas as pd

from scripts.tcga_util import (
    get_threshold_metrics,
    summarize_results,
    extract_coefficients,
    check_status,
    align_matrices,
    process_y_matrix,
    process_y_matrix_cancertype,
    train_model,
    build_top_feature_dictionary,
    get_feature,
    load_ensemble_dict,
)

# Load command arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "-a", "--algorithm", help="which algorithm to subset", default="vae"
)
parser.add_argument(
    "-u", "--use_all_features", action="store_true", help="boolean if all features"
)
args = parser.parse_args()

np.random.seed(123)

# Set arguments and constants
algorithm = args.algorithm
if algorithm not in ["pca", "ica", "nmf", "dae", "vae"]:
    algorithm = None
use_all_features = args.use_all_features
filter_prop = 0.05
filter_count = 15
folds = 5
max_iter = 100
alphas = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3]
l1_ratios = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4]
zs = [
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    12,
    14,
    16,
    18,
    20,
    25,
    30,
    35,
    40,
    45,
    50,
    60,
    70,
    80,
    90,
    100,
    125,
    150,
    200,
]

genes = ["TP53", "PTEN", "PIK3CA", "KRAS", "TTN"]
# cancertypes = ["KIRP", "OV", "UCEC", "LUAD", "BRCA"]

# Only use the first seed if the algorithm is None
if algorithm:
    seeds = ["165158", "451283", "486191", "908341", "978124"]
else:
    seeds = ["165158"]

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

# Load gene info
file = os.path.join("data", "top50_mutated_genes.tsv")
genes_df = pd.read_table(file)

# Load data to build y matrices
base_url = "https://github.com/greenelab/pancancer/raw"
commit = "2a0683b68017fb226f4053e63415e4356191734f"

# Load data
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

# Retrieve algorithm dictionary
algorithm_dict = load_ensemble_dict(
    zs=zs, seeds=seeds, algorithm=algorithm, use_all_features=use_all_features
)

# Perform classification on genes
for gene_idx, gene_series in genes_df.iterrows():

    gene_name = gene_series.gene
    classification = gene_series.classification

    if gene_name not in genes:
        continue

    # Create list to store gene specific results
    gene_auc_list = []
    gene_aupr_list = []
    gene_coef_list = []
    gene_metrics_list = []

    # Create directory for the gene
    gene_dir = os.path.join("results", "mutation_ensemble", gene_name)
    os.makedirs(gene_dir, exist_ok=True)

    # Get gene info to build y matrix
    gene_info = genes_df.query("gene == @gene_name")
    classification = gene_info.classification.values[0]

    # Check if gene has been processed already
    if algorithm:
        check_file = os.path.join(gene_dir, "{}_coefficients.tsv.gz".format(gene_name))
    else:
        update_name = "{}_ensemble_all_alg".format(gene_name)
        check_file = os.path.join(
            gene_dir, "{}_coefficients.tsv.gz".format(update_name)
        )

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

    # Now, perform all the analyses for each X matrix
    for signal in algorithm_dict.keys():
        k_dim_dict = algorithm_dict[signal]
        for k in k_dim_dict.keys():
            z_train_file = k_dim_dict[k]["train"]
            z_test_file = k_dim_dict[k]["test"]

            # Load and process data
            train_samples, x_train_df, y_train_df = align_matrices(
                x_file_or_df=z_train_file, y=y_df, algorithm=None
            )

            test_samples, x_test_df, y_test_df = align_matrices(
                x_file_or_df=z_test_file, y=y_df, algorithm=None
            )

            # Train the model
            print(
                "Training model... gene: {}, "
                "ensemble algorithm: {}, signal: {}, z_dim: {}, ".format(
                    gene_name, algorithm, signal, k
                )
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
            if algorithm:
                algorithm_assign = "{}_ensemble".format(algorithm)
            else:
                algorithm_assign = "all_ensemble"

            coef_df = extract_coefficients(
                cv_pipeline=cv_pipeline,
                feature_names=x_train_df.columns,
                signal=signal,
                z_dim=k,
                seed=",".join(seeds),
                algorithm=algorithm_assign,
            )

            coef_df = coef_df.assign(gene=gene_name)

            # Store all results
            train_metrics_, train_roc_df, train_pr_df = summarize_results(
                results=y_train_results,
                gene_or_cancertype=gene_name,
                signal=signal,
                z_dim=k,
                seed="ensemble",
                algorithm=algorithm_assign,
                data_type="train",
            )
            test_metrics_, test_roc_df, test_pr_df = summarize_results(
                results=y_test_results,
                gene_or_cancertype=gene_name,
                signal=signal,
                z_dim=k,
                seed="ensemble",
                algorithm=algorithm_assign,
                data_type="test",
            )
            cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
                results=y_cv_results,
                gene_or_cancertype=gene_name,
                signal=signal,
                z_dim=k,
                seed="ensemble",
                algorithm=algorithm_assign,
                data_type="cv",
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

    if not algorithm:
        gene_name = "{}_ensemble_all_alg".format(gene_name)

    file = os.path.join(gene_dir, "{}_auc_threshold_metrics.tsv.gz".format(gene_name))
    gene_auc_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(gene_dir, "{}_aupr_threshold_metrics.tsv.gz".format(gene_name))
    gene_aupr_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    gene_coef_df.to_csv(
        check_file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(gene_dir, "{}_classify_metrics.tsv.gz".format(gene_name))
    gene_metrics_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

# Provide one vs all classifications for all 33 different cancertypes in TCGA
for acronym in sample_freeze_df.DISEASE.unique():

    # if acronym not in cancertypes:
    #    continue

    # Create list to store cancer-type specific results
    cancertype_auc_list = []
    cancertype_aupr_list = []
    cancertype_coef_list = []
    cancertype_metrics_list = []

    # Create directory for the cancer-type
    cancertype_dir = os.path.join("results", "cancer-type_ensemble", acronym)
    os.makedirs(cancertype_dir, exist_ok=True)

    y_df, count_df = process_y_matrix_cancertype(
        acronym=acronym,
        sample_freeze=sample_freeze_df,
        mutation_burden=mut_burden_df,
        hyper_filter=5,
    )

    # Check if cancer-type has been processed already
    if algorithm:
        check_file = os.path.join(
            cancertype_dir, "{}_coefficients.tsv.gz".format(acronym)
        )
    else:
        update_name = "{}_ensemble_all_alg".format(acronym)
        check_file = os.path.join(
            cancertype_dir, "{}_coefficients.tsv.gz".format(update_name)
        )

    if check_status(check_file):
        continue

    # Now, perform all the analyses for each X matrix
    for signal in algorithm_dict.keys():
        k_dim_dict = algorithm_dict[signal]
        for k in k_dim_dict.keys():
            z_train_file = k_dim_dict[k]["train"]
            z_test_file = k_dim_dict[k]["test"]

            # Load and process data
            train_samples, x_train_df, y_train_df = align_matrices(
                x_file_or_df=z_train_file,
                y=y_df,
                add_cancertype_covariate=False,
                algorithm=None,
            )

            test_samples, x_test_df, y_test_df = align_matrices(
                x_file_or_df=z_test_file,
                y=y_df,
                add_cancertype_covariate=False,
                algorithm=None,
            )

            # Train the model
            print(
                "Training model... cancer-type: {}, "
                "ensemble algorithm: {}, signal: {}, z_dim: {}, ".format(
                    acronym, algorithm, signal, k
                )
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
            if algorithm:
                algorithm_assign = "{}_ensemble".format(algorithm)
            else:
                algorithm_assign = "all_ensemble"

            coef_df = extract_coefficients(
                cv_pipeline=cv_pipeline,
                feature_names=x_train_df.columns,
                signal=signal,
                z_dim=k,
                seed=",".join(seeds),
                algorithm=algorithm_assign,
            )

            coef_df = coef_df.assign(acronym=acronym)

            # Store all results
            train_metrics_, train_roc_df, train_pr_df = summarize_results(
                results=y_train_results,
                gene_or_cancertype=acronym,
                signal=signal,
                z_dim=k,
                seed="ensemble",
                algorithm=algorithm_assign,
                data_type="train",
            )
            test_metrics_, test_roc_df, test_pr_df = summarize_results(
                results=y_test_results,
                gene_or_cancertype=acronym,
                signal=signal,
                z_dim=k,
                seed="ensemble",
                algorithm=algorithm_assign,
                data_type="test",
            )
            cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
                results=y_cv_results,
                gene_or_cancertype=acronym,
                signal=signal,
                z_dim=k,
                seed="ensemble",
                algorithm=algorithm_assign,
                data_type="cv",
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

    if not algorithm:
        acronym = "{}_ensemble_all_alg".format(acronym)

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

    cancertype_coef_df.to_csv(
        check_file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(cancertype_dir, "{}_classify_metrics.tsv.gz".format(acronym))
    cancertype_metrics_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )
