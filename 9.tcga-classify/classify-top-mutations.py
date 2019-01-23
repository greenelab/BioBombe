"""
2018 Gregory Way
9.tcga-classify/classify-top-mutations.py

Use top compression features to predict mutation status. Compare top 200 to top 1 and
to random 200 features to performance with 200 features alone.

Usage:

    python classify-top-mutations.py

Output:
Gene specific DataFrames storing ROC, precision-recall, and classifier coefficients
for every compression model trained in their ability to predict mutations. The genes
used are the top 50 most mutated genes in TCGA PanCanAtlas. A gene was considered
mutated if a non-silent mutation was observed by the MC3 mutation calling effort. An
additional metrics file that stores all gene AUROC and AUPR is also saved.
"""

import os
import numpy as np
import pandas as pd

from scripts.tcga_util import (
    get_threshold_metrics,
    summarize_results,
    extract_coefficients,
    align_matrices,
    process_y_matrix,
    train_model,
    build_top_feature_dictionary,
)

np.random.seed(123)

# Load constants
filter_prop = 0.05
filter_count = 15
folds = 5
max_iter = 100
alphas = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3]
l1_ratios = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4]
genes = ["TP53", "PTEN", "PIK3CA", "KRAS", "TTN"]
algorithms = ["pca", "ica", "nmf", "dae", "vae", "all"]
num_features = [1, 200]

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

# Obtain a dictionary of file directories for loading each feature matrix (X)
x_matrix_dict = build_top_feature_dictionary(
    algorithms=algorithms, genes=genes, num_features=num_features
)

for gene_name in genes:

    # Create list to store gene specific results
    gene_auc_list = []
    gene_aupr_list = []
    gene_coef_list = []
    gene_metrics_list = []

    # Create directory for the gene
    gene_dir = os.path.join("results", "top_feature_classification", gene_name)
    os.makedirs(gene_dir, exist_ok=True)

    # Get gene info to build y matrix
    gene_info = genes_df.query("gene == @gene_name")
    classification = gene_info.classification.values[0]

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

    for algorithm in algorithms:
        for n in num_features:
            for randomized in (False, True):
                if randomized and n == 200:
                    train_df = x_matrix_dict[algorithm][gene_name][n]["randomized"][
                        "train"
                    ]
                    test_df = x_matrix_dict[algorithm][gene_name][n]["randomized"][
                        "test"
                    ]
                    signal_is = "randomized"
                elif randomized and n == 1:
                    continue
                else:
                    train_df = x_matrix_dict[algorithm][gene_name][n]["train"]
                    test_df = x_matrix_dict[algorithm][gene_name][n]["test"]
                    signal_is = "real"

                # Load and process data
                train_samples, x_train_df, y_train_df = align_matrices(
                    x_file_or_df=train_df, y=y_df
                )

                test_samples, x_test_df, y_test_df = align_matrices(
                    x_file_or_df=test_df, y=y_df
                )

                # Train the model
                print(
                    "Training model... gene: {}, "
                    "algorithm: {}, feature number: {}, randomized: {}".format(
                        gene_name, algorithm, n, randomized
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
                y_cv_results = get_threshold_metrics(
                    y_train_df.status, y_cv_df, drop=False
                )

                # Get coefficients
                coef_df = extract_coefficients(
                    cv_pipeline=cv_pipeline,
                    feature_names=x_train_df.columns,
                    signal=signal_is,
                    z_dim="top_features: {}".format(n),
                    seed="any",
                    algorithm=algorithm,
                )

                coef_df = coef_df.assign(gene=gene_name)

                # Store all results
                train_metrics_, train_roc_df, train_pr_df = summarize_results(
                    y_train_results,
                    gene_name,
                    signal_is,
                    "top_features: {}".format(n),
                    "any",
                    algorithm,
                    "train",
                )
                test_metrics_, test_roc_df, test_pr_df = summarize_results(
                    y_test_results,
                    gene_name,
                    signal_is,
                    "top_features: {}".format(n),
                    "any",
                    algorithm,
                    "test",
                )
                cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
                    y_cv_results,
                    gene_name,
                    signal_is,
                    "top_features: {}".format(n),
                    "any",
                    algorithm,
                    "cv",
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

    file = os.path.join(gene_dir, "{}_auc_threshold_metrics.tsv.gz".format(gene_name))
    gene_auc_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(gene_dir, "{}_aupr_threshold_metrics.tsv.gz".format(gene_name))
    gene_aupr_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(gene_dir, "{}_coefficients.tsv.gz".format(gene_name))
    gene_coef_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )

    file = os.path.join(gene_dir, "{}_classify_metrics.tsv.gz".format(gene_name))
    gene_metrics_df.to_csv(
        file, sep="\t", index=False, compression="gzip", float_format="%.5g"
    )
