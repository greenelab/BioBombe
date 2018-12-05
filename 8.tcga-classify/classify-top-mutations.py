"""
2018 Gregory Way
7.tcga-classify/classify-top-mutations.py

Predict gene mutation status in TCGA tumors based on compressed expression features by
various algorithms and dimensionalities.

Usage:

    python classify-top-mutations.py

Output:
Gene specific DataFrames storing ROC, precision-recall, and classifier coefficients
for every compression model trained in their ability to predict mutations. The genes
used are the top 50 most mutated genes in TCGA PanCanAtlas. An additional metrics file
that stores all gene AUROC and AUPR is also saved.
"""

import os
import glob
import numpy as np
import pandas as pd

from scripts.tcga_util import (
    get_threshold_metrics,
    summarize_results,
    extract_coefficients,
    align_matrices,
    process_y_matrix,
    train_model,
)

np.random.seed(123)

# Load constants
filter_prop = 0.05
filter_count = 15
folds = 5
algorithms = ["pca", "ica", "nmf", "dae", "vae"]
signals = ["signal", "shuffled"]
alphas = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3]
l1_ratios = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4]

# Load genes
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

# Track total metrics for each gene in one file
full_metrics_list = []

# Obtain a dictionary of file directories for loading each feature matrix (X)
z_matrix_dict = {}
for signal in signals:
    z_matrix_dict[signal] = {}

    if signal == "signal":
        results_dir = "TCGA_results"
    else:
        results_dir = "TCGA_shuffled_results"

    matrix_dir = os.path.join(
        "..", "2.ensemble-z-analysis", "results", results_dir, "ensemble_z_matrices"
    )

    for comp_dir in os.listdir(matrix_dir):
        matrix_comp_dir = os.path.join(matrix_dir, comp_dir)
        z_dim = comp_dir.split("_")[2]
        z_matrix_dict[signal][z_dim] = {}

        for z_file in glob.glob("{}/*_z_*".format(matrix_comp_dir)):
            seed = os.path.basename(z_file).split("_")[1]

            if seed not in z_matrix_dict[signal][z_dim].keys():
                z_matrix_dict[signal][z_dim][seed] = {}

            if "_test_" in z_file:
                z_matrix_dict[signal][z_dim][seed]["test"] = z_file
            else:
                z_matrix_dict[signal][z_dim][seed]["train"] = z_file

for gene_idx, gene_series in genes_df.iterrows():

    gene_name = gene_series.gene
    classification = gene_series.classification

    # Create list to store gene specific results
    gene_auc_list = []
    gene_aupr_list = []
    gene_coef_list = []

    # Create directory for the gene
    gene_dir = os.path.join("results", gene_name)
    os.makedirs(gene_dir, exist_ok=True)

    # Process the y matrix for the given gene or pathway
    y_mutation_df = mutation_df.loc[:, gene_name]

    # Include copy number gains for oncogenes and copy number loss for tumor suppressors
    if classification == "Oncogene":
        y_copy_number_df = copy_gain_df.loc[:, gene_name]
    elif classification == "TSG":
        y_copy_number_df = copy_loss_df.loc[:, gene_name]

    y_df = process_y_matrix(
        y_copy=y_copy_number_df,
        y_mutation=y_mutation_df,
        gene=gene_name,
        sample_freeze=sample_freeze_df,
        mutation_burden=mut_burden_df,
        filter_count=filter_count,
        filter_prop=filter_prop,
        output_directory=gene_dir,
        hyper_filter=5,
    )

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
                        x_file=z_train_file, y=y_df, algorithm=alg
                    )

                    test_samples, x_test_df, y_test_df = align_matrices(
                        x_file=z_test_file, y=y_df, algorithm=alg
                    )

                    # Train the model
                    print(
                        "Training model... gene: {}, "
                        "algorithm: {}, signal: {}, z_dim: {}, "
                        "seed: {}".format(gene_name, alg, signal, z_dim, seed)
                    )

                    # Fit the model
                    cv_pipeline, y_pred_train_df, y_pred_test_df, y_cv_df = train_model(
                        x_train=x_train_df,
                        x_test=x_test_df,
                        y_train=y_train_df,
                        alphas=alphas,
                        l1_ratios=l1_ratios,
                        n_folds=folds,
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
                        gene=gene_name,
                        signal=signal,
                        z_dim=z_dim,
                        seed=seed,
                        algorithm=alg,
                    )

                    # Store all results
                    train_metrics_, train_roc_df, train_pr_df = summarize_results(
                        y_train_results, gene_name, signal, z_dim, seed, alg, "train"
                    )
                    test_metrics_, test_roc_df, test_pr_df = summarize_results(
                        y_test_results, gene_name, signal, z_dim, seed, alg, "test"
                    )
                    cv_metrics_, cv_roc_df, cv_pr_df = summarize_results(
                        y_cv_results, gene_name, signal, z_dim, seed, alg, "cv"
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

                    gene_auc_df = pd.concat([train_roc_df, test_roc_df, cv_roc_df])
                    gene_auc_list.append(gene_auc_df)

                    gene_aupr_df = pd.concat([train_pr_df, test_pr_df, cv_pr_df])
                    gene_aupr_list.append(gene_aupr_df)

                    gene_coef_list.append(coef_df)

    gene_auc_df = pd.concat(gene_auc_list)
    gene_aupr_df = pd.concat(gene_aupr_list)
    gene_coef_df = pd.concat(gene_coef_list)

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

# Now, compile all results and write to file
final_metrics_df = pd.concat(full_metrics_list)

file = os.path.join("results", "all_gene_classify_metrics.tsv.gz")
final_metrics_df.to_csv(
    file, sep="\t", index=False, compression="gzip", float_format="%.5g"
)
