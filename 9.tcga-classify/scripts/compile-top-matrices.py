"""
2018 Gregory Way
9.tcga-classify/scripts/classify-top-mutations.py

Predict if specific genes are mutated across TCGA tumors based on compressed expression
features by various algorithms and dimensionalities.

Usage:

    python scripts/classify-top-mutations.py

Output:
Gene specific DataFrames storing ROC, precision-recall, and classifier coefficients
for every compression model trained in their ability to predict mutations. The genes
used are the top 50 most mutated genes in TCGA PanCanAtlas. A gene was considered
mutated if a non-silent mutation was observed by the MC3 mutation calling effort. An
additional metrics file that stores all gene AUROC and AUPR is also saved.
"""

import os
import pandas as pd

from tcga_util import build_good_feature_matrix

# Get top 200 features per algorithm
n = 200
genes = ["TP53", "PTEN", "PIK3CA", "KRAS", "TTN"]
algorithms = ["pca", "ica", "nmf", "dae", "vae", "all"]
base_path = os.path.join("results", "top_features_matrices")

for gene in genes:
    # Set file name
    file = os.path.join(
        "results", "mutation", gene, "{}_coefficients.tsv.gz".format(gene)
    )

    # Load file
    coef_df = pd.read_table(file)

    # Process file
    split_df = pd.DataFrame(
        coef_df.feature.str.split("_").values.tolist(),
        columns=["tissue", "feature_num"],
    )
    coef_df = (
        pd.concat([coef_df, split_df], axis="columns")
        .reset_index(drop=True)
        .drop("tissue", axis="columns")
    )

    coef_df = coef_df[~coef_df.feature_num.isna()].query("feature_num != 'mut'")

    # Loop through all algorithms and extract specific features
    for algorithm in algorithms:

        # Inform status
        print("now processing: {} {}".format(gene, algorithm))

        # If the algorithm is not all, subset by specific algorithm
        if algorithm != "all":
            coef_subset_df = coef_df.query("algorithm == @algorithm")
        else:
            coef_subset_df = coef_df.copy()

        coef_subset_df = coef_subset_df.sort_values(by="abs", ascending=False).iloc[
            range(0, n), :
        ]

        # Loop through each specific feature that was previously identified as important
        all_test_features_df, all_train_features_df = build_good_feature_matrix(
            coef_df=coef_subset_df
        )

        base_file = os.path.join(
            base_path, "top_model_algorithm_{}_gene_{}".format(algorithm, gene)
        )
        test_file = "{}_test.tsv.gz".format(base_file)
        train_file = "{}_train.tsv.gz".format(base_file)

        all_test_features_df.to_csv(test_file, sep="\t", compression="gzip")
        all_train_features_df.to_csv(train_file, sep="\t", compression="gzip")
