"""
2018 Gregory Way
9.tcga-classify/scripts/compile_top_feature_matrics.py

Scan through all prediction models and select the features with the top coefficients.
Then, find those features and construct representative matrices to use in prediction
tasks for the same mutations.

Usage:

    python scripts/compile_top_feature_matrices.py

    with the following optional arguments:

        --num_top_features <how many features to select>
        --random_features <if included, will select random features>

Output:
Several X matrices to be used in downstream classification tasks
"""

import os
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "-n", "--num_top_features", help="how many top features to use", default=200
)
parser.add_argument(
    "-r", "--random_features", action="store_true", help="select random features"
)
args = parser.parse_args()

# Load constants
num_top_features = int(args.num_top_features)
random_features = args.random_features
genes = ["TP53", "PTEN", "PIK3CA", "KRAS", "TTN"]
algorithms = ["pca", "ica", "nmf", "dae", "vae", "all"]
base_path = os.path.join("results", "top_feature_matrices")


def get_feature(z_dim, algorithm, seed, feature):
    """
    Load z matrix and extract specific feature scores
    """
    base_file = os.path.join(
        "..",
        "2.ensemble-z-analysis",
        "results",
        "TCGA_results",
        "ensemble_z_matrices",
        "tcga_components_{}".format(z_dim),
    )
    train_file = os.path.join(base_file, "model_{}_z_matrix.tsv.gz".format(seed))
    test_file = os.path.join(base_file, "model_{}_z_test_matrix.tsv.gz".format(seed))

    test_df = pd.read_table(test_file, index_col=0).loc[:, feature].sort_index()
    train_df = pd.read_table(train_file, index_col=0).loc[:, feature].sort_index()

    return test_df, train_df

np.random.seed(123)

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
        print(
            "Now processing: {} {} for top {} features. Random features? {}".format(
                gene, algorithm, num_top_features, random_features
            )
        )

        # If the algorithm is not all, subset by specific algorithm
        if algorithm != "all":
            coef_subset_df = coef_df.query("algorithm == @algorithm")
        else:
            coef_subset_df = coef_df.copy()

        if random_features:
            feature_idx = np.random.randint(
                0, high=coef_subset_df.shape[0], size=num_top_features
            )
        else:
            feature_idx = range(0, num_top_features)

        # Subset coefficient dataframe
        coef_subset_df = coef_subset_df.sort_values(by="abs", ascending=False).iloc[
            feature_idx, :
        ]

        # Loop through each specific feature that was previously identified as important
        all_test_features = []
        all_train_features = []
        for feature_idx, feature_row in coef_subset_df.iterrows():

            # Extract identifier information for the specific feature
            z_dim = feature_row.z_dim
            seed = feature_row.seed
            feature = feature_row.feature
            alg = feature_row.algorithm

            # Load the specific z matrix (training and testing) and pull out feature
            test_feature_df, train_feature_df = get_feature(z_dim, alg, seed, feature)

            # Append to growing list
            all_test_features.append(test_feature_df)
            all_train_features.append(train_feature_df)

        # Once the list is compiled, turn into DataFrame and write to file
        base_file = os.path.join(
            base_path,
            "top_model_algorithm_{}_gene_{}_numtopfeatures_{}".format(
                algorithm, gene, num_top_features
            ),
        )

        if random_features:
            base_file = "{}_randomized".format(base_file)

        all_test_features_df = pd.concat(all_test_features, axis="columns")
        all_train_features_df = pd.concat(all_train_features, axis="columns")

        test_file = "{}_test.tsv.gz".format(base_file)
        train_file = "{}_train.tsv.gz".format(base_file)

        all_test_features_df.to_csv(test_file, sep="\t", compression="gzip")
        all_train_features_df.to_csv(train_file, sep="\t", compression="gzip")
