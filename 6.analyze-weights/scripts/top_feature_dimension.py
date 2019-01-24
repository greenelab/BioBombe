"""
2019 Gregory Way
6.analyze-weights/scripts/top_feature_dimension.py

Extract the top feature (by absolute value z score) for the given geneset collection
and dataset. We extract top scores by algorithm.

Usage:

    python scripts/top_feature_dimension.py

Output:
A dataframe for easy plotting of the top features per collection, dataset and algorithm
"""

import os
import pandas as pd

path = os.path.join("results", "top_features")

top_results_tuples = [
    ("TCGA", "GpH"),
    ("TCGA", "GpXCELL"),
    ("TCGA", "GpC4CM"),
    ("TCGA", "GpC2CPREACTOME"),
    ("TCGA", "GpC3TFT"),
    ("TARGET", "GpH"),
    ("TARGET", "GpXCELL"),
    ("TARGET", "GpC4CM"),
    ("GTEX", "GpXCELL"),
]

full_top_results = []
for dataset, collection in top_results_tuples:
    match_str = "{}_{}".format(dataset, collection)
    top_results_files = [x for x in os.listdir(path) if match_str in x]

    for geneset_file in top_results_files:
        geneset_df = (
            pd.read_table(os.path.join(path, geneset_file))
            .sort_values(by=["abs_z_score", "z"], ascending=False)
            .reset_index()
            .rename({"index": "absolute_rank"}, axis="columns")
            .assign(collection=collection, dataset=dataset)
        )

        full_top_results.append(geneset_df)

full_top_results_df = pd.concat(full_top_results).reset_index(drop=True)

file = os.path.join("results", "all_top_z_dimensions.tsv")
full_top_results_df.to_csv(file, sep="\t", index=False)
