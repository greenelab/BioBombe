#!/bin/bash
#
# Gregory Way 2018
# Interpret compression
# 9.tcga-classify/classify_analysis.sh
#
# Perform a series of classification tasks on compressed gene expression data

# Analysis 1 - Calculate predictions on the top 50 most mutated genes in cancer
python classify-top-mutations.py

# Analysis 2 - Calculate predictions on cancer-type membership
python classify-cancer-types.py

# Analysis 3 - Calculate predictions with uncompressed RNAseq data
python classify-with-raw-expression.py

# Analysis 4 - Build the top feature matrices and perform classifications with features
bash build_top_feature_matrices.py
python classify-using-top-compression-features.py
