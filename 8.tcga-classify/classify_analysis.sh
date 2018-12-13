#!/bin/bash
#
# Gregory Way 2018
# Interpret compression
# 8.tcga-classify/classify_analysis.sh
#
# Perform a series of classification tasks on compressed gene expression data

# Analysis 1 - Calculate predictions on the top 50 most mutated genes in cancer
python classify-top-mutations.py

# Analysis 2 - Calculate predictions on cancer-type membership
python classify-cancer-types.py

# Analysis 3 - Calculate predictions with uncompressed RNAseq data
python classify-with-raw-expression.py
