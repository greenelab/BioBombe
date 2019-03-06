#!/bin/bash
#
# Gregory Way 2019
# BioBombe
# 10.gene-expression-signatures/signature_analysis.sh
#
# Identify sample sex and MYCN amplification signatures from compressed gene expression

# Step 0 - Download Validation Data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.download-validation-data.ipynb

# Step 1 - Detect the top differentiating features for sex and MYCN amplification
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.separate.ipynb

# Step 2 - Apply the MYCN signature to an external validation dataset
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute 2.apply-mycn-signature.ipynb

# Step 3 - Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 4 - Visualize the results
Rscript scripts/nbconverted/3.visualize-separation.r
