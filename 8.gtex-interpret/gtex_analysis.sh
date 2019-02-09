#!/bin/bash
#
# Gregory Way 2019
# BioBombe
# 8.gtex-interpret/gtex_analysis.sh
#
# Perform the full pipeline for the GTEx analysis and validation

# Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 0 - Apply BioBombe Network Projection to the VAE features of interest
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.process-gtex-blood-vae-feature.ipynb

# Step 1 - Download and process publicly available gene expression datasets
#      A - Neutrophils: Two Leukemia Cell Line Datasets (GSE103706)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1A.download-neutrophil-data.ipynb

#      B - Monocytes: Several isolated populations in stages of hematopoiesis (GSE24759)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1B.download-hematopoietic-data.ipynb

# Step 2 - Apply the signatures identified in step 0 to the public datasets
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.apply-signatures.ipynb

# Step 3 - Apply compressed feature signatures to neutrophil and monocyte data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.detect-signature-separation-ttest.ipynb

# Step 4 - Generate the final figure visualizing the applied signatures
Rscript --vanilla scripts/nbconverted/4.visualize-signatures.r
