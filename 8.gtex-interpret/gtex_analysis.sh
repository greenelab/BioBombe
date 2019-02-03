#!/bin/bash
#
# Gregory Way 2019
# BioBombe
# 8.gtex-interpret/gtex_analysis.sh
#
# Perform the full pipeline for the GTEx analysis and validation

# Step 0 - Apply BioBombe Network Projection to the VAE features of interest
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.process-gtex-blood-vae-feature.ipynb

# Step 1 - Visualize BioBombe application (this will generate a supplementary figure)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.visualize-gtex-blood-interpretation.ipynb

# Step 2 - Download and process publicly available gene expression datasets
#      A - Neutrophils: Two Leukemia Cell Line Datasets (GSE103706)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2A.download-neutrophil-data.ipynb

#      B - Monocytes: Several isolated populations in stages of hematopoiesis (GSE24759)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2A.download-hematopoietic-data.ipynb

# Step 3 - Apply the signatures identified in step 0 to the public datasets
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.apply-signatures.ipynb

# Step 4 - Generate the final figure visualizing the applied signatures
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.visualize-signatures.ipynb

# Step 5 - Apply compressed feature signatures to neutrophil and monocyte data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 5.detect-signature-separation.ipynb

# Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb
