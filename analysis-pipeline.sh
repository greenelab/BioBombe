#!/bin/bash
#
# BioBombe Analysis Pipeline
# Gregory Way, 2018
#
# This file will perform the entire BioBombe analysis pipeline. Modules 0, 1, 2, and 3
# are not performed in this script. Module 0 downloads and processes raw gene expression
# data from GTEX, TCGA, and TARGET. Module 1 performs a parameter sweep over the neural
# network models (Denoising Autoencoder and Variational Autoencoder). Module 2 actually
# performs the serial compression with increasing bottleneck layers. Module 3 builds
# the real and permuted networks used to interpret the compressed features.
#
# The first step in this pipeline is to download the precomputed data output from
# Module 2

conda activate biobombe

##############################
# Step 1: Get Data (Module 2)
##############################
cd 2.sequential-compression

python download-results-archive.py

##############################
# Step 2: Analyze Components (Module 4)
##############################
cd ../4.analyze-components

# Visualize sample reconstruction and correlation between input and reconstructed output
bash components-analysis.sh

##############################
# Step 3: Analyze Stability (Module 5)
##############################
cd ../5.analyze-stability

# Apply SVCCA to weight matrices and compare algorithm stability across dimensions
bash stability_analysis.sh

##############################
# Step 4: Analyze Weights (Module 6)
##############################
cd ../6.biobombe-projection

# Project the networks onto the compressed weight matrices to derive biological insight
python interpret-compression.py

##############################
# Step 5: Analyze Coverage of Weight Matrices (Module 7)
##############################
cd ../7.analyze-coverage

# Project the networks onto the compressed weight matrices to derive biological insight
bash run_coverage_analysis.sh

##############################
# Step 6: Perform GTEx module interpretation of increased blood correlation
##############################
cd ../8.gtex-interpret

# Determine the difference between VAE models k = 2 and k = 3 and apply features to
# external datasets
bash gtex_analysis.sh

##############################
# Step 7: Perform the TCGA classification analysis predicting cancer type and mutations
##############################
cd ../9.tcga-classify

# Train several models and visualize results
bash classify_analysis.sh
