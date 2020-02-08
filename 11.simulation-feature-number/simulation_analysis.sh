#!/bin/bash
#
# Gregory Way 2019
# BioBombe
# 11.simulation-feature-number/simulation_analysis.sh
#
# Generate simulated data and then apply BioBombe approach to determine
# if the number of the compression feature is assocated with signal importance

# First, convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 0 - Generate simulated data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.generate-simulated-data.ipynb

# Step 1 - Perform BioBombe approach
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.compression-simulation.ipynb

# Step 2 - Visualize simulation results
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=100000 \
        --execute 2.visualize-feature-importance.ipynb
