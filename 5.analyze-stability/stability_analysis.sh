#!/bin/bash
#
# Gregory Way 2018
# Interpret compression
# 5.analyze-stability/stability_analysis.sh
#
# Apply SVCCA to the compressed representations to determine model stability and
# feature representation variety across iterations, algorithms, and dimensions

# Analysis 1 - Calculate within z stability across datasets and dimensions
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute 1.stability-within-z.ipynb

# Analysis 2 - Calculate across z stability within datasets
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute 2.stability-across-z.ipynb

# Analysis 3 - Visualize the results
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=100000 \
        --execute 3.stability-visualize.ipynb
