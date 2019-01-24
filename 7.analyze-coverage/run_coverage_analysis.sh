#!/bin/bash
#
# Gregory Way 2019
# BioBombe
# 7.analyze-coverage
#
# Extract a series of results and generate figures desribing coverage analysis

# Step 1 - Extract coverage results for select collection and data set pairs
python coverage-analysis.py

# Step 2 - Visualize coverage per individual model, ensemble model, and all models
Rscript nbconverted/scripts/visualize_coverage.r
