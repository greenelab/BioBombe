#!/bin/bash
#
# Extract and visualize geneset coverage metrics
# Gregory Way, 2018

####################
# Analysis 1
# Apply the matrix interpretation approach to GTEX data using XCELL genesets
####################
python get-coverage.py \
        --dataset 'TCGA' \
        --geneset 'gpc6'
