#!/bin/bash
#
# Extract and visualize geneset coverage metrics
# Gregory Way, 2018

####################
# Analysis 1
# Extract coverage statistics for TCGA oncogenic genesets
####################
python get-coverage.py \
        --dataset 'TCGA' \
        --geneset 'gpc6'
