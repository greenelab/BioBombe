#!/bin/bash

# Analysis 1
# Apply the matrix interpretation approach to GTEX data using XCELL genesets
python geneset_tracking.py --dataset 'GTEX' --metaedge 'GpXCELL'
Rscript --vanilla visualize_genesets.R --dataset 'GTEX' --gmt_name 'xcell_all_entrez.gmt'

# Analysis 2
# Track Cancer Hallmarks across the TCGA dataset
python geneset_tracking.py --dataset 'TCGA' --metaedge 'GpH'
Rscript --vanilla visualize_genesets.R --dataset 'TCGA' --gmt_name 'h.all.v6.1.entrez.gmt'
