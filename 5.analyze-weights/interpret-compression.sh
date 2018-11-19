#!/bin/bash

####################
# Analysis 1
# Apply the matrix interpretation approach to GTEX data using XCELL genesets
####################
python geneset_tracking.py \
        --dataset 'GTEX' \
        --metaedge 'GpXCELL'
Rscript --vanilla visualize_genesets.R \
        --dataset 'GTEX' \
        --gmt_name 'xcell_all_entrez.gmt' \
        --metaedge 'GpXCELL' \
        --gene_set_dir 'results/gtex/gpxcell/signal'

# Apply to shuffled dataset to set null distribution
python geneset_tracking.py \
        --dataset 'GTEX' \
        --metaedge 'GpXCELL' \
        --shuffled
Rscript --vanilla visualize_genesets.R \
        --dataset 'GTEX' \
        --gmt_name 'xcell_all_entrez.gmt' \
        --metaedge 'GpXCELL' \
        --gene_set_dir 'results/gtex/gpxcell/shuffled' \
        --shuffled

###################
# Analysis 2
# Track Cancer Hallmarks across the TCGA dataset
###################
python geneset_tracking.py \
        --dataset 'TCGA' \
        --metaedge 'GpH'
Rscript --vanilla visualize_genesets.R \
        --dataset 'TCGA' \
        --gmt_name 'h.all.v6.1.entrez.gmt' \
        --metaedge 'GpH' \
        --gene_set_dir 'results/tcga/gph/signal'

# Apply to shuffled dataset to set null distribution
python geneset_tracking.py \
        --dataset 'TCGA' \
        --metaedge 'GpH' \
        --shuffled
Rscript --vanilla visualize_genesets.R \
        --dataset 'TCGA' \
        --gmt_name 'h.all.v6.1.entrez.gmt' \
        --gene_set_dir 'results/tcga/gph/shuffled' \
        --shuffled

###################
# Analysis 3
# Also track XCELL across TCGA
###################
python geneset_tracking.py \
        --dataset 'TCGA' \
        --metaedge 'GpXCELL'
Rscript --vanilla visualize_genesets.R \
        --dataset 'TCGA' \
        --gmt_name 'xcell_all_entrez.gmt' \
        --metaedge 'GpXCELL' \
        --gene_set_dir 'results/tcga/gpxcell/signal'

# Apply to shuffled dataset to set null distribution
python geneset_tracking.py \
        --dataset 'TCGA' \
        --metaedge 'GpXCELL' \
        --shuffled
Rscript --vanilla visualize_genesets.R \
        --dataset 'TCGA' \
        --gmt_name 'xcell_all_entrez.gmt' \
        --metaedge 'GpXCELL' \
        --gene_set_dir 'results/tcga/gpxcell/shuffled' \
        --shuffled

###################
# Analysis 4
# Track Cancer Hallmarks across the TARGET dataset
###################
python geneset_tracking.py \
        --dataset 'TARGET' \
        --metaedge 'GpH'
Rscript --vanilla visualize_genesets.R \
        --dataset 'TARGET' \
        --gmt_name 'h.all.v6.1.entrez.gmt' \
        --metaedge 'GpH' \
        --gene_set_dir 'results/target/gph/signal'

# Apply to shuffled dataset to set null distribution
python geneset_tracking.py \
        --dataset 'TARGET' \
        --metaedge 'GpH' \
        --shuffled
Rscript --vanilla visualize_genesets.R \
        --dataset 'TARGET' \
        --gmt_name 'h.all.v6.1.entrez.gmt' \
        --metaedge 'GpH' \
        --gene_set_dir 'results/target/gph/shuffled' \
        --shuffled

###################
# Analysis 5
# Also track Cell-types across TARGET
###################
python geneset_tracking.py \
        --dataset 'TARGET' \
        --metaedge 'GpXCELL'
Rscript --vanilla visualize_genesets.R \
        --dataset 'TARGET' \
        --gmt_name 'xcell_all_entrez.gmt' \
        --metaedge 'GpXCELL' \
        --gene_set_dir 'results/target/gpxcell/signal'

# Apply to shuffled dataset to set null distribution
python geneset_tracking.py \
        --dataset 'TARGET' \
        --metaedge 'GpXCELL' \
        --shuffled
Rscript --vanilla visualize_genesets.R \
        --dataset 'TARGET' \
        --gmt_name 'xcell_all_entrez.gmt' \
        --metaedge 'GpXCELL' \
        --gene_set_dir 'results/target/gpxcell/shuffled' \
        --shuffled
