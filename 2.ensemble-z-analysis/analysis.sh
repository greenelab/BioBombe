#!/bin/bash

# Fit several compression algorithms on various sets of input data with
# increasing dimensionality.

# Note that the optimal hyperparameters for each dimension are stored in the
# dataset specific configuration files that were previously determined to be
# optimal in module `1.initial-z-sweep`.

# Note this was run at the PMACS cluster at University of Pennsvylvania on a
# cluster of 8 NVIDIA GEFORCE GTX Ti GPUs.

# Set Constants
PARAM_FILE_PREFIX='config/z_parameter_sweep_'
PMACS_FILE='../config/pmacs_config.tsv'
PYTHON_PATH='python'
NUM_SEEDS=10
NUM_GENES=8000
DATASETS=( 'TCGA' 'TARGET' 'GTEX')

# Loop through each dataset and perform analysis
for dat in "${DATASETS[@]}"
do
    PARAM_FILE=$PARAM_FILE_PREFIX$dat'.tsv'
    OUT_DIR='results/'$dat'_results/'

    python scripts/z_sweep_jobs_submit.py \
          --param_config $PARAM_FILE \
          --pmacs_config $PMACS_FILE \
          --out_dir $OUT_DIR \
          --num_seeds $NUM_SEEDS \
          --subset_mad_genes $NUM_GENES \
          --components 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 45 50 60 70 80 90 100 125 150 200 \
          --dataset $dat \
          --local
done

# Loop through all datasets again with shuffled gene expression
for dat in "${DATASETS[@]}"
do
    PARAM_FILE=$PARAM_FILE_PREFIX$dat'.tsv'
    OUT_DIR='results/'$dat'_shuffled_results/'

    python scripts/z_sweep_jobs_submit.py \
          --param_config $PARAM_FILE \
          --pmacs_config $PMACS_FILE \
          --out_dir $OUT_DIR \
          --num_seeds $NUM_SEEDS \
          --subset_mad_genes $NUM_GENES \
          --components 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 45 50 60 70 80 90 100 125 150 200 \
          --dataset $dat \
          --local \
          --shuffle
done

