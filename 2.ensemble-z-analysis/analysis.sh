#!/bin/bash

param_file='config/full_z_sweep_params.tsv'
pmacs_config='../config/pmacs_config.tsv'
out_dir='results'
python_path='python'
num_seeds=10

# Perform ensemble z sweep over several latent space dimensionality estimates
python scripts/z_sweep_jobs_submit.py --local --param_config $param_file --pmacs_config $pmacs_config \
        --out_dir $out_dir --python_path $python_path --num_seeds $num_seeds \
        --components 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 45 50 60 70 80 90 100 125 150 200
