#!/bin/bash
#
# Gregory Way 2018
# 
# Run after analysis.sh to summarize the results of the parameter sweep
# of ADAGE and Tybalt models applied across three datasets

# TCGA
python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_adage_TCGA/' \
  --output_filename 'param_sweep_adage_TCGA_full-results.tsv'

python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_tybalt_TCGA/' \
  --output_filename 'param_sweep_tybalt_TCGA_full-results.tsv'

# GTEx
python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_adage_GTEX/' \
  --output_filename 'param_sweep_adage_GTEX_full-results.tsv'

python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_tybalt_GTEX/' \
  --output_filename 'param_sweep_tybalt_GTEX_full-results.tsv'

# TARGET
python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_adage_TARGET/' \
  --output_filename 'param_sweep_adage_TARGET_full-results.tsv'

python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_tybalt_TARGET/' \
  --output_filename 'param_sweep_tybalt_TARGET_full-results.tsv'
