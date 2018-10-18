"""
2018 Gregory Way
5.analyze-components/gsea.py

Submit jobs that perform GSEA on all compiled weight matrices

Usage:

    Run in command line:

            python gsea.py

     with optional arguments:
       --pmacs_config       filepath pointing to PMACS configuration file
       --python_path        absolute path of PMACS python in select environment
                              default: '~/.conda/envs/tybalt-gpu/bin/python'
       --local              if provided, sweep will be run locally instead

Output:

    Will submit jobs for each ensemble of compressed weight matrix (directory)
"""

import os
import sys
import argparse
import pandas as pd

sys.path.insert(0, '../scripts/util')
from bsub_helper import bsub_help

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--pmacs_config',
                    default='../config/pmacs_config.tsv',
                    help='location of the configuration file for PMACS')
parser.add_argument('-p', '--python_path', help='absolute path of python',
                    default='python')
parser.add_argument('-l', '--local', action='store_true',
                    help='decision to run models locally instead of on PMACS')
args = parser.parse_args()

pmacs_config_file = args.pmacs_config
python_path = args.python_path
local = args.local

# Load data
config_df = pd.read_table(pmacs_config_file, index_col=0)

# Retrieve PMACS configuration
queue = config_df.loc['queue']['assign']
num_gpus = config_df.loc['num_gpus']['assign']
num_gpus_shared = config_df.loc['num_gpus_shared']['assign']
walltime = config_df.loc['walltime']['assign']

# Load constants to be used in downstream analysis
algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']
distrib_methods = ['full', 'full_squared', 'pos_neg', 'pos_neg_high_weight']
num_perm = '20'

# First step is to compile all files and directories
base_file = os.path.join('..', '2.ensemble-z-analysis', 'results')
datasets = ['TCGA', 'TARGET', 'GTEX']
processing = ['results', 'shuffled_results']

# Set default parameter combination
conda = ['conda', 'activate', 'interpret-compression', '&&']

# Generate a dictionary that stores the locations of the directories and other
# pertinent information about the model of interest
all_file_dict = {}
for dataset in datasets:
    all_file_dict[dataset] = {}
    for signal in processing:
        all_file_dict[dataset][signal] = {}
        full_dir = os.path.join(base_file, '{}_{}'.format(dataset, signal),
                                'ensemble_z_matrices')
        for comp_file in os.listdir(full_dir):
            z_dir = os.path.join(full_dir, comp_file)
            z_dim = comp_file.split('_')[2]
            all_file_dict[dataset][signal][z_dim] = z_dir

# Now, loop through this dictionary and generate an iterable list
all_commands = []
for dataset in all_file_dict.keys():
    dataset_file_dict = all_file_dict[dataset]
    for signal in dataset_file_dict.keys():
        signal_file_dict = dataset_file_dict[signal]
        for z in signal_file_dict.keys():
            input_z_weight_dir = signal_file_dict[z]

            command = [python_path, 'scripts/submit_gsea.py',
                       '--input_weight_dir', input_z_weight_dir,
                       '--z_dim', z,
                       '--dataset_name', dataset,
                       '--num_perm', num_perm,
                       '--algorithms', ' '.join(algorithms),
                       '--distrib_methods', ' '.join(distrib_methods)]

            if signal != 'results':
                command += ['--shuffled']

            if not local:
                command = conda + command

            all_commands.append(command)

# Submit the jobs
if __name__ == "__main__":
    for command in all_commands:
        print(command)
        b = bsub_help(command=command,
                      queue=queue,
                      num_gpus=num_gpus,
                      num_gpus_shared=num_gpus_shared,
                      walltime=walltime,
                      local=local)
        b.submit_command()
