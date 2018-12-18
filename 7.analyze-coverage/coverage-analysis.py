"""
Gregory Way 2018
Interpret Compression
7.analyze-coverage/coverage-analysis.py

Track and visualize gene set activity across compressed features

Usage:

    python coverage-analysis.py

Output:
Several geneset enrichment scores and many figures describing enrichment across features
"""

import os
import subprocess

datasets = ['GTEX', 'TCGA', 'TARGET']
metaedges = ['GpXCELL', 'GpH', 'GpC1', 'GpC2CPG', 'GpC2PREACTOME', 'GpC3MIR', 'GpC3TFT',
             'GpC4CGN', 'GpC4CM', 'GpC5BP', 'GpC5CC', 'GpC5MF', 'GpC6', 'GpC7']

for dataset in datasets:
    for metaedge in metaedges:
        metaedge_lower = metaedge.lower()
        output_dir = os.path.join('results', dataset.lower(), metaedge_lower)

        coverage_command = ['python', 'get-coverage.py',
                            '--dataset', dataset,
                            '--geneset', metaedge_lower]
        visualize_command = ['Rscript', '--vanilla', 'visualize_coverage.R',
                             '--dataset', dataset,
                             '--geneset', metaedge_lower,
                             '--output_dir', output_dir]
        print(coverage_command)
        subprocess.call(coverage_command)
        print(visualize_command)
        subprocess.call(visualize_command)
