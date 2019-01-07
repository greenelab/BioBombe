"""
Gregory Way 2018
Interpret Compression
7.analyze-coverage/coverage-analysis.py

Track gene set coverage across compressed features

Usage:

    python coverage-analysis.py

Output:
Several geneset coverage scores summarizing enrichment across features
"""

import os
import subprocess

coverage_pairs = [('TCGA', 'GpC5BP'),
                  ('TCGA', 'GpC4CM'),
                  ('GTEX', 'GpXCELL')]

for dataset, metaedge in coverage_pairs:
    metaedge_lower = metaedge.lower()
    output_dir = os.path.join('results', dataset.lower(), metaedge_lower)

    coverage_command = ['python', 'get-coverage.py',
                        '--dataset', dataset,
                        '--geneset', metaedge_lower]
    print(coverage_command)
    subprocess.call(coverage_command)
