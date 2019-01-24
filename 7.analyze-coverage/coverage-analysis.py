"""
Gregory Way 2018
Interpret Compression
7.analyze-coverage/coverage-analysis.py

Track gene set coverage across compressed features. We consider a geneset as "covered"
by the module if the geneset is discovered as the top feature in the BioBombe network
projection applied to the specific dataset with given internal dimension. The gene set
"coverage" refers to the percentage of individual genesets in a given resource (e.g.
xCell genesets, or GO Biological Processes) that are captured by the specific model.

Usage:

    python coverage-analysis.py

Output:
Several geneset coverage scores summarizing enrichment across features
"""

import os
import subprocess

coverage_pairs = [
    ("TCGA", "GpH"),
    ("TCGA", "GpXCELL"),
    ("TCGA", "GpC4CM"),
    ("TCGA", "GpC2CPREACTOME"),
    ("TCGA", "GpC3TFT"),
    ("TARGET", "GpH"),
    ("TARGET", "GpXCELL"),
    ("TARGET", "GpC4CM"),
    ("GTEX", "GpXCELL"),
]

for dataset, metaedge in coverage_pairs:
    metaedge_lower = metaedge.lower()
    output_dir = os.path.join("results", dataset.lower(), metaedge_lower)

    coverage_command = [
        "python",
        "get-coverage.py",
        "--dataset",
        dataset,
        "--geneset",
        metaedge_lower,
    ]
    print(coverage_command)
    subprocess.call(coverage_command)
