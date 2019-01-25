"""
Gregory Way 2018
Interpret Compression
5.analyze-weights/interpret-compression.py

Track and visualize gene set activity across compressed features

Usage:

    python interpret-compression.py

Output:
Several geneset enrichment scores and many figures describing enrichment across features
"""

import os
import subprocess

dataset_collection_tuples = [
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

gmt_name_dict = {
    "gpc1": "c1.all.v6.1.entrez.gmt",
    "gpc2cpg": "c2.cgp.v6.1.entrez.gmt",
    "gpc2cpreactome": "c2.cp.reactome.v6.1.entrez.gmt",
    "gpc3mir": "c3.mir.v6.1.entrez.gmt",
    "gpc3tft": "c3.tft.v6.1.entrez.gmt",
    "gpc4cgn": "c4.cgn.v6.1.entrez.gmt",
    "gpc4cm": "c4.cm.v6.1.entrez.gmt",
    "gpc5bp": "c5.bp.v6.1.entrez.gmt",
    "gpc5cc": "c5.cc.v6.1.entrez.gmt",
    "gpc5mf": "c5.mf.v6.1.entrez.gmt",
    "gpc6": "c6.all.v6.1.entrez.gmt",
    "gpc7": "c7.all.v6.1.entrez.gmt",
    "gph": "h.all.v6.1.entrez.gmt",
    "gpxcell": "xcell_all_entrez.gmt",
}

for dataset, metaedge in dataset_collection_tuples:
    metaedge_lower = metaedge.lower()
    gene_set_dir = os.path.join("results", dataset.lower(), metaedge_lower, 'signal')
    gmt_name = gmt_name_dict[metaedge_lower]

    geneset_command = [
        "python",
        "geneset_tracking.py",
        "--dataset",
        dataset,
        "--metaedge",
        metaedge,
    ]
    visualize_command = [
        "Rscript",
        "--vanilla",
        "visualize_genesets.R",
        "--dataset",
        dataset,
        "--gmt_name",
        gmt_name,
        "--metaedge",
        metaedge,
        "--gene_set_dir",
        gene_set_dir,
        "--save_top_results",
    ]
    print(geneset_command)
    subprocess.call(geneset_command)
    print(visualize_command)
    subprocess.call(visualize_command)
