"""
2018 Gregory Way
2.ensemble-z-analysis/download-biobombe-archive.py

Download BioBombe Results
All BioBombe results have been preprocessed, and have been versioned/stored on Zenodo

To acquire the data from scratch, perform:

    conda activate biobombe
    ./analysis.sh

To retrieve the precompiled data, run this script.

    python download-biobombe-archive.py

The script will download, unzip the data, and place them in appropriate directories.
"""

import os
import subprocess

os.makedirs("results", exist_ok=True)

records = ["2222463", "2222469", "2110752", "2221216", "2300616", "2386816"]
filenames = [
    "TARGET_results.zip",
    "TARGET_shuffled_results.zip",
    "TCGA_results.zip",
    "TCGA_shuffled_results.zip",
    "GTEX_results.zip",
    "GTEX_shuffled_results.zip",
]

# Downloaf files
downloaded_files = []
for record, filename in zip(records, filenames):
    url = "https://zenodo.org/record/{}/files/{}".format(record, filename)
    wget_command = [
        "wget",
        "--timestamping",
        "--no-verbose",
        "--directory-prefix",
        "results",
        url,
    ]

    download_file = os.path.join("results", filename)
    downloaded_files.append(download_file)

    subprocess.call(wget_command)

# Unzip the downloaded files
for filename in downloaded_files:
    unzip_command = ["unzip", filename, "-d", "results"]

    subprocess.call(unzip_command)

# Check the integrity of the downloaded files
check_integrity = ["md5sum", "-c", "md5sums.txt"]
subprocess.call(check_integrity)
