# Downloading and Processing Gene Expression Data

**Gregory Way, 2018**

This module stores scripts to download and process gene expression data.
The processed files are tracked in this repository, so there is no need to rerun the downloading scripts.
All processed files will be used for either training or evaluation.

## RNAeq Data

### The Cancer Genome Atlas PanCanAtlas

This data was generated as a multicenter effort to profile over 10,000 tumors from 33 different cancer-types.
The list of data used as part of this effort is listed in the [Genomic Data Commons of The National Cancer Institute](https://gdc.cancer.gov/about-data/publications/pancanatlas).
We download, process, and train compression models using the `RNA (Final)` data listed there.

