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

### TARGET

Therapeutically Applicable Research to Generate Effective Treatments (TARGET) has profiled over 700 cases of pediatric cancer from 7 different cancer-types.
We access the TARGET data using [UCSC Xena](https://xenabrowser.net/datapages/?dataset=target_RSEM_isoform_fpkm&host=https%3A%2F%2Ftoil.xenahubs.net).
We use the RSEM FPKM RNAseq processed data.

### GTEx

The Genotype-Tissue Expression ([GTEx](https://www.gtexportal.org/home/documentationPage)) project measured gene expression on over 11,000 healthy samples.
These samples represent several different tissue types.
We use [version 7](https://www.gtexportal.org/home/datasets) of GTEx RNAseq data (TPM normalized).

