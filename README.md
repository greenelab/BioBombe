![logo](https://raw.githubusercontent.com/greenelab/BioBombe/master/logo.png)

# Sequential Compression of Gene Expression Data Across Latent Space Dimensions

**Gregory Way and Casey Greene 2018**

**University of Pennsylvania**

[![DOI](https://zenodo.org/badge/126377943.svg)](https://zenodo.org/badge/latestdoi/126377943)

The repository stores data and data processing modules to sequentially compress gene expression data.

Named after the [mechanical device developed by cryptologists in World War II](https://en.wikipedia.org/wiki/Bombe) to decipher secret messages sent by [Enigma machines](https://en.wikipedia.org/wiki/Enigma_machine), BioBombe is used to enhance biological signatures in gene expression data.
Inspired by the number-crunching knobs of [Alan Turing's](https://en.wikipedia.org/wiki/Alan_Turing) device, BioBombe sequentially compresses gene expression input across latent dimensionalities and deciphers the the biological signals embedded within compressed gene expression features.

In this repository, we sequentially compress three different gene expression data sets (TCGA, GTEx, and TARGET) across 28 different latent dimensions (_k_) using five different algorithms (PCA, ICA, NMF, DAE, and VAE).
We evaluate each algorithm and dimension using a variety of metrics.
Our goal is to construct reproducible gene expression signatures with unsupervised learning.

Our approach is outlined below:

![overview](https://raw.githubusercontent.com/greenelab/BioBombe/master/compression-overview.png)

## BioBombe Training Implementation

Our model implementation is described below.

![implementation](https://raw.githubusercontent.com/greenelab/BioBombe/master/biobombe-implementation.png)

## Analysis Modules

To reproduce the results and figures of the analysis, the modules should be run in order.

| Name | Description |
| :--- | :---------- |
| [0.expression-download](0.expression-download/) | Download and process gene expression data to run through pipeline |
| [1.initial-k-sweep](1.initial-k-sweep/) | Determine a set of optimal hyperparameters for Tybalt and ADAGE models across a representative range of k dimensions |
| [2.sequential-compression](2.sequential-compression/) | Train various algorithms to compress gene expression data across a large range of k dimensions |
| [3.build-hetnets](3.build-hetnets/) | Download, process, and integrate various curated gene sets into a single heterogeneous network |
| [4.analyze-components](4.analyze-components/) | Visualize the reconstruction and sample correlation results of the sequential compression analysis |
| [5.analyze-stability](5.analyze-stability/) | Determine how stable compression solutions are between and across algorithms, and across dimensions |
| [6.biobombe-projection](6.biobombe-projection/) | Apply BioBombe matrix interpretation analysis and overrepresentation analyses to assign biological knowledge to compression features |
| [7.analyze-coverage](7.analyze-coverage/) | Determine the coverage, or proportion, of enriched gene sets in compressed latent space features for all models and ensembles of models |
| [8.gtex-interpret](8.gtex-interpret/) | Interpret compressed features in the GTEX data |
| [9.tcga-classify](9.tcga-classify/) | Input compressed features from TCGA data into supervised machine learning classifiers to detect pathway aberration |
| [10.gene-expression-signatures](10.gene-expression-signatures/) | Identify gene expression signatures for sample sex in GTEx and TCGA data, and MYCN amplification in TARGET data |

## Algorithms

See [2.sequential-compression](2.sequential-compression/) for more details.

## Computational Environment

All processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.5.0
conda env create --force --file environment.yml

conda activate biobombe
```
