![logo](https://raw.githubusercontent.com/greenelab/BioBombe/master/logo.png)

# Biological Interpretation of Serially Compressed Gene Expression Features using Network Projection

**Gregory Way and Casey Greene 2018**

**University of Pennsylvania**

The repository stores data and data processing modules to enable compressed gene expression feature interpretation.

Named after the [mechanical device developed by cryptologists in World War II](https://en.wikipedia.org/wiki/Bombe) to decipher secret messages sent by [Enigma machines](https://en.wikipedia.org/wiki/Enigma_machine), BioBombe is used to decipher hidden signals in gene expression data.
Inspired by the number-crunching knobs of [Alan Turing's](https://en.wikipedia.org/wiki/Alan_Turing) device, BioBombe serially compresses gene expression input with increasing bottleneck dimensionality and deciphers the learned compressed features using biological network projection approaches.

In this repository, we serially compress three different gene expression datasets across different bottleneck dimensions (_z_) for five different algorithms.
We evaluate each algorithm and dimension on a variety of metrics.
Our goal is to construct reproducible gene expression features with unsupervised learning, and to begin interpreting what these compression features represent using network projection approaches.

Our approach is outlined below:

![overview](https://raw.githubusercontent.com/greenelab/BioBombe/master/compression-overview.png)

## BioBombe Training Implementation

To acquire our transformed data, we independently trained five different algorithms using three different datasets.
Our model implementation is described below.

![implementation](https://raw.githubusercontent.com/greenelab/BioBombe/master/biobombe-implementation.png)

## Modules

To reproduce the results of the analysis, the modules should be run in order.

| Name | Description |
| :--- | :---------- |
| [0.expression-download](0.expression-download/) | Download and process gene expression data to run through pipeline |
| [1.initial-z-sweep](1.initial-z-sweep/) | Determine a set of optimal hyperparameters for Tybalt and ADAGE models across a representative range of z dimensionality |
| [2.ensemble-z-analysis](2.ensemble-z-analysis/) | Train various algorithms to compress gene expression data across a large range of z dimensions |
| [3.build-hetnets](3.build-hetnets/) | Download, process, and integrate various curated gene sets into a single heterogeneous network |
| [4.analyze-components](4.analyze-components/) | Visualize the reconstruction and sample correlation results of the ensemble z analysis |
| [5.analyze-stability](5.analyze-stability/) | Determine how stable compression solutions are between and across algorithms, and across dimensions |
| [6.analyze-weights](6.analyze-weights/) | Apply our matrix interpretation analysis and overrepresentation analyses to assign biological knowledge to compression features |
| [7.analyze-coverage](7.analyze-coverage/) | Determine the coverage, or proportion, of enriched gene sets in compressed latent space features for all models and ensembles of models |
| [8.gtex-interpret](9.gtex-interpret/) | Interpret compressed features in the GTEX data |
| [9.tcga-classify](9.tcga-classify/) | Input compressed features from TCGA data into supervised machine learning classifiers to detect pathway aberration |

## Algorithms

See [2.ensemble-z-analysis](2.ensemble-z-analysis) for more details.

## Computational Environment

All processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.5.0
conda env create --force --file environment.yml

conda activate biobombe
```
