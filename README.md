# Interpreting Compressed Gene Expression Features 

**Gregory Way 2018**

**University of Pennsylvania**

The repository stores data and data processing modules to enable compressed gene expression feature interpretation.

## Modules 

To reproduce the results of the analysis, the modules should be run in order.

| Name | Description |
| :--- | :---------- |
| [0.expression-download](0.expression-download/) | Download and process gene expression data to run through pipeline |
| [1.initial-z-sweep](1.initial-z-sweep/) | Determine a set of optimal hyperparameters for Tybalt and ADAGE models across a representative range of z dimensionality |
| [2.ensemble-z-analysis](2.ensemble-z-analysis/) | Train various algorithms to compress gene expression data across a large range of z dimensions |
| [3.build-hetnets](3.build-hetnets/) | Download, process, and integrate various curated gene sets into a single heterogeneous network |
| [4.analyze-components](4.analyze-components/) | Visualize the reconstruction and sample correlation results of the ensemble z analysis |
| [5.analyze-weights](5.analyze-weights/) | Apply our matrix interpretation analysis and GSEA to assign biological knowledge to compression features |

## Algorithms

See [2.ensemble-z-analysis](2.ensemble-z-analysis) for more details.

## Computational Environment

All processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.5.0
conda env create --force --file environment.yml

conda activate interpret-compression
```
