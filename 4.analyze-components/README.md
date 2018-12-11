# Analyzing Compressed Gene Expression Sample Outputs

**Gregory Way 2018**

Samples within each dataset are compressed with dimensionality reduction algorithms.
When these samples pass through the bottleneck, they are variably activated by each latent space feature.
The activations can represent biological or technical patterns present in a specific quantity for the given sample and feature.

The samples are also reconstructed back into the original input gene expression dimensions.
In this module, we explore the ability of compression models to reconstruct and capture signals present in input samples.

## Reproduce Analysis

The following results were obtained by running the following

```bash
conda activate biobombe

# Navigate into this module folder
cd 4.analyze-components
./analysis.sh
```

## Reconstruction

We compute the reconstruction cost using binary cross entropy between input and reconstructed output.
The analysis is provided in [1.visualize-reconstruction.ipynb](1.visualize-reconstruction.ipynb).

The results are stored in the [results](results/) folder and the figures are stored in dataset specific folders in [figures](figures/).

![reconstruction](https://raw.githubusercontent.com/greenelab/BioBombe/master/4.analyze-components/figures/reconstruction_summary.png)

## Correlation

We also track the Pearson and Spearman correlations of input samples to reconstructed output samples.
This analysis is provided in [2.visualize-sample-correlation.ipynb](2.visualize-sample-correlation.ipynb).

The results are stored as gzipped files in the [results](results/) folder and the figures are stored in dataset specific folders in [figures](figures/).

![correlation](https://raw.githubusercontent.com/greenelab/BioBombe/master/4.analyze-components/figures/correlation_summary.png)
