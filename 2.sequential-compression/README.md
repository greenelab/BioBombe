# Sequential Compression Analysis

**Gregory Way, 2018**

Gene expression data compression reveals coordinated gene expression modules that describe important biology.

In the following analysis, we apply five different compression algorithms to three different gene expression datasets.
We sequentially compress input data across different bottleneck dimensions (_k_).
We save the population of all models, for each algorithm, across _k_ for downstream analyses.

## Algorithms

We compress gene expression data with the following algorithms:

| Algorithm | Implementation |
| :-------- | :------------- |
| Principal Components Analysis (PCA) | [sklearn](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html) |
| Independent Components Analysis (ICA) | [sklearn](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FastICA.html) |
| Non-Negative Matrix Factorization (NMF) | [sklearn](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html) |
| Analysis of Denoising Autoencoders for Gene Expression (ADAGE) | [tybalt.models.Adage](https://github.com/greenelab/tybalt/blob/master/tybalt/models.py#L284)
| Variational Autoencoder (VAE; Tybalt) | [tybalt.models.Tybalt](https://github.com/greenelab/tybalt/blob/master/tybalt/models.py#L25)

## Evaluation Metrics

We will evaluate the solutions across the ensemble population over all _k_ dimensions.
For each of the populations, we will also track performance of training and testing sets independently.

1. Reconstruction Cost - Measures the binary cross entropy of input data to reconstruction
2. Training History - For neural network models (ADAGE, Tybalt), save the training progress of each model
   * For Tybalt, the KL Divergence and Reconstruction Loss are saved separately
3. Correlation of input sample to reconstructed sample - Measure how well certain samples traverse through the bottleneck.
   * Calculate Pearson and Spearman correlations
   * May reveal certain biases in sample reconstruction efficiency across algorithms

The population of weight and latent space matrices are saved for alternative downstream analyses.

## Download Results

This module takes a long time to run.
For convenience, we include the option to download archived and versioned results from zenodo.

To acquire these results, perform the following:

```bash
conda activate biobombe

cd 2.sequential-compression
python download-biobombe-archive.py
```

## Reproduce Analysis

To rerun the analysis from scratch, perform the following:

```bash
conda activate biobombe

# Navigate into this module folder
cd 2.sequential-compression
./analysis.sh
```
