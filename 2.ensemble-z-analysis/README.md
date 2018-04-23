# Ensemble Z Analysis

**Gregory Way, 2018**

Gene expression data compression reveals coordinated gene expression modules that describe important biology.
In the following analysis, we apply a series of compression algorithms to gene expression datasets.
Through each series, we alter the dimensionality (z) of the bottleneck and compare mathematical performance.
We also save the population of all models, for each algorithm, across z for downstream analyses.

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

We will evaluate the solutions across the ensemble population over all z dimensions.

1. Reconstruction Cost - Measures the binary cross entropy of input data to reconstruction
2. Training History - For neural network models (ADAGE, Tybalt), save the training progress of each model
   * For Tybalt, the KL Divergence and Reconstruction Loss are saved separately
3. Correlation of input sample to reconstructed sample - Measure how well certain samples traverse through the bottleneck.
   * Calculate Pearson and Spearman correlations
   * May reveal certain biases in sample reconstruction efficiency across algorithms

The population of weight and z matrices are also saved for alternative downstream analyses.
These analyses include automatic interpretation of compressed features with HetMech.

