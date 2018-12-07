# Results for Initial Hyperparameter Sweep Considering Different Latent Dimensionality

**Gregory Way 2018**

## Latent Space Dimensionality

Compression algorithms reduce the dimensionality of input data by enforcing the number of dimensions to bottleneck.
A common problem is the decision of how many "useful" latent space features are present in data.
The solution is optimized differently for different problems or goals.
For example, when visualizing large differences between groups of data, a highly restrictive bottleneck, usually between 2 or 3 features, is required.
However, when the goal is to extract meaningful patterns in the data that may have more subtle relationships across samples, the recommendations are opaque.
As the bottleneck relaxes, the ability to explain the patterns decreases and the possibility of false positives increases.

In order to determine an optimal _range_ of compression dimensions, we propose the following.
We will first sweep over various different dimensions (results provided below) and perform several evals (to be described later).

Before sweeping over a large number of different dimensions, we perform a hyperparameter sweep of select dimensions.
In this sense, we want to minimize the effect of poor hyperparameter combinations across different dimensions contributing to performance differences.
In other words, we want to isolate the effect of changing dimensionality on the observed patterns and solutions.
Therefore, we perform a parameter sweep over several hyperparameters for the two unsupervised neural network models.
The models include a variational autoencoder (VAE; Tybalt) and a denoising autoencoder (DAE; ADAGE).

The full analysis is provided, with results visualized, in [`visualize-parameter-sweep.ipynb`](visualize-parameter-sweep.ipynb).

### Summary Figure

![summary sweep figure](https://raw.githubusercontent.com/greenelab/interpret-compression/master/1.initial-z-sweep/figures/z_dimension_sweep_summary.png)

### Datasets

We perform a hyperparamter grid search across three different training datasets.
The datasets include `TCGA`, `GTEx`, and `TARGET`.
For more details about these datasets, refer to [0.expression-download/README.md](../0.expression-download/README.md).

### Number of dimensions

Previously, we used a latent space dimensionality of `100` ([Way and Greene 2018](https://doi.org/10.1142/9789813235533_0008)).
Here, we sweep over dimensions: `5`, `25`, `50`, `75`, `100`, and `125`.

To reproduce the data for this analysis run the following command:

```bash
# From the top directory
conda activate interpret-compression

# Navigate into z-sweep directory
cd 1.initial-z-sweep
bash analysis.sh
```

## Parameter Sweep

We sweep over the following parameter combinations for `Tybalt` and `ADAGE` models in the `TCGA` dataset:

| Variable | Tybalt Values | ADAGE Values |
| :------- | :------------ | :----------- |
| Dimensionality | 5, 25, 50, 75, 100, 125 | 5, 25, 50, 75, 100, 125 |
| Learning Rate | 0.0005, 0.001, 0.0015, 0.002, 0.0025 | 0.00005, 0.00001, 0.0005, 0.001, 0.0015, 0.002 |
| Batch Size | 50, 100, 150 | 50, 100 |
| Epochs | 50, 100 | 100 |
| Kappa | 0, 0.5, 1 | |
| Sparsity | | 0, 0.000001, 0.001 |
| Noise | | 0, 0.1, 0.5 |
| Weights | | tied |

This resulted in the training of 540 Tybalt models and 648 ADAGE models.
Note that we have also tested `ADAGE` models with untied weights in the TCGA dataset (data not shown).
In this setting, performance was worse than with tied weight models.
For all downstream applications we use ADAGE models with _tied_ weights.

Our goal was to determine optimal hyperparameter combinations for both models across various bottleneck dimensionalities across each of the three datasets.

All other hyperparameter combinations used across models and datasets are located in the `config/` folder.

## Results

All results can be viewed in [`1.initial-z-sweep/visualize-parameter-sweep.ipynb`](1.initial-z-sweep/visualize-parameter-sweep.ipynb).

In that notebook, we report the results in a series of visualizations and tables for `Tybalt` and `ADAGE` across `TCGA`, `GTEx`, and `TARGET`.
Note that the notebook was converted to an R script with:

```bash
# Convert R notebook to R script for execution
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted visualize-parameter-sweep.ipynb
```

## Summary

Selection of hyperparameters across different latent space dimensionality operated as expected.
Loss was higher for lower dimensions and lower dimensions benefited the most from increased regularization and higher learning rates.
Nevertheless, we have obtained a broad set of optimal hyperparameters for use in a larger and more specific sweep of dimensionality for each of the three analyzed datasets.

## Selected Optimal Hyperparamters

The analysis allowed us to select optimal hyperparameters for each dataset and algorithm combination.
We report the results below:

### TCGA

#### VAE (Tybalt)

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate |
| :--------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 100 | 50 | 0.002 |
| 25 | 0 | 100 | 50 | 0.0015 |
| 50 | 0 | 100 | 100 | 0.0015 |
| 75 | 0 | 100 | 150  | 0.0015 |
| 100 | 0 | 100 | 150 | 0.001 |
| 125 | 0 | 100 | 150 | 0.0005 |

#### DAE (ADAGE)

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate |
| :--------- | :------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 0.0 | 100 | 50 | 0.0015 |
| 25 | 0 | 0.0 | 100 | 50 | 0.0015 |
| 50 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 75 | 0 | 0.0 | 100 | 50  | 0.0005 |
| 100 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 125 | 0 | 0.0 | 100 | 50 | 0.0005 |

### GTEx

#### VAE (Tybalt)

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate |
| :--------- | :---- | :----- | :--------- | :------------ |
| 5 | 0.5 | 100 | 100 | 0.0025 |
| 25 | 0.5 | 100 | 100 | 0.0025 |
| 50 | 0.5 | 100 | 100 | 0.002 |
| 75 | 0.5 | 100 | 50  | 0.002 |
| 100 | 0.5 | 100 | 50 | 0.0015 |
| 125 | 0.5 | 100 | 50 | 0.0015 |

#### DAE (ADAGE)

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate |
| :--------- | :------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 0.1 | 100 | 50 | 0.001 |
| 25 | 0 | 0.0 | 100 | 50 | 0.001 |
| 50 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 75 | 0 | 0.0 | 100 | 50  | 0.0005 |
| 100 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 125 | 0 | 0.0 | 100 | 50 | 0.0005 |

### TARGET

#### VAE (Tybalt)

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate |
| :--------- | :---- | :----- | :--------- | :------------ |
| 5 | 0.5 | 100 | 25 | 0.0015 |
| 25 | 0.5 | 100 | 25 | 0.0015 |
| 50 | 0.5 | 100 | 25 | 0.0015 |
| 75 | 0.5 | 100 | 25  | 0.0015 |
| 100 | 0.5 | 100 | 25 | 0.0015 |
| 125 | 0.5 | 100 | 25 | 0.0005 |

#### DAE (ADAGE)

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate |
| :--------- | :------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 25 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 50 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 75 | 0 | 0.1 | 100 | 50  | 0.0005 |
| 100 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 125 | 0 | 0.1 | 100 | 50 | 0.0005 |
