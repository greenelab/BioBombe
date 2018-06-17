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
Therefore, we perform a parameter sweep over several hyperparameters for both Tybalt and ADAGE models below.

The analysis is provided, with results visualized, in [`visualize-parameter-sweep.ipynb`](visualize-parameter-sweep.ipynb).

### Datasets

We perform a hyperparamter grid search across three different training datasets.
The datasets include `TCGA`, `GTEx`, and `TARGET`.
For more details about these datasets, refer to [0.expression-download/README.md](../0.expression-download/README.md).

### Number of dimensions

Previously, we used latent space dimensionality of `100` ([Way and Greene 2018](https://doi.org/10.1142/9789813235533_0008)).
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

We report the results in a series of visualizations and tables for `Tybalt` and `ADAGE` across `TCGA`, `GTEx`, and `TARGET` separately below.

First, in order to compile the results of the parameter sweep, run the following commands:

```bash
# Compile parameter sweep results
bash scripts/summarize_sweep.sh

# Convert R notebook to R script for execution
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts visualize-parameter-sweep.ipynb

# Visualize the results of the sweep for all models for all datasets
Rscript --vanilla scripts/visualize-parameter-sweep.r
```

### TCGA

#### Tybalt

Tybalt models showed variable performance, but they were generally stable across all hyperparameter combinations (**Figure 1 (TCGA, Tybalt)**).

![](figures/tcga_results/z_parameter_final_loss_tybalt_TCGA.png?raw=true)

**Figure 1 (TCGA, Tybalt).** The loss of validation sets at the end of training for all 540 TCGA Tybalt models.

Model performance (measured by observed validation loss) improved after increasing the capacity of the models from 5 to 125 dimensions.
However, the performance began to level off after around 50 dimensions.
All other hyperparameters had minor effects and varied with dimensionality, including `Kappa`.

Selecting a constant learning rate, batch size, and epochs for all models, we stratify the entire training process across all tested dimensionalities (**Figure 2 (TCGA, Tybalt)**).

![](figures/tcga_results/z_parameter_tybalt_training_TCGA.png?raw=true)

**Figure 2 (TCGA, Tybalt).** Validation loss across all training epochs for a fixed combination of hyperparameters.
The `learning rate` was set to 0.0005, `batch size` 50, and `epochs` 100.
We also show performance differences across `kappa`.

This analysis allowed us to select optimal models based on tested hyperparameters.
For `TCGA` `Tybalt` models, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate |
| :--------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 100 | 50 | 0.002 |
| 25 | 0 | 100 | 50 | 0.0015 |
| 50 | 0 | 100 | 100 | 0.0015 |
| 75 | 0 | 100 | 150  | 0.0015 |
| 100 | 0 | 100 | 150 | 0.001 |
| 125 | 0 | 100 | 150 | 0.0005 |

Training of models with optimal hyperparameters are shown in **Figure 3 (TCGA, Tybalt)**.

![](figures/tcga_results/z_parameter_best_model_tybalt_TCGA.png?raw=true)

**Figure 3 (TCGA, Tybalt).** Training optimal Tybalt models across different latent space dimensions.

### ADAGE

By constraining the compression and decompression networks to contain the same weights (tied weights), ADAGE models had variable performance across models.
ADAGE models failed to converge with low learning rates (**Figure 4 (TCGA, ADAGE)**).

![](figures/tcga_results/z_parameter_final_loss_adage_TCGA.png?raw=true)

**Figure 4 (TCGA, ADAGE).** The loss of validation sets at the end of training for 648 tied weight TCGA ADAGE models.

After removing models that did not converge we see a clearer picture (**Figure 5 (TCGA, ADAGE)**).

![](figures/tcga_results/z_parameter_final_loss_remove_converge_adage_TCGA.png?raw=true)

**Figure 5 (TCGA, ADAGE).** The loss of validation sets at the end of training all TCGA ADAGE models that appeared to converge.

It appears the models perform better without any induced sparsity.

This analysis allowed us to select optimal models based on tested hyperparameters.
For tied weights ADAGE, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate |
| :--------- | :------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 0.0 | 100 | 50 | 0.0015 |
| 25 | 0 | 0.0 | 100 | 50 | 0.0015 |
| 50 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 75 | 0 | 0.0 | 100 | 50  | 0.0005 |
| 100 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 125 | 0 | 0.0 | 100 | 50 | 0.0005 |

It appears that `learning rate` decreases for higher dimensional models, while `epochs` are globally optimal at 100; `batch size` at 50; and `noise` and `sparsity` at 0 (**Figure 6 (TCGA, ADAGE)**).
See https://github.com/greenelab/tybalt/issues/127 for more details about zero noise.

![](figures/tcga_results/z_parameter_best_model_adage_TCGA.png?raw=true)

**Figure 6 (TCGA, ADAGE).** Training optimal ADAGE models across different latent space dimensions.

### GTEx

#### Tybalt

Tybalt models in the GTEx dataset also showed variable performance, but they were generally stable across all hyperparameter combinations (**Figure 7 (GTEx, Tybalt)**).

![](figures/gtex_results/z_parameter_final_loss_tybalt_GTEX.png?raw=true)

**Figure 7 (GTEx, Tybalt).** The loss of validation sets at the end of training for all 540 GTEx Tybalt models.

Model performance (measured by observed validation loss) improved after increasing the capacity of the models from 5 to 125 dimensions.
However, the performance began to level off after around 50 dimensions.
All other hyperparameters had minor effects and varied with dimensionality, including `Kappa`.

This analysis allowed us to select optimal models based on tested hyperparameters.
For `GTEx` `Tybalt` models, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate |
| :--------- | :---- | :----- | :--------- | :------------ |
| 5 | 0.5 | 100 | 100 | 0.0025 |
| 25 | 0.5 | 100 | 100 | 0.0025 |
| 50 | 0.5 | 100 | 100 | 0.002 |
| 75 | 0.5 | 100 | 50  | 0.002 |
| 100 | 0.5 | 100 | 50 | 0.0015 |
| 125 | 0.5 | 100 | 50 | 0.0015 |

Training of models with optimal hyperparameters are shown in **Figure 8 (GTEx, Tybalt)**.

![](figures/gtex_results/z_parameter_best_model_tybalt_GTEX.png?raw=true)

**Figure 8 (GTEx, Tybalt).** Training optimal GTEx Tybalt models across different latent space dimensions.

### ADAGE

By constraining the compression and decompression networks to contain the same weights (tied weights), GTEx ADAGE models also had variable performance across models.
GTEx ADAGE models failed to converge with low learning rates (**Figure 9 (GTEx, ADAGE)**).

![](figures/gtex_results/z_parameter_final_loss_adage_GTEX.png?raw=true)

**Figure 9 (GTEx, ADAGE).** The loss of validation sets at the end of training for 144 tied weight GTEx ADAGE models.

After removing models that did not converge we see a clearer picture (**Figure 10 (GTEx, ADAGE)**).

![](figures/gtex_results/z_parameter_final_loss_remove_converge_adage_GTEX.png?raw=true)

**Figure 10 (GTEx, ADAGE).** The loss of validation sets at the end of training all GTEx ADAGE models that appeared to converge.

This analysis allowed us to select optimal models based on tested hyperparameters.
For tied weights ADAGE, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate |
| :--------- | :------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 0.1 | 100 | 50 | 0.001 |
| 25 | 0 | 0.0 | 100 | 50 | 0.001 |
| 50 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 75 | 0 | 0.0 | 100 | 50  | 0.0005 |
| 100 | 0 | 0.0 | 100 | 50 | 0.0005 |
| 125 | 0 | 0.0 | 100 | 50 | 0.0005 |

It appears that `learning rate` decreases for higher dimensional models, while `epochs` are globally optimal at 100; `batch size` at 50; and `sparsity` at 0 (**Figure 11 (GTEx, ADAGE)**).
We did observe added noise improving the lower dimensional models.

![](figures/gtex_results/z_parameter_best_model_adage_GTEX.png?raw=true)

**Figure 11 (GTEx, ADAGE).** Training optimal GTEx ADAGE models across different latent space dimensions.

### TARGET

#### Tybalt

Tybalt models in the TARGET dataset also showed variable performance, but they were generally stable across all hyperparameter combinations (**Figure 12 (TARGET, Tybalt)**).

![](figures/target_results/z_parameter_final_loss_tybalt_TARGET.png?raw=true)

**Figure 12 (TARGET, Tybalt).** The loss of validation sets at the end of training for all 540 TARGET Tybalt models.

Model performance (measured by observed validation loss) improved after increasing the capacity of the models from 5 to 125 dimensions.
However, the performance began to level off after around 50 dimensions.
All other hyperparameters had minor effects and varied with dimensionality, including `Kappa`.

This analysis allowed us to select optimal models based on tested hyperparameters.
For `TARGET` `Tybalt` models, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate |
| :--------- | :---- | :----- | :--------- | :------------ |
| 5 | 0.5 | 100 | 25 | 0.0015 |
| 25 | 0.5 | 100 | 25 | 0.0015 |
| 50 | 0.5 | 100 | 25 | 0.0015 |
| 75 | 0.5 | 100 | 25  | 0.0015 |
| 100 | 0.5 | 100 | 25 | 0.0015 |
| 125 | 0.5 | 100 | 25 | 0.0005 |

Training of models with optimal hyperparameters are shown in **Figure 13 (TARGET, Tybalt)**.

![](figures/target_results/z_parameter_best_model_tybalt_TARGET.png?raw=true)

**Figure 13 (TARGET, Tybalt).** Training optimal TARGET Tybalt models across different latent space dimensions.

### ADAGE

TARGET ADAGE models also had variable performance across models.
In contrast to TCGA, and GTEx, all TARGET models appeared to converge.
This may be because of the smaller, and less complex, data (**Figure 14 (TARGET, ADAGE)**).

![](figures/target_results/z_parameter_final_loss_adage_TARGET.png?raw=true)

**Figure 14 (TARGET, ADAGE).** The loss of validation sets at the end of training for 216 tied weight GTEx ADAGE models.

This analysis allowed us to select optimal models based on tested hyperparameters.
For TARGET ADAGE, the optimal hyperparameters across dimensionality estimates were globally consistent:

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate |
| :--------- | :------- | :---- | :----- | :--------- | :------------ |
| 5 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 25 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 50 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 75 | 0 | 0.1 | 100 | 50  | 0.0005 |
| 100 | 0 | 0.1 | 100 | 50 | 0.0005 |
| 125 | 0 | 0.1 | 100 | 50 | 0.0005 |

It appears that `learning rate` is globally optimum at 0.0005; `epochs` at at 100; `batch size` at 50; `noise` at 0.1; and `sparsity` at 0 (**Figure 15 (TARGET, ADAGE)**).

![](figures/target_results/z_parameter_best_model_adage_TARGET.png?raw=true)

**Figure 15 (TARGET, ADAGE).** Training optimal TARGET ADAGE models across different latent space dimensions.

## Summary

Selection of hyperparameters across different latent space dimensionality operated as expected.
Loss was higher for lower dimensions and lower dimensions benefited the most from increased regularization and higher learning rates.
Nevertheless, we have obtained a broad set of optimal hyperparameters for use in a larger and more specific sweep of dimensionality for each of the three analyzed datasets.
