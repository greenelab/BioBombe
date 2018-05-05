# Results for Initial Hyperparameter Sweep Considering Different Latent Dimensionality

**Gregory Way 2018**

## Latent Space Dimensionality

Compression algorithms reduce the dimensionality of input data by enforcing the number of dimensions to bottleneck.
A common problem is the decision of how many "useful" latent space features are present in data.
The solution is different for different problems or goals.
For example, when visualizing large differences between groups of data, a highly restrictive bottleneck, usually between 2 or 3 features, is required.
However, when the goal is to extract meaningful patterns in the data that may have more subtle relationships across samples, the recommendations are opaque.
As the bottleneck relaxes, the ability to explain the patterns decreases and the possibility of false positives increases.

In order to determine an optimal _range_ of compression dimensions, we propose the following.
We will first sweep over various different dimensions (results provided below) and perform several evals (to be described later).

Before sweeping over a large number of different dimensions, we perform a hyperparameter sweep of select dimensions.
In this sense, we want to minimize the effect of poor hyperparameter combinations across different dimensions contributing to performance differences.
In other words, we want to isolate the effect of changing dimensionality on the observed patterns and solutions.
Therefore, we perform a parameter sweep over several hyperparameters for both Tybalt and ADAGE models below.

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

We sweep over the following parameter combinations for Tybalt and ADAGE models:

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
Importantly, we also include results of a parameter sweep of 1,080 ADAGE models with _untied_ weights, which we ran previously.
For all downstream applications we use ADAGE models with _tied_ weights, but we also report the _untied_ results here.

Our goal was to determine optimal hyperparameter combinations for both models across various bottleneck dimensionalities.

## Results

We report the results in a series of visualizations and tables for Tybalt and ADAGE separately below.

In order to compile the results of the parameter sweep, run the following commands:

```bash
# Compile Tybalt parameter sweep results
python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_tybalt/' --output 'parameter_sweep_tybalt_full_results.tsv'

# Compile ADAGE parameter sweep results
python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_adage/' --output 'parameter_sweep_adage_tiedweights_full_results.tsv'

# Compile untied ADAGE parameter sweep results
python scripts/summarize_paramsweep.py --results_directory 'param_sweep/param_sweep_adage_untied' --output 'parameter_sweep_adage_full_results.tsv'

# Visualize the results of the sweep for all models
Rscript --vanilla scripts/param_sweep_latent_space_viz.R
```

### Tybalt

Tybalt models had variable performance across models, but was generally stable across all hyperparameter combinations (**Figure 1**).

![](figures/z_param_tybalt/z_parameter_tybalt.png?raw=true)

**Figure 1.** The loss of validation sets at the end of training for all 540 Tybalt models.

Model performance (measured by observed validation loss) improved after increasing the capacity of the models from 5 to 125 dimensions.
However, the performance began to level off after around 50 dimensions.
All other hyperparameters had minor effects and varied with dimensionality, including `Kappa`.

Selecting a constant learning rate, batch size, and epochs for all models, we stratify the entire training process across all tested dimensionalities (**Figure 2**).

![](figures/z_param_tybalt/z_parameter_tybalt_training.png?raw=true)

**Figure 2.** Validation loss across all training epochs for a fixed combination of hyperparameters.
The `learning rate` was set to 0.0005, `batch size` 50, and `epochs` 100.
We also show performance differences across `kappa`.

This analysis allowed us to select optimal models based on tested hyperparameters.
For Tybalt, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Kappa | Epochs | Batch Size | Learning Rate | End Loss |
| :--------- | :---- | :----- | :--------- | :------------ | :------- |
| 5 | 0.5 | 100 | 50 | 0.001 | 2525.4 |
| 25 | 0 | 100 | 50 | 0.001 | 2479.6 |
| 50 | 0 | 100 | 100 | 0.001 | 2465.5 |
| 75 | 0 | 100 | 150  | 0.001 | 2460.5 |
| 100 | 0 | 100 | 150 | 0.0005 | 2456.5 |
| 125 | 0 | 100 | 150 | 0.0005 | 2457.4 |

Generally, it appears that the optimal `learning rate` and `kappa` decreases while the `batch size` increases as the dimensionality increases.
Training of models with optimal hyperparameters are shown in **figure 3**.

![](figures/z_param_tybalt/z_parameter_tybalt_best.png?raw=true)

**Figure 3.** Training optimal Tybalt models across different latent space dimensions.

### ADAGE

### Untied Weights

With untied weights, ADAGE models had variable performance across models and failed to converge with high levels of sparsity (**Figure 4**).
High levels of sparsity fail _worse_ with increasing dimensionality.

![](figures/z_param_adage/z_parameter_adage.png?raw=true)

**Figure 4.** The loss of validation sets at the end of training for all 1,080 untied weight ADAGE models.

After removing `sparsity = 0.001`, we see a clearer picture (**Figure 5**).

![](figures/z_param_adage/z_parameter_adage_remove_sparsity.png?raw=true)

**Figure 5.** The loss of validation sets at the end of training for 720 untied weight ADAGE models.

A similar pattern appears where lower dimensionality benefits from increased sparsity.
ADAGE models are also generally stable, particularly at high dimensions.

It appears that `learning rate` is globally optimal at 0.0005; epochs at 100; batch size at 50; sparsity at 0; with decreasing noise for larger z dimensions.

![](figures/z_param_adage/z_parameter_adage_bes.png?raw=true)

**Figure 6.** Training optimal untied weight ADAGE models across different latent space dimensions.

### Tied Weights

By constrianing the compression and decompression networks to contain the same weights (tied weights), ADAGE models had variable performance across models.
ADAGE models failed to converge with low learning rates (**Figure 7**).

![](figures/z_param_adage_tied_weights/z_parameter_adage_tiedweights.png?raw=true)

**Figure 7.** The loss of validation sets at the end of training for 648 tied weight ADAGE models.

After removing `learning rate = 1e-05` and `learning_rate = 5e-05`, we see a clearer picture (**Figure 8**).

![](figures/z_param_adage_tied_weights/z_param_adage_remove_learningrate_tiedweights.png?raw=true)

**Figure 8.** The loss of validation sets at the end of training for 432 tied weight ADAGE models.

It appears the models perform better without any induced sparsity

This analysis allowed us to select optimal models based on tested hyperparameters.
For tied weights ADAGE, the optimal hyperparameters across dimensionality estimates are:

| Dimensions | Sparsity | Noise | Epochs | Batch Size | Learning Rate | End Loss |
| :--------- | :------- | :---- | :----- | :--------- | :------------ | :------- |
| 5 | 0 | 0.0 | 100 | 50 | 0.0015 | 0.0042 |
| 25 | 0 | 0.0 | 100 | 50 | 0.0015 | 0.0029 |
| 50 | 0.0 | 0 | 100 | 50 | 0.0005 | 0.0023 |
| 75 | 0 | 0.0 | 100 | 50  | 0.0005 | 0.0019 |
| 100 | 0 | 0.0 | 100 | 50 | 0.0005 | 0.0017 |
| 125 | 0 | 0.0 | 100 | 50 | 0.0005 | 0.0016 |

It appears that `learning rate` decreases for higher dimensional models, while epochs are globally optimal at 100; batch size at 50; and noise and sparsity at 0.
See https://github.com/greenelab/tybalt/issues/127 for more details about zero noise.

![](figures/z_param_adage_tied_weights/z_parameter_adage_best_tiedweights.png?raw=true)

**Figure 9.** Training optimal ADAGE models across different latent space dimensions.

## Summary

Selection of hyperparameters across different latent space dimensionality operated as expected.
In general, tied weight ADAGE models performed better than untied weight ADAGE models, and required less regularization.
We will use tied weight ADAGE in all downstream analyses.
Loss was higher for lower dimensions and lower dimensions benefited the most from increased regularization.
Nevertheless, we have obtained a broad set of optimal hyperparameters for use in a larger and more specific sweep of dimensionality.
