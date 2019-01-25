# TCGA BioBombe Application

**Gregory Way 2019**

## Predicting cancer type and gene alteration status using compressed features

In this module, we use the compressed latent features from 5 algorithms and across 28 different dimensions, to predict cancer type and gene alteration status in The Cancer Genome Atlas (TCGA) gene expression data.

We predict cancer types and gene alterations using a logistic regression algorithm penalized with elastic net.
We previously used a similar approach to predict Ras pathway activation in the same dataset ([Way et al. 2018](https://doi.org/10.1016/j.celrep.2018.03.046)).
Rather than using raw gene expression features, we predict cancer types and mutations using compressed gene expression features.

Our goal in this analysis was to determine the ability of the constructed features to learn signals representing cancer type and gene alterations.
Our aim was to observe at which dimension and at what strength the signals appeared in the compression models.
Furthermore, the compression algorithms aggregate groups of genes together that represent coordinated biological processes.
Observing how these compressed features participated in the models (by model coefficients) we can make statements about the compressed features biologically.

**Notes**

* Samples with a hypermutater phenotype were removed from the dataset (greater than 5 standard deviations over the mean of log10 mutation rate)
* All models are adjusted for log10 mutation rate
* Gene alteration predictions only:
  * These models are also adjusted by cancer type
  * To regulate balanced class sizes, only cancer types with greater than 15% and more than 5 samples in either positive or negative gene alteration class were used to train models
  * Gene alteration was defined as any non-silent mutation in the gene of interest, or a copy number amplification for an oncogene, or copy number deletion for a tumor suppressor gene.

### Predicting cancer types

TCGA has profiled 33 different cancer-types.
We constructed supervised machine learning tasks to distinguish each cancer-type from all others.

### Predicting gene alterations

We selected the [top 50 most mutated genes](https://github.com/greenelab/BioBombe/blob/master/9.tcga-classify/top-50-pancanatlas-mutations.ipynb) in TCGA to predict.
Using non silent mutation and copy number status, we predicted the alteration status of these genes in the included cancer types.

## Compression classification results

All of our results are deposited in a [publicly available resource](https://doi.org/10.5281/zenodo.2535759).

![TCGA_BioBombe_Results](https://raw.githubusercontent.com/greenelab/BioBombe/master/9.tcga-classify/figures/tcga_biobombe_main_figure.png)

We highlight predictions in 5 select cancer types and 5 important cancer genes.
**Panel A** describes performance for predicting cancer types, while **Panel B** displays gene alteration prediction performance.
For both panels, the red line represents permuted data input into the compression algorithms, and blue lines are real data.
The red and blue dotted lines represents prediction performance with raw gene expression data, and the grey dotted line represents a hypothetical guess.
Note that the permuted models include covariate information.

It appears that all of these cancer types can be predicted with 100% accuracy with raw gene expression features and at early to intermediate compression dimensions.
However, the gene alteration status are predicted with modest accuracy, and signal only starts to emerge at intermediate to late compression dimensions.
This indicates that cancer-type signatures are features learned at early dimensions and are some of the highest sources of variation in this dataset, while mutation signatures are learned in higher dimensions.

### Sparsity

We also track the sparsity, as the percentage of the model coefficients that are zero, for real and compressed predictions of gene alteration status (**Panel C Above**).
The grey square represents the models with raw gene expression.
In all cases, it appears that the models with low dimension are not performing as well as the models with high dimensions.
It also appears that denoising autoencoders (DAE) are the sparsest models.

### Exploring a Denoising Autoencoder Model

We followed up with one of the sparse DAE models used to predict TP53 inactivation.
The feature we followed up included 200 compression features and 22 covariates (log10 mutations and 21 cancer type dummy variables) (`seed 908341`) and had 80.6% of the features with zero coefficients (179 out of 222).
Therefore 43 coefficients were non zero.
The top ranked coefficients are shown in **Panel D Above**.

We applied the network projection interpretation approach to the top DAE coefficients used in predicting TP53 alterations.
The top ranked `GpH` gene sets for these features included:

| Full Feature | Absolute Value Z Score | Gene Set                        | raw score | z_score  | ML coefficient |
|--------------|-------------|--------------------------------------------|-----------|----------|----------------|
| dae_156      | 4.48129     | HALLMARK_INTERFERON_GAMMA_RESPONSE         | 1.46557   | 4.48129  | 0.23558  |
| dae_156      | 4.92812     | HALLMARK_COAGULATION                       | -2.63723  | -4.92812 | 0.23558  |
| dae_122      | 7.87064     | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 5.53091   | 7.87064  | 0.19832  |
| dae_122      | 4.64101     | HALLMARK_KRAS_SIGNALING_DN                 | -1.28863  | -4.64101 | 0.19832  |
| dae_7        | 5.68651     | HALLMARK_MYC_TARGETS_V1                    | 3.34075   | 5.68651  | 0.16581  |
| dae_7        | 11.0546     | HALLMARK_MYOGENESIS                        | -6.15819  | -11.0546 | 0.16581  |
| dae_186      | 4.65232     | HALLMARK_INTERFERON_GAMMA_RESPONSE         | 2.93655   | 4.65232  | 0.14635  |
| dae_186      | 5.05547     | HALLMARK_FATTY_ACID_METABOLISM             | -2.34884  | -5.05547 | 0.14635  |
| dae_143      | 10.7977     | HALLMARK_COMPLEMENT                        | 2.39042   | 10.7977  | 0.13902  |
| dae_143      | 3.2203      | HALLMARK_E2F_TARGETS                       | -2.33799  | -3.2203  | 0.13902  |
| dae_88       | 7.40728     | HALLMARK_ESTROGEN_RESPONSE_EARLY           | 3.25406   | 7.40728  | -0.15494 |
| dae_88       | 6.57371     | HALLMARK_XENOBIOTIC_METABOLISM             | -3.46226  | -6.57371 | -0.15494 |
| dae_130      | 4.19995     | HALLMARK_OXIDATIVE_PHOSPHORYLATION         | 0.536709  | 4.19995  | -0.16391 |
| dae_130      | 9.04101     | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | -11.3443  | -9.04101 | -0.16391 |
| dae_24       | 3.79301     | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 2.48863   | 3.79301  | -0.18192 |
| dae_24       | 6.65545     | HALLMARK_KRAS_SIGNALING_DN                 | -2.61153  | -6.65545 | -0.18192 |
| dae_46       | 4.03099     | HALLMARK_INTERFERON_ALPHA_RESPONSE         | 2.27472   | 4.03099  | -0.22864 |
| dae_46       | 7.15103     | HALLMARK_OXIDATIVE_PHOSPHORYLATION         | -3.15448  | -7.15103 | -0.22864 |

These genesets represent the highest scoring hallmark processes activated and deactivated in samples with TP53 activation.

#### Model Performance

Receiver operating characteristic (ROC) curves (**Panel E Above**) and precision recall curves (**Panel F Above**) are presented for training, testing, and cross validation (CV) data partitions for predictions using raw gene expression features and compressed gene expression features.

## Predicting Gene Mutation Status with Top Scoring Compressed Features

We were also interested in observing performance using hand picked features.

![TCGA_BioBombe_Supplementary_Results](https://raw.githubusercontent.com/greenelab/BioBombe/master/9.tcga-classify/figures/supplemental_tcga_top_feature_summary.png)

We compared AUROC (**Panel A Above**) and AUPR (**Panel B Above**) in cross validation (CV) training folds for predictions using k = 200 features to 200 of the top scoring features across all tested k dimensions.
The height of the bar represents performance in these top 200 features.
We also compared predictions using the top 1 feature and 200 randomly selected features.
Also shown are predictions using raw RNAseq data.

We determined top features based on coefficient weights in their respective prediction tasks.
Our hypothesis was that different compression algorithms would aggregate together various signals applicable to the gene alteration prediction task.

Overall, it does not appear that aggregating top features improves performance, especially in linear models.
Also, aggregating features together removes some prediction signal, as the raw data predictions have higher performance than all predictions with compressed features.


## Reproducible Analysis

To reproduce the results of the TCGA classification analysis perform the following:

```bash
# Activate computational environment
conda activate biobombe

# Perform the analysis
cd 9.tcga-classify
./classify_analysis.sh
```
