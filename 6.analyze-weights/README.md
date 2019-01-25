# BioBombe Network Projection

**Gregory Way 2019**

This is one of the most important modules.
Here, we perform the BioBombe Network Projection procedure on select gene set collections and gene expression gene sets.

## Gene set collections and data sets used

We performed the network projection on the following:

| Dataset | Collection |
| :------ | :--------- |
| TCGA | GpH |
| TCGA | GpXCELL |
| TCGA | GpC4CM |
| TCGA | GpC2CPREACTOME |
| TCGA | GpC3TFT |
| TARGET | GpH |
| TARGET | GpXCELL |
| TARGET | GpC4CM |
| GTEX | GpXCELL |

## Reproducible Analysis

To reproduce the network projection results, perform the following:

```bash
# Activate computational environment
conda activate biobombe

# Perform the analysis
cd 6.analyze-weights
python interpret-compression.py

# Visualize feature dimension and rank
Rscript --vanilla scripts/nbconverted visualize-top-feature-dimension-and-rank.r
```

## Results

The `interpret-compression.py` script performs the interpretation of the compressed gene expression features.
These results are used in several downstream applications.
Additionally, we visualize what dimensions specific gene set collection features are learned, and to what degree.

### Main Figures

![Density Figure](https://raw.githubusercontent.com/greenelab/BioBombe/master/6.analyze-weights/figures/top_feature_density.png)

Shown above are the latent dimensions (_k_) that have the highest enrichment for the particular gene set and data set.
We were curious about the dimensions of specific gene sets, and if there was any difference to what dimension they were discovered at to the highest enrichment.
We also show GTEx GpXCELL relationships for stacked and relative densities.
It is clear that the bulk of the gene sets for all algorithms (particularly NMF) are learned at high dimensions.

![Rank Figure](https://raw.githubusercontent.com/greenelab/BioBombe/master/6.analyze-weights/figures/top_feature_rank.png)

We were also interested in the distribution of gene set enrichment ranks across algorithms and dimensions.
The above plot shows all gene set ranked by BioBombe z score.
