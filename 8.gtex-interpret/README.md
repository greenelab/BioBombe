# GTEx BioBombe Application

**Gregory Way 2019**

## Exploring increased correlation in GTEx Blood with increased model capacity

Previously, we recognized a sharp increase in the ability of variational autoencoders (VAEs) to capture the correlational structure of blood between z = 2 and z = 3.
This increase in model capacity by one, increases the mean correlation of samples by nearly 0.3.
We do not observe this pattern for other algorithms at this change in dimension.

![sample-correlation_Blood_GTEX_signal_pearson.png](https://raw.githubusercontent.com/greenelab/BioBombe/master/4.analyze-components/figures/GTEX/sample-correlation/sample-type/sample-correlation_Blood_GTEX_signal_pearson.png)

We do, however, also notice large improvements in other algorithms, but at other dimensions.
For example, we observe a large increase in blood correlation for PCA and ICA between 3 and 4, for NMF between 6 and 7, and, surprisingly, a large decrease in performance for DAE between 8 and 9.
We chose to explore the correlational change for VAE models because there are fewer features to interrogate.
See panel D below:

![gtex_blood_correlation_increase](https://raw.githubusercontent.com/greenelab/BioBombe/master/4.analyze-components/figures/correlation_summary.png)

The samples are processed and scores applied in [0.process-gtex-blood-vae-feature.ipynb](0.process-gtex-blood-vae-feature.ipynb).

The relationship of these initial VAE features are visualized in [1.visualize-gtex-blood-interpretation.ipynb](1.visualize-gtex-blood-interpretation.ipynb).

## BioBombe network projection application to VAE features

We applied our network projection approach to two VAE models (z=2 and z=3) using a network built from xCell genesets.
Many cell-type signatures were implicated in both VAE models including skeletal muscle, neurons, keratinocytes, and sebocytes.
However, a neutrophil signature was extracted from the VAE z = 3 model and not the VAE z = 2 model (**Panel A Below**).

Furthermore, if we take the mean of all gene sets across all features in both models, we can identify additional geneset signatures which may also be contributing to the full model performance.
Using the approach, we implicate many of the same cell-type signatures, but we also unviel a set of monocyte signatures that are more enriched in VAE z = 3 than z = 2.
The signature `Monocytes_FANTOM_2` appears to have the lowest enrichment in VAE z = 2 and a relatively high enrichment in z = 3 (**Panel B Below**).

Therefore, we elected to follow up with two features implicated in the VAE model z = 3 that may be helping with the sharp increase in correlation in blood tissues with a single increase in model capacity.
The genesets we chose to follow up with are `Neutrophils_HPCA_2` and `Monocytes_FANTOM_2`.

 ![gtex_main_figure](https://raw.githubusercontent.com/greenelab/BioBombe/master/8.gtex-interpret/figures/gtex_biobombe_main_figure.png)

## Tracking Neutrophil and Monocyte signatures across algorithms and dimensions

Each algorithm at all dimensions received a score for all xCell genesets.
We tracked this score and can visualize the enrichment as z increases for `Neutrophils_HPCA_2` (**Panel C Above**) and `Monocytes_FANTOM_2` (**Panel D Above**).

Evidently, the scores seem to improve for all algorithms as the dimensions increase, but there are spikes at intermediate dimensions.
We selected the top scoring feature for both genesets and applied these (along with the origin z = 3 high scoring features) to external datasets.

The scores for the highest scoring `Neutrophils_HPCA_2` geneset feature was feature 10 in VAE z = 14.
The geneset features enriched in this feature were highly associated with the genesets enriched in VAE z = 3 (**Panel E Above**).
The gene weights themselves were also highly enriched (**Panel F Above**).

## Validating signatures derived from compression algorithms on external datasets

### Neutrophils

We downloaded processed gene expression data from [GSE103706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103706) ([Rincon et al. 2018](https://doi.org/10.1186/s12864-018-4957-6)).

In this dataset there are two acute myeloid leukemia (AML) cell lines; PLB-985 and HL-60, and a total of 14 samples.
The cell lines are exposed to two treatments - DMSO and DMSO+Nutridoma - plus replicates with no treatment applied.
The treatments are demonstrated to induce neutrophil differentiation in these cell lines.

We hypothesized that our constructed feature identified through our interpret compression approach would have higher activation patterns in the cell lines with induced neutrophil differentiation.

The data is downloaded and processed in [2A.download-neutrophil-data.ipynb](2A.download-neutrophil-data.ipynb).

### Monocytes

We downloaded processed gene expression data from [GSE24759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24759) ([Novershtern et al. 2011](https://doi.org/10.1016/j.cell.2011.01.004)).

In this dataset there are 211 samples consisting of 38 distinct hematopoietic states in various stages of differentiation.

We hypothesized that our constructed feature identified through our interpret compression approach would have higher activation patterns in Monocytes.

The data is downloaded and processed in [2B.download-hematopoietic-data.ipynb](2A.download-hematopoietic-data.ipynb).

### Application of Signatures

The BioBombe derived compression signatures are applied to the external datasets in [3.apply-signatures.ipynb](3.apply-signatures.ipynb).

## Signatures reveal consistent activation patterns in neutrophils and monocytes

### Neutrophils

The two VAE features (feature 0 in VAE z = 3 and feature 10 in VAE z = 14) were applied to GSE103706 (**Panel G Above**).
Feature 0 from VAE z = 3 and feature 10 from VAE z = 14 had the highest scores for the specific neutrophil geneset.
It does not appear that the z = 3 feature was able to robustly separate the the two treatments from the untreated cell line controls. However, the untreated controls were tending towards negative scores.
The z = 14 feature perfectly separated the untreated cell lines from the treated cell lines and therefore validate the interpretation approach.

### Monocytes

The two features (feature 2 in VAE z = 3 and feature 200 in NMF z = 200) were applied to GSE24759 (**Panel H Above**).
Both features showed the highest scores in isolated monocytes.
Indeed, it appears that, specifically, `Mono2` cells were particularly enriched.
Granulocytes also had high scores along this axis, indicating similarity between specific hematopoietic differentiation stages.

## Reproducible Analysis

To reproduce the results of the stability analysis perform the following:

```bash
# Activate computational environment
conda activate biobombe

# Perform the analysis
cd 8.gtex-interpret
./gtex_analysis.sh
```
