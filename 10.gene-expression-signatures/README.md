# Detecting Gene Expression Signatures

**Gregory Way 2019**

In this module we perform a series of t-tests to detect sample sex in GTEx samples and TCGA patients based on gene expression profiles.
We also test for signatures representing MYCN amplification status in neuroblastoma (NBL) tumors using TARGET data.

![signatures](https://raw.githubusercontent.com/greenelab/BioBombe/master/10.gene-expression-signatures/figures/full_separation_plot.png)

## Reproducible Analysis

To reproduce the results of the gene expression signature analysis perform the following:

```bash
# Activate computational environment
conda activate biobombe

# Perform the analysis
cd 10.gene-expression-signatures
./signature_analysis.sh
```
