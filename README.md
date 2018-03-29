# Interpreting Compressed Gene Expression Features 

**Gregory Way 2018**

**University of Pennsylvania**

The repository stores data and data processing scripts of [MSigDB](http://software.broadinstitute.org/gsea/msigdb/index.jsp) gene sets.
These data are used in methods to assign biology to compressed gene expression features.

## Data Access

Individual MSigDB gene sets (version 6.1) were downloaded from [GSEA downloads](http://software.broadinstitute.org/gsea/downloads.jsp).
Also download the full gene set: `msigdb.v6.1.symbols.gmt`.

See [download_msigdb.ipynb](download_msigdb.ipynb) for specific details.
 
The genesets consist of 8 different collections, with many have sub-collections:

| Name | Collection | License |
| :--: | :--------: | ------: |
| H | Hallmark gene sets | CC-BY 4.0 |
| C1 | Positional gene sets | CC-BY 4.0 |
| C2 | Curated gene sets | CC-BY 4.0 (except KEGG, BioCarta, AAAS/STKE) |
| C2.CPG | Chemical and genetic perturbations | CC-BY 4.0 |
| C2.CP.Reactome | Reactome | CC-BY 4.0 |
| C3 | Motif gene sets | CC-BY 4.0 |
| C3.MIR | microRNA targets | CC-BY 4.0 |
| C3.TFT | Transcription factor targets | CC-BY 4.0 |
| C4 | Computational gene sets | CC-BY 4.0 |
| C4.CGN | Cancer gene neighborhoods | CC-BY 4.0 |
| C4.CM | Cancer modules | CC-BY 4.0 |
| C5 | Gene Ontology (GO) terms | CC-BY 4.0 |
| C5.BP | GO biological processes | CC-BY 4.0 |
| C5.CC | GO cellular components | CC-BY 4.0 |
| C5.MF | GO molecular functions | CC-BY 4.0 |
| C6 | Oncogenic gene sets | CC-BY 4.0 |
| C7 | Immunologic gene sets | CC-BY 4.0 |

We do not include the KEGG, BioCarta, AAAS/STKE gene sets in `C2` for easy dissemination.
For full license terms visit the [MSigDB license page](http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp).

## Computational Environment

All processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.5.0
conda env create --force --file environment.yml

conda activate interpret-compression
```
