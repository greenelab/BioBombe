# Generating Hetergeneous Networks for Compression Interpretation

**Gregory Way 2018**

This module downloads and processes several gene sets and integrates these gene sets into a heterogeneous network (hetnet)([Himmelstein et al. 2017](https://doi.org/10.7554/eLife.26726 "Systematic integration of biomedical knowledge prioritizes drugs for repurposing")).

This network will be projected onto compressed gene expression features to enable biological interpretation.

## Data

The module stores data and data processing scripts of [MSigDB](http://software.broadinstitute.org/gsea/msigdb/index.jsp) and [xCell](https://doi.org/10.1186/s13059-017-1349-1) gene sets.

### Molecular Signatures Database (MSigDB)

Individual MSigDB gene sets (version 6.1) were downloaded from [GSEA downloads](http://software.broadinstitute.org/gsea/downloads.jsp).
We also download the full gene set: `msigdb.v6.1.entrez.gmt`.

See [download_msigdb.ipynb](download_msigdb.ipynb) for specific details.

The genesets consist of 8 different collections; many also have sub-collections:

| Name | Collection | License |
| :--: | :--------- | ------: |
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

We do not include the KEGG, BioCarta, AAAS/STKE gene sets in `C2`.
For full license terms visit the [MSigDB license page](http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp).

### xCell

We download and process the 489 gene signatures from [Arun et al. 2017](https://doi.org/10.1186/s13059-017-1349-1 "xCell: digitally portraying the tissue cellular heterogeneity landscape").
These 489 signatures represent 64 different human cell types including CD8+ T Cells, Neutrophils, Macrophages, etc.

See [process_xCell.ipynb](process_xCell.ipynb) for specific details.

## Integration

The gene sets are uniformly processed and linked together in a single hetnet.
The hetnet can be projected onto a compressed gene expression matrix to rapidly assign biological significance to compressed gene expression features.
