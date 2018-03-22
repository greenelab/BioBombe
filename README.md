# Interpreting Compressed Gene Expression Features 

**Gregory Way 2018**

**University of Pennsylvania**

The repository stores data and data processing scripts of [MSigDB](http://software.broadinstitute.org/gsea/msigdb/index.jsp) gene sets.
These data are used in methods to assign biology to compressed gene expression features.

## Data Access

Download individual MSigDB gene sets (version 6.1) from data from [GSEA downloads](http://software.broadinstitute.org/gsea/downloads.jsp).
Note that you will have to [register](http://software.broadinstitute.org/gsea/register.jsp?next=index.jsp) first.
Place these genesets in the `data/` folder.
The genesets consist of 8 different collections:

| Name | Collection | License |
| :--: | :--------: | ------: |
| H | Hallmark gene sets | CC-BY 4.0 |
| C1 | Positional gene sets | CC-BY 4.0 |
| C2 | Curated gene sets | CC-BY 4.0 (except KEGG, BioCarta, AAAS/STKE) |
| C3 | Motif gene sets | CC-BY 4.0 |
| C4 | Computational gene sets | CC-BY 4.0 |
| C5 | Gene Ontology terms | CC-BY 4.0 |
| C6 | Oncogenic gene sets | CC-BY 4.0 |
| C7 | Immunologic gene sets | CC-BY 4.0 |

We do not include the KEGG, BioCarta, AAAS/STKE gene sets in `C2` for easy dissemination.
For full license terms visit the [MSigDB license page](http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp).

