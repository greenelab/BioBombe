"""
2018 Gregory Way
5.analyze-components/scripts/submit_gsea.py

Perform GSEA on all compiled weight matrices in the given directory

Usage:

    Run in command line:

            python scripts/submit_gsea.py

     with required command arguments:

       --input_weight_dir   directory storing ensemble weight matrices
       --z_dim              the given z dimension (bottleneck layer size)
                            (options: TCGA, GTEX or TARGET)
       --dataset_name       a string describing the name of the dataset used
       --out_dir            filepath of where to save the results
       --num_perm           number of GSEA permutations
       --shuffled           boolean if the data was shuffled before fitting
       --algorithms         the algorithms compressed data to interpret
       --distrib_methods    various mechanisms that aggregate gene lists from
                            compressed features. These gene lists are input to
                            GSEA

    and optional command arguments:

        --gene_sets         A keyword from Enrichr or a file path to .gmt file
        --translate         if present, then the gene symbols in the weight
                            matrices will be translated

Output:

    GSEA results for each z dimension in `5.analyze-weights/results/gsea`
"""

import argparse
from latent import run_gsea_pipeline_command

parser = argparse.ArgumentParser()
parser.add_argument('-w', '--input_weight_dir',
                    help='location of the directory storing weight matrices')
parser.add_argument('-z', '--z_dim',
                    help='the bottleneck dimensionality to focus on')
parser.add_argument('-d', '--dataset_name', help='name of the dataset',
                    choices=['TCGA', 'TARGET', 'GTEX'])
parser.add_argument('-n', '--num_perm', help='number of GSEA permutations')
parser.add_argument('-s', '--shuffled', action='store_true',
                    help='if the signal is shuffled or real data')
parser.add_argument('-a', '--algorithms', help='the algorithms to consider')
parser.add_argument('-m', '--distrib_methods',
                    help='the distribution methods to perform')
parser.add_argument('-g', '--gene_sets', default='KEGG_2016', nargs='+',
                    help='gene sets to perform GSEA using')
parser.add_argument('-t', '--translate', action='store_true',
                    help='translate hugo symbol into entrez gene')
args = parser.parse_args()

# Set constants
input_weight_dir = args.input_weight_dir
z_dim = args.z_dim
dataset_name = args.dataset_name
num_perm = args.num_perm
shuffled = args.shuffled
algorithms = args.algorithms.split()
distrib_methods = args.distrib_methods.split()
gene_sets = args.gene_sets
translate = args.translate

# Perform the analysis
if __name__ == "__main__":
    run_gsea_pipeline_command(
        input_weight_dir=input_weight_dir,
        z_dim=z_dim,
        dataset_name=dataset_name,
        num_perm=num_perm,
        shuffled_true=shuffled,
        algorithms=algorithms,
        distrib_methods=distrib_methods,
        translate=translate,
        gene_sets=gene_sets
    )
