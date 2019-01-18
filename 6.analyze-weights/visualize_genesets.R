# Gregory Way 2018
# 6.analyze-weights/visualize_genesets.R
#
# Description:
# Visualizes the results of the matrix approach to gene set compressed feature
# intepretation.

# Usage:
# The script is run in the command line
#
#     Rscript --vanilla visualize_genesets.R
#
# and takes seven arguments as input:
#
#         --dataset             Name of the Dataset
#         --gmt_name            The name of the GMT data
#                                 (e.g. "xcell_all_entrez.gmt")
#         --metaedge            The name of the metaedge used (e.g. "GpH")
#         --gene_set_dir        The input directory of results
#         --shuffled            If the input data is shuffled
#         --save_results        If included, save intermediate results
#         --no_plot             If included, do not create plots
#
# Output:
# Significantly overrepresented pathways from a WebGestalt Analysis

library(ggplot2)
library(dplyr)

source(file.path("scripts", "utils.R"))

# Load in command arguments
option_list <- list(optparse::make_option(c("-d", "--dataset"),
                                          type = "character",
                                          help = "Name of the dataset"),
                    optparse::make_option(c("-g", "--gmt_name"),
                                          type = "character",
                                          help = "name of the gmt file"),
                    optparse::make_option(c("-m", "--metaedge"),
                                          type = "character",
                                          help = "name of the metaedge"),
                    optparse::make_option(c("-i", "--gene_set_dir"),
                                          type = "character",
                                          help = "input directory location"),
                    optparse::make_option(c("-s", "--shuffled"),
                                          type = "character",
                                          action = "store_true",
                                          default = FALSE,
                                          help = "if input data is shuffled"),
                    optparse::make_option(c("-r", "--save_results"),
                                          type = "character",
                                          action = "store_true",
                                          default = FALSE,
                                          help = "if top results are saved"),
                    optparse::make_option(c("-p", "--plot_results"),
                                          type = "character",
                                          action = "store_true",
                                          default = TRUE,
                                          help = "if plots are saved"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Load arguments
dataset <- opt$dataset
gmt_name <- opt$gmt_name
metaedge <- opt$metaedge
gene_set_dir <- opt$gene_set_dir
shuffled <- opt$shuffled
save_results <- opt$save_results
plot_results <- opt$plot_results

# Find the name of all of the gene sets of interest
xcell_file <- file.path("..", "3.build-hetnets", "data", gmt_name)

con <- file(xcell_file, open = "r")

genesets <- c()
while (length(geneset <- readLines(con, n = 1, warn = FALSE)) > 0) {
  geneset <- unlist(strsplit(geneset, "\t"))[1]
  genesets <- unique(c(geneset, genesets))
}

close(con)

# Process the top results
top_features_df <- get_top_biobombe_results(gene_set_dir = gene_set_dir)

# Save results to file
if (save_results) {
  out_file <- paste(dataset, metaedge, "top_biobombe_scores.tsv.gz", sep = "_")
  out_file <- file.path("results", "top_features", out_file)
  readr::write_tsv(x = top_features_df, path = out_file)
}

if (plot_results) {
  for (gene_set in genesets) {
    print(paste("curating results for:", gene_set))
    top_results_df <- plot_gene_set(gene_set = gene_set,
                                    full_results_df = top_features_df,
                                    metaedge = metaedge,
                                    dataset = dataset,
                                    show_plot = TRUE,
                                    shuffled = shuffled)
  }
}

