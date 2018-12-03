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
# and takes five positional arguments as input:
#
#         --dataset             Name of the Dataset
#         --gmt_name            The name of the GMT data
#                                 (e.g. "xcell_all_entrez.gmt")
#         --metaedge            The name of the metaedge used (e.g. "GpH")
#         --gene_set_dir        The input directory of results
#         --shuffled            If the input data is shuffled
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
                                          help = "if top results are saved"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Load arguments
dataset <- opt$dataset
gmt_name <- opt$gmt_name
metaedge <- opt$metaedge
gene_set_dir <- opt$gene_set_dir
shuffled <- opt$shuffled
save_results <- opt$save_results

# Find the name of all of the gene sets of interest
xcell_file <- file.path("..", "3.build-hetnets", "data", gmt_name)

con <- file(xcell_file, open = "r")

genesets <- c()
while (length(geneset <- readLines(con, n = 1, warn = FALSE)) > 0) {
  geneset <- unlist(strsplit(geneset, "\t"))[1]
  genesets <- unique(c(geneset, genesets))
}

close(con)

for (gene_set in genesets) {
  top_results_df <- plot_gene_set(gene_set = gene_set,
                                  gene_set_dir = gene_set_dir,
                                  metaedge = metaedge,
                                  dataset = dataset,
                                  show_plot = FALSE,
                                  shuffled = shuffled)

  if (save_results) {
    out_file <- paste(dataset, metaedge, paste0(gene_set, ".tsv"), sep = "_")
    out_file <- file.path("results", "top_features", out_file)
    readr::write_tsv(x = top_results_df, path = out_file)
  }
}
