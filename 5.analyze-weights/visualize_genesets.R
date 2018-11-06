# Gregory Way 2018
# visualize_genesets.R
#
# Description: 
# Visualizes the results of the matrix approach to gene set compressed feature
# intepretation.

# Usage:
# The script is run in the command line
#
#     Rscript --vanilla visualize_genesets.R
#
# and takes two positional arguments as input:
#
#         --dataset             Name of the Dataset
#         --gmt_name            The name of the GMT data
#                                 (e.g. "xcell_all_entrez.gmt")
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
                                          help = "name of the gmt file"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Load arguments
dataset <- opt$dataset
gmt_name <- opt$gmt_name

# Find the name of all of the gene sets of interest
xcell_file <- file.path("..", "3.build-hetnets", "data", gmt_name)

con <- file(xcell_file, open = "r")

genesets <- c()
while (length(geneset <- readLines(con, n = 1, warn = FALSE)) > 0) {
  geneset <- unlist(strsplit(geneset, "\t"))[1]
  geneset <- substr(geneset, 1, nchar(geneset) - 2)
  genesets <- unique(c(geneset, genesets))
} 

close(con)

# Plot results for all cell_types
options(repr.plot.width = 8, repr.plot.height = 4)

for (cell_type in genesets) {
  plot_cell_type(cell_type = cell_type, dataset = dataset, show_plot = FALSE)
}
