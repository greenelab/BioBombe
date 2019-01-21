# Interpreting Compression Features
# 6.analyze-weights/scripts/utils.R
# Gregory Way 2018
#
# Usage: Import only
#
#   source(file.path("scripts", "utils.R"))


get_biobombe_results <- function(gene_set_dir) {
  # Given a directory of BioBombe results, load all files and extract the top
  # scoring features across for each algorithm, dimension, and geneset
  #
  # Arguments
  # gene_set_dir - a string indicating where the results are stored
  #
  # Output:
  # A dataframe of all top results across all dimensions

  gene_set_results <- list.files(gene_set_dir, full.names = TRUE)

  # Load in the specific cell type
  gene_set_list <- list()
  for (gene_set_file in gene_set_results) {

    z_dim <- unlist(strsplit(basename(gene_set_file), "_"))[3]

    gene_set_df <- readr::read_tsv(
      gene_set_file,
      col_types = readr::cols(
        .default = readr::col_integer(),
        model_type = readr::col_character(),
        variable = readr::col_character(),
        value = readr::col_double(),
        z_score = readr::col_double(),
        algorithm = readr::col_character())
    ) %>%
      dplyr::mutate(abs_z_score = abs(z_score)) %>%
      dplyr::group_by(variable, algorithm, z) %>%
      dplyr::filter(abs_z_score == max(abs_z_score)) %>%
      dplyr::distinct(abs_z_score, algorithm, feature, z, .keep_all = TRUE)

    gene_set_list[[z_dim]] <- gene_set_df
  }
  
  # Combine results
  full_results_df <- dplyr::bind_rows(gene_set_list)

  # Create factors for plotting
  full_results_df$z <-
    factor(full_results_df$z,
           levels =
             sort(as.numeric(paste(unique(full_results_df$z))))
    )

  full_results_df$algorithm <-
    factor(full_results_df$algorithm,
           levels = c("pca", "ica", "nmf", "dae", "vae"))

  full_results_df <- full_results_df %>%
    dplyr::arrange(desc(abs_z_score)) %>%
    dplyr::ungroup()

  return(full_results_df)
}

extract_top_biobombe_results <- function(biobombe_df) {
  # Extract out the top BioBombe results for each algorithm independent of z
  #
  # Arguments:
  # biobombe_df - the data frame output of `get_biobombe_results()`
  #
  # Output:
  # Returns a processed data frame of top results

  top_results_df <- biobombe_df %>%
    dplyr::group_by(variable, algorithm) %>%
    dplyr::filter(abs_z_score == max(abs_z_score)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(relative_geneset_rank = order(abs_z_score,
                                                decreasing = TRUE)) %>%
    dplyr::arrange(z) %>%
    dplyr::distinct(model_type, variable, abs_z_score, algorithm,
                    .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(abs_z_score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(absolute_rank = order(abs_z_score, decreasing = TRUE))
  
  return(top_results_df)
}


plot_gene_set <- function(gene_set,
                          full_results_df,
                          metaedge,
                          dataset = "gtex",
                          show_plot = TRUE,
                          shuffled = FALSE,
                          return_top = FALSE,
                          return_plot = FALSE) {
  # Logic to plot cell-type activation across z and across algorithm
  #
  # Arguments:
  # gene_set - a string indicating the cell type of interest
  # full_results_df - a dataframe storing all top results of the specific
  #                   dataset and metaedge combination
  # metaedge - a string indicating the metaedge used to identify enrichment
  # dataset - a string indicating the dataset to focus on
  # show_plot - boolean to indicate if the plot should be printed
  # shuffled - boolean if the data were shuffled prior to analysis
  # return_top - boolean if the top results should be saved to a file
  # return_plot - boolean if the plot should be returned
  #
  # Output:
  # Saves plot to file

  # Subset results to the specific geneset
  top_results_df <- full_results_df %>%
    dplyr::filter(grepl(gene_set, variable, fixed = TRUE))

  # Plot and save to file
  p <- ggplot(top_results_df,
              aes(x = z,
                  y = abs_z_score,
                  color = algorithm,
                  group = algorithm)) +
    geom_point(size = 0.5) +
    geom_line(lwd = 0.2) +
    scale_color_manual(name = "Algorithm",
                       values = c("#e41a1c",
                                  "#377eb8",
                                  "#4daf4a",
                                  "#984ea3",
                                  "#ff7f00"),
                       labels = c("pca" = "PCA",
                                  "ica" = "ICA",
                                  "nmf" = "NMF",
                                  "dae" = "DAE",
                                  "vae" = "VAE")) +
    theme_bw() +
    ggtitle(gene_set) +
    ylab("Absolute Value Z Score") +
    xlab("Z Dimensions") +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size = 5),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))

  if (shuffled) {
    base_dir <- file.path("figures", dataset, 'shuffled', metaedge)
  } else {
    base_dir <- file.path("figures", dataset, 'signal', metaedge)
  }

  dir.create(base_dir,
             showWarnings = FALSE,
             recursive = TRUE)

  fig_file <- file.path(base_dir, paste0("gene_set_", gene_set, ".png"))
  ggsave(plot = p, fig_file, dpi = 300, height = 3, width = 4.5)

  if (show_plot) {
    print(p)
  }

  if (return_top) {
    return(top_results_df)
  }

  if (return_plot) {
    return(p)
  }

}
