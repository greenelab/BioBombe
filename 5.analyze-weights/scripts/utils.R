# Interpreting Compression Features
# 5.analyze-weights/scripts/utils.R
# Gregory Way 2018
#
# Usage: Import only
#
#   source(file.path("scripts", "utils.R"))


plot_cell_type <- function(cell_type, dataset = "gtex", show_plot = TRUE) {
  # Logic to plot cell-type activation across z and across algorithm
  #
  # Arguments:
  # cell_type - a string indicating the cell type of interest from XCELL
  # dataset - a string indicating the dataset to focus on
  # show_plot - boolean to indicate if the plot should be printed
  #
  # Output:
  # Saves plot to file

  cell_type_dir <- file.path("results", tolower(dataset))
  cell_type_results <- list.files(cell_type_dir, full.names = TRUE)

  # Load in the specific cell type
  cell_type_list <- list()
  for (cell_type_file in cell_type_results) {

    z_dim <- unlist(strsplit(basename(cell_type_file), "_"))[3]

    cell_df <- readr::read_tsv(
      cell_type_file,
      col_types = readr::cols(
        .default = readr::col_integer(),
        model_type = readr::col_character(),
        variable = readr::col_character(),
        value = readr::col_double(),
        z_score = readr::col_double(),
        algorithm = readr::col_character())
    )

    cell_df <- cell_df %>%
      dplyr::filter(grepl(cell_type, variable, fixed=TRUE))

    cell_type_list[[z_dim]] <- cell_df
  }

  # Combine results
  full_results_df <- dplyr::bind_rows(cell_type_list)

  # Create factors for plotting
  full_results_df$z <-
    factor(full_results_df$z,
           levels =
             sort(as.numeric(paste(unique(full_results_df$z))))
    )

  full_results_df$algorithm <-
    factor(full_results_df$algorithm,
           levels = c("pca", "ica", "nmf", "dae", "vae"))

  # Create absolute value z score for plotting
  full_results_df <- full_results_df %>%
    dplyr::mutate(abs_z_score = abs(z_score))

  # Subset to the top variables only
  top_results_df <- full_results_df %>%
    dplyr::group_by(variable, algorithm, z) %>%
    dplyr::filter(abs_z_score == max(abs_z_score))

  # Plot and save to file
  p <- ggplot(top_results_df,
              aes(x = z,
                  y = abs_z_score,
                  color = algorithm,
                  group = algorithm)) +
    geom_point(size = 0.5) +
    geom_line(lwd = 0.2) +
    facet_grid(~variable) +
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
    ylab("Absolute Value Z Score") +
    xlab("Z Dimensions") +
    theme(axis.text.x = element_text(angle = 90, size = 5),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))

  base_dir <- file.path("figures", dataset, "cell_type")
  dir.create(base_dir,
             showWarnings = FALSE,
             recursive = TRUE)
  fig_file <- file.path(base_dir, paste0("cell_type_", cell_type, ".png"))
  ggsave(plot = p, fig_file, dpi = 300, height = 3, width = 7)

  if (show_plot) {
    print(p)
  }

}