# Interpretation of Compression Models
#
# Gregory Way 2018
#
# This script stores several helper functions for plotting the ensemble z sweep
#
# Usage:
# source("scripts/util.R")

compile_reconstruction_data <- function(dataset_name, data_focus = "all") {
  # Crawl through given folder structure to obtain dataset specific
  # reconstruction results
  #
  # Arguments:
  # dataset_name - the name of the dataset of interest
  # data_focus - which data to identify and process (depends on plotting)
  #              (default = "all", can also use "vae")
  # Output:
  # A long dataframe storing all reconstruction results
  # Determine each component directory
  component_dir_list <- list()
  results_types <- c("signal", "shuffled")

  if (data_focus == "all") {
    search_pattern <- "reconstruction.tsv"
  } else {
    search_pattern <- "tybalt_training_hist.tsv"
  }

  for (result_type in results_types) {
    if (result_type == "signal") {
      results_suffix <- "_results"
    } else if (result_type == "shuffled") {
      results_suffix <- "_shuffled_results"
    }
    results_dir <- file.path("..",
                             "2.ensemble-z-analysis",
                             "results",
                             paste0(dataset_name, results_suffix),
                             "ensemble_z_results")
    component_dirs <- list.dirs(results_dir, recursive = FALSE)
    component_dir_list[[result_type]] <- component_dirs
  }
  
  # Load and process reconstruction costs across folders
  component_list <- list()
  for (component_dir in names(component_dir_list)) {
    # Retreive the appropriate directory (shuffled or signal)
    component_dirs <- component_dir_list[[component_dir]]
    
    # Loop over the results directory and store results
    for (component_path in component_dirs) {
      # Extract the bottleneck dimensionality estimate
      num_component <- unlist(strsplit(basename(component_path), "_"))[1]
      
      # Find the reconstruction cost file
      component_files <- list.files(path = component_path,
                                    full.names = TRUE,
                                    pattern = search_pattern)
      # Load and process the reconstruction file - convert to wide format
      if (data_focus == "all") {
        component_df <- readr::read_tsv(
          component_files,
          col_types = readr::cols(
            .default = readr::col_double(),
            seed = readr::col_integer(),
            shuffled = readr::col_character(),
            data_type = readr::col_character()
            )
          ) %>%
          dplyr::mutate(num_comp = num_component,
                        dataset_id = dataset_name) %>%
          reshape2::melt(
            id.vars = c("dataset_id", "num_comp", "data_type", "shuffled",
                        "seed"),
            variable.name = "algorithm",
            value.name = "reconstruction_cost"
          )
      } else {
        component_df <- readr::read_tsv(
          component_files,
          col_types = readr::cols(
            .default = readr::col_double(),
            seed = readr::col_integer(),
            shuffle = readr::col_character(),
            epoch = readr::col_integer()
          )
        ) %>%
          dplyr::mutate(num_comp = num_component,
                        dataset_id = dataset_name) %>%
          dplyr::group_by(seed) %>%
          dplyr::filter(row_number() == n()) %>%
          reshape2::melt(
            id.vars = c("epoch", "loss", "val_loss", "seed", "shuffle",
                        "num_comp", "dataset_id"),
            measure.vars = c("recon", "kl"),
            variable.name = "loss_type",
            value.name = "partial_loss")
        
        levels(component_df$loss_type) <- c("Reconstruction", "KL Divergence")
      }

      # Store in list for future concatenation
      list_name <- paste0(num_component, component_dir)
      component_list[[list_name]] <- component_df
    }
  }
  
  # Combine long dataframe of reconstruction results
  reconstruction_cost_df <- dplyr::bind_rows(component_list)
  
  # Make sure factors are in order
  reconstruction_cost_df$num_comp <-
    factor(reconstruction_cost_df$num_comp,
           levels =
             sort(as.numeric(paste(unique(reconstruction_cost_df$num_comp))))
           )
  
  if (data_focus == "all") {
    reconstruction_cost_df$algorithm <-
      factor(reconstruction_cost_df$algorithm,
             levels = c("pca", "ica", "nmf", "dae", "vae"))

    reconstruction_cost_df$reconstruction_cost <-
      as.numeric(paste(reconstruction_cost_df$reconstruction_cost))
  } else {
    reconstruction_cost_df$partial_loss <-
      as.numeric(paste(reconstruction_cost_df$partial_loss))
  }

  return(reconstruction_cost_df)
}

plot_reconstruction_loss <- function(data_df) {
  # Crawl through given folder structure to obtain dataset specific
  # reconstruction results
  #
  # Arguments:
  # data_df - the dataframe to be plotted
  #
  # Output:
  # A ggplot object to be saved and viewed
  
  options(repr.plot.width = 9, repr.plot.height = 4)
  p <- ggplot(data = data_df, aes(x = num_comp, y = reconstruction_cost)) +
    geom_point(aes(color = algorithm, shape = data_type, alpha = shuffled),
               size = 0.5) +
    scale_alpha_manual(values = c(0.75, 0.15),
                       labels = c("Real", "Shuffled"),
                       name = "Data") +
    scale_shape_manual(name = "Data Type",
                       values = c(16, 17),
                       labels = c("Testing", "Training")) +
    scale_color_manual(name = "Algorithm",
                       values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                                  "#ff7f00"),
                       labels = c("PCA", "ICA", "NMF", "DAE", "VAE")) +
    facet_grid(~ algorithm) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Reconstruction Cost") +
    ggtitle(paste(dataset, "Reconstruction Cost")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 4),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))
  
  return(p)
}

plot_vae_training <- function(data_df) {
  # Crawl through given folder structure to obtain dataset specific
  # VAE reconstruction results
  #
  # Arguments:
  # data_df - the dataframe to be plotted
  #
  # Output:
  # A ggplot object to be saved and viewed

  options(repr.plot.width = 4, repr.plot.height = 6)
  p <- ggplot(data = data_df,
              aes(x = num_comp, y = partial_loss)) +
    geom_boxplot(aes(color = loss_type, fill = shuffle),
                 outlier.size = 0.1,
                 lwd = 0.3) +
    scale_fill_manual(values = c("black", "grey"),
                      labels = c("Real", "Shuffled"),
                      name = "Data") +
    scale_color_manual(name = "Loss Type",
                       values = c("#e41a1c", "#377eb8"),
                       labels = c("Reconstruction", "KL Divergence")) +
    facet_wrap(~loss_type, scales = "free", nrow = 2) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Reconstruction Cost") +
    ggtitle(paste(dataset, "VAE Loss")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 4),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))
  
  return(p)
}
