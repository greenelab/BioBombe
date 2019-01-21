# Interpretation of Compression Models
#
# Gregory Way 2018
#
# This script stores several helper functions for plotting results of the
# ensemble z sweep
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
    
    reconstruction_cost_df$algorithm <- 
      reconstruction_cost_df$algorithm %>% dplyr::recode_factor("pca" = "PCA",
                                                                "ica" = "ICA",
                                                                "nmf" = "NMF",
                                                                "dae" = "DAE",
                                                                "vae" = "VAE")

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
  p <- ggplot(data = data_df,
              aes(x = num_comp,
                  y = reconstruction_cost)) +
    geom_point(aes(color = algorithm,
                   shape = data_type,
                   alpha = shuffled),
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
                       labels = c("pca" = "PCA",
                                  "ica" = "ICA",
                                  "nmf" = "NMF",
                                  "dae" = "DAE",
                                  "vae" = "VAE")) +
    facet_grid(~ algorithm) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Reconstruction Cost") +
    ggtitle(paste(dataset, "Reconstruction Cost")) +
    theme_bw() +
    theme(axis.title = element_text(size = 7.5),
          axis.text.y = element_text(size = 6.5),
          axis.text.x = element_text(angle = 90, size = 3.2),
          plot.title = element_text(size = 8.5, hjust = 0.5),
          legend.title = element_text(size = 7.5),
          legend.text = element_text(size = 6.5),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"),
          legend.key.size = unit(0.7, "lines"))
  
  return(p)
}

plot_vae_training <- function(data_df) {
  # Plot VAE reconstruction and KL divergence losses
  #
  # Arguments:
  # data_df - the dataframe to be plotted
  #
  # Output:
  # A ggplot object to be saved and viewed

  options(repr.plot.width = 4, repr.plot.height = 6)
  p <- ggplot(data = data_df,
              aes(x = num_comp, y = partial_loss)) +
    geom_boxplot(aes(color = shuffle),
                 outlier.size = 0.1,
                 lwd = 0.3) +
    scale_color_manual(values = c("#e41a1c", "#377eb8"),
                      labels = c("Real", "Shuffled"),
                      name = "Data") +
    facet_wrap(~ loss_type, scales = "free", nrow = 2) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Penalty") +
    ggtitle(paste(dataset, "VAE Loss")) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 90, size = 4.3),
          axis.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 8),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"),
          strip.text = element_text(size = 6),
          legend.key.size = unit(0.7, "lines"))
  
  return(p)
}

compile_sample_correlation <- function(dataset_name) {
  # Crawl through given folder structure to obtain dataset specific sample
  # correlation results
  #
  # Arguments:
  # dataset_name - the name of the dataset of interest
  #
  # Output:
  # A long dataframe storing all sample correlation results
  # Determine each component directory
  component_dir_list <- list()
  results_types <- c("signal", "shuffled")
  
  for (result_type in results_types) {
    if (result_type == "signal") {
      results_suffix <- "_results"
    } else if (result_type == "shuffled") {
      results_suffix <- "_shuffled_results"
    }
    results_dir <- file.path("..", "2.ensemble-z-analysis", "results",
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
                                    pattern = "_sample_corr.tsv.gz")
      
      # Load and process the reconstruction file - convert to wide format
      component_df <- readr::read_tsv(
        component_files,
        col_types = readr::cols(
          .default = readr::col_character(),
          seed = readr::col_integer(),
          correlation = readr::col_double()
        )
      ) %>%
        dplyr::mutate(num_comp = num_component,
                      dataset_id = dataset_name,
                      shuffled = component_dir)
      
      # Store in list for future concatenation
      list_name <- paste0(num_component, component_dir)
      component_list[[list_name]] <- component_df
    }
  }
  
  # Combine long dataframe of reconstruction results
  sample_corr_df <- dplyr::bind_rows(component_list)
  
  # Make sure factors are in order
  sample_corr_df$num_comp <-
    factor(sample_corr_df$num_comp,
           levels = sort(as.numeric(paste(unique(sample_corr_df$num_comp)))))
  
  sample_corr_df$algorithm <-
    factor(sample_corr_df$algorithm,
           levels = c("pca", "ica", "nmf", "dae", "vae"))
  
  sample_corr_df$correlation <-
    as.numeric(paste(sample_corr_df$correlation))
  
  sample_corr_df$data <-
    factor(sample_corr_df$data,
           levels = c("training", "testing"))
  
  return(sample_corr_df)
}

plot_sample_correlation <- function(data_df,
                                    dataset_name,
                                    data_type,
                                    correlation_type,
                                    use_theme,
                                    return_figures = FALSE) {
  # Helper function to plot sample correlation observations across bottleneck
  # dimensions, algorithms, sample types, and training iterations
  #
  # Arguments:
  # data_df - a dataframe with dataset specific results of sample correlations
  # dataset_name - a string indicating the name of the dataset of interest
  # data_type - either "signal" or "shuffled" (if the data was shuffled before)
  # correlation_type - either "pearson" or "spearman"
  # use_theme - a ggplot2 theme option
  # return_figures - boolean indicating if the figures should be output to list
  #
  # Output:
  # The function will save all plots to file, and, optionally, will output a
  # list storing all ggplots

  return_plot_list <- list()

  # Create a directory of where to save the figures
  sample_cor_dir <- file.path("figures", dataset_name, "sample-correlation")
  dir.create(sample_cor_dir,
             showWarnings = FALSE,
             recursive = TRUE)

  # Subset the data
  data_df <- data_df %>%
    dplyr::filter(shuffled == !!data_type,
                  cor_type == !!correlation_type) %>%
    dplyr::group_by(id, data, num_comp, algorithm, sample_type) %>%
    dplyr::summarise(median_corr = median(correlation, na.rm = TRUE))

  # Define a specific plot name
  plot_name <- paste(dataset_name, data_type, correlation_type, sep = "_")

  # Plot aggregate data
  return_plot_list[["agg_plot"]] <-
    ggplot(data = data_df,
           aes(x = num_comp, y = median_corr)) +
    geom_boxplot(aes(color = algorithm),
                 size = 0.1,
                 outlier.size = 0.05) +
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
    facet_grid(data ~ algorithm) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Sample Correlation") +
    ggtitle(paste0(dataset_name, " Sample Correlation (", data_type, " ",
                   correlation_type, ")")) +
    theme_bw() +
    use_theme

  # Save Figure
  figure_base <- file.path(sample_cor_dir,
                           paste0("full_sample-correlation_", plot_name))
  ggsave(plot = return_plot_list[["agg_plot"]],
         paste0(figure_base, ".pdf"),
         dpi = 500,
         height = 7,
         width = 10)
  ggsave(plot = return_plot_list[["agg_plot"]],
         paste0(figure_base, ".png"),
         height = 7,
         width = 10)

  # Plot full picture of cancer-type specific results
  return_plot_list[["cancer_type_full_plot"]] <-
    ggplot(data = data_df,
           aes(x = num_comp, y = median_corr)) +
    geom_boxplot(aes(color = algorithm, linetype = data),
                 size = 0.1,
                 outlier.size = 0.05) +
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
    scale_linetype_manual(name = "Data Split",
                          values = c("solid", "dashed"),
                          labels = c("training" = "Training",
                                     "testing" = "Testing")) +
    facet_grid(sample_type ~ algorithm) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Sample Correlation") +
    ggtitle(paste0(dataset_name, " Sample Correlation (", data_type, " ",
                   correlation_type, ")")) +
    theme_bw() +
    use_theme

  # Save Figure
  if (dataset_name %in% c("TCGA", "GTEX")) {
    plot_height <- 20
  } else {
    plot_height <- 7
  }

  figure_base <- file.path(sample_cor_dir,
                           paste0("sample-type_sample-correlation_", plot_name))
  ggsave(plot = return_plot_list[["cancer_type_full_plot"]],
         paste0(figure_base, ".pdf"),
         dpi = 500,
         height = plot_height,
         width = 10)
  ggsave(plot = return_plot_list[["cancer_type_full_plot"]],
         paste0(figure_base, ".png"),
         height = plot_height,
         width = 10)

  # Now, plot individual sample-types
  sample_type_specific_list <- list()

  # Create a directory of where to save the figures
  sample_type_cor_dir <- file.path("figures", dataset_name,
                                   "sample-correlation", "sample-type")
  dir.create(sample_type_cor_dir,
             showWarnings = FALSE,
             recursive = TRUE)

  for (sample_type in unique(data_df$sample_type)) {
    # Subset the sample correlation data frame to specific data and correlation
    # types
    sample_signal_data_df <- data_df %>%
      dplyr::filter(sample_type == !!sample_type)

    sample_type_specific_list[[sample_type]] <-
      ggplot(data = sample_signal_data_df,
             aes(x = num_comp, y = median_corr)) +
      geom_boxplot(aes(color = algorithm), size = 0.1, outlier.size = 0.05) +
      facet_wrap( ~ data, nrow = 2) +
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
      xlab("Latent Space Dimensions (z)") +
      ylab("Sample Correlation") +
      ggtitle(paste0(dataset_name, " Sample Correlation (", sample_type, ")\n(",
                     data_type, " ", correlation_type, ")")) +
      theme_bw() +
      use_theme

    figure_base <- file.path(sample_type_cor_dir,
                             paste0("sample-correlation_",
                                    sample_type, "_", plot_name))
    ggsave(plot = sample_type_specific_list[[sample_type]],
           paste0(figure_base, ".pdf"),
           dpi = 500,
           height = 5,
           width = 6)
    ggsave(plot = sample_type_specific_list[[sample_type]],
           paste0(figure_base, ".png"),
           height = 5,
           width = 6)
  }

  if (return_figures) {
    return_plot_list[["individual_cancertypes"]] <- sample_type_specific_list
    return(return_plot_list)
  }
}

subset_correlations <- function(df, cor_type, data_type, signal_type) {
  # Given a correlations dataframe, subset based on correlation, data
  # and signal types
  # 
  # Arguments:
  # df - the dataframe of interest to subset
  # cor_type - a string of correlation to subset ("pearson" or "spearman")
  # data_type - a string of type of data used ("training" or "testing")
  # signal_type - a string of the signal type ("signal" or "shuffled")
  #
  # Output:
  # a subset dataframe with correlations summarized across all 5 seeds (median)
  
  subset_df <- df %>%
    dplyr::filter(cor_type == !!cor_type,
                  data == !!data_type,
                  shuffled == !!signal_type) %>%
    dplyr::group_by(id, data, num_comp, algorithm, sample_type) %>%
    dplyr::summarise(median_corr = median(correlation, na.rm = TRUE))
  
  subset_df$num_comp <-
    factor(subset_df$num_comp,
           levels = sort(as.numeric(paste(unique(subset_df$num_comp)))))
  
  subset_df$algorithm <-
    factor(subset_df$algorithm,
           levels = c("pca", "ica", "nmf", "dae", "vae"))
  
  subset_df$correlation <-
    as.numeric(paste(subset_df$median_corr))
  
  return(subset_df)
}

plot_correlation_summary <- function(df,
                                     cor_type = "Pearson",
                                     ylimits = c(0, 1)) {
  # Plot a full distribution of correlation summary across dimensions
  #
  # Arguments:
  # df - a summary dataframe (output from `subset_correlations`)
  # cor_type - a string indicating which correlation type was used (for title)
  # ylimits - a vector of length two indicating where limits should be drawn
  #
  # Output:
  # a ggplot2 object plot describing sample correlations across models
  
  full_corr_gg <- 
    ggplot(data = df, aes(x = num_comp,
                          y = median_corr)) +
    geom_boxplot(aes(fill = algorithm),
                 size = 0.1,
                 outlier.size = 0.03,
                 outlier.color = "lightgrey") +
    scale_fill_manual(name = "Algorithm",
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
    xlab("Latent Space Dimensions (z)") +
    ylab(paste0("Sample Correlation (", cor_type, ")")) +
    ylim(ylimits) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 7),
          axis.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 11),
          legend.text = element_text(size = 7),
          legend.key.size = unit(1, "lines"),
          strip.text.y = element_text(size = 10))

  return(full_corr_gg)
}

process_capacity <- function(summary_df, select_sample_types) {
  # Given a dataset summary dataframe, calculate the difference in correlation
  # values as the model capacity increases
  #
  # Arguments:
  # summary_df - dataset specific dataframe of sample correlations
  # select_samples - a character vector of samples of choice to subset
  #
  # Output:
  # A dataframe storing how correlation increases with capacity + 1
  
  capacity_gain_df <- summary_df %>%
    dplyr::group_by(algorithm, sample_type, cor_type, shuffled, data) %>%
    dplyr::arrange(num_comp, .by_group = TRUE) %>%
    dplyr::mutate(capacity_diff = 
                    mean_cor - dplyr::lag(mean_cor, 
                                          default = dplyr::first(mean_cor)))
  
  capacity_gain_df <- capacity_gain_df %>%
    dplyr::filter(cor_type == "pearson", shuffled == "signal") %>%
    dplyr::filter(sample_type %in% select_sample_types)
  
  capacity_gain_df$sample_type <- (
    factor(capacity_gain_df$sample_type, levels = select_sample_types)
  )
  
  capacity_gain_df$num_comp <-
    factor(capacity_gain_df$num_comp,
           levels = as.character(
             sort(as.numeric(paste(unique(capacity_gain_df$num_comp))))
             ))
  
  capacity_gain_df$algorithm <- 
    capacity_gain_df$algorithm %>%
    dplyr::recode("pca" = "PCA",
                  "ica" = "ICA",
                  "nmf" = "NMF",
                  "dae" = "DAE",
                  "vae" = "VAE")
  
  capacity_gain_df$algorithm <-
    factor(capacity_gain_df$algorithm,
           levels = c("PCA", "ICA", "NMF", "DAE", "VAE"))
  
  capacity_gain_df <- capacity_gain_df %>%
    dplyr::mutate(group_var = paste0(algorithm, "_", data))

    factor(capacity_gain_df$algorithm,
           levels = c("PCA", "ICA", "NMF", "DAE", "VAE"))
  
  return(capacity_gain_df)
}

plot_capacity_difference <- function(capacity_df) {
  # Given a capacity dataframe, visualize sample increases with capacity
  #
  # Arguments:
  # capacity_df - dataframe storing specific sample correlation increases
  #
  # Output:
  # A ggplot object plot
  
  gg <- ggplot(capacity_df,
         aes(x = num_comp,
             y = capacity_diff,
             color = algorithm,
             linetype = data,
             group = group_var)) +
    stat_summary(fun.y = mean,
                 geom = "line",
                 size = 0.25,
                 aes(group = group_var,
                     color = algorithm,
                     linetype = data)) +
    facet_grid(sample_type ~ .) +
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
    scale_linetype_manual(name = "Data Split",
                          values = c("solid", "dotted"),
                          labels = c("training" = "Training",
                                     "testing" = "Testing")) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90,
                                     size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 9),
          plot.title = element_text(hjust = 0.5,
                                    size = 11),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size = unit(1, "lines"),
          strip.text.y = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    labs(x = "Latent Space Dimensions (z)",
         y = expression(paste("Correlation Gain ",
                              '(z'['i'], ' - ', 'z'['i - 1'], ')'))) +
    scale_x_discrete(expand = c(0, 0))
  
  return(gg)
}

process_summary_df <- function(summary_df, select_sample_types, cor_type,
                               data_type, signal_type) {
  # Given a dataset summary df, subset specific model types and sample types
  # values as the model capacity increases
  #
  # Arguments:
  # summary_df - dataset specific dataframe of sample correlations
  # select_samples - a character vector of samples of choice to subset
  # cor_type - a string of correlation to subset ("pearson" or "spearman")
  # data_type - a string of type of data used ("training" or "testing")
  # signal_type - a string of the signal type ("signal" or "shuffled")
  #
  # Output:
  # Subset data ready for plotting
  
  subset_summary_df <- summary_df %>%
    dplyr::filter(cor_type == !!cor_type,
                  shuffled == !!signal_type,
                  data == !!data_type)
  
  subset_summary_df$num_comp <-
    factor(subset_summary_df$num_comp,
           levels = as.character(sort(as.numeric(
             paste(unique(subset_summary_df$num_comp)))
             )))
  
  subset_summary_df$algorithm <- subset_summary_df$algorithm %>%
    dplyr::recode("pca" = "PCA",
                  "ica" = "ICA",
                  "nmf" = "NMF",
                  "dae" = "DAE",
                  "vae" = "VAE")
  
  subset_summary_df$algorithm <-
    factor(subset_summary_df$algorithm,
           levels = rev(c("PCA", "ICA", "NMF", "DAE", "VAE")))
  
  subset_summary_df <- subset_summary_df %>%
    dplyr::filter(sample_type %in% select_sample_types)
  
  subset_summary_df$sample_type <- factor(subset_summary_df$sample_type,
                                          levels = select_sample_types)
  
  return(subset_summary_df)
}

plot_subset_summary <- function(subset_summary_df, palette) {
  # Plot heatmap of summarized sample correlations
  #
  # Arguments:
  # subset_summary_df - a summarized and subset dataframe ready for plotting
  #
  # Output:
  # a ggplot object of median correlations across algorithms and dimensionality
  
  select_sampletype_gg <-
    ggplot(subset_summary_df, aes(x = num_comp,
                                  y = algorithm)) +
    geom_tile(aes(fill = mean_cor),
              colour = "white") +
    scale_fill_gradientn(name = "Mean Pearson\nCorrelation",
                         colours = palette(100),
                         values = scales::rescale(c(1, 0.85, 0.6))) +
    facet_grid(sample_type ~ .) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90,
                                     size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 9),
          plot.title = element_text(hjust = 0.5,
                                    size = 11),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size = unit(1, "lines"),
          strip.text.y = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    xlab("Latent Space Dimensions (z)") +
    ylab("Algorithm") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
  
  return(select_sampletype_gg)
}

save_png_pdf <- function(p, path_prefix, height, width) {
  # Save a pdf and png of a given plot and filename
  #
  # Arguments:
  # p - a ggplot object to save
  # path_prefix - directory and file name of where to save figures
  # height - integer indicating the height of the plot (in "mm")
  # width - integer indicating the width of the plot (in "mm")
  #
  # Output:
  # Will save a pdf and png of the plot of interest
  
  for (extension in c(".png", ".pdf")) {
    full_path <- paste0(path_prefix, extension)
    cowplot::save_plot(filename = full_path,
                       plot = p,
                       base_height = height,
                       base_width = width,
                       unit = "mm")
  }
}
