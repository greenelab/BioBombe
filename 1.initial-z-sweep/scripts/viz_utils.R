# Interpretation of Compressed Gene Expression Features
# Gregory Way 2018
#
# Functions to help with processing and visualizing parameter sweep results
#
# Usage: Import Only
#
#           source('scripts/viz_utils.R')

readParamSweep <- function(param_file, algorithm) {
    # Read and process parameter sweep file
    #
    # Arguments:
    # param_file - a string storing the location of a tab separated file
    # algorithm - a string indicating either 'adage' or 'tybalt'
    #
    # Output:
    # a list of length two with a full dataframe and melted across id_vars

    if (algorithm == 'adage') {
        id_vars <-  c("learning_rate", "batch_size", "epochs", "noise",
                      "train_epoch", "num_components", "sparsity", "seed")
    } else if (algorithm == 'tybalt') {
        id_vars <- c("learning_rate", "batch_size", "epochs", "kappa",
                     "train_epoch", "num_components")
    } else {
        stop('algorithm must be either "adage" or "tybalt"')
    }

    # Load full parameter sweep file
    full_df <- readr::read_tsv(param_file,
                               col_types = readr::cols(
                                   .default = readr::col_character(),
                                   train_epoch = readr::col_integer(),
                                   loss = readr::col_double(),
                                   val_loss = readr::col_double()))

    # Melt the parameter sweep results
    melt_df <- reshape2::melt(full_df,
                              id.vars = id_vars,
                              measure.vars = c("loss", "val_loss"),
                              variable.name = "loss_type",
                              value.name = "loss")

    # Select the lowest loss results across parameter combinations
    select_df <- melt_df %>% dplyr::filter(loss_type == "val_loss")
    if (algorithm == 'adage') {
        select_df <- select_df %>%
            dplyr::group_by(learning_rate, batch_size, epochs,
                            noise, sparsity, num_components)
    } else if (algorithm == 'tybalt') {
        select_df <- select_df %>%
            dplyr::group_by(learning_rate, batch_size, epochs,
                            kappa, num_components)
    }
    select_df <- select_df %>% dplyr::summarize(end_loss = tail(loss, 1))
    select_df <- recodeParamSweep(select_df, algorithm)

    # Obtain a dataframe with the best params based on final validation loss
    best_params <- select_df %>%
        dplyr::group_by(num_components) %>%
        dplyr::top_n(1, -end_loss) %>%
        dplyr::arrange(as.numeric(paste(num_components)))

    # Output training curves a single full model across dimensions
    one_model_df <- full_df %>%
      dplyr::filter(learning_rate == "0.0005",
                    batch_size == "50",
                    epochs == 100)
    one_model_df$num_components <-
      dplyr::recode_factor(one_model_df$num_components,
                           `5` = "Latent Dim: 5",
                           `25` = "Latent Dim: 25",
                           `50` = "Latent Dim: 50",
                           `75` = "Latent Dim: 75",
                           `100` = "Latent Dim: 100",
                           `125` = "Latent Dim: 125")

    return(list("full_df" = full_df, "melt_df" = melt_df,
                "select_df" = select_df, "one_model_df" = one_model_df,
                "best_params" = best_params))
}

recodeParamSweep <- function(df, algorithm) {
    # Recode several variables before plotting input
    if (algorithm == 'adage') {
        df$noise <-
            dplyr::recode_factor(df$noise,
                                 "0.0" = "Noise: 0",
                                 "0.1" = "Noise: 0.1",
                                 "0.5" = "Noise: 0.5")
        spars <- df$sparsity
        df$sparsity <- factor(spars, levels = c("0.0", "1e-06", "0.001"))

        lr <- as.numeric(paste(df$learning_rate))
        lr <- paste("Learn:", lr)
        lr_levels <- paste("Learn:",
                           unique(sort(as.numeric(gsub("Learn: ", '', lr)))))
        df$learning_rate <- factor(lr, levels = lr_levels)
    } else if (algorithm == 'tybalt') {
        # Order batch size and epoch variables
        df$batch_size <- factor(df$batch_size,
                                levels = sort(unique(df$batch_size)))
        df$epochs <- factor(df$epochs, levels = sort(unique(df$epochs)))

        # Recode batch size and learning rate variables for plotting
        df$batch_size <-
          dplyr::recode_factor(df$batch_size,
                               `25` = "Batch: 25",
                               `50` = "Batch: 50",
                               `100` = "Batch: 100",
                               `150` = "Batch: 150")

        df$learning_rate <-
          dplyr::recode_factor(df$learning_rate,
                               `0.0005` = "Learn: 0.0005",
                               `0.001` = "Learn: 0.001",
                               `0.0015` = "Learn: 0.0015",
                               `0.002` = "Learn: 0.002",
                               `0.0025` = "Learn: 0.0025")
    }
    return(df)
}

# Store base theme for plotting
base_theme <- theme(axis.text = element_text(size = rel(0.5)),
                    axis.title = element_text(size = rel(0.7)),
                    axis.text.x = element_text(angle = 45),
                    strip.text = element_text(size = rel(0.5)),
                    legend.text = element_text(size = rel(0.6)),
                    legend.title = element_text(size = rel(0.8)),
                    legend.key.height = unit(0.5, "line"))

plotFinalLoss <- function(select_df, algorithm, dataset) {
  # Plot the final validation loss at training end for each set of parameters
  #
  # Arguments:
  #   select_df - a processed dataframe of final validation loss
  #   algorithm - a string indicating the algorithm to be plotted
  #   dataset - a string indicating the name of the dataset

  # Set title
  title <- paste0(dataset, ' - ', algorithm)

  p <- ggplot(select_df, aes(x = as.numeric(paste(num_components)),
                             y = end_loss)) +
    scale_x_continuous(breaks = c(5, 25, 50, 75, 100, 125)) +
    xlab("Latent Space Dimensionality") +
    ylab("Final Validation Loss") +
    theme_bw() + base_theme
  if (algorithm == 'Tybalt') {
    p <- p +
      geom_point(aes(shape = epochs, color = kappa), size = 0.8, alpha = 0.7,
                 position = position_jitter(w = 5, h = 0)) +
      scale_color_brewer(palette = "Dark2") +
      facet_grid(batch_size ~ learning_rate)
  } else {
    if (dataset == 'TARGET') {
      p <- p +
      geom_point(aes(shape = epochs, size = sparsity, color = batch_size),
                 alpha = 0.7, position = position_jitter(w = 5, h = 0)) +
      scale_color_brewer(name = "Batch Size", palette = "Dark2") +
      scale_size_manual(values = c(0.8, 0.4), name = "Sparsity")
    } else {
      p <- p +
      geom_point(aes(shape = epochs, size = batch_size, color = sparsity),
                 alpha = 0.7, position = position_jitter(w = 5, h = 0)) +
      scale_color_brewer(name = "Sparsity", palette = "Dark2") +
      scale_size_manual(values = c(0.8, 0.4), name = "Batch Size")
    }
    p <- p + facet_grid(noise ~ learning_rate) +
    scale_shape_discrete(name = "Epochs")
  }

  p <- p + ggtitle(title)

  return(p)
}

plotOneModel <- function(one_model_df, algorithm, dataset) {
  # Plot a combination of training curves given a single model over dimensions
  #
  # Arguments:
  #   one_model_df - a processed dataframe of full training curves
  #   algorithm - a string indicating the algorithm to be plotted
  #   dataset - a string indicating the name of the dataset

  # Set title
  title <- paste0(dataset, ' - ', algorithm)

  p <- ggplot(one_model_df,
              aes(x = as.numeric(paste(train_epoch)),
                  y = val_loss)) +
    geom_line(aes(color = kappa), size = 0.2, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~ num_components) +
    ylab("Final Validation Loss") +
    xlab("Training Epochs") +
    theme_bw() +
    base_theme +
    theme(legend.key.height = unit(1, "line")) +
    ggtitle(title)
  return(p)
}

plotBestModel <- function(best_model_df, algorithm, dataset) {
  # Plot the training curves of the best models for each algorithm
  #
  # Arguments:
  #   best_model_df - full training curves of subset best models
  #   algorithm - a string indicating the algorithm to be plotted
  #   dataset - a string indicating the name of the dataset

  # Set title
  title <- paste0(dataset, ' - ', algorithm)

  p <- ggplot(best_model_df, aes(x = train_epoch, y = loss)) +
    geom_line(aes(color = num_components, linetype = loss_type), size = 0.2) +
    xlab("Training Epoch") +
    ylab("Loss") +
    scale_color_brewer(name = "Latent Dimensions", palette = "Dark2") +
    scale_linetype_manual(name = "Loss Type", values = c("solid", "dotted"),
                          labels = c("Train", "Validation")) +
    theme_bw() + base_theme +
    theme(axis.text.x = element_text(angle = 0)) +
    ggtitle(title)

  return(p)
}
