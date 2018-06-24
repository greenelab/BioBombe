
library(ggplot2)

set.seed(123)

`%>%` <- dplyr::`%>%`

source(file.path('scripts', 'viz_utils.R'))

# Set input file names
tcga_adage_file <- file.path("results", "param_sweep_adage_TCGA_full-results.tsv")
tcga_tybalt_file <- file.path("results", "param_sweep_tybalt_TCGA_full-results.tsv")

gtex_adage_file <- file.path("results", "param_sweep_adage_GTEX_full-results.tsv")
gtex_tybalt_file <- file.path("results", "param_sweep_tybalt_GTEX_full-results.tsv")

target_adage_file <- file.path("results", "param_sweep_adage_TARGET_full-results.tsv")
target_tybalt_file <- file.path("results", "param_sweep_tybalt_TARGET_full-results.tsv")

# Set output directories
tcga_fig_dir <- file.path("figures", "tcga_results")
gtex_fig_dir <- file.path("figures", "gtex_results")
target_fig_dir <- file.path("figures", "target_results")

# Load Data
tcga_tybalt <- readParamSweep(tcga_tybalt_file, algorithm = 'tybalt')

tybalt_param_z_png <- file.path(tcga_fig_dir, "z_parameter_final_loss_tybalt_TCGA.png")
tybalt_param_z_pdf <- file.path(tcga_fig_dir, "z_parameter_final_loss_tybalt_TCGA.pdf")

p <- plotFinalLoss(tcga_tybalt$select_df, algorithm = 'Tybalt', dataset = 'TCGA')
p

ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

tybalt_one_model_png <- file.path(tcga_fig_dir, "z_parameter_tybalt_training_TCGA.png")
tybalt_one_model_pdf <- file.path(tcga_fig_dir, "z_parameter_tybalt_training_TCGA.pdf")

p <- plotOneModel(tcga_tybalt$one_model_df, algorithm = 'Tybalt', dataset = 'TCGA')
p

ggsave(tybalt_one_model_png, plot = p, height = 2.5, width = 5)
ggsave(tybalt_one_model_pdf, plot = p, height = 2.5, width = 5)

tcga_tybalt_best_params <- tcga_tybalt$best_params
tcga_tybalt_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_tybalt_params_TCGA.tsv")
readr::write_tsv(tcga_tybalt_best_params, best_param_file)

tcga_tybalt_good_training_df <- tcga_tybalt$melt_df %>%
  dplyr::filter(batch_size == 50,
                epochs == 100,
                kappa == "0.0") %>%
  dplyr::filter(
    (learning_rate == "0.002" & num_components == 5) |
      (learning_rate == "0.0015" & num_components == 25) |
      (learning_rate == "0.0015" & num_components == 50) |
      (learning_rate == "0.0015" & num_components == 75) |
      (learning_rate == "0.001" & num_components == 100) |
      (learning_rate == "0.0005" & num_components == 125)
  ) %>%
  dplyr::mutate(
      num_components =
          factor(num_components,
                 levels = sort(as.numeric(paste(unique(num_components)))))
  )

tybalt_best_model_png <- file.path(tcga_fig_dir, "z_parameter_best_model_tybalt_TCGA.png")
tybalt_best_model_pdf <- file.path(tcga_fig_dir, "z_parameter_best_model_tybalt_TCGA.pdf")

p <- plotBestModel(tcga_tybalt_good_training_df, algorithm = 'Tybalt', dataset = 'TCGA')
p

ggsave(tybalt_best_model_png, plot = p, height = 2.5, width = 4)
ggsave(tybalt_best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
tcga_adage <- readParamSweep(tcga_adage_file, algorithm = 'adage')

# Specify that the model has tied weights
adage_param_z_png <- file.path(tcga_fig_dir, "z_parameter_final_loss_adage_TCGA.png")
adage_param_z_pdf <- file.path(tcga_fig_dir, "z_parameter_final_loss_adage_TCGA.pdf")

p <- plotFinalLoss(tcga_adage$select_df, algorithm = 'ADAGE', dataset = 'TCGA')
p

ggsave(adage_param_z_png, plot = p, height = 2.5, width = 5.5)
ggsave(adage_param_z_pdf, plot = p, height = 2.5, width = 5.5)

# Several hyperparameter combinations did not converge
# This was particularly a result of the low learning rates - filter and replot
tcga_adage$select_df <- tcga_adage$select_df %>%
    dplyr::filter(end_loss < 0.01, learning_rate != 'Learn: 1e-05')

adage_sub_png <- file.path(tcga_fig_dir, "z_parameter_final_loss_remove_converge_adage_TCGA.png")
adage_sub_pdf <- file.path(tcga_fig_dir, "z_parameter_final_loss_remove_converge_adage_TCGA.pdf")

p <- plotFinalLoss(tcga_adage$select_df, algorithm = 'ADAGE', dataset = 'TCGA')
p

ggsave(adage_sub_png, plot = p, height = 2.5, width = 5.5)
ggsave(adage_sub_pdf, plot = p, height = 2.5, width = 5.5)

tcga_adage_best_params <- tcga_adage$best_params
tcga_adage_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_adage_params_TCGA.tsv")
readr::write_tsv(tcga_adage_best_params, best_param_file)

adage_tied_good_training_df <- tcga_adage$melt_df %>%
  dplyr::filter(sparsity == "0.0",
                epochs == 100,
                batch_size == 50,
                noise == "0.0") %>%
  dplyr::filter(
    (num_components == 5 & learning_rate == "0.0015") |
      (num_components == 25 & learning_rate == "0.0015") |
      (num_components == 50 & learning_rate == "0.0005") |
      (num_components == 75 & learning_rate == "0.0005") |
      (num_components == 100 & learning_rate == "0.0005") |
      (num_components == 125 & learning_rate == "0.0005")
  ) %>%
  dplyr::mutate(
      num_components =
          factor(num_components,
                 levels = sort(as.numeric(paste(unique(num_components)))))
  )

best_model_png <- file.path(tcga_fig_dir, "z_parameter_best_model_adage_TCGA.png")
best_model_pdf <- file.path(tcga_fig_dir, "z_parameter_best_model_adage_TCGA.pdf")

p <- plotBestModel(adage_tied_good_training_df, algorithm = 'ADAGE', dataset = 'TCGA')
p

ggsave(best_model_png, plot = p, height = 2.5, width = 4)
ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
gtex_tybalt <- readParamSweep(gtex_tybalt_file, algorithm = 'tybalt')

tybalt_param_z_png <- file.path(gtex_fig_dir, "z_parameter_final_loss_tybalt_GTEX.png")
tybalt_param_z_pdf <- file.path(gtex_fig_dir, "z_parameter_final_loss_tybalt_GTEX.pdf")

p <- plotFinalLoss(gtex_tybalt$select_df, algorithm = 'Tybalt', dataset = 'GTEx')
p

ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

gtex_tybalt_best_params <- gtex_tybalt$best_params
gtex_tybalt_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_tybalt_params_GTEX.tsv")
readr::write_tsv(gtex_tybalt_best_params, best_param_file)

gtex_tybalt_good_training_df <- gtex_tybalt$melt_df %>%
  dplyr::filter(epochs == 100, kappa == 0.5) %>%
  dplyr::filter(
    (learning_rate == "0.0025" & batch_size == 100 & num_components == 5) |
      (learning_rate == "0.0025" & batch_size == 100 & num_components == 25) |
      (learning_rate == "0.002" & batch_size == 100 & num_components == 50) |
      (learning_rate == "0.002" & batch_size == 50 & num_components == 75) |
      (learning_rate == "0.0015" & batch_size == 50 & num_components == 100) |
      (learning_rate == "0.0015" & batch_size == 50 & num_components == 125)
  ) %>%
  dplyr::mutate(
      num_components =
          factor(num_components,
                 levels = sort(as.numeric(paste(unique(num_components)))))
  )

best_model_png <- file.path(gtex_fig_dir, "z_parameter_best_model_tybalt_GTEX.png")
best_model_pdf <- file.path(gtex_fig_dir, "z_parameter_best_model_tybalt_GTEX.pdf")

p <- plotBestModel(gtex_tybalt_good_training_df, algorithm = 'Tybalt', dataset = 'GTEx')
p

ggsave(best_model_png, plot = p, height = 2.5, width = 4)
ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
gtex_adage <- readParamSweep(gtex_adage_file, algorithm = 'adage')

adage_param_z_png <- file.path(gtex_fig_dir, "z_parameter_final_loss_adage_GTEX.png")
adage_param_z_pdf <- file.path(gtex_fig_dir, "z_parameter_final_loss_adage_GTEX.pdf")

p <- plotFinalLoss(gtex_adage$select_df, algorithm = 'ADAGE', dataset = 'GTEx')
p

ggsave(adage_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(adage_param_z_pdf, plot = p, height = 3, width = 5.5)

# Several hyperparameter combinations did not converge
# This was particularly a result of the low learning rates - filter and replot
gtex_adage$select_df <- gtex_adage$select_df %>%
    dplyr::filter(end_loss < 0.01, learning_rate != 'Learn: 1e-05')

adage_sub_png <- file.path(gtex_fig_dir, "z_parameter_final_loss_remove_converge_adage_GTEX.png")
adage_sub_pdf <- file.path(gtex_fig_dir, "z_parameter_final_loss_remove_converge_adage_GTEX.pdf")

p <- plotFinalLoss(gtex_adage$select_df, algorithm = 'ADAGE', dataset = 'GTEx')
p

ggsave(adage_sub_png, plot = p, height = 2.5, width = 5.5)
ggsave(adage_sub_pdf, plot = p, height = 2.5, width = 5.5)

gtex_adage_best_params <- gtex_adage$best_params
gtex_adage_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_adage_params_GTEX.tsv")
readr::write_tsv(gtex_adage_best_params, best_param_file)

gtex_adage_tied_good_training_df <- gtex_adage$melt_df %>%
  dplyr::filter(sparsity == "0.0",
                epochs == 100,
                batch_size == 50) %>%
  dplyr::filter(
    (num_components == 5 & learning_rate == "0.001" & noise == "0.1") |
      (num_components == 25 & learning_rate == "0.001" & noise == "0.0") |
      (num_components == 50 & learning_rate == "0.0005" & noise == "0.0") |
      (num_components == 75 & learning_rate == "0.0005" & noise == "0.0") |
      (num_components == 100 & learning_rate == "0.0005" & noise == "0.0") |
      (num_components == 125 & learning_rate == "0.0005" & noise == "0.0")
  ) %>%
  dplyr::mutate(
      num_components =
          factor(num_components,
                 levels = sort(as.numeric(paste(unique(num_components)))))
  )

best_model_png <- file.path(gtex_fig_dir, "z_parameter_best_model_adage_GTEX.png")
best_model_pdf <- file.path(gtex_fig_dir, "z_parameter_best_model_adage_GTEX.pdf")

p <- plotBestModel(gtex_adage_tied_good_training_df, dataset = 'GTEx', algorithm = 'ADAGE')
p

ggsave(best_model_png, plot = p, height = 2.5, width = 4)
ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
target_tybalt <- readParamSweep(target_tybalt_file, algorithm = 'tybalt')

tybalt_param_z_png <- file.path(target_fig_dir, "z_parameter_final_loss_tybalt_TARGET.png")
tybalt_param_z_pdf <- file.path(target_fig_dir, "z_parameter_final_loss_tybalt_TARGET.pdf")

p <- plotFinalLoss(target_tybalt$select_df, algorithm = 'Tybalt', dataset = 'TARGET')
p

ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

target_tybalt_best_params <- target_tybalt$best_params
target_tybalt_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_tybalt_params_TARGET.tsv")
readr::write_tsv(target_tybalt_best_params, best_param_file)

target_tybalt_good_training_df <- target_tybalt$melt_df %>%
  dplyr::filter(batch_size == 25,
                epochs == 100,
                kappa == '0.5') %>%
  dplyr::filter(
    (num_components == 5 & learning_rate == "0.0015") |
      (num_components == 25 & learning_rate == "0.0015") |
      (num_components == 50 & learning_rate == "0.0015") |
      (num_components == 75 & learning_rate == "0.0015") |
      (num_components == 100 & learning_rate == "0.0015") |
      (num_components == 125 & learning_rate == "0.0005")
  ) %>%
  dplyr::mutate(
      num_components =
          factor(num_components,
                 levels = sort(as.numeric(paste(unique(num_components)))))
  )

best_model_png <- file.path(target_fig_dir, "z_parameter_best_model_tybalt_TARGET.png")
best_model_pdf <- file.path(target_fig_dir, "z_parameter_best_model_tybalt_TARGET.pdf")

p <- plotBestModel(target_tybalt_good_training_df, algorithm = 'Tybalt', dataset = 'TARGET')
p

ggsave(best_model_png, plot = p, height = 2.5, width = 4)
ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
target_adage <- readParamSweep(target_adage_file, algorithm = 'adage')

adage_param_z_png <- file.path(target_fig_dir, "z_parameter_final_loss_adage_TARGET.png")
adage_param_z_pdf <- file.path(target_fig_dir, "z_parameter_final_loss_adage_TARGET.pdf")

p <- plotFinalLoss(target_adage$select_df, algorithm = 'ADAGE', dataset = 'TARGET')
p

ggsave(adage_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(adage_param_z_pdf, plot = p, height = 3, width = 5.5)

target_adage_best_params <- target_adage$best_params
target_adage_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_adage_params_TARGET.tsv")
readr::write_tsv(target_adage_best_params, best_param_file)

target_adage_tied_good_training_df <- target_adage$melt_df %>%
  dplyr::filter(sparsity == "0.0",
                epochs == 100,
                batch_size == 50,
                noise == "0.1",
                learning_rate == "0.0005") %>%
  dplyr::mutate(
      num_components =
          factor(num_components,
                 levels = sort(as.numeric(paste(unique(num_components)))))
  )

best_model_png <- file.path(target_fig_dir, "z_parameter_best_model_adage_TARGET.png")
best_model_pdf <- file.path(target_fig_dir, "z_parameter_best_model_adage_TARGET.pdf")

p <- plotBestModel(target_adage_tied_good_training_df, algorithm = 'ADAGE', dataset = 'TARGET')
p

ggsave(best_model_png, plot = p, height = 2.5, width = 4)
ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)
