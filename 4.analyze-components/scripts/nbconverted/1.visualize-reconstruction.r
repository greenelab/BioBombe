
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Load helper functions
source(file.path("scripts", "util.R"))

# Define the dataset to compile results for
dataset <- 'TARGET'
target_recon_cost_df <- compile_reconstruction_data(dataset)
recon_file <- file.path("results", paste0("reconstruction_", dataset, ".tsv"))

# Write results to file
readr::write_tsv(target_recon_cost_df, path = recon_file)

base_dir <- file.path("figures", "TARGET")

g <- plot_reconstruction_loss(target_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Define the dataset to compile results for
target_vae_recon_cost_df <- compile_reconstruction_data(dataset, data_focus = "vae")

g <- plot_vae_training(target_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".png")),
       plot = g, height = 6, width = 4)

g

# Define the dataset to compile results for
dataset <- 'TCGA'
tcga_recon_cost_df <- compile_reconstruction_data(dataset)
recon_file <- file.path("results", paste0("reconstruction_", dataset, ".tsv"))

# Write results to file
readr::write_tsv(tcga_recon_cost_df, path = recon_file)

base_dir <- file.path("figures", "TCGA")

g <- plot_reconstruction_loss(tcga_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Define the dataset to compile results for
tcga_vae_recon_cost_df <- compile_reconstruction_data(dataset, data_focus = "vae")

g <- plot_vae_training(tcga_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".png")),
       plot = g, height = 6, width = 4)

g

# Subset to iterations that may have converged
tcga_recon_cost_df <- tcga_recon_cost_df %>% dplyr::filter(reconstruction_cost < 4000)

g <- plot_reconstruction_loss(tcga_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Subset to testing non-shuffled data
tcga_recon_cost_df <- tcga_recon_cost_df %>%
    dplyr::filter(data_type == 'testing', shuffled == 'False')

g <- plot_reconstruction_loss(tcga_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Remove shuffled data and replot
tcga_vae_recon_cost_df <- tcga_vae_recon_cost_df %>% dplyr::filter(shuffle == "False")

g <- plot_vae_training(tcga_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".png")),
       plot = g, height = 6, width = 4)

g

# Define the dataset to compile results for
dataset <- 'GTEX'
gtex_recon_cost_df <- compile_reconstruction_data(dataset)
recon_file <- file.path("results", paste0("reconstruction_", dataset, ".tsv"))

# Write results to file
readr::write_tsv(gtex_recon_cost_df, path = recon_file)

base_dir <- file.path("figures", "GTEX")

g <- plot_reconstruction_loss(gtex_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Define the dataset to compile results for
gtex_vae_recon_cost_df <- compile_reconstruction_data(dataset, data_focus = "vae")

g <- plot_vae_training(gtex_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".png")),
       plot = g, height = 6, width = 4)

g

# Subset to iterations that may have converged
gtex_recon_cost_df <- gtex_recon_cost_df %>% dplyr::filter(reconstruction_cost < 5000)

g <- plot_reconstruction_loss(gtex_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Subset to testing non-shuffled data
gtex_recon_cost_df <- gtex_recon_cost_df %>%
    dplyr::filter(data_type == 'testing', shuffled == 'False')

g <- plot_reconstruction_loss(gtex_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".png")),
       plot = g, height = 3, width = 9)

g

# Remove shuffled data and replot
gtex_vae_recon_cost_df <- gtex_vae_recon_cost_df %>% dplyr::filter(shuffle == "False")

g <- plot_vae_training(gtex_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".pdf")),
       plot = g, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".png")),
       plot = g, height = 6, width = 4)

g
