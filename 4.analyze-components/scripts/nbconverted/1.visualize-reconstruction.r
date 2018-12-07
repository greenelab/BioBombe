
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# Load helper functions
source(file.path("scripts", "util.R"))

# Define the dataset to compile results for
dataset <- 'TARGET'
target_recon_cost_df <- compile_reconstruction_data(dataset)
recon_file <- file.path("results", paste0("reconstruction_", dataset, ".tsv"))

# Write results to file
readr::write_tsv(target_recon_cost_df, path = recon_file)

base_dir <- file.path("figures", "TARGET")

target_recon_gg <- plot_reconstruction_loss(target_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".pdf")),
       plot = target_recon_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".png")),
       plot = target_recon_gg, height = 3, width = 9)

target_recon_gg

# Define the dataset to compile results for
target_vae_recon_cost_df <- compile_reconstruction_data(dataset, data_focus = "vae")

target_vae_loss_gg <- plot_vae_training(target_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".pdf")),
       plot = target_vae_loss_gg, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".png")),
       plot = target_vae_loss_gg, height = 6, width = 4)

target_vae_loss_gg

# Define the dataset to compile results for
dataset <- 'TCGA'
tcga_recon_cost_df <- compile_reconstruction_data(dataset)
recon_file <- file.path("results", paste0("reconstruction_", dataset, ".tsv"))

# Write results to file
readr::write_tsv(tcga_recon_cost_df, path = recon_file)

base_dir <- file.path("figures", "TCGA")

tcga_recon_gg <- plot_reconstruction_loss(tcga_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".pdf")),
       plot = tcga_recon_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".png")),
       plot = tcga_recon_gg, height = 3, width = 9)

tcga_recon_gg

# Define the dataset to compile results for
tcga_vae_recon_cost_df <- compile_reconstruction_data(dataset, data_focus = "vae")

tcga_vae_loss_gg <- plot_vae_training(tcga_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".pdf")),
       plot = tcga_vae_loss_gg, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".png")),
       plot = tcga_vae_loss_gg, height = 6, width = 4)

tcga_vae_loss_gg

# Subset to iterations that may have converged
tcga_recon_cost_df <- tcga_recon_cost_df %>% dplyr::filter(reconstruction_cost < 4000)

tcga_recon_filter_gg <- plot_reconstruction_loss(tcga_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".pdf")),
       plot = tcga_recon_filter_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".png")),
       plot = tcga_recon_filter_gg, height = 3, width = 9)

tcga_recon_filter_gg

# Subset to testing non-shuffled data
tcga_recon_cost_df <- tcga_recon_cost_df %>%
    dplyr::filter(data_type == 'testing', shuffled == 'False')

tcga_recon_filter_test_gg <- plot_reconstruction_loss(tcga_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".pdf")),
       plot = tcga_recon_filter_test_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".png")),
       plot = tcga_recon_filter_test_gg, height = 3, width = 9)

tcga_recon_filter_test_gg

# Remove shuffled data and replot
tcga_vae_recon_cost_df <- tcga_vae_recon_cost_df %>% dplyr::filter(shuffle == "False")

tcga_vae_loss_filter_test_gg <- plot_vae_training(tcga_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".pdf")),
       plot = tcga_vae_loss_filter_test_gg, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".png")),
       plot = tcga_vae_loss_filter_test_gg, height = 6, width = 4)

tcga_vae_loss_filter_test_gg

# Define the dataset to compile results for
dataset <- 'GTEX'
gtex_recon_cost_df <- compile_reconstruction_data(dataset)
recon_file <- file.path("results", paste0("reconstruction_", dataset, ".tsv"))

# Write results to file
readr::write_tsv(gtex_recon_cost_df, path = recon_file)

base_dir <- file.path("figures", "GTEX")

gtex_recon_gg <- plot_reconstruction_loss(gtex_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".pdf")),
       plot = gtex_recon_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_", dataset, ".png")),
       plot = gtex_recon_gg, height = 3, width = 9)

gtex_recon_gg

# Define the dataset to compile results for
gtex_vae_recon_cost_df <- compile_reconstruction_data(dataset, data_focus = "vae")

gtex_vae_loss_gg <- plot_vae_training(gtex_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".pdf")),
       plot = gtex_vae_loss_gg, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_", dataset, ".png")),
       plot = gtex_vae_loss_gg, height = 6, width = 4)

gtex_vae_loss_gg

# Subset to iterations that may have converged
gtex_recon_cost_df <- gtex_recon_cost_df %>% dplyr::filter(reconstruction_cost < 5000)

gtex_recon_filter_gg <- plot_reconstruction_loss(gtex_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".pdf")),
       plot = gtex_recon_filter_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_", dataset, ".png")),
       plot = gtex_recon_filter_gg, height = 3, width = 9)

gtex_recon_filter_gg

# Subset to testing non-shuffled data
gtex_recon_cost_df <- gtex_recon_cost_df %>%
    dplyr::filter(data_type == 'testing', shuffled == 'False')

gtex_recon_filter_test_gg <- plot_reconstruction_loss(gtex_recon_cost_df)

ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".pdf")),
       plot = gtex_recon_filter_test_gg, dpi = 500, height = 3, width = 9)
ggsave(file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", dataset, ".png")),
       plot = gtex_recon_filter_test_gg, height = 3, width = 9)

gtex_recon_filter_test_gg

# Remove shuffled data and replot
gtex_vae_recon_cost_df <- gtex_vae_recon_cost_df %>% dplyr::filter(shuffle == "False")

gtex_vae_loss_filter_test_gg <- plot_vae_training(gtex_vae_recon_cost_df)

ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".pdf")),
       plot = gtex_vae_loss_filter_test_gg, dpi = 500, height = 6, width = 4)
ggsave(file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", dataset, ".png")),
       plot = gtex_vae_loss_filter_test_gg, height = 6, width = 4)

gtex_vae_loss_filter_test_gg

legend <- get_legend(target_recon_gg) 

main_plot <- (
    cowplot::plot_grid(
        gtex_recon_filter_test_gg + ggtitle('GTEX') + xlab('') +
            theme(plot.margin = margin(t = 0.5, r = 0.2, b = 0, l = 0.4),
                  legend.position = "none"),
        tcga_recon_filter_test_gg + ggtitle('TCGA') + xlab('') +
            theme(plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.4),
                  legend.position = "none"),
        target_recon_gg + ggtitle('TARGET') +
            theme(plot.margin = margin(t = 0, r = 0.2, b = 0.3, l = 0.4),
                  legend.position = "none"),
        labels = c("A", "B", "C"),
        ncol = 1,
        nrow = 3
    )
)

main_plot = cowplot::plot_grid(main_plot, legend, rel_widths = c(1, 0.15), ncol = 2)
main_plot

for(extension in c('.png', '.pdf')) {
    sup_file <- paste0("reconstruction_summary", extension)
    sup_file <- file.path("figures", sup_file)
    cowplot::save_plot(filename = sup_file,
                       plot = main_plot,
                       base_height = 6.5,
                       base_width = 8)
}
