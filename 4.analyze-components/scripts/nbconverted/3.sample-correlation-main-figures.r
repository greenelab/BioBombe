
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))

# Load helper functions
source(file.path("scripts", "util.R"))

# Create theme
correlation_theme <- theme(axis.text.x = element_text(angle = 90, size = 5),
                           plot.title = element_text(hjust = 0.5),
                           legend.text = element_text(size = 8),
                           legend.key.size = unit(0.7, 'lines'))

file <- file.path("results", "TCGA_sample_correlation_phenotype.tsv.gz")
tcga_full_df <- readr::read_tsv(file,
                                col_types = readr::cols(
                                    .default = readr::col_character(),
                                    correlation = readr::col_double(),
                                    seed = readr::col_integer(),
                                    num_comp = readr::col_integer()
                                )
                               )

print(dim(tcga_full_df))

head(tcga_full_df)

tcga_subset_df <- subset_correlations(df = tcga_full_df,
                                      cor_type = "pearson",
                                      data_type = "testing",
                                      signal_type = "signal")

print(dim(tcga_subset_df))

head(tcga_subset_df)

tcga_full_corr_gg <- plot_correlation_summary(df = tcga_subset_df)
tcga_full_corr_gg

select_cancer_types <- c("BRCA", "COAD", "LGG", "PCPG", "LAML")
mean_cor_palette <- colorRampPalette(rev(brewer.pal(9, "YlGn")))

file <- file.path("results", "TCGA_sample_correlation_phenotype_summary.tsv.gz")
tcga_summary_df <- readr::read_tsv(file,
                                   col_types = readr::cols(
                                       .default = readr::col_character(),
                                       num_comp = readr::col_integer(),
                                       mean_cor = readr::col_double(),
                                       var_cor = readr::col_double()
                                   )
                                  )

head(tcga_summary_df)

tcga_subset_summary_df <- process_summary_df(summary_df = tcga_summary_df,
                                             select_sample_types = select_cancer_types,
                                             cor_type = 'pearson',
                                             data_type = 'testing',
                                             signal_type = 'signal')

tcga_select_cancertype_gg <- plot_subset_summary(tcga_subset_summary_df, palette = mean_cor_palette)
tcga_select_cancertype_gg

tcga_capacity_gain_df <- process_capacity(tcga_summary_df, select_cancer_types)
tcga_select_cancertype_capactity_gain_gg <- plot_capacity_difference(tcga_capacity_gain_df)

tcga_select_cancertype_capactity_gain_gg

select_tissues <- c("Liver", "Pancreas", "Blood")

file <- file.path("results", "GTEX_sample_correlation_phenotype_summary.tsv.gz")
gtex_summary_df <- readr::read_tsv(file,
                                   col_types = readr::cols(
                                       .default = readr::col_character(),
                                       num_comp = readr::col_integer(),
                                       mean_cor = readr::col_double(),
                                       var_cor = readr::col_double()
                                   )
                                  )

head(gtex_summary_df)

gtex_capacity_gain_df <- process_capacity(gtex_summary_df, select_tissues)
gtex_select_tissuetype_capactity_gain_gg <- plot_capacity_difference(gtex_capacity_gain_df)

gtex_select_tissuetype_capactity_gain_gg

alg_legend <- cowplot::get_legend(gtex_select_tissuetype_capactity_gain_gg) 
cor_legend <- cowplot::get_legend(tcga_select_cancertype_gg) 

legend <- (
    cowplot::plot_grid(
        cor_legend,
        alg_legend,
        nrow = 2
    )
)

main_plot <- (
    cowplot::plot_grid(
        tcga_full_corr_gg + theme(legend.position = 'none'),
        tcga_select_cancertype_gg + theme(legend.position = "none"),
        tcga_select_cancertype_capactity_gain_gg + theme(legend.position = 'none'),
        gtex_select_tissuetype_capactity_gain_gg + theme(legend.position = 'none'),
        labels = c("A", "B", "C", "D"),
        ncol = 2,
        nrow = 2
    )
)

main_plot = cowplot::plot_grid(main_plot, legend, rel_widths = c(1, 0.15), ncol = 2)
main_plot

for(extension in c('.png', '.pdf')) {
    sup_file <- paste0("correlation_summary", extension)
    sup_file <- file.path("figures", sup_file)
    cowplot::save_plot(filename = sup_file,
                       plot = main_plot,
                       base_height = 7.5,
                       base_width = 9)
}

tcga_shuffled_subset_df <- subset_correlations(
    df = tcga_full_df,
    cor_type = "pearson",
    data_type = "testing",
    signal_type = "shuffled"
    )

print(dim(tcga_shuffled_subset_df))

head(tcga_shuffled_subset_df)

tcga_full_corr_shuffled_gg <- plot_correlation_summary(df = tcga_shuffled_subset_df, ylimits = c(-0.03, 0.09))
tcga_full_corr_shuffled_gg

tcga_subset_shuffled_summary_df <- 
    process_summary_df(summary_df = tcga_summary_df,
                       select_sample_types = select_cancer_types,
                       cor_type = 'pearson',
                       data_type = 'testing',
                       signal_type = 'shuffled')

tcga_select_shuffled_cancertype_gg <-
    plot_subset_summary(tcga_subset_shuffled_summary_df, palette = mean_cor_palette)
tcga_select_shuffled_cancertype_gg

file <- file.path("results", "GTEX_sample_correlation_phenotype.tsv.gz")
gtex_full_df <- readr::read_tsv(file,
                                col_types = readr::cols(
                                    .default = readr::col_character(),
                                    correlation = readr::col_double(),
                                    seed = readr::col_integer(),
                                    num_comp = readr::col_integer()
                                )
                               )

print(dim(gtex_full_df))

head(gtex_full_df)

gtex_subset_df <- subset_correlations(df = gtex_full_df,
                                      cor_type = "pearson",
                                      data_type = "testing",
                                      signal_type = "signal")

     
gtex_subset_shuffled_df <- subset_correlations(df = gtex_full_df,
                                               cor_type = "pearson",
                                               data_type = "testing",
                                               signal_type = "shuffled")

gtex_full_corr_gg <- plot_correlation_summary(df = gtex_subset_df,
                                              ylimits = c(-0.15, 1))
gtex_full_corr_shuffled_gg <- plot_correlation_summary(df = gtex_subset_shuffled_df,
                                                       ylimits = c(-0.01, 0.15))

gtex_full_cor_gg <- cowplot::plot_grid(
    gtex_full_corr_gg + theme(legend.position = 'none') + ylab('Real Data\nSample Correlations') + ggtitle('GTEX') + xlab(''),
    gtex_full_corr_shuffled_gg + theme(legend.position = "none") + ylab('Shuffled Data\nSample Correlations'),
    labels = c("", ""),
    ncol = 1,
    nrow = 2
)

gtex_full_cor_gg

file <- file.path("results", "TARGET_sample_correlation_phenotype.tsv.gz")
target_full_df <- readr::read_tsv(file,
                                  col_types = readr::cols(
                                      .default = readr::col_character(),
                                      correlation = readr::col_double(),
                                      seed = readr::col_integer(),
                                      num_comp = readr::col_integer()
                                  )
                                 )

print(dim(target_full_df))

head(target_full_df)

target_subset_df <- subset_correlations(df = target_full_df,
                                        cor_type = "pearson",
                                        data_type = "testing",
                                        signal_type = "signal")

     
target_subset_shuffled_df <- subset_correlations(df = target_full_df,
                                                 cor_type = "pearson",
                                                 data_type = "testing",
                                                 signal_type = "shuffled")

target_full_corr_gg <- plot_correlation_summary(df = target_subset_df,
                                               ylimits = c(-0.15, 1))
target_full_corr_shuffled_gg <- plot_correlation_summary(df = target_subset_shuffled_df,
                                                         ylimits =  c(-0.01, 0.15))

target_full_cor_gg <- cowplot::plot_grid(
    target_full_corr_gg + theme(legend.position = 'none') + ylab('') + xlab('') + ggtitle('TARGET'),
    target_full_corr_shuffled_gg + theme(legend.position = "none") + ylab(''),
    labels = c("", ""),
    ncol = 1,
    nrow = 2
)

target_full_cor_gg

alg_legend <- cowplot::get_legend(tcga_full_corr_shuffled_gg) 
cor_legend <- cowplot::get_legend(tcga_select_shuffled_cancertype_gg) 

legend <- (
    cowplot::plot_grid(
        cor_legend,
        alg_legend,
        nrow = 2
    )
)

main_plot <- (
    cowplot::plot_grid(
        tcga_full_corr_shuffled_gg + theme(legend.position = 'none') + ggtitle('TCGA') + ylab('Shuffled Data\nSample Correlations'),
        tcga_select_shuffled_cancertype_gg + theme(legend.position = "none") + ggtitle('TCGA'),
        gtex_full_cor_gg,
        target_full_cor_gg,
        labels = c("A", "B", "C", "D"),
        ncol = 2,
        nrow = 2,
        rel_heights = c(0.9, 1.1)
    )
)

main_plot = cowplot::plot_grid(main_plot, legend, rel_widths = c(1, 0.15), ncol = 2)
main_plot

for(extension in c('.png', '.pdf')) {
    sup_file <- paste0("supplemental_correlation_summary", extension)
    sup_file <- file.path("figures", sup_file)
    cowplot::save_plot(filename = sup_file,
                       plot = main_plot,
                       base_height = 8.5,
                       base_width = 8)
}
