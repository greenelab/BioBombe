
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Define the dataset to compile results for
svcca_file <- file.path("results", "svcca_mean_correlation_within_z.tsv.gz")
svcca_df <- readr::read_tsv(svcca_file,
                            col_types = readr::cols(
                                .default = readr::col_character(),
                                svcca_mean_similarity = readr::col_double()))
head(svcca_df, 2)

 # Make sure factors are in order
svcca_df$z_dim <- factor(
    svcca_df$z_dim,
    levels = sort(as.numeric(paste(unique(svcca_df$z_dim))))
)

algorithms <- c("pca", "ica", "nmf", "dae", "vae")

svcca_df$algorithm_1 <- factor(
    svcca_df$algorithm_1,
    levels = algorithms
)

svcca_df$algorithm_2 <- factor(
    svcca_df$algorithm_2,
    levels = algorithms
)

options(repr.plot.width = 15, repr.plot.height = 9)

for (dataset in c("TARGET", "TCGA", "GTEX")) {
    
    for (shuffled in c(TRUE, FALSE)) {
        
        if (shuffled) {
            shuff = "shuffled"
            svcca_subset = 0
        } else {
            shuff = "signal"
            svcca_subset = 0.3
        }
        
        svcca_subset_df <- svcca_df %>%
            dplyr::filter(shuffled == !!shuff,
                          dataset == !!dataset,
                          svcca_mean_similarity > !!svcca_subset)
        
        out_figure <- file.path("figures", "svcca", paste(dataset, shuff, sep = '_'))
        plot_title <- paste0("SVCCA Mean Correlations\n", dataset, ' - ', shuff)
        
        g <- ggplot(svcca_subset_df, aes(x = z_dim, y = svcca_mean_similarity, fill = algorithm_1)) +
                geom_boxplot(outlier.size = 0.1, lwd = 0.8) +
                facet_grid(algorithm_1 ~ algorithm_2) +
                xlab("Z Dimension") +
                ylab("SVCCA Mean Similarity") +
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
                ggtitle(plot_title) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, size = 5),
                      plot.title = element_text(hjust = 0.5, size = 18),
                      legend.text = element_text(size = 12),
                      legend.key.size = unit(1, "lines"),
                      strip.text.x = element_text(size = 20),
                      strip.text.y = element_text(size = 20),
                      strip.background = element_rect(colour = "black", fill = "lightgrey"))
        
        ggsave(plot = g, filename = paste0(out_figure, ".png"), height = 8, width = 12)
        ggsave(plot = g, filename = paste0(out_figure, ".pdf"), height = 8, width = 12)
        
        print(g)
    }
}
