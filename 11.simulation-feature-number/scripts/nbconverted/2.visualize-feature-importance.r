suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

data_file <- file.path("results", "compression_simulation_results.tsv")
data_df <- readr::read_tsv(data_file, col_types = readr::cols()) %>%
    reshape2::melt(id.vars = c("k", "algorithm", "compressed_feature"),
                   variable.name = "feature",
                   value.name = "feature_importance") %>%
    dplyr::mutate(abs_feature_importance = abs(feature_importance)) %>%
    dplyr::filter(feature %in% c("feature_1", "feature_2", "feature_3",
                                 "feature_5", "feature_6", "feature_7")) %>%
    dplyr::group_by(k, algorithm, feature) %>%
    dplyr::top_n(n = 1, wt = abs_feature_importance) %>%
    tidyr::separate(compressed_feature, into = c("alg", "feature_num"), sep = "_", remove = FALSE)

data_df <- data_df %>%
    dplyr::mutate(feature_num_recode = paste(as.numeric(paste(feature_num)) + 1))

data_df$feature_num <- factor(paste(data_df$feature_num),
                              levels = c("0", "1", "2", "3", "4", "5"))

data_df$algorithm <- factor(data_df$algorithm, levels = c("PCA", "ICA", "NMF", "DAE", "VAE"))

data_df$feature_block <- dplyr::recode(
    data_df$feature,
    feature_1 = "feature_block_1",
    feature_2 = "feature_block_1",
    feature_3 = "feature_block_1",
    feature_5 = "feature_block_2",
    feature_6 = "feature_block_2",
    feature_7 = "feature_block_2",
)

print(dim(data_df))
head(data_df, 10)

feature_colors = c(
    "feature_1" = "#9B4A34",
    "feature_2" = "#D59673",
    "feature_3" = "#DC9C2A",
    "feature_5" = "#291F46",
    "feature_6" = "#373443",
    "feature_7" = "#4D5575"
)

feature_labels = c(
    "feature_1" = "1",
    "feature_2" = "2",
    "feature_3" = "3",
    "feature_5" = "5",
    "feature_6" = "6",
    "feature_7" = "7"
)


feature_block_colors = c(
    "feature_block_1" = "#C77A39",
    "feature_block_2" = "#2A2930"
)

feature_block_labels = c(
    "feature_block_1" = "Block 1",
    "feature_block_2" = "Block 2"
)


simulation_theme <- theme_bw() + 
  theme(axis.text = element_text(size = 5.5),
        axis.title = element_text(size = 6.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(1, "lines"),
        strip.background = element_rect(colour = "black",
                                        fill = "#fdfff4"),
        strip.text.x = element_text(size = 5,
                                    margin = margin(t = 2,
                                                    b = 1.5,
                                                    l = 0,
                                                    r = 0)))

append_k <- function(string) paste0("k=", string)

plot_gg <- ggplot(data_df,
       aes(x = feature_num_recode, y = feature_importance)) +
    geom_bar(stat="identity", aes(fill=feature), position = position_dodge()) +
    facet_wrap(k~algorithm,
               scales = "free",
               ncol = 5,
               labeller = labeller(k = as_labeller(append_k))) +
    scale_fill_manual(values = feature_colors, labels = feature_labels, name = "Feature") +
    xlab("Compressed Feature Number (Top Scoring Feature)") +
    ylab("Feature Importance (Weight)") +
    simulation_theme

plot_gg

for(extension in c('.png', '.pdf')) {
    fig_file <- paste0("simulated_feature_number", extension)
    fig_file <- file.path("figures", fig_file)
    ggsave(filename = fig_file,
           plot = plot_gg,
           height = 150,
           width = 125,
           units = "mm",
           dpi = 400)
}

data_count_df <- data_df %>%
    dplyr::filter(k == 6) %>%
    dplyr::group_by(algorithm, feature_block, feature_num_recode) %>%
    dplyr::count() %>%
    tidyr::drop_na()

data_count_df

plot_count_gg <- ggplot(data_count_df,
       aes(x = feature_num_recode, y = n, fill = feature_block)) +
    geom_bar(stat="identity", position = position_dodge()) +
    facet_wrap(~algorithm,
               scales = "fixed",
               ncol = 5) +
    scale_fill_manual(values = feature_block_colors, labels = feature_block_labels, name = "") +
    xlab("Compressed Feature Number (Top Scoring Feature)") +
    ylab("Number of Times\nFeature Block Grouped Together") +
    simulation_theme

plot_count_gg

for(extension in c('.png', '.pdf')) {
    fig_file <- paste0("simulated_feature_number_feature_block_summary", extension)
    fig_file <- file.path("figures", fig_file)
    ggsave(filename = fig_file,
           plot = plot_count_gg,
           height = 50,
           width = 125,
           units = "mm",
           dpi = 400)
}
