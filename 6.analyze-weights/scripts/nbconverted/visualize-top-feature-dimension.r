
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

algorithms <- c("PCA", "ICA", "NMF", "DAE", "VAE")

collections <- c("GpXCELL", "GpC4CM",
                 "GpH", "GpC2CPG",
                 "GpC2CPREACTOME", "GpC3TFT")

base_dir <- file.path("results", "top_features")

all_top_list <- list()
list_idx <- 1
for (top_file in list.files(base_dir)) {
    if (grepl("top_biobombe_scores", top_file)) {
        dataset <- strsplit(top_file, '_')[[1]][1]
        collection <- strsplit(top_file, '_')[[1]][2]
        
        top_df <- readr::read_tsv(file.path(base_dir, top_file),
                                col_types = readr::cols(
                                    .default = readr::col_character(),
                                    value = readr::col_double(),
                                    z_score = readr::col_double(),
                                    feature = readr::col_integer(),
                                    z = readr::col_integer(),
                                    abs_z_score = readr::col_double(),
                                    relative_geneset_rank = readr::col_integer(),
                                    absolute_rank = readr::col_integer())) %>%
        dplyr::mutate(collection = collection,
                      dataset = dataset)
        
        all_top_list[[list_idx]] <- top_df
        list_idx <- list_idx + 1
    }
}

all_top_df <- dplyr::bind_rows(all_top_list) %>%
    dplyr::filter(collection %in% collections)

# Covert z score to p value and generate bonferonni adjusted alphas
all_top_df <- all_top_df %>%
    dplyr::mutate(p_val = 2 * pnorm(-abs_z_score)) %>%
    dplyr::group_by(z) %>%
    dplyr::mutate(bon_alpha = 0.05 / z) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(is_bon_filtered = bon_alpha < p_val)

# Process data for plotting
all_top_df$z <- factor(all_top_df$z,
                       levels = sort(as.numeric(paste(unique(all_top_df$z)))))

all_top_df <- all_top_df %>%
    dplyr::group_by(z, collection, dataset, algorithm) %>% 
    dplyr::mutate(num_instance = dplyr::n()) %>%
    dplyr::mutate(normalized_z = num_instance / as.numeric(paste(z))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(algorithm = toupper(algorithm))

all_top_df$algorithm <- factor(all_top_df$algorithm, levels = algorithms)
all_top_df$collection <- factor(all_top_df$collection, levels = collections)

print(dim(all_top_df))
head(all_top_df)

# How many genesets are filtered
filtered_genesets_df <- reshape2::melt(all_top_df,
                                       id.vars = c("algorithm",
                                                   "dataset",
                                                   "collection"),
                                       measure.vars = 'is_bon_filtered',
                                       value.name = "filtered") %>%
    dplyr::group_by(algorithm, dataset, collection) %>%
    dplyr::summarize(num_filtered = sum(filtered),
                     percent_filtered = sum(filtered) / n()) %>%
    dplyr::arrange(desc(num_filtered))

out_file <- file.path("results", "bonferonni_filtered_counts.tsv")
readr::write_tsv(filtered_genesets_df, out_file)

filtered_genesets_df

# Perform Bonferonni Correction
filtered_top_df <- all_top_df %>%
    dplyr::filter(!is_bon_filtered)

# Split TCGA, GTEx, and TARGET data
tcga_top_df <- filtered_top_df %>%
  dplyr::filter(dataset == "TCGA")

gtex_top_df <- filtered_top_df %>%
  dplyr::filter(dataset == "GTEX")

target_top_df <- filtered_top_df %>%
  dplyr::filter(dataset == "TARGET")

plot_top_feature_distribution <- function(top_df, density_how = "fill") {
    g <- ggplot(data = top_df,
                aes(x = z,
                    fill = collection,
                    group = collection)) +
    geom_density(alpha = 0.5,
                 position = density_how) +
    coord_flip() +
    scale_fill_manual(name = "",
                      values = c("GpC2CPG" = "#0F1FFF",
                                 "GpC2CPREACTOME" = "green",
                                 "GpC3TFT" = "#EDF67D",
                                 "GpH" = "#49DCB1",
                                 "GpC4CM" = "#7D82B8",
                                 "GpXCELL" = "#EF767A"),
                    labels = c("GpC2CPG" = "Perturbations (C2CGP)", 
                               "GpC2CPREACTOME" = "Pathways (REACTOME)",
                               "GpC3TFT" = "TF Targets (C3TFT)",
                               "GpH" = "Hallmark Gene Sets",
                               "GpXCELL" = "Cell Types (xCell)",
                               "GpC4CM" = "Cancer Modules (C4CM)")) +
    facet_wrap(dataset ~ algorithm, ncol = 5, scales = 'free_x') +
    xlab("Latent Dimensions") +
    ylab("Relative Density") +
    theme_bw() +
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"),
          strip.text.x = element_text(size = 6,
                                      margin = margin(t = 3,
                                                      b = 1.5,
                                                      l = 0,
                                                      r = 0)),
          axis.title = element_text(size = 9),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 5),
          legend.position = 'bottom',
          legend.text = element_text(size = 6),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0)) +
     guides(fill = guide_legend(reverse = T))
    
    return(g)
    
}

tcga_gg <- plot_top_feature_distribution(tcga_top_df)

tcga_gg

target_gg <- plot_top_feature_distribution(target_top_df, density_how = "fill")

target_gg

theme_gtex <- theme(strip.background = element_rect(colour = "black",
                                                    fill = "#fdfff4"),
                    strip.text.x = element_text(size = 6,
                                                margin = margin(t = 3,
                                                                b = 1.5,
                                                                l = 0,
                                                                r = 0)),
                    axis.title = element_text(size = 9),
                    axis.text.x = element_text(size = 8),
                    axis.text.y = element_text(size = 5),
                    legend.position = 'bottom',
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 6),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.box.margin = margin(-8, 0, 0, 0))

gtex_stack_gg <- ggplot(data = gtex_top_df,
                        aes(x = z, fill = algorithm, group = algorithm)) +
    geom_density(alpha = 0.5, position = "stack") +
    coord_flip() +
    facet_wrap(dataset ~ collection, ncol = 5) +
    xlab("Latent Dimensions") +
    ylab("Density (Stacked)") +
    scale_fill_manual(name = "Algorithm",
                     values = c("#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00")) +
    theme_bw() +
    theme_gtex +
    guides(fill = guide_legend(reverse = T))

gtex_stack_gg

gtex_fill_gg <- ggplot(data = gtex_top_df,
                       aes(x = z, fill = algorithm, group = algorithm)) +
    geom_density(alpha = 0.5, position = "fill") +
    coord_flip() +
    facet_wrap(dataset ~ collection, ncol = 5) +
    xlab("Latent Dimensions") +
    ylab("Relative Density") +
    scale_fill_manual(name = "Algorithm",
                     values = c("#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00")) +
    theme_bw() +
    theme_gtex +
    guides(fill = guide_legend(reverse = T))


gtex_fill_gg

algorithm_legend_gg <- cowplot::get_legend(gtex_fill_gg)

b_and_c_gg <- cowplot::plot_grid(
    gtex_stack_gg + theme(legend.position = "none"),
    gtex_fill_gg + theme(legend.position = "none"),
    ncol = 2,
    labels = c("c", "d")
)

b_and_c_gg <- cowplot::plot_grid(
    b_and_c_gg,
    algorithm_legend_gg,
    rel_heights = c(1, 0.15),
    nrow = 2
)

main_plot <- cowplot::plot_grid(
    tcga_gg,
    target_gg,
    b_and_c_gg,
    nrow = 3,
    rel_heights = c(1.05, 0.95, 1),
    labels = c("a", "b", "")
)

main_plot

for(extension in c('.png', '.pdf')) {
    gg_file <- paste0("top_feature_density", extension)
    gg_file <- file.path("figures", gg_file)
    cowplot::save_plot(filename = gg_file,
                       plot = main_plot,
                       base_height = 210,
                       base_width = 170,
                       units = "mm")
}
