
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("scripts", "viz_util.R"))

# Load all mutation classification results
mut_path <- file.path("results", "mutation")
full_mutation_df <- load_results(results_path = mut_path, process_output = TRUE)

# Load all cancertype classification results
cancertype_path <- file.path("results", "cancer-type")
full_cancertype_df <- load_results(results_path = cancertype_path, process_output = TRUE)

# Setup plotting logic for supplmental plots
genes <- unique(levels(full_mutation_df$gene_or_cancertype))
cancertypes <- as.character(unique(full_cancertype_df$gene_or_cancertype))

# Create mutation plots
gg_list <- list()
for (plot_idx in seq(1, 50, 10)) {
    end_idx <- plot_idx+9
    use_genes <- genes[plot_idx:end_idx]
    
    subset_df <- full_mutation_df %>%
        dplyr::filter(gene_or_cancertype %in% use_genes)
    
    gg_list[[use_genes[1]]] <- plot_mutation_figure(df = subset_df)
}

# Create cancertype plots
for (plot_idx in seq(1, 33, 11)) {
    end_idx <- plot_idx+10
    use_cancertypes <- cancertypes[plot_idx:end_idx]
    
    subset_df <- full_cancertype_df %>%
        dplyr::filter(gene_or_cancertype %in% use_cancertypes)
    subset_df$gene_or_cancertype <- as.character(subset_df$gene_or_cancertype)
    
    gg_list[[use_cancertypes[1]]] <- plot_mutation_figure(df = subset_df)
}

# Save a series of plots
for (g_idx in 1:length(names(gg_list))) {
    g_name <- names(gg_list)[g_idx]
    g <- gg_list[[g_name]] + theme(legend.position = "bottom")

    # Save Figure
    for (extension in c(".png", ".pdf")) {
      fig_file <- paste0("supplementary_figure_tcga_classify_auroc_plotindex_",
                         g_idx,
                         "_name_",
                         g_name,
                         extension)

      fig_file <- file.path("figures", fig_file)
      ggplot2::ggsave(filename = fig_file,
                      plot = g,
                      height = 150,
                      width = 170,
                      units = "mm")
    }
}

classifier_base_theme <-
    theme(strip.background = element_rect(colour = "black", fill = "#fdfff4"),
          strip.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          legend.position = "bottom",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3))

focus_cancertypes <- c("KIRP", "OV", "UCEC", "LUAD", "BRCA")
focus_cancertype_df <- full_cancertype_df %>%
  dplyr::filter(gene_or_cancertype %in% focus_cancertypes) %>%
  dplyr::mutate(gene_or_cancertype =
                  factor(gene_or_cancertype, levels = focus_cancertypes))

panel_a_gg <- plot_mutation_figure(df = focus_cancertype_df)
panel_a_gg <- panel_a_gg +
    classifier_base_theme +
    ylab("CV AUROC") +
    guides(color = guide_legend(keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))

panel_a_gg

focus_genes <- c("TP53", "PTEN", "PIK3CA", "KRAS", "TTN")
focus_mut_df <- full_mutation_df %>%
  dplyr::filter(gene_or_cancertype %in% focus_genes) %>%
  dplyr::mutate(gene_or_cancertype =
                  factor(gene_or_cancertype, levels = focus_genes))

panel_b_gg <- plot_mutation_figure(df = focus_mut_df)

panel_b_gg <- panel_b_gg +
    classifier_base_theme +
    ylab("CV AUROC") +
    guides(color = guide_legend(keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))

panel_b_gg

# Load Results
full_coef_results <- load_results(results_path = mut_path,
                                  file_string = "coefficients")

# Process Results data for plotting
coef_df <- full_coef_results[["metrics"]] %>%
  dplyr::filter(gene %in% focus_genes)

num_zero_df <- coef_df %>%
  dplyr::group_by(gene, signal, z_dim, seed, algorithm) %>%
  dplyr::summarize_at("weight", funs(sum(. == 0)))

denom_df <- coef_df %>%
  dplyr::group_by(gene, signal, z_dim, seed, algorithm) %>%
  dplyr::summarise(num_features = n())

num_zero_df <- num_zero_df %>%
  dplyr::full_join(denom_df, by = c('gene', 'signal', 'z_dim', 'seed',
                                    'algorithm')) %>%
  dplyr::mutate(percent_zero = weight / num_features)

cv_metrics_df <- full_mutation_df %>% dplyr::filter(data_type == "cv")

sparsity_metric_df <-
  dplyr::left_join(num_zero_df, cv_metrics_df,
                   by = c("gene" = "gene_or_cancertype", "signal", "seed",
                          "algorithm", "z_dim")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = factor(gene, levels = focus_genes)) %>%
  dplyr::mutate(signal = factor(signal, levels = c("signal", "shuffled")))

sparsity_metric_df <- sparsity_metric_df %>%
  dplyr::mutate(z_dim_shape = ifelse(as.numeric(paste(sparsity_metric_df$z_dim)) >= 20, 'z >= 20', 'z < 20'))

panel_c_gg <- ggplot(sparsity_metric_df,
                     aes(x = percent_zero,
                         y = auroc)) +
  geom_point(aes(color = algorithm,
                 shape = z_dim_shape),
             alpha = 0.4) +
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
  scale_shape_manual(values = c("+", "o")) +
  ylim(c(0.5, 1)) +
  ylab("CV AUROC") +
  xlab("Percent Zero Coefficients") + 
  facet_grid(signal ~ gene) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "#fdfff4"),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3)) +
  guides(color = guide_legend(order = 1,
                              keywidth = 0.1,
                              keyheight = 0.1,
                              default.unit = "inch",
                              override.aes = list(size = 1.4,
                                                  alpha = 1)),
         shape = guide_legend(order = 2,
                              keywidth = 0.1,
                              keyheight = 0.1,
                              default.unit = "inch",
                              title = "Dimension",
                              override.aes = list(size = 3,
                                                  alpha = 1)))

panel_c_gg

# Find Model with High Sparsity and High Performance
top_gene <- "TP53"
top_feature_search_df <- sparsity_metric_df %>%
  dplyr::filter(algorithm == 'DAE',
                gene == top_gene,
                auroc > 0.75,
                percent_zero > 0.8) %>% 
  dplyr::top_n(1, auroc)

# What feature is most explanatory in this model?
top_tp53_features <- coef_df %>%
    dplyr::filter(z_dim == top_feature_search_df$z_dim,
                  seed == top_feature_search_df$seed,
                  algorithm == 'DAE',
                  gene == "TP53",
                  signal == 'signal') %>%
    dplyr::arrange(weight)

top_tp53_features$ranked <- 1:nrow(top_tp53_features)

panel_d_gg <- ggplot(top_tp53_features,
                     aes(x = ranked, y = weight)) +
  geom_point(alpha = 0.3,
             size = 0.02) +
  xlab("Weight Rank") +
  ylab("Weight") +
  geom_text_repel(data = subset(top_tp53_features,
                                (weight > 0.15 | weight < -0.14)),
                  arrow = arrow(length = unit(0.02, 'npc')),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  box.padding = 0.68,
                  point.padding = 0.22,
                  size = 2.2,
                  fontface = 'italic',
                  aes(x = ranked, y = weight, label = feature)) +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7))

panel_d_gg

metric_col_type <- readr::cols(
    .default = readr::col_double(),
    predictor = readr::col_character(),
    signal = readr::col_character(),
    z_dim = readr::col_integer(),
    seed = readr::col_integer(),
    algorithm = readr::col_character(),
    data_type = readr::col_character()
)

# Find metrics for the specific model
top_model_path <- file.path("results", "mutation", top_gene)
auc_file <- file.path(top_model_path, paste0(top_gene, "_auc_threshold_metrics.tsv.gz"))

roc_df <- readr::read_tsv(auc_file,
                          col_types = metric_col_type) %>%
  dplyr::filter(z_dim == top_feature_search_df$z_dim,
                seed == top_feature_search_df$seed,
                algorithm == 'dae')

aupr_file <- file.path(top_model_path, paste0(top_gene, "_aupr_threshold_metrics.tsv.gz"))
pr_df <- readr::read_tsv(aupr_file,
                         col_types = metric_col_type) %>%
  dplyr::filter(z_dim == top_feature_search_df$z_dim,
                seed == top_feature_search_df$seed,
                algorithm == 'dae')

# Load Raw metrics
auc_raw_file <- file.path(top_model_path,
                          paste0(top_gene, "_raw_auc_threshold_metrics.tsv.gz"))
roc_raw_df <- readr::read_tsv(auc_raw_file, col_types = metric_col_type) 

aupr_raw_file <- file.path(top_model_path,
                           paste0(top_gene, "_raw_aupr_threshold_metrics.tsv.gz"))
pr_raw_df <- readr::read_tsv(aupr_raw_file, col_types = metric_col_type)

full_roc_df <- dplyr::bind_rows(roc_df, roc_raw_df)
full_roc_df$model_groups <- paste(full_roc_df$algorithm, full_roc_df$data_type)

full_pr_df <- dplyr::bind_rows(pr_df, pr_raw_df)
full_pr_df$model_groups <- paste(full_pr_df$algorithm, full_pr_df$data_type)

# Setup plotting variables
curve_colors <- c("#1b9e77",
                  "#d95f02",
                  "#7570b3",
                  "#737373",
                  "#bdbdbd",
                  "#d9d9d9")

curve_labels <- c("dae cv" = "DAE CV",
                  "dae test" = "DAE Test",
                  "dae train" = "DAE Train",
                  "raw cv" = "Raw CV",
                  "raw test" = "Raw Test",
                  "raw train" = "Raw Train")

curve_base_theme <-
    theme(axis.title = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3))

panel_e_gg <- ggplot(full_roc_df,
                     aes(x = fpr,
                         y = tpr,
                         color = model_groups)) +
  geom_step(aes(linetype = signal), size = 0.4) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = 1,
                   yend = 1),
               linetype = "dashed",
               color = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(name = "Models",
                     values = curve_colors,
                     labels = curve_labels) +
  scale_linetype_manual(name = "Data",
                        values = c("dashed",
                                   "solid"),
                     labels = c("signal" = "Signal",
                                "shuffled" = "Shuffled")) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_bw() +
  curve_base_theme +
  guides(color = guide_legend(order = 1,
                              keywidth = 0.1,
                              keyheight = 0.1,
                              default.unit = "inch",
                              override.aes = list(size = 0.8)),
         linetype = FALSE)


panel_e_gg

panel_f_gg <- ggplot(full_pr_df,
                     aes(x = recall,
                         y = precision,
                         color = model_groups)) +
  geom_step(aes(linetype = signal), size = 0.4) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(name = "Models",
                     values = curve_colors,
                     labels = curve_labels) +
  scale_linetype_manual(name = "Data",
                        values = c("dashed",
                                   "solid"),
                        labels = c("signal" = "Signal",
                                   "shuffled" = "Shuffled")) +
  xlab("Recall") +
  ylab("Precision") +
  theme_bw() +
  curve_base_theme +
  guides(color = guide_legend(order = 1,
                              keywidth = 0.1,
                              keyheight = 0.1,
                              default.unit = "inch",
                              override.aes = list(size = 0.8)),
         linetype = FALSE)

panel_f_gg

a_and_b_legend_gg <- cowplot::get_legend(panel_a_gg)

a_and_b_gg <- cowplot::plot_grid(
    panel_a_gg + theme(legend.position = "none"),
    panel_b_gg + theme(legend.position = "none"),
    rel_widths = c(1, 1),
    labels = c("A", "B"),
    ncol = 2
)

a_and_b_gg <- cowplot::plot_grid(
    a_and_b_gg,
    a_and_b_legend_gg,
    rel_heights = c(1, 0.04),
    nrow = 2
)

e_and_f_legend_gg <- cowplot::get_legend(panel_e_gg)

e_and_f_gg <- cowplot::plot_grid(
    panel_e_gg + theme(legend.position = "none"),
    panel_f_gg + theme(legend.position = "none"),
    rel_widths = c(1, 1),
    labels = c("E", "F"),
    ncol = 2,
    nrow = 1
)

e_and_f_gg <- cowplot::plot_grid(
    e_and_f_gg,
    e_and_f_legend_gg,
    rel_widths = c(1, 0.15),
    ncol = 2
)

d_e_and_f_gg <- cowplot::plot_grid(
    panel_d_gg,
    e_and_f_gg,
    rel_widths = c(0.4, 1),
    labels = c("D", ""),
    ncol = 2
)

full_gg <- cowplot::plot_grid(
    a_and_b_gg,
    panel_c_gg,
    d_e_and_f_gg,
    nrow = 3,
    labels = c("", "C", ""),
    rel_heights = c(1.2, 1, 0.7)
)

full_gg

for(extension in c('.png', '.pdf')) {
    gg_file <- paste0("tcga_biobombe_main_figure", extension)
    gg_file <- file.path("figures", gg_file)
    cowplot::save_plot(filename = gg_file,
                       plot = full_gg,
                       dpi = 300,
                       base_height = 190,
                       base_width = 170,
                       units = "mm")
}
