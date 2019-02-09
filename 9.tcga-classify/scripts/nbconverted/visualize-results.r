
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("scripts", "viz_util.R"))

# Load individual models
cancertype_path <- file.path("results", "cancer-type")

full_cancertype_info <- load_results(results_path = cancertype_path,
                                     file_string = "classify_metrics",
                                     process_output = FALSE)

full_cancertype_df <- process_results(df = full_cancertype_info$metrics,
                                      raw_df = full_cancertype_info$raw_metrics)

# Load ensemble models
cancertype_ensemble_path <- file.path("results", "cancer-type_ensemble")

full_cancertype_ensemble_df <- load_results(results_path = cancertype_ensemble_path,
                                            file_string = "classify_metrics",
                                            process_output = FALSE,
                                            process_ensemble = TRUE)

full_cancertype_ensemble_df <- dplyr::bind_rows(full_cancertype_ensemble_df)
full_cancertype_ensemble_df <- process_results(df = full_cancertype_ensemble_df,
                                               raw_df = full_cancertype_info$raw_metrics)

full_cancertype_ensemble_df$algorithm <- full_cancertype_ensemble_df$algorithm %>%
      dplyr::recode_factor("all_ensemble" = "Model Ensemble",
                           "vae_ensemble" = "VAE Ensemble")

# Load classifier with all features
cancertype_all_ensemble_path <- file.path("results", "cancer-type_ensemble_all")

full_cancertype_all_ensemble_df <-
    load_results(results_path = cancertype_all_ensemble_path,
                 file_string = "classify_metrics",
                 process_output = FALSE,
                 process_all_features = TRUE)

# Make sure to filter only cross validation data (other data is prefiltered)
full_cancertype_all_ensemble_df <- full_cancertype_all_ensemble_df %>%
    dplyr::filter(data_type == "cv")

# Load individual models
mut_path <- file.path("results", "mutation")

full_mutation_info <- load_results(results_path = mut_path,
                                   file_string = "classify_metrics",
                                   process_output = FALSE)

full_mutation_df <- process_results(df = full_mutation_info$metrics,
                                    raw_df = full_mutation_info$raw_metrics)

# Load ensemble models
mut_ensemble_path <- file.path("results", "mutation_ensemble")

full_mutation_ensemble_df <- load_results(results_path = mut_ensemble_path,
                                          file_string = "classify_metrics",
                                          process_output = FALSE,
                                          process_ensemble = TRUE)

full_mutation_ensemble_df <- dplyr::bind_rows(full_mutation_ensemble_df)
full_mutation_ensemble_df <- process_results(df = full_mutation_ensemble_df,
                                             raw_df = full_mutation_info$raw_metrics)

full_mutation_ensemble_df$algorithm <- full_mutation_ensemble_df$algorithm %>%
      dplyr::recode_factor("all_ensemble" = "Model Ensemble",
                           "vae_ensemble" = "VAE Ensemble")

# Load classifier with all features
mutation_all_ensemble_path <- file.path("results", "mutation_ensemble_all")

full_mutation_all_ensemble_df <-
    load_results(results_path = mutation_all_ensemble_path,
                 file_string = "classify_metrics",
                 process_output = FALSE,
                 process_all_features = TRUE)

# Make sure to filter only cross validation data (other data is prefiltered)
full_mutation_all_ensemble_df <- full_mutation_all_ensemble_df %>%
    dplyr::filter(data_type == "cv")

# Setup plotting logic for supplmental plots
cancertypes <- as.character(unique(full_cancertype_df$gene_or_cancertype))
genes <- unique(levels(full_mutation_df$gene_or_cancertype))

# Create cancertype plots
gg_list <- list()

n_cancertype_plots <- 33
n_cancertype_plots_per_page <- 11

for (plot_idx in seq(1, n_cancertype_plots, n_cancertype_plots_per_page)) {
    end_idx <- plot_idx + n_cancertype_plots_per_page - 1
    use_cancertypes <- cancertypes[plot_idx:end_idx]
    
    subset_df <- full_cancertype_df %>%
        dplyr::filter(gene_or_cancertype %in% use_cancertypes)
    subset_df$gene_or_cancertype <- as.character(subset_df$gene_or_cancertype)
    
    gg_list[[use_cancertypes[1]]] <- plot_mutation_figure(df = subset_df,
                                                          auroc_or_aupr = "aupr")
}

# Create mutation plots (append to existing gglist)
n_top_mutations <- 50
n_plots_per_page <- 10

for (plot_idx in seq(1, n_top_mutations, n_plots_per_page)) {
    end_idx <- plot_idx + n_plots_per_page - 1
    use_genes <- genes[plot_idx:end_idx]
    
    subset_df <- full_mutation_df %>%
        dplyr::filter(gene_or_cancertype %in% use_genes)
    
    gg_list[[use_genes[1]]] <- plot_mutation_figure(df = subset_df,
                                                    auroc_or_aupr = "aupr")
}

# Save a series of plots
sup_fig_height = 150
sup_fig_width = 170 

for (plot_idx in 1:length(gg_list)) {
    
    # Extract index name of list
    plot_name <- names(gg_list)[plot_idx]
    
    # Put legend on the bottom
    g <- gg_list[[plot_name]] +
        theme(legend.position = "bottom")

    # Save Figure
    for (extension in c(".png", ".pdf")) {

        # Save figure with the index and name of the plot - will present in this
        # order in the supplementary figure
        fig_file <- paste0("supplementary_figure_tcga_classify_auc_plotindex_",
                           plot_idx,
                           "_name_",
                           plot_name,
                           extension)
        
        fig_file <- file.path("figures", fig_file)
        
        ggplot2::ggsave(filename = fig_file,
                        plot = g,
                        height = sup_fig_height,
                        width = sup_fig_width,
                        units = "mm")
    }
}

# Generate base theme
classifier_base_theme <-
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"),
          strip.text.x = element_text(size = 6,
                                      margin = margin(
                                          t = 4,
                                          b = 3,
                                          l = 0,
                                          r = 0)
                                     ),
          strip.text.y = element_text(size = 6,
                                      margin = margin(
                                          t = 0,
                                          b = 0,
                                          l = 3,
                                          r = 4)
                                     ),
          axis.title = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          legend.position = "bottom",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3))

focus_cancertypes <- c("KIRP", "OV", "UCEC", "LUAD", "BRCA")
focus_genes <- c("TP53", "PTEN", "PIK3CA", "KRAS", "TTN")
algorithms <- c("pca", "ica", "nmf", "dae", "vae")
algorithm_refactor <- c("PCA", "ICA", "NMF", "DAE", "VAE", "VAE Ensemble", "Model Ensemble")
plot_algorithm_refactor <- c("PCA", "ICA", "NMF", "DAE", "VAE", "VAE Ens.", "Ensemble")

focus_cancertype_df <- full_cancertype_df %>%
    dplyr::filter(gene_or_cancertype %in% focus_cancertypes) %>%
    dplyr::mutate(gene_or_cancertype =
                  factor(gene_or_cancertype,
                         levels = focus_cancertypes))

focus_cancertype_ensemble_df <- full_cancertype_ensemble_df %>%
    dplyr::filter(gene_or_cancertype %in% focus_cancertypes) %>%
    dplyr::mutate(gene_or_cancertype =
                  factor(gene_or_cancertype,
                         levels = focus_cancertypes))

plot_ready_cancertype_df <- dplyr::bind_rows(
    focus_cancertype_df,
    focus_cancertype_ensemble_df
)

# Abbreviate and reorder factors for plotting
plot_ready_cancertype_df$algorithm <- plot_ready_cancertype_df$algorithm %>%
      dplyr::recode_factor("Model Ensemble" = "Ensemble",
                           "VAE Ensemble" = "VAE Ens.")

plot_ready_cancertype_df$algorithm <- factor(plot_ready_cancertype_df$algorithm,
                                             levels = plot_algorithm_refactor)


panel_a_gg <- plot_mutation_figure(df = plot_ready_cancertype_df,
                                   auroc_or_aupr = "aupr")
panel_a_gg <- panel_a_gg +
    classifier_base_theme +
    ylab("CV AUPR") +
    guides(color = guide_legend(keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))

panel_a_gg

focus_mut_df <- full_mutation_df %>%
  dplyr::filter(gene_or_cancertype %in% focus_genes) %>%
  dplyr::mutate(gene_or_cancertype =
                  factor(gene_or_cancertype,
                         levels = focus_genes))

plot_ready_mut_df <- dplyr::bind_rows(focus_mut_df,
                                      full_mutation_ensemble_df)

# Abbreviate and reorder factors for plotting
plot_ready_mut_df$algorithm <- plot_ready_mut_df$algorithm %>%
      dplyr::recode_factor("Model Ensemble" = "Ensemble",
                           "VAE Ensemble" = "VAE Ens.")

plot_ready_mut_df$algorithm <- factor(plot_ready_mut_df$algorithm,
                                      levels = plot_algorithm_refactor)

panel_b_gg <- plot_mutation_figure(df = plot_ready_mut_df,
                                   auroc_or_aupr = "aupr")

panel_b_gg <- panel_b_gg +
    classifier_base_theme +
    ylab("CV AUPR") +
    guides(color = guide_legend(keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))

panel_b_gg

cancertype_delta_auc_df <- process_delta_auc(focus_cancertype_df,
                                             auroc_or_aupr = "aupr",
                                             seed = "165158")

cancertype_delta_auc_ensemble_df <- process_delta_auc(focus_cancertype_ensemble_df,
                                                      auroc_or_aupr = "aupr",
                                                      seed = "ensemble")

cancertype_delta_auc_ensemble_all_features_df <-
    process_delta_auc(full_cancertype_all_ensemble_df,
                      auroc_or_aupr = "aupr",
                      seed = "ensemble_all_features")

line_plot_ready_cancertype_ensemble_df <- dplyr::bind_rows(
    cancertype_delta_auc_df,
    cancertype_delta_auc_ensemble_df
)

# Correct factor coercion error
line_plot_ready_cancertype_ensemble_df$algorithm <-
    factor(line_plot_ready_cancertype_ensemble_df$algorithm,
           levels = algorithm_refactor)

head(line_plot_ready_cancertype_ensemble_df)

panel_c_gg <- plot_delta_auc_simple(plot_df = line_plot_ready_cancertype_ensemble_df,
                                    auroc_or_aupr = "aupr",
                                    plot_ensemble = TRUE,
                                    plot_all_features = TRUE,
                                    all_feature_df = cancertype_delta_auc_ensemble_all_features_df,
                                    plot_title = "Cancer Type")

panel_c_gg

mutation_delta_auroc_df <- process_delta_auc(focus_mut_df,
                                             auroc_or_aupr = "aupr",
                                             seed = "165158")

mutation_delta_auroc_ensemble_df <- process_delta_auc(full_mutation_ensemble_df,
                                                      auroc_or_aupr = "aupr",
                                                      seed = "ensemble")

mutation_delta_auc_ensemble_all_features_df <-
    process_delta_auc(full_mutation_all_ensemble_df,
                      auroc_or_aupr = "aupr",
                      seed = "ensemble_all_features")

line_plot_ready_mutation_ensemble_df <- dplyr::bind_rows(
    mutation_delta_auroc_df,
    mutation_delta_auroc_ensemble_df
)

# Correct factor coercion error
line_plot_ready_mutation_ensemble_df$algorithm <-
    factor(line_plot_ready_mutation_ensemble_df$algorithm,
           levels = algorithm_refactor)

head(line_plot_ready_mutation_ensemble_df)

panel_d_gg <- plot_delta_auc_simple(plot_df = line_plot_ready_mutation_ensemble_df,
                                    auroc_or_aupr = "aupr",
                                    plot_ensemble = TRUE,
                                    plot_all_features = TRUE,
                                    all_feature_df = mutation_delta_auc_ensemble_all_features_df,
                                    plot_title = "Mutations")

panel_d_gg <- panel_d_gg + ylim(c(-0.04, 0.3))
panel_d_gg

# Load Results
full_coef_results <- load_results(results_path = mut_path,
                                  file_string = "coefficients")

full_coef_ensemble_results_df <- load_results(results_path = mut_ensemble_path,
                                              file_string = "coefficients",
                                              process_ensemble = TRUE,
                                              process_output = FALSE)

full_coef_ensemble_all_features_results_df <-
    load_results(results_path = mutation_all_ensemble_path,
                 file_string = "coefficients",
                 process_all_features = TRUE,
                 process_output = FALSE)

# Process full coefficient results
coef_df <- full_coef_results[["metrics"]]
raw_coef_df <- full_coef_results[["raw_metrics"]]

# Process ensemble results
full_coef_ensemble_results_df <- dplyr::bind_rows(full_coef_ensemble_results_df)

# Adjust factor
full_coef_ensemble_results_df$algorithm <- full_coef_ensemble_results_df$algorithm %>%
      dplyr::recode_factor("all_ensemble" = "Model Ensemble",
                           "vae_ensemble" = "VAE Ensemble")

# Process sparsity
sparsity_metric_df <- process_sparsity(coef_df = coef_df,
                                       mut_df = full_mutation_df,
                                       focus_genes = focus_genes)

raw_sparsity_metric_df <- process_sparsity(coef_df = raw_coef_df,
                                           mut_df = full_mutation_info$raw_metrics,
                                           focus_genes = focus_genes)

ensemble_sparsity_metric_df <- process_sparsity(coef_df = full_coef_ensemble_results_df,
                                                mut_df = full_mutation_ensemble_df,
                                                focus_genes = focus_genes,
                                                process_ensemble = TRUE)

ensemble_all_feature_sparsity_metric_df <-
    process_sparsity(coef_df = full_coef_ensemble_all_features_results_df,
                     mut_df = full_mutation_all_ensemble_df,
                     focus_genes = focus_genes,
                     process_all_features = TRUE)

# Setup Plotting logic
scale_colors <- c(
    "PCA" = "#e41a1c",
    "ICA" = "#377eb8",
    "NMF" = "#4daf4a",
    "DAE" = "#984ea3",
    "VAE" = "#ff7f00",
    "Model Ensemble" = "brown",
    "VAE Ensemble" = "grey50"
)

scale_labels <- c(
    "PCA" = "PCA",
    "ICA" = "ICA",
    "NMF" = "NMF",
    "DAE" = "DAE",
    "VAE" = "VAE",
    "Model Ensemble" = "Model Ensemble",
    "VAE Ensemble" = "VAE Ensemble"
)

plot_ready_sparsity_df <- dplyr::bind_rows(
    sparsity_metric_df,
    ensemble_sparsity_metric_df
)

# Correct factor coercion error
plot_ready_sparsity_df$algorithm <-
    factor(plot_ready_sparsity_df$algorithm,
           levels = algorithm_refactor)

head(plot_ready_sparsity_df)

panel_e_gg <- ggplot(plot_ready_sparsity_df,
                     aes(x = percent_zero,
                         y = aupr)) +
  geom_point(aes(color = algorithm,
                 shape = z_dim_shape),
             size = 0.9,
             alpha = 0.5) +
  geom_point(data = raw_sparsity_metric_df,
             aes(x = percent_zero,
                 y = aupr),
             alpha = 0.8,
             size = 1,
             color = "black",
             fill = "grey",
             pch = 22) +
  geom_point(data = ensemble_all_feature_sparsity_metric_df,
             aes(x = percent_zero,
                 y = aupr),
             alpha = 0.8,
             size = 1,
             color = "#176620",
             fill = "#05a818",
             pch = 23) +
  scale_color_manual(name = "Algorithm",
                     values = scale_colors,
                     labels = scale_labels) +
  scale_shape_manual(values = c("+", "o")) +
  ylim(c(0.2, 1)) +
  ylab("CV AUPR") +
  xlab("Percent Zero Coefficients") + 
  facet_grid(signal ~ gene) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "#fdfff4"),
        strip.text = element_text(size = 7),
        strip.text.x = element_text(size = 6,
                                      margin = margin(
                                          t = 4,
                                          b = 3,
                                          l = 0,
                                          r = 0)
                                     ),
          strip.text.y = element_text(size = 6,
                                      margin = margin(
                                          t = 0,
                                          b = 0,
                                          l = 3,
                                          r = 4)
                                     ),
        axis.title = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
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

panel_e_gg

# Create a dummy plot to use in cowplot legend for sparsity figure
custom_data <- as.data.frame(cbind(c(1, 2), c(3, 4), c("A", "B")))
colnames(custom_data) <- c('one', 'two', 'three')
custom_data

custom_gg <- ggplot(custom_data,
                    aes(x = one, y = two)) +
    geom_point(aes(shape = three)) +
    scale_shape_manual(name = "Other",
                       values = c(22, 23),
                       labels = c("Raw",
                                  "All Ensemble")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3)) +
  guides(shape = guide_legend(order = 2,
                              keywidth = 0.1,
                              keyheight = 0.1,
                              default.unit = "inch",
                              title = "Other",
                              override.aes = list(fill = c("grey",
                                                             "#05a818"),
                                                  color = c("black",
                                                            "#176620"),
                                                  size = 1.7,
                                                  alpha = 1)))
custom_gg

# Find Model with High Sparsity and High Performance
top_gene <- "TP53"

top_feature_search_df <- ensemble_all_feature_sparsity_metric_df %>%
  dplyr::filter(algorithm == 'all_feature_ensemble',
                gene == top_gene) %>% 
  dplyr::top_n(1, aupr)

outfile <- file.path("results", "top_model_ensemble_all_features_tp53_feature_for_followup.tsv")
readr::write_tsv(top_feature_search_df, outfile)

top_feature_search_df

# This is an ensemble model, what feature is most explanatory in this model?
full_coef_ensemble_all_features_results_df$seed <- "ensemble_all_features"

top_tp53_features <- full_coef_ensemble_all_features_results_df %>%
    dplyr::filter(z_dim == top_feature_search_df$z_dim,
                  seed == top_feature_search_df$seed,
                  algorithm == 'all_feature_ensemble',
                  gene == top_gene,
                  signal == 'signal') %>%
    dplyr::arrange(weight)

top_tp53_features$ranked <- 1:nrow(top_tp53_features)

# Process plotting labels
algorithm_assign <- c()
algorithm_label <- c()
all_ks <- c()
for (feature in strsplit(top_tp53_features$feature, "_")) {
    feature_id_unlist <- unlist(feature)
    algorithm <- feature_id_unlist[1]
    feature <- feature_id_unlist[2]
    k <- feature_id_unlist[4]
    algorithm_label <- c(algorithm_label, algorithm)
    algorithm_assign <- c(algorithm_assign,
                        paste0(toupper(algorithm),
                               "_",
                               k))
    all_ks <- c(all_ks, k)
}

top_tp53_features$algorithm_label <- algorithm_label
top_tp53_features$algorithm_assign <- algorithm_assign
top_tp53_features$k <- all_ks
top_tp53_features$pos_label <- "Positive"
top_tp53_features[top_tp53_features$weight < 0, "pos_label"] <- 'Negative'

weight_rank_gg <- ggplot(top_tp53_features,
                     aes(x = ranked, y = weight)) +
  geom_point(alpha = 0.3,
             size = 0.02) +
  xlab("Weight Rank") +
  ylab("Weight") +
  geom_text_repel(data = subset(top_tp53_features,
                                (weight > 0.05 | weight < -0.045)),
                  arrow = arrow(length = unit(0.02, 'npc')),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  box.padding = 0.4,
                  point.padding = 0.22,
                  size = 1.4,
                  fontface = 'italic',
                  aes(x = ranked,
                      y = weight,
                      label = algorithm_assign)) +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

weight_rank_gg

metric_col_type <- readr::cols(
    .default = readr::col_double(),
    predictor = readr::col_character(),
    signal = readr::col_character(),
    z_dim = readr::col_character(),
    seed = readr::col_character(),
    algorithm = readr::col_character(),
    data_type = readr::col_character()
)

# Find metrics for the specific model
top_model_path <- file.path("results",
                            "mutation_ensemble_all",
                            top_gene)

auc_file <- file.path(top_model_path,
                      paste0(top_gene,
                             "_ensemble_all_features_auc_threshold_metrics.tsv.gz"))

roc_df <- readr::read_tsv(auc_file,
                          col_types = metric_col_type) %>%
  dplyr::filter(z_dim == top_feature_search_df$z_dim,
                seed == top_feature_search_df$seed,
                algorithm == 'all_feature_ensemble')

aupr_file <- file.path(top_model_path,
                       paste0(top_gene,
                              "_ensemble_all_features_aupr_threshold_metrics.tsv.gz"))

pr_df <- readr::read_tsv(aupr_file,
                         col_types = metric_col_type) %>%
  dplyr::filter(z_dim == top_feature_search_df$z_dim,
                seed == top_feature_search_df$seed,
                algorithm == 'all_feature_ensemble')

# Load Raw metrics
raw_model_path <- file.path("results", "mutation", top_gene)

auc_raw_file <- file.path(raw_model_path,
                          paste0(top_gene, "_raw_auc_threshold_metrics.tsv.gz"))
roc_raw_df <- readr::read_tsv(auc_raw_file, col_types = metric_col_type) 

aupr_raw_file <- file.path(raw_model_path,
                           paste0(top_gene, "_raw_aupr_threshold_metrics.tsv.gz"))
pr_raw_df <- readr::read_tsv(aupr_raw_file, col_types = metric_col_type)

# Process auc data
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

curve_labels <- c("all_feature_ensemble cv" = "Ensemble CV",
                  "all_feature_ensemble test" = "Ensemble Test",
                  "all_feature_ensemble train" = "Ensemble Train",
                  "raw cv" = "Raw CV",
                  "raw test" = "Raw Test",
                  "raw train" = "Raw Train")

curve_base_theme <-
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3))

panel_f_gg <- ggplot(full_roc_df,
                     aes(x = fpr,
                         y = tpr,
                         color = model_groups)) +
    geom_step(aes(linetype = signal),
              alpha = 0.7,
              size = 0.25) +
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

panel_f_gg

panel_g_gg <- ggplot(full_pr_df,
                     aes(x = recall,
                         y = precision,
                         color = model_groups)) +
  geom_step(aes(linetype = signal),
            alpha = 0.7,
            size = 0.25) +
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

panel_g_gg

algorithm_contribution_df <- top_tp53_features %>%
    dplyr::group_by(algorithm_label, pos_label, k) %>%
    dplyr::summarize(algorithm_contribution = sum(abs)) %>%
    dplyr::arrange(desc(algorithm_contribution)) %>%
    dplyr::filter(algorithm_contribution != 0)

algorithm_contribution_df$algorithm_label <-
    factor(algorithm_contribution_df$algorithm_label,
           levels = c("pca", "ica", "nmf", "dae", "vae", "log10"))

algorithm_contribution_df$pos_label <-
    factor(algorithm_contribution_df$pos_label,
           levels = c("Negative", "Positive"))

algorithm_contribution_df[algorithm_contribution_df$pos_label != "Positive", "algorithm_contribution"] <-
    algorithm_contribution_df[algorithm_contribution_df$pos_label != "Positive", "algorithm_contribution"] * -1

algorithm_contribution_df$k <-
    factor(algorithm_contribution_df$k,
           levels = sort(as.numeric(paste(unique(algorithm_contribution_df$k)))))

head(algorithm_contribution_df)

# Setup Plotting logic
bar_colors <- c(
    "pca" = "#e41a1c",
    "ica" = "#377eb8",
    "nmf" = "#4daf4a",
    "dae" = "#984ea3",
    "vae" = "#ff7f00",
    "log10" = "grey75"
)

bar_labels <- c(
    "pca" = "PCA",
    "ica" = "ICA",
    "nmf" = "NMF",
    "dae" = "DAE",
    "vae" = "VAE",
    "log10" = "log10"
)

panel_h_gg <- ggplot(algorithm_contribution_df,
                     aes(x = k,
                         y = algorithm_contribution,
                         fill = algorithm_label)) +
    geom_bar(stat = "identity",
             position = "stack") +
    geom_hline(yintercept = 0,
               lwd = 0.2,
               color = "black", 
               linetype = "dashed") +
    scale_fill_manual(name = "",
                      labels = bar_labels,
                      values = bar_colors) +
    ylab("Weight Sum") +
    xlab("k Dimension") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6, angle = 90),
          axis.text.y = element_text(size = 6),
          legend.position = "right",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3)) +
    guides(fill = guide_legend(order = 1,
                               keywidth = 0.1,
                               keyheight = 0.1,
                               default.unit = "inch",
                               override.aes = list(size = 0.8)))

panel_h_gg

a_and_b_legend_gg <- cowplot::get_legend(panel_a_gg)

a_and_b_gg <- cowplot::plot_grid(
    panel_a_gg + theme(legend.position = "none"),
    panel_b_gg + theme(legend.position = "none"),
    rel_widths = c(1, 1),
    labels = c("a", "b"),
    ncol = 2
)

a_and_b_gg <- cowplot::plot_grid(
    a_and_b_gg,
    a_and_b_legend_gg,
    rel_heights = c(1, 0.04),
    nrow = 2
)

c_and_d_gg <- cowplot::plot_grid(
    panel_c_gg +
        theme(legend.position = "none",
              plot.title = element_text(margin = margin(b = 1,
                                                        unit = "pt")),
              plot.margin = unit(c(5.5, 5.5, 0, 5.5),
                                 "pt")) +
        xlab(""),
    panel_d_gg + theme(legend.position = "none",
                       plot.title = element_text(margin = margin(b = 1,
                                                                 unit = "pt")),
                       plot.margin = unit(c(3, 5.5, 5.5, 5.5),
                                          "pt")),
    nrow = 2,
    labels = c("c", "d")
)

panel_e_legend <- cowplot::get_legend(panel_e_gg +
                                     theme(legend.justification = c(0.54, 0.6)))
custom_legend <- cowplot::get_legend(custom_gg +
                                     theme(legend.justification = c(0.30, 0.55)))

panel_e_full_legend <- (
    cowplot::plot_grid(
        cowplot::ggdraw(),
        panel_e_legend,
        custom_legend,
        ncol = 1,
        nrow = 3,
        rel_heights = c(0.5, 1, 1)
    )
)

panel_e_with_legend_gg = (
    cowplot::plot_grid(
        panel_e_gg + theme(legend.position = "none"),
        panel_e_full_legend,
        rel_widths = c(1, 0.2),
        ncol = 2
    )
)

c_d_and_e_gg <- cowplot::plot_grid(
    c_and_d_gg,
    panel_e_with_legend_gg,
    rel_widths = c(0.4, 1),
    labels = c("", "e"),
    ncol = 2
)

f_and_g_legend_gg <- cowplot::get_legend(panel_f_gg +
                                         theme(legend.justification = c(0.5, 0.7)))

f_and_g_gg <- cowplot::plot_grid(
    panel_f_gg + theme(legend.position = "none"),
    panel_g_gg + theme(legend.position = "none"),
    rel_widths = c(1, 1),
    labels = c("f", "g"),
    ncol = 2,
    nrow = 1
)

f_and_g_gg <- cowplot::plot_grid(
    f_and_g_gg,
    f_and_g_legend_gg,
    rel_widths = c(1, 0.2),
    ncol = 2
)

f_g_and_h_gg <- cowplot::plot_grid(
    f_and_g_gg,
    panel_h_gg + theme(legend.position = "right"),
    rel_widths = c(1, 0.6),
    labels = c("", "h"),
    ncol = 2
)

full_gg <- cowplot::plot_grid(
    a_and_b_gg,
    c_d_and_e_gg,
    f_g_and_h_gg,
    nrow = 3,
    labels = c("", "", ""),
    rel_heights = c(1.5, 0.9, 0.6)
)

full_gg

for(extension in c('.png', '.pdf')) {
    gg_file <- paste0("tcga_biobombe_main_figure", extension)
    gg_file <- file.path("figures", gg_file)
    cowplot::save_plot(filename = gg_file,
                       plot = full_gg,
                       dpi = 300,
                       base_height = 200,
                       base_width = 170,
                       units = "mm")
}
