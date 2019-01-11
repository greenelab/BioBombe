
suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpmisc))

# Load custom plotting functions
util_file = file.path("..", "6.analyze-weights", "scripts", "utils.R")
source(util_file)

set.seed(123)

# Load and process data
results_file <- file.path('results', 'gtex_vae_example_interpret_compression.tsv')
interpret_data_df <- (
    readr::read_tsv(results_file,
                    col_types = readr::cols(.default = readr::col_character(),
                                            raw_score = readr::col_double(),
                                            z_score = readr::col_double()))
    )

interpret_data_df$full_feature <- factor(interpret_data_df$full_feature,
                                  levels = c("vae_0_two", "vae_1_two", "vae_0_three",
                                             "vae_1_three", "vae_2_three"))


# Setup plotting logic
vae_labels <- c("vae_0_two" = "z = 2 (0)",
                "vae_1_two" = "z = 2 (1)",
                "vae_0_three" = "z = 3 (0)",
                "vae_1_three" = "z = 3 (1)",
                "vae_2_three" = "z = 3 (2)")

vae_colors <- c("#c994c7", "#dd1c77", "#78c679", "#31a354", "#006837")

color_logic <- (interpret_data_df$z_score > 13 | interpret_data_df$z_score < -18)

# Plot
panel_a_gg <- ggplot(interpret_data_df,
                     aes(x = raw_score,
                         y = z_score)) +
    geom_point(aes(color = full_feature),
               size = 0.14,
               alpha = 0.65) +
    scale_color_manual(name = "",
                       values = vae_colors,
                       labels =  vae_labels) +
    geom_text_repel(data = subset(interpret_data_df, color_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    size = 1.5,
                    fontface = "italic",
                    box.padding = 0.25,
                    point.padding = 0.15,
                    aes(x = raw_score,
                        y = z_score,
                        label = variable)) +
    xlab("BioBombe Raw Score") +
    ylab("BioBombe Z Score") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = 'bottom',
          legend.text = element_text(size = 4.7),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0)) +
    guides(color = guide_legend(nrow = 1,
                                ncol = 5,
                                byrow = FALSE,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch",
                                override.aes = list(size = 1.4,
                                                    alpha = 1)))

panel_a_gg

# Load and process data
file <- file.path('results', 'gtex_vae_example_differentiating_features.tsv')
feature_info_df <- readr::read_tsv(file,
                                   col_types = readr::cols(.default = readr::col_double(),
                                                           variable = readr::col_character()))

# Setup plotting logic
color_logic = feature_info_df$abs_diff > 3.5 | feature_info_df$two > 13  | feature_info_df$three > 13

panel_b_gg <- ggplot(feature_info_df, aes(x = three, y = two)) +
    geom_point(alpha = 0.5,
               size = 0.6,
               shape = 16,
               color = ifelse(color_logic, "red", "grey50")) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    geom_text_repel(data = subset(feature_info_df, color_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    size = 1.5,
                    fontface = "italic",
                    box.padding = 0.25,
                    point.padding = 0.15,
                    aes(x = three,
                        y = two,
                        label = variable)) +
    xlab("Mean Z Score for VAE (z = 3)") +
    ylab("Mean Z Score for VAE (z = 2)") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
    
panel_b_gg

gene_set_dir <- file.path("..", "6.analyze-weights", "results", "gtex", "gpxcell", "signal")
metaedge <- "GpXCELL"
dataset <- "GTEX"

line_plot_theme <-
    theme(strip.background = element_rect(colour = "black", fill = "#fdfff4"),
          strip.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          axis.title.x = element_text(margin = margin(t = 0.1, r = 0, b = 0, l = 0, unit = 'cm')),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 4.7),
          legend.title = element_text(size = 6),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3))

panel_c_gg <- plot_gene_set(gene_set = "Neutrophils_HPCA_2",
                            gene_set_dir = gene_set_dir,
                            metaedge = metaedge,
                            dataset = dataset,
                            show_plot = FALSE,
                            shuffled = FALSE,
                            return_top = FALSE,
                            return_plot = TRUE)

panel_c_gg <- panel_c_gg + line_plot_theme

panel_c_gg

panel_d_gg <- plot_gene_set(gene_set = "Monocytes_FANTOM_2",
                            gene_set_dir = gene_set_dir,
                            metaedge = metaedge,
                            dataset = dataset,
                            show_plot = FALSE,
                            shuffled = FALSE,
                            return_top = FALSE,
                            return_plot = TRUE)

panel_d_gg <- panel_d_gg + line_plot_theme

panel_d_gg

# Load and process data
file <- file.path('results', 'latent_feature_enrichment_comparison_neutrophil_genesets.tsv')
geneset_weights_df <- (
    readr::read_tsv(file,
                    col_types = readr::cols(.default = readr::col_double(),
                                            model_type_z3 = readr::col_character(),
                                            variable = readr::col_character(),
                                            algorithm_z3 = readr::col_character(),
                                            model_type_z14 = readr::col_character(),
                                            algorithm_z14 = readr::col_character()))
    ) %>%
    dplyr::filter(feature_z3 == 0,
                  feature_z14 == 10)

head(geneset_weights_df, 3)

# Process plotting logic
color_logic <- ((geneset_weights_df$z_score_z3 < -9 | geneset_weights_df$z_score_z3 > 10) |
                (geneset_weights_df$z_score_z14 < -10 | geneset_weights_df$z_score_z14 > 10))

# Formula to plot linear model
formula <- x ~ y

# Generate plot
panel_e_gg <- ggplot(geneset_weights_df,
                     aes(x = z_score_z3, y = z_score_z14)) +
    geom_point(alpha = 0.5,
               size = 0.6,
               shape = 16,
               color = ifelse(color_logic, "red", "grey50")) +
    geom_smooth(method = "lm",
                se = TRUE,
                color = "grey20",
                alpha = 0.4,
                lwd = 0.3) +
    geom_hline(yintercept = 0,
               linetype = 'dashed',
               color = 'grey50',
               lwd = 0.3) +
    geom_vline(xintercept = 0,
               linetype = 'dashed',
               color = 'grey50',
               lwd = 0.3) +
    geom_text_repel(data = subset(geneset_weights_df, color_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    size = 1.5,
                    fontface = "italic",
                    box.padding = 0.25,
                    point.padding = 0.15,
                    aes(x = z_score_z3,
                        y = z_score_z14,
                        label = variable)) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 label.x.npc = 0.8,
                 label.y.npc = 0.88,
                 formula = formula,
                 parse = TRUE,
                 size = 2,
                 na.rm = TRUE,
                 rr.digits = 1) +
    stat_fit_glance(method = "lm",
                    geom = "text",
                    label.x.npc = 0.8,
                    label.y.npc = 0.97,
                    method.args = list(formula = formula),
                    size = 2,
                    aes(label = paste("p = ",
                                      signif(..p.value.., digits = 1),
                                      sep = ""))) +
    xlab("Z Score for VAE z = 3 (Feature 0)") +
    ylab("Z Score for VAE z = 14\n(Feature 10)") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))

panel_e_gg

# Load and Process Data
file <- file.path('results', 'latent_feature_enrichment_comparison_neutrophil_genes.tsv')
gene_weights_df <- (
    readr::read_tsv(file,
                    col_types = readr::cols(.default = readr::col_double(),
                                            classification = readr::col_character(),
                                            gene = readr::col_character(),
                                            gene_set = readr::col_character()))
    )

geneset_classes <- c('Neutrophils',
                     'Monocytes',
                     'Keratinocytes',
                     'Skeletal Muscle',
                     'Neurons',
                     'Other Geneset',
                     'No Geneset')

gene_weights_df$classification <- factor(gene_weights_df$classification,
                                         levels = geneset_classes,
                                         ordered = TRUE)

head(gene_weights_df, 3)

color_labels <- c("Neutrophils" = "#1b9e77",
                  "Keratinocytes" = "#6a3d9a",
                  "Neurons" = "#1f78b4",
                  "Skeletal Muscle" = "#ff7f00",
                  "Monocytes" = "#e7298a",
                  "Other Geneset" = "#CFDEDA",
                  "No Geneset" = "#F5B8D4")

other_points_logic <- gene_weights_df$classification %in% c('Other Geneset', 'No Geneset')

panel_f_gg <- ggplot(gene_weights_df,
                     aes(x = vae_0_3,
                         y = vae_10)) +
    geom_point(data = subset(gene_weights_df, other_points_logic),
               aes(color = classification),
               shape = 16,
               size = 0.02,
               alpha = 0.4) +
    geom_point(data = subset(gene_weights_df, !other_points_logic),
               aes(color = classification),
               shape = 16,
               size = 0.2,
               alpha = 0.6) +
    geom_smooth(method = "lm",
                se = TRUE,
                color = "grey20",
                alpha = 0.4,
                lwd = 0.3) +
    geom_hline(yintercept = 0,
               linetype = 'dashed',
               color = 'grey50',
               lwd = 0.3) +
    geom_vline(xintercept = 0,
               linetype = 'dashed',
               color = 'grey50',
               lwd = 0.3) +
    stat_poly_eq(data = gene_weights_df, 
                 aes(label = paste(..rr.label..)),
                 label.x.npc = 0.2,
                 label.y.npc = 0.58,
                 formula = formula,
                 parse = TRUE,
                 size = 2,
                 na.rm = TRUE,
                 rr.digits = 1) +
    stat_fit_glance(method = "lm",
                    geom = "text",
                    label.x.npc = 0.2,
                    label.y.npc = 0.67,
                    method.args = list(formula = formula),
                    size = 2,
                    aes(label = paste("p = ",
                                      signif(..p.value.., digits = 1),
                                      sep = ""))) +
    scale_color_manual(name = "Gene Set Class",
                       values = color_labels,
                       breaks = geneset_classes) +
    xlab("Z Score for VAE z = 3 (Feature 0)") +
    ylab("Z Score for VAE z = 14\n(Feature 10)") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 4.7),
          legend.title = element_text(size = 6),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.margin = margin(t = -3, r = 0, b = -3, l = -3)) +
    guides(color = guide_legend(keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch",
                                override.aes = list(size = 1.4,
                                                    alpha = 1)))

panel_f_gg

# Load and process data
file <- file.path('results', 'neutrophil_data_biobombe_results.tsv')
full_neutrophil_results_df <-
    readr::read_tsv(file,
                    col_types = readr::cols(.default = readr::col_double(),
                                            full_id = readr::col_character(),
                                            cell_line = readr::col_character(),
                                            treatment = readr::col_character(),
                                            day = readr::col_character()))

full_neutrophil_results_df <-
    reshape2::melt(full_neutrophil_results_df,
                   id.vars = c('full_id', 'cell_line', 'treatment', 'day'),
                   variable.vars = c('vae_10', 'vae_0'),
                   variable.name = 'feature',
                   value.name = 'score')

full_neutrophil_results_df$feature <-
    dplyr::recode_factor(full_neutrophil_results_df$feature,
                         'vae_0' = 'VAE z = 3 (Feature 0)',
                         'vae_10' = 'VAE z = 14 (Feature 10)',
                         .ordered = TRUE)

head(full_neutrophil_results_df, 3)

panel_g_gg <- ggplot(full_neutrophil_results_df,
                     aes(x = cell_line, y = score)) +
    geom_jitter(aes(color = treatment), width = 0.1, height = 0) +
    coord_flip() +
    facet_wrap(~feature, scales = "free") +
    scale_color_manual(name = "",
                       values = c("DMSO" = "#c40000",
                                  "DMSO+Nutridoma" = "#fca82a",
                                  "Not Differentiated" = "#695eff")) +
    xlab('') +
    ylab('BioBombe Score') +
    ggtitle("GSE103706") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 9),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"),
          strip.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = 'bottom',
          legend.text = element_text(size = 4.7),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0)) +
    guides(color = guide_legend(nrow = 1,
                                ncol = 5,
                                byrow = FALSE,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch",
                                override.aes = list(size = 1.4,
                                                    alpha = 1)))

panel_g_gg

# Load and process results
file <- file.path('results', 'hematopoietic_data_biobombe_results.tsv')
full_heme_results_df <-
    readr::read_tsv(file,
                    col_types = readr::cols(.default = readr::col_double(),
                                            cell_type = readr::col_character(),
                                            replicate = readr::col_character(),
                                            cell = readr::col_character(),
                                            cell_class = readr::col_character()))

full_heme_results_df <-
    reshape2::melt(full_heme_results_df,
                   id.vars = c('cell_type', 'replicate', 'cell', 'cell_class'),
                   variable.vars = c('vae_2', 'nmf_6'),
                   variable.name = 'feature',
                   value.name = 'score')

full_heme_results_df$feature <-
    dplyr::recode_factor(full_heme_results_df$feature,
                         'vae_2' = 'VAE z = 3 (Feature 2)',
                         'nmf_6' = 'NMF z = 200 (Feature 6)',
                         .ordered = TRUE)

head(full_heme_results_df, 3)

# Plot
n_colors <- length(unique(full_heme_results_df$cell_type))

panel_h_gg <- ggplot(full_heme_results_df,
                     aes(x = cell_type, y = score)) +
    geom_boxplot(aes(color = cell_type), outlier.alpha = 0, lwd = 0.2) +
    geom_jitter(aes(color = cell_type), width = 0.1, size = 0.5, alpha = 0.8) +
    facet_wrap(~ feature, scales = "free", ncol = 2) +
    coord_flip() +
    scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_colors)) +
    xlab('') +
    ylab('BioBombe Score') +
    ggtitle("GSE24759") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 9),
          strip.background = element_rect(colour = "black", fill = "#fdfff4"),
          strip.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 5),
          legend.position = 'none')

panel_h_gg

c_and_d_legend_gg <- cowplot::get_legend(panel_c_gg) 
c_and_d_gg <- cowplot::plot_grid(
    panel_c_gg + theme(legend.position = "none"),
    panel_d_gg + theme(legend.position = "none"),
    nrow = 2,
    labels = c("C", "D")
)

c_and_d_gg <- cowplot::plot_grid(
    c_and_d_gg,
    c_and_d_legend_gg,
    rel_widths = c(1, 0.15),
    ncol = 2
)

c_and_d_gg

a_b_c_and_d_gg <- cowplot::plot_grid(
    panel_a_gg,
    panel_b_gg,
    c_and_d_gg,
    ncol = 3,
    labels = c("A", "B", "")
)

a_b_c_and_d_gg

e_and_f_gg <- cowplot::plot_grid(
    panel_e_gg,
    panel_f_gg,
    ncol = 2,
    labels = c("E", "F"),
    rel_widths = c(0.7, 1)
)

e_and_f_gg

g_and_h_gg <- cowplot::plot_grid(
    panel_g_gg,
    panel_h_gg,
    ncol = 2,
    labels = c("G", "H"),
    rel_widths = c(0.7, 1)
)

g_and_h_gg

full_gg <- cowplot::plot_grid(
    a_b_c_and_d_gg,
    e_and_f_gg,
    g_and_h_gg,
    nrow = 3,
    rel_heights = c(1.15, 0.70, 1.05)
)

full_gg

for(extension in c('.png', '.pdf')) {
    gg_file <- paste0("gtex_biobombe_main_figure", extension)
    gg_file <- file.path("figures", gg_file)
    cowplot::save_plot(filename = gg_file,
                       plot = full_gg,
                       base_height = 9,
                       base_width = 8)
}
