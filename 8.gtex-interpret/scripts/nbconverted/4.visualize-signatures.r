
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpmisc))

# Load custom plotting functions
util_file = file.path("..", "6.biobombe-projection", "scripts", "utils.R")
source(util_file)

set.seed(123)

# Setup plotting logic
vae_labels <- c("vae_0_two" = "k = 2 (0)",
                "vae_1_two" = "k = 2 (1)",
                "vae_0_three" = "k = 3 (0)",
                "vae_1_three" = "k = 3 (1)",
                "vae_2_three" = "k = 3 (2)")

vae_colors <- c("#c994c7",
                "#dd1c77",
                "#78c679",
                "#31a354",
                "#006837")

results_file <- file.path('results', 'gtex_vae_example_overrepresentation.tsv')
overrep_data_df <- (
    readr::read_tsv(results_file,
                    col_types = readr::cols(.default = readr::col_integer(),
                                            variable = readr::col_character(),
                                            full_feature = readr::col_character(),
                                            odds = readr::col_double(),
                                            pval = readr::col_double(),
                                            tailed = readr::col_character()))
    )

overrep_data_df$full_feature <- factor(overrep_data_df$full_feature,
                                       levels = c("vae_0_two",
                                                  "vae_1_two",
                                                  "vae_0_three",
                                                  "vae_1_three",
                                                  "vae_2_three"))

overrep_data_df$neg_log10_p <- -log10(overrep_data_df$pval)
head(overrep_data_df %>% dplyr::arrange(desc(neg_log10_p)))

color_logic <- overrep_data_df$neg_log10_p > 10

# Create Panel A for a Supplementary Figure
sup_panel_a_gg <- ggplot(overrep_data_df,
                         aes(x = odds,
                             y = neg_log10_p,
                             shape = tailed)) +
    geom_point(aes(color = as.factor(full_feature)),
                   size = 3) +
    scale_color_manual(name = "VAE Features",
                       values = vae_colors,
                       labels =  vae_labels) +
    scale_shape_manual(name = "Tail",
                       values = c("-", "+"),
                       labels = c("pos" = "Positive",
                                  "neg" = "Negative")) +
    geom_text_repel(data = subset(overrep_data_df, color_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    box.padding = 0.25,
                    point.padding = 0.15,
                    aes(x = odds,
                        y = neg_log10_p,
                        label = variable)) +
    xlab("Odds Ratio") +
    ylab("-log 10 p value") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = 'right',
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 4.7),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0)) +
    guides(color = guide_legend(order = 1,
                                nrow = 5,
                                ncol = 1,
                                byrow = FALSE,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch",
                                override.aes = list(size = 1.4,
                                                    alpha = 1)),
           shape = guide_legend(order = 2,
                                nrow = 2,
                                ncol = 1,
                                byrow = FALSE,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))

sup_panel_a_gg

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

head(interpret_data_df)

combined_results_df <- overrep_data_df %>%
    dplyr::full_join(interpret_data_df,
                     by = c('variable', 'full_feature'))

head(combined_results_df)

color_logic <- (combined_results_df$z_score > 14 |
                combined_results_df$z_score < -20 |
                combined_results_df$neg_log10_p > 13)

sup_panel_b_gg <- ggplot(combined_results_df,
                         aes(x = z_score,
                             y = neg_log10_p,
                             shape = tailed)) +
    geom_point(aes(color = as.factor(full_feature)),
               size = 2) +
    scale_color_manual(name = "VAE Features",
                       values = vae_colors,
                       labels =  vae_labels) +
    scale_shape_manual(name = "Tail",
                       values = c("-", "+"),
                       labels = c("pos" = "Positive",
                                  "neg" = "Negative")) +
    geom_text_repel(data = subset(combined_results_df, color_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    box.padding = 0.25,
                    point.padding = 0.15,
                    aes(x = z_score,
                        y = neg_log10_p,
                        label = variable)) +
    xlab("BioBombe Z Score") +
    ylab("-log 10 p value (Overrepresentation)") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = 'right',
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 4.7),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0)) +
    guides(color = guide_legend(order = 1,
                                nrow = 5,
                                ncol = 1,
                                byrow = FALSE,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch",
                                override.aes = list(size = 1.4,
                                                    alpha = 1)),
           shape = guide_legend(order = 2,
                                nrow = 2,
                                ncol = 1,
                                byrow = FALSE,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))

sup_panel_b_gg

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

sup_panel_c_gg <- ggplot(geneset_weights_df,
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
                    size = 1.2,
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
    xlab("Z Score for VAE k = 3 (Feature 0)") +
    ylab("Z Score for VAE k = 14\n(Feature 10)") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))

sup_panel_c_gg

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

sup_panel_d_gg <- ggplot(gene_weights_df,
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
    xlab("Z Score for VAE k = 3 (Feature 0)") +
    ylab("Z Score for VAE k = 14\n(Feature 10)") +
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

sup_panel_d_gg

legend_gg <- cowplot::get_legend(sup_panel_b_gg)

sup_panel_a_and_b_gg <- cowplot::plot_grid(
    sup_panel_a_gg +
        theme(legend.position = 'none'),
    sup_panel_b_gg +
        theme(legend.position = 'none'),
    nrow = 2,
    labels = c("a", "b"),
    rel_heights = c(1, 1)
)

sup_panel_a_and_b_gg <- cowplot::plot_grid(
    sup_panel_a_and_b_gg,
    legend_gg,
    rel_widths = c(1, 0.2),
    ncol = 2
)

sup_panel_c_and_d_gg <- cowplot::plot_grid(
    sup_panel_c_gg,
    sup_panel_d_gg,
    nrow = 2,
    labels = c("c", "d"),
    rel_widths = c(1, 1)
)

sup_gg <- cowplot::plot_grid(
    sup_panel_a_and_b_gg,
    sup_panel_c_and_d_gg,
    ncol = 2,
    labels = c("", ""),
    rel_widths = c(1, 0.8)
)

sup_gg

for(extension in c('.png', '.pdf')) {
    gg_file <- paste0("gtex_biobombe_supplementary_figure", extension)
    gg_file <- file.path("figures", gg_file)
    cowplot::save_plot(filename = gg_file,
                       plot = sup_gg,
                       base_height = 130,
                       base_width = 170,
                       units = "mm")
}

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
vae_labels <- c("vae_0_two" = "k = 2 (0)",
                "vae_1_two" = "k = 2 (1)",
                "vae_0_three" = "k = 3 (0)",
                "vae_1_three" = "k = 3 (1)",
                "vae_2_three" = "k = 3 (2)")

vae_colors <- c("#c994c7", "#dd1c77", "#78c679", "#31a354", "#006837")

color_logic <- (interpret_data_df$z_score > 13 |
                interpret_data_df$z_score < -18 |
                interpret_data_df$raw_score < -22)

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
                    box.padding = 0.47,
                    point.padding = 0.23,
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
    xlab("Mean BioBombe Z Score (VAE k = 3)") +
    ylab("Mean BioBombe Z Score (VAE k = 2)") +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
    
panel_b_gg

gene_set_dir <- file.path("..", "6.biobombe-projection",
                          "results", "gtex",
                          "gpxcell", "signal")
metaedge <- "GpXCELL"
dataset <- "GTEX"

# Compile all results (use this for plotting later)
biobombe_results_df <- get_biobombe_results(gene_set_dir = gene_set_dir)

line_plot_theme <-
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"),
          strip.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          axis.title.x = element_text(margin = margin(t = 0.1,
                                                      r = 0,
                                                      b = 0,
                                                      l = 0,
                                                      unit = 'cm')),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 4.7),
          legend.title = element_text(size = 6),
          legend.margin = margin(t = 0,
                                 r = 0,
                                 b = 0,
                                 l = 0),
          legend.box.margin = margin(t = -3,
                                     r = 0,
                                     b = -3,
                                     l = -3))

panel_c_gg <- plot_gene_set(gene_set = "Neutrophils_HPCA_2",
                            full_results_df = biobombe_results_df,
                            metaedge = metaedge,
                            dataset = dataset,
                            show_plot = FALSE,
                            shuffled = FALSE,
                            return_top = FALSE,
                            return_plot = TRUE)

panel_c_gg <- panel_c_gg +
    line_plot_theme +
    ylab("BioBombe Score") +
    theme(title = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.title.x = element_text(size = 6))

panel_c_gg

panel_d_gg <- plot_gene_set(gene_set = "Monocytes_FANTOM_2",
                            full_results_df = biobombe_results_df,
                            metaedge = metaedge,
                            dataset = dataset,
                            show_plot = FALSE,
                            shuffled = FALSE,
                            return_top = FALSE,
                            return_plot = TRUE)

panel_d_gg <- panel_d_gg +
    line_plot_theme +
    ylab("BioBombe Score") +
    theme(title = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.title.x = element_text(size = 6))

panel_d_gg

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
                         'vae_0' = 'VAE k = 3 (Feature 0)',
                         'vae_10' = 'VAE k = 14 (Feature 10)',
                         .ordered = TRUE)

head(full_neutrophil_results_df, 3)

panel_e_gg <- ggplot(full_neutrophil_results_df,
                     aes(x = cell_line,
                         y = score)) +
    geom_jitter(aes(color = treatment),
                width = 0.1,
                height = 0,
                size = 0.8,
                alpha = 0.8) +
    coord_flip() +
    facet_wrap(~feature,
               scales = "free") +
    scale_color_manual(name = "",
                       values = c("DMSO" = "#c40000",
                                  "DMSO+Nutridoma" = "#fca82a",
                                  "Not Differentiated" = "#695eff")) +
    xlab('') +
    ylab('BioBombe Score') +
    ggtitle("GSE103706") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"),
          strip.text = element_text(size = 4.7),
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

panel_e_gg

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
                         'vae_2' = 'VAE k = 3 (Feature 2)',
                         'nmf_6' = 'NMF k = 200 (Feature 6)',
                         .ordered = TRUE)

head(full_heme_results_df, 3)

# Plot
n_colors <- length(unique(full_heme_results_df$cell_type))

y_axis_colors <- rep("black", n_colors)
y_axis_colors[unique(full_heme_results_df$cell_type) == 'MONO2'] <- "red"

panel_f_gg <- ggplot(full_heme_results_df,
                     aes(x = cell_type,
                         y = score)) +
    geom_boxplot(aes(color = cell_type),
                 outlier.alpha = 0,
                 lwd = 0.2) +
    geom_jitter(aes(color = cell_type),
                width = 0.1,
                size = 0.2,
                alpha = 0.8) +
    facet_wrap(~ feature, scales = "free", ncol = 2) +
    coord_flip() +
    scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_colors)) +
    xlab('') +
    ylab('BioBombe Score') +
    ggtitle("GSE24759") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"),
          strip.text = element_text(size = 4.7),
          axis.title = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 3.5,
                                     color = y_axis_colors),
          legend.position = 'none')

panel_f_gg

cell_types <- c(
  "Neutrophil", "Monocyte"
)

line_plot <- list()
for (cell_type in cell_types) {
    file <- file.path("results", paste0("all_", tolower(cell_type),
                                        "_top_scores_and_separation.tsv"))
    final_results_df <- readr::read_tsv(file,
                                        col_types = readr::cols(
                                            .default = readr::col_integer(),
                                            algorithm = readr::col_character(),
                                            t_stat = readr::col_double(),
                                            t_p = readr::col_double(),
                                            neg_log_p = readr::col_double(),
                                            variable = readr::col_character(),
                                            model_type = readr::col_character(),
                                            value = readr::col_double(),
                                            z_score = readr::col_double(),
                                            abs_z_score = readr::col_double()))
    final_results_df$algorithm <- 
      dplyr::recode(final_results_df$algorithm,
                    "pca" = "PCA", "ica" = "ICA", "nmf" = "NMF", "dae" = "DAE",
                    "vae" = "VAE")
  
    # Create factors for plotting
    final_results_df$z <-
      factor(final_results_df$z,
             levels =
               sort(as.numeric(paste(unique(final_results_df$z))))
      )
  
    final_results_df$algorithm <-
      factor(final_results_df$algorithm,
             levels = c("PCA", "ICA", "NMF", "DAE", "VAE"))
  
    final_results_df <- final_results_df %>% dplyr::arrange(z)

    # Generate plot title
    if (cell_type == "Neutrophil") {
        gse <- "GSE103706"
    } else {
        gse <- "GSE24759"
    }
    plot_title <- paste(cell_type, "-", gse)

    line_plot[[cell_type]] <- # Plot and save to file
      ggplot(final_results_df,
             aes(x = z,
                 y = neg_log_p,
                 color = algorithm,
                 group = algorithm)) +
      geom_point(size = 0.1) +
      geom_smooth(lwd = 0.2,
                  formula = y ~ x,
                  method = "loess",
                  alpha = 0.2,
                  aes(fill = algorithm)) +
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
      theme_bw() +
      ggtitle(plot_title) +
      ylab("-log 10 P") +
      xlab("k Dimensions") +
      theme(axis.title.x = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            axis.text.x = element_text(angle = 90,
                                     size = 5),
            axis.text.y = element_text(size = 6),
            plot.title = element_text(hjust = 0.5,
                                    size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 7),
            legend.key.size = unit(0.7, "lines"))
}

# Plot
line_legend <- cowplot::get_legend(line_plot[['Neutrophil']])
panel_g_gg <- cowplot::plot_grid(
    line_plot[['Neutrophil']] + theme(legend.position = "none"),
    line_plot[['Monocyte']] + theme(legend.position = "none"),
    labels = c("g", ""),
    ncol = 2
)

panel_g_gg <- cowplot::plot_grid(
    panel_g_gg,
    line_legend,
    rel_widths = c(1, 0.15),
    ncol = 2
)

panel_g_gg

myPalette <- colorRampPalette(rev(c("#43b8d8",
                                    "#9e5eed",
                                    "#ef07c4",
                                    "#ef0707")))

career_path_plot <- list()
for (cell_type in cell_types) {
    file <- file.path("results", paste0("all_", tolower(cell_type),
                                      "_top_scores_and_separation.tsv"))
    final_results_df <- readr::read_tsv(file,
                                        col_types = readr::cols(
                                            .default = readr::col_integer(),
                                            algorithm = readr::col_character(),
                                            t_stat = readr::col_double(),
                                            t_p = readr::col_double(),
                                            neg_log_p = readr::col_double(),
                                            variable = readr::col_character(),
                                            model_type = readr::col_character(),
                                            value = readr::col_double(),
                                            z_score = readr::col_double(),
                                            abs_z_score = readr::col_double()))
    final_results_df$algorithm <- 
    dplyr::recode(final_results_df$algorithm,
                  "pca" = "PCA", "ica" = "ICA", "nmf" = "NMF", "dae" = "DAE",
                  "vae" = "VAE")

    # Create factors for plotting
    final_results_df$z <-
    factor(final_results_df$z,
           levels =
             sort(as.numeric(paste(unique(final_results_df$z))))
    )

    final_results_df$algorithm <-
    factor(final_results_df$algorithm,
           levels = c("PCA", "ICA", "NMF", "DAE", "VAE"))

    final_results_df <- final_results_df %>% dplyr::arrange(z)

    career_path_plot[[cell_type]] <-
    ggplot(final_results_df,
                          aes(y = neg_log_p,
                              x = abs_z_score,
                              group = algorithm)) +
    geom_point(size = 0.1, alpha = 0.9) +
    geom_path(aes(color = z_dim),
              lwd = 0.2) +
    facet_wrap(~algorithm, ncol = 5) +
    scale_color_gradientn(name = "k Dimension",
                          colours = myPalette(300),
                          values = scales::rescale(c(200, 125, 25, 10, 0)),
                          limits = c(0, 200)) +
    theme_bw() +
    ylab("-log10 P") +
    xlab("BioBombe Score") +
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"),
          strip.text.x = element_text(size = 5,
                                      margin = margin(t = 2,
                                                      b = 1,
                                                      l = 0,
                                                      r = 0)),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(angle = 90,
                                     size = 5),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          legend.key.size = unit(0.7, "lines"))
}

career_legend <- cowplot::get_legend(career_path_plot[['Neutrophil']])

panel_h_gg <- cowplot::plot_grid(
    career_path_plot[['Neutrophil']] + theme(legend.position = "none"),
    career_path_plot[['Monocyte']] + theme(legend.position = "none"),
    labels = c("h", ""),
    ncol = 2
)

panel_h_gg <- cowplot::plot_grid(
    panel_h_gg,
    career_legend,
    rel_widths = c(1, 0.15),
    ncol = 2
)

panel_h_gg

c_and_d_legend_gg <- cowplot::get_legend(panel_c_gg) 
c_and_d_gg <- cowplot::plot_grid(
    panel_c_gg + theme(legend.position = "none"),
    panel_d_gg + theme(legend.position = "none"),
    nrow = 2,
    labels = c("c", "d")
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
    labels = c("a", "b", ""),
    rel_widths = c(0.95, 0.95, 1)
)

a_b_c_and_d_gg

e_and_f_gg <- cowplot::plot_grid(
    panel_e_gg,
    panel_f_gg,
    ncol = 2,
    labels = c("e", "f"),
    rel_widths = c(0.7, 1)
)

e_and_f_gg

g_and_h_gg <- cowplot::plot_grid(
    panel_g_gg,
    panel_h_gg,
    nrow = 2,
    rel_heights = c(1.1, 0.9)
)

g_and_h_gg

full_gg <- cowplot::plot_grid(
    a_b_c_and_d_gg,
    e_and_f_gg,
    g_and_h_gg,
    nrow = 3,
    rel_heights = c(1, 1.1, 1)
)

full_gg

for(extension in c('.png', '.pdf')) {
    gg_file <- paste0("gtex_biobombe_main_figure", extension)
    gg_file <- file.path("figures", gg_file)
    cowplot::save_plot(filename = gg_file,
                       plot = full_gg,
                       base_height = 200,
                       base_width = 170,
                       units = "mm")
}
