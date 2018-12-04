
library(ggplot2)
library(ggrepel)
suppressPackageStartupMessages(library(dplyr))

file <- file.path('results', 'latent_feature_enrichment_comparison_neutrophil_genesets.tsv')
geneset_weights_df <- (
    readr::read_tsv(file,
                    col_types = readr::cols(.default = readr::col_double(),
                                            model_type_vae = readr::col_character(),
                                            variable = readr::col_character(),
                                            algorithm_vae = readr::col_character(),
                                            model_type_dae = readr::col_character(),
                                            algorithm_dae = readr::col_character()))
    )

head(geneset_weights_df)

color_logic <- ((geneset_weights_df$z_score_vae < -16 | geneset_weights_df$z_score_vae > 18) |
                (geneset_weights_df$z_score_dae < -17))

ggplot(geneset_weights_df, aes(x = z_score_vae, y = z_score_dae)) +
  geom_point(size = 0.2) +
  geom_text_repel(data = subset(geneset_weights_df, color_logic),
                  arrow = arrow(length = unit(0.01, "npc")),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  color = 'red',
                  size = 2,
                  fontface = "italic",
                  point.padding = 0.4,
                  aes(x = z_score_vae, y = z_score_dae, label = variable)) +
  xlab("Z Score (VAE)") +
  ylab("Z Score (DAE)") +
  theme_bw()

fig_file <- file.path("figures", "feature_comparison_genesets.png")
ggsave(fig_file, height = 4.5, width = 6, dpi = 500)

file <- file.path('results', 'latent_feature_enrichment_comparison_neutrophil_genes.tsv')
gene_weights_df <- (
    readr::read_tsv(file,
                    col_types = readr::cols(.default = readr::col_double(),
                                            classification = readr::col_character(),
                                            gene = readr::col_character(),
                                            gene_set = readr::col_character()))
    )

geneset_classes <- c('Neutrophils', 'Hematopoetic', 'Skeletal Muscle',
                     'Hepatocytes', 'Other Geneset', 'No Geneset')
gene_weights_df$classification <- factor(gene_weights_df$classification, geneset_classes)

head(gene_weights_df, 3)

other_df <- gene_weights_df %>% dplyr::filter(classification %in% c('Other Geneset', 'No Geneset'))
neut_df <- gene_weights_df %>% dplyr::filter(!(classification %in% c('Other Geneset', 'No Geneset')))

ggplot(gene_weights_df, aes(x = vae_1_3, y = dae_2_6, color = classification)) +
  geom_point(data = other_df, size = 0.05, alpha = 0.5) +
  geom_point(data = neut_df, size = 0.1) +
  scale_color_manual(name = "Gene Set Class",
                     values = c("Neutrophils" = "#FFB413",
                                "Hematopoetic" = "#FF0800",
                                "Skeletal Muscle" = "#0080FF",
                                "Hepatocytes" = "#009805",
                                "Other Geneset" = "#CFDEDA",
                                "No Geneset" = "#F5B8D4"),
                     labels = c("Neutrophils" = "Neutrophils",
                                "Hematopoetic" = "Hematopoetic",
                                "Skeletal Muscle" = "Skeletal Muscle",
                                "Hepatocytes" = "Hepatocytes",
                                "Other Geneset" = "Other Geneset",
                                "No Geneset" = "No Geneset")) +
  xlab("VAE Feature 1 (z=3)") +
  ylab("DAE Feature 2 (z=6)") +
  theme_bw()

fig_file <- file.path("figures", "feature_comparison_genes.png")
ggsave(fig_file, height = 4.5, width = 6, dpi = 500)
