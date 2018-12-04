
library(ggplot2)
library(ggrepel)
suppressPackageStartupMessages(library(dplyr))

results_file <- file.path('results', 'gtex_vae_example_interpret_compression.tsv')
interpret_data_df <- (
    readr::read_tsv(results_file,
                    col_types = readr::cols(.default = readr::col_character(),
                                            value = readr::col_double(),
                                            z_score = readr::col_double()))
    )

interpret_data_df$index <- factor(interpret_data_df$index,
                                  levels = c("vae_0_two", "vae_1_two", "vae_0_three",
                                             "vae_1_three", "vae_2"))

head(interpret_data_df)

color_logic <- (interpret_data_df$z_score > 20 | interpret_data_df$z_score < -20)

ggplot(interpret_data_df, aes(x=value, y=z_score)) +
  geom_point(aes(color=index), size = 0.2) +
  scale_color_manual(name = "VAE Features",
                     values = c("#db70ff", "#9b41f4", "#9af441", "#29b239",
                                "#277544"),
                     labels = c("vae_0_two" = "Model 1 - Feature 1",
                                "vae_1_two" = "Model 1 - Feature 2",
                                "vae_0_three" = "Model 2 - Feature 1",
                                "vae_1_three" = "Model 2 - Feature 2",
                                "vae_2" = "Model 2 - Feature 3")) +
  geom_text_repel(data = subset(interpret_data_df, color_logic),
                  arrow = arrow(length = unit(0.01, "npc")),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 2,
                  fontface = "italic",
                  point.padding = 0.1,
                  aes(x = value, y = z_score, label = variable)) +
  xlab("Raw Score") +
  ylab("Z Score") +
  theme_bw()

fig_file <- file.path("figures", "interpret_compression_vae_neutrophils.png")
ggsave(fig_file, height = 4.5, width = 6, dpi = 500)

results_file <- file.path('results', 'gtex_vae_example_overrepresentation.tsv')
overrep_data_df <- (
    readr::read_tsv(results_file,
                    col_types = readr::cols(.default = readr::col_integer(),
                                            X1 = readr::col_character(),
                                            feature = readr::col_character(),
                                            odds = readr::col_double(),
                                            pval = readr::col_double(),
                                            tailed = readr::col_character()))
    )

overrep_data_df$feature <- factor(overrep_data_df$feature,
                                  levels = c("vae_0_two", "vae_1_two", "vae_0_three",
                                             "vae_1_three", "vae_2"))

overrep_data_df$neg_log10_p <- -log10(overrep_data_df$pval)
head(overrep_data_df %>% dplyr::arrange(desc(neg_log10_p)))

color_logic <- overrep_data_df$neg_log10_p > 11

ggplot(overrep_data_df, aes(x = odds, y = neg_log10_p, shape = tailed)) +
  geom_point(aes(color=as.factor(feature)), size = 3) +
  scale_color_manual(name = "VAE Features",
                     values = c("#db70ff", "#9b41f4", "#9af441", "#29b239",
                                "#277544"),
                     labels = c("vae_0_two" = "Model 1 - Feature 1",
                                "vae_1_two" = "Model 1 - Feature 2",
                                "vae_0_three" = "Model 2 - Feature 1",
                                "vae_1_three" = "Model 2 - Feature 2",
                                "vae_2" = "Model 2 - Feature 3")) +
  scale_shape_manual(name = "Tail",
                     values = c("-", "+"),
                     labels = c("pos" = "Positive",
                                "neg" = "Negative")) +
  geom_text_repel(data = subset(overrep_data_df, color_logic),
                  arrow = arrow(length = unit(0.01, "npc")),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 2,
                  fontface = "italic",
                  point.padding = 0.3,
                  aes(x = odds, y = neg_log10_p, label = X1)) +
  xlab("Odds Ratio") +
  ylab("-log 10 p value") +
  theme_bw()

fig_file <- file.path("figures", "interpret_compression_vae_neutrophils_tailed.png")
ggsave(fig_file, height = 4, width = 6, dpi = 500)

combined_results_df <- overrep_data_df %>%
    dplyr::full_join(interpret_data_df, by = c('X1' = 'variable', 'feature' = 'index'))

head(combined_results_df)

color_logic <- (combined_results_df$z_score > 20 | combined_results_df$z_score < -20)

ggplot(combined_results_df,
       aes(x = neg_log10_p,
           y = z_score,
           shape = tailed)) +
  geom_point(aes(color = as.factor(feature)), size = 2) +
  scale_color_manual(name = "VAE Features",
                     values = c("#db70ff", "#9b41f4", "#9af441", "#29b239",
                                "#277544"),
                     labels = c("vae_0_two" = "Model 1 - Feature 1",
                                "vae_1_two" = "Model 1 - Feature 2",
                                "vae_0_three" = "Model 2 - Feature 1",
                                "vae_1_three" = "Model 2 - Feature 2",
                                "vae_2" = "Model 2 - Feature 3")) +
  scale_shape_manual(name = "Tail",
                     values = c("-", "+"),
                     labels = c("pos" = "Positive",
                                "neg" = "Negative")) +
  geom_text_repel(data = subset(combined_results_df, color_logic),
                  arrow = arrow(length = unit(0.01, "npc")),
                  segment.size = 0.3,
                  segment.alpha = 0.6,
                  size = 1.4,
                  fontface = "italic",
                  point.padding = 0.3,
                  aes(x = neg_log10_p, y = z_score, label = X1)) +
  ylab("Z Score (Matrix)") +
  xlab("-log 10 p value (Overrepresentation)") +
  theme_bw()

fig_file <- file.path("figures", "full_comparison_vae_neutrophils_tailed.png")
ggsave(fig_file, height = 4, width = 6, dpi = 500)
