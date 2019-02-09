library(dplyr)
library(ggplot2)

file <- file.path("results", "nbl_mycn_separation_target_t_test.tsv")
full_results_df <- readr::read_tsv(file)

top_results_df <- full_results_df %>%
  dplyr::group_by(algorithm, z_dim, signal) %>%
  dplyr::filter(neg_log_p == max(neg_log_p))

# Create factors for plotting
top_results_df$z_dim <-
  factor(top_results_df$z_dim,
         levels =
           sort(as.numeric(paste(unique(top_results_df$z_dim))))
  )

top_results_df$algorithm <-
  factor(top_results_df$algorithm,
         levels = c("pca", "ica", "nmf", "dae", "vae"))


ggplot(top_results_df,
       aes(x = z_dim,
           y = neg_log_p,
           color = algorithm,
           group = algorithm)) +
  geom_point(size = 0.5) +
  facet_wrap(~ signal, scales = "free_y") +
  geom_line(lwd = 0.2) +
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
  theme_bw() +
  ylab("-log10 P Value") +
  xlab("k Dimensions") +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.7, "lines"))


best_result <- top_results_df %>%
  dplyr::arrange(desc(neg_log_p)) %>%
  dplyr::ungroup() %>%
  dplyr::top_n(n = 1, wt = neg_log_p)

best_k_dim <- paste(best_result$z_dim)
best_seed <- paste(best_result$seed)
best_feature <- paste0(paste(best_result$algorithm), "_", paste(best_result$feature_num))
feature_file <- file.path("..",
                          "2.ensemble-z-analysis",
                          "results",
                          "TARGET_results",
                          "ensemble_z_matrices",
                          paste0("target_components_", best_k_dim),
                          paste0("model_", best_seed, "_z_matrix.tsv.gz"))

top_feature_df <- readr::read_tsv(feature_file) %>%
  dplyr::select(sample_id, !!best_feature, "vae_107")

# Load TARGET phenotype data
file <- file.path("..", "0.expression-download", "data", "2017-09-30-TARGET update harmonized.txt")
nbl_pheno_df <- readr::read_tsv(file)

full_target_ids <- c()
for (target_id in strsplit(top_feature_df$sample_id, "-")) {
  target_id_unlist <- unlist(target_id)
  full_target_ids <- c(full_target_ids, target_id_unlist[3])
}
top_feature_df$target_id <- full_target_ids
top_feature_df <- top_feature_df %>% dplyr::inner_join(nbl_pheno_df, by = c("target_id" = "usi"))


ggplot(top_feature_df, aes(y = vae_111, fill = `MYCN status`)) +
    geom_boxplot(alpha = 0.5)


weight_file <- file.path("..",
                         "2.ensemble-z-analysis",
                         "results",
                         "TARGET_results",
                         "ensemble_z_matrices",
                         paste0("target_components_", best_k_dim),
                         paste0("model_", best_seed, "_weight_matrix.tsv.gz"))

top_weight_df <- readr::read_tsv(weight_file)


ggplot(top_weight_df, aes(x = vae_111)) + geom_density()

top_weight_df %>% dplyr::select(gene_id, vae_111) %>% dplyr::arrange(desc(vae_111))






# Load GTEx phenotype data
file <- file.path("..", "0.expression-download", "download", "GTEx_v7_Annotations_SubjectPhenotypesDS.txt")
gtex_pheno_df <- readr::read_tsv(file)

full_gtex_ids <- c()
for (gtex_id in strsplit(top_feature_df$sample_id, "-")) {
  gtex_id_unlist <- unlist(gtex_id)
  sample_id <- paste0(gtex_id_unlist[1], "-", gtex_id_unlist[2])
  full_gtex_ids <- c(full_gtex_ids, sample_id)
}

top_feature_df$gtex_id <- full_gtex_ids
top_feature_df <- top_feature_df %>% dplyr::left_join(gtex_pheno_df, by = c("gtex_id" = "SUBJID"))

top_feature_df$SEX <- as.factor(top_feature_df$SEX)


ggplot(top_feature_df, aes(y=vae_108, fill = SEX, group=SEX)) + geom_boxplot(alpha = 0.5)


weight_file <- file.path("..",
                          "2.ensemble-z-analysis",
                          "results",
                          "GTEX_results",
                          "ensemble_z_matrices",
                          paste0("gtex_components_", best_k_dim),
                          paste0("model_", best_seed, "_weight_matrix.tsv.gz"))

top_weight_df <- readr::read_tsv(weight_file)


ggplot(top_weight_df, aes(y = vae_108)) + geom_boxplot()

top_weight_df %>% dplyr::select(gene_id, vae_108) %>% dplyr::arrange(desc(vae_108))
