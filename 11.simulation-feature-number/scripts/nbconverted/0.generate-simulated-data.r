suppressPackageStartupMessages(library(dplyr))

set.seed(1234)

n = 10000
p = 10

cov_mat = diag(p)

random_off_diag_structure <- abs(rnorm(n = length(cov_mat[lower.tri(cov_mat)]), mean = 0, sd = 0))

cov_mat[lower.tri(cov_mat)] <- random_off_diag_structure

cov_mat[2, 1] <- 0.95
cov_mat[3, 2] <- 0.90
cov_mat[3, 1] <- 0.93

cov_mat[6, 5] <- 0.90
cov_mat[7, 6] <- 0.85
cov_mat[7, 5] <- 0.88

cov_mat[upper.tri(cov_mat)] <- t(cov_mat)[upper.tri(cov_mat)]

cov_mat

cov_mat %>%
    dplyr::as_tibble(.name_repair = "minimal")

feature_ids <- paste0("feature_", seq(1, nrow(cov_mat)))

cov_mat_df <- cov_mat %>%
    dplyr::as_tibble(.name_repair = "minimal")

colnames(cov_mat_df) <- feature_ids

cov_mat_df <- cov_mat_df %>%
    dplyr::mutate(feature_num = feature_ids) %>%
    dplyr::select(feature_num, dplyr::everything())

out_file <- file.path("data", "simulated_covariance_structure.tsv")
cov_mat_df %>% readr::write_tsv(out_file)

cov_mat_df

simulated_data <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cov_mat) 
colnames(simulated_data) <- paste0("feature_", 1:ncol(simulated_data))
simulated_data <- simulated_data %>% dplyr::as_tibble(.name_repair = "minimal")

print(dim(simulated_data))
head(simulated_data)

out_file <- file.path("data", "simulated_signal_n1000_p10.tsv")
simulated_data %>% readr::write_tsv(out_file)
