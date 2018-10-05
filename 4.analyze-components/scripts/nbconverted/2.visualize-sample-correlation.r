
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Load helper functions
source(file.path("scripts", "util.R"))

# Create theme
correlation_theme <- theme(axis.text.x = element_text(angle = 90, size = 5),
                           plot.title = element_text(hjust = 0.5),
                           legend.text = element_text(size = 8),
                           legend.key.size = unit(0.7, 'lines'))

# Create list for data storage
sample_correlation_list <- list()

dataset_name <- "TARGET"

# 1) Load phenotype data
target_file <- file.path("..", "0.expression-download", "download", "TARGET_phenotype.gz")
target_pheno_df <- readr::read_tsv(target_file,
                                   col_types = readr::cols(
                                       .default = readr::col_character()))

colnames(target_pheno_df)[2] <- 'sample_type'

# 2) Load sample correlation data
sample_correlation_list[[dataset_name]] <- compile_sample_correlation(dataset_name = dataset_name)

# 3) Merge phenotype and sample correlation results
sample_correlation_list[[dataset_name]] <-
    sample_correlation_list[[dataset_name]] %>%
    dplyr::full_join(target_pheno_df, by = c("id" = "sample_id")) %>%
    na.omit

head(sample_correlation_list[[dataset_name]], 2)

dataset_name <- "TCGA"

# 1) Load phenotype data
tcga_file <- file.path("..", "0.expression-download", "data", "tcga_sample_identifiers.tsv")
tcga_pheno_df <- readr::read_tsv(tcga_file,
                                 col_types = readr::cols(
                                        .default = readr::col_character()))

colnames(tcga_pheno_df)[2] <- 'sample_class'
colnames(tcga_pheno_df)[3] <- 'sample_type'

# 2) Load sample correlation data
sample_correlation_list[[dataset_name]] <- compile_sample_correlation(dataset_name = dataset_name)

# 3) Merge phenotype and sample correlation results
sample_correlation_list[[dataset_name]] <-
    sample_correlation_list[[dataset_name]] %>%
    dplyr::full_join(tcga_pheno_df, by = c("id" = "sample_id")) %>%
    na.omit

head(sample_correlation_list[[dataset_name]], 2)

dataset_name <- "GTEX"

# 1) Load phenotype data
gtex_file <- file.path("..", "0.expression-download", "download",
                       "GTEx_v7_Annotations_SampleAttributesDS.txt")
gtex_pheno_df <- readr::read_tsv(gtex_file,
                                 col_types = readr::cols(
                                        .default = readr::col_character()))

colnames(gtex_pheno_df)[1] <- 'sample_id'
colnames(gtex_pheno_df)[6] <- 'sample_type'

# Subset gtex phenotype file for plotting
gtex_pheno_df <- gtex_pheno_df[, c('sample_id', 'sample_type')]

# 2) Load sample correlation data
sample_correlation_list[[dataset_name]] <- compile_sample_correlation(dataset_name = dataset_name)

# 3) Merge phenotype and sample correlation results
sample_correlation_list[[dataset_name]] <-
    sample_correlation_list[[dataset_name]] %>%
    dplyr::full_join(gtex_pheno_df, by = c("id" = "sample_id")) %>%
    na.omit

head(sample_correlation_list[[dataset_name]], 2)

for (dataset_name in names(sample_correlation_list)) {
    # Extract out the specific dataset correlation data
    sample_corr_df <- sample_correlation_list[[dataset_name]]
    
    # Loop through the datatype
    for (data_type in c("signal", "shuffled")) {
        
        # Loop over the correlation type
        for (correlation_type in c("pearson", "spearman")) {
            
            print(paste("Processing... ",
                        dataset_name, data_type, correlation_type))

            # Execute the plotting logic
            plot_sample_correlation(dataset_name = dataset_name,
                                    data_df = sample_corr_df,
                                    data_type = data_type,
                                    correlation_type = correlation_type,
                                    use_theme = correlation_theme,
                                    return_figures = FALSE)
        }
    }
}

