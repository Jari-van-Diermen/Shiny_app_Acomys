### load libraries
library(dplyr)
library(readr)
library(purrr)
library(magrittr)

### set working directory
setwd(file.path("~", "Documents", "Github", "Shiny_app_Acomys"))

### Source split functions
source("./split_functions.R")

### load and rename aBSREL data

load(file = file.path("data_files", "aBSREL_analysis_pval_count.RData"))
aBSREL_data <- results_tibble_count_cor
rm(results_tibble_count_cor)

load(file = file.path("data_files", "aBSREL_analysis_pval_count_bias.RData"))
aBSREL_data_bias <- results_tibble_count_cor
rm(results_tibble_count_cor)

### Load aBSREL-estimated dN/dS data

dNdS_df <- read_delim(file.path("data_files", "aBSREL_dNdS_df.tsv"),
                      delim = "\t", escape_double = FALSE,
                      col_types = cols(rate_class_number = col_integer(),
                                       omega_3 = col_double(), omega_3_prop = col_double()),
                      trim_ws = TRUE)

### Load MEME data

# MEME data with substitution columns
load(file = file.path("data_files", "MEME_substitutions.RData"))
MEME_data <- MEME_simple_nsites_subs
rm(MEME_simple_nsites_subs)

# Create an input tree column
MEME_data <- MEME_data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(tree = misc_data %>%
                  magrittr::extract2("input") %>%
                  magrittr::extract2("trees") %>%
                  magrittr::extract2("0") %>%
                    paste0(";")
                  ) %>%
  dplyr::ungroup()

# Simplify the data by removing unneccesary columns
MEME_data <- MEME_data %>%
  dplyr::select(-branch_attributes, -substitutions, -post_prob_omega_class, -misc_data)

# EBF data calculated from the MEME data
# Object called 'MEME_EBF_simple'
load(file = file.path("data_files", "MEME_EBF_values_simple.RData"))

### Split the data in multiple files on a per row basis

## Prepare the EBF data

# Filter EBF-data for selected human and rodent assemblies

included_assemblies <- c("REFERENCE", "HLacoCah2", "HLpsaObe1",
                         "HLmerUng1", "HLratNor7", "rn6", "HLmusPah1",
                         "HLmusCar1", "mm10", "mm39", "HLmesAur2",
                         "mesAur1", "HLcriGri3", "HLsigHis1", "HLonyTor1",
                         "HLperManBai2", "HLondZib1", "HLellLut1")

MEME_EBF_simple <- MEME_EBF_simple %>%
  mutate(EBF_table = purrr:::map(EBF_table, ~ {
    EBF_df_filt <- .x %>%
      dplyr::filter(branch %in% included_assemblies)
    return(EBF_df_filt)
  }))

# Split EBF data
split_rows(MEME_EBF_simple, RData_basename = "EBF_data_",
           subdir_name = file.path("~", "Documents", "Github",
                                   "Shiny_app_Acomys",
                                   "Shiny_Acomys", "data",
                                   "EBF_data"))

# get separate files for MEME data rows and unimapped-EBF values
split_rows(MEME_data)

# Combine aBSREL_data and aBSREL_data_bias into one dataframe

for (gene in unique(aBSREL_data_bias$genename)) {
  
  # If no duplicates, replace aBSREL_data generow with aBSREL_data_bias_generow
  if (nrow(aBSREL_data[aBSREL_data$genename == gene,]) == 1 &&
      nrow(aBSREL_data_bias[aBSREL_data_bias$genename == gene,]) == 1) {
    
    aBSREL_data[aBSREL_data$genename == gene,] <- aBSREL_data_bias[aBSREL_data_bias$genename == gene,]
  }
}

# simplify the aBSREL RData
aBSREL_data <- aBSREL_data %>%
  select(-all_of(c("omega_rates", "tree", "TRUE_count", "FALSE_count")))

save(aBSREL_data, file = file.path("~", "Documents", "Github",
                                   "Shiny_app_Acomys",
                                   "Shiny_Acomys", "data",
                                   "aBSREL_data.RData"))

# filter dN/dS data for Acomys data and save dN/dS data
dNdS_df <- dNdS_df %>%
  dplyr::filter(species == "HLacoCah2")

write_delim(dNdS_df, file = file.path("~", "Documents", "Github",
                                      "Shiny_app_Acomys",
                                      "Shiny_Acomys", "data",
                                      "aBSREL_acomys_dNdS_df.tsv"),
            delim = "\t")
