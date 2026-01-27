### load libraries
library(dplyr)
library(readr)
library(purrr)
library(magrittr)

### set working directory
setwd(file.path("~", "Documents", "Github", "Shiny_app_Acomys"))

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
load(file = file.path("data_files", "MEME_substitutions_small.RData"))

# # EBF data calculated from the MEME data
# load(file = file.path("data_files", "MEME_EBF_values_simple.RData"))

# Uniprot-mapped EBF data from MEME results for bubbleplot
load(file = file.path("data_files", "EBF_data_uniprot_bubbleplot.RData"))
EBF_data_unimapped_bubble <- EBF_data_uniprot_noMSA
rm(EBF_data_uniprot_noMSA)

# Uniprot-mapped EBF data from MEME results for sign sites
load(file = file.path("data_files", "EBF_data_uniprot.RData"))
EBF_data_unimapped <- EBF_data_mapped_noMSA
rm(EBF_data_mapped_noMSA)

# rename data
MEME_data <- MEME_simple_nsites_subs_small
rm(MEME_simple_nsites_subs_small)

# # Combine EBF_data with MEME_data 
# MEME_data %<>%
#   left_join(MEME_EBF_simple)

### Split the data in multiple files on a per row basis

# Function that takes the complete MEME or aBSREL data object, and creates a
# new RData object per row
#
# args:
#   RData_object: An RData object to split into multiple files using the rows
#   column_name: the column name of the column that will be used for the
#   splitting procedure
#   RData_basename: The basename that the separate files should have
#   subdir_name: The subdirectory that is created to store the generated RData
#   files
split_rows <- function(RData_object, column_name = "genename", 
                       RData_basename = "MEME_data_",
                       subdir_name = file.path("~", "Documents", "Github",
                                               "Shiny_app_Acomys",
                                               "Shiny_Acomys", "data",
                                               "MEME_data")) {
  column_vector <- RData_object[[column_name]]
  
  dir.create(subdir_name)
  
  # Loop through the rows and save each one as a separate RData file
  for (row_id in seq_along(RData_object[[column_name]])) {
    
    RData_row <- RData_object[row_id,]
    RData_genename <- RData_row[[column_name]]
    
    save(RData_row, file = file.path(subdir_name, paste0(RData_basename, RData_genename, ".RData")))
    cat(paste0("Saved ", RData_genename, " row as a separate RData object\n\n"))
  }
}

# get separate files for MEME data rows and unimapped-EBF values
split_rows(MEME_data)

split_rows(EBF_data_unimapped_bubble, RData_basename = "EBF_data_up_bubble_",
           subdir_name = file.path("~", "Documents", "Github",
                                   "Shiny_app_Acomys",
                                   "Shiny_Acomys", "data",
                                   "EBF_data_uniprot_bubbleplot"))

split_rows(EBF_data_unimapped, RData_basename = "EBF_data_up_",
           subdir_name = file.path("~", "Documents", "Github",
                                   "Shiny_app_Acomys",
                                   "Shiny_Acomys", "data",
                                   "EBF_data_uniprot"))

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
