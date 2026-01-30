### load libraries
library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(Biostrings)

### set working directory
setwd(file.path("~", "Documents", "Github", "Shiny_app_Acomys"))

### Source split functions
source("./split_functions.R")

### load REF_seq data (df with MEME input and original TOGA alignments)

load(file = file.path("data_files", "MEME_results_REF_seq.RData"))
msa_df <- MEME_results_REF_seq
rm(MEME_results_REF_seq)

### Split REF_seq_rows
split_rows(msa_df,
           column_name = "genename",
           RData_basename = "MEME_REF_seq_",
           subdir_name = file.path("~", "Documents", "Github",
                                   "Shiny_app_Acomys",
                                   "Shiny_Acomys", "data",
                                   "REF_seq_data"))
