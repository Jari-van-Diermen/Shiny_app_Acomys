### load libraries
library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(Biostrings)

### set working directory
setwd(file.path("~", "Documents", "Github", "Shiny_app_Acomys"))

### load REF_seq data (df with MEME input and original TOGA alignments)
uniprot_seqs <- read_delim("data_files/uniprot_seqs.tsv",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)
# Keep identifier rows only
uniprot_seqs <- uniprot_seqs %>%
  dplyr::select(transcript_id, genename, uniprotswissprot)

# Store in R shiny data folder
write_path <- file.path("~", "Documents", "Github", "Shiny_app_Acomys",
                        "Shiny_Acomys", "data", "Uniprot", "uniprot_ids.tsv")

# Write uniprot data
dir.create(dirname(write_path))
write_delim(uniprot_seqs, file = write_path, delim = "\t")
