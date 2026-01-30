### MEME backend for 'Multiple sequence alignment' tabpanel

# Get MEME_results (that do not change with the buttons in the sidebar)
get_MEME_results_MSA <- reactive({
  
  # Stop if no gene selected
  req(input$MEMEGeneInput)
  
  MEME_results <- get_MEME_data_row() %>%
    use_series(csv) %>%
    extract2(1) %>%
    # add an additional site column
    add_column(TOGA_Human_site = seq(nrow(.))) %>%
    # Join with translation df
    dplyr::right_join(get_gene_site_translations(), by = join_by(TOGA_Human_site)) %>%
    dplyr::relocate(MSA_site, TOGA_Human_site, TOGA_Acomys_site, UP_Human_site)
  
  # Needed to retrieve all non TOGA_Human_sites again
  MEME_results <- get_gene_site_translations() %>%
    dplyr::left_join(MEME_results, by = join_by(MSA_site)) %>%
    dplyr::select(MSA_site, `p-value`)
  
  MEME_results
})

# Load the dataframe with MSAs for gene
get_MSA_data_row <- reactive({
  
  # Stop if no gene selected
  req(input$MEMEGeneInput)
  
  # load the required RData object
  load(file.path("data", "REF_seq_data",
                 paste0("MEME_REF_seq_", get_genename_MEME(), ".RData")))
  return(RData_row)
})


# Load the data with EBF values for gene
get_EBF_data_row <- reactive({
  
  # Stop if no gene selected
  req(input$MEMEGeneInput)
  
  load(file.path("data", "EBF_data",
                 paste0("EBF_data_", get_genename_MEME(), ".RData"))) 
  return(RData_row)
})

# Get EBF dataframe from loaded data
get_EBF_dataframe <- reactive({
  
  EBF_data <- get_EBF_data_row()
  
  # Get EBF table
  EBF_table <- EBF_data %>%
    dplyr::pull(EBF_table) %>%
    magrittr::extract2(1) %>%
    dplyr::rename(TOGA_Human_site = site) %>%
    # Add other reference frames.
    tidyr::pivot_wider(names_from = branch, values_from = EBF) %>%
    dplyr::right_join(get_gene_site_translations() %>%
                        dplyr::mutate(TOGA_Human_site = as.integer(TOGA_Human_site)),
                      by = join_by(TOGA_Human_site)) %>%
    dplyr::relocate(MSA_site, TOGA_Human_site, TOGA_Acomys_site, UP_Human_site, REFERENCE) %>%
    tidyr::pivot_longer(REFERENCE:last_col(), names_to = "branch", values_to = "EBF") %>%
    # Sort by MSA_site and branch
    dplyr::arrange(branch, MSA_site)
  
  EBF_table
})

# Construct the MSA visualization
create_MSA_viz <- observe({
  
  ## Stop if no gene selected
  req(input$MEMEGeneInput)
  
  # File paths
  final_pdf <- file.path("www", "MSA_pdfs", paste0(get_genename_MEME(), "_alignment.pdf"))
  tmp_pdf <- file.path("www", "MSA_pdfs", paste0(get_genename_MEME(), "_alignment.tmp.pdf"))
  
  # Check if file already exists
  if (file.exists(final_pdf)) {
    return(NULL)
  }
  
  ## Get the TOGA MSA sequence
  MSA_data <- get_MSA_data_row()
  orig_seqs_row_seq <- MSA_data$original_sequence$misc_data
  
  ## limit alignment to Humans and rodents
  # get Boolean vector describing which sequences to include in MSA 
  seq_names <- names(orig_seqs_row_seq)
  incl_bool <- sapply(seq_names, function(id) {
    if (any(str_detect(id, included_assemblies))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  # filter MSA using boolean vector
  orig_seqs_row_seq <- orig_seqs_row_seq[incl_bool]
  
  ## Change assembly names to more interpretable species names
  names(orig_seqs_row_seq) <- sapply(names(orig_seqs_row_seq), function(id) {
    if (id == "REFERENCE") {
      return("Homo_sapiens")
    } else {
      # Retrieve assembly name from id
      assembly_name <- stringr::str_match(id, "vs_([[:alnum:]]+)\\t.*")[2]
      
      # Use this assembly name to retrieve species name
      species_name <- name_conversion %>%
        filter(`Assembly name` == assembly_name) %>%
        use_series(Species_Tree)
      
      # Replace assembly name with species name in id
      id_species_name <- stringr::str_replace(id, "vs_[[:alnum:]]+(\\t.*)",
                                              paste0(species_name, "\\1"))
      
      return(id_species_name)
    }
  })
  
  ## Get sequences of the MSA in the desired order
  
  # Assuming the 'included_assemblies' vector is in the correct order,
  # get vector of species names in the desired order.
  ordered_species_names <- sapply(included_assemblies, function(id) {
    if (id == "REFERENCE") {
      return("Homo_sapiens")
    } else {
      species_name <- name_conversion %>%
        filter(`Assembly name` == id) %>%
        use_series(Species_Tree)
      return(species_name)
    }
  })
  
  # Remove duplicate species names
  ordered_species_names <- ordered_species_names[!duplicated(ordered_species_names)]
  
  # Retrieve vector of species names from the sequence ids
  msa_species_names <- sapply(names(orig_seqs_row_seq), function(id) {
    if (id == "Homo_sapiens") {
      return(id)
    } else {
      # Get species name from id
      species_name <- stringr::str_match(id, "(\\w+)\\t.*")[2]
      return(species_name)
    }
  })
  
  # Retrieve the correct sequence ordering based on the order of
  # 'ordered_species_names'
  MSA_names_ordered <- names(orig_seqs_row_seq)[order(match(msa_species_names, ordered_species_names))]
  # Reorder MSA
  orig_seqs_row_seq <- orig_seqs_row_seq[MSA_names_ordered]
  
  ## translate the sequence
  orig_seqs_row_seq_aa <- suppressWarnings(translateGappedAln(orig_seqs_row_seq))
  
  ## Remove Gap-only positions in the species-filtered MSA
  
  # Calculate consensus matrix for each position
  cm_aa <- consensusMatrix(orig_seqs_row_seq_aa)
  
  # Find all positions that are gap-only
  pos_gaponly <- which(cm_aa[rownames(cm_aa) == "-",] == length(orig_seqs_row_seq_aa))
  
  # Transform gap-only positions to an IntegerRanges object
  IR_gaponly <- IRanges::reduce(IRanges::IRanges(start = pos_gaponly, end = pos_gaponly))
  
  # Filter out all gap-only positions
  seq_aa_gapfiltered <- Biostrings::replaceAt(orig_seqs_row_seq_aa, at = IR_gaponly)
  
  ## Remove duplicate transcript names
  seq_aa_gapfiltered <- seq_aa_gapfiltered[!duplicated(names(seq_aa_gapfiltered))]
  
  ## Replace the tabs in the sequence names with a single space character
  names(seq_aa_gapfiltered) <- gsub("\t", " ", names(seq_aa_gapfiltered))
  
  ## Get MEME P-values for gene
  site_pvals <- get_MEME_results_MSA() %>%
    dplyr::arrange(MSA_site)
  
  ## Get EBF values for A. cahirinus
  EBF_data <- get_EBF_dataframe() %>%
    dplyr::filter(branch == "HLacoCah2") %>%
    dplyr::select(MSA_site, EBF)
  
  # Get max non-inf EBF value
  max_EBF <- max(EBF_data$EBF[!is.infinite(EBF_data$EBF)], na.rm = TRUE)
  
  ## Prepare annotation used for the ggmsa visualization
  
  anno_data <- EBF_data %>%
    # Transform infinite values to highest non-inf value
    dplyr::mutate(EBF = if_else(is.infinite(EBF), max_EBF, EBF)) %>%
    # Join with MEME P-value data
    dplyr::right_join(site_pvals, by = join_by(MSA_site))
  
  # Use tidy_msa() to get a sequence conservation table. This will be rebuild
  # into an annotation dataframe compatible with ggmsa
  pval_annotation <- bar_data(ggmsa::tidy_msa(seq_aa_gapfiltered))
  
  pval_annotation <- anno_data %>%
    dplyr::rename(pos = MSA_site, p_value = `p-value`) %>%
    dplyr::right_join(pval_annotation, by = join_by(pos)) %>%
    dplyr::arrange(pos) %>%
    dplyr::mutate("-log10(p-value)" = -log10(p_value)) %>%
    dplyr::mutate(sign = p_value <= 0.05)
  
  ## Generate annotated alignment with the help of the `ggmsa` package
  
  # Calculate the needed parameters
  length_msa <- width(seq_aa_gapfiltered)[1]
  n_of_rows <- floor(length_msa / 75)
  len_last_row <- length_msa %% 75
  number_of_seqs <- length(seq_aa_gapfiltered)
  
  current_gene <- input$MEMEGeneInput
  
  withProgress(
    message = "Generating multiple sequence alignment PDF",
    value = 0,
    {
      # Create a PDF file with the annotated alignment
      pdf(file = tmp_pdf, height = 1.5+(0.08*number_of_seqs), width = 9)
      
      ## Use while loop to plot 75-length sections on different PDF pages.
      i = 0
      while (i < n_of_rows) {
        
        # Increment progress bar
        incProgress(
          amount = 1 / n_of_rows,
          detail = paste("Processing page", i, "of", n_of_rows)
        )
        
        pos_start <- i*75 + 1
        pos_end <- i*75 + 75
        
        # Build seq_length_df
        seq_length_df <- build_seq_lengths_df(seq_aa_gapfiltered,
                                              pos_end = pos_end)
        # Generate annotated MSA
        anno_p <- 
          suppressWarnings(
            suppressMessages(
              plot_annotated_msa_EBF(seq_aa_gapfiltered,
                                     start = pos_start,
                                     end = pos_end,
                                     # +3 added so length labels fit
                                     plot_end = (pos_end + 3),
                                     anno_df = pval_annotation,
                                     seq_length_df = seq_length_df,
                                     x_breaks_interval = 5,
                                     rel_heights = c(10, 0.8, 0.8, 0.8))
            )
          )
        
        print(anno_p)
        
        i <- i + 1
      }
      ## Generate last section of alignment
      
      incProgress(
        amount = 1 / n_of_rows,
        detail = paste("Processing page", i, "of", n_of_rows)
      )
      
      # If alignment is shorter than 75.
      if (i == 0) {
        pos_start <- 1
        pos_end <- len_last_row
        
        seq_length_df <- build_seq_lengths_df(seq_aa_gapfiltered,
                                              pos_end = pos_end)
        
        anno_p <- 
          suppressWarnings(
            suppressMessages(
              plot_annotated_msa_EBF(seq_aa_gapfiltered,
                                     start = pos_start,
                                     end = pos_end,
                                     # +3 added so length labels fit
                                     plot_end = (pos_start + 75 + 3),
                                     x_breaks_interval = 5,
                                     anno_df = pval_annotation,
                                     seq_length_df = seq_length_df,
                                     rel_heights = c(10, 0.8, 0.8, 0.8))
            )
          )
      } else {
        # If last section of a longer alignment
        seq_length_df <- build_seq_lengths_df(seq_aa_gapfiltered,
                                              pos_end = (pos_end + len_last_row))
        
        anno_p <- 
          suppressWarnings(
            suppressMessages(
              plot_annotated_msa_EBF(seq_aa_gapfiltered,
                                     start = pos_end,
                                     end = (pos_end + len_last_row),
                                     # +3 added so length labels fit
                                     plot_end = (pos_end + 75 + 3),
                                     x_breaks_interval = 5,
                                     anno_df = pval_annotation,
                                     seq_length_df = seq_length_df,
                                     rel_heights = c(10, 0.8, 0.8, 0.8))
            )
          )
      }
      
      print(anno_p)
      
    }
  )
  
  dev.off()
  
  # Change tmp file to final pdf
  file.rename(tmp_pdf, final_pdf)
})


output$MultipleAlignmentText <- renderUI({
  if (get_genename_MEME() == "") {
    return("Please select a gene to display the multiple sequence alignment")
  }
  return(NULL)
})

# Render the multiple sequence alignment
output$MultipleAlignment <- renderUI({
  
  # Stop if no gene selected
  req(input$MEMEGeneInput)
  
  # create URL
  gene_pdf <- file.path("MSA_pdfs", paste0(get_genename_MEME(), "_alignment.pdf"))
  
  # create pdf iframe
  tags$iframe(style="height:900px; width:100%; scrolling=yes", 
              src = gene_pdf)
})
