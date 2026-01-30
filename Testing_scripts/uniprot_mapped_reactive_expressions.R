get_unimapped_MEME_results <- reactive({

  # Stop if no gene selected
  req(input$MEMEGeneInput)

  MEME_UP_results <- MEME_subs %>%
    dplyr::filter(genename == get_genename_MEME()) %>%
    #dplyr::filter(genename == "IL12RB1") %>%
    separate_longer_delim(c(sites, TOGA_sites, branch_subs, branch_subs_ancestral, `&alpha;`,
                            `&beta;<sup>-</sup>`, `p<sup>-</sup>`,
                            `&beta;<sup>+</sup>`, `p<sup>+</sup>`, `LRT`, `p-value`,
                            pval_fdr, `#_branches_under_selection`, `Total_branch_length`,
                            `MEME_LogL`, `FEL_LogL`, `Variation_p`),
                          delim = stringr::regex("[;,]")) %>%
    # rename columns
    dplyr::rename(UP_Human_site = sites,
                  TOGA_Human_site = TOGA_sites) %>%
    ## format some columns as doubles and integers
    mutate(UP_Human_site = as.integer(UP_Human_site),
           TOGA_Human_site = as.integer(TOGA_Human_site),
           `#_branches_under_selection` = as.integer(`#_branches_under_selection`),
           `&alpha;` = as.double(`&alpha;`),
           `&beta;<sup>-</sup>` = as.double(`&beta;<sup>-</sup>`),
           `p<sup>-</sup>` = as.double(`p<sup>-</sup>`),
           `&beta;<sup>+</sup>` = as.double(`&beta;<sup>+</sup>`),
           `p<sup>+</sup>` = as.double(`p<sup>+</sup>`),
           `LRT` = as.double(`LRT`),
           `p-value` = as.double(`p-value`),
           `pval_fdr` = as.double(`pval_fdr`),
           `Total_branch_length` = as.double(`Total_branch_length`),
           `MEME_LogL` = as.double(`MEME_LogL`),
           `FEL_LogL` = as.double(`FEL_LogL`),
           `Variation_p` = as.double(`Variation_p`)) %>%
    # Combine with translation table
    dplyr::left_join(get_gene_site_translations(), by = join_by(UP_Human_site, TOGA_Human_site)) %>%
    #dplyr::left_join(gene_trans_df, by = join_by(UP_Human_site, TOGA_Human_site)) %>%
    dplyr::select(-all_of(c("site_included",
                            "uniprot_gn_symbol", "ensembl_gene_id",
                            "pval_fdr", "n_of_sites", "genename",
                            "transcript_id"))) %>%
    dplyr::relocate(MSA_site, TOGA_Human_site,
                    TOGA_Acomys_site, UP_Human_site,
                    branch_subs, branch_subs_ancestral,
                    uniprotswissprot) %>%
    dplyr::arrange(TOGA_Human_site) %>%
    # Correct location site position of branch_subs_ancestral column
    dplyr::mutate(branch_subs_ancestral = purrr::map2_chr(TOGA_Acomys_site,
                                                          branch_subs_ancestral,
                                                          translate_branch_subs)) %>%
    # remove non-significant sites if input$MEMETableSignSites is TRUE
    {if (input$MEMETableSignSites) dplyr::filter(., `p-value` <= input$p_val_select) else .}
})

# Get length of unimapped results
get_unimapped_MEME_length <- reactive({

  # Stop if no gene selected
  req(input$MEMEGeneInput)

  nrow(get_unimapped_MEME_results())
})

# Rendering table for uniprot-mapped sites
output$MEME_table_unimapped <- gt::render_gt({

  if(nrow(get_unimapped_MEME_results()) == 0) {
    return(NULL)
  }

  get_unimapped_MEME_results() %>%
    dplyr::rename("Amino acid difference between Human and A. cahirinus sequence" = "branch_subs",
           "Amino acid substitution that took place in the A. cahirinus lineage (i.e. branch)" = "branch_subs_ancestral") %>%
  gt() %>%
    tab_spanner(
      label = html("<strong><em>Position</em></strong>"),
      columns = c("MSA_site", "TOGA_Human_site", "TOGA_Acomys_site", "UP_Human_site")
    ) %>%
    tab_spanner(
      label = html("<strong><em>MEME maximum likelyhood estimation output</em></strong>"),
      columns = c(`&alpha;`, `&beta;<sup>-</sup>`, `p<sup>-</sup>`,
                  `&beta;<sup>+</sup>`, `p<sup>+</sup>`,
                  LRT, `p-value`,
                  `#_branches_under_selection`, `Total_branch_length`,
                  `MEME_LogL`, `FEL_LogL`, `Variation_p`)) %>%
    cols_align(align = "left") %>%
    opt_stylize(style = 1) %>%
    cols_label(
      MSA_site = "MSA site",
      TOGA_Human_site = "TOGA Human site",
      TOGA_Acomys_site = "TOGA Acomys site",
      UP_Human_site = "UP Human site",
      `&alpha;` = gt::html("&alpha;"),
      `&beta;<sup>-</sup>` = gt::html("&beta;<sup>-</sup>"),
      `p<sup>-</sup>` = gt::html("p<sup>-</sup>"),
      `&beta;<sup>+</sup>` = gt::html("&beta;<sup>+</sup>"),
      `p<sup>+</sup>` = gt::html("p<sup>+</sup>")
    ) %>%
    data_color(
      columns = LRT,
      palette = "Blues") %>%
    data_color(
      columns = `p-value`,
      palette = "Greens",
      reverse = TRUE) %>%
    opt_interactive(use_compact_mode = TRUE,
                    use_resizers = TRUE,
                    use_page_size_select = TRUE,
                    page_size_values = c(10, 25, 50, get_unimapped_MEME_length()))

  })