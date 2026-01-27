### Function to create whitespace
generate_br_tags <- function(n) {
  tagList(replicate(n, br(), simplify = FALSE))
}

### Create function for generating vertical lines in plotly graphs
vline <- function(x = 0, color = "#92351E") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, width = 1)
  )
}

### Function for generating rectangles in plotly graphs
ly_rect <- function(x0, x1, y0, y1, fillcolor) {
  list(type = "rect",
       fillcolor = fillcolor, opacity = 0.8,
       x0 = x0, x1 = x1, xref = "x",
       y0 = y0, y1 = y1, yref = "y",
       layer = "below",
       line = list(color = fillcolor))
}

### Function for creating the interactive aBSREL STRING cluster-based
### scatterplot.
###   args:
###     dNdS_acomys_clust: The dN/dS and STRING clusters combined dataframe
###     pal: The color mappings to the data-values for the traces
###     dNdS_column: Either 'mean_dNdS' or 'max_omega'. Thus, this selects if
###                  you want to visualize the mean dN/dS or the highest
###                  estimated omega class per gene.
###     xaxis_title: The xaxis title that is selected.
###     hover_text: Changes the hovertext for the dN/dS values
Create_STRINGScatterplot <- function(dNdS_acomys_clust,
                                     pal,
                                     dNdS_column = "mean_dNdS",
                                     xaxis_title = "Mean dN/dS",
                                     hover_text = "Mean dN/dS") {
  dNdS_column <- sym(dNdS_column)
  
  p <- eval(expr({
    # Create ranges for the hidden traces
    rangex <- range(`$`(dNdS_acomys_clust, !!dNdS_column))
    rangey <- range(dNdS_acomys_clust$P_value_c)
    
    # Plot plotly, where a separate trace is generated for each STRING cluster
    # (the non-clustered genes too have a separate trace).
    p <- plot_ly()
    for (clust in c(seq(max(dNdS_acomys_clust$Cluster, na.rm = TRUE)), NA)) {
      if (is.na(clust)) {
        p <- add_trace(p,
                       data = dplyr::filter(dNdS_acomys_clust, is.na(Cluster)),
                       type = "scatter",
                       mode = "markers",
                       x = ~ !!dNdS_column,
                       y = ~P_value_c,
                       color = ~sign,
                       colors = pal,
                       legendgroup = "Not Clustered",
                       text = ~paste0("<b>gene:</b> ", genename, "<br>",
                                      "<b>Ensembl transcript ID:</b> ", transcript_id, "<br>",
                                      "<b>FDR:</b> ", P_value_c, "<br>",
                                      "<b>-log10(FDR):</b> ", -log10(P_value_c), "<br>",
                                      "<b>", hover_text, ":</b> ", !!dNdS_column, "<br>",
                                      "<b>log10(", hover_text, "):</b> ", log10(!!dNdS_column)))
      } else {
        p <- add_trace(p,
                       data = dplyr::filter(dNdS_acomys_clust, Cluster == clust),
                       type = "scatter",
                       mode = "markers",
                       x = ~ !!dNdS_column,
                       y = ~P_value_c,
                       color = ~sign,
                       colors = pal,
                       legendgroup = paste("Cluster", clust),
                       text = ~paste0("<b>gene:</b> ", genename, "<br>",
                                      "<b>Ensembl transcript ID:</b> ", transcript_id, "<br>",
                                      "<b>FDR:</b> ", P_value_c, "<br>",
                                      "<b>-log10(FDR):</b> ", -log10(P_value_c), "<br>",
                                      "<b>mean dN/dS:</b> ", !!dNdS_column, "<br>",
                                      "<b>log10(mean dN/dS):</b> ", log10(!!dNdS_column)))
      }
    }
    # Add styling using layout()
    p %>% layout(xaxis = list(title = xaxis_title, type = "log", tickformat = ".1e"),
                 yaxis = list(title = "Acomys cahirinus diversifying selection P-value", type = "log", tickformat = ".1e", autorange = "reversed"),
                 legend = list(title = list(text = "<b>STRING clusters</b>"))) %>%
      # Add additional hidden traces with the same ranges as the x and y-axes,
      # so that the plot does not resize when a cluster is selected/deselected
      add_trace(x = rangex[2], y = rangey[2], type = "scatter", mode = "markers",
                showlegend = FALSE, opacity = 0, hoverinfo='skip') %>%
      add_trace(x = rangex[1], y = rangey[1], type = "scatter", mode = "markers",
                showlegend = FALSE, opacity = 0, hoverinfo='skip') 
  }))
  return(p)
}

translate_branch_subs <- function(aco_site, branch_sub) {
  
  if ((is.na(branch_sub)) | (is.na(aco_site))) {
    return(NA_character_)
  }
  
  subst_split <- strsplit(
    gsub("([a-zA-Z-]+)(\\d+)([a-zA-Z-]+)", "\\1,\\2,\\3", branch_sub),
    ","
  )[[1]]
  
  return(paste0(subst_split[1], aco_site, subst_split[3]))
}

# Function to link uniprot sequence annotations to the unimapped MEME
# substitutions data
link_protein_features <- function(site_df, features_df,
                                  exclude_features = c("CHAIN")) {
  
  uni_subs_features <- site_df %>%
    mutate(features = map2(site, genename, ~ {
      site <- as.integer(.x)
      gene <- .y
      
      # filter features df for selected gene
      domain_df <- features_df %>%
        dplyr::filter(Gene == gene) %>%
        dplyr::filter(!(feat_type %in% exclude_features)) %>%
        # determine if a site falls in a feature/features
        mutate(in_feature = site >= as.integer(Start) & site <= as.integer(End)) %>%
        dplyr::filter(in_feature == TRUE) %>%
        dplyr::select(feat_name, feat_type, Start, End)
      
      if (nrow(domain_df) == 0) {
        return(NA)
      } else {
        return(domain_df)
      }
    }))
  
  # unnest the features
  uni_subs_features <- uni_subs_features %>%
    unnest(features)
  
  if (!("feat_name" %in% colnames(uni_subs_features))) {
    # When no sites were linked to Uniprot features, still create the feat_name,
    # feat_type, Start and End columns
    uni_subs_features$feat_type <- NA
    uni_subs_features$feat_name <- NA
    uni_subs_features$Start <- NA
    uni_subs_features$End <- NA
  }
  
  return(uni_subs_features)
}

identify_feature_sites <- function(site_df, features_df,
                                   selected_features = c("DOMAIN", "REGION",
                                                         "MOTIF")) {
  # function that creates a column that identifies if a site was present in the
  # selected features
  
  uni_subs_identified <- site_df %>%
    mutate(in_feature = map2_lgl(site, genename, ~ {
      site <- as.integer(.x)
      gene <- .y
      
      # filter features df for selected gene
      domain_df <- features_df %>%
        dplyr::filter(entryName == gene) %>%
        dplyr::filter(type %in% selected_features) %>%
        # determine if a site falls in the selected features
        mutate(in_feature = site >= as.integer(begin) & site <= as.integer(end))
      
      if (any(domain_df$in_feature)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }))
  
  return(uni_subs_identified)
}

### Function that takes a vector of Uniprot features and simplifies the
### labels. For example, 'Fibronectin type-III 4' becomes
### 'Fibronectin type-III'. The output is a vector with the simplified
### labels 
Translate_UP_features <- function(Uniprot_features) {
  case_when(
    str_detect(Uniprot_features, coll("Cadherin")) ~ "Cadherin",
    str_detect(Uniprot_features, coll("DRBM")) ~ "DRBM",
    str_detect(Uniprot_features, regex("EGF-like \\d+; calcium-binding")) ~ "EGF-like; calcium-binding",
    str_detect(Uniprot_features, coll("UvrD-like helicase ATP-binding")) ~ "UvrD-like helicase ATP-binding",
    str_detect(Uniprot_features, coll("Fibronectin type-III")) ~ "Fibronectin type-III",
    str_detect(Uniprot_features, coll("Ig-like C2-type")) ~ "Ig-like C2-type",
    str_detect(Uniprot_features, coll("Bromo")) ~ "Bromo",
    str_detect(Uniprot_features, coll("BRCT")) ~ "BRCT",
    str_detect(Uniprot_features, coll("Sushi")) ~ "Sushi",
    str_detect(Uniprot_features, coll("PAS")) ~ "PAS",
    str_detect(Uniprot_features, coll("EF-hand")) ~ "EF-hand",
    str_detect(Uniprot_features, coll("RRM")) ~ "RRM",
    str_detect(Uniprot_features, coll("WW")) ~ "WW",
    str_detect(Uniprot_features, coll("EGF-like")) ~ "EGF-like",
    str_detect(Uniprot_features, coll("TSP type-1")) ~ "TSP type-1",
    str_detect(Uniprot_features, coll("CUB")) ~ "CUB",
    str_detect(Uniprot_features, regex("EGF-like \\d+; incomplete")) ~ "EGF-like; incomplete",
    str_detect(Uniprot_features, coll("Beta/gamma crystallin 'Greek key'")) ~ "Beta/gamma crystallin 'Greek key'",
    str_detect(Uniprot_features, coll("Ig-like")) ~ "Ig-like",
    str_detect(Uniprot_features, coll("Ig-like")) ~ "Ig-like",
    Uniprot_features == "NONE" ~ "NONE",
    is.character(Uniprot_features) ~ Uniprot_features
  )
}

cbind_genelists <- function(genelist_paths, col_names) {
  # Function that reads the genelists and binds the columns in a single dataframe.
  #   args:
  #     genelist_paths: A vector of the filepaths of the genelists.
  #     col_names: A vector of the names that the columns should have in the
  #                combined datraframe. Should be in the same order as genelist_paths.
  require(readr)
  require(dplyr)
  
  dfs <- lapply(seq_along(genelist_paths), function(i) {
    read_csv(genelist_paths[i], col_names = col_names[i], show_col_types = FALSE)
  })
  return(do.call("cbind", dfs))
}