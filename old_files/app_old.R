library(shiny)
library(tidyverse)
library(magrittr)
library(DT)
library(tibble)
library(shinythemes)
library(heatmaply)
library(plotly)
library(gt)
library(ggprism)
library(scales)
library(MetBrewer)
library(ggrepel)

### load aBSREL data
load(file = "data/aBSREL_data.RData")

### load name conversion file
name_conversion <- read_csv("data/name_conversion.csv")

### Load the STRING cluster data
acomys_clusters <- read_csv("data/acomys_bias_one2one--clustered--infl4.csv")

### Load the aBSREL dN/dS dataframe
dNdS_df <- read_delim("data/aBSREL_acomys_dNdS_df.tsv",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

dNdS_df <- dNdS_df %>%
  # calculate and retrieve max omega values
  mutate(max_omega = if_else(rate_class_number == 1, omega_1,
                             if_else(rate_class_number == 2, omega_2,
                                     if_else(rate_class_number == 3, omega_3, NA)))) %>%
  # calculate and retrieve max omega proportions
  mutate(max_omega_prop = if_else(rate_class_number == 1, omega_1_prop,
                                  if_else(rate_class_number == 2, omega_2_prop,
                                          if_else(rate_class_number == 3, omega_3_prop, NA)))) %>%
  # calculate and retrieve mean omega values
  mutate(mean_dNdS = ifelse(rate_class_number == 1, omega_1,
                            if_else(rate_class_number == 2, 
                                    omega_1*omega_1_prop+omega_2*omega_2_prop,
                                    if_else(rate_class_number == 3, 
                                            omega_1*omega_1_prop+omega_2*omega_2_prop+omega_3*omega_3_prop, NA))))

### Load the functional enrichments
clusterprofiler_allgenes_GO_KEGG <- read_delim("data/functional_enrichments/clusterprofiler_results_GO_KEGG.tsv",
                                               delim = "\t", escape_double = FALSE,
                                               trim_ws = TRUE)

### Load the cluster-specific functional enrichments
enrichGO_aco <- read_delim("data/functional_enrichments_cluster/enrichGO_aco_all.tsv",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Acomys_reduced_GO_clusters <- read_delim("data/functional_enrichments_cluster/Acomys_clusters_GO_df.tsv",
                                         delim = "\t", escape_double = FALSE,
                                         trim_ws = TRUE)

### Load the uniprot-mapped substitutions
MEME_subs <- read_delim("data/Uniprot_mapping/MEME_subs_mapped_sitefiltered_noMSAs.tsv",
                        delim = "\t", escape_double = FALSE,
                        col_types = cols(n_of_sites = col_integer(),
                                         sites = col_character()), trim_ws = TRUE)

### Load the uniprot-mapped substitutions for the MSA genelist
MSA_genelist <- read_delim("data/Uniprot_mapping/MEME_subs_uniprot_nogaps_onlysubsites.tsv",
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                        col_select = "genename")$genename

### Load the uniprot sequence annotations
progo_out <- read_delim("data/Uniprot_annotations/progo_out_combined.tsv",
                        delim = "\t", escape_double = FALSE,
                        trim_ws = TRUE)
# Remove unneccesary line
progo_out <- progo_out %>%
  dplyr::filter(!(feat_name == "Processed zona pellucida sperm-binding protein 1" & is.na(End)))

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

### load the genenames and corresponding ensemble IDs for the significant and
### non-significant genes

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

### get function arguments (filepaths and column names)
genelist_path <- file.path("data", "genelists_bias_one2one")
genelist_sign_path <- c(file.path(genelist_path, "aBSREL_human_gene_ids_sign.txt"),
                        file.path(genelist_path, "aBSREL_human_genenames_sign.txt"),
                        file.path(genelist_path, "aBSREL_human_transcript_ids_sign.txt"))
genelist_path <- file.path("data", "genelists_one2one")
genelist_all_path <- c(file.path(genelist_path, "aBSREL_human_gene_ids.txt"),
                       file.path(genelist_path, "aBSREL_human_genenames.txt"),
                       file.path(genelist_path, "aBSREL_human_transcript_ids.txt"))
col_names <- c("ensemble", "names", "ensembleTrans")

### get genelist dataframes
sign_genes <- cbind_genelists(genelist_sign_path, col_names)
all_genes <- cbind_genelists(genelist_all_path, col_names)

### get genenames and transcript IDs for aBSREL
# all tested genes
genenames <- all_genes$names
transcript_ids <- all_genes$ensembleTrans
# sign genes
genenames_sign <- sign_genes$names
transcript_ids_sign <- sign_genes$ensembleTrans

### get genenames and transcript IDs for MEME results
MEME_genenames <- sign_genes$names
MEME_transcript_ids <- sign_genes$ensembleTrans

# Define UI ----
ui <- navbarPage(
  title = HTML(paste(em("Acomys cahirinus"), "study")),
  theme = shinytheme("flatly"),
  tabPanel("aBSREL results", 
           fluidRow(
             tabsetPanel(
               tabPanel("aBSREL LRT results",
                        column(3,
                               wellPanel(
                                 h3("Gene selection"),
                                 selectizeInput(inputId = "aBSRELGeneInput", "Select the gene for which you want to see the results", 
                                                choices = NULL),
                                 checkboxInput("aBSRELSignSwitch", HTML(paste("Only allow selecting genes with significant signs of episodic",
                                                                              "diversifying selection in", em("Acomys cahirinus"))), value = FALSE),
                                 checkboxInput("aBSRELTranscriptIDSelect", "Choose genes based on their transcript IDs", value = FALSE)),
                               wellPanel(
                                 h3("Table options"),
                                 checkboxInput("aBSRELTableSignBranches", "Only show significant Branches", value = FALSE)
                               )),
                        column(9,
                               div(p(strong("This page displays the aBSREL likelyhood ratio test (LRT) results for the selected gene.",
                                        "Every row represents an aBSREL LRT performed for a specific phylogenetic branch (i.e. species).",
                                        "The", em("Homo sapiens, A. cahirinus"), "and", em("M. musculus"), "branches were tested for",
                                        "positive selection in every analyzed gene, while for some genes an additional random selection",
                                        "of 3 rodent branches and 9 non-rodent branches were tested for positive selection.",
                                        "Nevertheless, only the A. cahirinus branch LRT results were utilized for the generation of",
                                        "the", em("A. cahirinus"), "positively selected genelist.")),
                                   p("The following columns are listed in the table below:"),
                                   tags$ul(
                                     tags$li(strong(em("Branch name:")), em("The name of the tested phylogenetic branch, which is a combination of the genome assembly name and the TOGA projection ID")),
                                     tags$li(strong(em("Genome assembly name:")), em("The identifier of the genome assembly used for this species (i.e. the genome assembly were TOGA identified the orthologs)")),
                                     tags$li(strong(em("TOGA projection ID:")), em("The transcript name that was given to the projected ortholog found in this genome assembly. This unique transcript identifier is made up of the human reference transcript ID, the genename and the chain ID (in the format ReferenceTranscript.GeneSymbol.ChainID)")),
                                     tags$li(strong(em("Likelyhood ratio test statistic:")), em("The LRT statistic that was used to calculate the uncorrected and FDR-adjusted p-values")),
                                     tags$li(strong(em("Uncorrected P-value:")), em("The raw (episodic) positive selection P-value for this branch, calculated from the LRT statistic. See the", a("aBSREL paper by Smith et al. (2015)", href = "https://pubmed.ncbi.nlm.nih.gov/25697341/"), "for the exact details on this calculation.")),
                                     tags$li(strong(em("FDR-adjusted P-value:")), em("The (episodic) positive selection P-value for this branch, adjusted for the number of genes analyzed by HyPhy aBSREL")),
                                     tags$li(strong(em("Uncorrected P-value significant (P-value <= 0.05):")), em("A logical value (i.e. true or false) indicating if the uncorrected P-value was below the significance threshold of 0.05")),
                                     tags$li(strong(em("FDR-adjusted P-value significant (P-value <= 0.05):")), em("A logical value (i.e. true or false) indicating if the FDR-adjusted P-value was below the significance threshold of 0.05")),
                                  )),
                               br(),
                               uiOutput("aBSRELTableDescription"),
                               br(),
                               DT::dataTableOutput("aBSREL_table"),
                        )),
               tabPanel("STRING interaction network",
                        column(3,
                               wellPanel(
                                 h4("Display options for STRING clusters"),
                                 checkboxInput("STRINGplotSwitch", HTML("Display all clusters with at least three genes"), value = FALSE),
                                 br(),
                                 sliderInput("STRINGClustTableSlider", "Cluster number:",
                                             min = 1, max = max(enrichGO_aco$`\`__mclCluster\``,
                                                                na.rm = TRUE),
                                             value = c(1, 8)),
                                 checkboxInput("STRINGDisplayNonClust", HTML("Display non-clustered genes in table"), value = FALSE),
                                 actionButton("STRINGClustTableAction", "Refresh")
                               )),
                        column(9,
                               div(p(strong("This page displays the STRING clusters that were identified in the", em("A. cahirinus"),
                                            "(episodic) positively selected genelist. The clusters were generated using the Markov cluster (MCL) algorithm with an inflation parameter of 4")),
                                   p(strong("We identified 8 STRING clusters with at least 5 genes and 64 STRING clusters with the minimum of 2 genes",
                                            "These can be explored using the interactive scatterplot and table below. Furthermore, the desired STRING",
                                            "clusters in the table can be selected using the range buttons on the left side of the page."))
                                   
                               ),
                               uiOutput("aBSRELSTRINGplot"),
                               br(),
                               p(strong("Double click a cluster in the cluster legend of the scatterplot to view that cluster and hide all others.",
                                        "Double click that cluster again in the cluster legend to visualize all clusters again.",
                                        "To view specific details about specific genes, hover over them in the scatterplot")),
                               fluidRow(
                                 column(12,
                                        checkboxInput("STRINGScatterdNdSSwitch", "Display scatterplot for the highest estimated dN/dS values (i.e. the highest estimated omega class) per gene", value = FALSE,
                                                      width = "100%"))),
                               plotlyOutput("aBSRELSTRINGScatterplot"),
                               br(),
                               div(p(strong("This table displays the HyPhy aBSREL output data for the", em("A. cahirinus"),
                                            "Positively selected genelist for each individual STRING cluster. Select the",
                                            "desired STRING clusters using the range buttons on the left side of the page")),
                                   p("The following columns are listed in the table below:"),
                                   tags$ul(
                                     tags$li(strong(em("STRING cluster:")), em("The STRING cluster in which the gene was clustered by the MCL algorithm.")),
                                     tags$li(strong(em("genename:")), em("The gene symbol")),
                                     tags$li(strong(em("Likelyhood ratio test statistic:")), em("The LRT statistic that was used to calculate the uncorrected and FDR-adjusted p-values")),
                                     tags$li(strong(em("Uncorrected P-value:")), em("The raw (episodic) positive selection P-value for this branch, calculated from the LRT statistic. See the", a("aBSREL paper by Smith et al. (2015)", href = "https://pubmed.ncbi.nlm.nih.gov/25697341/"), "for the exact details on this calculation.")),
                                     tags$li(strong(em("FDR-adjusted P-value:")), em("The (episodic) positive selection P-value for this branch, adjusted for the number of genes analyzed by HyPhy aBSREL")),
                                     tags$li(strong(em("Inferred branch length:")), em("The length of the A. cahirinus branch, as estimated in the aBSREL full adaptive model. This length is given as the number of substitutions per site.")),
                                     tags$li(strong(em("Inferred dN branch length:")), em("The non-synonynous component of the inferred A. cahirinus branch length. This length is given as the number of non-synonymous substitutions per site.")),
                                     tags$li(strong(em("Inferred dS branch length:")), em("The synonynous component of the inferred A. cahirinus branch length. This length is given as the number of synonymous substitutions per site.")),
                                     tags$li(strong(em("Mean dN/dS:")), em("The mean dN/dS ratio (i.e. the non-synonymous substitution rate divided by the synonymous substitution rate) estimated for the A. cahirinus branch for this gene. This mean dN/dS ratio is calculated by averaging the dN/dS ratios for each estimated omega class, taking into account the proportion of sites assigned to these omega classes.")),
                                     tags$li(strong(em("Maximum dN/dS:")), em("The maximum (max) dN/dS is the highest estimated dN/dS ratio (i.e. omega class) for the A. cahirinus branch for this gene")),
                                     tags$li(strong(em("proportion of sites with maximum dN/dS:")), em("The proportion of codon sites were the maximum dN/dS was estimated. In other words, the proportion of codon sites which were assigned to the omega class with the highest dN/dS ratio.")),
                                   )),
                               gt_output("STRINGClusterTable"))),
               tabPanel("Functional overrepresentation analysis",
                        column(3,
                               wellPanel(
                                 h4("Display options for functionel enrichment of full aBSREL genelist"),
                               ),
                               generate_br_tags(40),
                               wellPanel(
                                 h4("Display options for functionel enrichment of STRING clusters"),
                                 checkboxInput("FunEnrichReducedSwitch", HTML("Display all GO-terms per cluster instead of only the reduced GO-terms"), value = FALSE),
                                 h4("Select which clusters to display"),
                                 sliderInput("FunEnrichSelectClust", "Cluster number:",
                                             min = 1, max = max(enrichGO_aco$`\`__mclCluster\``,
                                                                na.rm = TRUE),
                                             value = c(1, 8)),
                                 actionButton("FunEnrichClustAction", "Refresh")
                               )),
                        column(9,
                               div(p(strong("This page contains the functional overrepresentation analysis results of both complete aBSREL positively selected genelist and the cluster-specific overrepresentation results",
                                  "The clusterProfiler (4.10.0) R-package was used to perform this functional over-representation analysis. Moreover, we used the rrvgo R package (1.15.1) to reduce the number of",
                                  "cluster-specific significant GO-terms based on their semantic similarity (Sayols, 2023)")),
                                  p("Just like the STRING interaction network table, the desired STRING clusters in the table",
                                    "and plot can be selected using the range buttons on the left side of the page.",
                                    "There is also the option to display all the GO-terms per cluster, and not just the reduced GO-terms.")),
                               hr(),
                               strong("Functional overrepresentation analysis of complete positively selected genelist"),
                               br(),
                               plotOutput("FunEnrichaBSRELAll", width = "100%", height = "300px"),
                               br(),
                               gt_output("FunEnrichaBSRELAllTable"),
                               hr(),
                               strong("Functional overrepresentation analysis of the individual positively selected STRING clusters"),
                               br(),
                               # The inline = TRUE is required so that the plot
                               # and table don't overlap
                               plotOutput("FunEnrichaBSRELClusters", inline = TRUE),
                               br(),
                               gt_output("FunEnrichaBSRELClustersTable")
                               ))
             )
           )
  ),
  tabPanel("MEME results",
                    fluidRow(
                      column(3,
                             wellPanel(
                               h3("Gene selection"),
                               selectizeInput(inputId = "MEMEGeneInput", "Select the gene for which you want to see the results", 
                                              choices = NULL),
                               checkboxInput("MEMETranscriptIDSelect", "Choose genes based on their transcript IDs", value = FALSE)),
                             wellPanel(
                               h3("Table options"),
                               checkboxInput("MEMETableSignSites", "Only show significant sites", value = FALSE),
                               numericInput("p_val_select", "p-value significance cut-off", value = 0.1)
                             )),
                      column(9,
                             tabsetPanel(
                               tabPanel("MEME LRT results",
                                        br(),
                                        uiOutput("MEMETableDescription"),
                                        br(),
                                        p(paste("This HyPhy MEME analysis was a follow-up analysis after",
                                                "the HyPhy aBSREL analysis. During the aBSREL analysis,",
                                                "genes were found that exhibited signs of episodic diversifying",
                                                "selection in Acomys cahirinus. However, HyPhy aBSREL",
                                                "detects episodic diversifying selection in specific branches of",
                                                "a phylogeny without informing us which sites are responsible.",
                                                "To get an indication of which codon sites could be under episodic",
                                                "diversifying selection in the genes found by aBSREL, HyPhy MEME",
                                                "was used. HyPhy MEME is a mixed effects model of evolution (MEME)",
                                                "that can detect episodic diversifying selection at an individual",
                                                "site (among all branches). Thus, HyPhy MEME might give us additional",
                                                "insights into which codon sites might be responsible for the episodic",
                                                "diversifying selection signal in the genelist found by aBSREL.",
                                                "A word of caution is advised however, as HyPhy MEME can only detect",
                                                "if a codon site shows signs of diversifying selection among all branches",
                                                "(i.e. MEME detects selection at an individual site, not an individual",
                                                "branch-site). So, HyPhy MEME predicting that a codon site shows",
                                                "significant diversifying selection does not guarantee that any individual",
                                                "branch (like the Acomys cahirinus branch) exhibits diversifying selection",
                                                "at that codon site. nevertheless, HyPhy MEME can give us an indication", 
                                                "which sites might be involved in the episodic diversifying selection signal",
                                                "in the Acomys cahirinus genes, even when this limitation is kept in mind.")),
                                        p(paste("Every Gene was individually analyzed using HyPhy MEME. The test results are",
                                                "described in the tables below, where the column names represent the following:")),
                                        tags$ul(
                                          tags$li(HTML(paste(strong(HTML("&alpha;")), ": Synonymous substitution rate at a site."))),
                                          tags$li(HTML(paste(strong(HTML("&beta;<sup>-</sup>")), ": Non-synonymous substitution rate at a site for the negative/neutral evolution component."))),
                                          tags$li(HTML(paste(strong(HTML("p<sup>-</sup>")), ": Mixture distribution weight allocated to &beta;<sup>-</sup>; loosely -- the proportion of the tree evolving neutrally or under negative selection."))),
                                          tags$li(HTML(paste(strong(HTML("&beta;<sup>+</sup>")), ": Non-synonymous substitution rate at a site for the positive/neutral evolution component."))),
                                          tags$li(HTML(paste(strong(HTML("p<sup>+</sup>")), ": Mixture distribution weight allocated to &beta;<sup>+</sup>; loosely -- the proportion of the tree evolving neutrally or under positive selection."))),
                                          tags$li(HTML(paste(strong("LRT"), ": Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;."))),
                                          tags$li(HTML(paste(strong("p-value"), ": Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;."))),
                                          tags$li(HTML(paste(strong("number of branches under selection"), ": The (very approximate and rough) estimate of how many branches may have been under selection at this site, i.e., had an empirical Bayes factor of 100 or more for the &beta;<sup>+</sup> rate."))),
                                          tags$li(HTML(paste(strong("Total branch length"), ": The total length of branches contributing to inference at this site, and used to scale dN-dS."))),
                                          tags$li(HTML(paste(strong("MEME LogL"), ": Site Log-likelihood under the MEME model."))),
                                          tags$li(HTML(paste(strong("FEL LogL"), ": Site Log-likelihood under the FEL model."))),
                                          tags$li(HTML(paste(strong("Variation p"), ": Asymptotic p-value for whether or not there is evidence of dN/dS variation across branches."))),
                                        ),
                                        hr(),
                                        p(strong(HTML(paste("This table and summary data represents the MEME likelyhood ratio test results for each site of the",
                                                "selected gene. Here, the positional site information is relative to the human gene",
                                                "that was used in the TOGA multiple sequence alignment (i.e. the multiple sequence",
                                                "alignment with mammalian orthologs of the selected gene, available at the",
                                                a("TOGA ortholog database", href="https://genome.senckenberg.de//download/TOGA/human_hg38_reference/MultipleCodonAlignments/")
                                                )))),
                                        br(),
                                        DT::dataTableOutput("MEME_table"),
                                        br(),
                                        uiOutput("MEMEResultDescription"),
                                        hr(),
                                        p(strong(HTML(paste("This table and summary data represents the positively-selected sites that could be mapped to a canonical Uniprot protein sequence.",
                                                            "Here, the positional site information is relative to the canonical protein sequence of the human Uniprot accession.",
                                                            "Furthermore, only the sites are displayed were an amino acid substitution took place in the A. cahirinus branch of the phylogeny."
                                        )))),
                                        br(),
                                        DT::dataTableOutput("MEME_table_unimapped"),
                                        br(),
                                        uiOutput("MEMEResultDescription_unimapped"),
                                        hr()
                                        ),
                               tabPanel("Multiple sequence alignment",
                                        h3("Multiple sequence alignment"),
                                        p(paste("This tab displays the protein multiple sequence alignment (MSA),",
                                                "which were created by translating the TOGA codon MSAs.")),
                                        p("Description of the MSAs:"),
                                        tags$ul(
                                          tags$li("Blue columns represent columns where over 50% of the species in the MSA had an identical amino acid."),
                                          tags$li("A dot represents a residue that was an 'NNN' codon (i.e. A codon that was masked in the original TOGA alignment) or a codon gap. Keep in mind that this MSA represents the MEME input sequence, where positions that represented gaps in the human sequence were removed. Thus, the sequences of the other mammals in this MSA can be missing insertions and could be shorter than the original TOGA ortholog."),
                                          tags$li("The red arrows point at columns that have evidence of episodic diversifying selection, according to HyPhy MEME.")),
                                        br(),
                                        radioButtons(
                                          "MSAradioselect",
                                          "Select MSA display option",
                                          choices = list("Identity/conservation" = 1, "Charge" = 2, "Structure" = 3),
                                          selected = 1
                                        ),
                                        uiOutput("MultipleAlignmentText"),
                                        uiOutput("MultipleAlignment"),
                                        hr()
                                        ),
                               tabPanel("Substitutions",
                                        div(h3("Predicted substitutions from the last common mammalian ancestor to", em("A. cahirinus"), "at all diversifying codon sites"),
                                            p("Of note: These site positions are relative to the human gene that was used in the TOGA MSA (i.e. the multiple sequence",
                                              "alignment with mammalian orthologs of the selected gene, available at the",
                                              a("TOGA ortholog database.", href="https://genome.senckenberg.de//download/TOGA/human_hg38_reference/MultipleCodonAlignments/"),
                                              "Thus, these site positions can not be directly used for locating substitutions on the Uniprot canonical sequence.")),
                                        hr(),
                                        DT::dataTableOutput("MEME_acomys_sub"),
                                        hr(),
                                        h3("All substitutions at diversifying codon sites"),
                                        div(p("This table displays the predicted substitutions (by HyPhy MEME) at diversifying sites of all tested species. Substitutions that appear to be indels are indicated with the term 'Gap' (e.g. C-333-Gap)."),
                                            p("Of note: These site positions are relative to the human gene that was used in the TOGA MSA (i.e. the multiple sequence",
                                              "alignment with mammalian orthologs of the selected gene, available at the",
                                              a("TOGA ortholog database.", href="https://genome.senckenberg.de//download/TOGA/human_hg38_reference/MultipleCodonAlignments/"),
                                              "Thus, these site positions can not be directly used for locating substitutions on the Uniprot canonical sequence.")),
                                        hr(),
                                        DT::dataTableOutput("MEME_branch_subs"),
                                        hr()
                               ),
                               tabPanel("EBF",
                                        uiOutput("EBFBubbleplotDesc"),
                                        plotlyOutput("EBF_bubbleplot"),
                                        br(),
                                        uiOutput("EBFTableDesc"),
                                        fluidRow(
                                          column(3,
                                                 checkboxInput("MEME_EBF_inf_remove", "Remove infinite EBF values", value = TRUE)),
                                          column(9,
                                                 checkboxInput("MEME_EBF_acomys_select", "Only select Acomys cahirinus EBF values", value = TRUE))),
                                        gt_output("EBF_table")
                                        ),
                               tabPanel("Uniprot sequence annotations",
                                        uiOutput("DomainPageDesc"),
                                        br(),
                                        p(HTML(paste0("For the genes where positively selected sites were mapped to ",
                                                 "their corresponding Uniprot accessions, this page displays ",
                                                 "these sites mapped to the canonical sequence isoform. ",
                                                 "These protein sequence annotations for the genes were retrieved by ",
                                                 "the EMBL-EBI Proteins REST API, more specifically the features service. ",
                                                 "The protein visualized below is annotated with the DOMAIN, ",
                                                 "REGION and MOTIF Uniprot annotation types, which were retrieved . Other types of ",
                                                 "Uniprot annotations can be found in the table below."))),
                                        br(),
                                        strong(paste("Select A. cahirinus empirical Bayes factor (EBF)",
                                                     "cutoff value. EBF values below this value will not",
                                                     "be displayed in the Uniprot sequence annotation plot",
                                                     "and table")),
                                        numericInput("UniDomainEBFSelect", NULL, value = 100),
                                        br(),
                                        plotlyOutput("UniDomainPlot", height = "150px"),
                                        br(),
                                        uiOutput("UniDomainSiteDesc"),
                                        gt_output("UniDomainSiteTable"),
                                        hr(),
                                        uiOutput("UniDomainFeatDesc"),
                                        gt_output("UniDomainFeatTable"),
                                        )
                             )
                      )
                    )
           )
  )

# Define server logic ----
server <- function(input, output, session) {
  
  observe({

    # freezes all reactive expressions until everything has been updated.
    # Without this, there is one second where it tries to treat gene symbols as
    # transcript IDs or vice versa, leading to errors.
    freezeReactiveValue(input, "aBSRELGeneInput")
    
    if (input$aBSRELTranscriptIDSelect) {
      if (input$aBSRELSignSwitch) {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = transcript_ids_sign, server = TRUE,
                             options = list(placeholder = "Transcript id"))
      } else {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = transcript_ids, server = TRUE,
                             options = list(placeholder = "Transcript id"))
      }
    } else {
      if (input$aBSRELSignSwitch) {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = genenames_sign, server = TRUE,
                             options = list(placeholder = "Gene name"))
      } else {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = genenames, server = TRUE,
                             options = list(placeholder = "Gene name"))
      }
    }

  })
  
  # Selection of STRINGp cluster plot
  output$aBSRELSTRINGplot <- renderUI({
    if (input$STRINGplotSwitch) {
      tags$img(src = "./images/stringdb_all_clusters_with_three_genes.svg", width = 510 * 1.6, height = 472 * 1.6)
    } else {
      tags$img(src = "./images/Acomys_stringdb_clusters.svg", width = 510 * 1.5, height = 472 * 1.5)
    }
  })
  
  # Generate table of aBSREL info for each STRING cluster
  STRINGClustTable <- eventReactive(input$STRINGClustTableAction, {
    
    # Filter for one2one significant A. cahirinus genes
    aBSREL_table <- aBSREL_data %>%
      dplyr::filter(genename %in% sign_genes$names) %>%
      # filter for the presence of orthologs from nearest phylogenetic neighbors
      # of A. cahirinus
      dplyr::filter(MRCA_present_aco == TRUE) %>%
      dplyr::select(-all_of(c("Batch", "!_present", "JSON_exists",
                              "MRCA_present_aco", "MRCA_present_mus",
                              "significant_prop", "significant_prop",
                              "branches_tested")))
    
    # Add ensemble gene IDs to dataframe
    aBSREL_table <- all_genes %>%
      as_tibble() %>%
      select(-names) %>%
      rename(gene_id = ensemble, transcript_id = ensembleTrans) %>%
      right_join(aBSREL_table)
    
    # Join the omega (i.e. dN/dS) data with the aBSREL LRT data
    aBSREL_table_omega <- dNdS_df %>%
      dplyr::filter(species == "HLacoCah2") %>%
      dplyr::select(-all_of(c("transcript_id", "P_value", "MRCA_present_aco", "MRCA_present_mus",
                              "species"))) %>%
      right_join(aBSREL_table)
    
    # Retrieve the likelyhood ratio test statistic from the 'test_results'
    # list-column
    aBSREL_table_omega_test <- aBSREL_table_omega %>%
      rowwise() %>%
      mutate(test_results = dplyr::filter(test_results, branch_assembly == "HLacoCah2")) %>%
      unnest(test_results) %>%
      ungroup() %>%
      dplyr::select(-all_of(c("branch_name", "branch_projection", "p_value_cor",
                              "significant", "significant_cor"))) %>%
      # give Columns more representative names
      dplyr::rename(p_value_uncorrected = p_value, FDR = P_value)
    
    # add species name to data
    aBSREL_table_omega_test <- name_conversion %>%
      dplyr::select(`Assembly name`, Species_Tree) %>%
      rename(branch_assembly = `Assembly name`) %>%
      right_join(aBSREL_table_omega_test) %>%
      dplyr::select(-branch_assembly)
    
    # Add STRINGdb clusters to data
    aBSREL_table_omega_test_clusters <- acomys_clusters %>%
      rename(Cluster = `__mclCluster`, gene_id = `query term`) %>%
      dplyr::select(gene_id, Cluster) %>%
      right_join(aBSREL_table_omega_test)
    
    # Filter for selected cluster range
    if (input$STRINGDisplayNonClust) {
      aBSREL_table_omega_test_clusters <- aBSREL_table_omega_test_clusters %>%
        dplyr::filter(Cluster %in% seq(input$STRINGClustTableSlider[1],
                                       input$STRINGClustTableSlider[2]) |
                        is.na(Cluster))
    } else {
      aBSREL_table_omega_test_clusters <- aBSREL_table_omega_test_clusters %>%
        dplyr::filter(Cluster %in% seq(input$STRINGClustTableSlider[1],
                                       input$STRINGClustTableSlider[2]))
    }
    
    # Generate table for likelyhood ratio test statistics and model estimations
    # per gene and cluster
    gt_table <- aBSREL_table_omega_test_clusters %>%
      # Make custom cluster labels
      dplyr::mutate(Cluster_label = if_else(is.na(Cluster), "Not clustered",
                                            paste0("Cluster ", Cluster)),
                    Cluster_label = factor(Cluster_label,
                                           levels = c(paste0("Cluster ",
                                                             seq(max(Cluster, na.rm = TRUE))),
                                                      "Not clustered"))) %>%
      dplyr::select(all_of(c("Cluster_label", "genename", "LRT",
                             "p_value_uncorrected", "FDR", "full_adapt",
                             "full_adapt_dN", "full_adapt_dS",
                             "mean_dNdS", "max_omega", "max_omega_prop"))) %>%
      # sort by FDR
      dplyr::arrange(Cluster_label, FDR) %>%
      # Start creating the table
      gt(groupname_col = "Cluster_label", row_group_as_column = TRUE) %>%
      tab_header(title = "HyPhy aBSREL output per STRING cluster") %>%
      tab_spanner(label = html("<em><strong>Likelyhood ratio test</em></strong>"),
                  columns = c(LRT, p_value_uncorrected, FDR)) %>%
      tab_spanner(label = html("<strong><em>Model estimation</em></strong>"),
                  columns = c(full_adapt, full_adapt_dN,
                              full_adapt_dS, mean_dNdS, max_omega,
                              max_omega_prop)) %>%
      fmt_number(decimals = 2, sep_mark = "") %>%
      cols_label(contains("p_value_uncorrected") ~ "Uncorrected P-value",
                 contains("LRT") ~ "likelyhood ratio test statistic",
                 contains("FDR") ~ "FDR-adjusted P-value",
                 contains("full_adapt") ~ "Inferred branch length",
                 contains("full_adapt_dN") ~ "Inferred dN branch length",
                 contains("full_adapt_dS") ~ "Inferred dS branch length",
                 contains("mean_dNdS") ~ "Mean dN/dS",
                 contains("max_omega") ~ "maximum dN/dS",
                 contains("max_omega_prop") ~ "proportion of sites with maximum dN/dS") %>%
      tab_stubhead(label = "STRING cluster") %>%
      cols_align(align = "left") %>%
      opt_stylize(style = 1) %>%
      data_color(
        columns = FDR,
        palette = "Greens",
        reverse = TRUE) %>%
      data_color(
        columns = mean_dNdS,
        palette = "Blues"
      ) %>%
      tab_style(
        locations = cells_body(columns = genename),
        style = cell_fill(color = "gray95"))
    
      # Reorder row groups
    if (input$STRINGDisplayNonClust) {
      gt_table %>%
        row_group_order(groups = c(paste0("Cluster ", seq(input$STRINGClustTableSlider[1],
                                                          input$STRINGClustTableSlider[2])),
                                   "Not clustered")) %>%
        opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
    } else {
      gt_table %>%
        row_group_order(groups = paste0("Cluster ", seq(input$STRINGClustTableSlider[1],
                                                        input$STRINGClustTableSlider[2]))) %>%
        opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
    }
  }, ignoreNULL = FALSE)
  
  output$STRINGClusterTable <- render_gt({
    STRINGClustTable()
    })
  
  output$aBSRELSTRINGScatterplot <- renderPlotly({
    # remove genes from the data where the closest neighbors of the A. cahirinus
    # branch had no annotated ortholog (that was classified as intact, partial
    # intact or uncertain loss)
    dNdS_df <- dNdS_df %>%
      filter(MRCA_present_aco == TRUE) %>%
      # Add a separate column where P-values are slightly shifted
      # (by half of the minimum non-0 p-value in the data).
      mutate(P_value_c = P_value + (min(P_value[P_value > 0])/2))
    
    # For the upcoming log10 transformations, deal with the zero values. This is
    # done by adding the smallest value that is not zero to every individual value
    dNdS_df <- dNdS_df %>%
      mutate(mean_dNdS = mean_dNdS + min(mean_dNdS[mean_dNdS > 0], na.rm = TRUE),
             max_omega = max_omega + min(max_omega[max_omega > 0], na.rm = TRUE),
             full_adapt = full_adapt + min(full_adapt[full_adapt > 0], na.rm = TRUE),
             full_adapt_dN = full_adapt_dN + min(full_adapt_dN[full_adapt_dN > 0], na.rm = TRUE),
             full_adapt_dS = full_adapt_dS + min(full_adapt_dS[full_adapt_dS > 0], na.rm = TRUE))
    
    dNdS_acomys <- dNdS_df %>%
      # Filter for Acomys dN/dS data
      dplyr::filter(species == "HLacoCah2")
    
    dNdS_acomys_clust <-  acomys_clusters %>%
      # Add STRING clusters to dN/dS data
      dplyr::rename(Cluster = `__mclCluster`, genename = `display name`) %>%
      dplyr::select(genename, Cluster) %>%
      right_join(dNdS_acomys) %>%
      # Rename data values for plotly
      mutate(sign = if_else(sign == TRUE & is.na(Cluster),
                            "Not clustered - Significant",
                            if_else(sign == FALSE & is.na(Cluster),
                                    "Not clustered - Not significant",
                                    if_else(!is.na(Cluster),
                                            paste0("Cluster ", Cluster, " - Significant"),
                                            paste0("Cluster ", Cluster, " - Not significant")))))
    
    # Create color mapping for data values
    pal <- c(rep("#d37750", max(dNdS_acomys_clust$Cluster, na.rm = TRUE)),
             "#d37750", "#f2c88f")
    pal <- setNames(pal, c(paste0("Cluster ", seq(max(dNdS_acomys_clust$Cluster, na.rm = TRUE)),
                                  " - Significant"), "Not clustered - Significant", "Not clustered - Not significant"))
    
    if (input$STRINGScatterdNdSSwitch) {
      Create_STRINGScatterplot(dNdS_acomys_clust,
                               pal,
                               dNdS_column = "max_omega",
                               xaxis_title = "dN/dS of highest estimated omega class",
                               hover_text = "dN/dS highest omega class")
    } else {
      Create_STRINGScatterplot(dNdS_acomys_clust,
                               pal)
    }
  })
  
  output$FunEnrichaBSRELAll <- renderPlot({
    
    ggplot(clusterprofiler_allgenes_GO_KEGG,
           aes(y = reorder(Description, p.adjust, decreasing = FALSE), 
                      x = p.adjust, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      geom_linerange(aes(xmin = 0, xmax = p.adjust), linewidth = 2) +
      ylab("GO Biological process") + xlab("False discovery rate") +
      theme_prism() +
      scale_size_continuous(name = "number of\ngenes",
                            range = c(3, 10)) +
      scale_y_discrete(labels = label_wrap(30)) +
      theme(axis.title.y = element_blank(), legend.title = element_text()) +
      scale_color_gradientn(colors = MetBrewer::met.brewer("OKeeffe2", type = "continuous"),
                            name = "number of\ngenes") +
      facet_grid(scales = "free_y", rows = vars(source),
                 space = "free_y")
  })
  
  output$FunEnrichaBSRELAllTable <- render_gt({
    gt(clusterprofiler_allgenes_GO_KEGG, groupname_col = "source",
       row_group_as_column = TRUE) %>%
      cols_align(align = "left") %>%
      opt_stylize(style = 1) %>%
      data_color(
        columns = Count,
        palette = "Blues") %>%
      data_color(
        columns = qvalue,
        palette = "Greens",
        reverse = TRUE) %>%
      opt_interactive(use_resizers = TRUE)
      
  })

  get_GOterm_plot_height <- eventReactive(input$FunEnrichClustAction, {
    # Get the number of GO-terms that will be plotted, which will
    # specify the height of the plot.
    if (input$FunEnrichReducedSwitch) {
      n_of_terms <- enrichGO_aco %>%
        dplyr::filter(Cluster %in% seq(input$FunEnrichSelectClust[1],
                                       input$FunEnrichSelectClust[2])) %>%
        nrow()
      20 * n_of_terms
    } else {
      # Get the number of reduced GO-terms that will be plotted, which will
      # specify the height of the plot.
      n_of_terms <- Acomys_reduced_GO_clusters %>%
        dplyr::filter(stringdb_cluster %in% seq(input$FunEnrichSelectClust[1],
                                                input$FunEnrichSelectClust[2])) %>%
        nrow()
      15 * n_of_terms
    }
  }, ignoreNULL = FALSE)
  
  FunEnrichClustPlot <- eventReactive(input$FunEnrichClustAction, {
    
    if (input$FunEnrichReducedSwitch) {
      ggplot(subset(enrichGO_aco, Cluster %in% seq(input$FunEnrichSelectClust[1],
                                                   input$FunEnrichSelectClust[2])),
             aes(y = reorder(Description, -log10(`p.adjust`), decreasing = FALSE), 
                 x = -log10(`p.adjust`), color = Count)) +
        geom_point(aes(size = Count)) +
        geom_linerange(aes(xmin = 0, xmax = -log10(`p.adjust`))) +
        ylab("GO-terms (Biological process)") + xlab("-log10(FDR)") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 12)) +
        scale_size_continuous(name = "number of\ngenes") +
        scale_y_discrete(labels = label_wrap(100)) +
        # reverse color palette using rev()
        scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Cross", type = "continuous")),
                              name = "number of\ngenes") +
        facet_grid(scales = "free_y", rows = vars(Cluster),
                   space = "free_y")
    } else {
      ggplot(subset(Acomys_reduced_GO_clusters,
                    stringdb_cluster %in% seq(input$FunEnrichSelectClust[1],
                                              input$FunEnrichSelectClust[2])),
             aes(y = reorder(term, score, decreasing = FALSE), 
                 x = score, color = n_of_genes)) +
        geom_point(aes(size = terms_per_parent)) +
        geom_linerange(aes(xmin = 0, xmax = score)) +
        ylab("Reduced GO-terms (Biological process)") + xlab("-log10(FDR)") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 12)) +
        scale_size_continuous(name = "number of\nGO terms") +
        scale_y_discrete(labels = label_wrap(100)) +
        # reverse color palette using rev()
        scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Cross", type = "continuous")),
                              name = "number of\ngenes") +
        facet_grid(scales = "free_y", rows = vars(stringdb_cluster),
                   space = "free_y")
    }
  }, ignoreNULL = FALSE)
    
    output$FunEnrichaBSRELClusters <- renderPlot({
      FunEnrichClustPlot()
    },
    width = 1400,
    height = function() {get_GOterm_plot_height()})
    
    FunEnrichClustTable <- eventReactive(input$FunEnrichClustAction, {
      
      # Render reduced GO-terms unless FunEnrichReducedSwitch is set to TRUE
      if (input$FunEnrichReducedSwitch) {
        enrichGO_aco %>%
          dplyr::select(-`\`__mclCluster\``) %>%
          dplyr::filter(Cluster %in% seq(input$FunEnrichSelectClust[1],
                                         input$FunEnrichSelectClust[2])) %>%
          dplyr::mutate(Cluster_label = if_else(is.na(Cluster), "Not clustered",
                                                paste0("Cluster ", Cluster)),
                        Cluster_label = factor(Cluster_label,
                                               levels = c(paste0("Cluster ",
                                                                 seq(max(Cluster, na.rm = TRUE))),
                                                          "Not clustered"))) %>%
          dplyr::select(-Cluster) %>%
          gt(groupname_col = "Cluster_label", row_group_as_column = TRUE) %>%
          opt_stylize(style = 1) %>%
          cols_align(align = "left") %>%
          data_color(
            columns = qvalue,
            palette = "Greens",
            reverse = TRUE) %>%
          data_color(
            columns = Count,
            palette = "Blues") %>%
          opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
      } else {
        Acomys_reduced_GO_clusters %>%
          dplyr::select(-all_of(c("cluster", "parent", "parentTerm", "termDispensability",
                                  "termUniquenessWithinCluster"))) %>%
          dplyr::filter(stringdb_cluster %in% seq(input$FunEnrichSelectClust[1],
                                                  input$FunEnrichSelectClust[2])) %>%
          dplyr::mutate(Cluster_label = if_else(is.na(stringdb_cluster), "Not clustered",
                                                paste0("Cluster ", stringdb_cluster)),
                        Cluster_label = factor(Cluster_label,
                                               levels = c(paste0("Cluster ",
                                                                 seq(max(stringdb_cluster, na.rm = TRUE))),
                                                          "Not clustered"))) %>%
          dplyr::select(-stringdb_cluster) %>%
          dplyr::rename("Representative GO of reduced GO-terms" = go,
                        "-log10(FDR) of representative GO" = score,
                        "representative GO-term size" = size,
                        "GO description" = term,
                        "Term Uniqueness" = termUniqueness,
                        "number of genes in reduced GO-terms" = n_of_genes,
                        "Gene IDs" = geneID,
                        "Number of reduced GO-terms" = terms_per_parent) %>%
          gt(groupname_col = "Cluster_label", row_group_as_column = TRUE) %>%
          opt_stylize(style = 1) %>%
          cols_align(align = "left") %>%
          data_color(
            columns = "-log10(FDR) of representative GO",
            palette = "Greens") %>%
          data_color(
            columns = "number of genes in reduced GO-terms",
            palette = "Blues") %>%
          opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
      }
    }, ignoreNULL = FALSE)
    
    output$FunEnrichaBSRELClustersTable <- render_gt({
      FunEnrichClustTable()
    })
    
  # retrieve which gene has been selected
  get_genename <- reactive({
    
    aBSREL_input <- input$aBSRELGeneInput
    
    # return nothing when no gene has been chosen
    if (aBSREL_input == "" ) {
      return("")
    }
    
    # select using transcript ID, if that option has been set
    if (input$aBSRELTranscriptIDSelect) {
      all_genes$names[all_genes$ensembleTrans == aBSREL_input]
    } else {
      all_genes$names[all_genes$names == aBSREL_input]
    }
  })

  output$aBSRELTableDescription <- renderUI({
    
    genename <- get_genename()
    
    if (genename == "") {
      # If no gene has been selected
      strong("Select a gene to display its results")
    } else {
      # If a gene has been selected
      strong(paste("aBSREL results of", genename))
    }
  })
  
  get_aBSREL_results <- reactive({
    
    genename_chosen <- get_genename()
    
    # Makes sure nothing gets generated when no gene has been selected
    if (genename_chosen == "") {
      return(NULL)
    }
    
    # generation of datatable
    aBSREL_data %>%
      filter(genename == genename_chosen) %>%
      use_series(test_results) %>%
      extract2(1) %>%
      rename("Branch name" = "branch_name",
             "Genome assembly name" = "branch_assembly",
             "TOGA projection ID" = "branch_projection",
             "Likelyhood ratio test statistic" = "LRT",
             "Uncorrected P-value" = "p_value", "FDR-adjusted P-value" = "p_value_cor",
             "Uncorrected P-value significant (P-value < 0.05)" = "significant",
             "FDR-adjusted P-value significant (P-value < 0.05)" = "significant_cor") %>%
      # remove non-significant branches of input$aBSRELTableSignBranches is TRUE
      {if (input$aBSRELTableSignBranches) 
        filter(., `FDR-adjusted P-value significant (P-value < 0.05)` == TRUE) else .} %>%
      # change assembly names to species names
      mutate(`Branch name` = map_chr(`Branch name`, ~ {
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else if (!(.x %in% name_conversion$`Assembly name`)) {
          return(.x)
        } else if (.x %in% name_conversion$`Assembly name`) {
          index <- which(.x == name_conversion$`Assembly name`)
          species_tree <- name_conversion$Species_Tree[index]
          return(species_tree)
        } else {
          cat("ERROR! non-recognized branch name!")
          return(NULL)
        }
      }))
  })
  
  get_aBSREL_results_length <- reactive({
    
    nrow(get_aBSREL_results())
  })
  
  output$aBSREL_table <- DT::renderDataTable({
    
    # Makes sure no table gets generated when no gene has been chosen
    if (get_genename() == "") {
      return(tibble())
    }
    
    get_aBSREL_results()
    
  }, escape = FALSE, rownames = FALSE, 
  options = list(scrollY = "500px", scrollX = TRUE, 
                 lengthMenu = c(50, 100, 250, get_aBSREL_results_length()))
  )
  
  ### Backend for MEME results page
  
  observe({
    
    freezeReactiveValue(input, "MEMEGeneInput")
    
    if (input$MEMETranscriptIDSelect) {
      
        updateSelectizeInput(session, "MEMEGeneInput", choices = MEME_transcript_ids, server = TRUE, 
                             options = list(placeholder = "Transcript id"))
    } else {
        updateSelectizeInput(session, "MEMEGeneInput", choices = MEME_genenames, server = TRUE, 
                             options = list(placeholder = "Gene name"))
    }
  })
  
  # retrieve which gene has been selected
  get_genename_MEME <- reactive({
    
    # return nothing when no gene has been chosen
    if (input$MEMEGeneInput == "") {
      return("")
    }
    
    if (input$MEMETranscriptIDSelect) {
      all_genes$names[all_genes$ensembleTrans == input$MEMEGeneInput]
    } else {
      all_genes$names[all_genes$names == input$MEMEGeneInput]
    }
  })
  
  ### MEME backend for 'test results' tabpanel
  
  output$MEMETableDescription <- renderUI({
    
    if (get_genename_MEME() == "") {
      # If no gene has been selected
      paragraph <- strong("Select a gene to display its results")
      paragraph
    } else {
      # If a gene has been selected
      paragraph <- paste0(strong(paste("MEME results of", get_genename_MEME())))
      HTML(paragraph)
    }
  })
  
  # extract the MEME_data row from the selected gene
  get_MEME_data_row <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    # load the required RData object
    load(file.path("data", "MEME_data",
                   paste0("MEME_data_", get_genename_MEME(), ".RData")))
    return(RData_row)
  })
  
  get_MEME_results <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }

    MEME_row <- get_MEME_data_row()
    
    MEME_row %>%
      use_series(csv) %>%
      extract2(1) %>%
      # add an additional site column
      add_column(site = seq(nrow(.))) %>%
      relocate(site) %>%
      # remove non-significant sites if input$MEMETableSignSites is TRUE
      {if (input$MEMETableSignSites) filter(., `p-value` <= input$p_val_select) else .} 
  })
  
  get_MEME_results_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_results())
  })
  
  get_unimapped_MEME_results <- reactive({
    
    # Makes sure no table gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_subs %>%
      dplyr::filter(genename == get_genename_MEME()) %>%
      separate_longer_delim(c(sites, branch_subs, branch_subs_ancestral, `&alpha;`,
                              `&beta;<sup>-</sup>`, `p<sup>-</sup>`,
                              `&beta;<sup>+</sup>`, `p<sup>+</sup>`, `LRT`, `p-value`,
                              pval_fdr, `#_branches_under_selection`, `Total_branch_length`,
                              `MEME_LogL`, `FEL_LogL`, `Variation_p`),
                            delim = stringr::regex("[;,]")) %>%
      # rename columns
      dplyr::rename(site = sites) %>%
      ## format some columns as doubles and integers
      mutate(site = as.integer(site),
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
      # remove non-significant sites if input$MEMETableSignSites is TRUE
      {if (input$MEMETableSignSites) filter(., `p-value` <= input$p_val_select) else .} 
  })
  
  output$MEME_table <- DT::renderDataTable({
    
    # Makes sure no table gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    get_MEME_results()
    
  # escape = FALSE allows HTML tags to be used in the table
  # options = list(scrollY = "500px", scrollX = TRUE) shortens the table and 
  # scrollX allows scrolling in the horizontal direction
  # lengthMenu picks the selectable lengths of the datatable
    
  }, escape = FALSE, rownames = FALSE, 
  options = list(scrollY = "500px", scrollX = TRUE, 
                 lengthMenu = c(50, 100, 250, get_MEME_results_length()))
  )
  
  # render block of summarizing text 
  output$MEMEResultDescription <- renderUI({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }

    # filter with p-value <= 0.1 if input$MEMETableSignSites is FALSE
    filtered_results <- get_MEME_results() %>%
      {if (input$MEMETableSignSites == FALSE) filter(., `p-value` <= 0.1) else .}
      
    # get number of significant sites
    n_of_sites <- filtered_results %>%  
      use_series(site) %>%
      length()
    
    # get significant site location
    sites <- filtered_results %>%
      use_series(site) %>%
      as.character() %>%
      paste(collapse = ", ")
    
    paragraph <- paste0(h4(HTML(paste(em(get_genename_MEME()), "result summary:"))), 
                        p(HTML(paste(get_genename_MEME(), "was found to have", strong(n_of_sites),
                                strong("codon sites"), "with significant signs (p <= ", input$p_val_select, ") of episodic diversifying",
                                "selection"))),
                        p(HTML(paste("These codon sites are", strong(sites), ".")))
                        )
    HTML(paragraph)
  })
  
  # Rendering table for uniprot-mapped sites
  output$MEME_table_unimapped <- DT::renderDataTable({
    
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(NULL)
    }
    
    get_unimapped_MEME_results() %>%
      dplyr::select(-all_of(c("site_included",
                              "uniprot_gn_symbol", "ensembl_gene_id",
                              "pval_fdr", "n_of_sites", "genename",
                              "transcript_id"))) %>%
      relocate(branch_subs, branch_subs_ancestral, .after = uniprotswissprot) %>%
      rename("Amino acid difference between Human and A. cahirinus sequence" = "branch_subs",
             "Amino acid substitution that took place in the A. cahirinus lineage (i.e. branch)" = "branch_subs_ancestral")
    }, escape = FALSE, rownames = FALSE, 
    options = list(scrollY = "500px", scrollX = TRUE, 
                   lengthMenu = c(25, 50, 100)))
  
  # render block of summarizing text for uniprot-mapped sites 
  output$MEMEResultDescription_unimapped <- renderUI({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(HTML(paste0(h4(HTML(paste("No sites could be mapped to the uniprot accession of gene",
                           em(get_genename_MEME())))))))
    }
    
    # filter with p-value <= 0.1 if input$MEMETableSignSites is FALSE
    filtered_results <- get_unimapped_MEME_results() %>%
      {if (input$MEMETableSignSites == FALSE) filter(., `p-value` <= 0.1) else .}
    
    # get number of significant sites
    n_of_sites <- filtered_results %>%  
      use_series(site) %>%
      length()
    
    # get significant site location
    sites <- filtered_results %>%
      use_series(site) %>%
      as.character() %>%
      paste(collapse = ", ")
    
    uniprot_acc <- filtered_results %>%
      use_series(uniprotswissprot) %>%
      .[1]
    
    paragraph <- paste0(h4(HTML(paste(em(get_genename_MEME()), "result summary:"))), 
                        p(HTML(paste(get_genename_MEME(), "was found to have", strong(n_of_sites),
                                     strong("codon sites"), "with significant signs (p <= ", input$p_val_select, ") of episodic diversifying",
                                     "selection that we unambigiously mapped to the canonical protein sequence of",
                                     get_genename_MEME(), paste0("(", uniprot_acc, ")")
                                     ))),
                        p(HTML(paste("These codon sites are", strong(sites), ".")))
    )
    HTML(paragraph)
  })
  
  ### MEME backend for 'Multiple sequence alignment' tabpanel
  
  output$MultipleAlignmentText <- renderUI({
    if (get_genename_MEME() == "") {
      return("Please select a gene to display the multiple sequence alignment")
    }
    if (!(get_genename_MEME() %in% MSA_genelist)) {
      return("This gene was not mapped to a uniprot entry, and consequently no MSA was generated for this gene. Please select a different gene to display the multiple sequence alignment")
    } else {
      return(NULL)
    }
  })
  
  # Render the multiple sequence alignment
  output$MultipleAlignment <- renderUI({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    if (!(get_genename_MEME() %in% MSA_genelist)) {
      return(NULL)
    }
    
    # create URL
    if (input$MSAradioselect == 1) {
      gene_URL <- paste0("MEME_protein_alignments/", get_genename_MEME(), "_alignment.pdf")
    } else if (input$MSAradioselect == 2) {
      gene_URL <- paste0("MEME_protein_alignments_charge/", get_genename_MEME(), "_alignment.pdf")
    } else if (input$MSAradioselect == 3) {
      gene_URL <- paste0("MEME_protein_alignments_structure/", get_genename_MEME(), "_alignment.pdf")
    }
    
    # create pdf iframe
    tags$iframe(style="height:900px; width:100%; scrolling=yes", 
                src = gene_URL)
  })
  
  ### MEME backend for 'Substitutions' tabpanel
  
  # Get acomys substitution data from MEME_data
  get_MEME_acomys_sub <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_row <- get_MEME_data_row()
    
    MEME_row <- MEME_row %>%
      use_series(substitutions) %>%
      extract2(1)
    
    HLacoCah2_colname <- colnames(MEME_row)[str_detect(colnames(MEME_row), "HLacoCah2")]
    
    MEME_row %>%
      rename("Acomys_cahirinus" = all_of(HLacoCah2_colname))
  })
  
  get_MEME_acomys_sub_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_acomys_sub())
  })
  
  output$MEME_acomys_sub <- DT::renderDataTable({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    get_MEME_acomys_sub()
  }, rownames = FALSE, options = list(scrollY = "300px", scrollX = TRUE, 
                                      lengthMenu = c(get_MEME_acomys_sub_length(), 20))
  )
  
  get_MEME_branch_subs <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_row <- get_MEME_data_row()
    
    MEME_row %>%
      use_series(branch_subs) %>%
      extract2(1) %>%
      select(!node) %>%
      rename("Branch name" = "label", "Substitutions" = "subs") %>%
      # change assembly names to species names
      mutate(`Branch name` = map_chr(`Branch name`, ~ {
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else if (!(.x %in% name_conversion$`Assembly name`)) {
          return(.x)
        } else if (.x %in% name_conversion$`Assembly name`) {
          index <- which(.x == name_conversion$`Assembly name`)
          species_tree <- name_conversion$Species_Tree[index]
          return(species_tree)
        } else {
          cat("ERROR! non-recognized branch name!")
          return(NULL)
        }
    }))
  })
  
  get_MEME_branch_subs_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_branch_subs())
  })
  
  output$MEME_branch_subs <- DT::renderDataTable({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    get_MEME_branch_subs()
  }, rownames = FALSE, options = list(scrollY = "500px", scrollX = TRUE, 
                                      lengthMenu = c(get_MEME_branch_subs_length(), 20))
  )
  
  ### MEME backend for 'EBF' tabpanel
  
  # Load and process the EBF data for the selected gene for the bubbleplot
  get_EBF_data_bubbleplot <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    if (!file.exists(file.path("data", "EBF_data_uniprot_bubbleplot",
                              paste0("EBF_data_up_bubble_", get_genename_MEME(), ".RData")))) {
      return("no_file")
    }
    
    sign_level <- 0.05
    
    ## load the required RData object
    load(file.path("data", "EBF_data_uniprot_bubbleplot",
                   paste0("EBF_data_up_bubble_", get_genename_MEME(), ".RData")))
    
    ## Process the EBF data to a longer format
    
    # Create semicolon separated site column
    RData_row <- RData_row %>%
      rowwise() %>%
      mutate(site = paste(seq(length(strsplit(`&alpha;`, ";")[[1]])), collapse = ";")) %>%
      ungroup()
    
    # get EBF table
    EBF_data_table <- RData_row %>%
      use_series(EBF_table) %>%
      extract2(1)
    
    # filter for selected assemblies
    EBF_data_table %<>% filter(branch %in% c("REFERENCE", "HLacoCah2", "HLpsaObe1",
                                        "HLmerUng1", "HLratNor7", "rn6", "HLmusPah1",
                                        "HLmusCar1", "mm10", "mm39", "HLmesAur2",
                                        "mesAur1", "HLcriGri3", "HLsigHis1",
                                        "HLonyTor1", "HLperManBai2", "HLondZib1",
                                        "HLellLut1"))
    
    # Add the MLE site information the the plot
    RData_row_long <- RData_row %>%
      dplyr::select(-any_of(c("EBF_table", "transcript_id", "uniprotswissprot",
                              "Sequence", "ensembl_gene_id", "genename"))) %>%
      separate_longer_delim(c(site, `&alpha;`, `&beta;<sup>-</sup>`,
                              `p<sup>-</sup>`, `&beta;<sup>+</sup>`, `p<sup>+</sup>`,
                              `LRT`, `p-value`, `#_branches_under_selection`, `Total_branch_length`,
                              `MEME_LogL`, `FEL_LogL`, `Variation_p`),
                            delim = ";") %>%
      # change character to doubles or integers
      mutate(site = as.integer(site),
             `&alpha;` = as.double(`&alpha;`),
             `&beta;<sup>-</sup>` = as.double(`&beta;<sup>-</sup>`),
             `p<sup>-</sup>` = as.double(`p<sup>-</sup>`),
             `&beta;<sup>+</sup>` = as.double(`&beta;<sup>+</sup>`),
             `p<sup>+</sup>` = as.double(`p<sup>+</sup>`),
             LRT = as.double(LRT),
             `p-value` = as.double(`p-value`),
             `#_branches_under_selection` = as.integer(`#_branches_under_selection`),
             Total_branch_length = as.double(Total_branch_length),
             MEME_LogL = as.double(MEME_LogL),
             FEL_LogL = as.double(FEL_LogL),
             Variation_p = as.double(Variation_p))
    
    EBF_data_table_MLE <- EBF_data_table %>%
      mutate(site = as.integer(site)) %>%
      left_join(RData_row_long)
    
    # Create column that indicates if site is significant
    EBF_data_table_MLE %<>%
      mutate(sign = `p-value` <= sign_level) %>%
      mutate(site_sign = if_else(as.integer(sign) == 1, site, NA_integer_))
    
    # Get sites that have EBF 100 in the A. cahirinus branch and that are
    # significant
    EBF100_sites <- EBF_data_table_MLE %>%
      filter(EBF >= 100 & sign == TRUE & branch == "HLacoCah2") %>%
      pull(site) %>%
      unique()
    
    # Create column that indicates which sites are both significant AND have
    # A. cahirinus branch EBF >= 100.
    EBF_data_table_MLE %<>%
      mutate(EBF_100 = EBF >= 100 & sign == TRUE,
             sign_EBF = site %in% EBF100_sites,
             sign_EBF_site = if_else(sign_EBF, site, NA_integer_))
    
    # Perform a log transformation, making the EBF data more normal, improving
    # the visualization.
    EBF_data_table_MLE %<>%
      mutate(EBF = log10(EBF))
    
    # replace assembly names to species_tree names
    EBF_data_table_MLE %<>%
      mutate(branch = map_chr(branch, ~ {
        
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else {
          index <- which(name_conversion$`Assembly name` == .x)
          return(name_conversion$Species_Tree[index])
        }
      }))
    return(EBF_data_table_MLE)
  })
  
  output$EBFBubbleplotDesc <- renderUI({
    if (is.null(get_EBF_data_bubbleplot())) {
      return("Select a gene to generate the EBF plot")
    }
    
    if (all(get_EBF_data_bubbleplot() == "no_file")) {
      return("Sites could not be mapped to a uniprot entry for this gene. Thus, the EBF values could also not be mapped. Please select a different gene")
    }
  })
  
  # Load and process the EBF data for the selected gene for the table
  get_EBF_data_table <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    if (!file.exists(file.path("data", "EBF_data_uniprot",
                               paste0("EBF_data_up_", get_genename_MEME(), ".RData")))) {
      return("no_file")
    }
    
    ## load the required RData object
    load(file.path("data", "EBF_data_uniprot",
                   paste0("EBF_data_up_", get_genename_MEME(), ".RData")))
    
    ## Process the EBF data to a longer format
    
    # Create semicolon separated site column
    RData_row <- RData_row %>%
      rowwise() %>%
      mutate(site = paste(seq(length(strsplit(`&alpha;`, ";")[[1]])), collapse = ";")) %>%
      ungroup()
    
    # get EBF table
    EBF_data_table <- RData_row %>%
      use_series(EBF_table) %>%
      extract2(1)
    
    # Filter EBF_table for leaf nodes only, using the name conversion df
    EBF_data_table %<>% filter(branch %in% c("REFERENCE", "HLacoCah2", "HLpsaObe1",
                                                "HLmerUng1", "HLratNor7", "rn6", "HLmusPah1",
                                                "HLmusCar1", "mm10", "mm39", "HLmesAur2",
                                                "mesAur1", "HLcriGri3", "HLsigHis1",
                                                "HLonyTor1", "HLperManBai2", "HLondZib1",
                                                "HLellLut1"))
    
    # Perform a log transformation, making the EBF data more normal, improving
    # the visualization.
    EBF_data_table %<>%
      mutate(log10_EBF = log10(EBF))
    
    # replace assembly names to species_tree names
    EBF_data_table %<>%
      mutate(branch = map_chr(branch, ~ {
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else {
          index <- which(name_conversion$`Assembly name` == .x)
          return(name_conversion$Species_Tree[index])
        }
      }))
    return(EBF_data_table)
  })
  
  output$EBFTableDesc <- renderUI({
    if (is.null(get_EBF_data_table())) {
      return("Select a gene to generate the EBF table")
    }
    
    if (all(get_EBF_data_table() == "no_file")) {
      return("Sites could not be mapped to a uniprot entry for this gene. Thus, the EBF values could also not be mapped. Please select a different gene")
    }
  })
  
  # render bubbleplot
  output$EBF_bubbleplot <- renderPlotly({
    
    if (is.null(get_EBF_data_bubbleplot())) {
      return(NULL)
    }
    
    if (all(get_EBF_data_bubbleplot() == "no_file")) {
      return(NULL)
    }
    
    # indicate order of assemblies
    assemblies_order <- c("REFERENCE", "HLacoCah2", "HLpsaObe1",
                          "HLmerUng1", "HLratNor7", "rn6", "HLmusPah1",
                          "HLmusCar1", "mm10", "mm39", "HLmesAur2",
                          "mesAur1", "HLcriGri3", "HLsigHis1", "HLonyTor1",
                          "HLperManBai2", "HLondZib1", "HLellLut1")
    
    species_order <- unique(
      sapply(assemblies_order, function(assembly) {
        if (assembly == "REFERENCE") {
          return("Homo_sapiens")
        } else {
          assembly_index <- which(name_conversion$`Assembly name` == assembly)
          return(name_conversion$Species_Tree[assembly_index])
        }
      }))
    
    # remove infinite rows
    EBF_table_MLE <- get_EBF_data_bubbleplot() %>%
      filter(!(is.infinite(EBF)))
    
    EBF_table_MLE <- EBF_table_MLE %>%
      mutate(EBF_100 = if_else(EBF_100 == FALSE, "MEME LRT p-value > 0.05\nor\nEBF < 100",
                               "MEME LRT p-value <= 0.05\nand\nEBF >= 100")) %>%
      filter(EBF >= 0)
    
    # Generate data that will be hidden in the plot, for the purpose of keeping
    # the plot range identical when hiding/selecting legend items
    rangex <- range(EBF_table_MLE$site)
    ranges <- expand(tibble(x = rep(rangex, length.out = length(species_order)),
                            y = species_order),
                     crossing(x, y))
    
    # Generate plotly
    plot_ly() %>%
      add_trace(data = filter(EBF_table_MLE, EBF_100 == "MEME LRT p-value > 0.05\nor\nEBF < 100"),
                name = "MEME LRT p-value > 0.05\nor\nEBF < 100",
                x = ~site,
                y = ~branch,
                color = ~EBF_100,
                colors = rev(as.character(met.brewer("OKeeffe2", type = "continuous", n = 2))),
                type = "scatter",
                mode = "markers",
                size = ~EBF,
                sizes = c(0.0001, 200),
                text = ~paste("<b>log10(EBF):</b> ", EBF)) %>%
      add_trace(data = filter(EBF_table_MLE, EBF_100 == "MEME LRT p-value <= 0.05\nand\nEBF >= 100"),
                name = "MEME LRT p-value <= 0.05\nand\nEBF >= 100",
                x = ~site,
                y = ~branch,
                color = ~EBF_100,
                colors = rev(as.character(met.brewer("OKeeffe2", type = "continuous", n = 2))),
                type = "scatter",
                mode = "markers",
                size = ~EBF,
                sizes = c(0.0001, 200),
                text = ~paste("<b>log10(EBF):</b> ", EBF)) %>%
      # Add an additional hidden trace with the same ranges as the x and y-axes,
      # so that the plot does not resize when a cluster is selected/deselected
      add_trace(data = ranges,
                x = ~x,
                y = ~y,
                type = "scatter",
                mode = "markers",
                showlegend = FALSE,
                opacity = 0,
                # Stops the data from showing up while hovering over the points
                hoverinfo='skip') %>%
      layout(shapes = lapply(unique(EBF_table_MLE$sign_EBF_site[!is.na(EBF_table_MLE$sign_EBF_site)]), function(i) {
        vline(i)
      }),
      yaxis = list(title = "", categoryorder = "array",
                   categoryarray = species_order),
      xaxis = list(title = "Codon site"),
      legend = list(itemsizing = "constant"))
  })
  
  MEME_EBF_inf_remover <- reactive({
    
    if (is.null(get_EBF_data_table()) | all(get_EBF_data_table() == "no_file")) {
      return(NULL)
    }
    get_EBF_data_table() %>%
      {if (input$MEME_EBF_inf_remove) filter(., !(is.infinite(EBF))) else .}
  })
  
  MEME_EBF_acomys_selecter <- reactive({
    
    if (is.null(get_EBF_data_table()) | all(get_EBF_data_table() == "no_file")) {
      return(NULL)
    }
    MEME_EBF_inf_remover() %>%
      {if (input$MEME_EBF_acomys_select) filter(., branch == "Acomys_cahirinus") else .}
  })
  
  # render datatable
  output$EBF_table <- render_gt({
    
    if (is.null(get_EBF_data_table()) | all(get_EBF_data_table() == "no_file")) {
      return(NULL)
    }
    gt(MEME_EBF_acomys_selecter(), groupname_col = "branch",
       row_group_as_column = TRUE) %>%
      cols_align(align = "left") %>%
      opt_stylize(style = 1) %>%
      data_color(
        columns = EBF,
        palette = "Blues") %>%
      data_color(
        columns = log10_EBF,
        palette = "Greens") %>%
      opt_interactive(use_resizers = TRUE)
  })
  
  ### MEME backend for 'Uniprot sequence annotations' tab
  
  output$DomainPageDesc <- renderUI({
    
    if (get_genename_MEME() == "") {
      return(h4("Select a gene to visualize the MEME-identified positively selected sites and their corresponding Uniprot annotations"))
    }
    
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(p("MEME-identified positively selected sites could not be mapped",
               "to a Uniprot sequence for gene", get_genename_MEME(),
               ". Please select a different gene in order to display positively",
               "selected sites and their corresponding Uniprot feature annotations"))
    }
    
    h4(paste0("Uniprot sequence annotations for ", get_genename_MEME()))
  })
  
  get_uniprot_anno_dfs <- reactive({
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    # Return NULL when gene was not mapped to Uniprot
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(NULL)
    }
    
    # Get unimapped MEME substitutions
    up_anno_subs_df <- get_unimapped_MEME_results()
    # Get unimapped MEME EBF values
    up_anno_EBF_df <- get_EBF_data_table()
    
    ### Filter the EBF_data to only include HLacoCah2 EBF values
    up_anno_EBF_df <- up_anno_EBF_df %>%
      dplyr::filter(branch == "Acomys_cahirinus")
    
    ### Add EBF (empirical bayes factor values to up_anno_subs_df)
    up_anno_subs_df <- up_anno_EBF_df %>%
      dplyr::select(-branch) %>%
      right_join(up_anno_subs_df)
    
    ### Filter for for the selected EBF cutoff
    if (is.na(input$UniDomainEBFSelect)) {
      up_anno_subs_df <- up_anno_subs_df %>%
        dplyr::filter(EBF >= 0)
    } else {
      up_anno_subs_df <- up_anno_subs_df %>%
        dplyr::filter(EBF >= input$UniDomainEBFSelect)
    }
    
    # Filter progo_out for selected gene
    progo_out <- progo_out %>%
      dplyr::filter(Gene == get_genename_MEME())
    
    ### Add feature length
    progo_out <- progo_out %>%
      mutate(length = End - Start)
    
    if (nrow(up_anno_subs_df) != 0) {
      ### Link the protein features to the sites
      up_anno_subs_df <- link_protein_features(up_anno_subs_df, progo_out)
    } else {
      up_anno_subs_df$feat_name <- NA
      up_anno_subs_df$feat_type <- NA
      up_anno_subs_df$Start <- NA
      up_anno_subs_df$End <- NA
    }
    
    return(list(up_anno_subs_df, progo_out))
    })
  
  output$UniDomainPlot <- renderPlotly({
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    # Return NULL when gene was not mapped to Uniprot
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(NULL)
    }
    
    ### Process progo_out to input for uniprot annotation plot
    
    # Filter and rename columns
    draw_features <- get_uniprot_anno_dfs()[[2]] %>%
      dplyr::select(any_of(c("Cluster", "feat_type", "feat_name", "Start", "End", "length",
                             "uniprotswissprot", "Gene"))) %>%
      dplyr::rename(type = "feat_type", description = "feat_name", begin = "Start",
                    end = "End", entryName = "Gene")
    
    # In the description column, replace NA with NONE
    draw_features <- draw_features %>%
      mutate(description = if_else(is.na(description), "NONE", description))
    
    ### Create plotly input dataframe with positively selected sites
    ### (annotated sites)
    
    if (nrow(get_uniprot_anno_dfs()[[1]]) != 0) {
      # Create custom ANNO data for the selected cluster
      ANNO_data <- get_uniprot_anno_dfs()[[1]] %>%
        dplyr::distinct(genename, site, .keep_all = TRUE) %>%
        # Get column that identifies if sites are in the features to be highlighted
        identify_feature_sites(draw_features,
                               selected_features = c("DOMAIN",
                                                     "REGION",
                                                     "MOTIF")) %>%
        dplyr::select(any_of(c("Cluster", "site", "feat_type", "feat_name", "Start", "End",
                               "uniprotswissprot", "genename", "in_feature"))) %>%
        dplyr::rename(type = "feat_type", description = "feat_name", begin = "Start",
                      end = "End", accession = "uniprotswissprot",
                      entryName = "genename") %>%
        # Modify data to ANNO data
        mutate(type = "ANNO",
               description = as.character(site),
               begin = as.integer(site),
               end = as.integer(site),
               length = 0)
      
      ANNO_data <- ANNO_data %>%
        mutate(in_feature = if_else(in_feature == TRUE, "Positively selected site\nLocated inside Domain, Region or Motif",
                                    "Positively selected site\nLocated outside Domain, Region or Motif"))
    }

    ### Separate CHAIN and SIGNAL features from DOMAIN, REGION and MOTIF features
    draw_chain <- draw_features %>%
      dplyr::filter(type == "CHAIN") %>%
      # For some genes multiple chains are given. For example, a chain for the
      # full-length protein and a chain for the C-terminally truncated protein.
      # In these cases we select the longest chain for visualization 
      dplyr::arrange(desc(length)) %>%
      dplyr::slice(1)
    
    draw_signal <- draw_features %>%
      dplyr::filter(type == "SIGNAL")
    
    draw_features <- draw_features %>%
      dplyr::filter(type %in% c("DOMAIN",
                         "REGION",
                         "MOTIF"))
    
    ### Determine the color mappings for the Uniprot annotations.
    ### We link a different color to each different Uniprot annotation,
    ### making sure that essentially identical annotations (Like Cadherin,
    ### Cadherin 1, Cadherin 2, Cadherin 3, ext.) get the same color.
    
    # Get and simplify the Uniprot feature names
    present_features <- unique(draw_features$description)
    present_features <- setNames(present_features,
                                        Translate_UP_features(present_features))
    
    # Get a color for each essentially different feature
    shape_pal <- rev(as.character(met.brewer("Manet",type = "continuous",
                                             n = length(unique(names(present_features))))))
    shape_pal <- setNames(shape_pal, unique(names(present_features)))
    
    # Map the colors determined for each feature to each separate Uniprot feature
    shape_pal <- sapply(names(present_features), function(feat) {
      unname(shape_pal[names(shape_pal) == feat])
    })
    names(shape_pal) <- present_features
    
    pal <- c("#551F00", "#32B2DA")
    pal <- setNames(pal, c("Positively selected site\nLocated inside Domain, Region or Motif",
                           "Positively selected site\nLocated outside Domain, Region or Motif"))
    
    ### Create plotly. The if-statement makes sure that the signal sequence is
    ### only drawn when the annotation is present
    if (nrow(draw_signal) == 0) {
      p <- plot_ly() %>%
        layout(shapes = c(list(ly_rect(x0 = draw_chain$begin[1],
                                       x1 = draw_chain$end[1],
                                       y0 = 0.90,
                                       y1 = 1.10,
                                       fillcolor = "gray")),
                          lapply(seq(nrow(draw_features)), function(i) {
                            ly_rect(x0 = draw_features$begin[i],
                                    x1 = draw_features$end[i],
                                    y0 = 0.75,
                                    y1 = 1.25,
                                    fillcolor = shape_pal[names(shape_pal) == draw_features$description[i]])
                          })),
               yaxis = list(visible = FALSE),
               xaxis = list(title = "<b>Amino acid position</b>",
                            showgrid = FALSE, zeroline = FALSE))
    } else {
      p <- plot_ly() %>%
        layout(shapes = c(list(ly_rect(x0 = draw_chain$begin[1],
                                       x1 = draw_chain$end[1],
                                       y0 = 0.90,
                                       y1 = 1.10,
                                       fillcolor = "gray"),
                               ly_rect(x0 = draw_signal$begin[1],
                                       x1 = draw_signal$end[1],
                                       y0 = 0.90,
                                       y1 = 1.10,
                                       fillcolor = "red")),
                          lapply(seq(nrow(draw_features)), function(i) {
                            ly_rect(x0 = draw_features$begin[i],
                                    x1 = draw_features$end[i],
                                    y0 = 0.75,
                                    y1 = 1.25,
                                    fillcolor = shape_pal[names(shape_pal) == draw_features$description[i]])
                          })),
               yaxis = list(visible = FALSE),
               xaxis = list(title = "<b>Amino acid position</b>",
                            showgrid = FALSE, zeroline = FALSE))
    }
    
    if (nrow(get_uniprot_anno_dfs()[[1]]) != 0) {
      # Draw sites on plotly
      p <- p %>% add_trace(data = ANNO_data,
                           type = "scatter",
                           mode = "markers",
                           color = ~in_feature,
                           colors = pal,
                           text = ~paste0("<b>Amino acid position:</b> ", begin, "<br>",
                                          "<b>Uniprot accession:</b> ", accession),
                           hoverinfo = c("text"),
                           opacity = 1,
                           x = ~begin,
                           y = 1.25)
    }
    
    for (i in seq(nrow(draw_features))) {
      p <- p %>% add_trace(type = "scatter",
                           mode = "markers",
                           x = c(draw_features$begin[i],
                                 draw_features$begin[i],
                                 draw_features$end[i],
                                 draw_features$end[i],
                                 draw_features$begin[i]),
                           y = c(0.75, 1.25, 1.25, 0.75, 0.75),
                           fill = "toself",
                           hoverlabel = list(bgcolor = shape_pal[names(shape_pal) == draw_features$description[i]]),
                           text = paste0("<b>Name:</b> ", draw_features$description[i], "<br>",
                                         "<b>Uniprot annotation type:</b> ", draw_features$type[i], "<br>",
                                         "<b>Start codon:</b>", draw_features$begin[i], "<br>",
                                         "<b>End codon:</b>", draw_features$end[i]),
                           hoverinfo = "text",
                           showlegend = FALSE,
                           opacity = 0,
                           name = draw_features$description[i])
    }
    
    return(p)
    })
  
  
  output$UniDomainSiteDesc <- renderUI({
    # Return NULL when no gene was selected or when gene was not mapped to Uniprot
    if ((get_genename_MEME() == "") || (nrow(get_unimapped_MEME_results()) == 0)) {
      return(NULL)
    }
    
    main_text <- strong(paste("Table describing MEME-identified positively selected",
                        "sites for gene", get_genename_MEME(), ", that are linked to their respective Uniprot sequence",
                        "annotations"))
    
    if (nrow(get_uniprot_anno_dfs()[[1]]) == 0) {
      return(div(p(main_text),
                 p(paste("No MEME-identified positively selected sites for gene",
                          get_genename_MEME(),
                          "remain under the current EBF and/or p-value settings")))
             )
    }
    div(p(main_text))
    
  })

  output$UniDomainSiteTable <- render_gt({
    # Return NULL when no gene is selected, when gene was not mapped to a
    # Uniprot sequence or when no positively selected sites remain under the
    # current p-value/EBF settings
    if ((get_genename_MEME() == "") ||
        (nrow(get_uniprot_anno_dfs()[[1]]) == 0) ||
        (nrow(get_unimapped_MEME_results()) == 0)) {
      return(NULL)
    }
    
    get_uniprot_anno_dfs()[[1]] %>%
      dplyr::select(-all_of(c("n_of_sites", "Total_branch_length", "pval_fdr", "uniprot_gn_symbol", "site_included"))) %>%
      dplyr::rename("log10(EBF)" = log10_EBF,
                    "Ensembl transcript ID" = transcript_id,
                    "gene" = genename,
                    "Amino acid difference between Human and A. cahirinus sequence" = branch_subs,
                    "Likelyhood ratio test statistic" = LRT,
                    "Ensembl gene ID" = ensembl_gene_id,
                    "Amino acid substitution that took place in the A. cahirinus lineage (i.e. branch)" = branch_subs_ancestral,
                    "Uniprot accession" = uniprotswissprot,
                    "Uniprot feature type" = feat_name,
                    "Uniprot feature name" = feat_type,
                    "Feature Start site" = Start,
                    "Feature End site" = End) %>%
    gt(groupname_col = "gene") %>%
      opt_stylize(style = 1) %>%
      cols_align(align = "left") %>%
      tab_spanner(label = html("<strong><em>Gene info</em></strong>"),
                  columns = c(site, gene, `Ensembl transcript ID`, `Ensembl gene ID`)) %>%
      tab_spanner(label = html("<strong><em>Uniprot feature annotation</em></strong>"),
                  columns = c(`Uniprot accession`, `Uniprot feature type`,
                              `Uniprot feature name`, `Feature Start site`,
                              `Feature End site`)) %>%
      tab_spanner(label = html("<strong><em>A. cahirinus EBF</em></strong>"),
                  columns = c(EBF, `log10(EBF)`)) %>%
      tab_spanner(label = html("<strong><em>MEME maximum likelyhood estimation output</em></strong>"),
                  columns = c(`&alpha;`, `&beta;<sup>-</sup>`, `p<sup>-</sup>`,
                              `&beta;<sup>+</sup>`, `p<sup>+</sup>`,
                              `Likelyhood ratio test statistic`, `p-value`,
                              `#_branches_under_selection`, MEME_LogL,
                              FEL_LogL, Variation_p)) %>%
      tab_spanner(label = html("<strong><em>MEME-predicted substitutions</em></strong>"),
                  columns = c(`Amino acid difference between Human and A. cahirinus sequence`,
                              `Amino acid substitution that took place in the A. cahirinus lineage (i.e. branch)`)) %>%
      cols_move(c(`Uniprot accession`, `Uniprot feature type`,
                  `Uniprot feature name`, `Feature Start site`,
                  `Feature End site`), after = `Ensembl gene ID`) %>%
      cols_move_to_end(c(`Amino acid difference between Human and A. cahirinus sequence`,
                         `Amino acid substitution that took place in the A. cahirinus lineage (i.e. branch)`)) %>%
      cols_label(`&alpha;` = gt::html("&alpha;"),
                 `&beta;<sup>-</sup>` = gt::html("&beta;<sup>-</sup>"),
                 `p<sup>-</sup>` = gt::html("p<sup>-</sup>"),
                 `&beta;<sup>+</sup>` = gt::html("&beta;<sup>+</sup>"),
                 `p<sup>+</sup>` = gt::html("p<sup>+</sup>")) %>%
      data_color(
        columns = `p-value`,
        palette = "Greens",
        reverse = TRUE) %>%
      data_color(
        columns = EBF,
        palette = "Blues") %>%
      opt_interactive(use_compact_mode = TRUE)
  })
  
  output$UniDomainFeatDesc <- renderUI({
    # Return NULL when no gene was selected or when gene was not mapped to Uniprot
    if ((get_genename_MEME() == "") || (nrow(get_unimapped_MEME_results()) == 0)) {
      return(NULL)
    }
    
    main_text <- div(strong(paste0("Table describing the Uniprot sequence annotations ",
                                  "(i.e. features) found for gene ", get_genename_MEME(),
                                  ". If more information is needed about the EMBL-EBI Proteins API,",
                                  " follow this link to the")),
                     strong(a("Proteins API webpage", href="https://www.ebi.ac.uk/proteins/api/doc/#/")),
                     br(),
                     br(),
                     p("Feature type CHAIN is used to draw the protein sequence",
                       "(gray) seen in the above visualization. Potential SIGNAL",
                       "peptides at the start of the protein sequence are",
                       "visualized with the color red"))
    
    return(main_text)
  })
  
  output$UniDomainFeatTable <- render_gt({
    # Return NULL when no gene was selected or when gene was not mapped to Uniprot
    if ((get_genename_MEME() == "") || (nrow(get_unimapped_MEME_results()) == 0)) {
      return(NULL)
    }
    
    get_uniprot_anno_dfs()[[2]] %>%
      dplyr::rename("Feature Start site" = Start,
                    "Feature End site" = End,
                    "Uniprot feature name" = feat_name,
                    "Uniprot feature type" = feat_type,
                    "Feature length" = length) %>%
      gt(groupname_col = "Gene") %>%
      opt_stylize(style = 1) %>%
      cols_align(align = "left") %>%
      data_color(
        columns = `Feature length`,
        palette = "Greens") %>%
      opt_interactive(use_compact_mode = TRUE)
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)