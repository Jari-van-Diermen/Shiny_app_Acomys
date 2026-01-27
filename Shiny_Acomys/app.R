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
library(Biostrings)
library(ggmsa)
library(bslib)
library(ape)
library(ggtree)

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
MEME_subs <- read_delim("data/uniprot_mapping_pairwise/MEME_subs_mapped_sitefiltered_noMSAs.tsv",
                        delim = "\t", escape_double = FALSE,
                        col_types = cols(n_of_sites = col_integer(),
                                         sites = col_character(),
                                         TOGA_sites = col_character()), trim_ws = TRUE)

### Load the uniprot sequence annotations
progo_out <- read_delim("data/Uniprot_annotations/progo_out_combined.tsv",
                        delim = "\t", escape_double = FALSE,
                        trim_ws = TRUE)

### Load the MEME site translation table
MEME_trans_df <- read_delim("data/translation_table/MEME_gene_site_translation_all_positions.tsv",
                            delim = "\t",
                            escape_double = FALSE,
                            trim_ws = TRUE)

MEME_trans_df <- MEME_trans_df %>%
  dplyr::rename(TOGA_Human_site = TOGA_orig_Human_site,
                TOGA_Acomys_site = TOGA_orig_Acomys_site,
                MSA_site = ggmsa_site) 

### Source functions
source("app_functions.R")

### load the genenames and corresponding ensemble IDs for the significant and
### non-significant genes

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

### create coord selection list for interactive switching
coord_options <- list("1" = c("MSA_site", "MSA (in the multiple sequence alignment tab),"),
                      "2" = c("TOGA_Human_site", "Human TOGA sequence,"),
                      "3" = c("TOGA_Acomys_site", "A. cahirinus TOGA sequence,"),
                      "4" = c("UP_Human_site", "Human canonical uniprot sequence,"))

# CSS override (fpr tables without scrollbars)
tags$head(
  tags$style(HTML("
      /* Disable DT internal scrolling for one table */
      .no-dt-scroll .dataTables_scrollBody {
        overflow-y: visible !important;
        max-height: none !important;
      }
    ")))

# Define UI ----
ui <- page_navbar(
  title = HTML("Title to be determined"),
  theme = bs_theme(version = 5, bootswatch = "lux"),
  nav_panel("aBSREL results",
            navset_underline(
              nav_panel("aBSREL LRT results",
                        layout_sidebar(
                          sidebar = sidebar(
                            width = 350,
                            h3("Gene selection"),
                            selectizeInput(inputId = "aBSRELGeneInput", "Select the gene for which you want to see the results", 
                                           choices = NULL),
                            checkboxInput("aBSRELSignSwitch", HTML(paste("Only allow selecting genes with significant signs of episodic",
                                                                         "diversifying selection in", em("Acomys cahirinus"))), value = FALSE),
                            checkboxInput("aBSRELTranscriptIDSelect", "Choose genes based on their transcript IDs", value = FALSE),
                            h3("Table options"),
                            checkboxInput("aBSRELTableSignBranches", "Only show significant Branches", value = FALSE)
                          ),
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
                          card(
                            full_screen = TRUE,
                            card_header(uiOutput("aBSRELTableDescription")),
                            card_body(
                              # Wrapped inside CSS class that disables scrolling
                              div(
                                class = "no-dt-scroll",
                                DT::dataTableOutput("aBSREL_table")
                              )
                            )
                          )
                        )),
              nav_panel("STRING interaction network",
                        layout_sidebar(
                          sidebar = sidebar(
                            h4("Display options for STRING clusters"),
                            checkboxInput("STRINGplotSwitch", HTML("Display all clusters with at least three genes"), value = FALSE),
                            br(),
                            sliderInput("STRINGClustTableSlider", "Cluster number:",
                                        min = 1, max = max(enrichGO_aco$`\`__mclCluster\``,
                                                           na.rm = TRUE),
                                        value = c(1, 8)),
                            checkboxInput("STRINGDisplayNonClust", HTML("Display non-clustered genes in table"), value = FALSE),
                            actionButton("STRINGClustTableAction", "Refresh")
                          ),
                          div(p(strong("This page displays the STRING clusters that were identified in the", em("A. cahirinus"),
                                       "(episodic) positively selected genelist. The clusters were generated using the Markov cluster (MCL) algorithm with an inflation parameter of 4")),
                              p(strong("We identified 8 STRING clusters with at least 5 genes and 64 STRING clusters with the minimum of 2 genes",
                                       "These can be explored using the interactive scatterplot and table below. Furthermore, the desired STRING",
                                       "clusters in the table can be selected using the range buttons on the left side of the page."))
                              
                          ),
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
              nav_panel("Functional overrepresentation analysis",
                        layout_sidebar(
                          sidebar = sidebar(
                            #generate_br_tags(40),
                            h4("Display options for functionel enrichment of STRING clusters"),
                            checkboxInput("FunEnrichReducedSwitch", HTML("Display all GO-terms per cluster instead of only the reduced GO-terms"), value = FALSE),
                            h4("Select which clusters to display"),
                            sliderInput("FunEnrichSelectClust", "Cluster number:",
                                        min = 1, max = max(enrichGO_aco$`\`__mclCluster\``,
                                                           na.rm = TRUE),
                                        value = c(1, 8)),
                            actionButton("FunEnrichClustAction", "Refresh")
                          ),
                          
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
                        )
              )
            )
  ),
  nav_panel("MEME results",
            layout_sidebar(
              sidebar = sidebar(
                width = 350,
                h3("Gene selection"),
                selectizeInput(inputId = "MEMEGeneInput", "Select the gene for which you want to see the results", 
                               choices = NULL),
                checkboxInput("MEMETranscriptIDSelect", "Choose genes based on their transcript IDs", value = FALSE),
                h3("Table options"),
                checkboxInput("MEMETableSignSites", "Only show significant sites", value = FALSE),
                numericInput("p_val_select", "p-value significance cut-off", value = 0.05),
                radioButtons(
                  "MEMEFilterSitesForSubs",
                  label = "Filter for sites with substitutions",
                  choices = list(
                    "All sites" = 1,
                    "Sites with substitution in Homo sapiens" = 2,
                    "Sites with substitution in A. cahirinus" = 3),
                  selected = 1
                )
              ),
              navset_underline(
                nav_panel("MEME LRT results",
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
                          card(
                            full_screen = TRUE,
                            card_header(uiOutput("MEMETableDescription")),
                            card_body(
                              p(HTML(paste("This table and summary data represents the MEME likelyhood ratio test results for each site of the",
                                           "selected gene. Here, the positional site information is relative to the human gene",
                                           "that was used in the TOGA multiple sequence alignment (i.e. the multiple sequence",
                                           "alignment with mammalian orthologs of the selected gene, available at the",
                                           a("TOGA ortholog database", href="https://genome.senckenberg.de//download/TOGA/human_hg38_reference/MultipleCodonAlignments/")
                              ))),
                              gt::gt_output("MEME_table")
                            )),
                          card(
                            full_screen = TRUE,
                            card_header(uiOutput("MEMEResultDescriptionTitle")),
                            card_body(
                              selectInput("MEMEResultCoordSelector", label = "Select reference frame", choices = list("MSA site" = 1,
                                                                                                                      "TOGA Human site" = 2,
                                                                                                                      "TOGA Acomys_site" = 3,
                                                                                                                      "UP Human site" = 4),
                                          selected = 2),
                              uiOutput("MEMEResultDescription"),
                            )),
                ),
                nav_panel("Multiple sequence alignment",
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
                nav_panel("Substitutions",
                          div(p("This page visualizes predicted substitutions at positively selected sites, identified using the HyPhy MEME framework across",
                                "the mammalian species included in the evolutionary model."),
                              p("These substitutions are displayed as a tree view for a selected site, or in the form of an interactive table where all substitutions are displayed among the mammalian species",
                                "The tree view is particularly usefull for identifying where a substitution took place according to the MEME evolutionary model."),
                              p("Substitutions are in the form of A11Q. Substitutions that appear to be indels are indicated with the term 'Gap' (e.g. C-333-Gap)."),
                              p("Here, the TOGA-Human and TOGA-Acomys reference frames refer to the positions in the TOGA-identified Human and", em("A. cahirinus"), "orthologs respectively.",
                                "These orthologs can be found at the", a("TOGA ortholog database.", href="https://genome.senckenberg.de//download/TOGA/"))),
                          card(
                            full_screen = TRUE,
                            card_header(h5("MEME-identified substitutions tree view")),
                            card_body(
                              layout_sidebar(
                                sidebar = sidebar(
                                  selectInput("MEMETreeCoordSelector", label = "Select reference frame", choices = list("MSA site" = 1,
                                                                                                                        "TOGA Human site" = 2,
                                                                                                                        "TOGA Acomys_site" = 3,
                                                                                                                        "UP Human site" = 4),
                                              selected = 2),
                                  selectizeInput(inputId = "SubSiteSelect", "Select site for phylogenetic tree annotation", choices = NULL),
                                ),
                                uiOutput("SubPhyloText"),
                                plotOutput("SubSitePhylo", width = "100%", height = "1400px"),
                              )
                            )
                          ),
                          card(
                            full_screen = TRUE,
                            card_header(h5("All substitutions at positively selected codon sites")),
                            card_body(
                              # Wrapped inside CSS class that disables scrolling
                              div(
                                class = "no-dt-scroll", 
                                DT::dataTableOutput("MEME_branch_subs"),
                              ),
                            )
                          )
                ),
                nav_panel("EBF",
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
                nav_panel("Uniprot sequence annotations",
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
  ),
  nav_spacer(),
  nav_item(bslib::input_dark_mode(mode = "light"))
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
      dplyr::rename(gene_id = ensemble, transcript_id = ensembleTrans) %>%
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
      dplyr::rename(branch_assembly = `Assembly name`) %>%
      right_join(aBSREL_table_omega_test) %>%
      dplyr::select(-branch_assembly)
    
    # Add STRINGdb clusters to data
    aBSREL_table_omega_test_clusters <- acomys_clusters %>%
      dplyr::rename(Cluster = `__mclCluster`, gene_id = `query term`) %>%
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
      dplyr::filter(MRCA_present_aco == TRUE) %>%
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
    
    # Stop if no gene selected
    req(input$aBSRELGeneInput)
    
    genename_chosen <- get_genename()
    
    # generation of datatable
    aBSREL_data %>%
      dplyr::filter(genename == genename_chosen) %>%
      magrittr::use_series(test_results) %>%
      magrittr::extract2(1) %>%
      dplyr::rename("Branch name" = "branch_name",
             "Genome assembly name" = "branch_assembly",
             "TOGA projection ID" = "branch_projection",
             "Likelyhood ratio test statistic" = "LRT",
             "Uncorrected P-value" = "p_value",
             "FDR-adjusted P-value" = "p_value_cor",
             "Uncorrected P-value significant (P-value < 0.05)" = "significant",
             "FDR-adjusted P-value significant (P-value < 0.05)" = "significant_cor") %>%
      # remove non-significant branches of input$aBSRELTableSignBranches is TRUE
      {if (input$aBSRELTableSignBranches) 
        dplyr::filter(., `FDR-adjusted P-value significant (P-value < 0.05)` == TRUE) else .} %>%
      # change assembly names to species names
      dplyr::mutate(`Branch name` = map_chr(`Branch name`, ~ {
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
    
    # Stop if no gene selected
    req(input$aBSRELGeneInput)
    
    get_aBSREL_results()
    
  }, escape = FALSE, rownames = FALSE, 
  options = list(scrollX = TRUE, 
                 lengthMenu = c(get_aBSREL_results_length(), 10, 20))
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
      paragraph <- h5("Select a gene to display its results")
      paragraph
    } else {
      # If a gene has been selected
      paragraph <- paste0(h5(paste("MEME results of", get_genename_MEME())))
      HTML(paragraph)
    }
  })
  
  # extract the MEME_data row from the selected gene
  get_MEME_data_row <- reactive({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    # load the required RData object
    load(file.path("data", "MEME_data",
                   paste0("MEME_data_", get_genename_MEME(), ".RData")))
    return(RData_row)
  })
  
  get_gene_site_translations <- reactive({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    # Prepare translation table for gene
    gene_trans_df <- MEME_trans_df %>%
      dplyr::filter(genename == get_genename_MEME()) %>%
      dplyr::select(-transcript_id, -genename)
    gene_trans_df
  })
  
  get_MEME_results <- reactive({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)

    MEME_row <- get_MEME_data_row()
  
    MEME_results <- MEME_row %>%
      use_series(csv) %>%
      extract2(1) %>%
      # add an additional site column
      add_column(TOGA_Human_site = seq(nrow(.))) %>%
      # Join with translation df
      dplyr::right_join(get_gene_site_translations(), by = join_by(TOGA_Human_site)) %>%
      dplyr::relocate(MSA_site, TOGA_Human_site, TOGA_Acomys_site, UP_Human_site)
    
    # Needed to retrieve all non TOGA_Human_sites again
    MEME_results <- get_gene_site_translations() %>%
      dplyr::select(MSA_site) %>%
      dplyr::left_join(MEME_results, by = join_by(MSA_site)) %>%
      # remove non-significant sites if input$MEMETableSignSites is TRUE
      {if (input$MEMETableSignSites) dplyr::filter(., `p-value` <= input$p_val_select) else .}
    
    # Combine with substitutions
    MEME_results_with_subs <- get_MEME_acomys_sub() %>%
      dplyr::mutate(branch_assembly = replace(branch_assembly, label == "REFERENCE", "REFERENCE"),
                    Species_Tree = replace(Species_Tree, label == "REFERENCE", "Homo_sapiens")) %>%
      dplyr::filter(Species_Tree %in% c("Homo_sapiens", "Acomys_cahirinus")) %>%
      # change substitution site if A. cahirinus
      dplyr::rowwise() %>%
      dplyr::mutate(subs = if_else(Species_Tree == "Acomys_cahirinus",
                                   translate_branch_subs(TOGA_Acomys_site, subs),
                                   subs)) %>%
      ungroup() %>%
      dplyr::select("MSA_site", "TOGA_Human_site", "TOGA_Acomys_site", "UP_Human_site", "subs", "Species_Tree") %>%
      tidyr::pivot_wider(values_from = "subs", names_from = "Species_Tree") %>%
      dplyr::right_join(MEME_results, by = join_by(MSA_site, TOGA_Human_site, TOGA_Acomys_site, UP_Human_site)) %>%
      dplyr::relocate(Homo_sapiens, Acomys_cahirinus, .after = last_col()) %>%
      dplyr::rename("Substitution in Homo sapiens branch" = Homo_sapiens,
                    "Substitution in A. cahirinus branch" = Acomys_cahirinus) %>%
      dplyr::arrange(MSA_site)
    
    # Filter for substitutions if selected
    if (input$MEMEFilterSitesForSubs == 2) {
      MEME_results_with_subs <- MEME_results_with_subs %>%
        dplyr::filter(!is.na(`Substitution in Homo sapiens branch`))
    }
    if (input$MEMEFilterSitesForSubs == 3) {
      MEME_results_with_subs <- MEME_results_with_subs %>%
        dplyr::filter(!is.na(`Substitution in A. cahirinus branch`))
    }
    MEME_results_with_subs
  })
  
  get_MEME_results_length <- reactive({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    nrow(get_MEME_results())
  })
  
  output$MEME_table <- gt::render_gt({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    get_MEME_results() %>%
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
                    `# branches under selection`, `Total branch length`,
                    `MEME LogL`, `FEL LogL`, `Variation p`)) %>%
      tab_spanner(
        label = html("<strong><em>Substitutions</em></strong>"),
        columns = c(`Substitution in Homo sapiens branch`,
                    `Substitution in A. cahirinus branch`)) %>%
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
                      page_size_values = c(10, 25, 50, 100, get_MEME_results_length()))
    
  })
  
  output$MEMEResultDescriptionTitle <- renderUI({
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      # If no gene has been selected
      paragraph <- h5("Select a gene to display its results")
      paragraph
    } else {
      HTML(paste0(h5(HTML(paste(get_genename_MEME(), "result summary:")))))
    }
  })
  
  # render block of summarizing text 
  output$MEMEResultDescription <- renderUI({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    # Get selected coordinate frame
    coord_sel <- coord_options[[input$MEMEResultCoordSelector]]
    
    # filter with p-value <= 0.1 if input$MEMETableSignSites is FALSE
    filtered_results <- get_MEME_results() %>%
      {if (input$MEMETableSignSites == FALSE) dplyr::filter(., `p-value` <= input$p_val_select) else .}
    
    # get number of significant sites
    n_of_sites <- filtered_results %>%  
      dplyr::pull(coord_sel[1]) %>%
      length()
    
    # get significant site location
    sites <- filtered_results %>%
      dplyr::filter(!is.na(.data[[coord_sel[1]]])) %>%
      dplyr::pull(coord_sel[1]) %>%
      as.character()
    
    # Check if any sites are remaining
    if (length(sites) == 0) {
      sites <- "No sites identified."
    }
    
    sites <- sites %>%
      paste(collapse = ", ")
    
    paragraph <- paste0(p(HTML(paste(get_genename_MEME(), "was found to have", strong(n_of_sites),
                                     strong("codon sites"), "with significant signs (p <= ", input$p_val_select, ") of episodic diversifying",
                                     "selection"))),
                        p(HTML(paste("Relative to the", coord_sel[2], "these codon sites are:", br(), br(), em(sites))))
    )
    HTML(paragraph)
  })
  
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
  
  # # Rendering table for uniprot-mapped sites
  # output$MEME_table_unimapped <- gt::render_gt({
  #   
  #   if(nrow(get_unimapped_MEME_results()) == 0) {
  #     return(NULL)
  #   }
  #   
  #   get_unimapped_MEME_results() %>%
  #     dplyr::rename("Amino acid difference between Human and A. cahirinus sequence" = "branch_subs",
  #            "Amino acid substitution that took place in the A. cahirinus lineage (i.e. branch)" = "branch_subs_ancestral") %>%
  #   gt() %>%
  #     tab_spanner(
  #       label = html("<strong><em>Position</em></strong>"),
  #       columns = c("MSA_site", "TOGA_Human_site", "TOGA_Acomys_site", "UP_Human_site")
  #     ) %>%
  #     tab_spanner(
  #       label = html("<strong><em>MEME maximum likelyhood estimation output</em></strong>"),
  #       columns = c(`&alpha;`, `&beta;<sup>-</sup>`, `p<sup>-</sup>`,
  #                   `&beta;<sup>+</sup>`, `p<sup>+</sup>`,
  #                   LRT, `p-value`,
  #                   `#_branches_under_selection`, `Total_branch_length`,
  #                   `MEME_LogL`, `FEL_LogL`, `Variation_p`)) %>%
  #     cols_align(align = "left") %>%
  #     opt_stylize(style = 1) %>%
  #     cols_label(
  #       MSA_site = "MSA site",
  #       TOGA_Human_site = "TOGA Human site",
  #       TOGA_Acomys_site = "TOGA Acomys site",
  #       UP_Human_site = "UP Human site",
  #       `&alpha;` = gt::html("&alpha;"),
  #       `&beta;<sup>-</sup>` = gt::html("&beta;<sup>-</sup>"),
  #       `p<sup>-</sup>` = gt::html("p<sup>-</sup>"),
  #       `&beta;<sup>+</sup>` = gt::html("&beta;<sup>+</sup>"),
  #       `p<sup>+</sup>` = gt::html("p<sup>+</sup>")
  #     ) %>%
  #     data_color(
  #       columns = LRT,
  #       palette = "Blues") %>%
  #     data_color(
  #       columns = `p-value`,
  #       palette = "Greens",
  #       reverse = TRUE) %>%
  #     opt_interactive(use_compact_mode = TRUE,
  #                     use_resizers = TRUE,
  #                     use_page_size_select = TRUE,
  #                     page_size_values = c(10, 25, 50, get_unimapped_MEME_length()))
  #   
  #   })
  
  ### MEME backend for 'Multiple sequence alignment' tabpanel
  
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
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    # Get selected coordinate frame
    coord_sel <- coord_options[[input$MEMETreeCoordSelector]]
    
    branch_substitutions <- get_MEME_data_row() %>%
      dplyr::pull(branch_subs) %>%
      magrittr::extract2(1)
    
    # Prepare the substitution dataframe for tree merging
    site_branch_subs <- branch_substitutions %>%
      #dplyr::filter(label %in% node_lab) %>%
      tidyr::separate_longer_delim(subs, delim = "/") %>%
      dplyr::mutate(subs = if_else(subs == "", NA_character_, subs)) %>%
      tidyr::separate_wider_regex(subs, c("(?:\\b(?:[ACDEFGHIKLMNPQRSTVWY]|Gap|InMut))?",
                                          "-?",
                                          TOGA_Human_site = "\\d+",
                                          "-?",
                                          "(?:(?:[ACDEFGHIKLMNPQRSTVWY]|Gap|InMut)\\b)?"),
                                  cols_remove = FALSE) %>%
      dplyr::mutate(TOGA_Human_site = as.integer(TOGA_Human_site)) %>%
      # Add other reference frames
      dplyr::left_join(get_gene_site_translations() %>%
                         dplyr::mutate(TOGA_Human_site = as.integer(TOGA_Human_site)) %>%
                         # `.&TOGA_orig_Human_site` are the sites in the
                         # branch_substitutions column
                         dplyr::filter((TOGA_Human_site %in% .$TOGA_Human_site) &
                                         (!is.na(TOGA_Human_site))),
                       by = join_by(TOGA_Human_site)) %>%
      # Add species names
      dplyr::left_join(name_conversion %>%
                         dplyr::select(`Assembly name`, Species_Tree) %>%
                         dplyr::rename(branch_assembly = `Assembly name`),
                       by = join_by(branch_assembly)) %>%
      # Change substitution to match coordinate frame
      dplyr::mutate(subs = purrr::map2_chr(.data[[coord_sel[1]]],
                                           subs,
                                           translate_branch_subs))
    
    site_branch_subs
  })
  
  get_MEME_acomys_sub_pval <- reactive({
    
    # Get selected coordinate frame
    coord_sel <- coord_options[[input$MEMETreeCoordSelector]]
    
    # Link to p-values of sites
    site_branch_subs <- get_MEME_results() %>%
      dplyr::select(MSA_site, `p-value`) %>%
      dplyr::right_join(get_MEME_acomys_sub(), by = join_by(MSA_site)) %>%
      # remove non-significant sites if input$MEMETableSignSites is TRUE
      {if (input$MEMETableSignSites) dplyr::filter(., `p-value` <= input$p_val_select) else .} %>%
      dplyr::mutate(`p-value` = round(`p-value`, digits = 4)) %>%
      # Add homo sapiens
      dplyr::mutate(Species_Tree = replace(Species_Tree, label == "REFERENCE", "Homo_sapiens"))
    
    # Filter for sites with substitutions, if selected
    if (input$MEMEFilterSitesForSubs == 2) {
      sites_to_keep <- site_branch_subs %>%
        dplyr::filter((Species_Tree == "Homo_sapiens") & (!is.na(subs))) %>%
        dplyr::pull(.data[[coord_sel[1]]])
      
      site_branch_subs <- site_branch_subs %>%
        dplyr::filter(.data[[coord_sel[1]]] %in% sites_to_keep)
    }
    if (input$MEMEFilterSitesForSubs == 3) {
      sites_to_keep <- site_branch_subs %>%
        dplyr::filter((Species_Tree == "Acomys_cahirinus") & (!is.na(subs))) %>%
        dplyr::pull(.data[[coord_sel[1]]])
      
      site_branch_subs <- site_branch_subs %>%
        dplyr::filter(.data[[coord_sel[1]]] %in% sites_to_keep)
    }
    site_branch_subs
  })
  
  
  get_site_options <- reactive({
    
    # Stop if no gene selected
    req(input$MEMEGeneInput)
    
    # Get selected coordinate frame
    coord_sel <- coord_options[[input$MEMETreeCoordSelector]]
    
    site_options <- get_MEME_acomys_sub_pval() %>%
      dplyr::pull(.data[[coord_sel[1]]]) %>%
      .[!is.na(.)] %>%
      unique() %>%
      .[order(.)] %>%
      as.character()
    
    site_options
  })
  
  # Get the site selection options
  observe({
  
    freezeReactiveValue(input, "SubSiteSelect")
    updateSelectizeInput(session, "SubSiteSelect", choices = get_site_options(),
                         server = TRUE, 
                         options = list(placeholder = "Site"))
  })
  
  # Get the display table
  get_sub_display_table <- reactive({
    
    req(input$MEMEGeneInput)
    
    # Makes sure nothing gets generated when no sites has been selected
    if (is.null(get_site_options())) {
      return(NULL)
    }
    
    site_display_table <- get_MEME_acomys_sub_pval() %>%
      dplyr::mutate(label = if_else(!is.na(Species_Tree), Species_Tree, label)) %>%
      dplyr::select(label, branch_assembly, branch_projection, MSA_site,
                    TOGA_Human_site, TOGA_Acomys_site, UP_Human_site, subs, `p-value`) %>%
      dplyr::rename(Assembly_name = branch_assembly,
                    TOGA_projection_name = branch_projection,
                    substitution = subs) %>%
      dplyr::arrange(label) %>%
      # Add homo sapiens
      dplyr::mutate(label = replace(label, label == "REFERENCE", "Homo_sapiens"),
                    Assembly_name = replace(Assembly_name, label == "Homo_sapiens", "REFERENCE")) %>%
      # Remove NA (i.e. no substitution) rows
      dplyr::filter(!is.na(substitution))
  })
  
  get_display_table_length <- reactive({
    
    req(input$MEMEGeneInput)
    
    # Makes sure nothing gets generated when no sites has been selected
    if (is.null(get_site_options())) {
      return(NULL)
    }
    
    nrow(get_sub_display_table())
  })
  
  # This prepares the substitution dataframe for use with the phylogenetic tree
  get_MEME_acomys_sub_for_phylo <- reactive({
    
    req(input$MEMEGeneInput)
    
    # Makes sure nothing gets generated when no sites has been selected
    if (length(get_site_options()) == 0) {
      return(NULL)
    }
    
    # Get selected coordinate frame
    coord_sel <- coord_options[[input$MEMETreeCoordSelector]]
    
    cols_to_remove <- c("TOGA_Human_site", "MSA_site", "TOGA_Acomys_site", "UP_Human_site")
    cols_to_remove <- cols_to_remove[cols_to_remove != coord_sel[1]]

    site_branch_subs <- get_MEME_acomys_sub_pval() %>%
      dplyr::select(!all_of(c("p-value", cols_to_remove))) %>%
      dplyr::distinct(label, subs, .keep_all = TRUE) %>%
      tidyr::pivot_wider(names_from = all_of(coord_sel[1]), values_from = subs) %>%
      # Join with tibble_tree
      {if ("branch.length" %in% colnames(.)) {
        dplyr::select(., -node, -branch.length, -branch_assembly, -branch_projection)
      } else {
        dplyr::select(., -node, -branch_assembly, -branch_projection)
      }
      }
    
    site_branch_subs
  })
  
  prepare_phylogeny <- reactive({
    
    req(input$MEMEGeneInput)
    
    # Makes sure nothing gets generated when no sites has been selected
    if (length(get_site_options()) == 0) {
      return(NULL)
    }
    
    # Get and prepare tree from MEME data
    gene_tree <- get_MEME_data_row() %>%
      dplyr::pull(tree)
    
    gene_tree <- ape::read.tree(text = gene_tree)
    tibble_tree <- as_tibble(gene_tree)
    
    # Pattern used to retrieve species names from label
    branch_pattern <- "vs_(.+?)_(.+?_.+?_.+?)$"
    
    # Get branch assembly (e.g. HLacoCah2) and branch projection 
    # (e.g. ENST00000011653_CD4_109) from branch names in separate columns.
    tibble_tree$branch_assembly <- sapply(tibble_tree$label, function(branch) {
      if (branch == "REFERENCE") {
        return("REFERENCE")
      }
      if (stringr::str_detect(branch, pattern = branch_pattern)) {
        return(stringr::str_match(branch, branch_pattern)[2])
      } else {
        return(NA)
      }
    })
    tibble_tree$branch_projection <- sapply(tibble_tree$label, function(branch) {
      if (branch == "REFERENCE") {
        return("REFERENCE")
      }
      if (stringr::str_detect(branch, pattern = branch_pattern)) {
        return(stringr::str_match(branch, branch_pattern)[3])
      } else {
        return(NA)
      }
    })
    
    # Merge tree with substitutions
    tibble_tree_anno <- tibble_tree %>%
      dplyr::left_join(get_MEME_acomys_sub_for_phylo(), by = join_by(label)) %>%
      dplyr::select(-Species_Tree) %>%
      dplyr::left_join(name_conversion %>%
                         dplyr::select(`Assembly name`, Species_Tree) %>%
                         dplyr::rename(branch_assembly = `Assembly name`),
                       by = join_by(branch_assembly)) %>%
      # Add homo sapiens
      dplyr::mutate(Species_Tree = replace(Species_Tree, branch_assembly == "REFERENCE", "Homo_sapiens"))
    
    # prepare order data
    orders <- c("Rodentia", "Lagomorpha", "Primates", "Scandentia", "Artiodactyla",
                "Carnivora", "Pholidota", "Perissodactyla", "Dermoptera", "Chiroptera",
                "Eulipotyphla", "Afrotheria", "Xenarthra")
    
    # filter out species not in the tree and add order information 
    species_orders <- name_conversion %>%
      dplyr::filter(Species_Tree %in% tibble_tree_anno$Species_Tree) %>%
      .[!(duplicated(.$Species_Tree)),] %>%
      dplyr::mutate(species_order = purrr::map_chr(`Taxonomic Lineage`, ~ {
        for (spec_order in orders) {
          if (stringr::str_detect(.x, spec_order)) {
            return(spec_order)
          }
        }    
      })) %>%
      dplyr::select(Species_Tree, species_order) %>%
      #manually add homo sapiens order
      dplyr::add_row(Species_Tree = "Homo_sapiens", species_order = "Primates")
    
    # Annotate tree with order information
    tibble_tree_anno <- tibble_tree_anno %>%
      left_join(species_orders, by = 'Species_Tree')
    
    clades <- lapply(orders, function(i) {
      selected_nodes <- tibble_tree_anno %>%
        filter(species_order == i) %>%
        select(node) %>%
        flatten_int()
      return(selected_nodes)
    })
    names(clades) <- orders
    
    # Remove zero length list elements
    clades <- clades[sapply(clades, length) != 0]
    
    tibble_tree_anno <- ggtree::groupOTU(tibble_tree_anno, clades)
    # remove unnecessary species_order column
    tibble_tree_anno <- tibble_tree_anno %>%
      select(!species_order)
    
    # Find node representing Acomys Cahirinus
    Acomys_node <- tibble_tree_anno %>%
      dplyr::filter(Species_Tree == "Acomys_cahirinus") %>%
      dplyr::pull(node)
    
    # Replace the labels with cleaner species names
    tibble_tree_nodes <- tibble_tree_anno %>%
      dplyr::mutate(label = if_else(!is.na(Species_Tree), Species_Tree, label))
    
    # Detect MRCA of the species orders
    grp_MRCA_nodes <- sapply(unique(tibble_tree_nodes$group), function(grp) {
      # "0" represents the root, which has no MCRA
      if (grp == "0") {
        return(NA)
      }
      grp_species <- tibble_tree_nodes$Species_Tree[tibble_tree_nodes$group == grp]
      ggtree::MRCA(tibble_tree_nodes, grp_species[!is.na(grp_species)])$node
    })
    grp_MRCA_nodes <- grp_MRCA_nodes[!is.na(grp_MRCA_nodes)]
    
    # Create column indicating whether a node is an MRCA node or an Acomys node
    tibble_tree_nodes <- tibble_tree_nodes %>%
      dplyr::mutate(is_MRCA = node %in% grp_MRCA_nodes,
                    is_aco = node == Acomys_node)
    
    Anno_tree <- treeio::as.treedata(tibble_tree_nodes)

    # Plot base tree with tippoints colored by species order/group
    p1 <- ggtree(Anno_tree, aes(x, y), open.angle = 10, size = 0.5) +
      geom_tippoint(mapping = aes(color = group, subset = !is_aco)) +
      #give acomys different label point
      geom_tippoint(mapping = aes(subset = is_aco),
                    shape = 17, size = 2, color = MetBrewer::met.brewer("Cross")[5])
    # Highlight the species orders/groups using geom_highlight
    p2 <- p1 + geom_hilight(mapping = aes(subset = is_MRCA, fill = group), alpha = 0.2) +
      guides(alpha = "none") +
      scale_fill_manual(name = "Order",
                        values = MetBrewer::met.brewer("Cross", type = "continuous", n = 14),
                        breaks = c("Carnivora", "Pholidota", "Perissodactyla",
                                   "Artiodactyla", "Chiroptera", "Eulipotyphla",
                                   "Rodentia", "Lagomorpha", "Primates", "Dermoptera",
                                   "Scandentia", "Afrotheria", "Xenarthra")) +
      scale_color_manual(name = "Order",
                         values = MetBrewer::met.brewer("Cross", type = "continuous", n = 14),
                         breaks = c("Carnivora", "Pholidota", "Perissodactyla",
                                    "Artiodactyla", "Chiroptera", "Eulipotyphla",
                                    "Rodentia", "Lagomorpha", "Primates", "Dermoptera",
                                    "Scandentia", "Afrotheria", "Xenarthra")) +
      # Add tip labels
      geom_tiplab(align = TRUE) +
      # move legend position
      theme(legend.position = "bottom") +
      # Adds 30% more space to the right
      hexpand(0.15, direction = 1)
    
    # Add branch length scale, if the tree contained branch lengths
    if ("branch.length" %in% colnames(tibble_tree_anno)) {
      p2 <- p2 + geom_treescale()
    }
    
    p2
  })
  
  # Add substitution annotation
  annotate_phylo <- reactive({
    
    req(input$SubSiteSelect)
    req(input$MEMEGeneInput)
    
    # Makes sure nothing gets generated when no sites has been selected
    if (length(get_site_options()) == 0) {
      return(NULL)
    }
    
    prepare_phylogeny() +
      geom_label(aes(label = .data[[input$SubSiteSelect]]), hjust = 1.2)
  })
  
  # Display message when no gene is selected
  output$SubPhyloText <- renderUI({
    if (get_genename_MEME() == "") {
      return("Please select a gene to display the phylogenetic tree with the substitution data.")
    } else if (is.null(get_site_options())) {
      return(NULL)
    } else if (length(get_site_options()) == 0) {
      # Get selected coordinate frame
      coord_sel <- coord_options[[input$MEMETreeCoordSelector]]
      return(paste0("For ", get_genename_MEME(), " in the ", coord_sel[1], " reference frame: No substitutions found"))
    } else {
      return(NULL)
    }
  })
  
  output$SubSitePhylo <- renderPlot({
    
    req(input$MEMEGeneInput)
    
    # Makes sure nothing gets generated when no sites has been selected
    if (length(get_site_options()) == 0) {
      return(NULL)
    }
    
    annotate_phylo()
  })
  
  output$MEME_branch_subs <- DT::renderDataTable({
    
    req(input$MEMEGeneInput)
    
    get_sub_display_table()
  }, rownames = FALSE, options = list(scrollX = TRUE,
                                      lengthMenu = c(20, 50, 100, get_display_table_length()))
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
    EBF_data_table %<>% dplyr::filter(branch %in% c("REFERENCE", "HLacoCah2", "HLpsaObe1",
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
      dplyr::filter(EBF >= 100 & sign == TRUE & branch == "HLacoCah2") %>%
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
    EBF_data_table %<>% dplyr::filter(branch %in% c("REFERENCE", "HLacoCah2", "HLpsaObe1",
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
      dplyr::filter(!(is.infinite(EBF)))
    
    EBF_table_MLE <- EBF_table_MLE %>%
      mutate(EBF_100 = if_else(EBF_100 == FALSE, "MEME LRT p-value > 0.05\nor\nEBF < 100",
                               "MEME LRT p-value <= 0.05\nand\nEBF >= 100")) %>%
      dplyr::filter(EBF >= 0)
    
    # Generate data that will be hidden in the plot, for the purpose of keeping
    # the plot range identical when hiding/selecting legend items
    rangex <- range(EBF_table_MLE$site)
    ranges <- expand(tibble(x = rep(rangex, length.out = length(species_order)),
                            y = species_order),
                     crossing(x, y))
    
    # Generate plotly
    plot_ly() %>%
      add_trace(data = dplyr::filter(EBF_table_MLE, EBF_100 == "MEME LRT p-value > 0.05\nor\nEBF < 100"),
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
      add_trace(data = dplyr::filter(EBF_table_MLE, EBF_100 == "MEME LRT p-value <= 0.05\nand\nEBF >= 100"),
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
      {if (input$MEME_EBF_inf_remove) dplyr::filter(., !(is.infinite(EBF))) else .}
  })
  
  MEME_EBF_acomys_selecter <- reactive({
    
    if (is.null(get_EBF_data_table()) | all(get_EBF_data_table() == "no_file")) {
      return(NULL)
    }
    MEME_EBF_inf_remover() %>%
      {if (input$MEME_EBF_acomys_select) dplyr::filter(., branch == "Acomys_cahirinus") else .}
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