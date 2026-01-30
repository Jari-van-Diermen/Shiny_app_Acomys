
### Load the MEME site translation table
MEME_trans_df <- read_delim("data/translation_table/MEME_gene_site_translation_all_positions.tsv",
                            delim = "\t",
                            escape_double = FALSE,
                            trim_ws = TRUE)

# load the required RData object
load(file.path("data", "MEME_data",
               paste0("MEME_data_", "LRRIQ1", ".RData")))

MEME_row <- RData_row
rm(RData_row)

MEME_row$tree

name_conversion$Species_Tree

gene_assembly_names <- MEME_row %>%
  use_series(branch_subs) %>%
  extract2(1) %>%
  dplyr::filter(!is.na(branch_assembly)) %>%
  dplyr::pull(branch_assembly) %>%
  unique()

gene_tree <- MEME_row %>%
  dplyr::pull(tree)

library(ape)
library(ggtree)

gene_tree <- ape::read.tree(text = gene_tree)
tibble_tree <- as_tibble(gene_tree)

branch_pattern <- "vs_(.+?)_(.+?_.+?_.+?)$"
species_branch <- "HLacoCah2"

# Get branch assembly (e.g. HLacoCah2) and branch projection 
# (e.g. ENST00000011653_CD4_109) from branch names in separate columns.
tibble_tree$branch_assembly <- sapply(tibble_tree$label, function(branch) {
  if (branch == "REFERENCE") {
    return("REFERENCE")
  }
  if (str_detect(branch, pattern = branch_pattern)) {
    return(str_match(branch, branch_pattern)[2])
  } else {
    return(NA)
  }
})
tibble_tree$branch_projection <- sapply(tibble_tree$label, function(branch) {
  if (branch == "REFERENCE") {
    return("REFERENCE")
  }
  if (str_detect(branch, pattern = branch_pattern)) {
    return(str_match(branch, branch_pattern)[3])
  } else {
    return(NA)
  }
})

# # get node ID of species_branch and root
# species_branch_node <- tibble_tree$node[which(tibble_tree$branch_assembly == species_branch)]
# root_node <- tibble_tree$node[which(tibble_tree$label == "")]
# # get nodepath from species_branch to root
# nodes <- ape::nodepath(gene_tree, from = root_node, to = species_branch_node)
# 
# # retrieve node labels
# node_lab <- ggtree::nodelab(gene_tree, nodes)
# node_lab <- replace(node_lab, node_lab == "", "root")
# 
# # get significant sites
# sign_sites <- MEME_row %>%
#   magrittr::use_series("sites") %>%
#   magrittr::extract2(1) %>%
#   as.character()

branch_substitutions <- MEME_row %>%
  dplyr::pull(branch_subs) %>%
  magrittr::extract2(1)

# Prepare the substitution dataframe for tree merging
site_branch_subs <- branch_substitutions %>%
  #dplyr::filter(label %in% node_lab) %>%
  tidyr::separate_longer_delim(subs, delim = "/") %>%
  dplyr::mutate(subs = if_else(subs == "", NA_character_, subs)) %>%
  tidyr::separate_wider_regex(subs, c("\\b[ACDEFGHIKLMNPQRSTVWY]-?",
                                      TOGA_orig_Human_site = "\\d+",
                                      "-?(?:[ACDEFGHIKLMNPQRSTVWY]|Gap)\\b"),
                              cols_remove = FALSE) %>%
  dplyr::mutate(TOGA_orig_Human_site = as.integer(TOGA_orig_Human_site)) %>%
  dplyr::left_join(MEME_trans_df %>%
                     dplyr::mutate(TOGA_orig_Human_site = as.integer(TOGA_orig_Human_site)) %>%
                     # `.&TOGA_orig_Human_site` are the sites in the
                     # branch_substitutions column
                     dplyr::filter((genename == "CD4") &
                                   (TOGA_orig_Human_site %in% .$TOGA_orig_Human_site) &
                                   (!is.na(TOGA_orig_Human_site))) %>%
                     dplyr::select(-transcript_id, -genename),
                   by = join_by(TOGA_orig_Human_site)) %>%
  # Rename site coordinate columns
  dplyr::rename(MSA_site = ggmsa_site,
                TOGA_Human_site = TOGA_orig_Human_site,
                TOGA_Acomys_site = TOGA_orig_Acomys_site) %>% 
  # Add species names
  dplyr::left_join(name_conversion %>%
                     dplyr::select(`Assembly name`, Species_Tree) %>%
                     dplyr::rename(branch_assembly = `Assembly name`),
                   by = join_by(branch_assembly)) %>%
  # Change substitution to match coordinate frame
  dplyr::mutate(subs = purrr::map2_chr(.data[["MSA_site"]],
                                       subs,
                                       translate_branch_subs))

  site_branch_subs %>%
    dplyr::pull(.data[["MSA_site"]]) %>%
    .[!is.na(.)] %>%
    unique() %>%
    .[order(.)] %>%
    as.character()


  # site_branch_subs can now be used as a table for display
  site_display_table <- site_branch_subs %>%
    dplyr::mutate(label = if_else(!is.na(Species_Tree), Species_Tree, label)) %>%
    dplyr::select(label, branch_assembly, branch_projection, MSA_site,
                  TOGA_Human_site, TOGA_Acomys_site, UP_Human_site, subs) %>%
    dplyr::rename(Assembly_name = branch_assembly,
                  TOGA_projection_name = branch_projection,
                  substitution = subs) %>%
    dplyr::arrange(label) %>%
    # Add homo sapiens
    dplyr::mutate(label = replace(label, label == "REFERENCE", "Homo_sapiens"),
                  Assembly_name = replace(Assembly_name, label == "Homo_sapiens", "REFERENCE")) %>%
    # Remove NA (i.e. no substitution) rows
    dplyr::filter(!is.na(substitution))

  # Then we can use the site_branch_subs for the tree tibble
  sel_ref_frame <- "TOGA_Human_site"
  cols_to_remove <- c("TOGA_Human_site", "MSA_site", "TOGA_Acomys_site", "UP_Human_site")
  cols_to_remove <- cols_to_remove[cols_to_remove != sel_ref_frame]
  site_branch_subs <- site_branch_subs %>%
    dplyr::select(!all_of(cols_to_remove)) %>%
    tidyr::pivot_wider(names_from = all_of(sel_ref_frame), values_from = subs) %>%
    # Join with tibble_tree
    dplyr::select(-node, -branch.length, -branch_assembly, -branch_projection, -`NA`)

# Merge tree with substitutions
  tibble_tree_anno <- tibble_tree %>%
    dplyr::left_join(site_branch_subs, by = join_by(label)) %>%
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
  if (grp == "0") {
    return(NA)
  }
  grp_species <- tibble_tree_nodes$Species_Tree[tibble_tree_nodes$group == grp]
  ggtree::MRCA(tibble_tree_nodes, grp_species[!is.na(grp_species)])$node
})
grp_MRCA_nodes <- grp_MRCA_nodes[!is.na(grp_MRCA_nodes)]

Anno_tree <- treeio::as.treedata(tibble_tree_nodes)

# Plot base tree with tippoints colored by species order/group
p1 <- ggtree(Anno_tree, aes(x, y), open.angle = 10, size = 0.5) +
  geom_tippoint(mapping = aes(color = group, subset = !(node %in% c(Acomys_node)))) +
  #give acomys different label point
  geom_tippoint(mapping = aes(subset = node %in% c(Acomys_node)),
                shape = 17, size = 2, color = MetBrewer::met.brewer("Cross")[5])
# Highlight the species orders/groups using geom_highlight
p2 <- p1 + geom_hilight(mapping = aes(subset = node %in% grp_MRCA_nodes,
                                      fill = group), alpha = 0.2) +
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
  # Add branch length scale
  geom_treescale() +
  # move legend position
  theme(legend.position = "bottom")
  
# Add substitution annotation
p3 <- p2 + geom_label(aes(label = `43`), hjust = 1.2)
p3
  

validate()
  
MEME_row$tree
  
  