### define functions needed for creating the annotated alignments using the
### ggmsa package

#' Split a gapped DNA alignment into codons
#'
#' Splits each sequence in a gapped DNA multiple sequence alignment into
#' codons (triplets), preserving gap characters.
#'
#' @param myAln A \code{\link[Biostrings]{DNAStringSet}} containing a gapped
#'   DNA alignment. All sequences are assumed to have equal width and to be
#'   in frame.
#'
#' @return
#' A list with one element per sequence. Each element is a character vector
#' of codons.
#'
#' @seealso
#' \code{\link[Biostrings]{DNAStringSet}},
#' \code{\link[Biostrings]{Views}}
getCodons <- function(myAln) {
  seqs <- as.character(myAln)
  len <- width(myAln)[1]
  starts <- seq(from=1, to=len, by=3)
  ends <- starts + 2
  myViews <- lapply(myAln, function(x) { 
    Views(x, starts, ends)
  })
  myCodons <- lapply(myViews, function(x) {
    as.character(DNAStringSet(x))
  })
  myCodons
}

#' Translate codons into an amino-acid sequence
#'
#' Translates a character vector of codons into an amino-acid sequence using
#' the standard genetic code, with special handling of gap codons and
#' untranslatable codons.
#'
#' @param myCodons Character vector of codons (length-3 strings).
#' @param unknownCodonTranslatesTo Single character used to replace codons
#'   that cannot be translated (e.g. due to frameshifts).
#'   Default is \code{"-"}.
#'
#' @return
#' A single character string representing the translated amino-acid sequence.
#'
#' @seealso
#' \code{\link[Biostrings]{GENETIC_CODE}}
translateCodons <- function(myCodons, unknownCodonTranslatesTo="-") {
  ## make new genetic code
  gapCodon <- "-"
  names(gapCodon) <- "---"
  my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
  
  ## translate the codons
  pep <- my_GENETIC_CODE[myCodons]
  
  ## check for codons that were not possible to translate, e.g. frameshift codons
  if (sum(is.na(pep))>0) {
    warning(paste0("\nwarning - there were codons I could not translate. Using this character ", unknownCodonTranslatesTo, "\n\n"))
    pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
  }
  
  ## prep for output
  pep <- paste(pep, collapse="")
  return(pep)
}

#' Translate a gapped DNA alignment into an amino-acid alignment
#'
#' Wrapper function that splits a gapped DNA alignment into codons and
#' translates each sequence into an amino-acid sequence.
#'
#' @param myAln A \code{\link[Biostrings]{DNAStringSet}} containing a gapped
#'   DNA alignment.
#' @param unknownCodonTranslatesTo Single character used to replace codons
#'   that cannot be translated. Default is \code{"-"}.
#'
#' @return
#' An \code{\link[Biostrings]{AAStringSet}} containing the translated
#' amino-acid alignment.
#'
#' @seealso
#' \code{\link{getCodons}},
#' \code{\link{translateCodons}},
#' \code{\link[Biostrings]{AAStringSet}}
translateGappedAln <- function(myAln, unknownCodonTranslatesTo="-") {
  myCodons <- getCodons(myAln)
  myAAaln <- AAStringSet(unlist(lapply(myCodons, translateCodons, unknownCodonTranslatesTo=unknownCodonTranslatesTo)))
  return(myAAaln)
}

#' Build mapping from ungapped sequence positions to gapped reference (MSA)
#'
#' @param ref_gapped Character string of the reference sequence with gaps ("-")
#'
#' @return data.frame with ungapped_position, ref_position (MSA column)
BuildMSAMapping <- function(ref_gapped) {
  
  ref_gapped <- as.character(ref_gapped)
  
  ref_chars <- strsplit(ref_gapped, "")[[1]]
  
  ungapped_pos <- 0
  
  mapping <- data.frame(
    ungapped_position = integer(0),
    ref_position = integer(0)
  )
  
  for (i in seq_along(ref_chars)) {
    if (ref_chars[i] != "-") {
      ungapped_pos <- ungapped_pos + 1
      
      mapping <- rbind(
        mapping,
        data.frame(
          ungapped_position = ungapped_pos,
          ref_position = i
        )
      )
    }
  }
  
  return(mapping)
}

bar_data <- function(tidy){
  character_position <- unique(tidy$position)
  conservation_score <- lapply(character_position, function(j) {
    cloumn_data <- tidy[tidy$position == j, ]
    character_frequency <- table(cloumn_data$character) %>% as.data.frame
    max_frequency <- character_frequency[character_frequency[2] ==
                                           max(character_frequency[2]),]
    max_frequency$Var1 <- as.character(max_frequency$Var1)
    if(nrow(max_frequency) == 1) {
      max_frequency <- max_frequency[1,]
    }else {
      max_frequency <- max_frequency[1,]
    }
  }) %>% do.call("rbind", .)
  conservation_score["pos"] <- character_position
  return(conservation_score)
}

plot_annotated_msa_EBF <- function(seq_aa,
                               start,
                               end,
                               plot_end,
                               anno_df,
                               seq_length_df = NULL,
                               colors_sign = c("FALSE" = "#B3B3B3",
                                               "TRUE" = "#FF7F0E"),
                               x_expand = 1.1,
                               x_breaks_interval = 5,
                               rel_heights = c(10, 1, 1, 1),
                               rel_widths = c(10, 1),
                               EBF_col_range = c(0, 1000)) {
  
  # Subset p-value annotation for the selected range
  anno_subset <- anno_df[(anno_df$pos >= start) & (anno_df$pos <= end),]
  
  p <- ggmsa(seq_aa, start = start, end = end, seq_name = TRUE) +
    theme(axis.text.y = element_text(size = 4),
          axis.text.x = element_text(size = 6),
          plot.margin = unit(c(5, 0, 0.5, 0), "pt")) +
    scale_x_continuous(expand = expansion(add = x_expand),
                       breaks = scales::breaks_width(width = x_breaks_interval)) +
    coord_cartesian(xlim = c(start, plot_end)) +
    geom_seqlogo()
  
  # Plot the sequence lengths if given
  if (!is.null(seq_length_df)) {
    p <- p + geom_text(data = seq_length_df, aes(x = x, y = y, label = seq_length),
                       hjust = 0,
                       vjust = 0.5,
                       size = 2)
  }
  
  pval_barplot <- ggplot(anno_subset, aes(x = pos, y = 1,
                                          fill = `-log10(p-value)`)) +
    geom_tile(width = 1, height = 1) +
    scale_x_continuous(expand = expansion(add = x_expand)) +
    coord_cartesian(xlim = c(start, plot_end)) +
    scale_fill_viridis_c(option = "viridis", direction = 1,
                         limits = c(0, 3),
                         oob = scales::squish,
                         name = "MEME site P-value\n[-log10(P-value)]",
                         na.value = "white") +
    ylab("MEME site P-value") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_text(size = 5, angle = 0,
                                      hjust = 1, vjust = 0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(0.3, "cm"),      # shrink key height
          legend.key.width  = unit(0.3, "cm"),      # shrink key width
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "pt"))
  
  
  legend <- cowplot::get_legend(
    pval_barplot +
      theme(
        legend.position = "left",
        legend.justification = "left",
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0)
      )
  )
  
  pval_barplot_nolegend <- pval_barplot +
    theme(legend.position = "none")
  
  sign_cb <- ggplot(anno_subset, aes(x = pos, y = 1, fill = as.factor(sign))) +
    geom_tile(width = 1, height = 1) +
    scale_x_continuous(expand = expansion(add = x_expand)) +
    coord_cartesian(xlim = c(start, plot_end)) +
    scale_fill_manual(values = colors_sign, name = "MEME site\nsignificance",
                      na.value = "white",
                      na.translate = FALSE) +
    ylab("MEME site Significance") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_text(size = 5, angle = 0,
                                      hjust = 1, vjust = 0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(0.3, "cm"),      # shrink key height
          legend.key.width  = unit(0.3, "cm"),      # shrink key width
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "pt"))
  
  legend_sign <- cowplot::get_legend(
    sign_cb +
      theme(
        legend.position = "left",
        legend.justification = "left",
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0)
      )
  )
  
  sign_cb_nolegend <- sign_cb +
    theme(legend.position = "none")
  
  ### new
  
  EBF_cb <- ggplot(anno_subset, aes(x = pos, y = 1, fill = EBF)) +
    geom_tile(width = 1, height = 1) +
    scale_x_continuous(expand = expansion(add = x_expand)) +
    coord_cartesian(xlim = c(start, plot_end)) +
    scale_fill_viridis_c(option = "magma", direction = 1,
                         limits = EBF_col_range,
                         oob = scales::squish,
                         name = "A. cahirinus EBF",
                         na.value = "white") +
    ylab("A. cahirinus EBF") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_text(size = 5, angle = 0,
                                      hjust = 1, vjust = 0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(0.3, "cm"),      # shrink key height
          legend.key.width  = unit(0.3, "cm"),      # shrink key width
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "pt"))
  
  legend_EBF <- cowplot::get_legend(
    EBF_cb +
      theme(
        legend.position = "left",
        legend.justification = "left",
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0)
      )
  )
  
  EBF_cb_nolegend <- EBF_cb +
    theme(legend.position = "none")
  
  ### end of new
  
  p_legends <- plot_grid(legend, legend_sign, legend_EBF,
                         ncol = 1)
  
  anno_p <- plot_grid(p,
                      pval_barplot_nolegend,
                      sign_cb_nolegend,
                      EBF_cb_nolegend,
                      ncol = 1,
                      rel_heights = rel_heights,
                      align = "vh",
                      axis = "lt")
  
  anno_p <- plot_grid(anno_p, p_legends,
                      ncol = 2,
                      rel_widths = rel_widths,
                      align = "vh",
                      axis = "lt") +
    # Add margins to final cowplot
    theme(plot.margin = unit(c(5, 5, 5, 5), "pt"))  # top, right, bottom, left
  
  return(anno_p)
}

build_seq_lengths_df <- function(seq, pos_end, gap_symbol = "-") {
  # Get ungapped lengths of each individual sequence until `pos_end`
  seq_lengths <- sapply(names(seq), function(id) {
    id_seq <- AAString(as.character(seq[id]),
                       start = 1, nchar = pos_end)
    id_seq_no_gap <- replaceAt(id_seq,at = matchPattern(gap_symbol, id_seq),
                               value = "")
    return(length(id_seq_no_gap))
  })
  
  seq_length_df <- data.frame(
    seq_name = names(seq_lengths),
    seq_length = seq_lengths,
    # position on x-axis: end of the alignment for each sequence
    x = pos_end + 1,       # slightly past the right edge
    y = rev(seq_along(seq_lengths))  # vertical position of each sequence in ggmsa
  )
  return(seq_length_df)
}
