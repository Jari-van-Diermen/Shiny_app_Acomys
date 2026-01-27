#' Split an RData object into multiple files by row
#'
#' Takes a complete MEME or aBSREL-like data object and creates a new RData file
#' for each row, based on a specified column.
#'
#' @param RData_object An RData object to split into multiple files using its rows.
#' @param column_name Character string giving the name of the column used for
#'   splitting the object.
#' @param RData_basename Character string giving the basename for the generated
#'   RData files.
#' @param subdir_name Character string giving the name of the subdirectory where
#'   the generated RData files will be stored.
#'
#' @return Invisibly returns a character vector with the paths to the generated
#'   RData files.
#'
#' @export
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
