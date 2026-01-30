library(aws.s3)

s3BucketName <- "acomys-positive-selection-rshiny"

Sys.setenv("AWS_ACCESS_KEY_ID" = "key",
           "AWS_SECRET_ACCESS_KEY" = "secret_key",
           "AWS_DEFAULT_REGION" = "region")

# Get object keys
object_keys <- aws.s3::get_bucket_df("acomys-positive-selection-rshiny")$Key

# Selected needed object keys
sel_object_keys <- object_keys[stringr::str_detect(object_keys, pattern = paste0("/", "ABI3BP", "/"))]
# Order selected objects
sel_object_keys <- sel_object_keys[
  order(as.integer(stringr::str_split_i(
    stringr::str_replace(sel_object_keys,
                         stringr::fixed(".png"),
                         ""),
    pattern = "_", i = -1)))]

tmp_dir_png <- "MEME_protein_alignments_png_tmp"
tmp_dir_png_gene <- file.path("www", tmp_dir_png, "ABI3BP")
final_dir_png <- "MEME_protein_alignments_png"
final_dir_symlink <- file.path("www", final_dir_png)

for (obj in sel_object_keys) {
  # Get path to temporary download directory
  obj_tmp <- stringr::str_replace(obj, stringr::fixed(final_dir_png), tmp_dir_png)
  
  # Save in temporary folder while retrieving from S3 bucket
  aws.s3::save_object(obj, s3BucketName, file = file.path("www", obj_tmp))
}

# Check the number of retrieved files
n_of_retrieved_objects <- length(list.files(tmp_dir_png_gene))
n_of_objects_to_retrieve <- length(sel_object_keys)

# Check if download succeeded
if (n_of_retrieved_objects == 0) {
  stop("ERROR, no objects found in S3 bucket for the selected gene. This should not happen.")
} else if (n_of_retrieved_objects != n_of_objects_to_retrieve) {
  stop("ERROR in retrieving MSAs from S3 bucket. Fewer files downloaded than expected")
} else {
  # if succeeded, symlink final directory to temporary directory
  from <- normalizePath(tmp_dir_png_gene)
  to <- file.path(normalizePath(final_dir_symlink), "ABI3BP")
  file.symlink(from, to)
}

# For accessing, R shiny already looks in the 'www' dir, so filepaths in
# 'sel_object_keys' are already correct

# get number of files
n_of_gene_MSAs <- length(sel_object_keys)
