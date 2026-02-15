#! /usr/bin/env Rscript
# workflow/scripts/attach_metadata.R
# Attach experimental metadata to merged Seurat object

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(readr)
    library(dplyr)
    library(qs)
})

# Set up logging
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

# Get input and output paths
merged_file <- snakemake@input[["merged"]]
metadata_file <- snakemake@input[["metadata"]]
output_file <- snakemake@output[["annotated"]]
join_column <- snakemake@params[["join_column"]]

message("Attaching metadata to merged Seurat object")
message("Input: ", merged_file)
message("Metadata: ", metadata_file)
message("Join column: ", join_column)

# Load merged object
message("Loading merged object...")
merged_object <- qread(merged_file)
message("Loaded object with ", ncol(merged_object), " cells")

# Load metadata
message("Loading metadata...")
meta <- read_csv(metadata_file, show_col_types = FALSE)
message("Loaded metadata with ", nrow(meta), " rows and ", ncol(meta), " columns")

# Rename the join column in metadata to match the Seurat object
# Typically metadata has 'sample' column and Seurat object has 'assignment'
if ("sample" %in% colnames(meta) && join_column == "assignment") {
    meta <- meta %>% rename(assignment = sample)
}

# Join metadata to cell metadata
message("Joining metadata on '", join_column, "' column...")
current_meta <- merged_object[[]]
new_meta <- left_join(current_meta, meta, by = join_column)

# Check for successful join
n_matched <- sum(!is.na(new_meta[[colnames(meta)[2]]]))
message("Matched metadata for ", n_matched, " / ", nrow(new_meta), " cells")

# Update metadata in Seurat object
merged_object[[]] <- new_meta

# Save annotated object
message("Saving annotated object to: ", output_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
qsave(merged_object, file = output_file)

message("Done attaching metadata")

sink(type = "message"); sink()
close(log_con)
