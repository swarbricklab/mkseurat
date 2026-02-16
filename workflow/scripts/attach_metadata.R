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
message("Metadata join column: ", join_column)

# Load merged object
message("Loading merged object...")
merged_object <- qread(merged_file)
message("Loaded object with ", ncol(merged_object), " cells")

# Load metadata
message("Loading metadata...")
meta <- read_csv(metadata_file, show_col_types = FALSE)
message("Loaded metadata with ", nrow(meta), " rows and ", ncol(meta), " columns")

# Rename the metadata join column to 'sample_id' to match the Seurat object
# The Seurat object always uses 'sample_id', but the metadata CSV may use a different column name
if (join_column != "sample_id") {
    if (!join_column %in% colnames(meta)) {
        stop("Join column '", join_column, "' not found in metadata. Available columns: ",
             paste(colnames(meta), collapse = ", "))
    }
    message("Renaming metadata column '", join_column, "' to 'sample_id' for join")
    meta <- meta %>% rename(sample_id = !!sym(join_column))
}

# Remove capture_id column if present (already stored in Seurat object as 'capture')
# and deduplicate to get one row per sample_id
if ("capture_id" %in% colnames(meta)) {
    message("Removing 'capture_id' column (already in Seurat as 'capture')")
    meta <- meta %>% select(-capture_id)
}
n_before <- nrow(meta)
meta <- meta %>% distinct()
message("Deduplicated metadata: ", n_before, " -> ", nrow(meta), " rows")

# Join metadata to cell metadata on 'sample_id'
message("Joining metadata on 'sample_id' column...")
current_meta <- merged_object[[]]
n_cells <- nrow(current_meta)
# Preserve cell barcodes (rownames) through the join
current_meta$cell_barcode <- rownames(current_meta)
new_meta <- left_join(current_meta, meta, by = "sample_id")

# Verify join didn't create duplicates
if (nrow(new_meta) != n_cells) {
    stop("Join created duplicate rows (", nrow(new_meta), " vs ", n_cells, 
         "). Check for duplicate sample_id values in metadata.")
}

rownames(new_meta) <- new_meta$cell_barcode
new_meta$cell_barcode <- NULL

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
