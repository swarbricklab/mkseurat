#! /usr/bin/env Rscript
# workflow/scripts/merge_captures.R
# Merge individual capture Seurat objects into a single object

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(qs)
})

# Set up logging
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

# Get input and output paths
rds_files <- snakemake@input[["rds_files"]]
output_file <- snakemake@output[["merged"]]

message("Merging ", length(rds_files), " capture objects...")
message("Output file: ", output_file)

# Load all objects
objects <- list()
for (i in seq_along(rds_files)) {
    message("Loading: ", rds_files[i])
    objects[[i]] <- readRDS(rds_files[i])
    message("  - Cells: ", ncol(objects[[i]]))
}

# Merge objects
if (length(objects) == 1) {
    message("Only one capture - no merging needed")
    merged_object <- objects[[1]]
} else {
    message("Merging captures...")
    merged_object <- merge(
        x = objects[[1]],
        y = objects[2:length(objects)]
    )
    
    # Join layers for RNA assay
    message("Joining RNA layers...")
    merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])
    
    # Join layers for AB assay if present
    if ("AB" %in% names(merged_object@assays)) {
        message("Joining AB layers...")
        merged_object[["AB"]] <- JoinLayers(merged_object[["AB"]])
    }
}

message("Merged object has ", ncol(merged_object), " cells total")

# Save merged object
message("Saving merged object to: ", output_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
qsave(merged_object, file = output_file)

message("Done merging captures")

sink(type = "message"); sink()
close(log_con)
