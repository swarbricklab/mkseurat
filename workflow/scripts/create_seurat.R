#! /usr/bin/env Rscript
# workflow/scripts/create_seurat.R
# Create a Seurat object from Cell Ranger output for a single capture

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(readr)
    library(dplyr)
    library(tibble)
})

# Set up logging
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

# Extract parameters
matrix_dir <- snakemake@input[["matrix_dir"]]
output_rds <- snakemake@output[["rds"]]
assignment_root <- snakemake@params[["assignment_root"]]
annotation_root <- snakemake@params[["annotation_root"]]
ambient_root <- snakemake@params[["ambient_root"]]
modality <- snakemake@params[["modality"]]

# Extract capture name from wildcard
capture <- snakemake@wildcards[["capture"]]

message("Creating Seurat object for capture: ", capture)
message("Matrix directory: ", matrix_dir)

# Helper function to load sample assignments
load_assignments <- function(capture, barcodes, assignment_root) {
    if (is.null(assignment_root) || assignment_root == "") {
        message("No assignment root configured. Using capture as sample ID.")
        return(data.frame(
            status = rep("singlet", length(barcodes)),
            assignment = rep(capture, length(barcodes)),
            row.names = barcodes
        ))
    }
    
    assignment_path <- file.path(assignment_root, capture, "cell_assignment.tsv")
    message("Looking for assignments at: ", assignment_path)
    
    if (!file.exists(assignment_path)) {
        message("No sample assignments found. Using capture as sample ID.")
        return(data.frame(
            status = rep("singlet", length(barcodes)),
            assignment = rep(capture, length(barcodes)),
            row.names = barcodes
        ))
    }
    
    assignments <- read_tsv(assignment_path, show_col_types = FALSE) %>%
        column_to_rownames("barcode")
    message("Loaded assignments for ", nrow(assignments), " cells")
    return(assignments)
}

# Helper function to load cell annotations
load_annotations <- function(capture, annotation_root) {
    if (is.null(annotation_root) || annotation_root == "") {
        message("No annotation root configured.")
        return(data.frame())
    }
    
    annotation_path <- file.path(annotation_root, capture, "cell_types.csv")
    message("Looking for annotations at: ", annotation_path)
    
    if (!file.exists(annotation_path)) {
        message("No cell annotations found.")
        return(data.frame())
    }
    
    annotations <- read_csv(annotation_path, show_col_types = FALSE) %>%
        column_to_rownames("barcode")
    message("Loaded annotations for ", nrow(annotations), " cells")
    return(annotations)
}

# Helper function to load ambient profiles
load_ambient <- function(capture, ambient_root) {
    if (is.null(ambient_root) || ambient_root == "") {
        message("No ambient root configured.")
        return(data.frame())
    }
    
    ambient_path <- file.path(ambient_root, capture, "ambient_summary.csv")
    message("Looking for ambient profiles at: ", ambient_path)
    
    if (!file.exists(ambient_path)) {
        message("No ambient profiles found.")
        return(data.frame())
    }
    
    ambient <- read_csv(ambient_path, show_col_types = FALSE) %>%
        column_to_rownames("barcode")
    message("Loaded ambient profiles for ", nrow(ambient), " cells")
    return(ambient)
}

# Read Cell Ranger output
message("Reading Cell Ranger output...")
counts <- Read10X(data.dir = matrix_dir)

# Handle multimodal data
if (is.list(counts)) {
    message("Detected multimodal data")
    available_assays <- names(counts)
    message("Available assays: ", paste(available_assays, collapse = ", "))
    
    # Determine which modality to use for main object
    if (modality == "auto" || modality == "Gene Expression") {
        if ("Gene Expression" %in% available_assays) {
            gex_counts <- counts[["Gene Expression"]]
        } else {
            stop("Gene Expression assay not found in multimodal data")
        }
    } else {
        gex_counts <- counts[[modality]]
    }
    
    # Create main object with gene expression
    message("Creating Seurat object with Gene Expression data...")
    obj <- CreateSeuratObject(counts = gex_counts)
    
    # Add antibody capture assay if present
    if ("Antibody Capture" %in% available_assays) {
        message("Adding Antibody Capture assay...")
        ab_counts <- counts[["Antibody Capture"]]
        # Filter to cells present in GEX
        common_cells <- intersect(colnames(obj), colnames(ab_counts))
        ab_assay <- CreateAssay5Object(counts = ab_counts[, common_cells])
        obj[["AB"]] <- ab_assay
        message("Added AB assay with ", nrow(ab_counts), " features")
    }
} else {
    # Unimodal data (gene expression only)
    message("Detected unimodal data (gene expression only)")
    obj <- CreateSeuratObject(counts = counts)
}

message("Created Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features")

# Get barcodes for metadata attachment
barcodes <- colnames(obj)

# Load and attach metadata
obj <- obj %>%
    AddMetaData(metadata = load_assignments(capture, barcodes, assignment_root)) %>%
    AddMetaData(metadata = load_annotations(capture, annotation_root)) %>%
    AddMetaData(metadata = load_ambient(capture, ambient_root))

# Add capture ID to cell barcodes (for uniqueness when merging)
prefixed_barcodes <- paste0(capture, "_", colnames(obj))
colnames(obj) <- prefixed_barcodes

# Add capture metadata
obj$capture <- capture

message("Saving Seurat object to: ", output_rds)
saveRDS(obj, file = output_rds)
message("Done creating Seurat object for capture: ", capture)

sink(type = "message"); sink()
close(log_con)
