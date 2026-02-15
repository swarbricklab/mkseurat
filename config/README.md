# Configuration

Configuring the workflow involves editing config files and preparing input data.

## Config file

To use this workflow as a module in a dataset project, copy the template config file:

```bash
cd top/of/project
mkdir -p config/mkseurat
cp modules/mkseurat/config/template.yaml config/mkseurat/config.yaml
```

Edit the config file to match your dataset structure.

### Required dependencies

| Key | Description |
|-----|-------------|
| `deps.cellranger` | Path to Cell Ranger filtered feature barcode matrices. Each capture should be a subdirectory containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| `deps.captures` | CSV file listing captures to process. Must have a `capture` column |

### Optional dependencies

| Key | Description |
|-----|-------------|
| `deps.demux` | Path to sample assignment files from SNP demux. Each capture directory should contain `cell_assignment.tsv` |
| `deps.annotation` | Path to cell type annotations. Each capture directory should contain `cell_types.csv` |
| `deps.ambient` | Path to ambient RNA profiles. Each capture directory should contain `ambient_summary.csv` |
| `deps.metadata` | Path to experimental metadata CSV. Will be joined to cell metadata |

### Outputs

| Key | Description |
|-----|-------------|
| `outs.results` | Directory for output files (merged Seurat object) |
| `outs.logs` | Directory for log files |

### Parameters

| Key | Description | Default |
|-----|-------------|---------|
| `params.modality` | How to handle multimodal data (`auto`, `Gene Expression`) | `auto` |
| `params.metadata_join_column` | Column to join metadata on | `assignment` |

## Captures CSV

The captures CSV file must have at minimum a `capture` column listing the capture IDs to process.
Capture IDs should match the subdirectory names in `deps.cellranger`.

Example:
```csv
capture
NC001
NC002
NC003
```

## Sample assignment files

If using SNP-based demultiplexing, each capture should have a `cell_assignment.tsv` file with columns:
- `barcode`: Cell barcode
- `status`: Assignment status (singlet, doublet, unassigned)
- `assignment`: Sample assignment

## Metadata file

Optional experimental metadata CSV should have a column matching `params.metadata_join_column` (default: `assignment` or `sample`).
Additional columns will be added to cell metadata in the Seurat object.
