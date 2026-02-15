# Installation

## Prerequisites

- Conda or Mamba
- Git (for submodule installation)
- DVC (for data tracking, optional)

## Snakemake Environment

Create the Snakemake environment:

```bash
conda env create -f env/snakemake_7.32.4.yaml
```

Activate the environment before running the workflow:

```bash
conda activate snakemake_7.32.4
```

## Installing as a Submodule

To use this workflow in a dataset project:

```bash
cd /path/to/dataset
git submodule add <repo-url> modules/mkseurat
git submodule update --init --recursive
```

## R Packages

The required R packages are installed automatically via conda when running rules.
The environment is defined in `workflow/envs/seurat.yaml`.

Key packages:
- Seurat 5.1
- SeuratObject
- qs (for fast serialization)
- tidyverse (readr, dplyr, tibble)
