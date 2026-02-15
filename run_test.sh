#!/bin/bash
# Run mkseurat workflow in standalone mode (for testing)

set -e
eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

workflow_profile="--workflow-profile profiles/workflow"

snakemake $workflow_profile \
    --snakefile workflow/Snakefile \
    --configfile config/test.yaml \
    "$@"
