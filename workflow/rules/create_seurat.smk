# workflow/rules/create_seurat.smk
# Rules for creating individual Seurat objects per capture

rule create_seurat_object:
    """
    Create a Seurat object from Cell Ranger filtered feature barcode matrix.
    Handles both unimodal (gene expression only) and multimodal (GEX + antibody capture) data.
    Attaches sample assignments, cell annotations, and ambient profiles if available.
    Prefixes barcodes with capture ID to ensure uniqueness when merging.
    """
    input:
        matrix_dir = lambda wc: Path(config['deps']['cellranger']) / wc.capture
    output:
        rds = temp(out_dir / "per_capture/{capture}.rds")
    log:
        log_dir / "create_seurat/{capture}.log"
    params:
        assignment_root = config['deps'].get('demux', ''),
        annotation_root = config['deps'].get('annotation', ''),
        ambient_root = config['deps'].get('ambient', ''),
        modality = config.get('params', {}).get('modality', 'auto')
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/create_seurat.R"
