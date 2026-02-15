# workflow/rules/metadata.smk
# Rules for attaching experimental metadata

rule attach_metadata:
    """
    Attach experimental metadata (e.g., sample information, clinical data)
    to the merged Seurat object. Joins metadata on sample/assignment column.
    """
    input:
        merged = out_dir / "merged.qs",
        metadata = config['deps']['metadata']
    output:
        annotated = out_dir / "merged_annotated.qs"
    log:
        log_dir / "attach_metadata.log"
    params:
        join_column = config.get('params', {}).get('metadata_join_column', 'assignment')
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/attach_metadata.R"
