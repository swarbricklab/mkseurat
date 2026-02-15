# workflow/rules/merge.smk
# Rules for merging per-capture Seurat objects

rule merge_captures:
    """
    Merge all per-capture Seurat objects into a single merged object.
    Uses Seurat's merge function and joins layers for proper integration.
    """
    input:
        rds_files = expand(out_dir / "per_capture/{capture}.rds", capture=captures)
    output:
        merged = out_dir / "merged.qs"
    log:
        log_dir / "merge_captures.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/merge_captures.R"
