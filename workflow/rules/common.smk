# workflow/rules/common.smk
# Common functions and utilities for the mkseurat workflow

def get_cellranger_path(wildcards):
    """Return the path to Cell Ranger output for a capture."""
    return Path(config['deps']['cellranger']) / wildcards.capture

def get_assignment_path(wildcards):
    """Return the path to sample assignment file for a capture."""
    if 'demux' in config['deps']:
        return Path(config['deps']['demux']) / wildcards.capture / "cell_assignment.tsv"
    return None

def get_annotation_path(wildcards):
    """Return the path to cell annotation file for a capture."""
    if 'annotation' in config['deps']:
        return Path(config['deps']['annotation']) / wildcards.capture / "cell_types.csv"
    return None

def get_ambient_path(wildcards):
    """Return the path to ambient profile file for a capture."""
    if 'ambient' in config['deps']:
        return Path(config['deps']['ambient']) / wildcards.capture / "ambient_summary.csv"
    return None

def has_optional_input(config_key):
    """Check if an optional input is configured."""
    return config_key in config.get('deps', {}) and config['deps'][config_key]
