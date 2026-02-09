"""
Shared patterns configuration. The patterns file path is set in main.py
so that generate_csv, plotting, and trendline all use the same file.
"""
import json
import os

# Set by main.py before importing generate_csv, plotting, trendline.
_patterns_file = None


def set_patterns_file(path):
    """Set the patterns JSON file path. Call from main.py before other steps."""
    global _patterns_file
    _patterns_file = path


def get_patterns_file():
    """Return the patterns file path set in main.py (via set_patterns_file)."""
    return _patterns_file


def load_patterns(patterns_file_path=None):
    """Load patterns and general_mutation_map from the configured patterns JSON."""
    path = patterns_file_path if patterns_file_path is not None else get_patterns_file()
    if path is None:
        raise FileNotFoundError(
            "Patterns file not set. In main.py call patterns_config.set_patterns_file(path) "
            "with your patterns JSON path (e.g. telomere_patterns_2x.json) before other steps."
        )
    with open(path, 'r') as f:
        data = json.load(f)
    version = data.get('version', 'unknown')
    return data['patterns'], data['general_mutation_map'], version


def get_patterns_version():
    """Load and return the version string from the configured patterns file (for graph titles)."""
    path = get_patterns_file()
    if path is None:
        return 'unknown'
    try:
        with open(path, 'r') as f:
            return json.load(f).get('version', 'unknown')
    except Exception:
        return 'unknown'
