# Telomere Mutational Signature Analysis

This repository analyzes telomeric mutation signatures from sequencing reads and generates summary CSV outputs plus visualization panels for age- and telomere-associated trends.

The current workflow is centered on `analysis/main.py`, which provides a command-line interface for:

- generating a versioned telomere analysis CSV from FASTQ/FASTA data,
- plotting from an existing CSV, or
- running both steps end-to-end.

## Requirements

- Python 3.10+ (the CLI uses modern type syntax)
- Install dependencies:

```bash
pip install numpy pandas matplotlib seaborn scipy HTSeq
```

## Repository Layout

```text
analysis/
  main.py                                     # CLI entrypoint (generate, plot, run)
  generate_csv.py                             # Sequence parsing + feature/CSV generation
  plotting.py                                 # Plotting pipelines (signatures, spearman, trendlines, curve fitting)
  usage.md                                    # More detailed command reference
  telomere_patterns_2x.json                   # 2x repeat pattern definitions (versioned)
  telomere_patterns_3x.json                   # 3x repeat pattern definitions (versioned)
  greider_methods_table_s2_outliers_removed.csv
  histograms/                                 # Histogram outputs
  plots/                                      # Per-sample mutational signature bar plots
  spearman's plots/                           # Spearman scatterplots + pairwise heatmap
  spearman_correlations/                      # Timestamped 2x2 Spearman panels
  trendlines/                                 # Timestamped 2x2 trendline panels
  curve_fitting_plots/                        # Curve-fitting outputs and summary CSV

greider_data_download/                        # Input FASTQ/FASTA(.gz) files
data/, old/                                   # Legacy and reference data/materials
```

## Inputs

1. Put sequence files in `greider_data_download/`.
   - Supported formats: `.fastq`, `.fastq.gz`, `.fasta`, `.fasta.gz`, `.fa`, `.fa.gz`, `.fas`, `.fas.gz`
2. Provide a patterns file:
   - `analysis/telomere_patterns_2x.json` or
   - `analysis/telomere_patterns_3x.json`
3. (Optional but recommended) Provide metadata CSV:
   - `analysis/greider_methods_table_s2_outliers_removed.csv`
   - Required columns: `fastq file name`, `Age (Years)`, `Mean Telomere Length (bps)`

## CLI Usage

Run commands from the repository root:

```bash
python analysis/main.py <command> [options]
```

### 1) Generate CSV only

```bash
python analysis/main.py generate \
  --patterns analysis/telomere_patterns_2x.json \
  --fastq-dir greider_data_download \
  --metadata analysis/greider_methods_table_s2_outliers_removed.csv
```

- Produces a versioned CSV such as `telomere_analysis_2x_repeat.csv` (or a custom file if `--csv-out` is provided).
- CSV includes raw mutation counts, strand-normalized per-1k rates, grouped mutation features, and engineered fields.

### 2) Plot from existing CSV

```bash
python analysis/main.py plot \
  --patterns analysis/telomere_patterns_2x.json \
  --csv analysis/telomere_analysis_2x_repeat.csv
```

Skip plot groups as needed:

- `--no-hist`
- `--no-signatures`
- `--no-spearman`
- `--no-pairwise`
- `--no-trendlines`
- `--no-curve`

### 3) Full pipeline (generate + plot)

```bash
python analysis/main.py run \
  --patterns analysis/telomere_patterns_2x.json \
  --fastq-dir greider_data_download \
  --metadata analysis/greider_methods_table_s2_outliers_removed.csv
```

## Outputs

- **Main CSV**: `telomere_analysis_<patterns-version>.csv`
- **Per-sample mutational signature plots**: `analysis/plots/`
- **Histogram outputs**: `analysis/histograms/`
- **Spearman outputs**: `analysis/spearman's plots/` and `analysis/spearman_correlations/`
- **Trendline outputs**: `analysis/trendlines/`
- **Curve-fitting outputs**: `analysis/curve_fitting_plots/` (including `curve_fitting_results.csv`)

## Notes

- Pattern version in JSON (for example, `2x_repeat`, `3x_repeat`) drives default output naming.
- Metadata is optional in CSV generation; if unavailable, `Age` and `Telomere_Length` fields are left blank.
- Some output files are timestamped to avoid overwriting previous runs.

## Troubleshooting

- **No sequence files found**: verify `--fastq-dir` path and file extensions.
- **Missing metadata values**: check sample-name normalization (metadata uses names mapped from underscores to dots).
- **Import errors**: reinstall dependencies in the active environment.
- **Unexpected CSV name resolution**: ensure `--patterns` points to the same JSON used to generate the CSV.
