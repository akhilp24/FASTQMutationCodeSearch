## Telomere Analysis CLI

This module lets you generate telomere mutation CSVs from FASTQ/FASTA files and run the plotting/analysis pipeline from the command line.

All commands are typically invoked from the repository root like:

```bash
python analysis/main.py <command> [options...]
```

Commands share these key concepts:

- **Patterns JSON**: telomere pattern definitions (e.g. `analysis/telomere_patterns_2x.json`).
- **Metadata CSV**: age and telomere length table (e.g. `analysis/greider_methods_table_s2_outliers_removed.csv`).
- **FASTQ directory**: folder containing `.fastq(.gz)` or `.fasta(.gz)` files.

---

## 1. Generate CSV only

Use `generate` to create a single `telomere_analysis_*.csv` file (containing both raw and normalized values) from raw sequence data.

```bash
python analysis/main.py generate \
  --patterns analysis/telomere_patterns_2x.json \
  [--metadata analysis/greider_methods_table_s2_outliers_removed.csv] \
  --fastq-dir greider_data_download \
  [--csv-out analysis/telomere_analysis_2x_repeat_feb4_2026.csv]
```

- **`--patterns`** (required): path to the telomere patterns JSON file.
- **`--metadata`** (optional): path to the age/length table CSV.
- **`--fastq-dir`** (required): directory containing FASTQ/FASTA input files.
- **`--csv-out`** (optional): explicit output CSV path.
  - If omitted, the code derives a filename from the patterns `version` field, e.g.  
    `telomere_analysis_2x_repeat_feb4_2026.csv`.

The generated CSV contains:

- `FileName`, `Age`, `Telomere_Length`, `Total_Reads`, `c_strand`, `g_strand`
- All raw mutation count columns (`*_mutations_*`)
- Per‑1k mutation rate columns (`*_per_1k`)
- Grouped per‑1k mutation rates, summed per‑1k features, and engineered features
  (transition/transversion composites, strand sums, telomere bins, etc.).

---

## 2. Run plots on an existing CSV

Use `plot` when you already have a normalized `telomere_analysis_*.csv` and just want to regenerate figures.

```bash
python analysis/main.py plot \
  --patterns analysis/telomere_patterns_2x.json \
  [--csv analysis/telomere_analysis_2x_repeat_feb4_2026.csv] \
  [--no-hist] [--no-signatures] [--no-spearman] \
  [--no-pairwise] [--no-trendlines] [--no-curve]
```

Arguments:

- **`--patterns`** (required): same patterns JSON used when generating the CSV.
  - Used for labelling plots and, if `--csv` is omitted, to infer the CSV filename.
- **`--csv`** (optional): explicit CSV path. If omitted, the CSV path is
  inferred from the patterns version (expects the single combined CSV naming).

Optional flags (to _skip_ certain plot groups):

- **`--no-hist`**: skip age‑binned mutation histograms and per‑file mutation bar plot.
  - Outputs (by default):
    - `analysis/histograms/histogram.png`
    - `analysis/histograms/mutations_per_file_histogram.png`
- **`--no-signatures`**: skip mutational signature barplots for each sample.
- **`--no-spearman`**: skip simple Spearman vs Age scatter plots.
- **`--no-pairwise`**: skip pairwise Spearman correlation heatmap.
- **`--no-trendlines`**: skip 2×2 trendline + R² grid and 2×2 Spearman grids.
- **`--no-curve`**: skip curve‑fitting analysis (telomere length vs age and key mutation metrics vs age).

When **not** skipped, plots are written into these locations:

- Mutational signatures: `analysis/plots/`
- Spearman scatter plots and composite score: `analysis/spearman's plots/`
- Spearman heatmaps: `analysis/spearman's plots/`
- Trendlines + 4‑panel Spearman grids: `analysis/trendlines/` and `analysis/spearman_correlations/`
- Curve fitting: `analysis/curve_fitting_plots/`

---

## 3. Full pipeline: generate CSV + run plots

Use `run` to go from raw FASTQ/FASTA to a CSV and all plots in a single command:

```bash
python analysis/main.py run \
  --patterns analysis/telomere_patterns_2x.json \
  [--metadata analysis/greider_methods_table_s2_outliers_removed.csv] \
  --fastq-dir greider_data_download \
  [--csv-out analysis/telomere_analysis_2x_repeat_feb4_2026.csv] \
  [--no-hist] [--no-signatures] [--no-spearman] \
  [--no-pairwise] [--no-trendlines] [--no-curve]
```

This:

1. Calls the **`generate`** step to build the combined CSV (raw + normalized).
2. Calls the **`plot`** step with the same patterns file and the freshly generated CSV.

Flags such as `--no-hist`, `--no-trendlines`, etc., behave the same as for `plot`.

---

## 4. Typical workflows

- **2× repeat analysis from scratch:**

  ```bash
  python analysis/main.py run \
    --patterns analysis/telomere_patterns_2x.json \
    --metadata analysis/greider_methods_table_s2_outliers_removed.csv \
    --fastq-dir greider_data_download
  ```

- **Rerun only plots on an existing 2× CSV:**

  ```bash
  python analysis/main.py plot \
    --patterns analysis/telomere_patterns_2x.json \
    --csv analysis/telomere_analysis_2x_repeat_feb4_2026_normalized.csv
  ```

- **Just generate a 3× repeat CSV (no plots):**

  ```bash
  python analysis/main.py generate \
    --patterns analysis/telomere_patterns_3x.json \
    --metadata analysis/greider_methods_table_s2_outliers_removed.csv \
    --fastq-dir greider_data_download \
    --csv-out analysis/telomere_analysis_3x_repeat_custom.csv
  ```

Refer to `analysis/generate_csv.py` and `analysis/plotting.py` if you need to customize which columns or variables are used in specific plots.
