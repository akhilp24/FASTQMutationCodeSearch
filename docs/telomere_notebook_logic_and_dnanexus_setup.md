# Telomere Notebook Logic and DNAnexus Setup

This note summarizes the logic in `analysis/telomere_analysis_ukb_cram_update_fixed.ipynb` and gives practical setup steps for running it on DNAnexus from this GitHub repository.

## Purpose of the notebook

The notebook is a two-stage telomere workflow:

1. **Generate analysis table** from UK Biobank CRAM input:
   - Detect telomeric reads on G and C strands
   - Count canonical and mutant telomeric patterns
   - Normalize mutation counts
   - Write a per-sample CSV and a human-readable summary text file
2. **Plot outputs** from the generated CSV:
   - Histograms / boxplots by age bins
   - Age trendline plots
   - Spearman correlation plots with age
   - Pairwise correlation heatmap
   - Per-sample mutational signature barplots

The notebook is designed for UKB-style CRAM data and supports either:
- a **single CRAM file**, or
- a **directory of CRAM files**.

## Pipeline logic (code walkthrough)

## 1) Environment and imports

Early cells install required Python packages (`pysam`, `pandas`, `numpy`, `seaborn`, `scipy`, `matplotlib`) and import analysis libraries.

Why this matters:
- `pysam` is required for reading CRAM/BAM records
- `pandas/numpy` handle tabular aggregation and normalization
- `scipy/seaborn/matplotlib` support downstream statistics + plotting

## 2) Configuration and metadata loading

The configuration cell defines key runtime parameters:

- `CRAM_INPUT`: path to one CRAM file or a directory of CRAMs
- `OUTPUT_DIR`: destination directory for CSV, plots, and summary text
- `REPEAT_K`: number of telomeric repeat units in a motif window (default shown: 3)
- `REFERENCE_FASTA`: optional CRAM reference for decoding
- `CSV_OUT`, `OUTPUT_TXT`: main output files

The notebook also loads UKB metadata (`/mnt/project/Telomere/Telomere_length_age.csv`) to attach covariates (especially age) to each sample.

## 3) Pattern-generation strategy (core analytical change)

The central logic improves telomere counting accuracy with a **phase-aware k-repeat approach**.

### Biological/computational issue addressed

Telomeric hexamers can appear in any rotational phase in sequencing reads. If only one phase is searched, most signal is missed.  
At the same time, searching many short mutant motifs inflates false positives in random genomic sequence.

### Implemented solution

The notebook:

- builds canonical G-strand (`GGGTTA`) and C-strand (`CCCTAA`) motifs repeated `k` times
- generates **all 6 rotational phases** for each canonical pattern
- creates mutation-specific variants (single-base substitutions in the second repeat)
- keeps quantification at **k-repeat length** (e.g., 18-mer for `k=3`) rather than single hexamer mutation counting

Net effect:
- better sensitivity to true telomeric reads (phase-aware)
- better specificity for mutations (longer k-mer matching)

## 4) CRAM scanning and per-sample counting

Generation cells then process reads from each CRAM:

1. iterate reads from `pysam`
2. pre-screen reads as potentially telomeric
3. match read sequence against canonical/mutant k-repeat phase patterns
4. accumulate per-sample counters, including:
   - total reads
   - G-strand and C-strand telomeric reads
   - mutation subtype counts
   - aggregate mutation totals
5. write per-sample records to the output CSV

A quick-check cell validates that CRAM paths and FASTA path exist before running heavy computation.

## 5) Summary text output

The `_write_output_summary(...)` function writes a plain-text report (`output.txt`) containing:

- run metadata (version, repeat-k, input mode)
- total reads and telomeric counts
- mutation totals and rates
- sample-level and cohort-level summaries

This provides a PI-readable QC and run summary without needing to parse the full CSV.

## 6) Plotting layer

After CSV generation, plotting functions consume the generated table and produce publication-ready visual summaries:

- age-binned mutation distribution plots
- age vs mutation trendlines
- Spearman correlation plots for numeric variables vs age
- pairwise correlation heatmap across mutation-derived variables
- per-sample mutational signature bars

Output directories are created inside `OUTPUT_DIR` (or configured subfolders).

## 7) End-to-end execution behavior

The notebook is organized as:

- **[1/2] Generate CSV** from CRAM input
- optional `dx upload output.txt --path Telomere/` convenience upload
- **[2/2] Run plots** from generated CSV

This means one run gives both quantitative table outputs and visual diagnostics.

## Expected inputs and outputs

### Inputs

- UKB CRAM data (`CRAM_INPUT`)
- optional CRAM reference FASTA (`REFERENCE_FASTA`)
- metadata table (`Telomere_length_age.csv`) for age annotation

### Outputs

- `telomere_analysis_ukb.csv` (main per-sample analysis table)
- `output.txt` (human-readable run summary)
- plot files in subdirectories (histograms, trendlines, Spearman plots, heatmap, signatures)

---

## DNAnexus setup instructions (GitHub pull workflow)

These steps describe how to load and run this notebook on DNAnexus by pulling this repository from GitHub.

## 1) Create/enter a DNAnexus environment

1. Open your DNAnexus project.
2. Launch a JupyterLab (or terminal-enabled) environment with permission to read the required UKB data paths.
3. Ensure `dx` authentication is active in that environment.

## 2) Pull this GitHub repository

In the DNAnexus terminal:

```bash
git clone <YOUR_REPO_URL>
cd GreiderCodeSearch
```

If already cloned:

```bash
cd GreiderCodeSearch
git pull
```

## 3) Open the notebook

Navigate to:

`analysis/telomere_analysis_ukb_cram_update_fixed.ipynb`

Run cells top-to-bottom the first time.

## 4) Update configuration paths for your DNAnexus project

In the UKB configuration cell, set:

- `CRAM_INPUT` to your target CRAM file or directory in mounted DNAnexus storage
- `REFERENCE_FASTA` to the correct FASTA path (or `None` if reference-free CRAMs are supported)
- `OUTPUT_DIR` to a writable location (commonly `/opt/notebooks` in DNAnexus Jupyter sessions)
- `REPEAT_K` to desired repeat size (default in notebook is 3)

Also confirm metadata path exists:

- `/mnt/project/Telomere/Telomere_length_age.csv`

If your project layout differs, update this path in the metadata-loading cell.

## 5) Execute analysis

1. Run package-install and import cells.
2. Run configuration + metadata cells.
3. Run quick-check cell to verify file existence.
4. Run generation cell (`[1/2]`) and wait for completion.
5. Run plotting cell (`[2/2]`) to generate visual outputs.

## 6) Upload/export outputs to DNAnexus project storage

Example upload command from notebook terminal:

```bash
dx upload /opt/notebooks/output.txt --path Telomere/
```

You can similarly upload:

- `/opt/notebooks/telomere_analysis_ukb.csv`
- plot directories under `/opt/notebooks/` (or your configured output path)

## 7) Reproducibility recommendations

For reproducible reruns and PI reporting:

- keep `REPEAT_K` fixed across comparisons unless intentionally testing sensitivity
- record commit hash of the repository version used
- preserve `output.txt` and the CSV together with generated plots for each run

---

## Short PI-facing interpretation

In plain terms: this notebook quantifies telomeric sequence content and mutation-like deviations from UKB CRAM reads using phase-aware telomere motifs, then relates these metrics to age with multiple statistical visualizations.  
The key methodological point is that it combines **all 6 motif phases** (sensitivity) with **k-repeat-length matching** (specificity), which is intended to improve biological signal recovery while reducing random-sequence false positives.
