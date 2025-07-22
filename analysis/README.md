# Telomere Mutational Signature Analysis

This repository contains code for analyzing mutational signatures in telomeric DNA sequences from FASTQ files. The analysis identifies and quantifies all possible single base substitutions (transitions and transversions) at specific positions within telomeric repeat sequences.

## Overview

The code performs comprehensive mutational signature analysis by:

- Counting all 6 possible single base substitutions (C>A, C>G, C>T, T>A, T>C, T>G)
- Analyzing mutations at 3 specific positions within telomeric repeat sequences
- Processing both G-strand and C-strand contexts
- Normalizing mutation rates per 1000 telomeric repeats
- Generating multiple types of visualizations including individual mutational signatures, trendlines, correlations, and age-group distributions

## Requirements

- Python 3.6+
- Required packages:
  - `matplotlib`
  - `pandas`
  - `numpy`
  - `seaborn`
  - `scipy`
  - Standard library modules: `gzip`, `collections`, `csv`, `os`, `glob`

## Input Data

### Required Files

1. **FASTQ Files**: Place your FASTQ files (`.fastq` or `.fastq.gz`) in the `greider_data_download/` directory
2. **Age Data**: `greider_methods_table_s2.csv` - Contains age and telomere length information for each sample
   - Must have columns: `fastq file name`, `Age (Years)`, `Mean Telomere Length (bps)`
   - Filename format should match: `JH47.F86.NB70` (dots instead of underscores)

### Data Format

The age data CSV should have this structure:

```csv
fastq file name,Age (Years),Mean Telomere Length (bps)
JH47.F86.NB70,3,5000
JH65.F86.NB70,90,3000
...
```

## Usage

### 1. Run the Main Analysis

```bash
python main.py
```

This will:

- Process all FASTQ files in the `greider_data_download/` directory
- Count telomeric repeat sequences and mutations
- Generate `telomere_analysis.csv` with detailed mutation counts and rates
- Print summary statistics to console

### 2. Generate Individual Mutational Signature Plots

```bash
python plotting.py telomere_analysis.csv
```

This will:

- Read the analysis results from `telomere_analysis.csv`
- Create individual mutational signature plots for each sample
- Save plots in the `plots/` directory with filenames matching the FASTQ files
- Each plot shows the percentage of each mutation type by position and strand context

### 3. Generate Age Correlation Analysis

```bash
python trendline.py
```

This will:

- Read the analysis results from `telomere_analysis.csv`
- Create two plots:
  - `trendline.png`: Scatter plots with trendlines showing mutation rates vs age
  - `spearman_correlation.png`: Scatter plots with Spearman correlation statistics
- Focus on total mutations and G>T mutations at positions 1, 2, and 3

### 4. Generate Age Group Distribution Analysis

```bash
python histogram.py
```

This will:

- Read the analysis results from `telomere_analysis.csv`
- Create `histogram.png` with box plots showing mutation rate distributions by age groups (10-year bins)
- Analyze the same variables as trendline.py but grouped by age ranges

## Output

### CSV Output (`telomere_analysis.csv`)

The analysis generates a comprehensive CSV with the following columns:

#### Basic Information

- `FileName`: Name of the processed FASTQ file
- `Age`: Age in years from the metadata
- `Telomere_Length`: Mean telomere length in base pairs
- `2x_cstrand`: Total count of C-strand telomeric repeats
- `2xg_strand`: Total count of G-strand telomeric repeats

#### Mutation Counts

For each mutation type (C>A, C>G, C>T, T>A, T>C, T>G) at positions 1, 2, and 3:

**G-strand mutations:**

- `G_A_g1`, `G_A_g2`, `G_A_g3`: G→A mutations at positions 1, 2, 3
- `G_C_g1`, `G_C_g2`, `G_C_g3`: G→C mutations at positions 1, 2, 3
- `G_T_g1`, `G_T_g2`, `G_T_g3`: G→T mutations at positions 1, 2, 3

**C-strand mutations:**

- `C_A_c1`, `C_A_c2`, `C_A_c3`: C→A mutations at positions 1, 2, 3
- `C_G_c1`, `C_G_c2`, `C_G_c3`: C→G mutations at positions 1, 2, 3
- `C_T_c1`, `C_T_c2`, `C_T_c3`: C→T mutations at positions 1, 2, 3

**T-strand mutations:**

- `T_A_t1`, `T_A_t2`, `T_A_t3`: T→A mutations at positions 1, 2, 3
- `T_C_t1`, `T_C_t2`, `T_C_t3`: T→C mutations at positions 1, 2, 3
- `T_G_t1`, `T_G_t2`, `T_G_t3`: T→G mutations at positions 1, 2, 3

#### Normalized Rates

For each mutation type, normalized rates per 1000 telomeric repeats:

- `*_per_1k`: Rate per 1000 repeats for each mutation type and position
- `total_mutations_over_total_g_per_1k`: Overall mutation rate per 1000 G-strand repeats

### Plot Outputs

#### Individual Mutational Signatures (`plots/` directory)

- One plot per sample showing percentage of each mutation type
- Color-coded bars following standard mutational signature conventions
- Includes both G-strand and C-strand contexts for each mutation type
- Shows mutations at all three positions within telomeric repeats

#### Age Correlation Analysis

- **`trendline.png`**: 2x2 subplot showing scatter plots with trendlines for:

  - Total mutations per 1000bp vs age
  - G>T mutations at position 1 vs age
  - G>T mutations at position 2 vs age
  - G>T mutations at position 3 vs age

- **`spearman_correlation.png`**: Same variables as trendline.png but with Spearman correlation statistics displayed in titles

#### Age Group Distribution Analysis

- **`histogram.png`**: 2x2 subplot showing box plots for the same variables as trendline.py, but grouped by 10-year age bins

## Telomeric Repeat Patterns

The analysis searches for these specific telomeric repeat patterns:

- **G-strand pattern**: `GGGTTAGGGTTA`
- **C-strand pattern**: `CCCTAACCCTAA`

Mutations are counted at the three G positions in the G-strand pattern and the three C positions in the C-strand pattern.

## Mutation Types

The analysis covers all possible single base substitutions:

**Transitions (purine↔purine, pyrimidine↔pyrimidine):**

- A↔G (purine transition)
- C↔T (pyrimidine transition)

**Transversions (purine↔pyrimidine):**

- A↔C, A↔T, G↔C, G↔T

## Example Workflow

```bash
# 1. Ensure your data is in place
ls greider_data_download/*.fastq*
ls greider_methods_table_s2.csv

# 2. Run the main analysis
python main.py

# 3. Generate individual mutational signature plots
python plotting.py telomere_analysis.csv

# 4. Generate age correlation analysis
python trendline.py

# 5. Generate age group distribution analysis
python histogram.py

# 6. View results
head telomere_analysis.csv
ls plots/
open trendline.png spearman_correlation.png histogram.png
```

## Notes

- The code processes both compressed (`.gz`) and uncompressed FASTQ files
- Mutation rates are normalized per 1000 telomeric repeats to account for different sequence depths
- The analysis includes both forward and reverse complement sequences
- All mutation types are counted regardless of frequency to provide complete mutational signatures
- The plotting scripts use seaborn for enhanced visual styling
- Individual mutational signature plots are saved in a `plots/` directory with filenames matching the original FASTQ files

## Troubleshooting

- **No FASTQ files found**: Ensure files are in the `greider_data_download/` directory
- **Age data not found**: Check that `greider_methods_table_s2.csv` exists and has correct column names
- **Filename mismatch**: Convert underscores to dots in the age data CSV to match FASTQ filenames
- **Missing dependencies**: Install required packages with `pip install matplotlib pandas numpy seaborn scipy`

## Next Steps

- Characterize super healthy individuals with low mutations/high age and see what their relative mutation rates are across positions
- Analyze the relationship between telomere length and mutation rates
- Investigate position-specific mutation patterns in different age groups
