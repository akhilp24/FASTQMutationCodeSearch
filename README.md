# Telomere Mutational Signature Analysis

This repository contains code for analyzing mutational signatures in telomeric DNA sequences from FASTQ files. The analysis identifies and quantifies all possible single base substitutions (transitions and transversions) at specific positions within telomeric repeat sequences.

## Overview

The code performs comprehensive mutational signature analysis by:

- Counting all 6 possible single base substitutions (C>A, C>G, C>T, T>A, T>C, T>G)
- Analyzing mutations at 3 specific positions within telomeric repeat sequences
- Processing both G-strand and C-strand contexts
- Normalizing mutation rates per 1000 telomeric repeats
- Generating mutational signature plots

## Requirements

- Python 3.6+
- Required packages:
  - `matplotlib`
  - `pandas`
  - `numpy`
  - Standard library modules: `gzip`, `collections`, `csv`, `os`, `glob`

## Input Data

### Required Files

1. **FASTQ Files**: Place your FASTQ files (`.fastq` or `.fastq.gz`) in the `greider_data/` directory
2. **Age Data**: `greider_methods_table_s2.csv` - Contains age information for each sample
   - Must have columns: `fastq file name`, `Age (Years)`
   - Filename format should match: `JH47.F86.NB70` (dots instead of underscores)

### Data Format

The age data CSV should have this structure:

```csv
fastq file name,Age (Years)
JH47.F86.NB70,3
JH65.F86.NB70,90
...
```

## Usage

### 1. Run the Main Analysis

```bash
python main.py
```

This will:

- Process all FASTQ files in the `greider_data/` directory
- Count telomeric repeat sequences and mutations
- Generate `telomere_analysis.csv` with detailed mutation counts and rates
- Print summary statistics to console

### 2. Generate Mutational Signature Plot

```bash
python plotting.py telomere_analysis.csv [output_filename.png]
```

This will:

- Read the analysis results from `telomere_analysis.csv`
- Create a mutational signature plot showing the percentage of each mutation type
- Save the plot as `mutational_signatures.png` (or custom filename)

## Output

### CSV Output (`telomere_analysis.csv`)

The analysis generates a comprehensive CSV with the following columns:

#### Basic Information

- `FileName`: Name of the processed FASTQ file
- `Age`: Age in years from the metadata
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

### Plot Output

The mutational signature plot shows:

- Percentage of each of the 6 mutation types (C>A, C>G, C>T, T>A, T>C, T>G)
- Color-coded bars following standard mutational signature conventions
- Aggregated data from all positions and both strand contexts

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
ls greider_data/*.fastq*
ls greider_methods_table_s2.csv

# 2. Run the analysis
python main.py

# 3. Generate the plot
python plotting.py telomere_analysis.csv

# 4. View results
head telomere_analysis.csv
open mutational_signatures.png
```

## Notes

- The code processes both compressed (`.gz`) and uncompressed FASTQ files
- Mutation rates are normalized per 1000 telomeric repeats to account for different sequence depths
- The analysis includes both forward and reverse complement sequences
- All mutation types are counted regardless of frequency to provide complete mutational signatures

## Troubleshooting

- **No FASTQ files found**: Ensure files are in the `greider_data/` directory
- **Age data not found**: Check that `greider_methods_table_s2.csv` exists and has correct column names
- **Filename mismatch**: Convert underscores to dots in the age data CSV to match FASTQ filenames
