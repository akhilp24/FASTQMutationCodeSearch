# Telomere Mutational Signature Analysis & Age Prediction

This repository provides tools for analyzing mutational signatures in telomeric DNA and predicting sample age from telomeric features. The codebase is organized into two main modules:

- **Analysis** (`/analysis`): Extracts, quantifies, and visualizes telomeric mutations from sequencing data.
- **Modeling** (`/modeling`): Builds and applies machine learning models to predict age from telomeric mutation features.

---

## Table of Contents

1. [Requirements](#requirements)
2. [Directory Structure](#directory-structure)
3. [Data Preparation](#data-preparation)
4. [Analysis Workflow](#analysis-workflow)
5. [Modeling Workflow](#modeling-workflow)
6. [Outputs](#outputs)
7. [Troubleshooting](#troubleshooting)
8. [Next Steps](#next-steps)

---

## Requirements

- Python 3.6+
- Install dependencies:
  ```bash
  pip install matplotlib pandas numpy seaborn scipy scikit-learn joblib
  ```

---

## Directory Structure

```
analysis/
  main.py                  # Main pipeline: processes FASTQ, outputs mutation CSV
  generate_csv.py          # (Advanced) Generate CSV from custom data
  plotting.py              # Individual mutational signature plots
  trendline.py             # Age correlation and trendline plots
  histogram.py             # Age group distribution plots
  scatterplot.py           # (Advanced) Custom scatterplots
  telomere_analysis.csv    # Main output: mutation counts/rates per sample
  telomere_analysis_3x_repeat.csv # (Optional) 3x repeat analysis
  greider_methods_table_s2_outliers_removed.csv    # Sample metadata (age, telomere length)
  plots/                   # Individual signature plots
  output/                  # Summary plots (trendline, histogram, etc.)
  resources/               # Documentation, calculation notes

modeling/
  features.py              # Feature engineering from mutation CSV
  train_model.py           # Train age prediction model
  predict_age.py           # Predict age for new samples
  age_predictor.joblib     # Saved model
  feature_importances.png  # Model feature importance plot
  predicted_vs_actual.png  # Model performance plot
```

---

## Data Preparation

1. **Sequencing Data**: Place your `.fastq` or `.fastq.gz` files in `greider_data_download/`.
2. **Sample Metadata**: Ensure `analysis/greider_methods_table_s2_outliers_removed.csv` exists with columns:
   - `fastq file name`, `Age (Years)`, `Mean Telomere Length (bps)`
   - Filenames should use dots (e.g., `JH47.F86.NB70`).

---

## Analysis Workflow

### 1. Run Main Analysis

Extracts telomeric repeat counts and mutation signatures from all FASTQ files.

```bash
cd analysis
python main.py
```

- **Input**: FASTQ files, `greider_methods_table_s2_outliers_removed.csv`
- **Output**: `telomere_analysis.csv` (mutation counts/rates per sample)

### 2. Generate Plots

#### a. Individual Mutational Signature Plots

```bash
python plotting.py telomere_analysis.csv
```

- Output: Bar plots per sample in `plots/`

#### b. Age Correlation & Trendline Plots

```bash
python trendline.py
```

- Output: `trendline.png`, `spearman_correlation.png` in `output/`

#### c. Age Group Distribution (Boxplots)

```bash
python histogram.py
```

- Output: `histogram.png` in `output/`

#### d. (Optional) Custom Scatterplots

```bash
python scatterplot.py
```

- Output: Custom scatterplots (see script for details)

---

## Modeling Workflow

### 1. Feature Engineering

Extracts features from the mutation CSV for modeling.

```bash
cd ../modeling
python features.py
```

- Output: Feature matrix (see script for details)

### 2. Train Age Prediction Model

```bash
python train_model.py
```

- Output: `age_predictor.joblib`, `feature_importances.png`, `predicted_vs_actual.png`

### 3. Predict Age for New Samples

```bash
python predict_age.py
```

- Output: Predicted ages for input samples

---

## Outputs

- **CSV**: `telomere_analysis.csv` — mutation counts and normalized rates per sample
- **Plots**:
  - Individual signature plots: `analysis/plots/`
  - Summary plots: `analysis/output/` (`trendline.png`, `spearman_correlation.png`, `histogram.png`)
  - Model plots: `modeling/feature_importances.png`, `modeling/predicted_vs_actual.png`
- **Model**: `modeling/age_predictor.joblib`

---

## Troubleshooting

- **No FASTQ files found**: Check `greider_data_download/`
- **Metadata missing**: Ensure `greider_methods_table_s2_outliers_removed.csv` is present and formatted correctly
- **Dependencies**: Install with `pip install -r requirements.txt` (if available)
- **Filename mismatches**: Use dots, not underscores, in sample names

---

## Next Steps

- Explore mutation patterns in healthy aging
- Analyze telomere length vs mutation rates
- Refine age prediction models with additional features

---

## Contact

For questions or contributions, please open an issue or pull request.
