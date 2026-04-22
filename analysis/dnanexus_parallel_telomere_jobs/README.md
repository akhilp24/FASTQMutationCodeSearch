# DNANexus Parallel UKB CRAM Telomere Runs

This folder runs the updated telomere CRAM logic from `analysis/telomere_analysis_ukb_cram_update_fixed.ipynb` in **parallel** on DNANexus, then creates a final combined CSV and plotting outputs.

## Files

- `submit_parallel_dnanexus_jobs.py` - parallel fan-out/fan-in launcher (IDs txt -> jobs -> combined CSV)
- `run_single_cram_analysis.py` - per-participant CRAM worker using phase-aware k-repeat analysis
- `plot_from_combined_csv.py` - plotting pipeline from combined CSV
- `participants_4_example.txt` - starter 4-ID input file

## Why this is robust for UKB scale

- Uses one DNANexus job per participant (embarrassingly parallel fan-out)
- No serial loop processing of CRAM content
- Per-sample CSVs are isolated and then merged only after all jobs finish
- Keeps final outputs in a timestamped run folder for reproducibility

## Prerequisites

- `dx` CLI installed and logged in (`dx whoami` works)
- Local Python environment with `pandas` available
- Access to project `project-J0J4VJ8J7pf8gJKjjZgYjkzJ`
- DNANexus metadata files used by default:
  - Reference FASTA: `file-G572Pj8JykJZ52BP4z6GqB21`
  - Participant metadata CSV (age/sex/smoking/telomere fields): `file-J4PfY50J7pfKqx2KzZpzy5zb`

## 1) Prepare participant IDs

Create a text file with one participant ID per line.

## 2) Submit parallel jobs (4-patient pilot)

From repo root:

```bash
python "analysis/dnanexus_parallel_telomere_jobs/submit_parallel_dnanexus_jobs.py" \
  --participant-ids "analysis/dnanexus_parallel_telomere_jobs/participants_4_example.txt" \
  --project "project-J0J4VJ8J7pf8gJKjjZgYjkzJ" \
  --cram-folder "/Bulk/GATK and GraphTyper WGS/Whole genome GATK CRAM files and indices [500k release]/10/" \
  --output-base-folder "/Telomere/parallel_ukb_runs" \
  --instance-type "mem2_hdd2_x4" \
  --k 3 \
  --combined-local-csv "analysis/dnanexus_parallel_telomere_jobs/combined_telomere_results.csv"
```

What this does:

1. Resolves each participant's CRAM in the UKB folder
2. Launches all participant jobs in parallel (not in series)
3. Waits for all jobs to complete
4. Downloads per-sample CSVs and writes one local combined CSV
5. Uploads combined CSV to `.../run_<timestamp>/combined/` on DNANexus

The launcher now creates the DNANexus output folders automatically before uploads/jobs.

### Why metadata CSV is needed

The CRAM contains sequencing reads, but not the phenotype fields needed for downstream analysis plots (age and related covariates).  
This pipeline uses your metadata file (`file-J4PfY50J7pfKqx2KzZpzy5zb`) to add:

- `Age` (from `Age at recruitment` when present)
- `Telomere_Length` proxy (auto-selected from available telomere-style columns)
- additional participant fields such as sex and smoking status

## 3) Generate plots from final combined CSV

```bash
python "analysis/dnanexus_parallel_telomere_jobs/plot_from_combined_csv.py" \
  --combined-csv "analysis/dnanexus_parallel_telomere_jobs/combined_telomere_results.csv" \
  --output-dir "analysis/dnanexus_parallel_telomere_jobs/plots"
```

This creates histogram, Spearman, heatmap, and trendline outputs compatible with your notebook-style downstream plotting workflow.

## Notes

- If a participant ID has no CRAM match in the target folder, the launcher fails fast.
- By default, IDs with no CRAM match are skipped and reported; use `--strict-missing` to fail fast.
- You can scale by increasing the number of IDs in the text file.
- For large runs, tune `--instance-type` for memory/throughput trade-offs.
