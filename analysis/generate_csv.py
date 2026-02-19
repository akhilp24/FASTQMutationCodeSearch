#!/usr/bin/env python3

import gzip
import json
from collections import defaultdict
import csv
import os
import numpy as np
import glob
import HTSeq


def load_patterns(patterns_file_path):
    """Load patterns and general_mutation_map from the given patterns JSON path."""
    with open(patterns_file_path, 'r') as f:
        data = json.load(f)
    version = data.get('version', 'unknown')
    return data['patterns'], data['general_mutation_map'], version


def read_sequence_file(file_path: str):
    """Read FASTQ or FASTA file and yield sequences."""
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    # Detect format by file extension
    is_fasta = any(ext in file_path.lower() for ext in ['.fasta', '.fa', '.fas'])
    
    with open_func(file_path, mode) as f:
        if is_fasta:
            # FASTA format
            current_sequence = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # New sequence header - yield previous sequence if exists
                    if current_sequence:
                        yield current_sequence
                    current_sequence = ""
                else:
                    # Sequence line - append to current sequence
                    current_sequence += line
            # Yield the last sequence
            if current_sequence:
                yield current_sequence
        else:
            # FASTQ format
            while True:
                header = f.readline().strip()
                if not header:
                    break
                sequence = f.readline().strip()
                _ = f.readline()  # + line
                _ = f.readline()  # quality line
                yield sequence

def count_patterns(sequence: str, pattern: str) -> int:
    return sequence.count(pattern)

def count_total_reads(file_path: str) -> int:
    """Count the total number of reads in a FASTQ or FASTA file."""
    count = 0
    for _ in HTSeq.FastqReader(file_path):
        count += 1
    return count

def load_age_data(metadata_file_path=None):
    """Load age data from greider_methods_table_s2_outliers_removed.csv."""
    age_data = {}
    
    # Try to find the metadata file
    if metadata_file_path is None:
        # Look in current directory first
        if os.path.exists('greider_methods_table_s2_outliers_removed.csv'):
            metadata_file_path = 'greider_methods_table_s2_outliers_removed.csv'
        # Look in analysis directory
        elif os.path.exists('../analysis/greider_methods_table_s2_outliers_removed.csv'):
            metadata_file_path = '../analysis/greider_methods_table_s2_outliers_removed.csv'
        # Look in parent directory
        elif os.path.exists('analysis/greider_methods_table_s2_outliers_removed.csv'):
            metadata_file_path = 'analysis/greider_methods_table_s2_outliers_removed.csv'
        else:
            raise FileNotFoundError("greider_methods_table_s2_outliers_removed.csv not found in current directory, analysis directory, or parent directory")
    
    with open(metadata_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            fastq_name = row['fastq file name'].replace('_', '.')
            age_data[fastq_name] = row['Age (Years)']
    return age_data

def load_length_data(metadata_file_path=None):
    """Load length data from greider_methods_table_s2_outliers_removed.csv."""
    length_data = {}
    
    # Try to find the metadata file
    if metadata_file_path is None:
        # Look in current directory first
        if os.path.exists('greider_methods_table_s2_outliers_removed.csv'):
            metadata_file_path = 'greider_methods_table_s2_outliers_removed.csv'
        # Look in analysis directory
        elif os.path.exists('../analysis/greider_methods_table_s2_outliers_removed.csv'):
            metadata_file_path = '../analysis/greider_methods_table_s2_outliers_removed.csv'
        # Look in parent directory
        elif os.path.exists('analysis/greider_methods_table_s2_outliers_removed.csv'):
            metadata_file_path = 'analysis/greider_methods_table_s2_outliers_removed.csv'
        else:
            raise FileNotFoundError("greider_methods_table_s2_outliers_removed.csv not found in current directory, analysis directory, or parent directory")
    
    with open(metadata_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            fastq_name = row['fastq file name'].replace('_', '.')
            length_data[fastq_name] = row['Mean Telomere Length (bps)']
    return length_data

def get_sequence_files(directory: str):
    """Get all FASTQ and FASTA files in the given directory."""
    sequence_files = []
    
    # FASTQ files
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fastq")))
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fastq.gz")))
    
    # FASTA files
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fasta")))
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fasta.gz")))
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fa")))
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fa.gz")))
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fas")))
    sequence_files.extend(glob.glob(os.path.join(directory, "*.fas.gz")))
    
    return sorted(sequence_files)  # Sort for consistent ordering

def generate_csv(
    data_dir: str,
    output_callback=None,
    metadata_file_path=None,
    patterns_file_path=None,
    output_csv_path=None,
):
    """
    Generate a single CSV file from sequence data containing both raw counts and
    normalized (per-1k) metrics.

    Args:
        data_dir: Directory containing sequence files.
        output_callback: Optional callback function to receive console output.
        metadata_file_path: Optional path to metadata file (age/length table). If not
            provided or not found, Age and Telomere_Length will be left empty.
        patterns_file_path: Path to patterns JSON (e.g. from main.py).
        output_csv_path: Optional path for output CSV. If None, writes to a file whose
            name encodes the patterns version, e.g. 'telomere_analysis_<version>.csv'.
            If provided, parent directory is created if needed.
    """
    if patterns_file_path is None:
        raise ValueError("patterns_file_path is required; pass it from main.py.")

    # Load patterns from reference file (also provides version for filenames).
    patterns, general_mutation_map, patterns_version = load_patterns(patterns_file_path)

    # Decide output path for the combined CSV (raw + normalized).
    if output_csv_path is None:
        safe_version = patterns_version.replace(' ', '_')
        csv_path = f"telomere_analysis_{safe_version}.csv"
    else:
        csv_path = output_csv_path
        parent = os.path.dirname(csv_path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    sequence_files = get_sequence_files(data_dir)
    
    if not sequence_files:
        message = f"No FASTQ or FASTA files found in {data_dir} directory"
        print(message)
        if output_callback:
            output_callback(message)
        return
    
    # Load age/length metadata if available; otherwise treat as optional.
    try:
        age_data = load_age_data(metadata_file_path)
        length_data = load_length_data(metadata_file_path)
    except FileNotFoundError:
        age_data = {}
        length_data = {}

    # --- Define headers for the combined CSV ---
    fieldnames = [
        'FileName', 'Age', 'Telomere_Length', 'Total_Reads',
        'c_strand', 'g_strand',
    ]

    mutation_keys = []
    for group in ['g_strand_mutations', 'c_strand_mutations']:
        for subkey in patterns[group].keys():
            mutation_keys.append(f"{group}_{subkey}")
    # Raw mutation counts
    fieldnames.extend(mutation_keys)

    # Normalized rate fields for each mutation key
    fieldnames.extend([f"{k}_per_1k" for k in mutation_keys])

    # --- Add general mutation per_1k columns ---
    general_mutation_headers = []
    for strand, mutmap in general_mutation_map.items():
        for mut in mutmap:
            general_mutation_headers.append(f"{strand}_{mut}_per_1k")
    fieldnames.extend(general_mutation_headers)

    # Add engineered features to the CSV header
    fieldnames.extend([
        'g_strand_mutations_sum_per_1k',
        'c_strand_mutations_sum_per_1k',
        'log_telomere_length',
        'telomere_length_bin',
        'mutation_rate_normalized_by_length',
    ])

    # Add summed per-1k columns for each mutation type per strand
    summed_per_1k_headers = []
    for strand, mutmap in general_mutation_map.items():
        for mut, subtypes in mutmap.items():
            summed_per_1k_headers.append(f"{strand}_{mut}_sum_per_1k")
    fieldnames.extend(summed_per_1k_headers)

    # Add total mutations field at the end (now strand-specific normalized)
    fieldnames.append('total_mutations_per_1k_strand_specific')

    # Add legacy-style total mutation rate fields without hard-coding repeat count in the name.
    # The repeat context (2x vs 3x) is conveyed by the output CSV filename/patterns version instead.
    fieldnames.extend([
        'total_mutations_over_total_g_strand_per_1k',
        'total_mutations_over_total_c_strand_per_1k',
    ])
        
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for file_path in sequence_files:
            counts = defaultdict(int)
            
            # Count total reads in the file
            total_reads = count_total_reads(file_path)
            
            for sequence in read_sequence_file(file_path):  # Changed function name
                # Count c-strand in forward direction only
                counts['c_strand'] += count_patterns(sequence, patterns['c_strand'])
                # Count g-strand in forward direction only
                counts['g_strand'] += count_patterns(sequence, patterns['g_strand'])
                # Count all mutation sub-patterns (frameshifts disabled)
                for group in ['g_strand_mutations', 'c_strand_mutations']:
                    for subkey, subpattern in patterns[group].items():
                        counts[f"{group}_{subkey}"] += count_patterns(sequence, subpattern)
            
            filename = os.path.basename(file_path)
            # Handle both FASTQ and FASTA extensions
            filename_base = filename
            for ext in ['.fastq.gz', '.fastq', '.fasta.gz', '.fasta', '.fa.gz', '.fa', '.fas.gz', '.fas']:
                if filename_base.endswith(ext):
                    filename_base = filename_base[:-len(ext)]
                    break
            
            age = age_data.get(filename_base, '')
            length = length_data.get(filename_base, '')
            g_strand_total = counts['g_strand']
            c_strand_total = counts['c_strand']

            # Calculate strand-specific normalization denominators
            g_strand_mutations_total = sum(counts[k] for k in counts if k.startswith('g_strand_mutations_'))
            c_strand_mutations_total = sum(counts[k] for k in counts if k.startswith('c_strand_mutations_'))
            
            g_strand_normalizer = g_strand_total + g_strand_mutations_total
            c_strand_normalizer = c_strand_total + c_strand_mutations_total

            # New strand-specific normalization function
            def per_1k_strand_specific(val, strand_normalizer):
                return (val / strand_normalizer) * 1000 if strand_normalizer > 0 else 0
            
            # Keep the old function for backward compatibility with strand-based normalization
            def per_1k(val, total):
                return (val / total) * 1000 if total > 0 else 0
                
            row = {
                'FileName': filename,
                'Age': age,
                'Telomere_Length': length,
                'Total_Reads': total_reads,
                'c_strand': c_strand_total,
                'g_strand': g_strand_total,
            }
            
            for k in mutation_keys:
                row[k] = counts.get(k, 0)
                # Use strand-specific normalization for mutations
                if k.startswith('g_strand_mutations_'):
                    row[f"{k}_per_1k"] = per_1k_strand_specific(counts.get(k, 0), g_strand_normalizer)
                elif k.startswith('c_strand_mutations_'):
                    row[f"{k}_per_1k"] = per_1k_strand_specific(counts.get(k, 0), c_strand_normalizer)
                else:
                    # Fallback for any other keys
                    row[f"{k}_per_1k"] = per_1k_strand_specific(counts.get(k, 0), g_strand_normalizer)

            # --- General mutation per_1k columns ---
            for strand, mutmap in general_mutation_map.items():
                for mut, subtypes in mutmap.items():
                    # Regular mutations only (frameshifts disabled)
                    total = sum(counts.get(f"{strand}_mutations_{subtype}", 0) for subtype in subtypes)
                    # Use strand-specific normalization
                    if strand == 'g_strand':
                        row[f"{strand}_{mut}_per_1k"] = per_1k_strand_specific(total, g_strand_normalizer)
                    else:  # c_strand
                        row[f"{strand}_{mut}_per_1k"] = per_1k_strand_specific(total, c_strand_normalizer)

            # --- Summed per-1k columns for each mutation type per strand ---
            for strand, mutmap in general_mutation_map.items():
                for mut, subtypes in mutmap.items():
                    # Regular mutations only (frameshifts disabled)
                    per_1k_sum = sum(
                        row.get(f"{strand}_mutations_{subtype}_per_1k", 0)
                        for subtype in subtypes
                    )
                    row[f"{strand}_{mut}_sum_per_1k"] = per_1k_sum
            
            # Total mutations (sum all mutation counts)
            total_mutations = sum(counts[k] for k in mutation_keys)
            # Use weighted average of strand-specific normalizations
            if g_strand_normalizer > 0 or c_strand_normalizer > 0:
                g_weight = g_strand_mutations_total / total_mutations if total_mutations > 0 else 0
                c_weight = c_strand_mutations_total / total_mutations if total_mutations > 0 else 0
                weighted_normalizer = (g_weight * g_strand_normalizer) + (c_weight * c_strand_normalizer)
                row['total_mutations_per_1k_strand_specific'] = per_1k_strand_specific(total_mutations, weighted_normalizer)
            else:
                row['total_mutations_per_1k_strand_specific'] = 0
            
            # Legacy-style fields for backward compatibility - now use strand-specific normalization
            # Names no longer encode repeat count; that is tracked at the file level.
            row['total_mutations_over_total_g_strand_per_1k'] = per_1k_strand_specific(total_mutations, g_strand_normalizer)
            row['total_mutations_over_total_c_strand_per_1k'] = per_1k_strand_specific(total_mutations, c_strand_normalizer)

            # --- New engineered features ---

            # Composite per strand
            row['g_strand_mutations_sum_per_1k'] = sum(row.get(k, 0) for k in row if k.startswith('g_strand_mutations') and k.endswith('_per_1k'))
            row['c_strand_mutations_sum_per_1k'] = sum(row.get(k, 0) for k in row if k.startswith('c_strand_mutations') and k.endswith('_per_1k'))

            # Log-transformed telomere length
            try:
                row['log_telomere_length'] = np.log1p(float(row['Telomere_Length'])) if row['Telomere_Length'] else 0
            except Exception:
                row['log_telomere_length'] = 0

            # Telomere length bin
            try:
                length_val = float(row['Telomere_Length'])
                if length_val < 5000:
                    row['telomere_length_bin'] = 'short'
                elif length_val < 8000:
                    row['telomere_length_bin'] = 'medium'
                else:
                    row['telomere_length_bin'] = 'long'
            except Exception:
                row['telomere_length_bin'] = 'unknown'

           

            # --- New composite scores for higher correlation ---
            
            # 1. Mutation Rate Normalized by Length: (Sum of G-strand mutations per 1k) / Telomere_Length
            try:
                telomere_length = float(row['Telomere_Length'])
                if telomere_length > 0:
                    row['mutation_rate_normalized_by_length'] = row['g_strand_mutations_sum_per_1k'] / telomere_length
                else:
                    row['mutation_rate_normalized_by_length'] = 0
            except Exception:
                row['mutation_rate_normalized_by_length'] = 0

           


           

            # Write combined row
            writer.writerow(row)
            

            messages = [
                f"\nProcessing {filename}:",
                f"Age: {age}",
                f"Telomere Length: {length}",
                f"Total Reads: {total_reads}",
                f"{patterns_version} c-strand total: {c_strand_total}",
                f"{patterns_version} g-strand total: {g_strand_total}",
                f"G-strand mutations total: {g_strand_mutations_total}",
                f"C-strand mutations total: {c_strand_mutations_total}",
                f"G-strand normalizer ({patterns_version} + mutations): {g_strand_normalizer}",
                f"C-strand normalizer ({patterns_version} + mutations): {c_strand_normalizer}",
                f"Total mutations found: {total_mutations}",
            ]
            
            if counts['g_strand'] == 0 and counts['c_strand'] == 0:
                messages.append(f"Warning: No telomere sequences found in {filename}")
            
            # Print to console and send to callback
            for message in messages:
                print(message)
                if output_callback:
                    output_callback(message)

if __name__ == "__main__":  
    generate_csv(".")
