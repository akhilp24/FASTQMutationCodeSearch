#!/usr/bin/env python3

import gzip
from collections import defaultdict
import csv
import os
import numpy as np
import glob

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

def load_age_data():
    """Load age data from greider_methods_table_s2.csv."""
    age_data = {}
    with open('greider_methods_table_s2.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            fastq_name = row['fastq file name'].replace('_', '.')
            age_data[fastq_name] = row['Age (Years)']
    return age_data

def load_length_data():
    """Load length data from greider_methods_table_s2.csv."""
    length_data = {}
    with open('greider_methods_table_s2.csv', 'r') as f:
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

def generate_csv(data_dir: str):
    sequence_files = get_sequence_files(data_dir)
    
    if not sequence_files:
        print(f"No FASTQ or FASTA files found in {data_dir} directory")
        return
    
    # Load age data
    age_data = load_age_data()
    length_data = load_length_data()
    
    patterns = {
        'c_strand': "CCCTAACCCTAA",
        'g_strand': "GGGTTAGGGTTA",
        'g_strand_mutations': {
            'G>A_g1': "GGGTTAAGGTTA",
            'G>A_g2': "GGGTTAGAGTTA",
            'G>A_g3': "GGGTTAGGATTA",
            'G>C_g1': "GGGTTACGGTTA",
            'G>C_g2': "GGGTTAGCGTTA",
            'G>C_g3': "GGGTTAGGCTTA",
            'G>T_g1': "GGGTTATGGTTA",
            'G>T_g2': "GGGTTAGTGTTA",
            'G>T_g3': "GGGTTAGGTTTA",
            'T>A_t1': "GGGTTAGGGATA",
            'T>A_t2': "GGGTTAGGGTAA",
            'T>C_t1': "GGGTCAGGGTTA",
            'T>C_t2': "GGGTTACGGTTA",
            'T>G_t1': "GGGTGAGGGTTA",
            'T>G_t2': "GGGTTAGGGGTA",
            'A>T_a1': "GGGTTAGGGTTT",
            'A>G_a1': "GGGTTAGGGTTG",
            'A>C_a1': "GGGTTAGGGTTC",
        },
        'c_strand_mutations': {
            'C>A_c1': "CCCTAAACCTAA",
            'C>A_c2': "CCCTAACACTAA",
            'C>A_c3': "CCCTAACCATAA",
            'C>G_c1': "CCCTAAGCCTAA",
            'C>G_c2': "CCCTAACGCTAA",
            'C>G_c3': "CCCTAACCGTAA",
            'C>T_c1': "CCCTAATCCTAA",
            'C>T_c2': "CCCTAACTCTAA",
            'C>T_c3': "CCCTAACCTTAA",
            'T>A_t1': "CCCTAACCCAAA",
            'T>C_t1': "CCCTAACCCCAA",
            'T>G_t1': "CCCTAACCCGAA",
            'A>T_a1': "CCCTAACCCTTA",
            'A>T_a2': "CCCTAACCCTAT",
            'A>G_a1': "CCCTAACCCTGA",
            'A>G_a2': "CCCTAACCCTAG",
            'A>C_a1': "CCCTAACCCTCA",
            'A>C_a2': "CCCTAACCCTAC",
        },
    }
    
    with open('telomere_analysis.csv', 'w', newline='') as csvfile:
        fieldnames = [
            'FileName', 'Age', 'Telomere_Length',
            'c_strand', 'g_strand',
        ]
        mutation_keys = []
        for group in ['g_strand_mutations', 'c_strand_mutations']:
            for subkey in patterns[group].keys():
                mutation_keys.append(f"{group}_{subkey}")
        fieldnames.extend(mutation_keys)

        # Add normalized rate fields for each mutation key
        fieldnames.extend([f"{k}_per_1k" for k in mutation_keys])

        # --- Add general mutation per_1k columns ---
        general_mutation_map = {
            'g_strand': {
                'G>A': ['G>A_g1', 'G>A_g2', 'G>A_g3'],
                'G>C': ['G>C_g1', 'G>C_g2', 'G>C_g3'],
                'G>T': ['G>T_g1', 'G>T_g2', 'G>T_g3'],
                'T>A': ['T>A_t1', 'T>A_t2'],
                'T>C': ['T>C_t1', 'T>C_t2'],
                'T>G': ['T>G_t1', 'T>G_t2'],
                'A>T': ['A>T_a1'],
                'A>G': ['A>G_a1'],
                'A>C': ['A>C_a1'],
            },
            'c_strand': {
                'C>A': ['C>A_c1', 'C>A_c2', 'C>A_c3'],
                'C>G': ['C>G_c1', 'C>G_c2', 'C>G_c3'],
                'C>T': ['C>T_c1', 'C>T_c2', 'C>T_c3'],
                'T>A': ['T>A_t1'],
                'T>C': ['T>C_t1'],
                'T>G': ['T>G_t1'],
                'A>T': ['A>T_a1', 'A>T_a2'],
                'A>G': ['A>G_a1', 'A>G_a2'],
                'A>C': ['A>C_a1', 'A>C_a2'],
            }
        }
        general_mutation_headers = []
        for strand, mutmap in general_mutation_map.items():
            for mut in mutmap:
                general_mutation_headers.append(f"{strand}_{mut}_per_1k")
        fieldnames.extend(general_mutation_headers)

        # Add new engineered features to the CSV header
        fieldnames.extend([
            'composite_transition_per_1k',
            'composite_transversion_per_1k',
            'mutation_ratio_CtoA_CtoT',
            'mutation_ratio_GtoA_GtoT',
            'g_strand_mutations_sum_per_1k',
            'c_strand_mutations_sum_per_1k',
            'log_telomere_length',
            'telomere_length_bin',
            'telomere_transition_interaction'
        ])

        # Add summed per-1k columns for each mutation type per strand
        summed_per_1k_headers = []
        for strand, mutmap in general_mutation_map.items():
            for mut, subtypes in mutmap.items():
                summed_per_1k_headers.append(f"{strand}_{mut}_sum_per_1k")
        fieldnames.extend(summed_per_1k_headers)

        # Add total mutations field at the end
        fieldnames.append('total_mutations_over_total_g_strand_2xrepeats_per_1k')
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for file_path in sequence_files:
            counts = defaultdict(int)
            
            for sequence in read_sequence_file(file_path):  # Changed function name
                # Count c-strand in forward direction only
                counts['c_strand'] += count_patterns(sequence, patterns['c_strand'])
                # Count g-strand in forward direction only
                counts['g_strand'] += count_patterns(sequence, patterns['g_strand'])
                # Count all mutation sub-patterns
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

            # Calculate per 1k rates for each mutation type
            def per_1k(val, total):
                return (val / total) * 1000 if total > 0 else 0
            row = {
                'FileName': filename,
                'Age': age,
                'Telomere_Length': length,
                'c_strand': c_strand_total,
                'g_strand': g_strand_total,
            }

            for k in mutation_keys:
                row[k] = counts.get(k, 0)
                if k.startswith('g_strand_mutations'):
                    norm_total = g_strand_total
                elif k.startswith('c_strand_mutations'):
                    norm_total = c_strand_total
                else:
                    norm_total = g_strand_total  # fallback, should not occur
                row[f"{k}_per_1k"] = per_1k(counts.get(k, 0), norm_total)

            # --- General mutation per_1k columns ---
            for strand, mutmap in general_mutation_map.items():
                strand_total = g_strand_total if strand == 'g_strand' else c_strand_total
                for mut, subtypes in mutmap.items():
                    total = sum(counts.get(f"{strand}_mutations_{subtype}", 0) for subtype in subtypes)
                    row[f"{strand}_{mut}_per_1k"] = per_1k(total, strand_total)

            # --- Summed per-1k columns for each mutation type per strand ---
            for strand, mutmap in general_mutation_map.items():
                for mut, subtypes in mutmap.items():
                    per_1k_sum = sum(row.get(f"{strand}_mutations_{subtype}_per_1k", 0) for subtype in subtypes)
                    row[f"{strand}_{mut}_sum_per_1k"] = per_1k_sum
            
            # Total mutations (sum all mutation counts)
            total_mutations = sum(counts[k] for k in mutation_keys)
            row['total_mutations_over_total_g_strand_2xrepeats_per_1k'] = per_1k(total_mutations, g_strand_total)

            # --- New engineered features ---
           

            # Composite scores (per_1k only)
            transition_keys = [
                'g_strand_mutations_G>A_g1_per_1k', 'g_strand_mutations_G>A_g2_per_1k', 'g_strand_mutations_G>A_g3_per_1k',
                'c_strand_mutations_C>T_c1_per_1k', 'c_strand_mutations_C>T_c2_per_1k', 'c_strand_mutations_C>T_c3_per_1k',
                'g_strand_mutations_A>G_a1_per_1k', 'c_strand_mutations_A>G_a1_per_1k', 'c_strand_mutations_A>G_a2_per_1k',
                'g_strand_mutations_T>C_t1_per_1k', 'g_strand_mutations_T>C_t2_per_1k'
            ]
            transversion_keys = [
                'g_strand_mutations_G>C_g1_per_1k', 'g_strand_mutations_G>C_g2_per_1k', 'g_strand_mutations_G>C_g3_per_1k',
                'g_strand_mutations_G>T_g1_per_1k', 'g_strand_mutations_G>T_g2_per_1k', 'g_strand_mutations_G>T_g3_per_1k',
                'g_strand_mutations_T>A_t1_per_1k', 'g_strand_mutations_T>A_t2_per_1k',
                'g_strand_mutations_T>G_t1_per_1k', 'g_strand_mutations_T>G_t2_per_1k',
                'g_strand_mutations_A>T_a1_per_1k', 'g_strand_mutations_A>C_a1_per_1k',
                'c_strand_mutations_C>A_c1_per_1k', 'c_strand_mutations_C>A_c2_per_1k', 'c_strand_mutations_C>A_c3_per_1k',
                'c_strand_mutations_C>G_c1_per_1k', 'c_strand_mutations_C>G_c2_per_1k', 'c_strand_mutations_C>G_c3_per_1k',
                'c_strand_mutations_T>A_t1_per_1k', 'c_strand_mutations_T>G_t1_per_1k',
                'c_strand_mutations_A>T_a1_per_1k', 'c_strand_mutations_A>T_a2_per_1k',
                'c_strand_mutations_A>C_a1_per_1k', 'c_strand_mutations_A>C_a2_per_1k'
            ]

            row['composite_transition_per_1k'] = sum(row.get(k, 0) for k in transition_keys)
            row['composite_transversion_per_1k'] = sum(row.get(k, 0) for k in transversion_keys)

            # Mutation ratios
            def safe_ratio(a, b):
                return a / b if b != 0 else 0

            row['mutation_ratio_CtoA_CtoT'] = safe_ratio(
                sum(row.get(k, 0) for k in ['c_strand_mutations_C>A_c1_per_1k', 'c_strand_mutations_C>A_c2_per_1k', 'c_strand_mutations_C>A_c3_per_1k']),
                sum(row.get(k, 0) for k in ['c_strand_mutations_C>T_c1_per_1k', 'c_strand_mutations_C>T_c2_per_1k', 'c_strand_mutations_C>T_c3_per_1k']) + 1e-6
            )
            row['mutation_ratio_GtoA_GtoT'] = safe_ratio(
                sum(row.get(k, 0) for k in ['g_strand_mutations_G>A_g1_per_1k', 'g_strand_mutations_G>A_g2_per_1k', 'g_strand_mutations_G>A_g3_per_1k']),
                sum(row.get(k, 0) for k in ['g_strand_mutations_G>T_g1_per_1k', 'g_strand_mutations_G>T_g2_per_1k', 'g_strand_mutations_G>T_g3_per_1k']) + 1e-6
            )

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

            # Interaction term (telomere length × composite_transition_per_1k)
            try:
                row['telomere_transition_interaction'] = float(row['Telomere_Length']) * row['composite_transition_per_1k']
            except Exception:
                row['telomere_transition_interaction'] = 0

            writer.writerow(row)
            
            print(f"\nProcessing {filename}:")
            print(f"Age: {age}")
            print(f"Telomere Length: {length}")
            print(f"2x cstrand total: {c_strand_total}")
            print(f"2x g strand total: {g_strand_total}")
            if counts['g_strand'] == 0 and counts['c_strand'] == 0:
                print(f"Example sequence from {filename}: {sequence}")

if __name__ == "__main__":
    generate_csv()
