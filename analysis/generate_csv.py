#!/usr/bin/env python3

import gzip
from collections import defaultdict
import csv
import os
import numpy as np
import glob
import HTSeq

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
        # 'g_strand_frameshifts': {
        #     'del_firstG': "GGGTTAGGTTAG",
        #     'del_firstT': "GGGTTAGGGTAG",
        #     'del_secondT': "GGGTTAGGGTAG",
        #     'del_A': "GGGTTAGGGTTG",
        #     'ins_G_pos1': "GGGTTAGGGGTTA",
        #     'ins_G_pos4': "GGGTTAGGGTTA",
        #     'ins_T_pos5': "GGGTTATGGGTTA",
        #     'ins_A_pos6': "GGGTTAAGGGTTA",
        #     'ins_G_pos7': "GGGTTAGGGGTTA",
        #     'ins_T_pos10': "GGGTTAGGGTTTA",
        #     'ins_A_pos12': "GGGTTAGGGTTAA",
        #     'fs_GGTTA_to_GGTAG': "GGGTTAGGTAG",
        #     'fs_slip_rep1': "GGGTTGGGTTA",
        #     'fs_slip_rep2': "GGGTTAGGTTA",
        # },
        # 'c_strand_frameshifts': {
        #     'del_firstC': "ACCTAACCCTAA",
        #     'del_T': "CCCAACCCTAA",
        #     'del_firstA': "CCCTACCCTAA",
        #     'del_secondA': "CCCTACCTAA",
        #     'del_midC': "CCCTAACCTAA",
        #     'del_lastT': "CCCTAACCCAA",
        #     'del_finalA': "CCCTAACCCTA",
        #     'ins_C_pos1': "CCCCTAACCCTAA",
        #     'ins_T_pos4': "CCCTAACCCTAA",
        #     'ins_A_pos5': "CCCTAAACCCTAA",
        #     'ins_A_pos6': "CCCTAAACCCTAA",
        #     'ins_C_pos7': "CCCTAACCCCTAA",
        #     'ins_T_pos10': "CCCTAACCCTTAA",
        #     'ins_A_pos11': "CCCTAACCCTAAA",
        #     'fs_delA_insT': "CCCTAACCCTATA",
        #     'fs_shift_end': "CCCTAACCCTAG",
        #     'fs_shift_mid': "CCCTAACGCTAA",
        #     'fs_dup_partial': "CCCTAACCCCTAAA",
        #     'fs_TAA_to_TAG': "CCCTAACCCTAG",
        #     'fs_slip_rep1': "CCCTACCCTAA",
        #     'fs_slip_rep2': "CCCTAACCTAA",
        # },
    }
    
    with open('telomere_analysis.csv', 'w', newline='') as csvfile:
        fieldnames = [
            'FileName', 'Age', 'Telomere_Length', 'Total_Reads',
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
            # 'mutation_ratio_CtoA_CtoT',
            # 'mutation_ratio_GtoA_GtoT',
            'g_strand_mutations_sum_per_1k',
            'c_strand_mutations_sum_per_1k',
            'log_telomere_length',
            'telomere_length_bin',
            'telomere_transition_interaction',
            'mutation_rate_normalized_by_length',
            'log_telomere_tg_composite',
            # 'ratio_ga_g3_ct_c3',
            'composite_score'
        ])

        # Add summed per-1k columns for each mutation type per strand
        summed_per_1k_headers = []
        for strand, mutmap in general_mutation_map.items():
            for mut, subtypes in mutmap.items():
                summed_per_1k_headers.append(f"{strand}_{mut}_sum_per_1k")
        fieldnames.extend(summed_per_1k_headers)

        # Add total mutations field at the end (now strand-specific normalized)
        fieldnames.append('total_mutations_per_1k_strand_specific')

        # Add legacy field names for backward compatibility 
        fieldnames.extend([
            'total_mutations_over_total_g_strand_2xrepeats_per_1k',
            'total_mutations_over_total_c_strand_2xrepeats_per_1k'
        ])
        
        # Add frameshift-specific fields (disabled)
        # fieldnames.extend([
        #     'total_g_strand_frameshifts_per_1k',
        #     'total_c_strand_frameshifts_per_1k',
        #     'total_all_frameshifts_per_1k',
        #     'total_g_strand_deletions_per_1k',
        #     'total_c_strand_deletions_per_1k', 
        #     'total_g_strand_insertions_per_1k',
        #     'total_c_strand_insertions_per_1k'
        # ])

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
                    per_1k_sum = sum(row.get(f"{strand}_mutations_{subtype}_per_1k", 0) for subtype in subtypes)
                    row[f"{strand}_{mut}_sum_per_1k"] = per_1k_sum
            
            # Total mutations (sum all mutation counts) - now normalized by strand-specific normalizers
            total_mutations = sum(counts[k] for k in mutation_keys)
            # Use weighted average of strand-specific normalizations
            if g_strand_normalizer > 0 or c_strand_normalizer > 0:
                g_weight = g_strand_mutations_total / total_mutations if total_mutations > 0 else 0
                c_weight = c_strand_mutations_total / total_mutations if total_mutations > 0 else 0
                weighted_normalizer = (g_weight * g_strand_normalizer) + (c_weight * c_strand_normalizer)
                row['total_mutations_per_1k_strand_specific'] = per_1k_strand_specific(total_mutations, weighted_normalizer)
            else:
                row['total_mutations_per_1k_strand_specific'] = 0
            
            # Legacy fields for backward compatibility - now use strand-specific normalization
            row['total_mutations_over_total_g_strand_2xrepeats_per_1k'] = per_1k_strand_specific(total_mutations, g_strand_normalizer)
            row['total_mutations_over_total_c_strand_2xrepeats_per_1k'] = per_1k_strand_specific(total_mutations, c_strand_normalizer)

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


            # Composite per strand (frameshifts disabled)
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

            # 2. Log-Transformed Composite: log(Telomere_Length) * (G-strand mutations T>G t1 per 1k)
            try:
                telomere_length = float(row['Telomere_Length'])
                tg_t1_per_1k = row.get('g_strand_mutations_T>G_t1_per_1k', 0)
                if telomere_length > 0:
                    row['log_telomere_tg_composite'] = np.log(telomere_length) * tg_t1_per_1k
                else:
                    row['log_telomere_tg_composite'] = 0
            except Exception:
                row['log_telomere_tg_composite'] = 0


            # 4. Linear Combination (Weighted Sum): -0.6 * Telomere_Length + 0.4 * (Total mutations per 1k strand-specific)
            try:
                telomere_length = float(row['Telomere_Length'])
                total_mutations_per_1k = row.get('total_mutations_per_1k_strand_specific', 0)
                row['composite_score'] = -0.6 * telomere_length + 0.4 * total_mutations_per_1k
            except Exception:
                row['composite_score'] = 0

            writer.writerow(row)
            
            print(f"\nProcessing {filename}:")
            print(f"Age: {age}")
            print(f"Telomere Length: {length}")
            print(f"Total Reads: {total_reads}")
            print(f"2x c-strand total: {c_strand_total}")
            print(f"2x g-strand total: {g_strand_total}")
            print(f"G-strand mutations total: {g_strand_mutations_total}")
            print(f"C-strand mutations total: {c_strand_mutations_total}")
            print(f"G-strand normalizer (2x repeats + mutations): {g_strand_normalizer}")
            print(f"C-strand normalizer (2x repeats + mutations): {c_strand_normalizer}")
            print(f"Total mutations found: {total_mutations}")
            # print(f"Frameshifts - G-strand: {row.get('total_g_strand_frameshifts_per_1k', 0):.2f} per 1k reads")
            # print(f"Frameshifts - C-strand: {row.get('total_c_strand_frameshifts_per_1k', 0):.2f} per 1k reads")
            if counts['g_strand'] == 0 and counts['c_strand'] == 0:
                print(f"Warning: No telomere sequences found in {filename}")

if __name__ == "__main__":  
    generate_csv()
