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

def count_total_reads(file_path: str) -> int:
    """Count the total number of reads in a FASTQ or FASTA file."""
    read_count = 0
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    # Detect format by file extension
    is_fasta = any(ext in file_path.lower() for ext in ['.fasta', '.fa', '.fas'])
    
    with open_func(file_path, mode) as f:
        if is_fasta:
            # FASTA format - count header lines starting with '>'
            for line in f:
                if line.startswith('>'):
                    read_count += 1
        else:
            # FASTQ format - count every 4th line (headers start with '@')
            line_count = 0
            for line in f:
                line_count += 1
                if line_count % 4 == 1:  # Every 4th line starting from 1 is a header
                    read_count += 1
    
    return read_count

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
        'g_strand_frameshifts': {
            # Deletions (various deletion patterns from GGGTTAGGGTTA)
            'del_firstG': "GGTTAGGGTTA",      # Delete first G
            'del_firstT': "GGGTAGGGTTA",      # Delete first T  
            'del_secondT': "GGGTGGGTTA",      # Delete second T (from TTA)
            'del_A': "GGGTTGGGTTA",           # Delete A
            'del_midG': "GGGTTAGGTTA",        # Delete one G from GGG in middle
            'del_lastT': "GGGTTAGGGTA",       # Delete one T from TTA at end
            'del_finalA': "GGGTTAGGGTT",      # Delete final A
            
            # Insertions (13bp variants from 12bp standard)
            'ins_G_pos1': "GGGGTTTAGGGTTA",   # Insert G at position 1
            'ins_G_pos4': "GGGTTGAGGGTTA",    # Insert G after TTA
            'ins_T_pos5': "GGGTTATGGGTTA",    # Insert T after A
            'ins_A_pos6': "GGGTTAAGGGTTA",    # Insert A after first repeat
            'ins_G_pos7': "GGGTTAGGGGTTA",    # Insert G in second repeat
            'ins_T_pos10': "GGGTTAGGGTTTA",   # Insert T in TTA
            'ins_A_pos12': "GGGTTAGGGTTAA",   # Insert A at end
            
            # Compound frameshift: examples like GGGTTAGGGTTA -> GGGTTAGGTTAG
            'fs_delT_insG': "GGGTTAGGTTAG",    # Delete T, insert G (your example)
            'fs_shift_end': "GGGTTAGGGTAG",    # Shift at end: TTA -> TAG  
            'fs_shift_mid': "GGGTTAGAGGTTA",   # Shift in middle
            'fs_dup_partial': "GGGTTAGGGGGTTTA", # Partial duplication
            
            # Other common frameshift patterns
            'fs_GGTTA_to_GGTAG': "GGGTTAGGTAG",  # TTA->TAG shift
            'fs_slip_rep1': "GGGTTGGGTTA",       # Slippage in first repeat
            'fs_slip_rep2': "GGGTTAGGTTA",       # Slippage in second repeat
        },
        'c_strand_frameshifts': {
            # Deletions for C-strand (CCCTAACCCTAA -> variants)
            'del_firstC': "CCTAACCCTAA",      # Delete first C
            'del_T': "CCCAACCCTAA",           # Delete T
            'del_firstA': "CCCTACCCTAA",      # Delete first A
            'del_secondA': "CCCTACCTAA",      # Delete second A (from AA in middle)
            'del_midC': "CCCTAACCTAA",        # Delete one C from CCC in middle
            'del_lastT': "CCCTAACCCAA",       # Delete T in second repeat  
            'del_finalA': "CCCTAACCCTA",      # Delete final A
            
            # Insertions for C-strand
            'ins_C_pos1': "CCCCTAACCCTAA",    # Insert C at start
            'ins_T_pos4': "CCCTAACCCTAA",     # Insert T after CCT
            'ins_A_pos5': "CCCTAAACCCTAA",    # Insert A after TA
            'ins_A_pos6': "CCCTAAACCCTAA",    # Insert A after first repeat
            'ins_C_pos7': "CCCTAACCCCTAA",    # Insert C in second repeat
            'ins_T_pos10': "CCCTAACCCTTAA",   # Insert T in second repeat
            'ins_A_pos11': "CCCTAACCCTAAA",   # Insert A at end
            
            # Compound frameshifts for C-strand
            'fs_delA_insT': "CCCTAACCCTATA",   # Similar pattern to G-strand example
            'fs_shift_end': "CCCTAACCCTAG",    # AA -> AG shift
            'fs_shift_mid': "CCCTAACGCTAA",    # Shift in middle
            'fs_dup_partial': "CCCTAACCCCTAAA", # Partial duplication
            
            # Other C-strand frameshift patterns
            'fs_TAA_to_TAG': "CCCTAACCCTAG",   # TAA->TAG shift
            'fs_slip_rep1': "CCCTACCCTAA",     # Slippage in first repeat
            'fs_slip_rep2': "CCCTAACCTAA",     # Slippage in second repeat
        },
    }
    
    with open('telomere_analysis.csv', 'w', newline='') as csvfile:
        fieldnames = [
            'FileName', 'Age', 'Telomere_Length', 'Total_Reads',
            'c_strand', 'g_strand',
        ]
        mutation_keys = []
        for group in ['g_strand_mutations', 'c_strand_mutations', 'g_strand_frameshifts', 'c_strand_frameshifts']:
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
                # Frameshift categories
                'deletions': ['del_firstG', 'del_firstT', 'del_secondT', 'del_A', 'del_midG', 'del_lastT', 'del_finalA'],
                'insertions': ['ins_G_pos1', 'ins_G_pos4', 'ins_T_pos5', 'ins_A_pos6', 
                              'ins_G_pos7', 'ins_T_pos10', 'ins_A_pos12'],
                'compound_fs': ['fs_delT_insG', 'fs_shift_end', 'fs_shift_mid', 'fs_dup_partial',
                               'fs_GGTTA_to_GGTAG', 'fs_slip_rep1', 'fs_slip_rep2'],
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
                # Frameshift categories
                'deletions': ['del_firstC', 'del_T', 'del_firstA', 'del_secondA', 'del_midC', 'del_lastT', 'del_finalA'],
                'insertions': ['ins_C_pos1', 'ins_T_pos4', 'ins_A_pos5', 'ins_A_pos6',
                              'ins_C_pos7', 'ins_T_pos10', 'ins_A_pos11'],
                'compound_fs': ['fs_delA_insT', 'fs_shift_end', 'fs_shift_mid', 'fs_dup_partial',
                               'fs_TAA_to_TAG', 'fs_slip_rep1', 'fs_slip_rep2'],
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
            'telomere_transition_interaction',
            'mutation_rate_normalized_by_length',
            'log_telomere_tg_composite',
            'ratio_ga_g3_ct_c3',
            'composite_score'
        ])

        # Add summed per-1k columns for each mutation type per strand
        summed_per_1k_headers = []
        for strand, mutmap in general_mutation_map.items():
            for mut, subtypes in mutmap.items():
                summed_per_1k_headers.append(f"{strand}_{mut}_sum_per_1k")
        fieldnames.extend(summed_per_1k_headers)

        # Add total mutations field at the end (now read-based)
        fieldnames.append('total_mutations_per_1k_reads')

        # Add legacy field names for backward compatibility 
        fieldnames.extend([
            'total_mutations_over_total_g_strand_2xrepeats_per_1k',
            'total_mutations_over_total_c_strand_2xrepeats_per_1k'
        ])
        
        # Add frameshift-specific fields
        fieldnames.extend([
            'total_g_strand_frameshifts_per_1k',
            'total_c_strand_frameshifts_per_1k',
            'total_all_frameshifts_per_1k',
            'total_g_strand_deletions_per_1k',
            'total_c_strand_deletions_per_1k', 
            'total_g_strand_insertions_per_1k',
            'total_c_strand_insertions_per_1k'
        ])

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
                # Count all mutation sub-patterns including frameshifts
                for group in ['g_strand_mutations', 'c_strand_mutations', 'g_strand_frameshifts', 'c_strand_frameshifts']:
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

            # Calculate per 1k rates for each mutation type (now based on total reads)
            def per_1k_reads(val, total_reads):
                return (val / total_reads) * 1000 if total_reads > 0 else 0
            
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
                # Use read-based normalization for all mutations
                row[f"{k}_per_1k"] = per_1k_reads(counts.get(k, 0), total_reads)

            # --- General mutation per_1k columns ---
            for strand, mutmap in general_mutation_map.items():
                for mut, subtypes in mutmap.items():
                    # Handle both regular mutations and frameshifts
                    if mut in ['deletions', 'insertions', 'compound_fs']:
                        # Frameshift mutations
                        total = sum(counts.get(f"{strand}_frameshifts_{subtype}", 0) for subtype in subtypes)
                    else:
                        # Regular mutations
                        total = sum(counts.get(f"{strand}_mutations_{subtype}", 0) for subtype in subtypes)
                    # Use read-based normalization
                    row[f"{strand}_{mut}_per_1k"] = per_1k_reads(total, total_reads)

            # --- Summed per-1k columns for each mutation type per strand ---
            for strand, mutmap in general_mutation_map.items():
                for mut, subtypes in mutmap.items():
                    # Handle both regular mutations and frameshifts
                    if mut in ['deletions', 'insertions', 'compound_fs']:
                        # Frameshift mutations
                        per_1k_sum = sum(row.get(f"{strand}_frameshifts_{subtype}_per_1k", 0) for subtype in subtypes)
                    else:
                        # Regular mutations
                        per_1k_sum = sum(row.get(f"{strand}_mutations_{subtype}_per_1k", 0) for subtype in subtypes)
                    row[f"{strand}_{mut}_sum_per_1k"] = per_1k_sum
            
            # Total mutations (sum all mutation counts) - now normalized by total reads
            total_mutations = sum(counts[k] for k in mutation_keys)
            row['total_mutations_per_1k_reads'] = per_1k_reads(total_mutations, total_reads)
            
            # Legacy fields for backward compatibility
            row['total_mutations_over_total_g_strand_2xrepeats_per_1k'] = per_1k_reads(total_mutations, total_reads)
            row['total_mutations_over_total_c_strand_2xrepeats_per_1k'] = per_1k_reads(total_mutations, total_reads)

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

            # Composite per strand (including frameshifts)
            row['g_strand_mutations_sum_per_1k'] = sum(row.get(k, 0) for k in row if (k.startswith('g_strand_mutations') or k.startswith('g_strand_frameshifts')) and k.endswith('_per_1k'))
            row['c_strand_mutations_sum_per_1k'] = sum(row.get(k, 0) for k in row if (k.startswith('c_strand_mutations') or k.startswith('c_strand_frameshifts')) and k.endswith('_per_1k'))

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

            # 3. Ratio of Specific Mutations: (G-strand mutations G>A g3 per 1k) / (C-strand mutations C>T c3 per 1k)
            try:
                ga_g3_per_1k = row.get('g_strand_mutations_G>A_g3_per_1k', 0)
                ct_c3_per_1k = row.get('c_strand_mutations_C>T_c3_per_1k', 0)
                if ct_c3_per_1k > 0:
                    row['ratio_ga_g3_ct_c3'] = ga_g3_per_1k / ct_c3_per_1k
                else:
                    row['ratio_ga_g3_ct_c3'] = 0
            except Exception:
                row['ratio_ga_g3_ct_c3'] = 0

            # 4. Linear Combination (Weighted Sum): -0.6 * Telomere_Length + 0.4 * (Total mutations per 1k reads)
            try:
                telomere_length = float(row['Telomere_Length'])
                total_mutations_per_1k = row.get('total_mutations_per_1k_reads', 0)
                row['composite_score'] = -0.6 * telomere_length + 0.4 * total_mutations_per_1k
            except Exception:
                row['composite_score'] = 0

            # --- Calculate frameshift totals ---
            
            # Total frameshifts per strand
            row['total_g_strand_frameshifts_per_1k'] = sum(row.get(k, 0) for k in row if k.startswith('g_strand_frameshifts') and k.endswith('_per_1k'))
            row['total_c_strand_frameshifts_per_1k'] = sum(row.get(k, 0) for k in row if k.startswith('c_strand_frameshifts') and k.endswith('_per_1k'))
            row['total_all_frameshifts_per_1k'] = row['total_g_strand_frameshifts_per_1k'] + row['total_c_strand_frameshifts_per_1k']
            
            # Deletion totals
            g_deletions = ['g_strand_frameshifts_del_firstG_per_1k', 'g_strand_frameshifts_del_firstT_per_1k', 
                          'g_strand_frameshifts_del_secondT_per_1k', 'g_strand_frameshifts_del_A_per_1k', 
                          'g_strand_frameshifts_del_midG_per_1k', 'g_strand_frameshifts_del_lastT_per_1k', 
                          'g_strand_frameshifts_del_finalA_per_1k']
            row['total_g_strand_deletions_per_1k'] = sum(row.get(k, 0) for k in g_deletions)
            
            c_deletions = ['c_strand_frameshifts_del_firstC_per_1k', 'c_strand_frameshifts_del_T_per_1k',
                          'c_strand_frameshifts_del_firstA_per_1k', 'c_strand_frameshifts_del_secondA_per_1k',
                          'c_strand_frameshifts_del_midC_per_1k', 'c_strand_frameshifts_del_lastT_per_1k',
                          'c_strand_frameshifts_del_finalA_per_1k']
            row['total_c_strand_deletions_per_1k'] = sum(row.get(k, 0) for k in c_deletions)
            
            # Insertion totals
            g_insertions = ['g_strand_frameshifts_ins_G_pos1_per_1k', 'g_strand_frameshifts_ins_G_pos4_per_1k',
                           'g_strand_frameshifts_ins_T_pos5_per_1k', 'g_strand_frameshifts_ins_A_pos6_per_1k',
                           'g_strand_frameshifts_ins_G_pos7_per_1k', 'g_strand_frameshifts_ins_T_pos10_per_1k',
                           'g_strand_frameshifts_ins_A_pos12_per_1k']
            row['total_g_strand_insertions_per_1k'] = sum(row.get(k, 0) for k in g_insertions)
            
            c_insertions = ['c_strand_frameshifts_ins_C_pos1_per_1k', 'c_strand_frameshifts_ins_T_pos4_per_1k',
                           'c_strand_frameshifts_ins_A_pos5_per_1k', 'c_strand_frameshifts_ins_A_pos6_per_1k',
                           'c_strand_frameshifts_ins_C_pos7_per_1k', 'c_strand_frameshifts_ins_T_pos10_per_1k',
                           'c_strand_frameshifts_ins_A_pos11_per_1k']
            row['total_c_strand_insertions_per_1k'] = sum(row.get(k, 0) for k in c_insertions)

            writer.writerow(row)
            
            print(f"\nProcessing {filename}:")
            print(f"Age: {age}")
            print(f"Telomere Length: {length}")
            print(f"Total Reads: {total_reads}")
            print(f"2x c-strand total: {c_strand_total}")
            print(f"2x g-strand total: {g_strand_total}")
            print(f"Total mutations found: {total_mutations}")
            print(f"Frameshifts - G-strand: {row.get('total_g_strand_frameshifts_per_1k', 0):.2f} per 1k reads")
            print(f"Frameshifts - C-strand: {row.get('total_c_strand_frameshifts_per_1k', 0):.2f} per 1k reads")
            if counts['g_strand'] == 0 and counts['c_strand'] == 0:
                print(f"Warning: No telomere sequences found in {filename}")

if __name__ == "__main__":
    data_directory = "../data/tesla_minion1"  # Default data directory
    generate_csv(data_directory)
