#!/usr/bin/env python3

import gzip
from collections import defaultdict
import csv
import os
import glob

def read_fastq(file_path: str):
    """Read FASTQ file and yield sequences."""
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    with open_func(file_path, mode) as f:
        while True:
            header = f.readline().strip()
            if not header:  # End of file
                break
            sequence = f.readline().strip()
            _ = f.readline()  
            _ = f.readline() 
            yield sequence

def get_reverse_complement(sequence: str) -> str:
    """Get the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(sequence))

def count_patterns(sequence: str, pattern: str) -> int:
    """Count occurrences of a pattern in a sequence."""
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

def get_fastq_files(directory: str):
    """Get all FASTQ files in the given directory."""
    # Look for both .fastq and .fastq.gz files
    fastq_files = glob.glob(os.path.join(directory, "*.fastq"))
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fastq.gz")))
    return sorted(fastq_files)  # Sort for consistent ordering

def main():
    fastq_files = get_fastq_files("greider_data_download")
    
    if not fastq_files:
        print("No FASTQ files found in greider_data directory")
        return
    
    # Load age data
    age_data = load_age_data()
    length_data = load_length_data()
    
    # Expanded patterns for all 6 mutation types at 3 positions each (G-strand context)
    patterns = {
        'c_strand': "CCCTAACCCTAA",  # C-strand pattern
        'g_strand': "GGGTTAGGGTTA",  # G-strand pattern
        # G-strand mutations (positions 1, 2, 3)
        'G>A_g1': "GGGTTAAGGTTA",
        'G>A_g2': "GGGTTAGAGTTA",
        'G>A_g3': "GGGTTAGGATTA",
        'G>C_g1': "GGGTTACGGTTA",
        'G>C_g2': "GGGTTAGCGTTA",
        'G>C_g3': "GGGTTAGGCTTA",
        'G>T_g1': "GGGTTATGGTTA",
        'G>T_g2': "GGGTTAGTGTTA",
        'G>T_g3': "GGGTTAGGTTTA",
        # C-strand mutations (positions 1, 2, 3)
        'C>A_c1': "CCATAACCCTAA",
        'C>A_c2': "CCCTAAACCTAA",
        'C>A_c3': "CCCTAACCATAA",
        'C>G_c1': "CCGTAACCCTAA",
        'C>G_c2': "CCCTAGGCCTAA",
        'C>G_c3': "CCCTAACCCTGA",
        'C>T_c1': "CCTTAACCCTAA",
        'C>T_c2': "CCCTATTCCTAA",
        'C>T_c3': "CCCTAACCCTTA",
        # T-strand mutations (positions 1, 2, 3)
        'T>A_t1': "GGGTAAGGGTTA",
        'T>A_t2': "GGGTTAAAGTTA",
        'T>A_t3': "GGGTTAGGGAAA",
        'T>C_t1': "GGGTCAGGGTTA",
        'T>C_t2': "GGGTTACGGTTA",
        'T>C_t3': "GGGTTAGGGTCA",
        'T>G_t1': "GGGTGAGGGTTA",
        'T>G_t2': "GGGTTAGGGGTA",
        'T>G_t3': "GGGTTAGGGGTA",
    }
    
    # Create CSV file
    with open('telomere_analysis.csv', 'w', newline='') as csvfile:
        # Define CSV headers
        fieldnames = [
            'FileName', 'Age', 'Telomere_Length', '2x_cstrand', '2xg_strand',
            # G-strand mutation counts
            'G_A_g1', 'G_A_g2', 'G_A_g3',
            'G_C_g1', 'G_C_g2', 'G_C_g3',
            'G_T_g1', 'G_T_g2', 'G_T_g3',
            # C-strand mutation counts
            'C_A_c1', 'C_A_c2', 'C_A_c3',
            'C_G_c1', 'C_G_c2', 'C_G_c3',
            'C_T_c1', 'C_T_c2', 'C_T_c3',
            # T-strand mutation counts
            'T_A_t1', 'T_A_t2', 'T_A_t3',
            'T_C_t1', 'T_C_t2', 'T_C_t3',
            'T_G_t1', 'T_G_t2', 'T_G_t3',
            # Normalized rates (per 1k g_strand)
            'G_A_g1_per_1k', 'G_A_g2_per_1k', 'G_A_g3_per_1k',
            'G_C_g1_per_1k', 'G_C_g2_per_1k', 'G_C_g3_per_1k',
            'G_T_g1_per_1k', 'G_T_g2_per_1k', 'G_T_g3_per_1k',
            # Normalized rates (per 1k c_strand)
            'C_A_c1_per_1k', 'C_A_c2_per_1k', 'C_A_c3_per_1k',
            'C_G_c1_per_1k', 'C_G_c2_per_1k', 'C_G_c3_per_1k',
            'C_T_c1_per_1k', 'C_T_c2_per_1k', 'C_T_c3_per_1k',
            # Normalized rates (per 1k t_strand, if needed)
            'T_A_t1_per_1k', 'T_A_t2_per_1k', 'T_A_t3_per_1k',
            'T_C_t1_per_1k', 'T_C_t2_per_1k', 'T_C_t3_per_1k',
            'T_G_t1_per_1k', 'T_G_t2_per_1k', 'T_G_t3_per_1k',
            # Totals
            'total_mutations_over_total_g_per_1k'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for file_path in fastq_files:
            counts = defaultdict(int)
            
            for sequence in read_fastq(file_path):
                # Count c-strand in forward direction only
                counts['c_strand'] += count_patterns(sequence, patterns['c_strand'])
                # Count g-strand in forward direction only
                counts['g_strand'] += count_patterns(sequence, patterns['g_strand'])
                # Count all mutation patterns
                for name, pattern in patterns.items():
                    if name not in ['c_strand', 'g_strand']:
                        counts[name] += count_patterns(sequence, pattern)
            
            filename = os.path.basename(file_path)
            filename_base = filename.replace('.fastq', '').replace('.gz', '')
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
                '2x_cstrand': c_strand_total,
                '2xg_strand': g_strand_total,
            }
            # G-strand
            for mut in ['G>A', 'G>C', 'G>T']:
                for i in range(1, 4):
                    key = f'{mut}_g{i}'
                    row[f'{mut.replace(">", "_")}_g{i}'] = counts.get(key, 0)
                    row[f'{mut.replace(">", "_")}_g{i}_per_1k'] = per_1k(counts.get(key, 0), g_strand_total)
            # C-strand
            for mut in ['C>A', 'C>G', 'C>T']:
                for i in range(1, 4):
                    key = f'{mut}_c{i}'
                    row[f'{mut.replace(">", "_")}_c{i}'] = counts.get(key, 0)
                    row[f'{mut.replace(">", "_")}_c{i}_per_1k'] = per_1k(counts.get(key, 0), c_strand_total)
            # T-strand (if needed, else can be omitted)
            for mut in ['T>A', 'T>C', 'T>G']:
                for i in range(1, 4):
                    key = f'{mut}_t{i}'
                    row[f'{mut.replace(">", "_")}_t{i}'] = counts.get(key, 0)
                    row[f'{mut.replace(">", "_")}_t{i}_per_1k'] = per_1k(counts.get(key, 0), g_strand_total)  # or t_strand_total if available
            
            # Total mutations (sum all mutation counts)
            total_mutations = sum(counts[k] for k in counts if k not in ['c_strand', 'g_strand'])
            row['total_mutations_over_total_g_strand_2xrepeats_per_1k'] = per_1k(total_mutations, g_strand_total)
            writer.writerow(row)
            
            # Print to console as well
            print(f"\nProcessing {filename}:")
            print(f"Age: {age}")
            print(f"Telomere Length: {length}")
            print(f"2x cstrand total: {c_strand_total}")
            print(f"2x g strand total: {g_strand_total}")

if __name__ == "__main__":
    main()
