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
            _ = f.readline()  # Skip the '+' line
            _ = f.readline()  # Skip the quality line
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
            # Convert the filename format to match our FASTQ files
            # e.g., "JH47.F86_NB70" -> "JH47.F86.NB70"
            fastq_name = row['fastq file name'].replace('_', '.')
            age_data[fastq_name] = row['Age (Years)']
    return age_data

def get_fastq_files(directory: str):
    """Get all FASTQ files in the given directory."""
    # Look for both .fastq and .fastq.gz files
    fastq_files = glob.glob(os.path.join(directory, "*.fastq"))
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fastq.gz")))
    return sorted(fastq_files)  # Sort for consistent ordering

def main():
    # Get all FASTQ files from the greider_data directory
    fastq_files = get_fastq_files("greider_data")
    
    if not fastq_files:
        print("No FASTQ files found in greider_data directory")
        return
    
    # Load age data
    age_data = load_age_data()
    
    patterns = {
        'c_strand': "CCCTAACCCTAA",  # C-strand pattern
        'g_strand': "GGGTTAGGGTTA",  # G-strand pattern
        'G>T_g1': "GGGTTATGGTTA",
        'G>T_g2': "GGGTTAGTGTTA",
        'G>T_g3': "GGGTTAGGTTTA",
        'G>A_g1': "GGGTTAAGGTTA",
        'G>A_g2': "GGGTTAGAGTTA",
        'G>A_g3': "GGGTTAGGATTA"
    }
    
    # Create CSV file
    with open('telomere_analysis.csv', 'w', newline='') as csvfile:
        # Define CSV headers
        fieldnames = [
            'FileName', 'Age', '2x_cstrand', '2xg_strand', 
            'G_T_g1', 'G_T_g2', 'G_T_g3', 
            'G_A_g1', 'G_A_g2', 'G_A_g3',
            'G_T_g1_over_total_g_per_1k',
            'G_T_g2_over_total_g_per_1k',
            'G_T_g3_over_total_g_per_1k',
            'G_A_g1_over_total_g_per_1k',
            'G_A_g2_over_total_g_per_1k',
            'G_A_g3_over_total_g_per_1k',
            'pos1_mutation_rate_per_1k',
            'pos2_mutation_rate_per_1k',
            'pos3_mutation_rate_per_1k',
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
                
                # Count mutations in both directions
                for name, pattern in patterns.items():
                    if name not in ['c_strand', 'g_strand']:
                        counts[name] += count_patterns(sequence, pattern)
            
            # Get just the filename from the path
            filename = os.path.basename(file_path)
            # Remove .fastq or .fastq.gz extension for matching
            filename_base = filename.replace('.fastq', '').replace('.gz', '')
            
            # Get age from the data table
            age = age_data.get(filename_base, '')
            
            # Calculate per 1k rates
            g_strand_total = counts['g_strand']
            if g_strand_total > 0:  # Avoid division by zero
                g_t_g1_rate = (counts['G>T_g1'] / g_strand_total) * 1000
                g_t_g2_rate = (counts['G>T_g2'] / g_strand_total) * 1000
                g_t_g3_rate = (counts['G>T_g3'] / g_strand_total) * 1000
                
                g_a_g1_rate = (counts['G>A_g1'] / g_strand_total) * 1000
                g_a_g2_rate = (counts['G>A_g2'] / g_strand_total) * 1000
                g_a_g3_rate = (counts['G>A_g3'] / g_strand_total) * 1000
                
                pos1_rate = ((counts['G>A_g1'] + counts['G>T_g1']) / g_strand_total) * 1000
                pos2_rate = ((counts['G>A_g2'] + counts['G>T_g2']) / g_strand_total) * 1000
                pos3_rate = ((counts['G>A_g3'] + counts['G>T_g3']) / g_strand_total) * 1000
                
                total_mutations = sum(counts[p] for p in patterns if p not in ['c_strand', 'g_strand'])
                total_mutation_rate = (total_mutations / g_strand_total) * 1000
            else:
                g_t_g1_rate = g_t_g2_rate = g_t_g3_rate = 0
                g_a_g1_rate = g_a_g2_rate = g_a_g3_rate = 0
                pos1_rate = pos2_rate = pos3_rate = 0
                total_mutation_rate = 0
            
            # Write row to CSV
            writer.writerow({
                'FileName': filename,
                'Age': age,
                '2x_cstrand': counts['c_strand'],
                '2xg_strand': counts['g_strand'],
                'G_T_g1': counts['G>T_g1'],
                'G_T_g2': counts['G>T_g2'],
                'G_T_g3': counts['G>T_g3'],
                'G_A_g1': counts['G>A_g1'],
                'G_A_g2': counts['G>A_g2'],
                'G_A_g3': counts['G>A_g3'],
                'G_T_g1_over_total_g_per_1k': g_t_g1_rate,
                'G_T_g2_over_total_g_per_1k': g_t_g2_rate,
                'G_T_g3_over_total_g_per_1k': g_t_g3_rate,
                'G_A_g1_over_total_g_per_1k': g_a_g1_rate,
                'G_A_g2_over_total_g_per_1k': g_a_g2_rate,
                'G_A_g3_over_total_g_per_1k': g_a_g3_rate,
                'pos1_mutation_rate_per_1k': pos1_rate,
                'pos2_mutation_rate_per_1k': pos2_rate,
                'pos3_mutation_rate_per_1k': pos3_rate,
                'total_mutations_over_total_g_per_1k': total_mutation_rate
            })
            
            # Print to console as well
            print(f"\nProcessing {filename}:")
            print(f"Age: {age}")
            print(f"2x cstrand total: {counts['c_strand']}")
            print(f"From strand:")
            print(f"  c1: {counts['c_strand']}")
            print(f"  c2: {counts['c_strand']}")
            print(f"  c3: {counts['c_strand']}")
            
            print(f"\n2x g strand total: {counts['g_strand']}")
            
            print("\nG>T:")
            print(f"  g1: {counts['G>T_g1']}")
            print(f"  g2: {counts['G>T_g2']}")
            print(f"  g3: {counts['G>T_g3']}")
            
            print("\nG>A:")
            print(f"  g1: {counts['G>A_g1']}")
            print(f"  g2: {counts['G>A_g2']}")
            print(f"  g3: {counts['G>A_g3']}")

if __name__ == "__main__":
    main()
