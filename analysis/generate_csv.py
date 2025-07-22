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
            if not header:
                break
            sequence = f.readline().strip()
            _ = f.readline()  
            _ = f.readline() 
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

def get_fastq_files(directory: str):
    """Get all FASTQ files in the given directory."""
    # Look for both .fastq and .fastq.gz files
    fastq_files = glob.glob(os.path.join(directory, "*.fastq"))
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fastq.gz")))
    return sorted(fastq_files)  # Sort for consistent ordering

def generate_csv(data_dir: str):
    fastq_files = get_fastq_files(data_dir)
    
    if not fastq_files:
        print(f"No FASTQ files found in {data_dir} directory")
        return
    
    # Load age data
    age_data = load_age_data()
    length_data = load_length_data()
    
    patterns = {
        'c_strand': "CCCTAACCCTAACCCTAA",
        'g_strand': "GGGTTAGGGTTAGGGTTA",
        'g_strand_mutations': {
        'G>A_g1': "GGGTTAAGGTTAGGGTTA",
        'G>A_g2': "GGGTTAGAGTTAGGGTTA",
        'G>A_g3': "GGGTTAGGATTAGGGTTA",
        'G>C_g1': "GGGTTACGGTTAGGGTTA",
        'G>C_g2': "GGGTTAGCGTTAGGGTTA",
        'G>C_g3': "GGGTTAGGCTTAGGGTTA",
        'G>T_g1': "GGGTTATGGTTAGGGTTA",
        'G>T_g2': "GGGTTAGTGTTAGGGTTA",
        'G>T_g3': "GGGTTAGGTTTAGGGTTA",
        'T>A_t1': "GGGTTAGGGATAGGGTTA",
        'T>A_t2': "GGGTTAGGGTAAGGGTTA",
        'T>C_t1': "GGGTCAGGGTTAGGGTTA",
        'T>C_t2': "GGGTTACGGTTAGGGTTA",
        'T>G_t1': "GGGTGAGGGTTAGGGTTA",
        'T>G_t2': "GGGTTAGGGGTAGGGTTA",
        'A>T_a1': "GGGTTAGGGTTTGGGTTA",
        'A>G_a1': "GGGTTAGGGTTGGGGTTA",
        'A>C_a1': "GGGTTAGGGTTCGGGTTA",
        },
        'c_strand_mutations': {
        'C>A_c1': "CCCTAAACCTAACCCTAA",
        'C>A_c2': "CCCTAACACTAACCCTAA",
        'C>A_c3': "CCCTAACCATAACCCTAA",
        'C>G_c1': "CCCTAAGCCTAACCCTAA",
        'C>G_c2': "CCCTAACGCTAACCCTAA",
        'C>G_c3': "CCCTAACCGTAACCCTAA",
        'C>T_c1': "CCCTAATCCTAACCCTAA",
        'C>T_c2': "CCCTAACTCTAACCCTAA",
        'C>T_c3': "CCCTAACCTTAACCCTAA",
        'T>A_t1': "CCCTAACCCAAACCCTAA",
        'T>C_t1': "CCCTAACCCCAACCCTAA",
        'T>G_t1': "CCCTAACCCGAACCCTAA",
        'A>T_a1': "CCCTAACCCTTACCCTAA",
        'A>T_a2': "CCCTAACCCTATCCCTAA",
        'A>G_a1': "CCCTAACCCTGACCCTAA",
        'A>G_a2': "CCCTAACCCTAGCCCTAA",
        'A>C_a1': "CCCTAACCCTCACCCTAA",
        'A>C_a2': "CCCTAACCCTACCCCTAA",
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

        # Add total mutations field at the end
        fieldnames.append('total_mutations_over_total_g_strand_3xrepeats_per_1k')
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for file_path in fastq_files:
            counts = defaultdict(int)
            
            for sequence in read_fastq(file_path):
                # Count c-strand in forward direction only
                counts['c_strand'] += count_patterns(sequence, patterns['c_strand'])
                # Count g-strand in forward direction only
                counts['g_strand'] += count_patterns(sequence, patterns['g_strand'])
                # Count all mutation sub-patterns
                for group in ['g_strand_mutations', 'c_strand_mutations']:
                    for subkey, subpattern in patterns[group].items():
                        counts[f"{group}_{subkey}"] += count_patterns(sequence, subpattern)
            
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
            
            # Total mutations (sum all mutation counts)
            total_mutations = sum(counts[k] for k in mutation_keys)
            row['total_mutations_over_total_g_strand_3xrepeats_per_1k'] = per_1k(total_mutations, g_strand_total)
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
