import gzip
import csv
import os

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

def read_fastq(file_path):
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

def count_patterns(sequence, pattern):
    return sequence.count(pattern)

def extract_features_from_fastq(file_path):
    from collections import defaultdict
    counts = defaultdict(int)
    for sequence in read_fastq(file_path):
        counts['c_strand'] += count_patterns(sequence, patterns['c_strand'])
        counts['g_strand'] += count_patterns(sequence, patterns['g_strand'])
        for group in ['g_strand_mutations', 'c_strand_mutations']:
            for subkey, subpattern in patterns[group].items():
                counts[f"{group}_{subkey}"] += count_patterns(sequence, subpattern)
    g_strand_total = counts['g_strand']
    c_strand_total = counts['c_strand']
    mutation_keys = []
    for group in ['g_strand_mutations', 'c_strand_mutations']:
        for subkey in patterns[group].keys():
            mutation_keys.append(f"{group}_{subkey}")
    def per_1k(val, total):
        return (val / total) * 1000 if total > 0 else 0
    features = {
        'c_strand': c_strand_total,
        'g_strand': g_strand_total,
    }
    for k in mutation_keys:
        features[k] = counts.get(k, 0)
        if k.startswith('g_strand_mutations'):
            norm_total = g_strand_total
        elif k.startswith('c_strand_mutations'):
            norm_total = c_strand_total
        else:
            norm_total = g_strand_total
        features[f"{k}_per_1k"] = per_1k(counts.get(k, 0), norm_total)
    total_mutations = sum(counts[k] for k in mutation_keys)
    features['total_mutations_over_total_g_strand_2xrepeats_per_1k'] = per_1k(total_mutations, g_strand_total)

    # Add Telomere_Length from greider_methods_table_s2.csv
    # Extract filename base (without path, .fastq, or .gz)
    filename = os.path.basename(file_path)
    filename_base = filename.replace('.fastq', '').replace('.gz', '')
    telomere_length = ''
    try:
        with open('/Users/akhilpeddikuppa/FieldLab/GreiderCodeSearch/analysis/greider_methods_table_s2.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                fastq_name = row['fastq file name'].replace('_', '.')
                if fastq_name == filename_base:
                    telomere_length = row['Mean Telomere Length (bps)']
                    break
    except Exception as e:
        telomere_length = ''
    features['Telomere_Length'] = telomere_length
    return features
