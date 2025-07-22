import gzip

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
    features['total_mutations_over_total_g_strand_3xrepeats_per_1k'] = per_1k(total_mutations, g_strand_total)
    return features
