#!/usr/bin/env python3

import gzip
from collections import defaultdict
import csv
import os
import time
import numpy as np
import glob
from datetime import datetime

from pattern_generator import generate_patterns

_ANALYSIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Telomeric pre-screen: skip full pattern counting on reads that don't even
# contain a canonical telomere hexamer. Cuts runtime substantially on large
# FASTQ files where most reads aren't telomeric.
_TELOMERE_ANCHORS = ('GGGTTA', 'CCCTAA')


def _is_telomeric(sequence: str) -> bool:
    """Return True if the sequence contains any canonical telomere hexamer."""
    return any(anchor in sequence for anchor in _TELOMERE_ANCHORS)


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


def load_age_data(metadata_file_path=None):
    """Load age data from greider_methods_table_s2_outliers_removed.csv (returns {} if not found)."""
    age_data = {}
    if metadata_file_path is None:
        for candidate in [
            'greider_methods_table_s2_outliers_removed.csv',
            '../analysis/greider_methods_table_s2_outliers_removed.csv',
            'analysis/greider_methods_table_s2_outliers_removed.csv',
        ]:
            if os.path.exists(candidate):
                metadata_file_path = candidate
                break
        else:
            return age_data

    with open(metadata_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            fastq_name = row['fastq file name'].replace('_', '.')
            age_data[fastq_name] = row['Age (Years)']
    return age_data


def load_length_data(metadata_file_path=None):
    """Load length data from greider_methods_table_s2_outliers_removed.csv (returns {} if not found)."""
    length_data = {}
    if metadata_file_path is None:
        for candidate in [
            'greider_methods_table_s2_outliers_removed.csv',
            '../analysis/greider_methods_table_s2_outliers_removed.csv',
            'analysis/greider_methods_table_s2_outliers_removed.csv',
        ]:
            if os.path.exists(candidate):
                metadata_file_path = candidate
                break
        else:
            return length_data

    with open(metadata_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            fastq_name = row['fastq file name'].replace('_', '.')
            length_data[fastq_name] = row['Mean Telomere Length (bps)']
    return length_data


def get_sequence_files(directory: str):
    """Get all FASTQ and FASTA files in the given directory."""
    sequence_files = []
    for ext in ['*.fastq', '*.fastq.gz', '*.fasta', '*.fasta.gz', '*.fa', '*.fa.gz', '*.fas', '*.fas.gz']:
        sequence_files.extend(glob.glob(os.path.join(directory, ext)))
    return sorted(sequence_files)  # Sort for consistent ordering


def _write_output_summary(txt_path, summary_rows, version, k, data_dir, csv_path,
                           file_type, patterns, general_mutation_map):
    """Write a human-readable plain-text summary of the analysis run."""
    total_reads = sum(r['total_reads'] for r in summary_rows)
    total_c = sum(r['c_strand'] for r in summary_rows)
    total_g = sum(r['g_strand'] for r in summary_rows)
    total_telomere = total_c + total_g
    total_mutations = sum(r['total_mutations'] for r in summary_rows)

    with open(txt_path, 'w') as f:
        f.write("=" * 72 + "\n")
        f.write("TELOMERE ANALYSIS — OUTPUT SUMMARY\n")
        f.write("=" * 72 + "\n")
        f.write(f"Run date/time       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Pattern version     : {version}  (k={k} repeats of GGGTTA/CCCTAA)\n")
        f.write(f"Input directory     : {data_dir}\n")
        f.write(f"Input file type     : {file_type}\n")
        f.write(f"Output CSV          : {csv_path}\n")
        f.write(f"Output TXT          : {txt_path}\n")
        f.write("\n")
        f.write("CANONICAL PATTERNS\n")
        f.write(f"  G-strand ({k}x): {patterns['g_strand']}\n")
        f.write(f"  C-strand ({k}x): {patterns['c_strand']}\n")
        f.write(f"  Mutations evaluated in the 2nd repeat (positions 6-11)\n")
        f.write(f"  G-strand variants : {len(patterns['g_strand_mutations'])}\n")
        f.write(f"  C-strand variants : {len(patterns['c_strand_mutations'])}\n")
        f.write("\n")
        f.write("SINGLE-BASE MODIFICATION VARIANTS — G-strand (2nd repeat)\n")
        for subkey, seq in patterns['g_strand_mutations'].items():
            f.write(f"  {subkey:<12}: {seq}\n")
        f.write("\n")
        f.write("SINGLE-BASE MODIFICATION VARIANTS — C-strand (2nd repeat)\n")
        for subkey, seq in patterns['c_strand_mutations'].items():
            f.write(f"  {subkey:<12}: {seq}\n")
        f.write("\n")
        f.write("=" * 72 + "\n")
        f.write("RUN SUMMARY\n")
        f.write("=" * 72 + "\n")
        f.write(f"Samples processed       : {len(summary_rows)}\n")
        f.write(f"Total reads             : {total_reads:,}\n")
        f.write(f"Total telomere reads    : {total_telomere:,}  (c-strand + g-strand hits)\n")
        f.write(f"  C-strand hits         : {total_c:,}\n")
        f.write(f"  G-strand hits         : {total_g:,}\n")
        f.write(f"Total mutations found   : {total_mutations:,}\n")
        if total_reads > 0:
            f.write(f"Telomere read fraction  : {total_telomere / total_reads:.6f}\n")
        f.write("\n")
        f.write("PER-SAMPLE RESULTS\n")
        f.write("-" * 72 + "\n")
        f.write(f"{'Sample':<22} {'Total Reads':>13} {'C-strand':>10} {'G-strand':>10} {'Mutations':>10}\n")
        f.write("-" * 72 + "\n")
        for r in summary_rows:
            f.write(
                f"{r['sample_id']:<22} {r['total_reads']:>13,} "
                f"{r['c_strand']:>10} {r['g_strand']:>10} {r['total_mutations']:>10}\n"
            )
        f.write("-" * 72 + "\n")
        f.write(
            f"{'TOTAL':<22} {total_reads:>13,} {total_c:>10} {total_g:>10} {total_mutations:>10}\n"
        )
        f.write("\n")
        f.write("MUTATION TYPE GROUPINGS\n")
        f.write("-" * 72 + "\n")
        for strand, mutmap in general_mutation_map.items():
            f.write(f"\n{strand}:\n")
            for mut_type, subkeys in mutmap.items():
                f.write(f"  {mut_type}: {', '.join(subkeys)}\n")
        f.write("\n")
        f.write("=" * 72 + "\n")
        f.write("END OF SUMMARY\n")
        f.write("=" * 72 + "\n")


def generate_csv(
    data_dir: str,
    k: int = 3,
    output_callback=None,
    metadata_file_path=None,
    output_csv_path=None,
    output_txt_path=None,
):
    """
    Generate a single CSV file from FASTQ/FASTA telomere data containing both
    raw counts and normalized (per-1k) metrics.

    Patterns are generated programmatically from the repeat count k via
    pattern_generator.generate_patterns() (no JSON patterns file needed).
    Reads are pre-screened for telomeric anchors before full pattern counting.

    Args:
        data_dir: Directory containing sequence files.
        k: Number of GGGTTA/CCCTAA repeat units (must be >= 2).
        output_callback: Optional callback function to receive console output.
        metadata_file_path: Optional path to metadata file (age/length table). If not
            provided or not found, Age and Telomere_Length will be left empty.
        output_csv_path: Optional path for output CSV. If None, writes to a file whose
            name encodes the pattern version, e.g. 'telomere_analysis_<k>x_repeat.csv'.
        output_txt_path: Optional path for a human-readable run summary. If given (and
            at least one sample was processed), a summary is written there.

    Returns:
        Path to the CSV that was written, or None if no sequence files were found.
    """
    patterns, general_mutation_map, version = generate_patterns(k)

    if output_csv_path is None:
        safe_version = version.replace(' ', '_')
        csv_path = os.path.join(_ANALYSIS_DIR, f"telomere_analysis_{safe_version}.csv")
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
        return None

    msg = f"Found {len(sequence_files)} FASTQ/FASTA file(s) in {data_dir}"
    print(msg)
    if output_callback:
        output_callback(msg)

    age_data = load_age_data(metadata_file_path)
    length_data = load_length_data(metadata_file_path)

    mutation_keys = []
    for group in ['g_strand_mutations', 'c_strand_mutations']:
        for subkey in patterns[group].keys():
            mutation_keys.append(f"{group}_{subkey}")

    fieldnames = ['FileName', 'Age', 'Telomere_Length', 'Total_Reads', 'c_strand', 'g_strand']
    fieldnames.extend(mutation_keys)
    fieldnames.extend([f"{mk}_per_1k" for mk in mutation_keys])
    general_mutation_headers = []
    for strand, mutmap in general_mutation_map.items():
        for mut in mutmap:
            general_mutation_headers.append(f"{strand}_{mut}_per_1k")
    fieldnames.extend(general_mutation_headers)
    fieldnames.extend([
        'g_strand_mutations_sum_per_1k', 'c_strand_mutations_sum_per_1k',
        'log_telomere_length', 'telomere_length_bin', 'mutation_rate_normalized_by_length',
    ])
    for strand, mutmap in general_mutation_map.items():
        for mut in mutmap:
            fieldnames.append(f"{strand}_{mut}_sum_per_1k")
    fieldnames.append('total_mutations_per_1k_strand_specific')
    fieldnames.extend([
        'total_mutations_over_total_g_strand_per_1k',
        'total_mutations_over_total_c_strand_per_1k',
    ])

    summary_rows = []

    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for file_path in sequence_files:
            filename = os.path.basename(file_path)
            sample_id = filename
            for ext in ['.fastq.gz', '.fastq', '.fasta.gz', '.fasta', '.fa.gz', '.fa', '.fas.gz', '.fas']:
                if sample_id.endswith(ext):
                    sample_id = sample_id[:-len(ext)]
                    break

            counts = defaultdict(int)
            total_reads = 0
            telomeric_reads = 0
            _PROGRESS_INTERVAL = 1_000_000
            _loop_start = time.time()

            for sequence in read_sequence_file(file_path):
                total_reads += 1
                if not _is_telomeric(sequence):
                    continue
                telomeric_reads += 1
                counts['c_strand'] += count_patterns(sequence, patterns['c_strand'])
                counts['g_strand'] += count_patterns(sequence, patterns['g_strand'])
                for group in ['g_strand_mutations', 'c_strand_mutations']:
                    for subkey, subpattern in patterns[group].items():
                        counts[f"{group}_{subkey}"] += count_patterns(sequence, subpattern)
                if total_reads % _PROGRESS_INTERVAL == 0:
                    _elapsed = time.time() - _loop_start
                    _rate = total_reads / _elapsed if _elapsed > 0 else 0
                    print(
                        f'  [progress] {total_reads:>12,} reads scanned  '
                        f'| {telomeric_reads:,} telomeric  '
                        f'| {_rate:,.0f} reads/sec  '
                        f'| elapsed {int(_elapsed // 60)}m {int(_elapsed % 60):02d}s',
                        flush=True
                    )

            age = age_data.get(sample_id, '')
            length = length_data.get(sample_id, '')
            g_strand_total = counts['g_strand']
            c_strand_total = counts['c_strand']
            g_strand_mutations_total = sum(counts[mk] for mk in counts if mk.startswith('g_strand_mutations_'))
            c_strand_mutations_total = sum(counts[mk] for mk in counts if mk.startswith('c_strand_mutations_'))
            g_strand_normalizer = g_strand_total + g_strand_mutations_total
            c_strand_normalizer = c_strand_total + c_strand_mutations_total

            def per_1k_ss(val, norm):
                return (val / norm) * 1000 if norm > 0 else 0

            row = {
                'FileName': filename, 'Age': age, 'Telomere_Length': length,
                'Total_Reads': total_reads,
                'c_strand': c_strand_total, 'g_strand': g_strand_total,
            }
            for mk in mutation_keys:
                row[mk] = counts.get(mk, 0)
                norm = g_strand_normalizer if mk.startswith('g_strand_mutations_') else c_strand_normalizer
                row[f"{mk}_per_1k"] = per_1k_ss(counts.get(mk, 0), norm)
            for strand, mutmap in general_mutation_map.items():
                norm = g_strand_normalizer if strand == 'g_strand' else c_strand_normalizer
                for mut, subtypes in mutmap.items():
                    total = sum(counts.get(f"{strand}_mutations_{st}", 0) for st in subtypes)
                    row[f"{strand}_{mut}_per_1k"] = per_1k_ss(total, norm)
            for strand, mutmap in general_mutation_map.items():
                for mut, subtypes in mutmap.items():
                    row[f"{strand}_{mut}_sum_per_1k"] = sum(
                        row.get(f"{strand}_mutations_{st}_per_1k", 0) for st in subtypes
                    )
            total_mutations = sum(counts[mk] for mk in mutation_keys)
            if total_mutations > 0 and (g_strand_normalizer > 0 or c_strand_normalizer > 0):
                g_w = g_strand_mutations_total / total_mutations
                c_w = c_strand_mutations_total / total_mutations
                w_norm = g_w * g_strand_normalizer + c_w * c_strand_normalizer
                row['total_mutations_per_1k_strand_specific'] = per_1k_ss(total_mutations, w_norm)
            else:
                row['total_mutations_per_1k_strand_specific'] = 0
            row['total_mutations_over_total_g_strand_per_1k'] = per_1k_ss(total_mutations, g_strand_normalizer)
            row['total_mutations_over_total_c_strand_per_1k'] = per_1k_ss(total_mutations, c_strand_normalizer)
            row['g_strand_mutations_sum_per_1k'] = sum(
                v for rk, v in row.items() if rk.startswith('g_strand_mutations') and rk.endswith('_per_1k')
            )
            row['c_strand_mutations_sum_per_1k'] = sum(
                v for rk, v in row.items() if rk.startswith('c_strand_mutations') and rk.endswith('_per_1k')
            )
            try:
                row['log_telomere_length'] = np.log1p(float(length)) if length else 0
            except Exception:
                row['log_telomere_length'] = 0
            try:
                lv = float(length)
                row['telomere_length_bin'] = 'short' if lv < 5000 else ('medium' if lv < 8000 else 'long')
            except Exception:
                row['telomere_length_bin'] = 'unknown'
            try:
                lv = float(length)
                row['mutation_rate_normalized_by_length'] = row['g_strand_mutations_sum_per_1k'] / lv if lv > 0 else 0
            except Exception:
                row['mutation_rate_normalized_by_length'] = 0

            writer.writerow(row)

            messages = [
                f"\nSample  : {sample_id}  ({filename})",
                f"  Age              : {age}",
                f"  Telomere Length  : {length}",
                f"  Total Reads      : {total_reads:,}",
                f"  Telomeric Reads  : {telomeric_reads:,}",
                f"  {version} c-strand: {c_strand_total}",
                f"  {version} g-strand: {g_strand_total}",
                f"  G-strand mutations: {g_strand_mutations_total}",
                f"  C-strand mutations: {c_strand_mutations_total}",
                f"  G-strand normalizer: {g_strand_normalizer}",
                f"  C-strand normalizer: {c_strand_normalizer}",
                f"  Total mutations  : {total_mutations}",
            ]
            if g_strand_total == 0 and c_strand_total == 0:
                messages.append(f"  WARNING: No telomere sequences found in {sample_id}")
            for message in messages:
                print(message)
                if output_callback:
                    output_callback(message)

            summary_rows.append({
                'sample_id': sample_id, 'file': filename,
                'total_reads': total_reads,
                'c_strand': c_strand_total, 'g_strand': g_strand_total,
                'g_strand_mutations': g_strand_mutations_total,
                'c_strand_mutations': c_strand_mutations_total,
                'total_mutations': total_mutations,
                'age': age, 'telomere_length': length,
            })

    if output_txt_path and summary_rows:
        _write_output_summary(
            output_txt_path, summary_rows, version, k, data_dir, csv_path,
            'FASTQ/FASTA', patterns, general_mutation_map,
        )
        print(f"\nSummary written to: {output_txt_path}")

    return csv_path


if __name__ == "__main__":
    generate_csv(".")
