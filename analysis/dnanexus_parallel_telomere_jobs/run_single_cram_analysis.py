#!/usr/bin/env python3
"""Run telomere analysis on one CRAM and emit one-row CSV.

This script ports the phase-aware k-repeat logic from
analysis/telomere_analysis_ukb_cram_update_fixed.ipynb to a CLI job worker.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import pysam

G_REPEAT_UNIT = "GGGTTA"
C_REPEAT_UNIT = "CCCTAA"

G_STRAND_MUT_SPEC = {
    "G>A_g1": (0, "A"), "G>A_g2": (1, "A"), "G>A_g3": (2, "A"),
    "G>C_g1": (0, "C"), "G>C_g2": (1, "C"), "G>C_g3": (2, "C"),
    "G>T_g1": (0, "T"), "G>T_g2": (1, "T"), "G>T_g3": (2, "T"),
    "T>A_t1": (3, "A"), "T>A_t2": (4, "A"),
    "T>C_t1": (3, "C"), "T>C_t2": (4, "C"),
    "T>G_t1": (3, "G"), "T>G_t2": (4, "G"),
    "A>T_a1": (5, "T"), "A>G_a1": (5, "G"), "A>C_a1": (5, "C"),
}

C_STRAND_MUT_SPEC = {
    "C>A_c1": (0, "A"), "C>A_c2": (1, "A"), "C>A_c3": (2, "A"),
    "C>G_c1": (0, "G"), "C>G_c2": (1, "G"), "C>G_c3": (2, "G"),
    "C>T_c1": (0, "T"), "C>T_c2": (1, "T"), "C>T_c3": (2, "T"),
    "T>A_t1": (3, "A"), "T>C_t1": (3, "C"), "T>G_t1": (3, "G"),
    "A>T_a1": (4, "T"), "A>T_a2": (5, "T"),
    "A>G_a1": (4, "G"), "A>G_a2": (5, "G"),
    "A>C_a1": (4, "C"), "A>C_a2": (5, "C"),
}

GENERAL_MUTATION_MAP = {
    "g_strand": {
        "G>A": ["G>A_g1", "G>A_g2", "G>A_g3"],
        "G>C": ["G>C_g1", "G>C_g2", "G>C_g3"],
        "G>T": ["G>T_g1", "G>T_g2", "G>T_g3"],
        "T>A": ["T>A_t1", "T>A_t2"], "T>C": ["T>C_t1", "T>C_t2"], "T>G": ["T>G_t1", "T>G_t2"],
        "A>T": ["A>T_a1"], "A>G": ["A>G_a1"], "A>C": ["A>C_a1"],
    },
    "c_strand": {
        "C>A": ["C>A_c1", "C>A_c2", "C>A_c3"],
        "C>G": ["C>G_c1", "C>G_c2", "C>G_c3"],
        "C>T": ["C>T_c1", "C>T_c2", "C>T_c3"],
        "T>A": ["T>A_t1"], "T>C": ["T>C_t1"], "T>G": ["T>G_t1"],
        "A>T": ["A>T_a1", "A>T_a2"], "A>G": ["A>G_a1", "A>G_a2"], "A>C": ["A>C_a1", "A>C_a2"],
    },
}


def _hexamer_rotations(hexamer: str) -> List[str]:
    return [hexamer[i:] + hexamer[:i] for i in range(len(hexamer))]


def _make_kmer_phases(canonical_hex: str, k: int, mut_pos: int | None = None, mut_base: str | None = None) -> List[str]:
    long_seq = list(canonical_hex * (k + 1))
    if mut_pos is not None and mut_base is not None:
        long_seq[len(canonical_hex) + mut_pos] = mut_base
    long_seq = "".join(long_seq)
    return [long_seq[phase: phase + len(canonical_hex) * k] for phase in range(len(canonical_hex))]


G_HEX_ROTATIONS = _hexamer_rotations(G_REPEAT_UNIT)
C_HEX_ROTATIONS = _hexamer_rotations(C_REPEAT_UNIT)
TELOMERE_ANCHORS = frozenset(G_HEX_ROTATIONS + C_HEX_ROTATIONS)


def generate_patterns(k: int) -> Tuple[Dict[str, object], Dict[str, Dict[str, List[str]]], str]:
    if k < 2:
        raise ValueError(f"k must be >= 2; got {k}")
    g_can = _make_kmer_phases(G_REPEAT_UNIT, k)
    c_can = _make_kmer_phases(C_REPEAT_UNIT, k)
    g_mut = {name: _make_kmer_phases(G_REPEAT_UNIT, k, pos, base) for name, (pos, base) in G_STRAND_MUT_SPEC.items()}
    c_mut = {name: _make_kmer_phases(C_REPEAT_UNIT, k, pos, base) for name, (pos, base) in C_STRAND_MUT_SPEC.items()}
    return {
        "g_strand": g_can,
        "c_strand": c_can,
        "g_strand_mutations": g_mut,
        "c_strand_mutations": c_mut,
    }, GENERAL_MUTATION_MAP, f"{k}x_repeat"


def _is_telomeric(sequence: str) -> bool:
    return any(anchor in sequence for anchor in TELOMERE_ANCHORS)


def _count_canonical_hexamers(sequence: str) -> Tuple[int, int]:
    return sum(sequence.count(r) for r in G_HEX_ROTATIONS), sum(sequence.count(r) for r in C_HEX_ROTATIONS)


def _count_pattern_phases(sequence: str, phase_list: List[str]) -> int:
    return sum(sequence.count(p) for p in phase_list)


def _iter_cram_sequences(cram_path: str, reference_fasta: str) -> Iterable[str]:
    with pysam.AlignmentFile(cram_path, "rc", reference_filename=reference_fasta) as af:
        for read in af.fetch(until_eof=True):
            if read.query_sequence:
                yield read.query_sequence


def _pick_first_existing(df: pd.DataFrame, candidates: List[str]) -> str | None:
    by_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in by_lower:
            return by_lower[cand.lower()]
    return None


def _load_metadata_map(metadata_csv: str) -> Dict[str, Dict[str, object]]:
    df = pd.read_csv(metadata_csv)
    id_col = next((c for c in df.columns if c.lower() in {"participant id", "participant_id", "eid", "id"}), None)
    age_col = _pick_first_existing(df, ["Age at recruitment", "Age"])
    tel_col = _pick_first_existing(
        df,
        [
            "Telomere_Length",
            "Mean Telomere Length (bps)",
            "Adjusted T/S ratio | Instance 0",
            "Unadjusted T/S ratio | Instance 0",
            "Z-adjusted T/S log | Instance 0",
        ],
    )
    if id_col is None:
        raise ValueError(f"Could not find participant ID column in {metadata_csv}.")
    mapping: Dict[str, Dict[str, object]] = {}
    for _, row in df.iterrows():
        pid = str(row[id_col]).strip()
        mapping[pid] = {
            "Age": row[age_col] if age_col else "",
            "Telomere_Length": row[tel_col] if tel_col else "",
            "all_fields": {col: row[col] for col in df.columns if col != id_col},
        }
    return mapping


def _safe_float(v: object) -> float | None:
    try:
        if v == "" or pd.isna(v):
            return None
        return float(v)
    except Exception:
        return None


def analyze_single_cram(cram_path: str, participant_id: str, reference_fasta: str, metadata_csv: str, k: int) -> Dict[str, object]:
    patterns, general_mutation_map, _ = generate_patterns(k)
    mutation_keys: List[str] = []
    for group in ["g_strand_mutations", "c_strand_mutations"]:
        for subkey in patterns[group].keys():
            mutation_keys.append(f"{group}_{subkey}")

    metadata = _load_metadata_map(metadata_csv).get(participant_id, {"Age": "", "Telomere_Length": "", "all_fields": {}})
    counts: Dict[str, int] = defaultdict(int)
    total_reads = 0
    telomeric_reads = 0
    telomeric_hex_threshold = 4

    for seq in _iter_cram_sequences(cram_path, reference_fasta):
        total_reads += 1
        if not _is_telomeric(seq):
            continue
        g_hex, c_hex = _count_canonical_hexamers(seq)
        if (g_hex + c_hex) < telomeric_hex_threshold:
            continue
        telomeric_reads += 1
        counts["g_strand"] += _count_pattern_phases(seq, patterns["g_strand"])
        counts["c_strand"] += _count_pattern_phases(seq, patterns["c_strand"])
        for group in ["g_strand_mutations", "c_strand_mutations"]:
            for subkey, phase_list in patterns[group].items():
                counts[f"{group}_{subkey}"] += _count_pattern_phases(seq, phase_list)

    g_strand_total = counts["g_strand"]
    c_strand_total = counts["c_strand"]
    g_strand_mut_total = sum(counts[k] for k in counts if k.startswith("g_strand_mutations_"))
    c_strand_mut_total = sum(counts[k] for k in counts if k.startswith("c_strand_mutations_"))
    g_norm = g_strand_total + g_strand_mut_total
    c_norm = c_strand_total + c_strand_mut_total

    def per_1k_ss(val: float, norm: float) -> float:
        return (val / norm) * 1000.0 if norm > 0 else 0.0

    row: Dict[str, object] = {
        "FileName": participant_id,
        "Age": metadata["Age"],
        "Telomere_Length": metadata["Telomere_Length"],
        "Total_Reads": total_reads,
        "Telomeric_Reads": telomeric_reads,
        "c_strand": c_strand_total,
        "g_strand": g_strand_total,
    }
    row.update(metadata.get("all_fields", {}))

    for mk in mutation_keys:
        val = counts.get(mk, 0)
        row[mk] = val
        norm = g_norm if mk.startswith("g_strand_mutations_") else c_norm
        row[f"{mk}_per_1k"] = per_1k_ss(val, norm)

    for strand, mutmap in general_mutation_map.items():
        norm = g_norm if strand == "g_strand" else c_norm
        for mut, subtypes in mutmap.items():
            total = sum(counts.get(f"{strand}_mutations_{subtype}", 0) for subtype in subtypes)
            row[f"{strand}_{mut}_per_1k"] = per_1k_ss(total, norm)

    for strand, mutmap in general_mutation_map.items():
        for mut, subtypes in mutmap.items():
            row[f"{strand}_{mut}_sum_per_1k"] = sum(row.get(f"{strand}_mutations_{subtype}_per_1k", 0) for subtype in subtypes)

    total_mutations = sum(counts[mk] for mk in mutation_keys)
    if total_mutations > 0 and (g_norm > 0 or c_norm > 0):
        g_w = g_strand_mut_total / total_mutations
        c_w = c_strand_mut_total / total_mutations
        weighted_norm = g_w * g_norm + c_w * c_norm
        row["total_mutations_per_1k_strand_specific"] = per_1k_ss(total_mutations, weighted_norm)
    else:
        row["total_mutations_per_1k_strand_specific"] = 0

    row["total_mutations_over_total_g_strand_per_1k"] = per_1k_ss(total_mutations, g_norm)
    row["total_mutations_over_total_c_strand_per_1k"] = per_1k_ss(total_mutations, c_norm)
    row["g_strand_mutations_sum_per_1k"] = sum(
        v for k, v in row.items() if k.startswith("g_strand_mutations") and k.endswith("_per_1k")
    )
    row["c_strand_mutations_sum_per_1k"] = sum(
        v for k, v in row.items() if k.startswith("c_strand_mutations") and k.endswith("_per_1k")
    )

    length_f = _safe_float(row["Telomere_Length"])
    row["log_telomere_length"] = float(np.log1p(length_f)) if length_f is not None else 0.0
    if length_f is None:
        row["telomere_length_bin"] = "unknown"
        row["mutation_rate_normalized_by_length"] = 0.0
    else:
        row["telomere_length_bin"] = "short" if length_f < 5000 else ("medium" if length_f < 8000 else "long")
        row["mutation_rate_normalized_by_length"] = row["g_strand_mutations_sum_per_1k"] / length_f if length_f > 0 else 0.0

    return row


def write_single_row_csv(row: Dict[str, object], output_csv: str) -> None:
    parent = os.path.dirname(output_csv)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run telomere analysis for one participant CRAM.")
    parser.add_argument("--cram", required=True, help="Local CRAM path")
    parser.add_argument("--participant-id", required=True, help="Participant ID for metadata lookup")
    parser.add_argument("--reference-fasta", required=True, help="Local reference FASTA path")
    parser.add_argument("--metadata-csv", required=True, help="Local UKB metadata CSV (participant ID, age, sex, smoking, ratios, etc.)")
    parser.add_argument("--output-csv", required=True, help="Output CSV path (single row)")
    parser.add_argument("--k", type=int, default=3, help="Repeat count for k-mer patterns (default: 3)")
    args = parser.parse_args()

    row = analyze_single_cram(
        cram_path=args.cram,
        participant_id=str(args.participant_id),
        reference_fasta=args.reference_fasta,
        metadata_csv=args.metadata_csv,
        k=args.k,
    )
    write_single_row_csv(row, args.output_csv)
    print(f"Wrote single-sample output: {args.output_csv}")


if __name__ == "__main__":
    main()
