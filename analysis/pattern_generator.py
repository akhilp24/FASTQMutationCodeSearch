"""
Telomere pattern generator for k-repeat analysis.

Instead of loading patterns from a JSON file, this module generates telomere
repeat patterns programmatically based on a given k value (number of GGGTTA repeats).

The telomere G-strand repeat unit is GGGTTA and C-strand is CCCTAA.
Single-base modifications are always evaluated in the SECOND repeat within
a k-repeat pattern, since that is the biologically relevant site for
identifying single base modification variants in telomeric reads.
"""

G_REPEAT_UNIT = "GGGTTA"
C_REPEAT_UNIT = "CCCTAA"

# Mutation specification: {subkey: (position_within_repeat_unit, original_base, new_base)}
# G-strand: GGGTTA  →  G=pos0, G=pos1, G=pos2, T=pos3, T=pos4, A=pos5
_G_STRAND_MUT_SPEC = {
    "G>A_g1": (0, "G", "A"),
    "G>A_g2": (1, "G", "A"),
    "G>A_g3": (2, "G", "A"),
    "G>C_g1": (0, "G", "C"),
    "G>C_g2": (1, "G", "C"),
    "G>C_g3": (2, "G", "C"),
    "G>T_g1": (0, "G", "T"),
    "G>T_g2": (1, "G", "T"),
    "G>T_g3": (2, "G", "T"),
    "T>A_t1": (3, "T", "A"),
    "T>A_t2": (4, "T", "A"),
    "T>C_t1": (3, "T", "C"),
    "T>C_t2": (4, "T", "C"),
    "T>G_t1": (3, "T", "G"),
    "T>G_t2": (4, "T", "G"),
    "A>T_a1": (5, "A", "T"),
    "A>G_a1": (5, "A", "G"),
    "A>C_a1": (5, "A", "C"),
}

# C-strand: CCCTAA  →  C=pos0, C=pos1, C=pos2, T=pos3, A=pos4, A=pos5
_C_STRAND_MUT_SPEC = {
    "C>A_c1": (0, "C", "A"),
    "C>A_c2": (1, "C", "A"),
    "C>A_c3": (2, "C", "A"),
    "C>G_c1": (0, "C", "G"),
    "C>G_c2": (1, "C", "G"),
    "C>G_c3": (2, "C", "G"),
    "C>T_c1": (0, "C", "T"),
    "C>T_c2": (1, "C", "T"),
    "C>T_c3": (2, "C", "T"),
    "T>A_t1": (3, "T", "A"),
    "T>C_t1": (3, "T", "C"),
    "T>G_t1": (3, "T", "G"),
    "A>T_a1": (4, "A", "T"),
    "A>T_a2": (5, "A", "T"),
    "A>G_a1": (4, "A", "G"),
    "A>G_a2": (5, "A", "G"),
    "A>C_a1": (4, "A", "C"),
    "A>C_a2": (5, "A", "C"),
}

# Maps broad mutation categories to the specific subkeys generated above.
# This is fixed regardless of k.
GENERAL_MUTATION_MAP = {
    "g_strand": {
        "G>A": ["G>A_g1", "G>A_g2", "G>A_g3"],
        "G>C": ["G>C_g1", "G>C_g2", "G>C_g3"],
        "G>T": ["G>T_g1", "G>T_g2", "G>T_g3"],
        "T>A": ["T>A_t1", "T>A_t2"],
        "T>C": ["T>C_t1", "T>C_t2"],
        "T>G": ["T>G_t1", "T>G_t2"],
        "A>T": ["A>T_a1"],
        "A>G": ["A>G_a1"],
        "A>C": ["A>C_a1"],
    },
    "c_strand": {
        "C>A": ["C>A_c1", "C>A_c2", "C>A_c3"],
        "C>G": ["C>G_c1", "C>G_c2", "C>G_c3"],
        "C>T": ["C>T_c1", "C>T_c2", "C>T_c3"],
        "T>A": ["T>A_t1"],
        "T>C": ["T>C_t1"],
        "T>G": ["T>G_t1"],
        "A>T": ["A>T_a1", "A>T_a2"],
        "A>G": ["A>G_a1", "A>G_a2"],
        "A>C": ["A>C_a1", "A>C_a2"],
    },
}


def generate_patterns(k: int):
    """
    Generate telomere patterns for k repeats of the GGGTTA/CCCTAA unit.

    Single-base modifications are evaluated in the SECOND repeat only.
    For example, k=3 gives canonical G-strand GGGTTAGGGTTAGGGTTA, and all
    variant patterns mutate exactly one base within the middle GGGTTA block.

    Args:
        k: Number of telomere repeat units (must be >= 2 so the second repeat exists).

    Returns:
        Tuple of (patterns, general_mutation_map, version) where:
          - patterns: dict with keys 'g_strand', 'c_strand',
            'g_strand_mutations' (dict subkey→sequence),
            'c_strand_mutations' (dict subkey→sequence)
          - general_mutation_map: dict grouping subkeys by mutation category
          - version: label string, e.g. '3x_repeat'
    """
    if k < 2:
        raise ValueError(
            f"k must be >= 2 (the second repeat must exist for mutation analysis); got k={k}"
        )

    g_canonical = G_REPEAT_UNIT * k
    c_canonical = C_REPEAT_UNIT * k

    # The second repeat starts at index 6 (one repeat_unit length into the sequence).
    offset = len(G_REPEAT_UNIT)  # == 6

    def _apply(canonical: str, pos_in_repeat: int, new_base: str) -> str:
        idx = offset + pos_in_repeat
        return canonical[:idx] + new_base + canonical[idx + 1:]

    g_mutations = {
        subkey: _apply(g_canonical, pos, new_base)
        for subkey, (pos, _orig, new_base) in _G_STRAND_MUT_SPEC.items()
    }

    c_mutations = {
        subkey: _apply(c_canonical, pos, new_base)
        for subkey, (pos, _orig, new_base) in _C_STRAND_MUT_SPEC.items()
    }

    patterns = {
        "g_strand": g_canonical,
        "c_strand": c_canonical,
        "g_strand_mutations": g_mutations,
        "c_strand_mutations": c_mutations,
    }

    return patterns, GENERAL_MUTATION_MAP, f"{k}x_repeat"


def verify_patterns(k: int, verbose: bool = False) -> bool:
    """
    Run a quick sanity check that generated patterns are the correct length
    and differ from canonical by exactly one base in the second repeat.

    Returns True if all checks pass.
    """
    patterns, _, version = generate_patterns(k)
    canonical_g = patterns["g_strand"]
    canonical_c = patterns["c_strand"]
    expected_len = len(G_REPEAT_UNIT) * k

    ok = True
    for subkey, seq in patterns["g_strand_mutations"].items():
        if len(seq) != expected_len:
            print(f"FAIL {subkey}: length {len(seq)} != {expected_len}")
            ok = False
            continue
        diffs = [i for i in range(expected_len) if seq[i] != canonical_g[i]]
        if len(diffs) != 1:
            print(f"FAIL {subkey}: {len(diffs)} differences from canonical (expected 1): positions {diffs}")
            ok = False
        elif verbose:
            pos = diffs[0]
            print(f"OK   {subkey}: pos {pos} ({canonical_g[pos]}->{seq[pos]}), in repeat {pos // 6 + 1}")

    for subkey, seq in patterns["c_strand_mutations"].items():
        if len(seq) != expected_len:
            print(f"FAIL {subkey}: length {len(seq)} != {expected_len}")
            ok = False
            continue
        diffs = [i for i in range(expected_len) if seq[i] != canonical_c[i]]
        if len(diffs) != 1:
            print(f"FAIL {subkey}: {len(diffs)} differences from canonical (expected 1): positions {diffs}")
            ok = False
        elif verbose:
            pos = diffs[0]
            print(f"OK   {subkey}: pos {pos} ({canonical_c[pos]}->{seq[pos]}), in repeat {pos // 6 + 1}")

    if ok and verbose:
        print(f"All patterns valid for k={k} ({version})")
    return ok


if __name__ == "__main__":
    for k in [2, 3, 4, 5, 7]:
        ok = verify_patterns(k, verbose=True)
        assert ok, f"Pattern verification failed for k={k}"
    print("All verifications passed.")
