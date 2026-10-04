"""Alphabet, IUPAC handling and input validation helpers."""

from __future__ import annotations

MIN_REFS = 3
MAX_REFS = 20
MIN_QUERIES = 1
MAX_QUERIES = 5
MAX_COLUMNS = 500
MAX_PENDANT = 2.0

# IUPAC degenerate base codes -> set of compatible state indices (A,C,G,T).
STATE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}

IUPAC_STATES: dict[str, frozenset[int]] = {
    "A": frozenset({0}),
    "C": frozenset({1}),
    "G": frozenset({2}),
    "T": frozenset({3}),
    "R": frozenset({0, 2}),
    "Y": frozenset({1, 3}),
    "S": frozenset({1, 2}),
    "W": frozenset({0, 3}),
    "K": frozenset({2, 3}),
    "M": frozenset({0, 1}),
    "B": frozenset({1, 2, 3}),
    "D": frozenset({0, 2, 3}),
    "H": frozenset({0, 1, 3}),
    "V": frozenset({0, 1, 2}),
    "N": frozenset({0, 1, 2, 3}),
    "-": frozenset(),  # gap: treated as missing / unknown
    "?": frozenset(),
}

# Codes that carry no nucleotide evidence at all.
NON_EVIDENCE = frozenset({"-", "?", "N"})


class InputError(ValueError):
    """Validation failure that must reject the whole submission (HTTP 422)."""


def leaf_vector(code: str) -> list[float]:
    """Observed-state vector for a leaf at one column.

    Ambiguity codes give 1.0 to every compatible state; gaps and fully
    unknown codes give 1.0 to all states (no evidence).
    """
    states = IUPAC_STATES[code]
    if not states:
        return [1.0, 1.0, 1.0, 1.0]
    return [1.0 if i in states else 0.0 for i in range(4)]


def validate_alignment(records, kind, minimum, maximum):
    """Validate one FASTA record set; returns the column count."""
    if not (minimum <= len(records) <= maximum):
        raise InputError(
            f"{kind}: expected {minimum}..{maximum} sequences, got {len(records)}"
        )
    width = None
    for seq_id, seq in records:
        if not seq:
            raise InputError(f"{kind} '{seq_id}': empty sequence")
        if width is None:
            width = len(seq)
        elif len(seq) != width:
            raise InputError(
                f"{kind} '{seq_id}': length {len(seq)} != {width}; "
                "all sequences must be aligned to equal length"
            )
        for col, ch in enumerate(seq, start=1):
            if ch not in IUPAC_STATES:
                raise InputError(
                    f"{kind} '{seq_id}' column {col}: illegal character {ch!r}"
                )
    assert width is not None
    if width > MAX_COLUMNS:
        raise InputError(f"{kind}: {width} columns exceeds limit {MAX_COLUMNS}")
    return width


def check_unique_ids(refs, queries):
    seen = set()
    for seq_id, _ in list(refs) + list(queries):
        if seq_id in seen:
            raise InputError(f"duplicate sequence id '{seq_id}'")
        seen.add(seq_id)
