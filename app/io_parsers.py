"""FASTA and Newick parsing built on Biopython."""

from __future__ import annotations

import math
from io import StringIO

from Bio import Phylo

from .schemas import InputError


def parse_fasta(text, kind):
    records = []
    seq_id = None
    chunks = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id is not None:
                records.append((seq_id, "".join(chunks)))
            seq_id = line[1:].strip().split()[0] if line[1:].strip() else ""
            if not seq_id:
                raise InputError(f"{kind} FASTA line {lineno}: empty identifier")
            chunks = []
        else:
            if seq_id is None:
                raise InputError(
                    f"{kind} FASTA line {lineno}: sequence data before first header"
                )
            chunks.append(line.upper())
    if seq_id is not None:
        records.append((seq_id, "".join(chunks)))
    if not records:
        raise InputError(f"{kind} FASTA: no records found")
    return records


def parse_newick(text):
    """Parse the reference Newick string and validate branch lengths."""
    try:
        tree = Phylo.read(StringIO(text.strip()), "newick")
    except Exception as exc:  # Biopython raises several parser errors
        raise InputError(f"newick: parse error: {exc}") from exc
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        bl = clade.branch_length
        if bl is None:
            raise InputError(
                f"newick: branch above '{clade.name or 'internal node'}' "
                "has no length; all branches need positive finite lengths"
            )
        if not math.isfinite(bl) or bl <= 0.0:
            raise InputError(
                f"newick: branch above '{clade.name or 'internal node'}' "
                f"has non-positive or non-finite length {bl!r}"
            )
    return tree.root
