"""Per-query placement: enumerate every branch, optimize insertion point
and pendant length by maximizing the JC69 log-likelihood."""

from __future__ import annotations

import math

from .schemas import MAX_PENDANT, NON_EVIDENCE, leaf_vector

_GOLD = 0.6180339887498949


def _golden_maximize(func, lo, hi, iterations=40):
    a, b = lo, hi
    c = b - _GOLD * (b - a)
    d = a + _GOLD * (b - a)
    fc, fd = func(c), func(d)
    for _ in range(iterations):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + _GOLD * (b - a)
            fd = func(d)
        else:
            b, d, fd = d, c, fc
            c = b - _GOLD * (b - a)
            fc = func(c)
    x = 0.5 * (a + b)
    return x, func(x)


class PlacementResult:
    def __init__(self, edge_id, distal_length, pendant_length, log_likelihood):
        self.edge_id = edge_id
        self.distal_length = distal_length
        self.pendant_length = pendant_length
        self.log_likelihood = log_likelihood
        self.weight = 0.0  # filled in after normalisation across all branches


def query_has_evidence(seq):
    """A query is placeable only if at least one column carries nucleotide
    evidence (not a gap, not a fully ambiguous code)."""
    return any(ch not in NON_EVIDENCE for ch in seq)


def optimize_on_edge(engine, edge, query_cols):
    """Maximize log-likelihood over (distal_len in [0, L], pendant in [0, 2])."""
    L = edge.length

    def best_pendant(x):
        y, ll = _golden_maximize(
            lambda t: engine.grafted_log_likelihood(edge, x, t, query_cols),
            0.0, MAX_PENDANT, iterations=30)
        return y, ll

    # coarse scan over the split point, then golden-section refinement
    grid = [L * i / 8.0 for i in range(9)]
    scored = [(best_pendant(x)[1], x) for x in grid]
    best_ll, best_x = max(scored, key=lambda t: t[0])
    idx = grid.index(best_x)
    lo = grid[idx - 1] if idx > 0 else grid[0]
    hi = grid[idx + 1] if idx < len(grid) - 1 else grid[-1]
    x, ll = _golden_maximize(lambda t: best_pendant(t)[1], lo, hi, iterations=25)
    y, ll = best_pendant(x)
    return x, y, ll


def place_query(engine, seq_id, seq):
    """Place one query on every branch; returns (results, reason)."""
    if not query_has_evidence(seq):
        return None, (
            "query carries no nucleotide evidence (only gaps or fully "
            "ambiguous codes); likelihood is flat, placement impossible"
        )
    query_cols = [leaf_vector(ch) for ch in seq]
    results = []
    for edge in engine.tree.edges:
        distal, pendant, ll = optimize_on_edge(engine, edge, query_cols)
        results.append(PlacementResult(edge.eid, distal, pendant, ll))
    # relative likelihood weights normalised across all branches
    best = max(r.log_likelihood for r in results)
    rel = [math.exp(r.log_likelihood - best) for r in results]
    total = sum(rel)
    for r, w in zip(results, rel):
        r.weight = w / total
    # sort by likelihood descending, ties broken by stable edge number
    results.sort(key=lambda r: (-r.log_likelihood, r.edge_id))
    return results, None
