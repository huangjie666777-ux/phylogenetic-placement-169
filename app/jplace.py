"""jplace version 3 rendering for placement results."""

from __future__ import annotations

FIELDS = ["edge_num", "likelihood", "like_weight_ratio",
          "distal_length", "pendant_length"]


def build_jplace(tree, placements, program="edna-fixed-tree-placer"):
    """placements: list of (query_id, [PlacementResult, ...])."""
    entries = []
    for query_id, results in placements:
        p = [[r.edge_id, r.log_likelihood, r.weight,
              r.distal_length, r.pendant_length] for r in results]
        entries.append({"p": p, "n": [query_id]})
    return {
        "tree": tree.to_newick(annotate_edges=True),
        "placements": entries,
        "metadata": {
            "program": program,
            "model": "JC69",
            "note": ("like_weight_ratio is the relative likelihood weight "
                     "normalised across all branches; it is relative "
                     "support, not a probability of correct species "
                     "assignment."),
        },
        "version": 3,
        "fields": FIELDS,
    }
