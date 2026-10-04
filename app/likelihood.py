"""JC69 likelihood of a fixed alignment on the reference tree.

Felsenstein's pruning algorithm with per-site scaling so long alignments
cannot underflow.  Ancestral states are summed over exactly; no Hamming
distance or nearest-neighbour shortcuts are used anywhere.
"""

from __future__ import annotations

import math

from .schemas import leaf_vector


def jc69_probs(t):
    """Return (p_same, p_diff) for JC69 at branch length t."""
    if t <= 0.0:
        return 1.0, 0.0
    e = math.exp(-4.0 * t / 3.0)
    return 0.25 + 0.75 * e, 0.25 - 0.25 * e


class LikelihoodEngine:
    """Caches per-column leaf vectors and directed partial likelihoods."""

    def __init__(self, tree, ref_records):
        self.tree = tree
        self.n_cols = len(ref_records[0][1])
        # leaf column vectors: node nid -> list over columns of [4] floats
        self.leaf_cols = {}
        for seq_id, seq in ref_records:
            node = tree.leaf(seq_id)
            self.leaf_cols[node.nid] = [leaf_vector(ch) for ch in seq]
        # directed partials cache: (nid, parent_nid) -> (cols, log_scale, col_sums)
        self._partials_cache = {}

    def _propagate(self, cols, length):
        """Propagate scaled partial vectors across a branch of given length."""
        p_same, p_diff = jc69_probs(length)
        out = []
        for cv in cols:
            total = cv[0] + cv[1] + cv[2] + cv[3]
            out.append([cv[s] * p_same + (total - cv[s]) * p_diff for s in range(4)])
        return out

    def upward_partials(self, node, parent):
        """Scaled partial likelihoods at node for the subtree away from
        parent.  Returns (cols x 4 scaled vectors, summed log scale)."""
        key = (node.nid, parent.nid if parent else -1)
        if key in self._partials_cache:
            return self._partials_cache[key]
        if node.is_leaf:
            result = (self.leaf_cols[node.nid], 0.0)
            self._partials_cache[key] = result
            return result
        child_parts = []
        for edge in node.edges:
            other = edge.other(node)
            if other is parent:
                continue
            vec, log_scale = self.upward_partials(other, node)
            child_parts.append((self._propagate(vec, edge.length), log_scale))
        (v1, s1), (v2, s2) = child_parts
        cols = []
        log_scale = s1 + s2
        for c in range(self.n_cols):
            a, b = v1[c], v2[c]
            col = [a[0] * b[0], a[1] * b[1], a[2] * b[2], a[3] * b[3]]
            m = max(col)
            if m > 0.0:
                inv = 1.0 / m
                col = [x * inv for x in col]
                log_scale += math.log(m)
            cols.append(col)
        result = (cols, log_scale)
        self._partials_cache[key] = result
        return result

    def grafted_log_likelihood(self, edge, distal_len, pendant_len, query_cols):
        """Log-likelihood of the whole alignment with the query grafted onto
        the given edge.  The edge (a--b, length L) is split: distal_len on
        the side of edge.b, L - distal_len on the side of edge.a; the new
        leaf hangs off the inserted node with branch pendant_len."""
        if distal_len < 0.0 or distal_len > edge.length or pendant_len < 0.0:
            return -math.inf
        part_u, scale_u = self.upward_partials(edge.a, edge.b)
        part_v, scale_v = self.upward_partials(edge.b, edge.a)
        prox_same, prox_diff = jc69_probs(edge.length - distal_len)
        dist_same, dist_diff = jc69_probs(distal_len)
        pend_same, pend_diff = jc69_probs(pendant_len)
        log_l = scale_u + scale_v
        for c in range(self.n_cols):
            au = part_u[c]
            av = part_v[c]
            qv = query_cols[c]
            sum_u = au[0] + au[1] + au[2] + au[3]
            sum_v = av[0] + av[1] + av[2] + av[3]
            sum_q = qv[0] + qv[1] + qv[2] + qv[3]
            site = 0.0
            for s in range(4):
                pu = au[s] * prox_same + (sum_u - au[s]) * prox_diff
                pv = av[s] * dist_same + (sum_v - av[s]) * dist_diff
                pq = qv[s] * pend_same + (sum_q - qv[s]) * pend_diff
                site += pu * pv * pq
            site *= 0.25  # JC69 equilibrium frequencies at the inserted node
            if site <= 0.0:
                return -math.inf
            log_l += math.log(site)
        return log_l
