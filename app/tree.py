"""Unrooted binary reference tree: structure, validation, edge numbering."""

from __future__ import annotations

from collections import deque

from .schemas import InputError


class Node:
    __slots__ = ("nid", "name", "edges")

    def __init__(self, nid, name=None):
        self.nid = nid
        self.name = name          # leaf name, None for internal nodes
        self.edges = []           # list of Edge

    @property
    def is_leaf(self):
        return self.name is not None


class Edge:
    __slots__ = ("eid", "a", "b", "length")

    def __init__(self, eid, a, b, length):
        self.eid = eid
        self.a = a
        self.b = b
        self.length = length

    def other(self, node):
        return self.b if node is self.a else self.a


class ReferenceTree:
    """Fixed unrooted binary tree with original branch lengths.

    A rooted-binary Newick (root of degree 2) is accepted: the degree-2
    root is suppressed and its two child branches merged, preserving the
    unrooted topology and total branch lengths exactly.
    """

    def __init__(self, root_clade, ref_ids):
        self.nodes = []
        self.edges = []
        self._build(root_clade)
        self._validate_binary()
        self._validate_leaves(ref_ids)
        self._number_edges()
        self._by_name = {n.name: n for n in self.nodes if n.is_leaf}

    # -- construction ----------------------------------------------------
    def _make_node(self, clade):
        node = Node(len(self.nodes), clade.name if clade.is_terminal() else None)
        self.nodes.append(node)
        return node

    def _connect(self, clade, parent_node, length):
        node = self._make_node(clade)
        edge = Edge(len(self.edges), parent_node, node, float(length))
        self.edges.append(edge)
        parent_node.edges.append(edge)
        node.edges.append(edge)
        for child in clade.clades:
            self._connect(child, node, child.branch_length)
        return node

    def _build(self, root_clade):
        if len(root_clade.clades) == 2:
            # suppress the degree-2 root, merging both child branch lengths
            c1, c2 = root_clade.clades
            n1 = self._make_node(c1)
            n2 = self._make_node(c2)
            merged = float(c1.branch_length) + float(c2.branch_length)
            edge = Edge(len(self.edges), n1, n2, merged)
            self.edges.append(edge)
            n1.edges.append(edge)
            n2.edges.append(edge)
            for child in c1.clades:
                self._connect(child, n1, child.branch_length)
            for child in c2.clades:
                self._connect(child, n2, child.branch_length)
        else:
            root = self._make_node(root_clade)
            for child in root_clade.clades:
                self._connect(child, root, child.branch_length)

    def _validate_binary(self):
        n_leaves = sum(1 for n in self.nodes if n.is_leaf)
        for node in self.nodes:
            degree = len(node.edges)
            if node.is_leaf:
                if degree != 1:
                    raise InputError(f"newick: leaf '{node.name}' has degree {degree}")
            elif degree != 3:
                raise InputError(
                    "newick: tree must be an unrooted binary tree "
                    "(every internal node degree 3); found internal node "
                    f"with degree {degree}"
                )
        expected_edges = 2 * n_leaves - 3
        if len(self.edges) != expected_edges:
            raise InputError(
                f"newick: expected {expected_edges} branches for {n_leaves} "
                f"leaves in an unrooted binary tree, got {len(self.edges)}"
            )

    def _validate_leaves(self, ref_ids):
        leaf_names = sorted(n.name for n in self.nodes if n.is_leaf)
        if leaf_names != sorted(ref_ids):
            missing = sorted(set(ref_ids) - set(leaf_names))
            extra = sorted(set(leaf_names) - set(ref_ids))
            raise InputError(
                f"newick: tree leaves do not match reference ids; "
                f"missing on tree: {missing or '[]'}, not in alignment: {extra or '[]'}"
            )

    def _number_edges(self):
        """Stable edge numbers: BFS from the leaf with the smallest id."""
        start = min((n for n in self.nodes if n.is_leaf), key=lambda n: n.name)
        order = []
        seen_nodes = {start.nid}
        seen_edges = set()
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for edge in sorted(node.edges, key=lambda e: e.eid):
                if edge.eid in seen_edges:
                    continue
                seen_edges.add(edge.eid)
                order.append(edge)
                other = edge.other(node)
                if other.nid not in seen_nodes:
                    seen_nodes.add(other.nid)
                    queue.append(other)
        for new_eid, edge in enumerate(order):
            edge.eid = new_eid
        self.edges.sort(key=lambda e: e.eid)

    # -- accessors -------------------------------------------------------
    def leaf(self, name):
        return self._by_name[name]

    def to_newick(self, annotate_edges=True):
        """Serialize the unrooted tree (trifurcating root at the internal
        node with the smallest id), with {edge_num} tags on each branch."""
        root = min((n for n in self.nodes if not n.is_leaf), key=lambda n: n.nid)

        def render(node, parent):
            if node.is_leaf:
                out = node.name
            else:
                parts = []
                for edge in node.edges:
                    other = edge.other(node)
                    if other is parent:
                        continue
                    parts.append(render(other, node))
                out = "(" + ",".join(parts) + ")"
            if parent is not None:
                edge = next(e for e in node.edges if e.other(node) is parent)
                out += f":{edge.length:.10g}"
                if annotate_edges:
                    out += "{" + str(edge.eid) + "}"
            return out

        return render(root, None) + ";"
