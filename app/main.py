"""FastAPI entrypoint: validate, place, and deliver results."""

from __future__ import annotations

import json

from fastapi import FastAPI
from fastapi.responses import JSONResponse, Response
from pydantic import BaseModel, Field

from . import io_parsers
from .jplace import build_jplace
from .likelihood import LikelihoodEngine
from .optimize import place_query
from .schemas import (InputError, MAX_QUERIES, MAX_REFS, MIN_QUERIES,
                      MIN_REFS, check_unique_ids, validate_alignment)
from .tree import ReferenceTree

app = FastAPI(title="eDNA fixed-tree placement backend", version="1.0.0")


class PlacementRequest(BaseModel):
    reference_fasta: str = Field(..., description="Aligned reference FASTA")
    query_fasta: str = Field(..., description="Aligned query FASTA")
    newick: str = Field(..., description="Unrooted binary Newick with branch lengths")


def _run_placement(req: PlacementRequest):
    refs = io_parsers.parse_fasta(req.reference_fasta, "reference")
    queries = io_parsers.parse_fasta(req.query_fasta, "query")
    n_cols = validate_alignment(refs, "reference", MIN_REFS, MAX_REFS)
    q_cols = validate_alignment(queries, "query", MIN_QUERIES, MAX_QUERIES)
    if q_cols != n_cols:
        raise InputError(
            f"query alignment has {q_cols} columns but reference has "
            f"{n_cols}; all sequences must share the same alignment length"
        )
    check_unique_ids(refs, queries)
    root_clade = io_parsers.parse_newick(req.newick)
    tree = ReferenceTree(root_clade, [seq_id for seq_id, _ in refs])
    engine = LikelihoodEngine(tree, refs)

    placed, failed = [], []
    for seq_id, seq in queries:
        results, reason = place_query(engine, seq_id, seq)
        if results is None:
            failed.append({"id": seq_id, "reason": reason})
        else:
            placed.append((seq_id, results))
    return tree, n_cols, placed, failed


@app.exception_handler(InputError)
async def input_error_handler(_, exc: InputError):
    return JSONResponse(status_code=422, content={"detail": str(exc)})


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.post("/api/place")
def place(req: PlacementRequest):
    tree, n_cols, placed, failed = _run_placement(req)
    return {
        "n_reference": len(tree._by_name),
        "n_columns": n_cols,
        "n_edges": len(tree.edges),
        "placements": [
            {
                "id": seq_id,
                "results": [
                    {
                        "edge_num": r.edge_id,
                        "log_likelihood": r.log_likelihood,
                        "like_weight_ratio": r.weight,
                        "distal_length": r.distal_length,
                        "pendant_length": r.pendant_length,
                    }
                    for r in results
                ],
            }
            for seq_id, results in placed
        ],
        "unplaceable": failed,
        "weight_note": ("like_weight_ratio is relative support normalised "
                        "across all branches, not a probability of correct "
                        "species assignment."),
    }


@app.post("/api/place/jplace")
def place_jplace(req: PlacementRequest):
    tree, _, placed, failed = _run_placement(req)
    doc = build_jplace(tree, placed)
    if failed:
        doc["metadata"]["unplaceable"] = failed
    payload = json.dumps(doc, indent=2)
    return Response(
        content=payload,
        media_type="application/json",
        headers={"Content-Disposition": "attachment; filename=placements.jplace"},
    )
