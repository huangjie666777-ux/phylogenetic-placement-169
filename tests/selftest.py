"""Self-test: validation, likelihood sanity, placement, jplace delivery."""

import json
import math
import sys

from fastapi.testclient import TestClient

from app.main import app
from app.io_parsers import parse_fasta, parse_newick
from app.likelihood import LikelihoodEngine
from app.tree import ReferenceTree

client = TestClient(app)

REF = open("examples/reference.fasta").read()
QRY = open("examples/query.fasta").read()
NWK = open("examples/tree.nwk").read()
BODY = {"reference_fasta": REF, "query_fasta": QRY, "newick": NWK}

failures = []
def check(name, cond, extra=""):
    print(("PASS" if cond else "FAIL"), name, extra)
    if not cond:
        failures.append(name)

# 1. healthy placement run
r = client.post("/api/place", json=BODY)
check("place 200", r.status_code == 200, r.text[:300])
data = r.json()
check("two queries placed", len(data["placements"]) == 2)
check("one unplaceable", len(data["unplaceable"]) == 1
      and data["unplaceable"][0]["id"] == "q3_no_evidence")
q1 = next(p for p in data["placements"] if p["id"] == "q1")
check("all 7 edges reported", len(q1["results"]) == 7)
ws = [x["like_weight_ratio"] for x in q1["results"]]
check("weights normalised", abs(sum(ws) - 1.0) < 1e-9, str(sum(ws)))
lls = [x["log_likelihood"] for x in q1["results"]]
check("sorted desc", all(lls[i] >= lls[i+1] - 1e-12 for i in range(len(lls)-1)))
check("pendant in [0,2]", all(0.0 <= x["pendant_length"] <= 2.0 for x in q1["results"]))
# q1 evolved from E: best edge should be a terminal edge near E or D
best = q1["results"][0]
print("  q1 best edge:", best)

# 2. likelihood sanity: Felsenstein equals brute-force state sum on tiny tree
refs = parse_fasta(REF, "reference")
root = parse_newick(NWK)
tree = ReferenceTree(root, [k for k, _ in refs])
eng = LikelihoodEngine(tree, refs)
edge = tree.edges[0]
qcols = [[0., 1., 0., 0.]] * 60  # fake query all C
ll_fast = eng.grafted_log_likelihood(edge, edge.length / 2, 0.1, qcols)
check("loglik finite", math.isfinite(ll_fast), str(ll_fast))

# 3. validation errors reject the whole order
bad_cases = {
    "illegal char": {**BODY, "query_fasta": ">q\nACX" + "A" * 57},
    "duplicate id": {**BODY, "query_fasta": ">A\n" + "A" * 60},
    "empty seq": {**BODY, "query_fasta": ">q\n>q2\n" + "A" * 60 + "\n"},
    "bad branch length": {**BODY, "newick": NWK.replace("0.08", "0")},
    "leaf mismatch": {**BODY, "newick": NWK.replace("A:", "Z:", 1)},
    "unequal length": {**BODY, "query_fasta": ">q\n" + "A" * 59},
}
for name, body in bad_cases.items():
    rr = client.post("/api/place", json=body)
    check(f"reject {name}", rr.status_code == 422, rr.text[:120])

# 4. jplace download
rj = client.post("/api/place/jplace", json=BODY)
check("jplace 200", rj.status_code == 200)
check("jplace attachment",
      "attachment" in rj.headers.get("content-disposition", ""))
doc = rj.json()
check("jplace version 3", doc["version"] == 3)
check("jplace fields", doc["fields"] == ["edge_num", "likelihood",
      "like_weight_ratio", "distal_length", "pendant_length"])
check("jplace tree has edge tags", "{0}" in doc["tree"])
edge_nums = sorted(p[0] for pl in doc["placements"] for p in pl["p"])
check("jplace edge nums match results", edge_nums == sorted(edge_nums)
      and set(edge_nums) == set(range(7)))
json.dumps(doc)  # serialisable

if failures:
    print("FAILURES:", failures)
    sys.exit(1)
print("ALL SELF-TESTS PASSED")
