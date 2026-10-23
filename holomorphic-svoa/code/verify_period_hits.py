# Copied from claude/kspan/verify_period_hits.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Verify period-ring hits of mixed_products.py (--periods 3) for the order-3 / order-6 rows and write a
certificate.  Uses persistent/mixed_rings/run.py:verify (rebuilt generator, self-duality, rank, generator-only
banded DP for the index, GJF shift g, fixed subcode of h = g^{m/3}, contracted one-period code and its
balance by exhaustive enumeration when dim <= 30, else by the banded DP).

usage: python verify_period_hits.py results/kspan/mixed_k15_run1.json results/kspan/k15_certificate.json
"""
import json, sys, time
from math import gcd
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mixed_rings import verify

src, dst = sys.argv[1], sys.argv[2]
data = json.load(open(src)); r = int(data.get("r", 3)); periods = int(data["periods"]); k = int(data["k"]); m = gcd(k, 24)
out = []
for i, h in enumerate(data["hits"]):
    t = time.time()
    rep, G = verify(h["cells"] * periods, r, h["kind"] == "cyc", periods, m)
    rep["seconds"] = round(time.time() - t, 1)
    print(i, h["kind"], h["index"], json.dumps({kk: v for kk, v in rep.items() if kk != "balance_exhaustive"}), rep.get("balance_exhaustive"), flush=True)
    out.append({"hit": i, "kind": h["kind"], "widths": h["widths"], "periods": periods, "r": r, "claimed_index": h["index"],
                "cells": h["cells"], "generator": G.tolist(), "verification": rep})
ok = [o for o in out if o["verification"]["ok"]]
Path(dst).write_text(json.dumps({
    "what": f"self-dual ternary codes of length {12*k} (c = {12*k}, k = {k}) with index +-24, a monomial automorphism g of order {m} "
            f"(shift of the period ring) whose element h = g^(m/3) has cycle shape 3^({4*k}) and balanced fixed maximal words; "
            f"GJF Theorem 5.6 gives orbifold index 24/{m} = {24//m}, the TMF minimum at c = {12*k}",
    "construction": "mixed-width period ring: one period = cells of widths " + ",".join(map(str, data["widths"])) +
                    f" (common delay count r = {r}), repeated {periods} times; code = {{(a | aA)}}, A block-bidiagonal; "
                    "negacyclic: g = z^2 (shift by two periods with wrap sign), cyclic: g = z",
    "found_by": f"kspan/mixed_products.py --periods {periods} ({data['tried']:,} products)",
    "verified_by": "kspan/verify_period_hits.py via persistent/mixed_rings/run.py:verify",
    "n_hits": len(out), "n_verified": len(ok), "hits": out}, indent=0))
print("verified", len(ok), "of", len(out), "->", dst)
