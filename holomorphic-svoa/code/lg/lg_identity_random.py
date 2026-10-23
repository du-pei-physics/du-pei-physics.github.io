#!/usr/bin/env python3
# Copied from GappedCFT/code/lg_identity_random.py (a companion project of the scft-code-search workspace) for the holomorphic-svoa
# note; only pointers to the working note were changed (see tools/package_lg.py).
"""Test I_LG = Index(C)^2 on random self-dual ternary codes of length 12 and 24 (random neighbor walks).
A neighbor of C: pick an isotropic v not in C (v.v = 0), C' = (C cap v^perp) + <v>. Seeded, exact."""
import json, random, sys
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
from lg_count_dense import span, rref3, rank3, is_self_dual, analyse, paley_S

def neighbor(G, rng):
    k, n = G.shape
    C = span(G)
    while True:
        v = np.array([rng.randrange(3) for _ in range(n)], dtype=np.int64)
        if (v @ v) % 3 or not v.any(): continue
        if rank3(np.vstack([G, v])) == k: continue           # v in C
        break
    D = C[(C @ v) % 3 == 0]
    B = rref3(D); B = B[(B != 0).any(1)]
    H = rref3(np.vstack([B, v])); H = H[(H != 0).any(1)]
    assert H.shape[0] == k and is_self_dual(H)
    return H

rng = random.Random(20261021)
golay = np.array([[1 if i == r else 0 for i in range(6)] + [int(ch) for ch in s] for r, s in
                  enumerate(["011111", "101221", "110122", "121012", "122101", "112210"])], dtype=np.int64)
P24 = np.concatenate([np.eye(12, dtype=np.int64), paley_S()], axis=1)
results = []
for start, steps, reps in [(golay, 3, 15), (P24, 4, 15)]:
    for r in range(reps):
        G = start
        for _ in range(steps):
            G = neighbor(G, rng)
        res = analyse(f"n{G.shape[1]}_walk{r}", G); results.append(res)
        print(json.dumps({k: res[k] for k in ("code", "d", "full_support_words", "Index(C)", "I_LG", "match")}))
print("all match:", all(r["match"] for r in results), "distinct (Index, I_LG):", sorted({(r["Index(C)"], r["I_LG"]) for r in results}))
