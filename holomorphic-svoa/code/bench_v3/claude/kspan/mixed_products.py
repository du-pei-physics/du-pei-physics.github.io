"""Product rings with cells of mixed widths (12 and 6) and a common number r of delayed
channels, for the rows that need no symmetry (m = 1, 2) or only a period structure (m = 3, 6).

Legacy cell form per cell j: Q_j orthogonal (q_j x q_j), S_j anti-orthogonal, the first r outputs
of a Q_j are delayed into cell j+1; the check block j+1 is (p, u_{j+1}) S_{j+1}.  Tables are
3^r x 3^r for every width, so cells of different widths chain freely (this is how Codex's k = 11
code is built: five 12-cells and one 6-cell).  Index = Tr(J T_1 ... T_N) (negacyclic) or
Tr(T_1 ... T_N) (cyclic); for a period-p ring with 3 periods the GJF conditions are
Tr(J M) = 0, |Tr(J M^3)| = 24 with M the period product.

usage: python mixed_products.py --k 13 --widths 12,12,12,12,12,12,6 --minutes 10 --seed 1 --out results/kspan/mixed_k13.json
       python mixed_products.py --k 15 --widths 12,12,6 --periods 3 ...   (balanced period ring)
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cells_general import table, anti_bases, hadamard_Q, random_reflection, is_orth, is_anti


def random_cell(rng, q, r, bases, metric2=0.0, s_refl=0):
    """Random legacy-form cell.  metric2 (0/1): change the delay space to one of
    discriminant 2 (Codex's 'second discriminant': Q -> Q R, S -> R^{-1} S with R mixing one delayed
    and one kept coordinate by [[1,1],[1,2]]; R^T R is block-diagonal for the split, so A A^T = -I
    still holds although Q, S are no longer (anti-)orthogonal).  s_refl: up to this many random
    reflections applied to S on the right (keeps S anti-orthogonal)."""
    name = list(bases)[rng.integers(len(bases))]
    S = bases[name]
    if rng.random() < 0.5:
        perm = rng.permutation(q); sg = rng.integers(1, 3, q); S = (S[:, perm] * sg[None, :]) % 3
    for _ in range(int(rng.integers(0, s_refl + 1))):
        S = (S @ random_reflection(rng, q)) % 3
    Q = hadamard_Q(q, rng)
    for _ in range(int(rng.integers(0, 3))):
        Q = (Q @ random_reflection(rng, q)) % 3
    assert is_orth(Q) and is_anti(S)
    desc = {"q": q, "base": name}
    if metric2:
        # the delay-space Gram must be the SAME for every cell of a ring (the check block of cell j+1
        # receives cell j's delayed outputs), so metric 2 is a ring-level choice: every cell gets the
        # norm-2 direction at delayed index 0 (metric2 = 1), with its own kept partner j
        i, j = 0, int(rng.integers(r, q))
        R = np.eye(q, dtype=np.int64); R[np.ix_([i, j], [i, j])] = [[1, 1], [1, 2]]
        Ri = np.eye(q, dtype=np.int64); Ri[np.ix_([i, j], [i, j])] = [[2, 2], [2, 1]]
        assert not ((R @ Ri) % 3 - np.eye(q, dtype=np.int64)).any()
        Q = (Q @ R) % 3; S = (Ri @ S) % 3; desc["metric2"] = [i, j]
    T, neg = table(Q, S, r)
    desc.update({"Q": Q.tolist(), "S": S.tolist()})
    return T, neg, desc


def trJ(M, neg):
    return int(sum(int(M[neg[i], i]) for i in range(M.shape[0])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True); ap.add_argument("--widths", required=True)
    ap.add_argument("--periods", type=int, default=1); ap.add_argument("--r", type=int, default=3)
    ap.add_argument("--pool", type=int, default=40); ap.add_argument("--minutes", type=float, default=10)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--out", required=True)
    ap.add_argument("--refresh", type=int, default=0, help="regenerate the cell pools every this many products (0 = never)")
    ap.add_argument("--metric2", type=int, default=0, help="1 = all cells use a discriminant-2 delay space"); ap.add_argument("--s-refl", type=int, default=0)
    a = ap.parse_args()
    widths = [int(w) for w in a.widths.split(",")]
    assert sum(widths) * a.periods == 6 * a.k, "sum of widths times periods must be 6k"
    rng = np.random.default_rng(a.seed)
    pools = {q: [] for q in set(widths)}
    bases = {q: anti_bases(q) for q in set(widths)}
    for q in pools:
        for _ in range(a.pool):
            pools[q].append(random_cell(rng, q, a.r, bases[q], a.metric2, a.s_refl))
    neg = pools[widths[0]][0][1]
    n = 3 ** a.r
    t_end = time.time() + 60 * a.minutes
    tried = 0; hits = []; small = {}; v2min = 99; balanced = 0
    while time.time() < t_end:
        if a.refresh and tried and tried % a.refresh == 0:
            pools = {q: [random_cell(rng, q, a.r, bases[q], a.metric2, a.s_refl) for _ in range(a.pool)] for q in pools}
        picks = [(q, int(rng.integers(len(pools[q])))) for q in widths]
        M = np.eye(n, dtype=object)
        for q, i in picks:
            M = M.dot(pools[q][i][0])
        tried += 1
        if a.periods == 1:
            vals = {"nega": trJ(M, neg), "cyc": int(np.trace(M))}
        else:
            M3 = M.dot(M).dot(M)
            vals = {}
            if trJ(M, neg) == 0:
                vals["nega"] = trJ(M3, neg); balanced += 1
            if int(np.trace(M)) == 0:
                vals["cyc"] = int(np.trace(M3))
        for kind, v in vals.items():
            av = abs(v)
            if av and av <= 240:
                small[str(av)] = small.get(str(av), 0) + 1
            if av == 24:
                hits.append({"k": a.k, "kind": kind, "widths": widths, "periods": a.periods, "index": v,
                             "cells": [pools[q][i][2] for q, i in picks]})
                Path(a.out).write_text(json.dumps({"k": a.k, "tried": tried, "hits": hits, "small": small}, indent=1))
    Path(a.out).write_text(json.dumps({"k": a.k, "widths": widths, "periods": a.periods, "tried": tried, "balanced": balanced,
                                       "hits": hits, "small_values": dict(sorted(small.items(), key=lambda kv: int(kv[0])))}, indent=1))
    print(json.dumps({"k": a.k, "tried": tried, "balanced": balanced, "hits": len(hits), "small": dict(sorted(small.items(), key=lambda kv: int(kv[0]))[:8])}))


if __name__ == "__main__":
    main()
