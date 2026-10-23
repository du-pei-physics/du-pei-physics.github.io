#!/usr/bin/env python3
# Copied from GappedCFT/code/pointwise_identity.py (a companion project of the scft-code-search workspace) for the holomorphic-svoa
# note; only pointers to the working note were changed (see tools/package_lg.py).
"""Pointwise and charge-resolved form of the index identity (../../notes/lg_elliptic_genus.md section 2).

For a self-dual ternary code C of length n and a full-support codeword w, let
    N(w) = sum_{c in C, c <= w} (-1)^{n - wt(c)},      c <= w meaning c_i in {0, w_i},
    sign(w) = (-1)^{#{i : w_i = 2}},  and  w^ in {+1,-1}^n the lift of w (1 -> +1, 2 -> -1).
Claims checked here (exact integers):
  (1) N(w) = sign(w) Index(C) for every full-support w            (pointwise identity; proof in the public note);
  (2) sum_w s^{sum w^} N(w) = Index(C) sum_w sign(w) s^{-sum w^}  (the q^0 term of the elliptic genus of LG/C with
      left U(1) fugacity s = y^{1/6}, against Index(C) times the Ramond ground-state character of V_Lambda(C)).
Codes: the ternary Golay code, the Pless symmetry code P24 and the extended QR code of length 24 (as in
lg_count_dense.py), and 20 random self-dual codes of length 12 and 24 from the neighbor walk of lg_identity_random.py.

usage: pointwise_identity.py
"""
import os, sys, itertools, random
from collections import Counter
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from lg_count_dense import paley_S, qr24, is_self_dual
from lg_identity_random import neighbor


def codewords(G):
    G = np.array(G, dtype=np.int64) % 3
    k, n = G.shape
    msgs = np.array(list(itertools.product(range(3), repeat=k)), dtype=np.int64)
    return (msgs @ G) % 3

def check(name, G):
    C = codewords(G); n = C.shape[1]
    full = C[(C != 0).all(axis=1)]
    sgn = lambda W: (-1) ** (W == 2).sum(axis=-1)
    index = int(sgn(full).sum())
    wt = (C != 0).sum(axis=1)
    lg, cft, ok = Counter(), Counter(), 0
    for w in full:
        le = ((C == 0) | (C == w)).all(axis=1)              # c <= w
        N = int(((-1) ** (n - wt[le])).sum())
        q = int(np.where(w == 1, 1, -1).sum())
        lg[q] += N; cft[-q] += index * int(sgn(w))
        ok += (N == index * int(sgn(w)))
    return n, len(full), index, ok, lg == cft

if __name__ == "__main__":
    golay = [[1 if i == r else 0 for i in range(6)] + [int(ch) for ch in s] for r, s in
             enumerate(["011111", "101221", "110122", "121012", "122101", "112210"])]
    P24 = np.concatenate([np.eye(12, dtype=np.int64), paley_S()], axis=1)
    codes = [("Golay", np.array(golay, dtype=np.int64)), ("Pless P24", P24), ("extended QR24", qr24())]
    rng = random.Random(20261022)
    for start, steps, reps in [(np.array(golay, dtype=np.int64), 3, 10), (P24, 4, 10)]:
        for r in range(reps):
            G = start
            for _ in range(steps): G = neighbor(G, rng)
            codes.append((f"walk n{G.shape[1]} #{r}", G))
    allok = True
    for name, G in codes:
        assert is_self_dual(G)
        n, nf, index, ok, charge = check(name, G)
        allok &= (ok == nf and charge)
        print(f"{name}: n = {n}, full-support words {nf}, Index {index}; pointwise identity {ok}/{nf}; "
              f"q^0 charge-resolved match: {charge}")
    print("all codes pass:", allok)
