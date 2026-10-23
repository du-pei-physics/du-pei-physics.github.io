#!/usr/bin/env python3
# Copied from GappedCFT/code/quintuple_check.py (a companion project of the scft-code-search workspace) for the holomorphic-svoa
# note; only pointers to the working note were changed (see tools/package_lg.py).
"""Exact check of the quintuple product identity in the form used in ../../notes/lg_elliptic_genus.md section 3:
    sum_n q^{n(3n+1)/2} (x^{3n} - x^{-3n-1}) = prod_{n>=1} (1-q^n)(1-x q^n)(1-x^{-1} q^{n-1})(1-x^2 q^{2n-1})(1-x^{-2} q^{2n-1}),
as Laurent series in x and power series in q up to q^QM (integers).

usage: quintuple_check.py
"""
# as exact Laurent series in x, power series in q, up to q^QM.
QM = 14
def mul(A, B):
    out = {}
    for (qa, xa), ca in A.items():
        for (qb, xb), cb in B.items():
            q = qa + qb
            if q > QM: continue
            out[(q, xa + xb)] = out.get((q, xa + xb), 0) + ca * cb
    return {k: v for k, v in out.items() if v}
def lin(qe, xe):  # (1 - x^xe q^qe)
    return {(0, 0): 1, (qe, xe): -1} if qe <= QM else {(0, 0): 1}
P = {(0, 0): 1}
for n in range(1, QM + 2):
    for f in (lin(n, 0), lin(n, 1), lin(n - 1, -1), lin(2 * n - 1, 2), lin(2 * n - 1, -2)):
        P = mul(P, f)
L = {}
for n in range(-10, 11):
    e = n * (3 * n + 1) // 2
    if e > QM: continue
    L[(e, 3 * n)] = L.get((e, 3 * n), 0) + 1
    L[(e, -3 * n - 1)] = L.get((e, -3 * n - 1), 0) - 1
L = {k: v for k, v in L.items() if v}
print("quintuple product identity holds up to q^%d:" % QM, P == L, "| terms:", len(P))
