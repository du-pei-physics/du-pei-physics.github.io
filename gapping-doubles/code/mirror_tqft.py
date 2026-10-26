#!/usr/bin/env python3
"""Vacuum bookkeeping of the massive phase of (sum Phi_i^3)/C under the duality line N (Section 4 of the note).

Per factor the (1,1) Morse deformation h = Re z^3 + (3/2)|z|^2 has one Z3-fixed minimum (even) and three Z3-permuted
saddles (odd): the shape of the regular TY(Z3) module, with N exchanging the fixed vacuum and the permuted ones.
After the C-orbifold the vacua are labelled by S = set of factors in the minimum, a C-orbit of saddle labels on S^c
(3^{n-|S|} / |C / C_S| orbits, C_S = words supported in S) and a character of the stabilizer C_S; parity (-1)^{n-|S|}.
  signed count  I(C) = sum_S (-3)^{n-|S|} |C_S|^2 / |C|           (= Index(C)^2; holomorphic-svoa note, Sec. L)
  N maps S-type vacua to S^c-type vacua.  For a self-dual code this is a bijection of the counts,
  3^{n-|S|}|C_S|^2 = 3^{|S|}|C_{S^c}|^2, because dim C_{S^c} - dim C_S = n/2 - |S|.
  The one-sided N-gauging is consistent only for 8 | n (the chiral -1 of the code theory has twisted energy n/16);
  for the examples here (n = 4, 12) the halved numbers are formal complement-orbit counts, not vacuum counts of a
  gauged theory.
Exact integers; enumerates all 2^n subsets (n <= 16 here: tetracode, Golay); the identities are proved for all n.
usage: mirror_tqft.py
"""
import itertools
import numpy as np


def rref3(M):
    M = M.copy() % 3; r = 0
    for c in range(M.shape[1]):
        piv = next((i for i in range(r, M.shape[0]) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]; M[r] = (M[r] * M[r, c]) % 3          # inverse of 1 is 1, of 2 is 2
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % 3
        r += 1
        if r == M.shape[0]: break
    return M


def rank3(M):
    R = rref3(M); return int((R != 0).any(1).sum())


def span(G):
    k = G.shape[0]
    co = np.array(list(itertools.product(range(3), repeat=k)), dtype=np.int64)
    return (co @ G) % 3


def is_self_dual(G):
    return not ((G @ G.T) % 3).any() and 2 * G.shape[0] == G.shape[1] and rank3(G) == G.shape[0]


def min_wt(C):
    w = (C != 0).sum(1); return int(w[w > 0].min())


def golay():
    """Extended ternary Golay code as the extended QR code of length 11 + 1."""
    p = 11
    Q = {pow(a, 2, p) for a in range(1, p)}
    N = set(range(1, p)) - Q
    for const in range(3):
        for base in (Q, N):
            w = np.zeros(p, dtype=np.int64); w[list(base)] = 1; w[0] = const
            G = np.array([np.roll(w, s) for s in range(p)]) % 3
            if rank3(G) == 6:
                G = rref3(G)[:6]
                ext = np.concatenate([G, (-G.sum(1, keepdims=True)) % 3], axis=1)
                if is_self_dual(ext) and min_wt(span(ext)) == 6:
                    return ext
    raise RuntimeError


def tetracode():
    return np.array([[1, 0, 1, 1], [0, 1, 1, 2]], dtype=np.int64)


def signed_index(C):
    """Index(C) = sum over full-support words of prod_i (+1 if c_i = 1, -1 if c_i = 2)."""
    tot = 0
    for c in C:
        if np.all(c != 0): tot += int(np.prod(np.where(c == 1, 1, -1)))
    return tot


def analyse(name, G):
    n = G.shape[1]; k = G.shape[0]
    C = span(G)
    assert is_self_dual(G) and len(C) == 3 ** k
    idx = signed_index(C)
    supp = [frozenset(np.nonzero(c)[0].tolist()) for c in C]
    dimCS = {}
    for r in range(n + 1):
        for S in itertools.combinations(range(n), r):
            Sf = frozenset(S)
            cnt = sum(1 for s in supp if s <= Sf)        # |C_S|
            dimCS[Sf] = cnt
    I = 0
    bij = True
    for Sf, cnt in dimCS.items():
        I += (-3) ** (n - len(Sf)) * cnt * cnt
        Sc = frozenset(range(n)) - Sf
        if 3 ** (n - len(Sf)) * cnt ** 2 != 3 ** len(Sf) * dimCS[Sc] ** 2: bij = False
    assert I % (3 ** k) == 0
    I //= 3 ** k
    nS = {len(Sf): 0 for Sf in dimCS}
    for Sf, cnt in dimCS.items(): nS[len(Sf)] += 3 ** (n - len(Sf)) * cnt * cnt   # |C| x number of S-type vacua
    gauged = "one-sided gauging consistent (8 | n)" if n % 8 == 0 else "one-sided gauging anomalous (8 does not divide n): formal count"
    print(f"{name}: n = {n}, Index(C) = {idx}, signed vacuum count I(C) = {I} (Index^2 = {idx * idx}), "
          f"N-bijection S <-> S^c for all S: {bij}; complement-orbit signed count I/2 = {I // 2 if I % 2 == 0 else I / 2} [{gauged}]")
    print("   |C| x (vacua with |S| = r), r = 0..n:", [nS[r] for r in range(n + 1)])
    tot = sum(nS.values()); assert tot % (3 ** k) == 0
    print(f"   number of vacua (unsigned) = {tot // 3 ** k}; complement-orbit count {tot // 3 ** k // 2}")
    return I == idx * idx and bij


if __name__ == "__main__":
    ok = analyse("tetracode", tetracode()) and analyse("Golay", golay())
    print("all checks pass:", ok)
