#!/usr/bin/env python3
# Copied from GappedCFT/code/lg_count_dense.py (a companion project of the scft-code-search workspace) for the holomorphic-svoa
# note; only pointers to the working note were changed (see tools/package_lg.py).
"""I_LG = sum_{x in C full support} sum_{c in C, c <= x} (-1)^{wt c} versus Index(C)^2, by enumeration (n <= 24).

Codes: the Pless symmetry code P24 = (I_12 | S), with S the skew Paley conference matrix of order 12 (the paper's S),
and the extended ternary quadratic-residue code of length 24. Both are extremal [24, 12, 9]. Exact integer arithmetic.
"""
import itertools, json
import numpy as np

def paley_S():
    q = 11
    chi = lambda a: 0 if a % q == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)
    S = np.zeros((12, 12), dtype=np.int64)            # index 0 = infinity, 1..11 = residues 1..10, 0
    pts = [None] + list(range(1, 11)) + [0]
    for i, x in enumerate(pts):
        for j, y in enumerate(pts):
            if x is None and y is None: S[i, j] = 0
            elif x is None: S[i, j] = 1
            elif y is None: S[i, j] = -1
            else: S[i, j] = chi(y - x)
    return S % 3

def qr24():
    p = 23
    Q = {pow(a, 2, p) for a in range(1, p)}
    # generator of the length-23 QR code: idempotent-free route via the cyclic code spanned by shifts of the
    # indicator-based word; simplest: g(x) = prod_{r in Q} (x - alpha^r) over F_3^11 -- use the polynomial
    # factorization of x^23 - 1 over F_3 by brute-force search of the degree-11 factor whose roots are alpha^Q.
    # Instead, build the code as the span of the cyclic shifts of e = sum_{r in N} x^r (+ constant) and test.
    N = set(range(1, p)) - Q
    for const in range(3):
        for base in (Q, N):
            w = np.zeros(p, dtype=np.int64); w[list(base)] = 1; w[0] = const
            G = np.array([np.roll(w, s) for s in range(p)]) % 3
            r = rank3(G)
            if r == 12:
                G = rref3(G)[:12]
                ext = np.concatenate([G, (-G.sum(1, keepdims=True)) % 3], axis=1)
                C = span(ext)
                if is_self_dual(ext) and min_wt(C) == 9:
                    return ext
    raise RuntimeError("QR construction failed")

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
    return not ((G @ G.T) % 3).any() and 2 * G.shape[0] == G.shape[1]

def min_wt(C):
    w = (C != 0).sum(1); return int(w[w > 0].min())

def analyse(name, G):
    C = span(G)
    n = G.shape[1]
    full = C[(C != 0).all(1)]
    index = int(((-1) ** (full == 2).sum(1)).sum())
    bits = 1 << np.arange(n, dtype=np.int64)
    ones = ((C == 1) * bits).sum(1); twos = ((C == 2) * bits).sum(1)
    sgn = (-1) ** (C != 0).sum(1)
    total = 0
    for x in full:
        xo = int(((x == 1) * bits).sum()); xt = int(((x == 2) * bits).sum())
        ok = ((ones & xt) == 0) & ((twos & xo) == 0)                      # c_i in {0, x_i}
        total += int(sgn[ok].sum())
    return {"code": name, "n": n, "self_dual": bool(is_self_dual(G)), "d": min_wt(C), "full_support_words": int(len(full)),
            "Index(C)": index, "Index^2": index * index, "I_LG": total, "match": total == index * index}

if __name__ == "__main__":
    S = paley_S()
    P24 = np.concatenate([np.eye(12, dtype=np.int64), S], axis=1)
    for name, G in [("Pless P24", P24), ("extended QR24", qr24())]:
        print(json.dumps(analyse(name, G)))
