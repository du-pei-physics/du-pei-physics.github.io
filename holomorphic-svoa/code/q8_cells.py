# Copied from claude/k8/q8_cells.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Quaternion-equivariant transfer cells: a route to k = 8 (c = 96, index 3) that avoids Z_8.

Idea.  GJF's index formula for a (non-anomalous) orbifold by a finite group G of code
automorphisms is Index(V//G) = (1/|G|) sum over commuting pairs (a,b) of Z(a,b).  If every
commuting pair generates a CYCLIC subgroup, every Z(a,b) is a modular image of a twining trace
Z(1,a') = signed number of a'-fixed maximal words, and these vanish whenever -1 is a power of a'.
In Q_8 = {+-1, +-i, +-j, +-k} every commuting pair lies in a cyclic subgroup and every element
other than 1 has -1 as a power.  So a self-dual code of length 96 with Index 24 and a monomial
Q_8 action (i^2 = j^2 = k^2 = -1 = the scalar -I) would give orbifold index 24/8 = 3 -- without
any element of order 8, hence outside the 16-divisibility phenomenon of Z_8 actions.

Monomial Q_8 modules over F_3 with -1 acting as -I: the 4-dimensional one, Ind_{<-1>}^{Q_8}(sign),
coordinates (1, i, j, k), left multiplication; it is irr (+) irr for the 2-dimensional irreducible
rep (Q_8 < SL(2,3)).  A 12-coordinate cell is three such blocks: V = irr (x) F_3^6, and the
Q_8-commutant is M_6(F_3).  The standard form on V is omega (x) Omega with Omega symplectic on
F_3^6, so an equivariant anti-orthogonal S is 1 (x) s with Omega(ms, m's) = -Omega(m, m')
(a conformal symplectic map of multiplier -1), an equivariant orthogonal Q is 1 (x) q with
q in Sp(6,3), and an invariant delay space is irr (x) Y with Y a symplectic plane (r = 4).

This file builds such cells, their local tables (projection form, 81 states), checks everything
on expanded generators, and surveys the 2-adic valuation of ring indices (identical cells and
4-cell products), looking for |index| = 24.
usage: python q8_cells.py --selftest
       python q8_cells.py --survey --minutes 10 --seed 1 --out results/k8/q8_survey.json
"""
import argparse, itertools, json, sys, time
from math import gcd
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from f3lib import rref, rank, nullspace, is_self_dual, preserves, monomial_matrix
from cells_general import traces

# ---------------------------------------------------------------- Q_8 on four coordinates (1, i, j, k)
def signed_perm(images):
    """images[a] = (b, sign): e_a -> sign * e_b  (row-vector convention: x -> x M)."""
    M = np.zeros((4, 4), dtype=np.int64)
    for a, (b, s) in enumerate(images):
        M[a, b] = s % 3
    return M

SI = signed_perm([(1, 1), (0, -1), (3, 1), (2, -1)])     # i: 1->i, i->-1, j->k, k->-j
SJ = signed_perm([(2, 1), (3, -1), (0, -1), (1, 1)])     # j: 1->j, i->-k, j->-1, k->i
SK = (SI @ SJ) % 3                                       # k = i j
I4 = np.eye(4, dtype=np.int64)
assert np.array_equal((SI @ SI) % 3, (2 * I4) % 3) and np.array_equal((SJ @ SJ) % 3, (2 * I4) % 3)
assert np.array_equal((SK @ SK) % 3, (2 * I4) % 3) and np.array_equal((SJ @ SI) % 3, (2 * SK) % 3)


def block_diag(*ms):
    n = sum(m.shape[0] for m in ms); out = np.zeros((n, n), dtype=np.int64); i = 0
    for m in ms:
        out[i:i + m.shape[0], i:i + m.shape[0]] = m; i += m.shape[0]
    return out


Q = 12
GI, GJ, GK = block_diag(SI, SI, SI), block_diag(SJ, SJ, SJ), block_diag(SK, SK, SK)
GENS = [GI, GJ]


def commutant_basis():
    """Basis of {M : G M = M G for G in GENS} (12x12 matrices as 144-vectors)."""
    rows = []
    for G in GENS:
        # (G M - M G)_{ab} = sum_c G_ac M_cb - M_ac G_cb : linear in vec(M)
        A = np.zeros((Q * Q, Q * Q), dtype=np.int64)
        for a in range(Q):
            for b in range(Q):
                for c in range(Q):
                    A[a * Q + b, c * Q + b] += G[a, c]
                    A[a * Q + b, a * Q + c] -= G[c, b]
        rows.append(A % 3)
    return nullspace(np.vstack(rows) % 3)           # each row = vec(M)


COMM = commutant_basis()
assert COMM.shape[0] == 36, COMM.shape


def rand_commutant(rng):
    return ((rng.integers(0, 3, 36) @ COMM) % 3).reshape(Q, Q)


def is_anti(S):
    return np.array_equal((S @ S.T) % 3, (2 * np.eye(Q, dtype=np.int64)) % 3)


def is_orth(M):
    return np.array_equal((M @ M.T) % 3, np.eye(Q, dtype=np.int64))


def random_anti_equivariant(rng, seed_S):
    """S = seed_S . (random equivariant orthogonal): anti-orthogonal and equivariant."""
    return (seed_S @ random_orth_equivariant(rng)) % 3


def random_orth_equivariant(rng, steps=6):
    """Random equivariant orthogonal matrix: product of equivariant 'reflections' R = I - 2 P
    where P is an equivariant orthogonal projection onto an invariant non-degenerate subspace.
    Simplest: P = projection onto irr (x) <m> for a non-isotropic... irr(x)<m> is isotropic
    (Omega(m,m)=0), so use planes: P onto irr (x) Y, Y a random symplectic plane (dim U = 4)."""
    M = np.eye(Q, dtype=np.int64)
    for _ in range(steps):
        B = random_invariant_plane(rng)
        if B is None:
            continue
        P = projection(B)
        M = (M @ ((np.eye(Q, dtype=np.int64) - 2 * P) % 3)) % 3
    return M


def projection(B):
    gram = (B @ B.T) % 3
    n = gram.shape[0]
    r, piv = rref(np.hstack([gram, np.eye(n, dtype=np.int64)]))
    inv = r[:, n:] % 3
    return (B.T @ inv @ B) % 3


def invariant_span(v):
    """F_3-span of the Q_8-orbit of the row vector v (an invariant subspace)."""
    rows = [v]
    for G in (GI, GJ, GK):
        rows.append((v @ G) % 3)
    return rref(np.array(rows) % 3)[0]


def random_invariant_plane(rng, tries=200):
    """A 4-dimensional non-degenerate invariant subspace U = irr (x) Y: the Q_8-orbit span of a
    generic vector v = x_1 (x) m_1 + x_2 (x) m_2 is irr (x) <m_1, m_2>, of dimension 4; keep it
    when the standard form restricted to it is non-degenerate (Y a symplectic plane)."""
    for _ in range(tries):
        U = invariant_span(rng.integers(0, 3, Q))
        if U.shape[0] == 4 and rank((U @ U.T) % 3) == 4:
            return U
    return None


def seed_anti():
    """One equivariant anti-orthogonal matrix, block-diagonal: on a single 4-block W4 = irr (x) F_3^2
    the commutant is M_2(F_3) (81 elements) and contains multiplier -1 maps (1 (x) diag(1,-1) in a
    symplectic basis); find one by exhaustion of the 81 block-commutant elements."""
    rows = []
    for G in (SI, SJ):
        A = np.zeros((16, 16), dtype=np.int64)
        for a_ in range(4):
            for b in range(4):
                for c in range(4):
                    A[a_ * 4 + b, c * 4 + b] += G[a_, c]
                    A[a_ * 4 + b, a_ * 4 + c] -= G[c, b]
        rows.append(A % 3)
    cb = nullspace(np.vstack(rows) % 3)
    assert cb.shape[0] == 4, cb.shape
    for coeff in itertools.product(range(3), repeat=4):
        m = ((np.array(coeff) @ cb) % 3).reshape(4, 4)
        if np.array_equal((m @ m.T) % 3, (2 * I4) % 3):
            S = block_diag(m, m, m)
            assert is_anti(S) and all(np.array_equal((G @ S) % 3, (S @ G) % 3) for G in GENS)
            return S
    raise RuntimeError("no anti-orthogonal element in the block commutant")


# ---------------------------------------------------------------- tables in projection form
A_IN = np.array(list(itertools.product((1, 2), repeat=Q)), dtype=np.int64)
PHI_A = (1 - 2 * (np.sum(A_IN == 2, axis=1) % 2)).astype(np.int64)


def phi_rows(w):
    return np.where(np.all(w != 0, axis=1), 1 - 2 * (np.sum(w == 2, axis=1) % 2), 0)


def table(B, S):
    """Cell (U = row(B), S): T[c, c'] = sum_{a: a Pi = c' B} phi(a) phi((c B + a (I - Pi)) S)."""
    r = B.shape[0]; Pi = projection(B); n = 3 ** r
    gram = (B @ B.T) % 3
    rr, piv = rref(np.hstack([gram, np.eye(r, dtype=np.int64)])); ginv = rr[:, r:] % 3
    coords = lambda X: ((X @ B.T) % 3 @ ginv) % 3          # U-coordinates of vectors in U
    aPi = (A_IN @ Pi) % 3
    cprime = coords(aPi)
    pw = 3 ** np.arange(r)
    vidx = cprime @ pw
    rest = (A_IN - aPi) % 3
    Uout = (rest @ S) % 3
    states = np.array(list(itertools.product(range(3), repeat=r)), dtype=np.int64)[:, ::-1]
    T = np.zeros((n, n), dtype=np.int64)
    for pi, c in enumerate(states):
        out = (Uout + ((c @ B) % 3 @ S) % 3) % 3
        w = phi_rows(out) * PHI_A
        T[pi] = np.bincount(vidx, weights=w, minlength=n)
    neg = ((-states) % 3) @ pw
    return T, neg


def ring_generator(B, cells, cyclic):
    q = Q; L = len(cells); d = q * L
    Pi = projection(B); I = np.eye(q, dtype=np.int64)
    A = np.zeros((d, d), dtype=np.int64)
    for j in range(L):
        A[q*j:q*j+q, q*j:q*j+q] = ((I - Pi) @ cells[j]) % 3
        k = (j + 1) % L
        blk = (Pi @ cells[k]) % 3
        A[q*j:q*j+q, q*k:q*k+q] = blk if (k or cyclic) else (-blk) % 3
    return np.hstack([np.eye(d, dtype=np.int64), A]) % 3


def index_from_generator(G, q, L):
    d = q * L; A = G[:, d:] % 3
    E = [A[q*((j-1) % L):q*((j-1) % L)+q, q*j:q*j+q] for j in range(L)]
    keys = []
    for j in range(L):
        ks = sorted({tuple(x) for x in (A_IN @ E[j]) % 3}); keys.append((ks, {k: i for i, k in enumerate(ks)}))
    P = None
    for j in range(L):
        kin, _ = keys[j]; kout, posout = keys[(j+1) % L]
        out = np.array([posout[tuple(x)] for x in (A_IN @ E[(j+1) % L]) % 3])
        base = (A_IN @ A[q*j:q*j+q, q*j:q*j+q]) % 3
        T = np.zeros((len(kin), len(kout)), dtype=object)
        for i, s in enumerate(kin):
            w = phi_rows((base + np.array(s)) % 3) * PHI_A
            for kk in range(len(kout)):
                T[i, kk] = int(w[out == kk].sum())
        P = T if P is None else P.dot(T)
    return int(np.trace(P))


def v2(n):
    n = abs(int(n)); k = 0
    if n == 0:
        return 99
    while n % 2 == 0:
        n //= 2; k += 1
    return k


def trJ(M, neg):
    return int(sum(int(M[neg[i], i]) for i in range(M.shape[0])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--survey", action="store_true")
    ap.add_argument("--minutes", type=float, default=10); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    S0 = seed_anti()
    print("equivariant anti-orthogonal seed found; commutant dim", COMM.shape[0], flush=True)
    if a.selftest:
        B = random_invariant_plane(rng)
        cells = [random_anti_equivariant(rng, S0) for _ in range(4)]
        for S in cells:
            assert is_anti(S) and all(np.array_equal((G @ S) % 3, (S @ G) % 3) for G in GENS)
        assert all(np.array_equal(rref(np.vstack([B, (B @ G) % 3]))[0], B) for G in (GI, GJ, GK))
        for cyclic in (False, True):
            G = ring_generator(B, cells, cyclic)
            sd = is_self_dual(G)
            d = Q * 4
            H = block_diag(*([GI] * 8)); Hj = block_diag(*([GJ] * 8))     # internal Q_8 on all 8 blocks
            inv = preserves(G, H) and preserves(G, Hj)
            idx_dp = index_from_generator(G, Q, 4)
            Ts = [table(B, S) for S in cells]; neg = Ts[0][1]
            M = Ts[0][0].astype(object)
            for T, _ in Ts[1:]:
                M = M.dot(T.astype(object))
            idx_tab = int(np.trace(M)) if cyclic else trJ(M, neg)
            print(f"cyclic={cyclic}: self_dual={sd} Q8_invariant={inv} index DP={idx_dp} tables={idx_tab}")
        return
    if a.survey:
        t_end = time.time() + 60 * a.minutes
        stat = {"identical_nega": {}, "identical_cyc": {}, "prod4_nega": {}, "prod4_cyc": {}}
        pools = 0; hits = []; nprod = 0
        while time.time() < t_end:
            B = random_invariant_plane(rng)
            if B is None:
                continue
            pool = []; pool_S = []
            for _ in range(24):
                S = random_anti_equivariant(rng, S0)
                T, neg = table(B, S)
                pool.append(T); pool_S.append(S)
            pools += 1
            for T in pool[:6]:
                I, C = traces(T, neg, 8, 96)
                for L in range(1, 9):
                    for key, v in (("identical_nega", I[L]), ("identical_cyc", C[L])):
                        st = stat[key].setdefault(str(L), {"n": 0, "min_v2": 99, "gcd": 0, "min_abs": None})
                        st["n"] += 1; st["min_v2"] = min(st["min_v2"], v2(v)); st["gcd"] = gcd(st["gcd"], abs(v))
                        if v and (st["min_abs"] is None or abs(v) < st["min_abs"]): st["min_abs"] = abs(v)
            P = len(pool)
            for _ in range(4000):
                i1, i2, i3, i4 = rng.integers(0, P, 4)
                M = pool[i1].astype(object).dot(pool[i2]).dot(pool[i3]).dot(pool[i4])
                nprod += 1
                for key, v in (("prod4_nega", trJ(M, neg)), ("prod4_cyc", int(np.trace(M)))):
                    st = stat[key].setdefault("4", {"n": 0, "min_v2": 99, "gcd": 0, "min_abs": None, "v2_hist": {}})
                    st["n"] += 1; st["min_v2"] = min(st["min_v2"], v2(v)); st["gcd"] = gcd(st["gcd"], abs(v))
                    st["v2_hist"][str(v2(v))] = st["v2_hist"].get(str(v2(v)), 0) + 1
                    if v and (st["min_abs"] is None or abs(v) < st["min_abs"]): st["min_abs"] = abs(v)
                    if abs(v) == 24:
                        hits.append({"kind": key, "B": B.tolist(), "cells": [pool_S[i].tolist() for i in (i1, i2, i3, i4)], "index": v})
        out = {"pools": pools, "products": nprod, "stats": stat, "hits": hits}
        print(json.dumps({k: v for k, v in out.items() if k != "hits"}, indent=1)); print("hits:", len(hits))
        if a.out:
            Path(a.out).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
