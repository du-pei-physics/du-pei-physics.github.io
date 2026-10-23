"""Transfer cells, general form, with exact traces.  Self-contained (numpy only).

A cell is a pair (B, S) over F_3 with q = 12:
  * S  (q x q) with S S^T = -I;
  * B  (r x q) a basis of a non-degenerate subspace U of F_3^q (Gram G = B B^T invertible).
With Pi = B^T G^{-1} B (orthogonal projection onto U, acting on row vectors) put
    A(z) = (I + (z - 1) Pi) S   over F_3[z]/(z^L - eps),  eps = -1 (negacyclic) or +1 (cyclic).
A(z) A(z^{-1})^T = -I, so [I | A(z)] is self-dual of length 24 L.
The legacy form (Q, S0, r) of ../transfer_cells.py is B = Q[:, :r]^T, S = Q S0.
New compared with the legacy form: U need not have an orthonormal basis (non-square
discriminant is allowed).

States are coordinates c in F_3^r of u = c B in U.  Local table
    T[c, c'] = sum_{a in {1,2}^q, a Pi = c' B} phi(a) phi((c B + a (I - Pi)) S),
    Index(C_L^{eps=-1}) = Tr(J T^L),   Index(C_L^{eps=+1}) = Tr(T^L).
"""
from __future__ import annotations

import itertools
from functools import lru_cache

import numpy as np

Q = 12


@lru_cache(None)
def _inputs(q):
    a = np.array(list(itertools.product((1, 2), repeat=q)), dtype=np.int64)
    return a, 1 - 2 * (np.sum(a == 2, axis=1) % 2)


A_IN, PHI_A = _inputs(Q)


# ----------------------------------------------------------------------------- F_3 basics
def inv3(M):
    M = np.asarray(M, dtype=np.int64) % 3
    n = M.shape[0]
    a = np.hstack([M, np.eye(n, dtype=np.int64)])
    for c in range(n):
        p = c + int(np.nonzero(a[c:, c])[0][0])
        a[[c, p]] = a[[p, c]]
        a[c] = (a[c] * a[c, c]) % 3                 # x^{-1} = x in F_3^*
        for r_ in range(n):
            if r_ != c and a[r_, c]:
                a[r_] = (a[r_] - a[r_, c] * a[c]) % 3
    return a[:, n:]


def rank3(M):
    a = np.asarray(M, dtype=np.int64).copy() % 3
    r = 0
    for c in range(a.shape[1]):
        nz = np.nonzero(a[r:, c])[0]
        if not len(nz):
            continue
        p = r + nz[0]
        a[[r, p]] = a[[p, r]]
        a[r] = (a[r] * a[r, c]) % 3
        for i in range(a.shape[0]):
            if i != r and a[i, c]:
                a[i] = (a[i] - a[i, c] * a[r]) % 3
        r += 1
        if r == a.shape[0]:
            break
    return r


def reflection(v):
    v = np.asarray(v, dtype=np.int64) % 3
    nn = int(v @ v) % 3
    assert nn, "isotropic"
    return (np.eye(len(v), dtype=np.int64) - 2 * nn * np.outer(v, v)) % 3


def A0(q=Q):
    a = np.zeros((q, q), dtype=np.int64)
    for i in range(0, q, 2):
        a[i:i + 2, i:i + 2] = [[1, 1], [1, 2]]
    return a


def phi_rows(X):
    return np.where(np.all(X != 0, axis=-1), 1 - 2 * (np.sum(X == 2, axis=-1) % 2), 0)


@lru_cache(None)
def states(r):
    """All c in F_3^r, row i has index sum_j c_j 3^j."""
    if r == 0:
        return np.zeros((1, 0), dtype=np.int64)
    return np.array(list(itertools.product(range(3), repeat=r)), dtype=np.int64)[:, ::-1].copy()


def check_cell(B, S):
    I = np.eye(S.shape[0], dtype=np.int64)
    G = (B @ B.T) % 3
    return bool(np.array_equal((S @ S.T) % 3, (2 * I) % 3)) and rank3(G) == B.shape[0]


# ----------------------------------------------------------------------------- transfer matrix
def transfer(B, S):
    """Returns T (3^r x 3^r int64), neg (index of -c), norm (c G c^T mod 3)."""
    B = np.asarray(B, dtype=np.int64) % 3
    q = S.shape[0]
    a_in, phi_a = _inputs(q)
    r = B.shape[0]
    G = (B @ B.T) % 3
    Ginv = inv3(G)
    C = states(r)
    pw = 3 ** np.arange(r)
    cp = (a_in @ B.T % 3) @ Ginv % 3                     # coordinates of a Pi
    Y = ((a_in - cp @ B) % 3) @ S % 3                    # a (I - Pi) S
    CB = (C @ B % 3) @ S % 3                             # (c B) S
    out = cp @ pw
    n = 3 ** r
    T = np.zeros((n, n), dtype=np.int64)
    chunk = max(1, 2 ** 21 // (len(a_in) * q))
    for s in range(0, n, chunk):
        X = (Y[None, :, :] + CB[s:s + chunk, None, :]) % 3
        w = phi_rows(X) * phi_a[None, :]
        for k in range(w.shape[0]):
            T[s + k] = np.bincount(out, weights=w[k], minlength=n)
    neg = ((-C) % 3) @ pw
    norm = np.einsum('ij,jk,ik->i', C, G, C) % 3
    return T, neg, norm


def blocks(T, neg, norm, split_norm=True):
    """Integer blocks of T on (norm class) x (J = +-1), in the basis e_c +- e_{-c}.
    (If the cell size q is not divisible by 3, T shifts norm classes: use split_norm=False.)
    Returns list of (nclass, jsign, M)."""
    n = T.shape[0]
    out = []
    if not split_norm:
        norm = np.zeros(n, dtype=np.int64)
    for nc in range(3):
        idx = [i for i in range(n) if norm[i] == nc]
        reps = [i for i in idx if i <= neg[i]]               # 0 is its own negative
        plus = reps
        minus = [i for i in reps if neg[i] != i]
        if plus:
            Mp = T[np.ix_(plus, plus)] + T[np.ix_([neg[i] for i in plus], plus)]
            for a, i in enumerate(plus):
                if neg[i] == i:
                    Mp[a] = T[i, plus]
            out.append((nc, 1, Mp))
        if minus:
            Mm = T[np.ix_(minus, minus)] - T[np.ix_([neg[i] for i in minus], minus)]
            out.append((nc, -1, Mm))
    return out


# ----------------------------------------------------------------------------- exact traces
def _primes(k, top=2 ** 26):
    ps, x = [], top
    while len(ps) < k:
        x -= 1
        if x % 2 and all(x % d for d in range(3, int(x ** 0.5) + 1, 2)):
            ps.append(x)
    return ps


PRIMES = _primes(24)


def power_traces(mats, Lmax, signs=None):
    """Exact sum_k sign_k Tr(M_k^L) for L = 1..Lmax (Python ints), by CRT over ~2^26 primes."""
    if signs is None:
        signs = [1] * len(mats)
    nmax = max(m.shape[0] for m in mats)
    rho = max(float(np.abs(m).sum(axis=1).max()) for m in mats) or 1.0
    bits = np.log2(sum(m.shape[0] for m in mats)) + Lmax * np.log2(rho) + 2
    k = int(bits // 25) + 2
    assert k <= len(PRIMES) and nmax * (2 ** 26) ** 2 < 2 ** 63
    ps = PRIMES[:k]
    res = np.zeros((Lmax, k), dtype=object)
    for j, p in enumerate(ps):
        acc = [0] * Lmax
        for m, s in zip(mats, signs):
            mp = m % p
            P = mp.copy()
            for L in range(1, Lmax + 1):
                if L > 1:
                    P = (P @ mp) % p
                acc[L - 1] += s * int(np.trace(P) % p)
        for L in range(Lmax):
            res[L, j] = acc[L] % p
    M = 1
    for p in ps:
        M *= p
    out = []
    for L in range(Lmax):
        x = 0
        for j, p in enumerate(ps):
            m = M // p
            x = (x + int(res[L, j]) * m * pow(m, -1, p)) % M
        out.append(x - M if x > M // 2 else x)
    return out


def index_traces(T, neg, norm, Lmax, split_norm=True):
    """(I_L negacyclic = Tr J T^L, Tr T^L cyclic) for L = 1..Lmax, exact."""
    bl = blocks(T, neg, norm, split_norm)
    mats = [m for _, _, m in bl]
    sg = [j for _, j, _ in bl]
    neg_tr = power_traces(mats, Lmax, sg)
    cyc_tr = power_traces(mats, Lmax, [1] * len(mats))
    return neg_tr, cyc_tr


# ----------------------------------------------------------------------------- random cells
def rand_vec(rng, q=Q, norm=None, weight=None):
    assert weight is None or norm is not None or weight % 3, "weight divisible by 3 is always isotropic"
    while True:
        if weight is None:
            v = rng.integers(0, 3, q)
        else:
            v = np.zeros(q, dtype=np.int64)
            v[rng.choice(q, weight, replace=False)] = rng.integers(1, 3, weight)
        nn = int(v @ v) % 3
        if (norm is None and nn) or (norm is not None and nn == norm):
            return v


def random_orth(rng, depth, base=None, q=Q):
    M = np.eye(q, dtype=np.int64) if base is None else base.copy()
    for _ in range(depth):
        M = (reflection(rand_vec(rng, q)) @ M) % 3
    return M


def random_subspace(rng, r, q=Q):
    """Random non-degenerate r-dim subspace (any discriminant)."""
    while True:
        B = rng.integers(0, 3, (r, q))
        if rank3((B @ B.T) % 3) == r:
            return B


def legacy_to_general(Qm, S0, r):
    Qm = np.asarray(Qm, dtype=np.int64) % 3
    return Qm[:, :r].T.copy(), (Qm @ np.asarray(S0)) % 3


def generator(B, S, L, eps=-1):
    """[I | A] with A the qL x qL expansion of A(z) = (I + (z-1)Pi) S over z^L = eps."""
    B = np.asarray(B, dtype=np.int64) % 3
    Q = S.shape[0]
    G = (B @ B.T) % 3
    Pi = (B.T @ inv3(G) % 3) @ B % 3
    I = np.eye(Q, dtype=np.int64)
    A_same = ((I - Pi) % 3) @ S % 3
    A_next = Pi @ S % 3
    A = np.zeros((Q * L, Q * L), dtype=np.int64)
    for j in range(L):
        A[Q * j:Q * j + Q, Q * j:Q * j + Q] = (A[Q * j:Q * j + Q, Q * j:Q * j + Q] + A_same) % 3
        k = (j + 1) % L
        blk = A_next if k != 0 else (eps * A_next) % 3
        A[Q * j:Q * j + Q, Q * k:Q * k + Q] = (A[Q * j:Q * j + Q, Q * k:Q * k + Q] + blk) % 3
    return np.hstack([np.eye(Q * L, dtype=np.int64), A]) % 3


def generator_ring(B, Slist, eps=-1):
    """Ring with common U = rowspace(B) and per-cell matrices S_0..S_{L-1}:
    check block j = (a_j (I - Pi) + a_{j-1} Pi) S_j, with sign eps on the wrap L-1 -> 0."""
    B = np.asarray(B, dtype=np.int64) % 3
    L = len(Slist)
    Q = Slist[0].shape[0]
    G = (B @ B.T) % 3
    Pi = (B.T @ inv3(G) % 3) @ B % 3
    I = np.eye(Q, dtype=np.int64)
    A = np.zeros((Q * L, Q * L), dtype=np.int64)
    for j, Sj in enumerate(Slist):
        A[Q * j:Q * j + Q, Q * j:Q * j + Q] = (A[Q * j:Q * j + Q, Q * j:Q * j + Q] + ((I - Pi) % 3) @ Sj) % 3
        i = (j - 1) % L
        blk = (Pi @ Sj) % 3 if j != 0 else (eps * (Pi @ Sj)) % 3
        A[Q * i:Q * i + Q, Q * j:Q * j + Q] = (A[Q * i:Q * i + Q, Q * j:Q * j + Q] + blk) % 3
    return np.hstack([np.eye(Q * L, dtype=np.int64), A]) % 3
