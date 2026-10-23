# Copied from claude/census/hpcell.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Fast exact local tables for Hadamard-Paley transfer cells with an arbitrary delayed set.

Cell: Q = (H4 + H4 + H4) P_pi D_s  (three F_3-Fourier blocks placed on the 12 Paley points,
      with column signs s), S = skew Paley conference matrix C_12 on P^1(F_11)
      (coordinate 0 = inf, coordinate k = field element k mod 11 for k = 1..11),
      R = delayed set (any subset of the 12 positions).
Configuration = (blocks, parities, R):
  * blocks: three disjoint 4-subsets of {0..11} (positions = Paley points);
  * parities: one bit per block = parity of the number of -1 column signs on that block
    (sign patterns with an even number of -1's are +-characters of Z_2^2 and are absorbed
    by signed row permutations of H4, so only the parity matters);
  * R: delayed positions.
Row monomials of Q only change T by a global sign, so (|I_L|, nmax, |spectrum|) depend only
on the configuration.  T[p, v] = sum_{a in {1,2}^12, (aQ)_R = v} phi(a) phi(x S),
x = aQ with x_R := p.   Index_L = Tr(J T^L).
"""
import itertools
from functools import lru_cache
import numpy as np

H4 = np.array([[1, 1, 1, 1], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1]], dtype=np.int64)
HB = np.zeros((12, 12), dtype=np.int64)
for _k in range(3):
    HB[4 * _k:4 * _k + 4, 4 * _k:4 * _k + 4] = H4
PW = 3 ** np.arange(12)
A_IN = np.array(list(itertools.product((1, 2), repeat=12)), dtype=np.int64)   # 4096 x 12
PHI_A = (1 - 2 * (np.sum(A_IN == 2, axis=1) % 2)).astype(np.int64)
W_ALL = np.array(list(itertools.product(range(3), repeat=12)), dtype=np.int64)[:, ::-1]  # idx = sum w_i 3^i


def legendre(a, p=11):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def paley12():
    S = np.zeros((12, 12), dtype=np.int64)
    for i in range(12):
        for j in range(12):
            S[i, j] = 0 if i == j == 0 else (1 if i == 0 else (-1 if j == 0 else legendre(j - i))) % 3
    assert np.array_equal((S @ S.T) % 3, (2 * np.eye(12, dtype=np.int64)) % 3)
    return S


PALEY = paley12()


def phi_vec(X):
    """phi of each row: 0 if a zero entry, else (-1)^{#2's}."""
    return np.where(np.all(X != 0, axis=1), 1 - 2 * (np.sum(X == 2, axis=1) % 2), 0).astype(np.int64)


@lru_cache(maxsize=4)
def phi_out_for(S_key):
    S = np.array(S_key, dtype=np.int64)
    return phi_vec((W_ALL @ S) % 3)                                      # 3^12 entries


def matrix_Q(blocks, parities):
    """Q = HB P_pi D_s with pi = concatenation of the blocks and one sign flip per odd block."""
    pi = [p for b in blocks for p in b]
    P = np.zeros((12, 12), dtype=np.int64); P[np.arange(12), pi] = 1
    Q = (HB @ P) % 3
    for b, par in zip(blocks, parities):
        if par:
            Q[:, b[0]] = (-Q[:, b[0]]) % 3
    return Q


def tables(Q, R, S=PALEY, want_abs=True):
    """Signed (and unsigned) local tables for delayed set R (sorted positions)."""
    R = sorted(R); r = len(R); n = 3 ** r
    rest = [j for j in range(12) if j not in R]
    phi_out = phi_out_for(tuple(map(tuple, S.tolist())))
    y = (A_IN @ Q) % 3
    vidx = y[:, R] @ (3 ** np.arange(r))
    tail = y[:, rest] @ PW[rest]
    states = np.array(list(itertools.product(range(3), repeat=r)), dtype=np.int64)[:, ::-1]  # p idx = sum p_k 3^k
    off = states @ PW[R]
    neg = ((-states) % 3) @ (3 ** np.arange(r))
    idx = tail[None, :] + off[:, None]                                   # n x 4096
    w = phi_out[idx] * PHI_A[None, :]
    T = np.zeros((n, n), dtype=np.int64); Ta = None
    for p in range(n):
        T[p] = np.bincount(vidx, weights=w[p], minlength=n)
    if want_abs:
        Ta = np.zeros((n, n), dtype=np.int64)
        wa = np.abs(w)
        for p in range(n):
            Ta[p] = np.bincount(vidx, weights=wa[p], minlength=n)
    return T, Ta, neg


def tr_J(M, neg):
    return int(sum(int(M[neg[i], i]) for i in range(M.shape[0])))


def tr_J_powers(T, neg, Lmax=6):
    """Exact Tr(J T^L), L = 1..Lmax, with Python integers."""
    M = T.astype(object); P = M.copy(); out = []
    for L in range(1, Lmax + 1):
        if L > 1:
            P = P.dot(M)
        out.append(tr_J(P, neg))
    return out


def evaluate(blocks, parities, R, S=PALEY, Lmax=6, want_abs=True):
    """Quick seed test: I1 == 0 and |I3| == 24.  Returns (res, T, Ta, neg); Ta only if want_abs or hit."""
    Q = matrix_Q(blocks, parities)
    T, Ta, neg = tables(Q, R, S, want_abs=want_abs)
    I1 = tr_J(T, neg)
    T2 = T @ T; T3 = T2 @ T                                               # int64 safe: |T|<=4096, n<=243
    I3 = tr_J(T3, neg)
    I2 = tr_J(T2, neg)
    res = {"I1": I1, "I2": I2, "I3": I3}
    if I1 == 0 and abs(I3) == 24:
        if Ta is None:
            T, Ta, neg = tables(Q, R, S, want_abs=True)
        ev = np.linalg.eigvals(T.astype(float))
        res.update({"hit": True, "nmax": tr_J(Ta @ Ta @ Ta, neg),
                    "nz_eigs": int(np.sum(np.abs(ev) > 1e-6)),
                    "IL": tr_J_powers(T, neg, Lmax),
                    "eigs": sorted([(round(float(e.real), 6), round(float(e.imag), 6)) for e in ev if abs(e) > 1e-6])})
    return res, T, Ta, neg


def sector_spectra(T, neg):
    """Eigenvalues of T restricted to the J=+1 and J=-1 sectors (J f(v) = f(-v))."""
    n = T.shape[0]
    J = np.zeros((n, n)); J[np.arange(n), neg] = 1
    Tf = T.astype(float)
    out = {}
    for s, name in ((1, "J+"), (-1, "J-")):
        Pj = (np.eye(n) + s * J) / 2
        ev = np.linalg.eigvals(Pj @ Tf @ Pj)
        out[name] = sorted([complex(round(e.real, 6), round(e.imag, 6)) for e in ev if abs(e) > 1e-6], key=lambda z: (z.real, z.imag))
    return out
