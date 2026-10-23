#!/usr/bin/env python3
# Copied from claude/census_general/k6_general.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Complete NS/Ramond characters of T = V_{Lambda(C)}//Z_6 for any k = 6 code given by its generator, 22 Oct 2026.

Generalizes the census scripts of the holomorphic-svoa note (code/census/fc_seed_census.py, fc_seed_level1.py,
fc_seed_level2.py), which build the code from a Paley cell (Q, R) of an anti-periodic chain, to an arbitrary generator
(I_36 | A) of a chain of three 12-cells, anti-periodic ("nega", g = shift) or periodic ("cyc", g = -shift).  In both
cases g is a product of 24 signed 3-cycles with sign product -1 and g^3 = -1, which is all the formulas of
field_content.md use (traces tr g^j = 72, 0, 0, -72, 0, 0 on the oscillators; only g^2, g^4 fix momenta; twisted ground
energies 8/3, 19/6, 9/2), together with three per-code hypotheses that are checked here:
  * balance: the g^2-fixed full-support words are half even, half odd;
  * Cbar (the g^2-fixed subcode contracted to one coordinate per 3-cycle) is self-dual of dimension 12;
  * Cbar has no words of weight 1, 2, 4, 5 (so L = Lambda^{g^2} has no vectors of norm 1, 2, 4, 5).
Formulas (copied from the census scripts):
  n = (N + 8F)/6;  dim T_1 = A_6/6;  dim T_{3/2} = (N_3 + 2 A_3(Cbar))/6;  dim T_2 = (5256 + N_4 + 72 A_6)/6;
  dim T_{5/2} = (N_5 + 72 N_3)/6, with N_3 = 144 + A_9 + 6 A_6, N_4 = A_12 + 9 A_9 + 15 A_6,
  N_5 = A_15 + 12 A_12 + 36 A_9 + 152 A_6;  the degree-6 polynomial a in t = (theta_3/eta)^12 follows.
  Per-code check: D_1 = (144 N + 8 A_69)/6 + X(Cbar) + 864,288,768 (twisted sectors) must equal -4096 a_1.
A_w (w <= 15) and A_69, N come from weight- and zero-graded transfer traces over the generator (exact integers).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import itertools, sys
from pathlib import Path
import numpy as np
CODE = Path(__file__).resolve().parents[1]
for p in (CODE, CODE / "census"):
    sys.path.insert(0, str(p))
from f3lib import monomial_matrix, fixed_subcode, rref, is_self_dual, all_codewords   # noqa: E402
from mixed_rings import shift_monomial                                             # noqa: E402
from fc_seed_level1 import theta_twisted_level1, K_G                               # noqa: E402
from fc_seed_level2 import polynomial                                              # noqa: E402
import fc_closed_form as CF                                                        # noqa: E402
import verify as V                                                                 # noqa: E402

W = 12
ALL = np.array(list(itertools.product(range(3), repeat=W)), dtype=np.int64)
NZ = np.count_nonzero(ALL, axis=1)
CONTRACT = list(range(12)) + list(range(36, 48))          # one coordinate per 3-cycle of g^2 (first period)


def symmetry(kind):
    perm, signs = shift_monomial([W], 3, kind == "cyc", 1)
    if kind == "cyc":
        signs = (3 - signs) % 3                            # g = -shift: flip every sign (1 <-> 2)
    return monomial_matrix(perm.tolist(), signs.tolist())


def fmul(X, Y):
    """(X @ Y) mod 3 through float64 BLAS: exact, the entries are at most 12 * 2 * 2."""
    return np.rint(X.astype(np.float64) @ Y.astype(np.float64)).astype(np.int64) % 3


def graded(G, stat, cap):
    """out[k] = number of codewords with statistic k (weight or #zeros), k <= cap, by a transfer trace over the
    three blocks of (I_36 | A)."""
    A = G[:, 36:]
    blk = lambda i, j: A[12 * i:12 * i + 12, 12 * j:12 * j + 12]
    V_ = ALL[NZ >= 12 - cap] if stat == "zeros" else ALL
    st_a = (12 - np.count_nonzero(V_, axis=1)) if stat == "zeros" else np.count_nonzero(V_, axis=1)
    outs, states = [], []
    for i in range(3):
        o = fmul(V_, blk(i, (i + 1) % 3))
        code = o @ (3 ** np.arange(W))
        u, inv = np.unique(code, return_inverse=True)
        reps = np.zeros((len(u), W), dtype=np.int64); reps[inv] = o
        outs.append(inv); states.append(reps)
    Ms = []
    for i in range(3):
        inc = states[(i - 1) % 3]; base = fmul(V_, blk(i, i)); nout = len(states[i])
        # per parity coordinate j and value v of the incoming state: is (base_j + v) nonzero
        nzc = [[((base[:, j] + v) % 3 != 0).astype(np.int16) for v in range(3)] for j in range(W)]
        M = np.zeros((cap + 1, len(inc), nout), dtype=np.int64)
        for si, s in enumerate(inc):
            wp = nzc[0][s[0]].copy()
            for j in range(1, W):
                wp += nzc[j][s[j]]
            k = st_a + ((12 - wp) if stat == "zeros" else wp); ok = k <= cap
            M[:, si, :] = np.bincount(k[ok] * nout + outs[i][ok], minlength=(cap + 1) * nout).reshape(cap + 1, nout)
        Ms.append(M)
    # int64 is exact: entries of each M[k] are at most 3^12, so a triple product is below 3^36 * 27^2 < 2^63
    out = [0] * (cap + 1)
    nz = [[k for k in range(cap + 1) if M[k].any()] for M in Ms]
    for k0 in nz[0]:
        for k1 in nz[1]:
            if k0 + k1 > cap:
                continue
            P = Ms[0][k0] @ Ms[1][k1]
            for k2 in nz[2]:
                if k0 + k1 + k2 <= cap:
                    out[k0 + k1 + k2] += int(np.einsum("ij,ji->", P, Ms[2][k2]))
    return out


def analyse(G, kind):
    G = np.asarray(G, dtype=np.int64) % 3
    assert G.shape == (36, 72) and np.array_equal(G[:, :36], np.eye(36, dtype=np.int64)) and is_self_dual(G)
    g = symmetry(kind)
    assert not ((G @ g @ G.T) % 3).any(), "g does not preserve C"
    g3 = np.linalg.matrix_power(g, 3) % 3
    assert np.array_equal(g3, (2 * np.eye(72, dtype=np.int64)) % 3), "g^3 != -1"
    h = (g @ g) % 3
    I, N = V.band_index(G, [12, 12, 12])
    Fx = fixed_subcode(G, h); Fc = rref(Fx[:, CONTRACT] % 3)[0]
    cb = np.zeros(25, dtype=np.int64); e = o = 0
    for wds in all_codewords(Fc):
        np.add.at(cb, np.count_nonzero(wds, axis=1), 1)
        full = np.all(wds != 0, axis=1); par = np.sum(wds == 2, axis=1) % 2
        e += int((full & (par == 0)).sum()); o += int((full & (par == 1)).sum())
    cb = [int(x) for x in cb]; F = e + o
    hyp = {"fixed_dim_12": int(Fx.shape[0]) == 12, "Cbar_self_dual": bool(Fc.shape[0] == 12 and is_self_dual(Fc)),
           "balanced": e == o, "Cbar_no_wt_1245": all(cb[w] == 0 for w in (1, 2, 4, 5)), "abs_index_24": abs(I) == 24}
    rec = {"kind": kind, "I": int(I), "N": int(N), "F": F, "F_even_odd": [e, o], "Cbar_enum": cb, "hypotheses": hyp}
    if not all(hyp.values()):
        rec["ok"] = False
        return rec
    wt = graded(G, "weight", 15); z = graded(G, "zeros", 3)
    assert wt[0] == 1 and all(wt[w] == 0 for w in range(1, 16) if w % 3) and z[1] == z[2] == 0
    assert z[0] == N, (z[0], N)
    A6, A9, A12, A15, A69 = wt[6], wt[9], wt[12], wt[15], z[3]
    n6 = N + 8 * F; assert n6 % 6 == 0; n = n6 // 6
    N3 = 144 + A9 + 6 * A6; N4 = A12 + 9 * A9 + 15 * A6; N5 = A15 + 12 * A12 + 36 * A9 + 152 * A6
    vals = [A6, N3 + 2 * cb[3], 5256 + N4 + 72 * A6, N5 + 72 * N3]
    assert all(x % 6 == 0 for x in vals), vals
    T1, T32, T2, T52 = (x // 6 for x in vals)
    a = polynomial(n, T1, T32, T2, T52)
    NS = CF.eval_poly(a); RR = CF.eval_poly(a, "R")
    X, _ = theta_twisted_level1(cb)
    D1 = (144 * N + 8 * A69) // 6 + X + K_G
    derived = [NS.get(0, 0), NS.get(1, 0), NS.get(2, 0), RR.get(4, 0), RR.get(6, 0), RR.get(8, 0)]
    rec.update({"A": {"6": A6, "9": A9, "12": A12, "15": A15, "69": A69}, "n": n, "T1": T1, "T32": T32, "T2": T2,
                "T52": T52, "a": a, "D1": D1, "D1_check": D1 == -4096 * a[1], "derived_nonneg": all(x >= 0 for x in derived)})
    rec["ok"] = bool(rec["D1_check"] and rec["derived_nonneg"])
    return rec
