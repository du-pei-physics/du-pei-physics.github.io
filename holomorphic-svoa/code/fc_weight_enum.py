# Copied from claude/paper/fc_weight_enum.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Exact Hamming weight enumerator of the canonical seed code C (length 72) by a transfer-matrix DP.

C = {(a | aA)}, a = (a0,a1,a2) in F_3^36, A = [[D,E,0],[0,D,E],[-E,0,D]] (canonical_seed.json), E = Q_R S_R of rank 3.
Check block j is c_j = a_j D + a_{j-1} E (j=1,2), c_0 = a_0 D - a_2 E.  The state passed from period j to j+1 is
s_j = a_j Q_R in F_3^3 (27 states).  M_y[t, s](x) = sum_{a: a Q_R = s} x^{wt(a)+wt(aD + t S_R)} y^{#2's},
W_y(x) = sum_s (M_y^3)[-s, s].  y = +1 gives the weight enumerator, y = -1 the signed one (sign = (-1)^{#2's});
its x^72 coefficient is the index.  All arithmetic in int64; every coefficient is bounded by 3^36 < 2^63.

Also: the contracted g^2-fixed code Cbar (self-dual [24,12]) and its weight enumerator, by enumeration of 3^12 words.
Writes paper/fc_weight_enum.json.
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from f3lib import is_self_dual, rank, monomial_matrix, fixed_subcode, rref, all_codewords

seed = json.load(open(HERE.parent / "certificates" / "original" / "canonical_seed.json"))
Q = np.array(seed["Q"]) % 3; S = np.array(seed["S"]) % 3; R = seed["R"]; D = np.array(seed["D"]) % 3; E = np.array(seed["E"]) % 3
G = np.array(seed["generator"]) % 3
QR = Q[:, R]; SR = S[R, :]
assert not (((QR @ SR) % 3 - E) % 3).any()
A = np.block([[D, E, 0 * D], [0 * D, D, E], [(-E) % 3, 0 * D, D]]) % 3
assert np.array_equal(G, np.hstack([np.eye(36, dtype=np.int64), A]) % 3)

allA = np.array(list(itertools.product(range(3), repeat=12)), dtype=np.int64)          # 3^12 x 12
aD = (allA @ D) % 3; sidx = ((allA @ QR) % 3) @ np.array([9, 3, 1])
wa = np.count_nonzero(allA, axis=1); ta = np.count_nonzero(allA == 2, axis=1)
states = np.array(list(itertools.product(range(3), repeat=3)), dtype=np.int64)
neg = [int(((-s) % 3) @ np.array([9, 3, 1])) for s in states]


def transfer(y):
    M = np.zeros((27, 27, 73), dtype=np.int64)
    for ti, t in enumerate(states):
        v = (aD + t @ SR) % 3
        w = wa + np.count_nonzero(v, axis=1); par = (ta + np.count_nonzero(v == 2, axis=1)) % 2
        wt = np.ones_like(w) if y == 1 else 1 - 2 * par
        np.add.at(M[ti], (sidx, w), wt)
    return M


def polymatmul(X, Y):
    Z = np.zeros((27, 27, 73), dtype=np.int64)
    for i in range(73):
        Xi = X[:, :, i]
        if not Xi.any():
            continue
        for j in range(73 - i):
            if Y[:, :, j].any():
                Z[:, :, i + j] += Xi @ Y[:, :, j]
    return Z


out = {}
for y, name in ((1, "plain"), (-1, "signed")):
    M = transfer(y); M3 = polymatmul(polymatmul(M, M), M)
    W = sum(M3[neg[s], s] for s in range(27))
    out[name] = [int(c) for c in W]
plain, signed = out["plain"], out["signed"]
assert sum(plain) == 3 ** 36 and plain[0] == 1 and signed[0] == 1
assert plain[72] == 83000 and signed[72] == -24
assert all(plain[w] == 0 for w in range(73) if w % 3)
# MacWilliams self-check (C self-dual): W(x) = 3^-36 sum_w A_w (1+2x)^(72-w) (1-x)^w, exact with sympy
import sympy as sp
x = sp.symbols("x")
Wx = sum(c * x ** w for w, c in enumerate(plain))
mw = sp.expand(sum(c * (1 + 2 * x) ** (72 - w) * (1 - x) ** w for w, c in enumerate(plain)) / sp.Integer(3) ** 36)
assert sp.expand(mw - Wx) == 0
out["macwilliams_ok"] = True
out["even_odd_maximal"] = [(plain[72] + signed[72]) // 2, (plain[72] - signed[72]) // 2]
out["min_distance"] = min(w for w in range(1, 73) if plain[w])

# --- the symmetry g and the contracted g^2-fixed code
pass
from mixed_rings import shift_monomial
perm, signs = shift_monomial([12], 3, False, 1)
g = monomial_matrix(perm.tolist(), signs.tolist()); h = (g @ g) % 3
F = fixed_subcode(G, h); cols = list(range(12)) + list(range(36, 48))
Fc = rref(F[:, cols] % 3)[0]
assert Fc.shape[0] == 12 and is_self_dual(Fc)
wbar = np.zeros(25, dtype=np.int64); sbar = np.zeros(25, dtype=np.int64)
for words in all_codewords(Fc):
    w = np.count_nonzero(words, axis=1); p = np.count_nonzero(words == 2, axis=1) % 2
    np.add.at(wbar, w, 1); np.add.at(sbar, w, 1 - 2 * p)
out["Cbar_generator"] = Fc.tolist(); out["Cbar_weight_enumerator"] = [int(c) for c in wbar]
out["Cbar_signed_enumerator"] = [int(c) for c in sbar]
out["g_perm"] = perm.tolist(); out["g_signs"] = signs.tolist()
(HERE / "out" / "fc_weight_enum.json").write_text(json.dumps(out))
print("A_w (w: count):", {w: c for w, c in enumerate(plain) if c})
print("signed:", {w: c for w, c in enumerate(signed) if c})
print("d =", out["min_distance"], "even/odd maximal", out["even_odd_maximal"])
print("Cbar:", {w: c for w, c in enumerate(wbar) if c}, "signed", {w: c for w, c in enumerate(sbar) if c})
