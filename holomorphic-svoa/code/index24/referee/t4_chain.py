"""L4: explicit chain between two self-dual codes via a basis e_i of C/K and the
dual basis f_i of C'/K (K = C cap C'), switching e_i -> f_i one at a time.
Checks every intermediate code is self-dual, consecutive ones meet in codim 1,
and Index alternates sign mod 8 along the chain."""
import sys
import numpy as np
from f3 import *

n = int(sys.argv[1])
trials = int(sys.argv[2])
rng = np.random.default_rng(4)
m = n // 2
lens = []
for t in range(trials):
    C = rref(random_monomial(random_selfdual_AAT(n, rng), rng))[0]
    if t % 2 == 0:
        Cp = rref(random_monomial(tetra_sum(n), rng))[0]  # tetracode-sum endpoint
    else:
        Cp = rref(random_monomial(random_selfdual_AAT(n, rng), rng))[0]
    K = dual(np.vstack([C, Cp]), n)          # (C + C')^perp = C cap C'
    K = rref(K)[0] if K.shape[0] else np.zeros((0, n), dtype=np.int64)
    r = m - K.shape[0]
    # complete K to C and to C'
    def complete(Kb, Cb):
        basis = list(Kb)
        ext = []
        for row in Cb:
            if rank(np.array(basis + [row])) > len(basis):
                basis.append(row)
                ext.append(row)
        return np.array(ext).reshape(-1, n)
    E = complete(K, C)
    F = complete(K, Cp)
    assert E.shape[0] == r == F.shape[0]
    if r == 0:
        lens.append(0)
        continue
    Pm = (E @ F.T) % 3
    Pinv = inv_mod3(Pm)                        # pairing C/K x C'/K is perfect
    F2 = (Pinv.T @ F) % 3                      # e_i . f2_j = delta_ij
    assert ((E @ F2.T) % 3 == np.eye(r, dtype=np.int64)).all()
    chain = []
    for k in range(r + 1):
        rows = [K] if K.shape[0] else []
        rows += [F2[:k], E[k:]]
        Ck = np.vstack(rows) % 3
        assert is_self_dual(Ck, n)
        chain.append(Ck)
    assert rank(np.vstack([chain[0], C])) == m and rank(np.vstack([chain[-1], Cp])) == m
    idx = [index_fast(Ck) for Ck in chain]
    for a, b in zip(chain, chain[1:]):
        inter = m - (rank(np.vstack([a, b])) - m)
        assert inter == m - 1
    for x, y in zip(idx, idx[1:]):
        assert (x + y) % 8 == 0
    lens.append(r)
print(f"n={n}: {trials} chains ok; chain lengths seen {sorted(set(lens))}")
