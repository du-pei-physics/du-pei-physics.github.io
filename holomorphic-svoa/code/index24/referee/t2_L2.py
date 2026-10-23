"""L2 on random neighbour pairs; coset structure of D^perp/D; n = 4, 8 behaviour."""
import sys
import numpy as np
from collections import Counter
from f3 import *

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
n = int(sys.argv[1])
m = n // 2
npairs = {4: 200, 8: 200, 12: 300, 16: 200, 20: 100, 24: 150, 36: 60}[n]


def coset_analysis(C, Cp, D):
    """Brute-force D^perp, split into the 9 cosets of D, report q and full-support counts."""
    Dp = dual(D, n)
    assert Dp.shape[0] == m + 1
    # D^perp = C + C'; pick e in C\D, f in C'\D
    X = all_words(Dp)
    # coordinates of each word in D^perp/D: solve via the projections onto e, f
    e = next(r for r in C if rank(np.vstack([D, r])) > D.shape[0])
    f = next(r for r in Cp if rank(np.vstack([D, r])) > D.shape[0])
    # linear functionals on D^perp vanishing on D: x -> x.f picks the e-coefficient (times e.f), etc.
    ef = int(e @ f) % 3
    assert ef != 0
    a = (X @ f) % 3  # = alpha * (e.f) for x = alpha e + beta f + d
    b = (X @ e) % 3
    q = np.sum(X * X, axis=1) % 3
    vals = chi_prod(X)
    full = np.all(X != 0, axis=1)
    out = {}
    for A in range(3):
        for B in range(3):
            sel = (a == A) & (b == B)
            qs = set(q[sel].tolist())
            out[(A, B)] = (qs, int(full[sel].sum()), int(vals[sel].sum()))
    return out, index_brute(Dp)


bad = 0
cnt = Counter()
G = tetra_sum(n)
for step in range(30):  # burn in away from the tetracode sum
    G, _ = random_neighbor(G, rng)
    G = random_monomial(G, rng)
did_cosets = 0
for t in range(npairs):
    G = rref(random_monomial(G, rng))[0]
    Gp, D = random_neighbor(G, rng)
    assert is_self_dual(G, n) and is_self_dual(Gp, n)
    # D = C cap C' independently: (C + C')^perp
    K = dual(np.vstack([G, Gp]), n)
    assert K.shape[0] == m - 1 and rank(np.vstack([K, D])) == m - 1
    iC, iCp, iD = index_fast(G), index_fast(Gp), index_fast(D)
    iDp = index_fast(dual(D, n))
    lhs, rhs = iC + iCp, 4 * iD
    ok = (lhs == rhs)
    l1 = (iDp == 3 * iD)
    incl = (iDp == iC + iCp - iD)
    cnt[(ok, l1, incl)] += 1
    cnt[("iD even", iD % 2 == 0)] += 1
    cnt[("C mod 8", iC % 8)] += 1
    cnt[("C mod 3", iC % 3)] += 1
    cnt[("|C|", abs(iC))] += 1
    cnt[("|D|", abs(iD))] += 1
    if not ok:
        bad += 1
        if bad <= 3:
            print("L2 FAIL", n, iC, iCp, iD, iDp)
    if n <= 24 and did_cosets < (3 if n == 24 else 10):
        out, iDp_b = coset_analysis(G, Gp, D)
        assert iDp_b == iDp
        did_cosets += 1
        if did_cosets == 1 or (n in (4, 8) and iD != 0):
            print(f"n={n}: cosets (alpha,beta): (q-values, #full-support, index contribution)")
            for k_, v in out.items():
                print("   ", k_, v)
            print("    Index C, C', D, D^perp =", iC, iCp, iD, iDp)
    G = Gp
print(f"n={n}: pairs={npairs}  (L2 ok, L1 for D, inclusion-exclusion) counts:",
      {k: v for k, v in cnt.items() if isinstance(k[0], bool)})
print("   ", sorted((k, v) for k, v in cnt.items() if not isinstance(k[0], bool)))
