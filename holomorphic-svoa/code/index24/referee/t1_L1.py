"""L1 (Poisson) and L3, brute force; base case; fast-vs-brute index cross-check."""
import numpy as np
from fractions import Fraction
from f3 import *

rng = np.random.default_rng(1)

# chi-hat check, exact in Q(w) via complex numbers
w = np.exp(2j * np.pi / 3)
chi = {0: 0, 1: 1, 2: -1}
for t in range(3):
    ch = sum(chi[s] * w ** (s * t) for s in range(3))
    assert abs(ch - 1j * np.sqrt(3) * chi[t]) < 1e-12
print("chi-hat(t) = i sqrt3 chi(t): ok")

# tetracode words
T = tetracode()
W = all_words(T)
print("tetracode weights:", sorted(set((W != 0).sum(1).tolist())), "self-dual:", is_self_dual(T, 4))

# fast vs brute
for trial in range(300):
    n = int(rng.integers(2, 13))
    k = int(rng.integers(1, min(n, 9) + 1))
    G = rng.integers(0, 3, (k, n))
    assert index_fast(G) == index_brute(G), (G,)
print("index_fast == index_brute on 300 random codes: ok")

# L1: Index(W^perp) vs 3^{n/2-k} Index(W) (exact rational), brute force both sides.
stats = {}
for trial in range(1200):
    n = int(rng.integers(3, 17))
    kmax = min(n - 1, 11)
    kmin = max(1, n - 11)
    if kmin > kmax:
        continue
    k = int(rng.integers(kmin, kmax + 1))
    G = rng.integers(0, 3, (k, n))
    # occasionally make codes with many full-support words: add all-ones-ish rows
    if trial % 3 == 0:
        G[0] = rng.choice([1, 2], n)
    R, piv = rref(G)
    k = R.shape[0]
    H = dual(R, n)
    if 3 ** k > 3 ** 11 or 3 ** (n - k) > 3 ** 11:
        continue
    iW = index_brute(R)
    iWp = index_brute(H)
    # predicted with i^n factor
    lhs = Fraction(iWp)
    rhs_abs = Fraction(3) ** (Fraction(n, 2) - k) if n % 2 == 0 else None
    key = n % 4
    if n % 2 == 1:
        ok = (iW == 0 and iWp == 0)
        res = "both zero" if ok else f"FAIL {iW} {iWp}"
    else:
        sign = 1 if n % 4 == 0 else -1
        pred = sign * Fraction(3) ** (n // 2 - k) * iW
        ok = (lhs == pred)
        # also check whether the 4|n formula (sign +1) would hold
        pred_plus = Fraction(3) ** (n // 2 - k) * iW
        res = (ok, lhs == pred_plus)
    s = stats.setdefault(key, {"n": 0, "ok": 0, "plus_ok": 0, "nonzero": 0})
    s["n"] += 1
    if n % 2 == 1:
        s["ok"] += ok
    else:
        s["ok"] += ok
        s["plus_ok"] += (lhs == pred_plus)
    s["nonzero"] += (iW != 0)
    # L3
    if n % 2 == 0:
        assert iW % 2 == 0 and iWp % 2 == 0
for key in sorted(stats):
    print(f"n mod 4 = {key}: {stats[key]}  (ok = matches i^n-signed formula; plus_ok = matches formula with i^n=+1)")
print("L3 (Index even for n even): ok on all tested codes")
