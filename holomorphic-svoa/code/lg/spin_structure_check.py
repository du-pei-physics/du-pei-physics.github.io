#!/usr/bin/env python3
# Copied from GappedCFT/code/spin_structure_check.py (a companion project of the scft-code-search workspace) for the holomorphic-svoa
# note; only pointers to the working note were changed (see tools/package_lg.py).
"""Spin-structure-resolved comparison of the orbifold M^n / C with the double V_Lambda(C) (x) Vbar (../../notes/lg_elliptic_genus.md section 5).

M is the k = 1 N = 2 minimal model (the IR of W = Phi^3), written as the diagonal theory of V_{sqrt3 Z}. A state of one
field has left and right labels x_L, x_R in Z (momenta X = x for NS, X = x + 3/2 for R; weights X^2/6 + oscillators),
with x_L = x_R mod 3 in the untwisted sector. Data used (one field, Z3 twist a in space, b in time):
  - a-twisted sector: x_L - x_R = a mod 3;
  - Z3 phase (time twist b): NS: omega^{b x_L}; R: omega^{b (x_L + a)}, the phase fixed by the LG elliptic genus
    (../../notes/lg_elliptic_genus.md section 3: F_{a,b} = - sum_u iota_{u-a} omega^{b(u+a)} psi_u); the constant kappa below tests other R phases;
  - fermion parity (for the spin structures with (-1)^F inserted): (-1)^{x_L + x_R}; no discrete torsion.
The orbifold weight of a lattice point (x_L, x_R) in (Z^n)^2 is
    w_M = (1/|C|) sum_{c in C} [x_L - x_R = c] sum_{c' in C} prod_i phase_i(c_i, c'_i) * sign,
and the double's weight is w_V = [x_L mod 3 in C] [x_R mod 3 in C] * sign. The theta functions (all four spin
structures; the eta factors are common) agree iff the weights agree pointwise. Checked on all points of a small box
(tetracode) and on 20,000 random points per spin structure (Golay), exactly in Z[omega]; also shows that a constant extra
R-sector phase omega^{kappa b} (kappa = 1, 2) changes the R projection when the all-ones word is not in C (it is then
excluded by T-covariance, not by the elliptic genus, which need not change; see ../../notes/lg_elliptic_genus.md), and that any
per-field bilinear phase omega^{a b} drops out (it multiplies to omega^{c . c'} = 1 on a self-orthogonal code).

usage: spin_structure_check.py
"""
import itertools, random

def codewords(G):
    k, n = len(G), len(G[0])
    return [tuple(sum(a * r[j] for a, r in zip(co, G)) % 3 for j in range(n)) for co in itertools.product(range(3), repeat=k)]

def zsum_omega(exps):
    """sum_j omega^{e_j} as an exact element of Z[omega]: counts of exponents 0, 1, 2 -> (n0 - n2, n1 - n2)."""
    n = [0, 0, 0]
    for e in exps: n[e % 3] += 1
    return (n[0] - n[2], n[1] - n[2])

def weight_M(C, Cset, xL, xR, sector, kappa=0, torsion=False):
    n = len(xL)
    diff = tuple((a - b) % 3 for a, b in zip(xL, xR))
    if diff not in Cset: return (0, 0)
    a = diff
    exps = []
    for cp in C:
        e = 0
        for i in range(n):
            if sector == "NS": e += cp[i] * xL[i]
            else: e += cp[i] * (xL[i] + a[i] + kappa)
            if torsion: e += a[i] * cp[i]
        exps.append(e)
    s = zsum_omega(exps)
    assert s[0] % len(C) == 0 and s[1] % len(C) == 0
    return (s[0] // len(C), s[1] // len(C))

def weight_V(Cset, xL, xR):
    return (1, 0) if (tuple(x % 3 for x in xL) in Cset and tuple(x % 3 for x in xR) in Cset) else (0, 0)

def check(name, G, points, kappas=(0,), torsion=False):
    C = codewords(G); Cset = set(C)
    results = {}
    for sector in ("NS", "R"):
        for kappa in (kappas if sector == "R" else (0,)):
            bad = 0
            for xL, xR in points:
                wm = weight_M(C, Cset, xL, xR, sector, kappa, torsion)
                wv = weight_V(Cset, xL, xR)
                bad += (wm != wv)          # the sign (-1)^F and the momenta are common to both sides
            results[(sector, kappa)] = bad
    return results

if __name__ == "__main__":
    tetra = [[1, 0, 1, 1], [0, 1, 1, 2]]
    golay = [[1 if i == r else 0 for i in range(6)] + [int(ch) for ch in s] for r, s in
             enumerate(["011111", "101221", "110122", "121012", "122101", "112210"])]
    box = list(itertools.product(range(-2, 3), repeat=4))
    pts_t = [(xL, xR) for xL in box for xR in box]
    rng = random.Random(20261023)
    Cg = codewords(golay)
    pts_g = []
    for _ in range(20000):                                   # bias toward the orbifold's support: half from C x C classes
        if rng.random() < 0.5:
            u, w = rng.choice(Cg), rng.choice(Cg)
            xL = tuple(ui + 3 * rng.randint(-1, 1) for ui in u); xR = tuple(wi + 3 * rng.randint(-1, 1) for wi in w)
        else:
            xL = tuple(rng.randint(-3, 3) for _ in range(12)); xR = tuple(rng.randint(-3, 3) for _ in range(12))
        pts_g.append((xL, xR))
    for name, G, pts in (("tetracode (all points of the box [-2,2]^8)", tetra, pts_t),
                         ("Golay (20,000 sample points)", golay, pts_g)):
        r = check(name, G, pts, kappas=(0, 1, 2))
        rt = check(name, G, pts, kappas=(0,), torsion=True)
        print(f"{name}: mismatches NS: {r[('NS', 0)]}, R (LG phase): {r[('R', 0)]}, "
              f"R with extra omega^(kappa b): kappa=1 -> {r[('R', 1)]}, kappa=2 -> {r[('R', 2)]}; "
              f"with per-field bilinear torsion omega^(ab): NS {rt[('NS', 0)]}, R {rt[('R', 0)]}  (of {len(pts)})")
