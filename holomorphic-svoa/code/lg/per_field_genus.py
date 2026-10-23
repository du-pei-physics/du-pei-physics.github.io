#!/usr/bin/env python3
# Copied from GappedCFT/code/per_field_genus.py (a companion project of the scft-code-search workspace) for the holomorphic-svoa
# note; only pointers to the working note were changed (see tools/package_lg.py).
"""Per-field form of the elliptic genus identity LG/C = V (x) Vbar (../../notes/lg_elliptic_genus.md sections 3-4).

LG side (one field, W = Phi^3, left U(1) fugacity y, Z3 flavor twist a in space and b in time):
    F_{a,b}(tau, z) = y^{-a/3} theta1(tau, -2z/3 + (a tau + b)/3) / theta1(tau, z/3 + (a tau + b)/3),
and its charge projections G_{a,k} = (1/3) sum_b omega^{-b k} F_{a,b}. The LG/C genus is
    Z_LG/C = (1/|C|) sum_{c, c' in C} prod_i F_{c_i, c'_i} = sum_{c in C, k in C} prod_i G_{c_i, k_i}   (C = C^perp).
Lattice side (one copy of V_{sqrt3 Z}, Ramond sector, label u in Z3 = coset of the frame coordinate):
    psi_u(tau, z) = eta^{-1} sum_{x = u mod 3} (-1)^x y^{(x + 3/2)/3} q^{(x + 3/2)^2 / 6}
(the per-field fermion-parity factor e^{i pi (x + 3/2)} = (-1)^x (-i); the constants (-i) multiply to 1 over n = 0 mod 4
coordinates and are dropped), and iota_w = the y = 1, q^0 coefficient of psi_w (the per-field Witten index of a
right-moving label: iota_0 = 0). Then Index(C) psi_V = sum_{u in C, w in C} prod_i psi_{u_i} iota_{w_i}.
The script searches GL2(F3) for a relabeling (a, k) -> (u, w) with G_{a,k} = kappa psi_u iota_w for all nine (a, k)
(kappa a common constant), comparing exact series up to q^QMAX. If one exists, the two elliptic genera agree for every
self-dual ternary code, because the relabeling maps C x C onto itself.

Series: dict {(eq, ey): coefficient}, eq = exponent of q in units of 1/24, ey = exponent of y in units of 1/6;
coefficients in Q(omega) as pairs (p, r) = p + r omega with Fractions.

usage: per_field_genus.py
"""
from fractions import Fraction as Fr
import itertools

QMAX = 2                      # compare up to q^QMAX
EQ = 24 * QMAX                # in units of 1/24
YC = 6                        # y-exponent unit 1/6

# ---- Q(omega) arithmetic, omega^2 = -1 - omega
def cmul(x, y):
    a, b = x; c, d = y
    # (a + b w)(c + d w) = ac + (ad + bc) w + bd w^2 = (ac - bd) + (ad + bc - bd) w
    return (a * c - b * d, a * d + b * c - b * d)
def cadd(x, y): return (x[0] + y[0], x[1] + y[1])
ZERO, ONE = (Fr(0), Fr(0)), (Fr(1), Fr(0))
def wpow(k):
    k %= 3
    return [ONE, (Fr(0), Fr(1)), (Fr(-1), Fr(-1))][k]
def cscale(x, s): return (x[0] * s, x[1] * s)
def is_zero(x): return x[0] == 0 and x[1] == 0

def smul(A, B):
    out = {}
    for (q1, y1), c1 in A.items():
        for (q2, y2), c2 in B.items():
            q = q1 + q2
            if q > EQ: continue
            key = (q, y1 + y2)
            out[key] = cadd(out.get(key, ZERO), cmul(c1, c2))
    return {k: v for k, v in out.items() if not is_zero(v)}
def sadd(A, B):
    out = dict(A)
    for k, v in B.items(): out[k] = cadd(out.get(k, ZERO), v)
    return {k: v for k, v in out.items() if not is_zero(v)}
def sscale(A, c): return {k: cmul(v, c) for k, v in A.items() if not is_zero(cmul(v, c))}
def mono(eq, ey, c=ONE): return {(eq, ey): c}

def geom(term):
    """1/(1 - term) for a monomial term with positive q-exponent: sum_j term^j."""
    (eq, ey), c = next(iter(term.items()))
    assert eq > 0
    out, p = {(0, 0): ONE}, {(0, 0): ONE}
    while True:
        p = smul(p, term)
        if not p: break
        out = sadd(out, p)
    return out

def F_ab(a, b):
    """y^{-a/3} theta1(-2z/3 + u)/theta1(z/3 + u), u = (a tau + b)/3, via the product formula.
    y_N = y^{-2/3} q^{a/3} omega^b, y_D = y^{1/3} q^{a/3} omega^b; ratio = y^{1/2} (y_N - 1)/(y_D - 1) prod_n (...)."""
    yN = (8 * a, -4, wpow(b)); yD = (8 * a, 2, wpow(b))      # (eq in 1/24, ey in 1/6, coeff)
    pref = mono(0, 3 - 2 * a)                                 # y^{1/2} y^{-a/3}
    if a == 0:
        # (y_N - 1)/(y_D - 1) is a Laurent polynomial in t = y^{1/3}: solve by exact division of t^2 (y_N - 1) by (y_D - 1)
        # numerator: omega^b t^{-2} - 1 = t^{-2}(omega^b - t^2); denominator: omega^b t - 1
        # (omega^b - t^2)/(omega^b t - 1): for b = 0: (1 - t^2)/(t - 1) = -(1 + t); b != 0: -(omega^{2b} t + omega^b)
        if b == 0: poly = {(0, 0): (Fr(-1), Fr(0)), (0, 2): (Fr(-1), Fr(0))}
        else: poly = {(0, 2): cscale(wpow(2 * b), Fr(-1)), (0, 0): cscale(wpow(b), Fr(-1))}
        lead = smul(mono(0, -4), poly)                        # t^{-2} * poly
    else:
        # (y_N - 1)/(y_D - 1) = (1 - y_N)/(1 - y_D)
        lead = smul(sadd(mono(0, 0), mono(yN[0], yN[1], cscale(yN[2], Fr(-1)))), geom(mono(yD[0], yD[1], yD[2])))
    out = smul(pref, lead)
    for nq in range(1, QMAX + 1):
        e = 24 * nq
        for (eq, ey, c) in (yN,):
            out = smul(out, sadd(mono(0, 0), mono(e + eq, ey, cscale(c, Fr(-1)))))            # (1 - q^n y_N)
            cinv = wpow(-b)
            if e - eq > 0 or True:
                out = smul(out, sadd(mono(0, 0), mono(e - eq, -ey, cscale(cinv, Fr(-1)))))   # (1 - q^n / y_N)
        out = smul(out, geom(mono(e + yD[0], yD[1], yD[2])))                                  # 1/(1 - q^n y_D)
        out = smul(out, geom(mono(e - yD[0], -yD[1], wpow(-b))))                              # 1/(1 - q^n / y_D)
    return out

def G_ak(a, k):
    tot = {}
    for b in range(3):
        tot = sadd(tot, sscale(F_ab(a, b), wpow(-b * k)))
    return {key: cscale(v, Fr(1, 3)) for key, v in tot.items()}

def eta_inv():
    """q^{-1/24} prod (1 - q^n)^{-1}, as a series in units of 1/24 (shifted later)."""
    out = {(0, 0): ONE}
    for n in range(1, QMAX + 2):
        out = smul(out, geom(mono(24 * n, 0)))
    return out

def psi(u):
    """eta^{-1} sum_{x = u mod 3} (-1)^x y^{(x+3/2)/3} q^{(x+3/2)^2/6}; exponents: q in 1/24, y in 1/6."""
    out = {}
    for x in range(-12, 13):
        if x % 3 != u: continue
        h = Fr(2 * x + 3, 2) ** 2 / 6 - Fr(1, 24)            # with the eta^{-1} q^{-1/24} shift
        eq = h * 24
        if eq > EQ: continue
        ey = Fr(2 * x + 3, 2) / 3 * YC
        assert eq.denominator == 1 and ey.denominator == 1
        out[(int(eq), int(ey))] = cadd(out.get((int(eq), int(ey)), ZERO), ((-1) ** (x % 2) * ONE[0], Fr(0)))
    return smul(out, eta_inv())

def index_of(s):
    """y = 1, q^0 coefficient (lowest level 0)."""
    tot = ZERO
    for (eq, ey), c in s.items():
        if eq == 0: tot = cadd(tot, c)
    return tot

def proportional(A, B):
    """A = kappa B for a constant kappa in Q(omega)? returns kappa or None (both empty -> 'zero')."""
    if not A and not B: return "zero"
    if not A or not B or set(A) != set(B): return None
    k0 = next(iter(B)); b0 = B[k0]
    # kappa = A[k0]/B[k0]; division in Q(omega): x / y = x * conj(y) / |y|^2, conj(omega) = omega^2
    def conj(y): return (y[0] - y[1], -y[1])
    def norm(y): return y[0] * y[0] - y[0] * y[1] + y[1] * y[1]
    kappa = cscale(cmul(A[k0], conj(b0)), 1 / norm(b0))
    return kappa if all(A[k] == cmul(kappa, B[k]) for k in A) else None


if __name__ == "__main__":
    G = {(a, k): G_ak(a, k) for a in range(3) for k in range(3)}
    P = {u: psi(u) for u in range(3)}
    iota = {w: index_of(P[w]) for w in range(3)}
    print("per-field right indices iota_w:", {w: (str(v[0]), str(v[1])) for w, v in iota.items()})
    print("zero pattern of G_{a,k} (a, k):", [ak for ak, g in G.items() if not g])
    found = []
    for (p, q, r, s) in itertools.product(range(3), repeat=4):
        if (p * s - q * r) % 3 == 0: continue
        kappas = set(); ok = True
        for (a, k), g in G.items():
            u, w = (p * a + q * k) % 3, (r * a + s * k) % 3
            target = sscale(P[u], iota[w]) if not is_zero(iota[w]) else {}
            kap = proportional(g, target)
            if kap is None: ok = False; break
            if kap != "zero": kappas.add(kap)
        if ok and len(kappas) <= 1:
            found.append(((p, q, r, s), kappas))
    print(f"relabelings L = [[p, q], [r, s]] in GL2(F3) with G_(a,k) = kappa psi_u iota_w up to q^{QMAX}:")
    for L, kap in found:
        print("   ", L, "kappa =", [(str(k[0]), str(k[1])) for k in kap])
    print("found:", len(found))
