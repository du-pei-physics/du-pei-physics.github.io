#!/usr/bin/env python3
"""Z3 self-duality of the c = 1 N = 2 minimal model at level 1 (the IR of the LG model W = Phi^3), Section 4.

The fermionic theory is the diagonal U(1)_3 theory with odd self-dual Lorentzian momentum lattice
    Gamma = { ((3a + j)/sqrt3, (3b + j)/sqrt3) : a, b in Z, j in Z/3 },   p.p' = p_L p'_L - p_R p'_R,
characteristic vector w = (sqrt3, sqrt3) (spin structures).  The LG Z3, Phi -> omega Phi, acts by e^{2 pi i v.p} with
v = (1/sqrt3, sqrt3): v.p = j/3 mod 1, and 3 v.v = -8 is even, so v generates an honest Z3 in every sector (the
representative (1/sqrt3, 0) gives the same untwisted action but 3 v.v = 1 odd, i.e. g^3 = (-1)^F on twisted sectors).
Torus partition functions in the NS-NS spin structure with Z3 backgrounds (a, b) (eta factors omitted; they are the same
on both sides):
    Z(a, b) = e^{-i pi (a v).(b v)} sum_{p in Gamma + a v} e^{2 pi i (b v).p} q^{p_L^2/2} qbar^{p_R^2/2}.
Check: the Z3-gauged theory, Zhat(r, t) = (1/3) sum_{a,b} omega^{t a - r b} Z(a, b), equals Z(r, t) for all nine
backgrounds (exact self-duality in NS-NS), at several values of tau.  The other even spin structures are SL(2,Z) images
of NS-NS; the periodic-periodic sector is the protected table Z(a, b) = +-(3 delta_{a0} delta_{b0} - 1), whose Fourier
transform is minus itself (self-dual only up to Arf).
Exact form: with Theta_s(tau) = sum_m q^{(3m+s)^2/6} (so Theta_{s+3} = Theta_s = Theta_{-s}),
    Z(a, b) = sum_{j in Z/3} omega^{b(j+2a)} Theta_{j+a} conj(Theta_j),
and summing over b in Zhat(r, t) sets j = r + a; Theta_{r+2a} = Theta_{2r+a} and the substitution j -> j - r give
Zhat(r, t) = Z(r, t) identically in tau.  The script checks (1) the lattice sums against this closed form numerically,
(2) the nine NS identities and the nine periodic-periodic identities exactly, as polynomials in Theta_0, Theta_1 and
their conjugates over Z[omega], and (3) the lattice-sum identity numerically at four values of tau.
usage: u13_duality.py
"""
import cmath, itertools, math
import numpy as np

r3 = 3 ** 0.5
v = np.array([1 / r3, r3])
assert abs(3 * (v[0] ** 2 - v[1] ** 2) + 8) < 1e-12          # 3 v.v = -8
M = 9
G = np.array([((3 * a + j) / r3, (3 * b + j) / r3) for a, b in itertools.product(range(-M, M + 1), repeat=2) for j in range(3)])
om = cmath.exp(2j * math.pi / 3)


def Z(a, b, tau):
    s, t = a * v, b * v
    P = G + s
    ph = np.exp(2j * np.pi * (P[:, 0] * t[0] - P[:, 1] * t[1]))
    q = np.exp(1j * np.pi * tau * P[:, 0] ** 2) * np.exp(-1j * np.pi * np.conj(tau) * P[:, 1] ** 2)
    return np.exp(-1j * np.pi * (s[0] * t[0] - s[1] * t[1])) * np.sum(ph * q)


def theta(s, tau):
    m = np.arange(-M, M + 1)
    return np.sum(np.exp(1j * np.pi * tau * (3 * m + s) ** 2 / 3))


# (1) closed form against the lattice sums
cf = 0.0
for tau in (0.31 + 0.87j, -0.42 + 1.3j):
    for a, b in itertools.product(range(3), repeat=2):
        closed = sum(om ** (b * (j + 2 * a)) * theta(j + a, tau) * np.conj(theta(j, tau)) for j in range(3))
        cf = max(cf, abs(closed - Z(a, b, tau)))
print("closed form Z(a,b) = sum_j omega^{b(j+2a)} Theta_{j+a} conj(Theta_j) matches the lattice sums to",
      f"{cf:.1e}", "->", "pass" if cf < 1e-12 else "FAIL")


# (2) exact identities over Z[omega]; x + y omega is the pair (x, y), omega^2 = -1 - omega
def zmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0] - u[1] * v[1])


def wpow(n):
    return ((1, 0), (0, 1), (-1, -1))[n % 3]


def red(s):                                      # Theta_s depends on s only through {s, -s} mod 3
    return min(s % 3, -s % 3)


def add(P, key, c):
    x = P.get(key, (0, 0))
    P[key] = (x[0] + c[0], x[1] + c[1])


def clean(P):
    return {k: v for k, v in P.items() if v != (0, 0)}


def Zx(a, b):                                    # monomial (s, s') stands for Theta_s conj(Theta_s')
    P = {}
    for j in range(3):
        add(P, (red(j + a), red(j)), wpow(b * (j + 2 * a)))
    return clean(P)


ns_ok = rr_ok = 0
for r, t in itertools.product(range(3), repeat=2):
    lhs = {}
    for a, b in itertools.product(range(3), repeat=2):
        for k, c in Zx(a, b).items():
            add(lhs, k, zmul(wpow(t * a - r * b), c))
    ns_ok += clean(lhs) == {k: (3 * c[0], 3 * c[1]) for k, c in Zx(r, t).items()}
    I = lambda a, b: 3 * (a == 0 and b == 0) - 1
    rr = (0, 0)
    for a, b in itertools.product(range(3), repeat=2):
        w = wpow(t * a - r * b)
        rr = (rr[0] + I(a, b) * w[0], rr[1] + I(a, b) * w[1])
    rr_ok += rr == (-3 * I(r, t), 0)
print(f"exact over Z[omega]: NS identities Zhat = Z {ns_ok}/9, periodic-periodic identities Ihat = -I {rr_ok}/9 ->",
      "pass" if ns_ok == rr_ok == 9 else "FAIL")

# (3) the lattice-sum identity at four values of tau
worst = 0.0
for tau in (0.31 + 0.87j, -0.42 + 1.3j, 0.05 + 0.62j, 0.49 + 0.9j):
    tab = {(a, b): Z(a, b, tau) for a in range(3) for b in range(3)}
    err = max(abs(sum(om ** (t * a - r * b) * tab[(a, b)] for a in range(3) for b in range(3)) / 3 - tab[(r, t)])
              for r in range(3) for t in range(3))
    worst = max(worst, err)
    print(f"tau = {tau}: max over the nine backgrounds of |Zhat(r,t) - Z(r,t)| = {err:.1e};  Z(0,0) = {tab[(0, 0)].real:.6f}")
print("NS-NS Z3 self-duality holds to", f"{worst:.1e}", "->", "pass" if worst < 1e-12 else "FAIL")
