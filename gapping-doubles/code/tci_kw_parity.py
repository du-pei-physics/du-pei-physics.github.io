#!/usr/bin/env python3
"""Kramers-Wannier parity of the tricritical Ising primaries (exact, sympy).

M(4,5): h_{r,s} = ((5r - 4s)^2 - 1)/80, 1 <= r <= 3, 1 <= s <= 4; the Z2 line is phi_{3,1} (h = 3/2), the
Kramers-Wannier line N is phi_{2,1} (h = 7/16), N x N = 1 + phi_{3,1}.  A Verlinde line L_a acting on the primary
phi_b has eigenvalue S_{ab}/S_{1b}; phi_b commutes with N iff the eigenvalue equals d_N = S_{N1}/S_{11} = sqrt 2
(then the deformation by phi_b preserves N), and anticommutes iff it equals -sqrt 2.
  S_{(r,s),(rho,sigma)} = 2 sqrt(2/(p p')) (-1)^{1 + s rho + r sigma} sin(pi p' r rho / p) sin(pi p s sigma / p'),
  p = 4, p' = 5 (Di Francesco-Mathieu-Senechal (10.134)).
usage: tci_kw_parity.py
"""
import sympy as sp
p, pp = 4, 5
def S(r, s, rho, sig):
    return 2 * sp.sqrt(sp.Rational(2, p * pp)) * (-1) ** (1 + s * rho + r * sig) * sp.sin(sp.pi * pp * r * rho / p) * sp.sin(sp.pi * p * s * sig / pp)
def h(r, s): return sp.Rational((pp * r - p * s) ** 2 - (pp - p) ** 2, 4 * p * pp)
N = (2, 1); one = (1, 1)
dN = sp.nsimplify(sp.simplify(S(*N, *one) / S(*one, *one)))
print("d_N =", dN, " h_N =", h(*N), " h_eta =", h(3, 1))
for b, name in [((1, 2), "epsilon = phi_{1,2}"), ((1, 3), "epsilon' = phi_{1,3}"), ((1, 4), "phi_{1,4} = phi_{3,1} (eta)"), ((2, 2), "sigma = phi_{2,2}"), ((2, 1), "sigma' = phi_{2,1}")]:
    ev = sp.nsimplify(sp.simplify(S(*N, *b) / S(*one, *b)))
    verdict = "commutes with N (KW-even)" if ev == dN else ("anticommutes with N (KW-odd)" if ev == -dN else "neither (not a symmetric operator)")
    print(f"h = {h(*b)!s:>5}  {name:28s} eigenvalue of N: {ev!s:>8}  -> {verdict}")
