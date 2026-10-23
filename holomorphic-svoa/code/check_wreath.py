# Copied from claude/check_wreath.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Check: Arf (x) Sym2( Sym2(T (x) Arf) (x) Arf )  ==  D4-orbifold of T^4 with a
specific twist, as index polynomials in I = Index(T); and compare the two
D4 twists found by perm_orbifolds.py.  Exact rational arithmetic."""
import json
from fractions import Fraction
from pathlib import Path
import sympy
from perm_orbifolds import (closure, compose, orbits, hom_to_z2,
                            lift_commutator_sign)

I = sympy.symbols('I')
f2 = lambda x: sympy.expand((x**2 + 3*x) / 2)          # Koszul Sym^2 index map
chain = sympy.expand(-f2(-f2(-I)))                    # Arf.Sym2(Arf.Sym2(Arf.T)) with final Arf
# D4 = <(01), (23), (02)(13)> acting on 4 copies, blocks {0,1},{2,3}
n = 4
D4 = list(closure([(1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)], n))
pairs = [(g, h) for g in D4 for h in D4 if compose(g, h) == compose(h, g)]
swaps_blocks = lambda g: int(g[0] in (2, 3))
def poly(eps):
    return sympy.expand(sum(sympy.Rational(eps(g, h), 8) * I**orbits([g, h], n) for g, h in pairs))
chi_out = lambda g, h: -1 if (swaps_blocks(g) or swaps_blocks(h)) else 1
d4_chain = sympy.expand(-poly(chi_out).subs(I, -I))   # seed Arf, chi_outer twist, final Arf
res = {"chain_polynomial": str(chain), "D4_chi_outer_polynomial": str(d4_chain),
       "equal": sympy.simplify(chain - d4_chain) == 0,
       "chain_at_I=4": int(chain.subs(I, 4))}
# compare all D4 twists: sign functions on commuting pairs
homs = hom_to_z2(D4, n)
tor = {pr: lift_commutator_sign(*pr) for pr in pairs}
table = {}
for ci, chi in enumerate(homs):
    for t in (0, 1):
        eps = lambda g, h, chi=chi, t=t: (-1 if (chi[g] or chi[h]) else 1) * (tor[(g, h)] if t else 1)
        p = poly(eps)
        table[f"chi{ci}" + ("+tor" if t else "")] = {
            "is_chi_outer": all(chi[g] == swaps_blocks(g) for g in D4),
            "poly": str(p), "at_4": int(p.subs(I, 4)), "at_-4": int(p.subs(I, -4))}
res["D4_twists"] = table
# Arf on all four seeds = sign (-1)^{orb}; is it a (chi, tor) twist?
arf_sign = {pr: (-1) ** (4 - orbits(list(pr), n)) for pr in pairs}
match = []
for ci, chi in enumerate(homs):
    for t in (0, 1):
        if all(arf_sign[pr] == (-1 if (chi[pr[0]] or chi[pr[1]]) else 1) * (tor[pr] if t else 1) for pr in pairs):
            match.append(f"chi{ci}" + ("+tor" if t else ""))
res["seed_Arf_equals_twist_on_RR_tori"] = match
(Path(__file__).resolve().parent / "out" / "wreath_check.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
