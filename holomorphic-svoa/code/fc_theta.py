# Copied from claude/paper/fc_theta.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Theta series of Lambda(C) (norm <= 8), its fixed sublattice L = Lambda^{g^2} = sqrt3 Lambda(Cbar), the dual L*,
and the shifted cosets, from the weight enumerators (fc_weight_enum.json).  Prints counts by norm."""
import json
from fractions import Fraction
import fc_characters as fc
WE = fc.WE; U = fc.U


def by_norm(s, unit):     # exponent e (1/72 units of q^{v^2/2}) -> norm v^2 = 2 e / 72
    return {str(Fraction(2 * e, U)): c[0] for e, c in sorted(s.d.items()) if Fraction(2 * e, U) <= 8}


one = fc.one; par = fc.par; cs = fc.coord_series; te = fc.theta_from_enum
out = {
 "Lambda": by_norm(te(WE["plain"], 72, cs(0, {0}, one, 12), cs(0, {1}, one, 12)), 12),
 "Lambda+chi/2": by_norm(te(WE["plain"], 72, cs(3, {0}, one, 12), cs(3, {1}, one, 12)), 12),
 "L=Lambda^{g^2}": by_norm(te(WE["Cbar_weight_enumerator"], 24, cs(0, {0}, one, 36), cs(0, {1}, one, 36)), 36),
 "L*": by_norm(te(WE["Cbar_weight_enumerator"], 24, cs(0, {0}, one, 4), cs(0, {1}, one, 4)), 4),
 "L*+chi_L/2": by_norm(te(WE["Cbar_weight_enumerator"], 24, cs(3, {0}, one, 4), cs(3, {1}, one, 4)), 4),
}
json.dump(out, open(fc.HERE / "out" / "fc_theta.json", "w"), indent=1)
for k, v in out.items():
    print(k, v)
