# Copied from claude/paper/fc_sym2.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed, and one comment was reworded for publication
# (see tools/build_code_folder.py).
"""Low-weight content of T' = Sym^2(T) (c=144) and P = Sym^2(T' x Arf) (c=288) from the characters of T.

Untwisted sector of Sym^2(X): graded-symmetric part of X x X (swap with the Koszul sign).
  NS:  dim_h = 1/2 [ (chi*chi)_h + chi^AP(q^2)_h ],   super_h = 1/2 [ (chiAP*chiAP)_h + chi^AP(q^2)_h ]
       (Tr(sigma) = Tr((-1)^F sigma) = sum_a (-1)^{|a|} q^{2 h_a} on X x X).
  The sigma-twisted NS sector starts at c_X/12 (= 6 for X = T, 12 for X = T'), so below that weight the NS
  content is untwisted only.
Ramond ground states (L_0 = c/24):  untwisted: even = S^2 B + L^2 F, odd = B x F (B, F = even/odd ground states of X);
  twisted: one copy of the ground states of X with the same parity (Z^{PP}(sigma,1) = Z^{PP}(sigma,sigma) = I(X)).
Arf: X x Arf has the same NS sector and Ramond ground states with parities exchanged.
"""
import json
from fractions import Fraction
from pathlib import Path
HERE = Path(__file__).resolve().parent
ch = json.load(open(HERE / "out" / "fc_characters.json"))


def ser(d):
    return {Fraction(h): n for h, n in d.items()}


def mul(a, b, hmax):
    out = {}
    for h1, n1 in a.items():
        for h2, n2 in b.items():
            if h1 + h2 <= hmax:
                out[h1 + h2] = out.get(h1 + h2, 0) + n1 * n2
    return out


def sym2_ns(chi, chiAP, hmax):
    sq = mul(chi, chi, hmax); sqAP = mul(chiAP, chiAP, hmax)
    dbl = {2 * h: n for h, n in chiAP.items() if 2 * h <= hmax}
    hs = sorted(set(sq) | set(dbl))
    dim = {}; sup = {}
    for h in hs:
        x = sq.get(h, 0) + dbl.get(h, 0); y = sqAP.get(h, 0) + dbl.get(h, 0)
        assert x % 2 == 0 and y % 2 == 0
        dim[h] = x // 2; sup[h] = y // 2
    return dim, sup


def sym2_R(B, F):
    ue = B * (B + 1) // 2 + F * (F - 1) // 2; uo = B * F
    return {"untwisted": [ue, uo], "twisted": [B, F], "total": [ue + B, uo + F]}


def fmt(d):
    return {str(h): n for h, n in sorted(d.items()) if n}


HMAX = Fraction(5)
T, TAP = ser(ch["T_AA"]), ser(ch["T_AP"])
# R ground states of T: untwisted Z6-invariants (super -4) + 40 in each of the g^2, g^4 twisted sectors (super 0)
nR = int(ch["T_PA"]["3"]); sR = int(ch["T_PP"]["3"])
BT, FT = (nR + sR) // 2, (nR - sR) // 2
Tp, TpAP = sym2_ns(T, TAP, HMAX)                  # valid below 6
Pn, PnAP = sym2_ns(Tp, TpAP, HMAX)                # valid below 12
RTp = sym2_R(BT, FT)
Bp, Fp = RTp["total"]
RP = sym2_R(Fp, Bp)                               # Arf exchanges the parities of T' before squaring
res = {"T": {"NS": fmt(T), "R_ground_even_odd": [BT, FT], "index": BT - FT},
       "Tprime": {"NS": fmt(Tp), "NS_super": fmt(TpAP), "R_ground": RTp, "index": Bp - Fp},
       "P": {"NS": fmt(Pn), "NS_super": fmt(PnAP), "R_ground": RP, "index": RP["total"][0] - RP["total"][1]}}
I = BT - FT
assert res["Tprime"]["index"] == (I * I + 3 * I) // 2 == 2
J = -res["Tprime"]["index"]
assert res["P"]["index"] == (J * J + 3 * J) // 2 == -1
# D4 index polynomial: composing I -> (I^2+3I)/2, I -> -I, I -> (I^2+3I)/2 gives (I^4 + 6 I^3 + 3 I^2 - 18 I)/8 in I = I(T);
# the form (I^4 - 6 I^3 + 3 I^2 + 18 I)/8 of the original question is the same polynomial in -I(T) = +4 (NOTES 2.3 uses the index +4 seed U).
pd = Fraction(I ** 4 + 6 * I ** 3 + 3 * I ** 2 - 18 * I, 8); pb = Fraction((-I) ** 4 - 6 * (-I) ** 3 + 3 * (-I) ** 2 + 18 * (-I), 8)
assert pd == pb == -1
res["P"]["D4_index_polynomial_check"] = "(I^4+6I^3+3I^2-18I)/8 at I(T)=-4 = (I^4-6I^3+3I^2+18I)/8 at I=+4 = -1"
(HERE / "out" / "fc_sym2.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
