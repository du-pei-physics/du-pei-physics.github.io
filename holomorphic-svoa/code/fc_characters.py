# Copied from claude/paper/fc_characters.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Exact low-weight characters of V = V_{Lambda(C)} (c=72), T = V//Z_6, T' = Sym^2 T (c=144), P = Sym^2(T' x Arf) (c=288).

Conventions (see field_content.md, Section 2):
 * Series in z = q^{1/72}; a series stores exponents of q^{L_0 - c/24} in units of 1/72.  Coefficients lie in Z[w],
   w = exp(2 pi i/3), stored as integer pairs (a, b) = a + b w.
 * Z^{(a,b)}_u = Tr_{g^a-twisted sector, spin u}( g^b q^{L_0-c/24} ), u in {AA, AP, PA, PP} = (space, time);
   A = antiperiodic.  AA = NS without (-1)^F, AP = NS with (-1)^F, PA = R without, PP = R with (-1)^F.
   Spin labels as vectors of F_2^2: AA=(1,1), AP=(1,0), PA=(0,1), PP=(0,0); modular group acts on the right on both
   the holonomy pair (a,b) and the spin label: Z^{x}_u(gamma tau) = Z^{x gamma}_{u gamma}(tau).
 * Inputs (untwisted twining functions, exact):
     (0,0): Theta_Lambda / eta^72 with (-1)^F = e^{pi i <v,chi>};  R: Theta_{Lambda+chi/2}, (-1)^F = (-1)^{#2's} on ground states.
     (0,j), j=1,5:  Delta(3tau)/Delta(6tau)  (Frame shape 6^24 3^-24, no fixed momenta; vanishes in R).
     (0,3):        Delta(tau)^3/Delta(2tau)^3 (g^3 = -1; vanishes in R).
     (0,2),(0,4):  Theta_{L,(.)}/Delta(3tau), L = Lambda^{g^2} = sqrt3 Lambda(Cbar); standard lift (trivial on e^v, v in L,
                   and on the g^2-fixed Ramond states: GJF, assumption A1 in the write-up).
   Twisted sectors by S, T and, for pairs generating Z_6, by the Hermite-normal-form transformation of Delta(k gamma tau).
"""
import itertools, json, math, sys
from fractions import Fraction
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
U = 72                    # exponent unit: q^(1/72)
TGT = 2 * U               # report exponents (L_0 - c/24) <= 2, i.e. L_0 <= 5 for c = 72
CUT = TGT + 4 * U         # per-coordinate theta cutoff (the deepest denominator is q^3)
WE = json.load(open(HERE / "out" / "fc_weight_enum.json"))


# ---------------------------------------------------------------- Z[w] truncated Laurent series in q^(1/72)
INF = 10 ** 9


class S:
    """Sparse series sum c_e q^(e/72), exact for e <= prec."""
    def __init__(self, d=None, prec=INF):
        self.prec = prec
        self.d = {} if d is None else {e: c for e, c in d.items() if e <= prec and (c[0] or c[1])}

    @staticmethod
    def mono(e, c=(1, 0)):
        return S({e: c})

    def __add__(self, o):
        d = dict(self.d)
        for e, (a, b) in o.d.items():
            x, y = d.get(e, (0, 0)); d[e] = (x + a, y + b)
        return S(d, min(self.prec, o.prec))

    def __neg__(self):
        return S({e: (-a, -b) for e, (a, b) in self.d.items()}, self.prec)

    def __sub__(self, o):
        return self + (-o)

    def lead(self):
        return min(self.d) if self.d else INF

    def __mul__(self, o):
        if isinstance(o, int):
            return S({e: (o * a, o * b) for e, (a, b) in self.d.items()}, self.prec)
        prec = min(self.prec + o.lead(), o.prec + self.lead())
        d = {}
        for e1, (a, b) in self.d.items():
            for e2, (c, dd) in o.d.items():
                e = e1 + e2
                if e > prec:
                    continue
                x, y = d.get(e, (0, 0))
                d[e] = (x + a * c - b * dd, y + a * dd + b * c - b * dd)
        return S(d, prec)

    def inv(self):
        e0 = self.lead(); c0 = self.d[e0]; a, b = c0
        assert a * a - a * b + b * b == 1, "leading coefficient must be a unit"
        ic = (a - b, -b)                                     # (a + b w)^-1
        K = self.prec - 2 * e0                               # result exact up to exponent K
        f = {e - e0: c for e, c in self.d.items() if e > e0}
        g = {}
        for k in range(0, K + e0 + 1):
            s = (0, 0)
            for j, fj in f.items():
                if j <= k and (k - j) in g:
                    t = mul_c(fj, g[k - j]); s = (s[0] + t[0], s[1] + t[1])
            val = mul_c(ic, ((1 if k == 0 else 0) - s[0], -s[1]))
            if val != (0, 0):
                g[k] = val
        return S({k - e0: c for k, c in g.items()}, K)

    def tshift(self, k=1):
        """tau -> tau + k: multiply coefficient of q^(e/72) by exp(2 pi i k e/72); requires a sixth root of unity."""
        d = {}
        for e, c in self.d.items():
            assert e > TGT or (k * e) % 12 == 0, ("T-phase not in mu_6", e)
            if e <= TGT:
                d[e] = mul_c(c, zeta6((k * e) // 12))
        return S(d, min(self.prec, TGT))

    def iszero(self, upto):
        assert self.prec >= upto, (self.prec, upto)
        return all(e > upto for e in self.d)


def mul_c(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


def zeta6(k):
    """exp(2 pi i k/6) in Z[w]."""
    return [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)][k % 6]


def zetaD(D, m):
    assert 6 % D == 0
    return zeta6((6 // D) * m)


def powS(f, n):
    r = S.mono(0)
    for _ in range(n):
        r = r * f
    return r


# ---------------------------------------------------------------- eta products
NT = 60
TAU = [0] * (NT + 2)              # Ramanujan tau via prod (1-q^n)^24
_p = [1] + [0] * (NT + 1)
for n in range(1, NT + 2):
    for _ in range(24):
        _p = [_p[i] - (_p[i - n] if i >= n else 0) for i in range(NT + 2)]
for m in range(1, NT + 2):
    TAU[m] = _p[m - 1]


def Delta(A, B, D):
    """Delta((A tau + B)/D) = sum_m tau(m) zeta_D^{Bm} q^{Am/D}."""
    d = {}
    for m in range(1, NT + 2):
        e = U * A * m // D
        assert (U * A * m) % D == 0
        d[e] = mul_c(zetaD(D, B * m), (TAU[m], 0))
    return S(d, U * A * (NT + 2) // D - 1)


def hnf(N):
    """N = M' U with M' in SL2(Z), U = [[A,B],[0,D]], D > 0, 0 <= B < D."""
    (a, b), (c, d) = N
    g = math.gcd(a, c)
    r, s = -c // g, a // g
    # find p,q with p s - q r = 1
    def egcd(x, y):
        if y == 0:
            return (x, 1, 0)
        gg, u, v = egcd(y, x % y); return (gg, v, u - (x // y) * v)
    gg, p, q = egcd(s, -r)            # p s + q (-r) = gg = +-1
    if gg < 0:
        p, q = -p, -q
    Minv = [[p, q], [r, s]]
    A_ = p * a + q * c; B_ = p * b + q * d; C_ = r * a + s * c; D_ = r * b + s * d
    assert C_ == 0 and A_ * D_ == a * d - b * c
    if D_ < 0:
        A_, B_, D_ = -A_, -B_, -D_
    # right-multiplying U by [[1,t],[0,1]] is not allowed; instead adjust M': U -> [[1,t],[0,1]] U changes B by t D
    t = (B_ // D_); B_ -= t * D_
    return A_, B_, D_


def eta_ratio_at(gamma, ks_num, ks_den, power=1):
    """prod_{k in ks_num} Delta(k gamma tau) / prod_{k in ks_den} Delta(k gamma tau), each to `power`;
    Delta(k gamma tau) = ((c tau + d)/D_k)^12 Delta(U_k tau); the (c tau + d)^12 factors cancel (equal counts)."""
    (al, be), (c, d) = gamma
    assert al * d - be * c == 1 and len(ks_num) == len(ks_den)
    num = S.mono(0); den = S.mono(0); fac = Fraction(1)
    for k in ks_num:
        A, B, D = hnf([[k * al, k * be], [c, d]]); num = num * Delta(A, B, D); fac /= Fraction(D) ** 12
    for k in ks_den:
        A, B, D = hnf([[k * al, k * be], [c, d]]); den = den * Delta(A, B, D); fac *= Fraction(D) ** 12
    r = num * den.inv()
    r = powS(r, power)
    fac = fac ** power
    assert fac.denominator == 1
    return r * int(fac)


# ---------------------------------------------------------------- lattice theta series from weight enumerators
def coord_series(shift2, residues, sign, unit):
    """sum over integers x with x = r mod 3 (r in residues) of sign(x) q^{unit * (2x + shift2)^2 / 4}, in 1/72 units.
    shift2 = 0 for no shift (x), 3 for x + 3/2."""
    d = {}
    for x in range(-40, 41):
        if x % 3 not in residues:
            continue
        e = Fraction(unit * (2 * x + shift2) ** 2, 4)
        assert e.denominator == 1
        e = int(e)
        if e > CUT:
            continue
        a, _ = d.get(e, (0, 0)); d[e] = (a + sign(x), 0)
    return S(d, CUT)


def theta_from_enum(enum, n, f0, f1):
    """sum_w enum[w] f0^(n-w) f1^w (truncated)."""
    tot = S(); p0 = [S.mono(0)]; p1 = [S.mono(0)]
    for _ in range(n):
        p0.append(p0[-1] * f0); p1.append(p1[-1] * f1)
    for w, A in enumerate(enum):
        if A:
            tot = tot + p0[n - w] * p1[w] * int(A)
    return tot


one = lambda x: 1
par = lambda x: -1 if x % 2 else 1          # (-1)^x


def main():
    out = {}
    Aw = WE["plain"]; Asw = WE["signed"]; Cb = WE["Cbar_weight_enumerator"]; Cbs = WE["Cbar_signed_enumerator"]
    # Lambda: norm v^2 = sum x^2 / 3; q^{v^2/2} = q^{x^2/6} -> 12 x^2 units of 1/72; unit*(2x)^2/4 = unit*x^2 => unit 12
    th_L_AA = theta_from_enum(Aw, 72, coord_series(0, {0}, one, 12), coord_series(0, {1}, one, 12))
    th_L_AP = theta_from_enum(Aw, 72, coord_series(0, {0}, par, 12), coord_series(0, {1}, par, 12))
    # R coset Lambda + chi/2: coordinates (x + 3/2)/sqrt3 -> q^{(x+3/2)^2/6} = 12 (2x+3)^2/4 units
    th_R_PA = theta_from_enum(Aw, 72, coord_series(3, {0}, one, 12), coord_series(3, {1}, one, 12))
    eta72inv = Delta(1, 0, 1).inv(); eta72inv = powS(eta72inv, 3)          # 1/eta^72 = Delta^-3
    Z = {}
    Z[(0, 0)] = {"AA": th_L_AA * eta72inv, "AP": th_L_AP * eta72inv, "PA": th_R_PA * eta72inv,
                 "PP": S.mono(0, (-24, 0))}   # Witten index of V: constant -24 (Index(C), computed by DP)
    # (0,2), (0,4): L = {b in Z^24 : b mod 3 in Cbar}, v^2 = sum b^2 -> q^{b^2/2} = 36 b^2 units  (unit 36)
    LAA = theta_from_enum(Cb, 24, coord_series(0, {0}, one, 36), coord_series(0, {1}, one, 36))
    LAP = theta_from_enum(Cb, 24, coord_series(0, {0}, par, 36), coord_series(0, {1}, par, 36))
    LPA = theta_from_enum(Cb, 24, coord_series(3, {0}, one, 36), coord_series(3, {1}, one, 36))
    # PP: sign prod (-1)^{b+1} per coordinate; residue-2 coordinates give minus the residue-1 series -> signed enumerator
    sg = lambda x: -1 if (x + 1) % 2 else 1
    h1 = coord_series(3, {1}, sg, 36); h2 = coord_series(3, {2}, sg, 36)
    assert (h1 + h2).iszero(CUT)
    LPP = theta_from_enum(Cbs, 24, coord_series(3, {0}, sg, 36), h1)
    D3inv = Delta(3, 0, 1).inv()
    for j in (2, 4):
        Z[(0, j)] = {"AA": LAA * D3inv, "AP": LAP * D3inv, "PA": LPA * D3inv, "PP": LPP * D3inv}
    out["check_PP_02_vanishes"] = (LPP * D3inv).iszero(TGT)
    # (2,0),(4,0) by S: L* = L/3 -> w^2/2 = b^2/18 = 4 b^2 units (unit 4); L* + (1/2)1 -> (2b+3)^2/72 units
    Ls_AA = theta_from_enum(Cb, 24, coord_series(0, {0}, one, 4), coord_series(0, {1}, one, 4))
    Ls_PA = theta_from_enum(Cb, 24, coord_series(3, {0}, one, 4), coord_series(3, {1}, one, 4))
    Ls_AP = theta_from_enum(Cb, 24, coord_series(0, {0}, par, 4), coord_series(0, {1}, par, 4))   # (-1)^{sum b}
    D13inv = Delta(1, 0, 3).inv()                       # 1/Delta(tau/3) = 1/eta(tau/3)^24
    for a in (2, 4):
        Z[(a, 0)] = {"AA": Ls_AA * D13inv, "PA": Ls_PA * D13inv, "AP": Ls_AP * D13inv, "PP": S()}
    # T-images: Z^{(a,b)}_u (tau+1) = Z^{(a,a+b)}_{uT}(tau); uT: AA<->AP, PA->PA, PP->PP
    swapT = {"AA": "AP", "AP": "AA", "PA": "PA", "PP": "PP"}
    for a in (2, 4):
        prev = (a, 0)
        for _ in range(2):
            nxt = (a, (prev[1] + a) % 6)
            Z[nxt] = {swapT[u]: Z[prev][u].tshift() for u in Z[prev]}
            prev = nxt
        # level matching: one more T returns to (a,0)
        back = {swapT[u]: Z[prev][u].tshift() for u in Z[prev]}
        out[f"level_matching_{a}"] = all((back[u] - Z[(a, 0)][u]).iszero(TGT) for u in back)
    # (0,3), (3,0), (3,3)
    F3 = powS(Delta(1, 0, 1) * Delta(2, 0, 1).inv(), 3)
    Z[(0, 3)] = {"AA": F3, "AP": F3, "PA": S(), "PP": S()}
    G3 = powS(Delta(1, 0, 1) * Delta(1, 0, 2).inv(), 3) * (2 ** 36)
    Z[(3, 0)] = {"AA": G3, "PA": G3, "AP": S(), "PP": S()}
    Z[(3, 3)] = {swapT[u]: Z[(3, 0)][u].tshift() for u in Z[(3, 0)]}
    back = {swapT[u]: Z[(3, 3)][u].tshift() for u in Z[(3, 3)]}
    out["level_matching_3"] = all((back[u] - Z[(3, 0)][u]).iszero(TGT) for u in back)
    # order-6 pairs: Z^{(a,b)}_v(tau) = Z^{(0,1)}_{v gamma^-1}(gamma tau), gamma second row (a,b)
    spins = {"AA": (1, 1), "AP": (1, 0), "PA": (0, 1), "PP": (0, 0)}; inv_sp = {v: k for k, v in spins.items()}
    consistency = []
    for a in range(6):
        for b in range(6):
            if math.gcd(math.gcd(a, b), 6) != 1:
                continue
            found = []
            for c in (a, a + 6, a - 6, a + 12):
                for d in (b, b + 6, b - 6, b + 12):
                    if math.gcd(c, d) != 1:
                        continue
                    for al in range(-30, 31):
                        for be in range(-30, 31):
                            if al * d - be * c == 1:
                                found.append(((al, be), (c, d)))
            found = sorted(found, key=lambda gm: sum(abs(x) for r in gm for x in r))
            assert found, (a, b)
            # also lifts of (a,b) mod 6 with entries shifted by 6 are equivalent; keep first two distinct gammas
            res = []
            for gam in found[:3]:
                (al, be), (c, d) = gam
                f = eta_ratio_at(gam, [3], [6], power=24 // 24)
                f = powS(f, 1)
                r = {}
                gin = ((d, -be), (-c, al))     # gamma^{-1}
                for v, vec in spins.items():
                    w = ((vec[0] * gin[0][0] + vec[1] * gin[1][0]) % 2, (vec[0] * gin[0][1] + vec[1] * gin[1][1]) % 2)
                    r[v] = S() if inv_sp[w] in ("PA", "PP") else f
                res.append(r)
            consistency.append(all((res[0][u] - r[u]).iszero(TGT) for r in res for u in r))
            Z[(a, b)] = res[0]
    out["order6_gamma_independence"] = all(consistency)
    # Delta(3 gamma tau)/Delta(6 gamma tau) is Delta^1 ratio = eta(3.)^24/eta(6.)^24 = 1/eta_g. sanity: (0,1) itself
    assert (Z[(0, 1)]["AA"] - Delta(3, 0, 1) * Delta(6, 0, 1).inv()).iszero(TGT)
    return Z, out


def project(Z, a, u):
    tot = S()
    for b in range(6):
        tot = tot + Z[(a, b)][u]
    d = {}
    for e, (x, y) in tot.d.items():
        assert y == 0 and x % 6 == 0, ("non-integral projection", a, u, e, x, y)
        d[e] = (x // 6, 0)
    return S(d)


def table(s, shift, maxw=None):
    """dims by weight h = shift + exponent, for exponents <= TGT (asserted exact)."""
    assert s.prec >= TGT, s.prec
    for e, c in s.d.items():
        assert e > TGT or c[1] == 0, ("non-rational coefficient", e, c)
    return {str(Fraction(e, U) + shift): c[0] for e, c in sorted(s.d.items()) if e <= TGT}


if __name__ == "__main__":
    Z, out = main()
    res = {"checks": out}
    res["V"] = {u: table(Z[(0, 0)][u], 3, 6) for u in ("AA", "AP", "PA")}
    sectors = {}
    for a in range(6):
        sectors[a] = {u: project(Z, a, u) for u in ("AA", "AP", "PA", "PP")}
        res[f"T_sector_{a}"] = {u: table(sectors[a][u], 3, 6) for u in ("AA", "AP", "PA", "PP")}
        res[f"twisted_raw_{a}"] = {u: table(S({e: c for e, c in Z[(a, 0)][u].d.items()}), 3, 6) for u in ("AA", "PA")} if a else None
    for u in ("AA", "AP", "PA", "PP"):
        tot = S()
        for a in range(6):
            tot = tot + sectors[a][u]
        res[f"T_{u}"] = table(tot, 3, 6)
    # nonnegativity and spin-statistics in NS: weights in Z bosonic, Z+1/2 fermionic: AP coefficient = (-1)^{2h} AA
    ok = True
    for h, n in res["T_AA"].items():
        hf = Fraction(h); ok &= n >= 0 and (2 * hf).denominator == 1
        ok &= res["T_AP"].get(h, 0) == (n if (2 * hf) % 2 == 0 else -n)
    res["checks"]["T_NS_spin_statistics"] = ok
    res["checks"]["T_R_nonneg_integral_weights"] = all(n >= 0 and Fraction(h).denominator == 1 for h, n in res["T_PA"].items())
    # the other spin-orbifold: Kitaev twist for Z_2 = <g^3>, sign (-1)^{Arf(u + abar) + Arf(u)}, abar = (a, b) mod 2
    spins = {"AA": (1, 1), "AP": (1, 0), "PA": (0, 1), "PP": (0, 0)}
    Zk = {}; kit_nonzero = []
    for (a, b), fam in Z.items():
        Zk[(a, b)] = {}
        for u, f in fam.items():
            v = spins[u]; w = ((v[0] + a) % 2, (v[1] + b) % 2)
            flip = (w == (0, 0)) != (v == (0, 0))
            Zk[(a, b)][u] = -f if flip else f
            if flip:
                kit_nonzero.append((a, b, u)) if not f.iszero(TGT) else None
    for u in ("AA", "AP", "PA", "PP"):
        tot = S()
        for a in range(6):
            tot = tot + project(Zk, a, u)
        res[f"Ttilde_{u}"] = table(tot, 3, 6)
    res["checks"]["kitaev_sign_only_on_vanishing_functions"] = not kit_nonzero
    res["Ttilde_sector_1_AA"] = table(project(Zk, 1, "AA"), 3)
    res["Ttilde_sector_3_AA"] = table(project(Zk, 3, "AA"), 3)
    res["Ttilde_sector_2_PA_PP"] = [table(project(Zk, 2, "PA"), 3), table(project(Zk, 2, "PP"), 3)]
    res["Ttilde_sector_0_PA_PP"] = [table(project(Zk, 0, "PA"), 3), table(project(Zk, 0, "PP"), 3)]
    (HERE / "out" / "fc_characters.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res["checks"]))
    for k in ("V", "T_AA", "T_PA", "T_PP", "Ttilde_AA", "Ttilde_PA", "Ttilde_PP", "Ttilde_sector_1_AA", "Ttilde_sector_3_AA", "Ttilde_sector_2_PA_PP", "Ttilde_sector_0_PA_PP"):
        print(k, res[k] if k != "V" else {u: dict(list(res[k][u].items())[:8]) for u in res[k]})
    for a in range(6):
        print("sector", a, {u: dict(list(res[f"T_sector_{a}"][u].items())[:6]) for u in ("AA", "PA", "PP")})
