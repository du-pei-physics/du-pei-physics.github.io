# Copied from claude/paper/fc_closed_form.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Closed forms for the NS characters of T (c=72), T' (c=144), P (c=288).

Fact used (Gamma_theta argument, field_content.md Section 6): for a holomorphic N=1 SVOA X with 12 | c, the NS
character Z_X(tau) = Tr_NS q^{L_0 - c/24} is invariant under Gamma_theta = <S, T^2>, holomorphic on H, regular at
the cusp 1 (there it becomes the Ramond character, which starts at q^0), and has a pole of order c/12 in q^{1/2} at
infinity.  Hence Z_X = sum_{j=0}^{c/12} a_j t^j with the Hauptmodul
    t = (theta_3/eta)^12 = eta(tau)^48 / (eta(tau/2) eta(2 tau))^24 = q^{-1/2} + 24 + 276 q^{1/2} + ...
and then  Z~_X(tau) = Z_X(tau+1) = sum a_j (-1)^j t~^j,  t~ = (theta_4/eta)^12,
          Z^R_X(tau) = sum a_j (-1)^j (theta_2/eta)^{12 j},  (theta_2/eta)^12 = 2^12 eta(2tau)^24/eta(tau)^24.
In particular a_0 = Z^R_X at q^0 = number of Ramond ground states.

Series are dicts {k: c} meaning c q^{k/2}, truncated at k <= N.
"""
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent
N = 60                                     # keep q-exponents <= 30 (units of q^{1/2})


def trunc(d, n=N):
    return {k: v for k, v in d.items() if k <= n and v}


def mul(a, b, n=N):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= n:
                out[i + j] = out.get(i + j, 0) + x * y
    return trunc(out, n)


def add(*ss, coef=None):
    out = {}
    for i, s in enumerate(ss):
        c = 1 if coef is None else coef[i]
        for k, v in s.items():
            out[k] = out.get(k, 0) + c * v
    return trunc(out)


def power(a, e, n=N):
    r = {0: 1}
    for _ in range(e):
        r = mul(r, a, n)
    return r


def inv_unit(a, n=N):
    """1/a for a = 1 + O(s) with integer coefficients."""
    assert a.get(0) == 1 and min(a) == 0
    g = {0: 1}
    for k in range(1, n + 1):
        v = -sum(a.get(j, 0) * g.get(k - j, 0) for j in range(1, k + 1))
        if v:
            g[k] = v
    return g


M = N + 40                                 # inner precision (the t-powers carry negative exponents)
theta3 = trunc({n * n: (1 if n == 0 else 2) for n in range(0, 20)}, M)            # sum s^{n^2}, s = q^{1/2}
eta_q = {0: 1}                                                                    # prod (1 - q^n) in s
for n in range(1, M // 2 + 1):
    eta_q = mul(eta_q, {0: 1, 2 * n: -1}, M)
t_reg = mul(power(theta3, 12, M), inv_unit(power(eta_q, 12, M), M), M)           # s * t
eta2q = {0: 1}
for n in range(1, M // 4 + 1):
    eta2q = mul(eta2q, {0: 1, 4 * n: -1}, M)
r_reg = mul(power(eta2q, 24, M), inv_unit(power(eta_q, 24, M), M), M)            # (theta2/eta)^12 = 2^12 s^2 r_reg


def tpow(j, sign=1):
    """(sign-twisted) t^j as a Laurent series: t^j = s^{-j} t_reg^j; t(tau+1): s -> -s."""
    p = power(t_reg, j, M)
    out = {k - j: v for k, v in p.items()}
    if sign == -1:
        out = {k: v * (-1) ** (k % 2) for k, v in out.items()}
    return trunc(out)


TP = [tpow(j) for j in range(25)]


def fit(Z, deg):
    """Z = sum_{j<=deg} a_j t^j; Z given with exponents >= -deg (units q^{1/2}).  Solve top-down."""
    rest = dict(Z); a = [0] * (deg + 1)
    for j in range(deg, -1, -1):
        a[j] = rest.get(-j, 0)
        rest = add(rest, TP[j], coef=[1, -a[j]])
    return a, rest                         # rest must vanish (to precision) if Z is of the claimed form


def eval_poly(a, kind="NS"):
    tot = {}
    for j, c in enumerate(a):
        if not c:
            continue
        if kind == "NS":
            term = TP[j]
        elif kind == "AP":                 # Z(tau+1): s -> -s
            term = {k: v * (-1) ** (k % 2) for k, v in TP[j].items()}
        else:                              # Ramond: (-1)^j (2^12 q r)^j, q = s^2
            term = {k + 2 * j: v * (-(2 ** 12)) ** j for k, v in power(r_reg, j).items()}
        tot = add(tot, term, coef=[1, c])
    return tot


def s_sub(Z, factor):
    """tau -> factor * tau on a series in s = q^{1/2}: factor 2 -> s^2, factor 1/2 -> s^{1/2} (needs even keys)."""
    if factor == 2:
        return trunc({2 * k: v for k, v in Z.items()})
    assert all(k % 2 == 0 for k in Z if k <= 2 * N)
    return trunc({k // 2: v for k, v in Z.items()})


def sym2(Z, ZAP, ZR, I, arf):
    """NS and R characters of Sym^2(X) or Sym^2(X x Arf); see field_content.md Section 6.2."""
    eps = -1 if arf else 1
    ns = add(mul(Z, Z), s_sub(ZAP, 2), s_sub(ZR, 0.5), {0: eps * I})
    zr_half = s_sub(ZR, 0.5)
    zr_half_shift = {k: v * (-1) ** (k % 2) for k, v in zr_half.items()}     # Z^R((tau+1)/2)
    r = add(mul(ZR, ZR), {0: eps * I}, zr_half, zr_half_shift)
    for d in (ns, r):
        for k in d:
            assert d[k] % 2 == 0, ("odd coefficient", k, d[k])
    return {k: v // 2 for k, v in ns.items()}, {k: v // 2 for k, v in r.items()}


def weights(Z, c, upto):
    return {f"{(k + c // 12) / 2:g}": v for k, v in sorted(Z.items()) if (k + c // 12) / 2 <= upto}


if __name__ == "__main__":
    from fractions import Fraction
    ch = json.load(open(HERE / "out" / "fc_characters.json"))
    ZT_data = {int(2 * (Fraction(h) - 3)): n for h, n in ch["T_AA"].items()}      # exponents in s, weights <= 5
    res = {"t_expansion": {f"{k / 2:g}": v for k, v in sorted(TP[1].items()) if k <= 8}}
    # --- T: fit from weights 0..3 only, compare with weights 7/2..5 and with the R character
    lowT = {k: v for k, v in ZT_data.items() if k <= 0}
    aT, rest = fit(lowT, 6)
    ZT = eval_poly(aT)
    assert all(ZT.get(k, 0) == v for k, v in ZT_data.items()), "T: polynomial does not reproduce weights 7/2..5"
    ZRT = eval_poly(aT, "R"); ZAPT = eval_poly(aT, "AP")
    RT_data = {2 * (int(h) - 3): n for h, n in ch["T_PA"].items()}
    assert all(ZRT.get(k, 0) == v for k, v in RT_data.items()), "T: R character mismatch"
    assert aT[0] == 13940
    res["T"] = {"a_j": aT, "checks": "fit from h<=3 reproduces h=7/2..5 and the R character 13940 + 2579935232 q + ..."}
    # --- T' = Sym^2 T (plain, I(T) = -4): build from the closed Sym^2 formula, then fit a degree-12 polynomial
    ZTp, ZRTp = sym2(ZT, ZAPT, ZRT, -4, arf=False)
    aTp, rest = fit(ZTp, 12)
    assert not trunc(rest, N - 30), "T': Sym^2 character is not a polynomial in t of degree 12"
    assert all(eval_poly(aTp, "R").get(k, 0) == v for k, v in ZRTp.items() if k <= N - 30), "T': R mismatch"
    assert aTp[0] == 48587870 + 48587868
    res["Tprime"] = {"a_j": aTp}
    # --- P = Sym^2(T' x Arf), I(T') = 2
    ZAPTp = eval_poly(aTp, "AP"); ZRTp2 = eval_poly(aTp, "R")
    ZP, ZRP = sym2(eval_poly(aTp), ZAPTp, ZRTp2, 2, arf=True)
    aP, rest = fit(ZP, 24)
    assert not trunc(rest, N - 36), "P: not a polynomial of degree 24"
    assert all(eval_poly(aP, "R").get(k, 0) == v for k, v in ZRP.items() if k <= N - 36), "P: R mismatch"
    assert aP[0] == 2360781062549029 + 2360781062549030
    res["P"] = {"a_j": aP}
    # Faber form: Z = sum_{k} d_k F_k(t) + d_0 with d_k = dim at weight c/24 - k/2 (polar part), F_k = q^{-k/2} + O(q^{1/2})
    for name, a, c in (("T", aT, 72), ("Tprime", aTp, 144), ("P", aP, 288)):
        Z = eval_poly(a)
        res[name]["NS_dims"] = weights(Z, c, c / 24 + 2)
        res[name]["polar_part_and_constant"] = {f"{(k + c // 12) / 2:g}": v for k, v in sorted(Z.items()) if k <= 0}
        res[name]["R_character"] = {f"q^{k // 2}": v for k, v in sorted(eval_poly(a, "R").items()) if k <= 6}
    # Faber polynomials F_k (monic in t0 = t - 24, F_k = q^{-k/2} + O(q^{1/2})), k <= 6
    t0 = add(TP[1], {0: -24}); P0 = [{0: 1}] + [power(t0, j) for j in range(1, 25)]
    fab = {}; Fser = {}
    for k in range(1, 25):
        rest = P0[k]; co = {k: 1}
        for j in range(k - 1, -1, -1):
            cj = rest.get(-j, 0); co[j] = -cj
            rest = add(rest, P0[j], coef=[1, -cj])
        assert all(e > 0 or (e == -k and rest[e] == 1) for e in rest), k
        Fser[k] = rest
        fab[k] = " + ".join(f"({co[j]}) t0^{j}" for j in range(k, -1, -1) if co[j])
    res["Faber_polynomials_in_t0"] = {k: v for k, v in fab.items() if k <= 6}
    # check: Z = sum_{polar weights} dim F_k + dim at weight c/24
    for name, a, c in (("T", aT, 72), ("Tprime", aTp, 144), ("P", aP, 288)):
        Z = eval_poly(a); recon = {0: Z.get(0, 0)}
        for k in range(1, c // 12 + 1):
            if Z.get(-k, 0):
                recon = add(recon, Fser[k], coef=[1, Z[-k]])
        assert not trunc(add(recon, Z, coef=[1, -1]), N - 2 * (c // 12)), name
    res["Faber_form_checked"] = True
    # Ramond partition function (no (-1)^F): Z^R(tau) = Z~_NS(-1/tau) = P_X(t(1 - 1/tau)) = P_X(-u),
    # u = (theta2/eta)^12 = 4096 q prod(1+q^n)^24 (R character of 24 free fermions); b_j = (-1)^j a_j.
    for name, a, I in (("T", aT, -4), ("Tprime", aTp, 2), ("P", aP, -1)):
        b = [(-1) ** j * x for j, x in enumerate(a)]
        assert all(x > 0 for x in b)
        ZR = eval_poly(a, "R")
        res[name]["R_poly_in_u"] = b
        res[name]["R_even_odd_by_level"] = {f"q^{k // 2}": [(v + (I if k == 0 else 0)) // 2, (v - (I if k == 0 else 0)) // 2]
                                           for k, v in sorted(ZR.items()) if k <= 8}
    res["u_expansion"] = {f"q^{k // 2 + 1}": v * 4096 for k, v in sorted(power(r_reg, 1).items()) if k <= 10}
    (HERE / "out" / "fc_closed_form.json").write_text(json.dumps(res, indent=1))
    for k in ("T", "Tprime", "P"):
        print(k, "a_j =", res[k]["a_j"])
        print("   polar+const:", res[k]["polar_part_and_constant"])
    print("t =", res["t_expansion"])
    print("Faber:", fab)
