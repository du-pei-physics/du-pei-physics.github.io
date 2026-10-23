# Copied from claude/paper/power_operations_check.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Exact q-expansion checks for claude/paper/power_operations.md (S3).

All series are in t = q^{1/2} with exact Fraction coefficients and tracked precision.
Model (S2): f = a(q) beta^w in pi_{2w} K((q)), w = 12k, a = Z * Delta^k.
  P_2^string(f) sector characters (Ganter math/0701565 Def. 5.2, 5.9-5.10, Rem. 5.15):
    (1,1): a(q)^2            (1,s): 2^w a(q^2)   [2^w = P_2(beta^w) at the swap, i.e. psi^2(beta^w)]
    (s,1): a(q^{1/2})        (s,s): a(-q^{1/2})  [coefficient of q^{n/2} times (-1)^n]
  Stringy (Witten) Thom class of V_w = C^w (x) rho_2 built from mode products with zeta-regularised zero point
  (q^{+1/12} per complex dim for periodic modes, q^{-1/24} for antiperiodic modes; Ganter §3.2 principle).
  Sym^2_eps(f) := (1/2) sum_{(g,h)} eps(g,h) * P_2(f)(g,h) / u(V_w)(g,h) * u(C^{2w}),  u(C^{2w}) = Delta^{2k}.
Run:  .venv/bin/python claude/paper/power_operations_check.py
"""
from fractions import Fraction as F
import sys


class S:
    """Laurent series in t: coefficients dict {exp: coeff}, known for exponents < prec."""
    def __init__(self, c, prec):
        self.prec = prec
        self.c = {e: F(v) for e, v in c.items() if v != 0 and e < prec}

    def val(self):
        return min(self.c) if self.c else self.prec

    def __add__(self, o):
        p = min(self.prec, o.prec); d = dict(self.c)
        for e, v in o.c.items(): d[e] = d.get(e, 0) + v
        return S(d, p)

    def scale(self, s):
        return S({e: s * v for e, v in self.c.items()}, self.prec)

    def __sub__(self, o):
        return self + o.scale(-1)

    def __mul__(self, o):
        if not isinstance(o, S):
            return self.scale(o)
        p = min(self.prec + o.val(), o.prec + self.val()); d = {}
        for e1, v1 in self.c.items():
            for e2, v2 in o.c.items():
                e = e1 + e2
                if e < p: d[e] = d.get(e, 0) + v1 * v2
        return S(d, p)

    def inv(self):
        v = self.val(); a0 = self.c[v]; n = self.prec - v  # relative precision
        # 1/(a0 t^v (1 + r)) ; Newton-free: recursive coefficients
        b = [F(0)] * n; b[0] = 1 / a0
        for i in range(1, n):
            s = sum(self.c.get(v + j, 0) * b[i - j] for j in range(1, i + 1))
            b[i] = -s / a0
        return S({i - v: b[i] for i in range(n)}, n - v)

    def __rmul__(self, o):
        return self.scale(o)

    def __pow__(self, k):
        r = S({0: 1}, 10**9)
        for _ in range(k): r = r * self
        return r

    def subs_pow(self, m):          # t -> t^m  (q -> q^m)
        return S({e * m: v for e, v in self.c.items()}, self.prec * m - (m - 1))

    def halve(self):                # q -> q^{1/2}: requires even t-exponents
        assert all(e % 2 == 0 for e in self.c)
        return S({e // 2: v for e, v in self.c.items()}, (self.prec + 1) // 2)

    def sign_flip(self):            # t -> -t
        return S({e: (-1) ** (e % 2) * v for e, v in self.c.items()}, self.prec)

    def eq(self, o):
        p = min(self.prec, o.prec)
        return all(self.c.get(e, 0) == o.c.get(e, 0) for e in range(min(self.val(), o.val()), p))


def prod(M, factors, lead=0):
    """t^lead * prod over (e, c, m) of (1 - c t^e)^m, truncated at t^M."""
    r = S({lead: 1}, M)
    for e, c, m in factors:
        if e >= M: continue
        base = S({0: 1, e: -c}, M)
        for _ in range(m): r = r * base
    return r


def check(M_q):
    M = 2 * M_q                              # precision in t
    out = []
    Delta = prod(M, [(2 * n, 1, 24) for n in range(1, M_q + 1)], lead=2)        # q prod (1-q^n)^24
    E4 = S({2 * n: 240 * sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, M_q)} | {0: 1}, M)

    def P2(a, w):
        return {(1, 1): a * a, (1, 's'): a.subs_pow(2).scale(2 ** w),
                ('s', 1): a.halve(), ('s', 's'): a.halve().sign_flip()}

    def thom_Vw(k):
        w = 12 * k
        per = [(2 * n, 1, 2 * w) for n in range(1, M_q + 1)]                     # (1-q^n)^{2w}
        u11 = prod(M, per + per, lead=2 * (2 * w) // 12)                          # q^{2w/12}
        u1s = prod(M, per + [(2 * n, -1, 2 * w) for n in range(1, M_q + 1)], lead=2 * (2 * w) // 12).scale(2 ** w)
        # twisted: periodic part C^w: q^{w/12}; antiperiodic part C^w (x) sgn: q^{-w/24} prod_{j odd}(1 -+ q^{j/2})^{2w}
        lead_tw = 2 * w // 12 - 2 * w // 24                                      # in t-units: 2*(w/12 - w/24)
        us1 = prod(M, per + [(j, 1, 2 * w) for j in range(1, M, 2)], lead=lead_tw)
        uss = prod(M, per + [(j, -1, 2 * w) for j in range(1, M, 2)], lead=lead_tw).scale((-1) ** k)
        return {(1, 1): u11, (1, 's'): u1s, ('s', 1): us1, ('s', 's'): uss}

    def sym2(a, k, eps):
        w = 12 * k; P = P2(a, w); U = thom_Vw(k); tot = S({}, 10**9)
        for key in P:
            tot = tot + (P[key] * U[key].inv()).scale(eps[key])
        return (tot * Delta ** (2 * k)).scale(F(1, 2))

    plain = {(1, 1): 1, (1, 's'): 1, ('s', 1): 1, ('s', 's'): 1}
    arf = {(1, 1): 1, (1, 's'): -1, ('s', 1): -1, ('s', 's'): -1}
    borel = {(1, 1): 1, (1, 's'): 1, ('s', 1): 0, ('s', 's'): 0}

    # (H) H_infinity consistency: P_2^string(u(C^w)) = u(V_w) sector by sector, u(C^w) = Delta^k
    for k in (1, 2, 6):
        w = 12 * k; P = P2(Delta ** k, w); U = thom_Vw(k)
        uCw = prod(M, [(2 * n, 1, 2 * w) for n in range(1, M_q + 1)], lead=2 * w // 12)
        ok = uCw.eq(Delta ** k) and all(P[s].eq(U[s]) for s in P)
        out.append(f"(H) k={k}: u(C^w)=Delta^k and P2(u(C^w)) = u(V_w) in all four sectors: {ok}")

    # (A) index arithmetic on f = I Delta^k
    for (I, k, eps, name, target) in [(-4, 6, plain, 'plain', 2 * Delta ** 12),
                                     (2, 12, arf, 'Arf', -1 * Delta ** 24),
                                     (-2, 12, plain, 'plain', -1 * Delta ** 24),
                                     (-2, 12, arf, 'Arf', 5 * Delta ** 24),
                                     (-4, 6, borel, 'untwisted-only (Borel/K(1))', 6 * Delta ** 12)]:
        r = sym2(Delta ** k * I, k, eps)
        out.append(f"(A) Sym2_{name}({I} Delta^{k}) = {target.c[target.val()]} Delta^{2*k}: {r.eq(target)}"
                   f"  (checked to q^{r.prec//2}; leading q^{F(r.val(),2)}, coeff {r.c.get(r.val())})")

    # (W) wrong-degree variant: twisted terms without the Thom ratio give a(tau_H) of weight 12k, not 24k
    k = 6; a = Delta ** k * -4; P = P2(a, 72)
    naive = (P[(1, 1)] + P[(1, 's')].scale(F(1, 2 ** 72)) + P[('s', 1)] + P[('s', 's')]).scale(F(1, 2))
    out.append(f"(W) naive (1/2)[f^2 + sum_H a(tau_H)]: leading t-exponent {naive.val()} (q^{F(naive.val(), 2)}), "
               f"half-integral q-powers present: {any(e % 2 for e in naive.c)}")

    # (M) modularity + integrality for general f = I Delta^k + c E4^{3k}: result in M_{24k} with Z-coefficients
    for k, I, c in [(1, -24, 1), (2, 7, 3), (6, -4, 1)]:
        a = Delta ** k * I + (E4 ** (3 * k)) * c
        for name, eps in (('plain', plain), ('Arf', arf)):
            r = sym2(a, k, eps)
            # basis E4^{6k-3m} Delta^m, m = 0..2k (dim M_{24k} = 2k+1); solve triangular system on q^0..q^{2k}
            basis = [(E4 ** (6 * k - 3 * m)) * (Delta ** m) for m in range(2 * k + 1)]
            res = r; coeffs = []
            for m in range(2 * k + 1):
                cm = res.c.get(2 * m, F(0)); coeffs.append(cm); res = res - basis[m].scale(cm)
            zero = all(v == 0 for e, v in res.c.items() if e < res.prec)
            integral = all(v.denominator == 1 for v in r.c.values())
            out.append(f"(M) k={k}, f={I}Delta^{k}+{c}E4^{3*k}, {name}: in M_{24*k} (checked to q^{res.prec//2}): {zero};"
                       f" integral q-coefficients: {integral}")
    # (M') negative control: a sign pattern that is not SL(2,Z)-invariant must fail modularity
    k, I, c = 1, -24, 1
    a = Delta ** k * I + (E4 ** (3 * k)) * c
    for pat in ({(1, 1): 1, (1, 's'): -1, ('s', 1): 1, ('s', 's'): 1},
                {(1, 1): 1, (1, 's'): -1, ('s', 1): -1, ('s', 's'): 1}):
        r = sym2(a, k, pat)
        basis = [(E4 ** (6 * k - 3 * m)) * (Delta ** m) for m in range(2 * k + 1)]
        res = r
        for m in range(2 * k + 1):
            res = res - basis[m].scale(res.c.get(2 * m, F(0)))
        zero = all(v == 0 for e, v in res.c.items() if e < res.prec)
        out.append(f"(M') k=1, signs {[pat[x] for x in pat]}: in M_24: {zero}; half-integral q-powers: "
                   f"{any(e % 2 for e in r.c)}")
    return out


if __name__ == '__main__':
    M_q = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    for line in check(M_q):
        print(line)
