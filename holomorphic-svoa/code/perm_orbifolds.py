# Copied from claude/perm_orbifolds.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Witten indices of permutation orbifolds T^{(x)n} // G of a holomorphic SCFT T.

Model (standard; see AKL arXiv:2210.14923 and the notes):

  Index(T^n // G) = (1/|G|) sum_{g h = h g}  eps(g,h) * I^{orb<g,h>}

where I = Index(T) and orb<g,h> is the number of orbits of <g,h> on the n
copies.  Each orbit is a connected covering torus; for the RR (odd) spin
structure every cover again carries the RR structure and contributes I.
With the Koszul sign convention for permuting fermionic states eps = 1.
The allowed modifications ('twists', i.e. stacking an invertible
G x Z2^F-symmetric phase before gauging) are, on an RR torus:

  * Arf on the seed:            I -> -I
  * Kitaev twist chi in Hom(G,Z2): eps *= -1 if (chi(g),chi(h)) != (0,0)
    [ (-1)^{Arf(s+chi(a)) + Arf(s)} ];  for G = S2 this equals Arf on the seed
  * discrete torsion from H^2(S_n,U(1)) restricted to G:
    eps *= commutator sign of lifts of g,h in the Pin double cover of S_n.

Anomaly: by AKL App. A the S_n action on T^n is non-anomalous for all n iff
c(T) is in 12Z, which holds for every seed considered (c = 12k), hence for
every G <= S_n.  TMF lower bound on |I| for a seed at c = 12k: 24/gcd(k,24).

Outputs results/perm_orbifolds.json and prints the solutions of
  Index = +/-1 at c = 288,  +/-2 at c = 144,  +/-4 at c = 72.
"""
import itertools
import json
from fractions import Fraction
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------- permutations
def compose(p, q):            # (p q)(i) = p(q(i))
    return tuple(p[i] for i in q)


def inverse(p):
    out = [0] * len(p)
    for i, x in enumerate(p):
        out[x] = i
    return tuple(out)


def closure(gens, n):
    e = tuple(range(n))
    elems = {e}
    frontier = [e]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = compose(g, a)
                if b not in elems:
                    elems.add(b)
                    new.append(b)
        frontier = new
    return frozenset(elems)


def cycle_type(p):
    seen, out = set(), []
    for i in range(len(p)):
        if i not in seen:
            l, j = 0, i
            while j not in seen:
                seen.add(j); j = p[j]; l += 1
            out.append(l)
    return tuple(sorted(out))


def sign(p):
    return (-1) ** (len(p) - len(cycle_type(p)))


def orbits(gens, n):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for g in gens:
        for i in range(n):
            a, b = find(i), find(g[i])
            if a != b:
                parent[a] = b
    return len({find(i) for i in range(n)})


def invariant(G, n):
    orb = orbits(list(G), n)
    return (len(G), orb, tuple(sorted(cycle_type(g) for g in G)))


def subgroups(n):
    """All subgroups of S_n up to conjugacy (n <= 6).  Every subgroup K is
    <H', g> for a maximal subgroup H' < K, so joining each class
    representative with every element reaches all classes."""
    sym = list(itertools.permutations(range(n)))
    e = tuple(range(n))
    seen = set()
    reps = []

    def register(G, gens):
        if G in seen:
            return False
        for s_ in sym:
            si = inverse(s_)
            seen.add(frozenset(compose(compose(s_, g), si) for g in G))
        reps.append((G, gens))
        return True
    triv = frozenset([e])
    register(triv, [])
    frontier = [(triv, [])]
    while frontier:
        new = []
        for H, gens in frontier:
            for g in sym:
                if g in H:
                    continue
                G = closure(gens + [g], n)
                if register(G, gens + [g]):
                    new.append((G, gens + [g]))
        frontier = new
    return [G for G, _ in reps]


# ---------------------------------------------------------------- Pin lifts
def clifford_mul(a, b):
    """Product in Cl(n) with e_i^2 = +1; elements are dicts mask->Fraction."""
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            # sign from reordering: count pairs (i in ma, j in mb, i > j)
            s = 0
            x = ma
            while x:
                i = x.bit_length() - 1
                s += bin(mb & ((1 << i) - 1)).count("1")
                x &= ~(1 << i)
            c = ca * cb * (-1 if s % 2 else 1)
            m = ma ^ mb
            out[m] = out.get(m, 0) + c
    return {m: c for m, c in out.items() if c != 0}


def pin_lift(p):
    """Lift of a permutation: product of (e_i - e_j) over a transposition
    decomposition (unnormalised; signs/scales cancel in commutators)."""
    n = len(p)
    elem = {0: Fraction(1)}
    q = list(p)
    for i in range(n):
        while q[i] != i:
            j = q[i]
            # transposition (i j) applied on the left reduces q
            elem = clifford_mul({1 << i: Fraction(1), 1 << j: Fraction(-1)}, elem)
            q[i], q[j] = q[j], q[i]
    return elem


_LIFTS = {}


def lift_commutator_sign(g, h):
    for x in (g, h):
        if x not in _LIFTS:
            _LIFTS[x] = pin_lift(x)
    a, b = _LIFTS[g], _LIFTS[h]
    ab, ba = clifford_mul(a, b), clifford_mul(b, a)
    if ab == ba:
        return 1
    neg = {m: -c for m, c in ba.items()}
    assert ab == neg, "lifts neither commute nor anticommute"
    return -1


# ---------------------------------------------------------------- homomorphisms
def hom_to_z2(G, n):
    G = list(G)
    gens = []
    H = frozenset([tuple(range(n))])
    for g in G:
        if g not in H:
            gens.append(g)
            H = closure(gens, n)
        if len(H) == len(G):
            break
    homs = []
    for vals in itertools.product((0, 1), repeat=len(gens)):
        chi = {tuple(range(n)): 0}
        frontier = [tuple(range(n))]
        ok = True
        while frontier and ok:
            new = []
            for a in frontier:
                for g, v in zip(gens, vals):
                    b = compose(g, a)
                    val = (chi[a] + v) % 2
                    if b in chi:
                        if chi[b] != val:
                            ok = False
                            break
                    else:
                        chi[b] = val
                        new.append(b)
                if not ok:
                    break
            frontier = new
        if ok:
            homs.append(chi)
    return homs


# ---------------------------------------------------------------- index polynomials
def index_polynomials(G, n):
    """{twist_name: {orbit_count: Fraction coefficient}} (seed sign handled later)."""
    G = list(G)
    pairs = [(g, h) for g in G for h in G if compose(g, h) == compose(h, g)]
    orb = {pr: orbits(list(pr), n) for pr in pairs}
    tor = {pr: lift_commutator_sign(*pr) for pr in pairs}
    homs = hom_to_z2(G, n)
    polys = {}
    for ci, chi in enumerate(homs):
        for t in (0, 1):
            poly = {}
            for pr in pairs:
                eps = 1
                if chi[pr[0]] or chi[pr[1]]:
                    eps = -eps
                if t:
                    eps *= tor[pr]
                poly[orb[pr]] = poly.get(orb[pr], 0) + Fraction(eps, len(G))
            name = f"chi{ci}" + ("+tor" if t else "")
            if t and all(tor[pr] == 1 for pr in pairs):
                continue                     # torsion restricts trivially
            polys[name] = poly
    return polys


def evaluate(poly, I):
    return sum(c * I ** k for k, c in poly.items())


def describe(G, n):
    order = len(G)
    orb = orbits(list(G), n)
    abelian = all(compose(a, b) == compose(b, a) for a in G for b in G)
    return {"order": order, "orbits": orb, "abelian": abelian,
            "generators_cycle_types": sorted({cycle_type(g) for g in G})}


def tmf_min(c):
    k = c // 12
    return 24 // gcd(k, 24)


def main():
    out = {"single_step": [], "towers": []}
    targets = {288: 1, 144: 2, 72: 4}
    library = {}
    for n in range(2, 7):
        subs = subgroups(n)
        library[n] = []
        for G in subs:
            d = describe(G, n)
            for name, poly in index_polynomials(G, n).items():
                library[n].append((d, name, poly))
        print(f"n={n}: {len(subs)} subgroup classes, {len(library[n])} (G,twist) pairs")
    # -------- single-step:  T^n // G with T at c = ct/n
    for ct, it in targets.items():
        for n in library:
            if ct % n or (ct // n) % 12:
                continue
            c0 = ct // n
            m0 = tmf_min(c0)
            for d, name, poly in library[n]:
                if d["orbits"] != 1:
                    continue                 # intransitive: index factorises
                for mult in (1, 2, 3):
                    for s in (1, -1):
                        I = s * mult * m0
                        val = evaluate(poly, I)
                        if val.denominator != 1:
                            out.setdefault("non_integral", []).append(
                                [ct, n, d, name, I, str(val)])
                        if abs(val) == it:
                            out["single_step"].append({"c_target": ct, "n": n,
                                "seed_c": c0, "seed_index": I, "G": d,
                                "twist": name, "index": int(val)})
    # -------- towers of symmetric / cyclic steps from small seeds to c=288
    def sym_coeffs(I, n):
        coeffs = [1] + [0] * n
        for m in range(1, n + 1):
            # exact binomial series of (1-p^m)^{-I}
            f = [1]
            num = Fraction(1)
            for j in range(1, n // m + 1):
                num = num * (I + j - 1) / j
                f.append(num)
            new = [0] * (n + 1)
            for i, a in enumerate(coeffs):
                for j, b in enumerate(f):
                    if i + j * m <= n:
                        new[i + j * m] += a * b
            coeffs = new
        return int(coeffs[n])

    def jordan2(d):
        r = d * d
        x, p = d, 2
        while x > 1:
            if x % p == 0:
                r = r * (p * p - 1) // (p * p)
                while x % p == 0:
                    x //= p
            p += 1
        return r

    def cyc(I, n):
        tot = sum(jordan2(d) * I ** (n // d) for d in range(1, n + 1) if n % d == 0)
        assert tot % n == 0
        return tot // n
    steps = []
    for n in (2, 3, 4, 6, 8, 12, 24):
        steps.append((f"S{n}", n, lambda I, n=n: sym_coeffs(I, n)))
        steps.append((f"Z{n}", n, lambda I, n=n: cyc(I, n)))
    # also every transitive (G,twist) from the n<=6 library
    for n in library:
        for d, name, poly in library[n]:
            if d["orbits"] == 1:
                steps.append((f"G{n}[|G|={d['order']},{name}]", n,
                              lambda I, poly=poly: evaluate(poly, I)))
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def reach(c, I):
        """Tuple of step-sequences taking (c, I) to |I| = 1 at c = 288."""
        if c == 288:
            return ((),) if abs(I) == 1 else ()
        if abs(I) > 10 ** 30:
            return ()          # beyond the growth threshold (see notes)
        found = []
        for name, n, f in steps:
            if (c * n) > 288 or 288 % (c * n):
                continue
            for s in (1, -1):
                v = f(s * I)
                if isinstance(v, Fraction):
                    if v.denominator != 1:
                        continue
                    v = int(v)
                for tail in reach(c * n, v):
                    found.append(((("Arf." if s < 0 else "") + name, c * n, v),) + tail)
                    if len(found) > 50:
                        return tuple(found)
        return tuple(found)
    sols = []
    for c0 in (12, 24, 36, 48, 72, 96, 144):
        m0 = tmf_min(c0)
        for mult in (1, 2):
            for s in (1, -1):
                for tail in reach(c0, s * mult * m0):
                    sols.append([("seed", c0, s * mult * m0)] + list(tail))
    distinct = set()
    for p in sols:
        key = tuple((c, v) for _, c, v in p)
        if key not in distinct:
            distinct.add(key)
            out["towers"].append(p)
    (HERE / "out").mkdir(exist_ok=True)
    (HERE / "out" / "perm_orbifolds.json").write_text(
        json.dumps(out, indent=1, default=str))
    print("single-step solutions:")
    for s in out["single_step"]:
        print(" ", s)
    print("distinct tower index sequences reaching |I|=1 at c=288:")
    for p in out["towers"]:
        print("  ", " -> ".join(f"{nm}:(c={c},I={v})" for nm, c, v in p))
    print("non-integral evaluations:", len(out.get("non_integral", [])))


if __name__ == "__main__":
    main()
