# Copied from claude/census/symmetry.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Monomial automorphisms of the skew Paley matrix C_12 and orbit representatives of signed
Hadamard block partitions.

Positions: 0 = inf, k = field element k mod 11 (k = 1..11; position 11 = element 0).
PSL(2,11) acts on P^1(F_11) by Mobius maps.  For each sigma we solve for signs d (rows),
e (columns) with  S[sigma i, sigma j] = d_i e_j S[i, j]  (i != j).  Then (M_sigma, N_sigma)
with M[i, sigma i] = d_i, N[j, sigma j] = e_j satisfy M S = S N, and the cell
(Q, S, R) is equivalent to (Q M, S, sigma R).  Column signs of Q M: s'_{sigma i} = d_i s_i.
A configuration (blocks, parities, R) therefore maps to
   blocks -> sigma(blocks),  parity(b) -> parity(b) + parity(d|_b),  R -> sigma(R).
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from hpcell import PALEY, evaluate

p = 11
INF = 0
def pos_of(x):            # field element / INF -> position
    return 0 if x == "inf" else (11 if x % p == 0 else x % p)
def elem_of(k):           # position -> field element or "inf"
    return "inf" if k == 0 else (k % p)

def mobius(M, k):
    a, b, c, d = M
    x = elem_of(k)
    if x == "inf":
        return pos_of("inf" if c % p == 0 else (a * pow(c, -1, p)) % p)
    den = (c * x + d) % p
    return pos_of("inf" if den == 0 else ((a * x + b) * pow(den, -1, p)) % p)

def group_elements():
    """PSL(2,11) as position permutations (660 distinct), each with its sign vector d."""
    seen = {}
    for a, b, c, d in itertools.product(range(p), repeat=4):
        if (a * d - b * c) % p != 1:
            continue
        perm = tuple(mobius((a, b, c, d), k) for k in range(12))
        if perm in seen:
            continue
        sg = signs_for(perm)
        assert sg is not None, ("no monomial lift", (a, b, c, d))
        seen[perm] = sg
    return seen

def signs_for(perm):
    """Find d, e in {+-1}^12 with S[perm i, perm j] = d_i e_j S[i,j] (as +-1 for i != j)."""
    S = np.where(PALEY == 2, -1, PALEY)
    d = np.zeros(12, dtype=int); e = np.zeros(12, dtype=int)
    d[0] = 1
    for j in range(1, 12):
        e[j] = S[perm[0], perm[j]] * S[0, j]           # e_j = S'/S with d_0 = 1
    for i in range(1, 12):
        d[i] = S[perm[i], perm[1 if i != 1 else 2]] * S[i, 1 if i != 1 else 2] * e[1 if i != 1 else 2]
    # e_0 from any i != 0
    e[0] = S[perm[1], perm[0]] * S[1, 0] * d[1]
    for i in range(12):
        for j in range(12):
            if i != j and S[perm[i], perm[j]] != d[i] * e[j] * S[i, j]:
                return None
    return d.tolist(), e.tolist()

def act(perm, d, blocks, parities, R):
    nb = []
    for b, par in zip(blocks, parities):
        nb.append((tuple(sorted(perm[k] for k in b)), (par + sum(1 for k in b if d[k] == -1)) % 2))
    nb.sort()
    return tuple(nb), tuple(sorted(perm[k] for k in R))

def canon(blocks, parities):
    nb = sorted((tuple(sorted(b)), par % 2) for b, par in zip(blocks, parities))
    return tuple(nb)

def all_configs():
    pts = list(range(12))
    out = []
    for b1 in itertools.combinations(pts, 4):
        if b1[0] != 0:
            continue
        rest = [k for k in pts if k not in b1]
        for b2 in itertools.combinations(rest, 4):
            if b2[0] != rest[0]:
                continue
            b3 = tuple(k for k in rest if k not in b2)
            for pars in itertools.product((0, 1), repeat=3):
                out.append(canon((b1, b2, b3), pars))
    return out

def orbit_reps(G, configs):
    seen = set(); reps = []
    gens = list(G.items())
    for c in configs:
        if c in seen:
            continue
        orbit = {c}; frontier = [c]
        while frontier:
            x = frontier.pop()
            blocks = [b for b, _ in x]; pars = [q for _, q in x]
            for perm, (d, e) in gens:
                y, _ = act(perm, d, blocks, pars, ())
                if y not in orbit:
                    orbit.add(y); frontier.append(y)
        seen |= orbit
        stab = [perm for perm, (d, e) in gens if act(perm, d, [b for b, _ in c], [q for _, q in c], ())[0] == c]
        reps.append({"config": [[list(b), q] for b, q in c], "orbit_size": len(orbit), "stabilizer_order": len(stab)})
    return reps

if __name__ == "__main__":
    G = group_elements()
    print("group elements with monomial lift:", len(G))
    configs = all_configs()
    print("signed partitions:", len(configs))
    reps = orbit_reps(G, configs)
    print("orbits:", len(reps), "sum of orbit sizes:", sum(r["orbit_size"] for r in reps))
    from collections import Counter
    print("stabilizer orders:", Counter(r["stabilizer_order"] for r in reps))
    (Path(__file__).resolve().parent.parent / "out" / "census" / "reps.json").write_text(json.dumps({"group_order": len(G), "n_configs": len(configs), "reps": reps}, indent=1))

    # empirical check: invariants constant along orbits (random configs, random group elements, random R)
    rng = np.random.default_rng(1); items = list(G.items()); bad = 0
    for t in range(40):
        c = configs[rng.integers(len(configs))]
        blocks = [list(b) for b, _ in c]; pars = [q for _, q in c]
        r = int(rng.choice([2, 3, 4, 5])); R = sorted(rng.choice(12, r, replace=False).tolist())
        perm, (d, e) = items[rng.integers(len(items))]
        c2, R2 = act(perm, d, blocks, pars, R)
        o1, T1, Ta1, n1 = evaluate(blocks, pars, R)
        o2, T2, Ta2, n2 = evaluate([list(b) for b, _ in c2], [q for _, q in c2], list(R2))
        inv1 = (abs(o1["I1"]), abs(o1["I2"]), abs(o1["I3"]), sorted(np.round(np.abs(np.linalg.eigvals(T1.astype(float))), 5)))
        inv2 = (abs(o2["I1"]), abs(o2["I2"]), abs(o2["I3"]), sorted(np.round(np.abs(np.linalg.eigvals(T2.astype(float))), 5)))
        if inv1 != inv2:
            bad += 1; print("MISMATCH", c, R, perm, inv1[:3], inv2[:3])
    print("orbit-invariance mismatches:", bad, "/ 40")
    # control: a sign flip that is NOT in the group must change things sometimes
    diff = 0
    for t in range(20):
        c = configs[rng.integers(len(configs))]
        blocks = [list(b) for b, _ in c]; pars = [q for _, q in c]
        R = sorted(rng.choice(12, 3, replace=False).tolist())
        o1, T1, _, _ = evaluate(blocks, pars, R)
        o2, T2, _, _ = evaluate(blocks, [(pars[0] + 1) % 2] + pars[1:], R)
        if (o1["I1"], o1["I2"], o1["I3"]) != (o2["I1"], o2["I2"], o2["I3"]):
            diff += 1
    print("parity flip changes (I1,I2,I3) in", diff, "/ 20 cases (so parities are genuine data)")
