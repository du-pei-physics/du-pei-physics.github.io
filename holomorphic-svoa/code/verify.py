#!/usr/bin/env python3
"""Independent verifier for the certificates of the holomorphic-SVOA note.

    python3 verify.py                      # every certificates/k*.json next to this folder
    python3 verify.py ../certificates/k06_claude_canonical.json ...
    python3 verify.py --selftest           # band DP versus brute force on random small chains
    python3 verify.py --json out.json      # also write the full reports

Only numpy is needed.  Everything is exact: arithmetic over F_3 in int64, signed counts as exact
integers (Chinese remaindering over several primes for the transfer products).

For each certificate (a generator matrix in systematic form, called R in the paper and G in this code, plus a
signed-permutation action):
  0. the stored rowspace hash (SHA-256 of the RREF over F_3) is recomputed; distinct hashes = distinct codes;
  1. self-duality over F_3: G G^T = 0 and rank G = n/2;
  2. the code index  Index(C) = sum over codewords of full support of (-1)^{#coordinates equal to 2},
     and the number of such codewords (= Ramond ground states of the lattice SVOA, even + odd),
     computed from the generator alone by a banded transfer DP (method below);
  3. the symmetry: each generator preserves C (x in C iff x G^T = 0, as C is self-dual), its order and
     cycle type; for a cyclic group Z_m (m = gcd(k,24)) the hypotheses of Gaiotto-Johnson-Freyd
     Thm 5.6 as used here: g^(m/2) = -1 if m is even; if 3 | m, h = g^(m/3) has n/3 three-cycles and
     the codewords of full support fixed by h are balanced (signed count 0), computed by contracting
     the fixed subcode to one coordinate per 3-cycle; for Q8: the quaternion relations, |group| = 8,
     every commuting pair generates a cyclic subgroup and every element != 1 has -1 as a power;
  4. orbifold index = Index(C) / |G| and comparison with the TMF minimum 24/gcd(24,k).
The report keeps three things apart: `finite_checks_pass` (items 0-3: exact finite facts about the code and the
symmetry), `orbifold_index_status` (what the orbifold index rests on beyond them: GJF Thm 5.6 for cyclic groups, cited,
with its finite hypotheses checked; for Q8 "orbifold_index_conditional_on_anomaly_lemma", the non-anomaly argument of the
note, which no finite check can replace), and the derived rows (cited constructions + index arithmetic).
Finally the rows without a new code are recomputed by index arithmetic, and the twisted D4 orbifold of four copies of
the k = 6 theory is evaluated independently as a sum over the 40 commuting pairs of D4 (sector by sector).

Method for 2.  Write G = (I | A).  A codeword of full support is x = (a | aA) with a in {1,2}^{n/2}.
Split the message coordinates into consecutive blocks R_1..R_L and the parity coordinates into
blocks C_1..C_L of the same sizes ("band widths", a hint stored in the certificate).  The DP is
valid as soon as A[R_i, C_j] = 0 unless j in {i, i+1 mod L}; this zero pattern is CHECKED on the
generator, so a wrong hint can only make the check fail, never give a wrong number.  Then
    Index = Tr(T_1 T_2 ... T_L),
    T_i[u, v] = sum over a in {1,2}^{R_i} with a A[R_i,C_{i+1}] = v of phi(a) sigma(u + a A[R_i,C_i]),
with phi(a) = (-1)^{#2 in a}, sigma(b) = (-1)^{#2 in b} if b has no zero entry and 0 otherwise, and u
running over the vectors handed over from block i-1 (the boundary spins of the chain).  Nothing about
cells, Q, S or boundary coordinates is used.  Codes of dimension <= 24 (the contracted fixed codes at k = 6, 9) are instead enumerated by
brute force over an information set.
"""
import argparse, glob, hashlib, itertools, json, os, sys, time
from math import gcd

import numpy as np

# ----------------------------------------------------------------------------------- F_3 linear algebra

def rref(M):
    """Reduced row echelon form over F_3. Returns (R, pivots) with zero rows removed."""
    R = np.array(M, dtype=np.int64) % 3
    rows, cols = R.shape
    piv, r = [], 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(R[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            R[[r, p]] = R[[p, r]]
        if R[r, c] == 2:
            R[r] = (2 * R[r]) % 3
        others = np.nonzero(R[:, c])[0]
        others = others[others != r]
        if len(others):
            R[others] = (R[others] - np.outer(R[others, c], R[r])) % 3
        piv.append(c)
        r += 1
    return R[:r], piv


def rank3(M):
    return len(rref(M)[1])


def left_nullspace(K):
    """Basis (rows) of {a : a K = 0} over F_3."""
    Kt = np.array(K, dtype=np.int64).T % 3          # solve Kt a^T = 0
    R, piv = rref(Kt)
    n = Kt.shape[1]
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64)
        v[f] = 1
        for i, p in enumerate(piv):
            v[p] = (-R[i, f]) % 3
        basis.append(v)
    return np.array(basis, dtype=np.int64).reshape(len(basis), n)

# ----------------------------------------------------------------------------------- signed permutations
# A signed permutation (perm, signs) maps x to y with y[perm[i]] = signs[i] x[i].

def sp_apply(rows, perm, signs):
    rows = np.asarray(rows, dtype=np.int64)
    y = np.zeros_like(rows)
    y[:, perm] = (rows * np.asarray(signs, dtype=np.int64)) % 3
    return y


def sp_compose(a, b):
    """first a, then b."""
    pa, sa = a
    pb, sb = b
    perm = pb[pa]
    signs = sa * sb[pa]
    return perm, signs


def sp_identity(n):
    return np.arange(n), np.ones(n, dtype=np.int64)


def sp_minus(n):
    return np.arange(n), -np.ones(n, dtype=np.int64)


def sp_eq(a, b):
    return np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def sp_power(a, e):
    out = sp_identity(len(a[0]))
    for _ in range(e):
        out = sp_compose(out, a)
    return out


def sp_order(a, cap=10000):
    n = len(a[0]); cur = a; o = 1
    while not sp_eq(cur, sp_identity(n)):
        cur = sp_compose(cur, a); o += 1
        if o > cap:
            return None
    return o


def sp_cycles(a):
    perm = a[0]; n = len(perm); seen = np.zeros(n, bool); cyc = []
    for i in range(n):
        if not seen[i]:
            c = []; j = i
            while not seen[j]:
                seen[j] = True; c.append(j); j = perm[j]
            cyc.append(c)
    return cyc


def cycle_type(a):
    from collections import Counter
    cnt = Counter(len(c) for c in sp_cycles(a))
    return " ".join(f"{l}^{m}" for l, m in sorted(cnt.items()))


def sp_key(a):
    return (tuple(int(x) for x in a[0]), tuple(int(x) for x in a[1]))

# ----------------------------------------------------------------------------------- exact integer matrix products (CRT)

def _primes_below(start, count):
    out, p = [], start
    while len(out) < count:
        p -= 1
        if p % 2 and all(p % q for q in range(3, int(p ** 0.5) + 1, 2)):
            out.append(p)
    return out


PRIMES = _primes_below(1 << 26, 12)


def trace_of_product_exact(mats, bound_bits):
    """Tr(M_1 ... M_L) for integer matrices given as int64 arrays, exactly, by CRT."""
    need, prod_p, used = bound_bits + 2, 1, []
    for p in PRIMES:
        used.append(p); prod_p *= p
        if prod_p.bit_length() > need:
            break
    else:
        raise ValueError("bound too large for the prime table")
    residues = []
    for M, Mn in zip(mats, mats[1:]):
        assert M.dtype == np.int64 and Mn.dtype == np.int64 and M.shape[1] == Mn.shape[0], "incompatible factors"
    for p in used:
        # int64 safety: each entry of a product is a sum of `inner` terms < (p-1)^2, so inner*(p-1)^2 must stay < 2^63
        for M in mats[1:]:
            assert M.shape[0] * (p - 1) ** 2 < 2 ** 63, f"int64 overflow risk: inner dimension {M.shape[0]}, prime {p}"
        P = mats[0] % p
        for M in mats[1:]:
            P = (P @ (M % p)) % p
        residues.append(int(np.trace(P)) % p)
    x, mod = 0, 1                                   # CRT, then symmetric lift
    for r, p in zip(residues, used):
        t = ((r - x) * pow(mod, -1, p)) % p
        x += mod * t; mod *= p
    return x - mod if x > mod // 2 else x

# ----------------------------------------------------------------------------------- index computations

SIGN = np.array([0, 1, -1], dtype=np.int64)         # weight of a coordinate value 0, 1, 2


def all_pm(w):
    """all a in {1,2}^w as an (2^w, w) int8 array"""
    idx = np.arange(1 << w, dtype=np.int64)[:, None]
    return (1 + ((idx >> np.arange(w)) & 1)).astype(np.int8)


def band_partition_ok(A, widths):
    off = np.concatenate([[0], np.cumsum(widths)])
    L = len(widths)
    if off[-1] != A.shape[0] or A.shape[0] != A.shape[1] or L < 3:
        return False
    for i in range(L):
        for j in range(L):
            if j not in (i, (i + 1) % L) and A[off[i]:off[i + 1], off[j]:off[j + 1]].any():
                return False
    return True


def band_index(G, widths):
    """(Index, number of full-support codewords) by the banded transfer DP; None if not applicable."""
    G = np.asarray(G, dtype=np.int64) % 3
    d, n = G.shape
    if n != 2 * d or not np.array_equal(G[:, :d], np.eye(d, dtype=np.int64)):
        return None
    A = G[:, d:]
    if not band_partition_ok(A, widths):
        return None
    L = len(widths)
    # guards: block widths are enumerated (2^w words; float64 sums of +-1 are exact below 2^53) and boundary vectors
    # are encoded as base-3 integers in int64 (3^w < 2^63 for w <= 39)
    assert max(widths) <= 26, f"block of width {max(widths)}: 2^w enumeration too large for this verifier"
    off = np.concatenate([[0], np.cumsum(widths)])
    blk = lambda i, j: A[off[i]:off[i + 1], off[j]:off[j + 1]]
    ins = [all_pm(w).astype(np.int64) for w in widths]
    phis = [SIGN[a].prod(axis=1) for a in ins]
    # states handed from block i-1 into block i: the distinct values of a_{i-1} A[R_{i-1}, C_i]
    out_codes, states = [], []
    for i in range(L):
        j = (i + 1) % L
        v = (ins[i] @ blk(i, j)) % 3
        code = v @ (3 ** np.arange(v.shape[1], dtype=np.int64))
        uniq, inv = np.unique(code, return_inverse=True)
        out_codes.append(inv)
        reps = np.zeros((len(uniq), v.shape[1]), dtype=np.int64)
        reps[inv] = v
        states.append(reps)                          # states[i] = vectors entering block i+1
    Ts, Us = [], []
    for i in range(L):
        incoming = states[(i - 1) % L]
        base = (ins[i] @ blk(i, i)) % 3
        nout = len(states[i])
        T = np.zeros((len(incoming), nout), dtype=np.int64)
        U = np.zeros((len(incoming), nout), dtype=np.int64)
        for si, s in enumerate(incoming):
            vals = (base + s) % 3
            sig = SIGN[vals].prod(axis=1)
            T[si] = np.rint(np.bincount(out_codes[i], weights=(sig * phis[i]).astype(np.float64), minlength=nout)).astype(np.int64)
            U[si] = np.rint(np.bincount(out_codes[i], weights=(sig != 0).astype(np.float64), minlength=nout)).astype(np.int64)
        Ts.append(T); Us.append(U)
    bound = int(np.log2(max(len(s) for s in states))) + 1 + int(sum(widths))
    return trace_of_product_exact(Ts, bound), trace_of_product_exact(Us, bound)


def brute_index(G, max_dim=24):
    """(Index, number of full-support words) by enumerating an information set."""
    R, piv = rref(G)
    d = len(piv)
    if d > max_dim:
        return None
    if d == 0:
        return 0, 0
    idx, tot = 0, 0
    chunk = 1 << min(d, 16)
    hi_bits = d - min(d, 16)
    lo = all_pm(min(d, 16)).astype(np.int64)
    lo_part = (lo @ R[:min(d, 16)]) % 3
    for h in range(1 << hi_bits):
        hv = np.array([1 + ((h >> b) & 1) for b in range(hi_bits)], dtype=np.int64)
        shift = (hv @ R[min(d, 16):]) % 3 if hi_bits else 0
        X = (lo_part + shift) % 3
        sig = SIGN[X].prod(axis=1)
        idx += int(sig.sum()); tot += int((sig != 0).sum())
    return idx, tot


def code_index(G, widths):
    r = band_index(G, widths)
    if r is not None:
        return r + ("banded transfer DP",)
    r = brute_index(G)
    if r is not None:
        return r + ("brute force over an information set",)
    return None

# ----------------------------------------------------------------------------------- symmetry checks

def preserves(G, a):
    return not (sp_apply(G, *a) @ G.T % 3).any()


def fixed_subcode(G, a):
    n = G.shape[1]
    K = (sp_apply(G, *a) - G) % 3                  # a-th row: g(x) - x for x = row
    N = left_nullspace(K)
    if N.shape[0] == 0:
        return np.zeros((0, n), dtype=np.int64)
    return rref((N @ G) % 3)[0]


def balance_of_fixed_words(G, h, period_widths):
    """Signed count (and number) of full-support codewords fixed by h, via contraction to one
    coordinate per cycle.  Returns dict."""
    n = G.shape[1]
    F = fixed_subcode(G, h)
    rep = {"fixed_dim": int(F.shape[0])}
    cyc = sp_cycles(h)
    perm, signs = h
    eps_bits, reps, dead = 0, [], False
    for c in cyc:
        c0 = min(c)
        order = [c0]
        while perm[order[-1]] != c0:
            order.append(perm[order[-1]])
        # x_{perm[i]} = signs[i] x_i along the cycle; nonzero fixed values need product of signs = 1
        if len(order) % 2 == 0:
            rep["note"] = "cycle of even length: sign bookkeeping not implemented"
            return rep
        cum = 1; minus = 0
        for t, i in enumerate(order):
            if t > 0:
                minus += cum == -1
            cum *= signs[i]
        if cum != 1:
            dead = True
        eps_bits += minus
        reps.append(c0)
    if dead:
        rep.update({"signed_count": 0, "count": 0, "note": "some cycle has sign product -1: no fixed full-support word"})
        return rep
    reps = sorted(reps)
    Fc = rref(F[:, reps])[0]
    rep["contracted_length"] = len(reps)
    rep["contracted_dim"] = int(Fc.shape[0])
    rep["contracted_self_dual"] = bool(2 * Fc.shape[0] == len(reps) and not ((Fc @ Fc.T) % 3).any())
    res = brute_index(Fc)
    method = "brute force"
    if res is None:
        res = band_index(Fc, period_widths) if period_widths else None
        method = "banded transfer DP"
    if res is None:
        rep["note"] = "contracted code neither small nor banded: balance not computed"
        return rep
    eps = -1 if eps_bits % 2 else 1
    rep.update({"signed_count": eps * res[0], "count": res[1],
                "even": (res[1] + eps * res[0]) // 2, "odd": (res[1] - eps * res[0]) // 2, "method": method})
    return rep


def group_closure(gens, n):
    elems = {sp_key(sp_identity(n)): sp_identity(n)}
    frontier = list(elems.values())
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = sp_compose(x, g)
                k = sp_key(y)
                if k not in elems:
                    elems[k] = y; new.append(y)
        frontier = new
        if len(elems) > 1000:
            break
    return list(elems.values())


def check_q8(G, gens):
    n = G.shape[1]
    i, j = gens["i"], gens["j"]
    k = sp_compose(i, j)
    m1 = sp_minus(n)
    out = {"i,j preserve C": preserves(G, i) and preserves(G, j),
           "i^2=j^2=k^2=-1": all(sp_eq(sp_compose(x, x), m1) for x in (i, j, k)),
           "ji=-k": sp_eq(sp_compose(j, i), sp_compose(k, m1))}
    grp = group_closure([i, j], n)
    out["order"] = len(grp)
    ident = sp_identity(n)
    cyclic_pairs = True
    for a in grp:
        for b in grp:
            if sp_eq(sp_compose(a, b), sp_compose(b, a)):
                sub = group_closure([a, b], n)
                if not any(sp_order(x) == len(sub) for x in sub):
                    cyclic_pairs = False
    out["commuting pairs cyclic"] = cyclic_pairs
    out["-1 is a power of every element != 1"] = all(
        sp_eq(x, ident) or any(sp_eq(sp_power(x, e), m1) for e in range(1, 5)) for x in grp)
    out["i cycle type"] = cycle_type(i)
    out["ok"] = all(v for kk, v in out.items() if kk not in ("order", "i cycle type")) and out["order"] == 8
    return out


def check_cyclic(G, gen, m, period_widths):
    n = G.shape[1]
    out = {"g preserves C": preserves(G, gen)}
    o = sp_order(gen)
    out["order"] = o
    out["cycle type"] = cycle_type(gen)
    ok = out["g preserves C"] and o == m
    if m % 2 == 0:
        out[f"g^{m // 2} = -1"] = sp_eq(sp_power(gen, m // 2), sp_minus(n))
        ok = ok and out[f"g^{m // 2} = -1"]
    if m % 3 == 0:
        h = sp_power(gen, m // 3)
        ct = cycle_type(h)
        out["h = g^(m/3) cycle type"] = ct
        three = ct == f"3^{n // 3}"
        bal = balance_of_fixed_words(G, h, period_widths)
        out["fixed words of h"] = bal
        ok = ok and three and bal.get("signed_count") == 0
    out["ok"] = bool(ok)
    return out

# ----------------------------------------------------------------------------------- driver

def verify_certificate(path):
    t0 = time.process_time()
    c = json.load(open(path))
    G = np.array([[int(ch) for ch in row] for row in c["generator"]], dtype=np.int64)
    d, n = G.shape
    k = c["k"]
    rep = {"file": os.path.basename(path), "k": k, "length": n, "found_by": c["found_by"]["team"].split()[0]}
    rep["length = 12k"] = n == 12 * k
    rep["rowspace_sha256"] = rowspace_sha256(G)
    rep["rowspace hash = stored"] = rep["rowspace_sha256"] == c.get("rowspace_sha256")
    rep["self_dual"] = bool(2 * d == n and rank3(G) == d and not ((G @ G.T) % 3).any())
    ci = code_index(G, c.get("band_widths", []))
    if ci is None:
        rep["index"] = None; rep["ok"] = False
        return rep
    rep["index"], rep["full_support_words"], rep["index_method"] = ci
    rep["even_odd"] = [(ci[1] + ci[0]) // 2, (ci[1] - ci[0]) // 2]
    rep["index = claimed"] = rep["index"] == c["claimed_code_index"]
    rep["|index| = 24"] = abs(rep["index"]) == 24
    sym = c["symmetry"]
    m = gcd(24, k)
    gens = {name: (np.array(g["permutation"], dtype=np.int64), np.array(g["signs"], dtype=np.int64))
            for name, g in sym.get("generators", {}).items()}
    nb = len(c.get("band_widths", []))
    period_widths = c["band_widths"][: nb // 3] if nb % 3 == 0 and nb else None
    if sym["group"] == "Q8":
        s = check_q8(G, gens)
        order = 8
        rep["symmetry"] = "Q8, i: " + s["i cycle type"]
    elif m == 1:
        s = {"ok": True}
        order = 1
        rep["symmetry"] = "none needed (m = 1)"
    else:
        s = check_cyclic(G, gens["g"], m, period_widths)
        order = m
        desc = f"Z{m}, g: {s['cycle type']}"
        if m % 3 == 0:
            fw = s["fixed words of h"]
            desc += f"; h fixed words {fw.get('even')}/{fw.get('odd')}"
        rep["symmetry"] = desc
    rep["symmetry_checks"] = s
    rep["orbifold_order"] = order
    rep["orbifold_index"] = rep["index"] // order if rep["index"] % order == 0 else None
    rep["tmf_min"] = 24 // m
    rep["orbifold_index = claimed"] = rep["orbifold_index"] == c["claimed_orbifold_index"]
    # what the machine checks prove: finite facts about the code and the symmetry (exact computations) ...
    rep["finite_checks_pass"] = bool(rep["length = 12k"] and rep["rowspace hash = stored"] and rep["self_dual"]
                                     and rep["index = claimed"] and rep["|index| = 24"] and s["ok"])
    # ... and on what the orbifold index rests beyond them
    if order == 1:
        rep["orbifold_index_status"] = "no orbifold (m = 1): the index of the lattice SCFT itself"
    elif sym["group"] == "Q8":
        rep["orbifold_index_status"] = "orbifold_index_conditional_on_anomaly_lemma (Q8 non-anomalous: details F.3)"
    else:
        rep["orbifold_index_status"] = "GJF Thm 5.6 (cited), its finite hypotheses checked here"
    rep["ok"] = bool(rep["finite_checks_pass"] and rep["orbifold_index = claimed"]
                     and abs(rep["orbifold_index"]) == rep["tmf_min"])
    rep["cpu_seconds"] = round(time.process_time() - t0, 1)
    return rep


def rowspace_sha256(G):
    """SHA-256 of the RREF of G over F_3, rows written over {0,1,2} and joined by newlines."""
    R, _ = rref(G)
    return hashlib.sha256("\n".join("".join(str(int(v)) for v in row) for row in R).encode()).hexdigest()


def d4_sectors(I, arf=True):
    """Twisted D4 orbifold of four copies of a theory of index I: (1/8) sum over commuting pairs (g,h) of
    eps(g,h) I^{#orbits of <g,h>}, split by the conjugacy class of g (the twisted sector).  D4 = Z2 wr Z2 acts on the
    corners 0..3 with inner swaps a = (01), b = (23) and outer swap c = (02)(13); with arf=True,
    eps(g,h) = -1 exactly when (g,h) has nontrivial image in D4/<a,b> (the paper's twist), else eps = 1."""
    from fractions import Fraction
    def comp(p, q):
        return tuple(q[p[i]] for i in range(4))
    a, b, c = (1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)
    G = {(0, 1, 2, 3)}
    frontier = list(G)
    while frontier:
        new = []
        for x in frontier:
            for y in (a, b, c):
                z = comp(x, y)
                if z not in G:
                    G.add(z); new.append(z)
        frontier = new
    K = {(0, 1, 2, 3), a, b, comp(a, b)}

    def orbits(gens):
        seen, n = set(), 0
        for s in range(4):
            if s not in seen:
                n += 1; stack = [s]
                while stack:
                    x = stack.pop()
                    if x not in seen:
                        seen.add(x); stack.extend(g[x] for g in gens)
        return n

    def cls(g):
        inv = {x: tuple(sorted(range(4), key=lambda i: x[i])) for x in G}
        return frozenset(comp(comp(inv[x], g), x) for x in G)

    names = {"untwisted": (0, 1, 2, 3), "inner swaps": a, "both inner swaps": comp(a, b), "outer swaps": c,
             "rotations": comp(a, c)}
    out = {}
    for name, rep in names.items():
        C = cls(rep)
        tot = Fraction(0)
        for g in C:
            for h in G:
                if comp(g, h) == comp(h, g):
                    eps = -1 if arf and (g not in K or h not in K) else 1
                    tot += eps * Fraction(I) ** orbits([g, h])
        out[name] = tot / 8
    assert len(G) == 8 and sum(len(cls(r)) for r in names.values()) == 8
    return out


def derived_rows(results):
    """Index arithmetic for the rows obtained without a new code (symmetric squares, tensor products)."""
    best = {}
    for r in results:
        if r["ok"]:
            best.setdefault(r["k"], r["orbifold_index"])
    lines = []
    I6 = best.get(6)
    if I6 is not None:
        I12 = I6 * (I6 + 3) // 2                     # Sym^2 T, untwisted swap
        I24 = (-I12) * ((-I12) + 3) // 2              # Sym^2 (T' x Arf)
        lines.append((12, I12, f"T' = Sym^2 T6: I(I+3)/2 at I = {I6}"))
        lines.append((24, I24, f"T'' = Sym^2(T' x Arf): I(I+3)/2 at I = {-I12}  [= twisted D4 orbifold of T6^4]"))
        best[12], best[24] = I12, I24
        sec = d4_sectors(I6)
        plain = d4_sectors(I6, arf=False)
        poly = (I6 ** 4 + 6 * I6 ** 3 + 3 * I6 ** 2 - 18 * I6) // 8
        d4 = sum(sec.values())
        lines.append((24, int(d4), "D4 check at I(T6) = %d: sectors %s, total %s; without the Arf phase %s, total %s; "
                      "(I^4+6I^3+3I^2-18I)/8 = %d; agrees with the Sym^2 chain: %s"
                      % (I6, ", ".join(f"{k} {v}" for k, v in sec.items()), d4,
                         ", ".join(str(v) for v in plain.values()), sum(plain.values()), poly, d4 == I24 == poly)))
    for k, (a, b) in {17: (8, 9), 19: (16, 3), 20: (8, 12), 22: (16, 6), 23: (8, 15)}.items():
        ia, ib = best.get(a), best.get(b)
        if b == 3:
            ib = 8                                    # GJF Example 5.7: c = 36, index 8
        if ia is not None and ib is not None:
            lines.append((k, ia * ib, f"T{a} (x) T{b}: {ia} * {ib}"))
    out = []
    for k, idx, how in sorted(lines):
        tm = 24 // gcd(24, k)
        out.append(f"{k:>3} | {12 * k:>4} | cited+arithmetic | {idx:>4} | {tm:>3} | {'matches' if abs(idx) == tm else 'FAIL'} | {how}")
    return out


def selftest(trials=20, seed=1):
    """Random small self-dual chains (cells of width 2, 4 and 6): band DP == brute force."""
    rng = np.random.default_rng(seed)
    H4 = np.array([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 2, 2], [1, 2, 2, 1]])   # H4 H4^T = 4 = 1
    S2 = np.array([[1, 1], [1, 2]])                                          # S2 S2^T = 2 = -1
    chi5 = {0: 0, 1: 1, 4: 1, 2: 2, 3: 2}                                    # Legendre symbol mod 5 (2 = -1)
    C6 = np.zeros((6, 6), dtype=np.int64)                                    # Paley conference matrix, C6 C6^T = 5 = -1
    C6[0, 1:] = 1; C6[1:, 0] = 1
    for x in range(5):
        for y in range(5):
            C6[1 + x, 1 + y] = chi5[(y - x) % 5]

    def sperm(w):
        P = np.zeros((w, w), dtype=np.int64)
        P[np.arange(w), rng.permutation(w)] = rng.choice([1, 2], size=w)
        return P

    def orth(w):
        Q = sperm(w)
        if w < 4:
            return Q
        B = np.eye(w, dtype=np.int64)
        B[:4, :4] = H4
        for _ in range(3):
            Q = (Q @ B @ sperm(w)) % 3
        return Q

    def anti(w):
        if w == 6:
            base = C6
        else:
            base = np.zeros((w, w), dtype=np.int64)
            for b in range(0, w, 2):
                base[b:b + 2, b:b + 2] = S2
        return (orth(w) @ base @ orth(w)) % 3

    bad = 0
    for t in range(trials):
        # total message length 12 or 24 (code length 24 or 48): other lengths have index 0
        choices = [[4, 4, 4], [2, 4, 6], [6, 2, 4], [2, 2, 4, 4], [6, 6, 6, 6], [4, 6, 4, 6, 4], [6, 6, 4, 4, 4]]
        widths = choices[t % len(choices)] if t < 2 * len(choices) else choices[int(rng.integers(0, 4))]
        r = int(rng.integers(1, min(min(widths), 4)))
        cells = [(orth(w), anti(w)) for w in widths]
        L = len(widths); d = sum(widths)
        off = np.concatenate([[0], np.cumsum(widths)])
        A = np.zeros((d, d), dtype=np.int64)
        wrap = 1 if rng.random() < 0.5 else 2
        for j, (Q, S) in enumerate(cells):
            Sn = cells[(j + 1) % L][1]
            A[off[j]:off[j + 1], off[j]:off[j + 1]] = (Q[:, r:] @ S[r:, :]) % 3
            nxt = (Q[:, :r] @ Sn[:r, :]) % 3
            if j == L - 1:
                nxt = (wrap * nxt) % 3
            A[off[j]:off[j + 1], off[(j + 1) % L]:off[(j + 1) % L + 1]] = nxt
        G = np.hstack([np.eye(d, dtype=np.int64), A]) % 3
        sd = not ((G @ G.T) % 3).any()
        b1, b2 = band_index(G, widths), brute_index(G)
        ok = sd and b1 == b2
        bad += not ok
        print(f"selftest {t:2d}: widths {widths} r={r} {'nega' if wrap == 2 else 'cyc '} self-dual {sd}  "
              f"band (index, words) {b1}  brute {b2}  {'ok' if ok else 'MISMATCH'}")
    print("selftest:", "PASS" if bad == 0 else f"{bad} FAILURES")
    return bad == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--json")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)
    here = os.path.dirname(os.path.abspath(__file__))
    files = args.files or sorted(glob.glob(os.path.join(here, "..", "certificates", "k*.json")))
    print("numpy", np.__version__, "| python", sys.version.split()[0])
    print("finite checks = code and symmetry facts verified exactly here (self-duality, rowspace hash, index, symmetry,")
    print("balance, quaternion relations).  basis = what the orbifold index rests on beyond them: 'GJF 5.6' (cited theorem,")
    print("hypotheses checked), 'Q8 lemma' (conditional on the anomaly lemma of details F.3), 'm = 1' (no orbifold).")
    hdr = ("  k |    c | team   | length | self-dual | code index | R ground states (even/odd) | symmetry | orbifold index | basis"
           " | TMF min | finite checks | method | CPU s")
    print(hdr); print("-" * len(hdr))
    results = []
    for f in files:
        r = verify_certificate(f)
        results.append(r)
        st = r.get("orbifold_index_status", "")
        basis = "Q8 lemma" if st.startswith("orbifold_index_conditional") else ("m = 1" if st.startswith("no orbifold") else "GJF 5.6")
        print(f"{r['k']:>3} | {12 * r['k']:>4} | {r['found_by']:<6} | {r['length']:>6} | {str(r['self_dual']):<9} | {r.get('index')!s:>10} | "
              f"{r.get('even_odd', ['?', '?'])[0]}/{r.get('even_odd', ['?', '?'])[1]} | {r.get('symmetry', '?')} | {r.get('orbifold_index')!s:>4} | "
              f"{basis} | {r.get('tmf_min')!s:>3} | {'pass' if r.get('finite_checks_pass') else 'FAIL'} | {r.get('index_method', '-')} | "
              f"{r.get('cpu_seconds')}   [{r['file']}]")
    print()
    print("Rows without a new code: cited constructions + index arithmetic (no finite certificate; the theories exist by")
    print("the cited constructions: symmetric squares, tensor products, GJF's k = 3 theory).")
    print("  k |    c | type             | index | TMF min | arithmetic | how")
    for line in derived_rows(results):
        print(line)
    nbad = sum(not r["ok"] for r in results)
    print()
    distinct = len({r.get("rowspace_sha256") for r in results})
    nfin = sum(bool(r.get("finite_checks_pass")) for r in results)
    from collections import Counter
    basis = Counter(("Q8 lemma" if r.get("orbifold_index_status", "").startswith("orbifold_index_conditional") else
                     "m = 1" if r.get("orbifold_index_status", "").startswith("no orbifold") else "GJF 5.6") for r in results)
    print(f"{len(results)} certificates ({distinct} distinct codes by rowspace): finite checks pass for {nfin}; "
          f"orbifold index = TMF minimum for {len(results) - nbad} "
          f"(basis: {basis['GJF 5.6']} GJF Thm 5.6, {basis['Q8 lemma']} conditional on the Q8 anomaly lemma, {basis['m = 1']} m = 1); "
          f"{nbad} failed")
    if args.json:
        json.dump(results, open(args.json, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
