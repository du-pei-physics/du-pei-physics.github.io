"""Product rings for the minimal-index targets at c = 12k, k = 8..23
(see NOTE_kspan.md).  Length 12k = 2qN: N cells of size q with a common delay space U, period p,
M = T_{S_1} ... T_{S_p}, L = N/p periods.  m = gcd(k, 24), target index 24/m.

Conditions (negacyclic shift z of order 2N; cyclic z of order N):
  A  (m = 1, 2):  |Tr(J M)| = 24  (p = N, one period)  -- or cyclic |Tr M| = 24.  g = -1 only.
  B  (m = 3, 6):  Tr(J M) = 0 and |Tr(J M^3)| = 24 (3 periods); cyclic: Tr M = 0, |Tr M^3| = 24.
                  m = 3: g = z^{2p} (nega) / z^p (cyc), order 3.  m = 6: g = z^p (nega) / -z^p (cyc), order 6.
  C  (m = 4):     |Tr(J M^2)| = 24 (2 periods), g = z^p of order 4, g^2 = -1; negacyclic only.
Per family: pools of distinct tables for random U, random p-tuples screened mod two primes on
the (norm x J) blocks, exact check of survivors, inline verification of hits on the expanded
generator (verify_hit.verify_kspan), and exact small-value statistics on a sample.
"""
import json, time
from collections import Counter
import numpy as np
import core
from bases import BASES12, conf6x2
from search_p4 import PR, trpow_mod, prod_mod

PRS = (1048573, 1048571, 1048549)                 # three primes < 2^20; CRT exact for |value| < 2^40
MODP = PRS[0] * PRS[1] * PRS[2]
BIG = 1 << 40

A06 = core.A0(6); A06D = A06.copy(); A06D[:, 0] = (-A06D[:, 0]) % 3
BASES6 = {"A0": A06, "A0D": A06D, "C6": conf6x2()[:6, :6]}

FAMILIES = {
    "k8_sigma_p2":  dict(k=8,  q=12, N=4,  p=2,  m=8, cond="D", sigma=True),   # g = h z^2, order 8, |Tr J M^2| = 24
    "k16_sigma_p4": dict(k=16, q=12, N=8,  p=4,  m=8, cond="D", sigma=True),   # g = h z^4, order 8, |Tr J M^2| = 24
    "k10_q12":    dict(k=10, q=12, N=5,  p=5,  m=2, cond="A"),
    "k18_q12_p3": dict(k=18, q=12, N=9,  p=3,  m=6, cond="B"),
    "k14_q12":    dict(k=14, q=12, N=7,  p=7,  m=2, cond="A"),
    "k20_q12_p5": dict(k=20, q=12, N=10, p=5,  m=4, cond="C"),
    "k22_q12":    dict(k=22, q=12, N=11, p=11, m=2, cond="A"),
    "k9_q6_p3":   dict(k=9,  q=6,  N=9,  p=3,  m=3, cond="B"),
    "k15_q6_p5":  dict(k=15, q=6,  N=15, p=5,  m=3, cond="B"),
    "k21_q6_p7":  dict(k=21, q=6,  N=21, p=7,  m=3, cond="B"),
    "k11_q6":     dict(k=11, q=6,  N=11, p=11, m=1, cond="A"),
    "k13_q6":     dict(k=13, q=6,  N=13, p=13, m=1, cond="A"),
    "k17_q6":     dict(k=17, q=6,  N=17, p=17, m=1, cond="A"),
    "k19_q6":     dict(k=19, q=6,  N=19, p=19, m=1, cond="A"),
    "k23_q6":     dict(k=23, q=6,  N=23, p=23, m=1, cond="A"),
}
ORDER = ["k8_sigma_p2", "k10_q12", "k18_q12_p3", "k16_sigma_p4", "k14_q12", "k20_q12_p5", "k22_q12", "k9_q6_p3",
         "k11_q6", "k13_q6", "k15_q6_p5", "k17_q6", "k19_q6", "k21_q6_p7", "k23_q6"]
COND_L = {"A": (None, 1), "B": (1, 3), "C": (None, 2), "D": (None, 2)}        # (balance power, index power)


def twists_for(fam):
    return ("nega", "cyc") if fam["cond"] in ("A", "B", "D") else ("nega",)


def build_pool_sigma(rng, N):
    """sigma-equivariant cells with sigma-invariant U (search_sigma.py), for the m = 8 families."""
    import search_sigma as ss
    B = ss.random_invariant_U(rng)
    if B is None:
        return None
    P = ss.equivariant_pool(rng, B, N)
    if P is None:
        return None
    specs, stacks, signs = P
    return {"B": B, "r": B.shape[0], "specs": specs, "stacks": stacks, "signs": signs, "neg": core.transfer(B, specs[0])[1], "sigma": ss.SIGMA}


def build_pool(rng, fam, N, dmax=6):
    if fam.get("sigma"):
        return build_pool_sigma(rng, N)
    q = fam["q"]
    bases = BASES12 if q == 12 else BASES6
    names = list(bases)
    r = int(rng.choice([2, 3, 4, 5] if q == 12 else [2, 3, 4]))
    B = core.random_subspace(rng, r, q)
    specs, blks, shapes, seen, neg0 = [], [], None, set(), None
    for _ in range(10 * N):
        if len(specs) >= N:
            break
        S = core.random_orth(rng, int(rng.integers(0, dmax + 1)), bases[names[int(rng.integers(len(names)))]], q)
        if rng.integers(2):
            S[:, 0] = (-S[:, 0]) % 3
        T, neg, norm = core.transfer(B, S)
        if not T.any():
            continue
        key = T.tobytes()
        if key in seen:
            continue
        seen.add(key); neg0 = neg
        bl = core.blocks(T, neg, norm)
        shapes = [(nc, j, mm.shape[0]) for nc, j, mm in bl]
        specs.append(S); blks.append([mm for _, _, mm in bl])
    if len(specs) < 8:
        return None
    stacks = [np.array([b[k] for b in blks], dtype=np.int64) for k in range(len(shapes))]
    return {"B": B, "r": r, "specs": specs, "stacks": stacks, "signs": [j for _, j, _ in shapes], "neg": neg0}


def tr_mod(stacks, signs, idx, L, p):
    tot = np.zeros(idx.shape[0], dtype=np.int64)
    for st, sg in zip(stacks, signs):
        tot = (tot + sg * trpow_mod(prod_mod(st, idx, p), L, p)) % p
    return tot


def tr_crt(stacks, signs, idx, L):
    """Exact values (Python ints) for a batch idx (K, p) when |value| < 2^40, else None."""
    res = [tr_mod(stacks, signs, idx, L, p) for p in PRS]
    out = []
    for k in range(idx.shape[0]):
        x = 0
        for p, rr in zip(PRS, res):
            mm = MODP // p
            x = (x + int(rr[k]) * mm * pow(mm, -1, p)) % MODP
        x = x - MODP if x > MODP // 2 else x
        out.append(x if abs(x) < BIG else None)
    return out


def tr_exact(stacks, signs, cols, L):
    tot = 0
    for st, sg in zip(stacks, signs):
        M = np.eye(st.shape[1], dtype=object)
        for c in cols:
            M = M.dot(st[int(c)].astype(object))
        X = np.eye(st.shape[1], dtype=object)
        for _ in range(L):
            X = X.dot(M)
        tot += sg * int(np.trace(X))
    return tot


def new_stats():
    return {"pools": 0, "cells": 0, "tried": 0, "balanced": 0, "exact_sampled": 0, "min_abs_index": None,
            "small_values": {}, "hits": 0, "verified_hits": 0}


def run_family(name, rng, seconds, stats, pool_size=150, batch=50000, nstat=12, on_hit=None, verify=True, max_hits=None):
    """Sample products for `seconds`; update `stats` (dict per twist) in place; return list of hits."""
    fam = FAMILIES[name]; p = fam["p"]; Lb, Li = COND_L[fam["cond"]]
    t_end = time.time() + seconds; hits = []
    while time.time() < t_end:
        P = build_pool(rng, fam, pool_size)
        if P is None:
            continue
        Np = len(P["specs"])
        for tw in twists_for(fam):
            st = stats.setdefault(tw, new_stats()); st["pools"] += 1; st["cells"] += Np
        for _round in range(6):
            if time.time() > t_end:
                break
            idx = rng.integers(0, Np, (batch, p))
            for tw in twists_for(fam):
                signs = P["signs"] if tw == "nega" else [1] * len(P["signs"])
                st = stats[tw]; st["tried"] += len(idx)
                sel = idx
                if Lb is not None:
                    b0 = tr_mod(P["stacks"], signs, sel, Lb, PR[0]); sel = sel[b0 == 0]
                    if len(sel):
                        sel = sel[tr_mod(P["stacks"], signs, sel, Lb, PR[1]) == 0]
                    st["balanced"] += len(sel)
                if not len(sel):
                    continue
                # small-value statistics on a sample of the (balanced) tuples, exact via CRT (< 2^40)
                samp = sel[rng.permutation(len(sel))[:nstat]]
                vals = tr_crt(P["stacks"], signs, samp, Li)
                for v in vals:
                    st["exact_sampled"] += 1
                    if v is None or v == 0:
                        continue
                    v = abs(v)
                    if st["min_abs_index"] is None or v < st["min_abs_index"]:
                        st["min_abs_index"] = v
                    if v <= 240:
                        st["small_values"][str(v)] = st["small_values"].get(str(v), 0) + 1
                v0 = tr_mod(P["stacks"], signs, sel, Li, PR[0])
                cand = sel[(v0 == 24) | (v0 == PR[0] - 24)]
                if not len(cand):
                    continue
                v1 = tr_mod(P["stacks"], signs, cand, Li, PR[1])
                cand = cand[(v1 == 24) | (v1 == PR[1] - 24)]
                if len(cand):
                    v2 = tr_mod(P["stacks"], signs, cand, Li, PRS[2])
                    cand = cand[(v2 == 24) | (v2 == PRS[2] - 24)]
                for row in cand:
                    if Lb is not None and tr_exact(P["stacks"], signs, row, Lb) != 0:
                        continue
                    v = tr_exact(P["stacks"], signs, row, Li)
                    if abs(v) != 24:
                        continue
                    st["hits"] += 1
                    hit = {"family": name, "k": fam["k"], "q": fam["q"], "N": fam["N"], "p": p, "m": fam["m"], "twist": tw,
                           "r": P["r"], "B": P["B"].tolist(), "S_period": [P["specs"][int(c)].tolist() for c in row],
                           "index": v, "traces_1..3": [tr_exact(P["stacks"], signs, row, L) for L in (1, 2, 3)]}
                    if fam.get("sigma"):
                        hit["sigma"] = P["sigma"].tolist()
                    if verify:
                        try:
                            from verify_hit import verify_kspan, verify_sigma_order8
                            if fam.get("sigma"):
                                rep = verify_sigma_order8(P["B"], [P["specs"][int(c)] for c in row], P["sigma"], fam["N"] // p, tw)
                            else:
                                rep = verify_kspan(P["B"], [P["specs"][int(c)] for c in row], fam["N"] // p, fam["m"], tw)
                            hit["verification"] = rep
                            if rep.get("ok"):
                                st["verified_hits"] += 1
                        except Exception as e:          # verification failure is recorded, never silent
                            hit["verification"] = {"error": repr(e)}
                    hits.append(hit)
                    if on_hit:
                        on_hit(hit)
                    if max_hits is not None and len(hits) >= max_hits:
                        return hits
    return hits
