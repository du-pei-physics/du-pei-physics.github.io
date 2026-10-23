"""Q_8-equivariant product rings (parent's claude/k8/q8_cells.py cells) in the k24_transfer engine.

k = 16: rings of EIGHT different Q_8-equivariant 12-cells with a common invariant U (r = 4),
|Tr(J T_1...T_8)| = 24 (negacyclic) or |Tr(T_1...T_8)| = 24 (cyclic); orbifold by the internal
Q_8 (i^2 = j^2 = k^2 = -1) would have index 24/8 = 3 if non-anomalous.  Also records v_2 statistics.
Verification (verify_q8_ring): self-dual, i and j preserve the code (so Q_8 does), no fixed
codewords, index by (1) block products of core.transfer tables and (2) the parent's independent
generator-based DP q8_cells.index_from_generator.
usage: python q8_rings.py --minutes 5 --ncell 8 --seed 1 --out results/q8_k16_probe.json
"""
import argparse, json, sys, time
from collections import Counter
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "k8")); sys.path.insert(0, str(HERE.parent))
import core, kspan_rings as kr
import q8_cells as q8
from f3lib import is_self_dual, preserves, fixed_subcode

SEED_S = q8.seed_anti()
HBLK = lambda G, L: q8.block_diag(*([G] * (2 * L)))


def v2(x):
    k = 0
    while x % 2 == 0:
        x //= 2; k += 1
    return k


def q8_pool(rng, N):
    B = q8.random_invariant_plane(rng)
    if B is None:
        return None
    specs, blks, shapes, seen = [], [], None, set()
    for _ in range(6 * N):
        if len(specs) >= N:
            break
        S = q8.random_anti_equivariant(rng, SEED_S)
        T, neg, norm = core.transfer(B, S)
        if not T.any() or T.tobytes() in seen:
            continue
        seen.add(T.tobytes())
        bl = core.blocks(T, neg, norm); shapes = [(a, b, m.shape[0]) for a, b, m in bl]
        specs.append(S); blks.append([m for _, _, m in bl])
    if len(specs) < 8:
        return None
    return {"B": B, "r": B.shape[0], "specs": specs, "stacks": [np.array([b[k] for b in blks], dtype=np.int64) for k in range(len(shapes))],
            "signs": [j for _, j, _ in shapes]}


def verify_q8_ring(B, cells, twist):
    L = len(cells); cyclic = twist == "cyc"
    G = q8.ring_generator(B, cells, cyclic)
    HI, HJ = HBLK(q8.GI, L), HBLK(q8.GJ, L)
    rep = {"length": int(G.shape[1]), "k": int(G.shape[1] // 12), "twist": twist,
           "cells_anti_orthogonal": bool(all(q8.is_anti(S) for S in cells)),
           "cells_equivariant": bool(all(np.array_equal((Gq @ S) % 3, (S @ Gq) % 3) for S in cells for Gq in (q8.GI, q8.GJ))),
           "U_invariant": bool(core.rank3(np.vstack([B, (B @ q8.GI) % 3, (B @ q8.GJ) % 3])) == B.shape[0]),
           "self_dual": bool(is_self_dual(G)), "Q8_preserves": bool(preserves(G, HI) and preserves(G, HJ)),
           "fixed_subcode_dims_i_j": [int(fixed_subcode(G, H).shape[0]) for H in (HI, HJ)]}
    # method 1: my tables, exact block products
    Ts = [core.transfer(B, S)[0] for S in cells]; neg = core.transfer(B, cells[0])[1]
    from search_composite import exact_trJ
    from verify_hit import _tr_plain
    rep["index_tables"] = exact_trJ(Ts, neg, 1) if not cyclic else _tr_plain(Ts, 1)
    # method 2: parent's generator-based DP (independent implementation)
    rep["index_generator_DP"] = q8.index_from_generator(G, q8.Q, L)
    rep["ok"] = bool(rep["cells_anti_orthogonal"] and rep["cells_equivariant"] and rep["U_invariant"] and rep["self_dual"] and rep["Q8_preserves"]
                     and rep["fixed_subcode_dims_i_j"] == [0, 0] and rep["index_tables"] == rep["index_generator_DP"] and abs(rep["index_tables"]) == 24)
    return rep


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--minutes", type=float, default=5); ap.add_argument("--ncell", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--pool", type=int, default=120); ap.add_argument("--out", required=True)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed); t_end = time.time() + 60 * a.minutes
    # sanity: my tables == parent's tables
    P = q8_pool(rng, 10); assert P is not None
    T1, neg1 = q8.table(P["B"], P["specs"][0]); T2, neg2, _ = core.transfer(P["B"], P["specs"][0])
    assert np.array_equal(T1, T2) and np.array_equal(neg1, neg2), "table mismatch with q8_cells"
    st = {tw: {"tried": 0, "v2": Counter(), "small": Counter(), "min": None, "hits": 0} for tw in ("nega", "cyc")}
    hits = []; pools = 0
    while time.time() < t_end:
        P = q8_pool(rng, a.pool)
        if P is None:
            continue
        pools += 1; N = len(P["specs"])
        for _ in range(4):
            idx = rng.integers(0, N, (20000, a.ncell))
            for tw in ("nega", "cyc"):
                signs = P["signs"] if tw == "nega" else [1] * len(P["signs"])
                s = st[tw]; s["tried"] += len(idx)
                vals = kr.tr_crt(P["stacks"], signs, idx[:3000], 1)
                for v in vals:
                    if v: s["v2"][v2(abs(v))] += 1
                    if v and abs(v) <= 240: s["small"][abs(v)] += 1
                    if v and (s["min"] is None or abs(v) < s["min"]): s["min"] = abs(v)
                v0 = kr.tr_mod(P["stacks"], signs, idx, 1, kr.PR[0]); cand = idx[(v0 == 24) | (v0 == kr.PR[0] - 24)]
                if len(cand):
                    v1 = kr.tr_mod(P["stacks"], signs, cand, 1, kr.PR[1]); cand = cand[(v1 == 24) | (v1 == kr.PR[1] - 24)]
                for row in cand:
                    v = kr.tr_exact(P["stacks"], signs, row, 1)
                    if abs(v) != 24: continue
                    s["hits"] += 1
                    cells = [P["specs"][int(c)] for c in row]
                    rep = verify_q8_ring(P["B"], cells, tw)
                    hits.append({"k": 6 * a.ncell // 3 * 1, "ncell": a.ncell, "twist": tw, "B": P["B"].tolist(), "cells": [S.tolist() for S in cells], "index": v, "verification": rep})
                    Path(a.out).write_text(json.dumps({"partial": True, "hits": hits}, indent=1, default=str))
    summ = {"ncell": a.ncell, "pools": pools, "per_twist": {tw: {"tried": s["tried"], "v2_hist": dict(sorted(s["v2"].items())), "min_abs": s["min"],
            "small_values": dict(sorted(s["small"].items())), "hits": s["hits"]} for tw, s in st.items()},
            "verified_hits": sum(1 for h in hits if h["verification"]["ok"])}
    Path(a.out).write_text(json.dumps({"summary": summ, "hits": hits}, indent=1, default=str))
    print(json.dumps(summ))


if __name__ == "__main__":
    main()
