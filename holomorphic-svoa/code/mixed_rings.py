# Copied from claude/persistent/mixed_rings/run.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Persistent, checkpointed search for mixed-width period rings (k = 15, 18, 21: the rows that need an
order-3 or order-6 shift symmetry with balanced fixed words).  Same arithmetic as kspan/mixed_products.py:
cells of widths 12 and 6 with a common delay count r = 3, M = T_1 ... T_p for one period, three periods;
conditions  Tr(J M) = 0 (nega) / Tr(M) = 0 (cyc)  and  |Tr(J M^3)| = 24 / |Tr(M^3)| = 24.
Hits are verified inline from the rebuilt generator (self-duality, rank, generator-only banded DP for the
index, fixed subcode of the GJF shift g and its balance when the contracted code is small enough) and
appended to hits.jsonl; stop.flag honoured; state in state.json; resume continues the seed sequence.

usage (from this folder):
  nohup nice -n 19 ../../../.venv/bin/python run.py --name mixed_k15 --k 15 --widths 12,12,6 --seed 1 > ../mixed_k15/log.txt 2>&1 &
  --k 18 --widths 12,12,12 ;  --k 21 --widths 12,12,12,6
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from worker_common import Worker
PRIMES = (1048573.0, 1048571.0)   # two primes below 2^20 for the fast screen
from mixed_products import random_cell, trJ
from cells_general import anti_bases
from verify_mixed_hit import generator_from_cells, banded_index
from f3lib import is_self_dual, rank, monomial_matrix, preserves, order_of, fixed_subcode, rref, index_exhaustive


def shift_monomial(widths, periods, cyclic, step):
    """Monomial map of the ring code: shift by `step` periods on both halves (info | parity);
    negacyclic wrap contributes a sign -1 per wrap.  Returns (perm, signs) in f3lib convention."""
    per = sum(widths); d = per * periods; n = 2 * d
    perm = np.zeros(n, dtype=np.int64); signs = np.ones(n, dtype=np.int64)
    for half in (0, d):
        for j in range(periods):
            for t in range(per):
                src = half + j * per + t; dst_j = (j + step) % periods; dst = half + dst_j * per + t
                perm[src] = dst
                if not cyclic and (j + step) >= periods:
                    signs[src] = 2
    return perm, signs


def verify(cells, r, cyclic, periods, m):
    G, widths = generator_from_cells(cells, r, cyclic)
    rep = {"length": int(G.shape[1]), "widths": widths, "rank": int(rank(G)), "self_dual": bool(is_self_dual(G))}
    rep["index_banded_DP"] = banded_index(G, widths)
    # GJF shift: nega g = z^{2L/m} with L = periods = 3 -> step 2 if m = 3, step 1 if m = 6; cyc g = z (m=3) or -z (m=6)
    step = {3: 2, 6: 1}[m] if not cyclic else 1
    wp = widths[: len(widths) // periods]          # one period's widths (widths lists all periods)
    perm, signs = shift_monomial(wp, periods, cyclic, step)
    if cyclic and m == 6:
        signs = 3 - signs  # -z
    g = monomial_matrix(perm.tolist(), signs.tolist())
    rep["g_preserves"] = bool(preserves(G, g)); rep["g_order"] = int(order_of(g))
    # GJF balance condition concerns h = g^{m/3} (order 3; for m = 6 the order-6 g itself fixes nothing)
    h = g if m == 3 else (g @ g) % 3
    rep["h_order"] = int(order_of(h))
    F = fixed_subcode(G, h); rep["fixed_dim"] = int(F.shape[0])
    per = sum(widths) // periods          # widths is the full list (all periods); one period's span
    cols = list(range(0, per)) + list(range(G.shape[0], G.shape[0] + per))
    Fc = rref(F[:, cols] % 3)[0]
    rep["contracted_rank"] = int(Fc.shape[0]); rep["contracted_self_dual"] = bool(is_self_dual(Fc))
    if Fc.shape[0] <= 30:
        res = index_exhaustive(Fc, workers=1, prefix_bits=3, tag="mixedbal")
        rep["balance_exhaustive"] = res; rep["balanced"] = (res["index"] == 0) if isinstance(res, dict) else None
    else:
        try:
            rep["balance_banded_DP"] = banded_index(Fc, widths[: len(widths) // periods]); rep["balanced"] = rep["balance_banded_DP"] == 0
        except AssertionError:
            rep["balanced"] = None; rep["balance_note"] = "contracted generator not in banded form; verify by hand"
    rep["ok"] = bool(rep["self_dual"] and rep["rank"] == G.shape[0] and abs(rep["index_banded_DP"]) == 24 and rep["g_preserves"] and rep["balanced"])
    return rep, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True); ap.add_argument("--k", type=int, required=True); ap.add_argument("--widths", required=True)
    ap.add_argument("--periods", type=int, default=3); ap.add_argument("--r", type=int, default=3); ap.add_argument("--pool", type=int, default=60)
    ap.add_argument("--refresh", type=int, default=300000); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--max-minutes", type=float, default=None)
    ap.add_argument("--exclude-bases", default="", help="comma-separated base names to drop from the cell pools, e.g. A0")
    ap.add_argument("--metric2", type=int, default=0, help="1 = every cell uses a discriminant-2 delay space (ring-level choice)")
    ap.add_argument("--hits-max", type=int, default=0, help="store at most this many hit records (0 = all); further hits are only counted")
    ap.add_argument("--fast", type=int, default=0, help="1 = screen products modulo two primes with float matmul before exact evaluation (for r = 4)")
    ap.add_argument("--near-abs", type=int, default=72, help="archive exact descriptors of balanced products with 0 < |index| <= this (near hits)")
    ap.add_argument("--near-max", type=int, default=300, help="at most this many near hits per run")
    ap.add_argument("--s-refl", type=int, default=0, help="up to this many random reflections on S")
    a = ap.parse_args()
    widths = [int(w) for w in a.widths.split(",")]
    assert sum(widths) * a.periods == 6 * a.k
    from math import gcd
    m = gcd(a.k, 24)
    assert (a.periods == 3 and m in (3, 6)) or (a.periods == 1 and m in (1, 2)), "periods 3 for m in {3,6}, periods 1 for m in {1,2}"
    prev_stats = None
    try:
        sp = HERE / "out" / a.name / "state.json"
        prev_stats = json.loads(sp.read_text()).get("stats") if sp.exists() else None
    except Exception:
        prev_stats = None
    w = Worker(HERE / "out" / a.name, a.name)
    rng = np.random.default_rng(a.seed + 1000 * w.state.get("previous_runs", 0))
    excl = set(b for b in a.exclude_bases.split(",") if b)
    bases = {q: {k: v for k, v in anti_bases(q).items() if k not in excl} for q in set(widths)}
    assert all(bases[q] for q in bases), "a width has no bases left"
    t_end = time.time() + 60 * a.max_minutes if a.max_minutes else None
    st = prev_stats if prev_stats else w.state.get("stats", {"tried": 0, "balanced": {"nega": 0, "cyc": 0}, "min_abs": {"nega": None, "cyc": None},
                                "small": {"nega": {}, "cyc": {}}, "hits": {"nega": 0, "cyc": 0}, "verified": 0})
    cpu_prev = float(st.get("cpu_total", 0.0)); st.setdefault("cpu_total", cpu_prev)
    pools = None; n = 3 ** a.r; since_refresh = 0; refresh_at = 0; seen = set(); last_summary = time.time()
    while not w.should_stop() and (t_end is None or time.time() < t_end):
        w.wait_for_disk()
        if pools is None or since_refresh >= refresh_at:
            pools = {q: [random_cell(rng, q, a.r, bases[q], a.metric2, a.s_refl) for _ in range(a.pool)] for q in set(widths)}
            neg = pools[widths[0]][0][1]; since_refresh = 0; seen = set()
            # de-duplication (18 Oct): a pool of P cells gives only P^len(widths) distinct candidates; refresh after
            # at most half of them so that the sampler does not keep re-testing the same codes
            refresh_at = max(1, min(a.refresh, (a.pool ** len(widths)) // 2))
            if a.fast:
                pool_f = [{q: [np.asarray(np.mod(T, P), dtype=np.float64) for (T, _, _) in pools[q]] for q in pools} for P in PRIMES]
        picks = [(q, int(rng.integers(len(pools[q])))) for q in widths]
        key = tuple(i for _, i in picks)
        if key in seen:
            st["dup"] = st.get("dup", 0) + 1
            continue
        seen.add(key)
        st["tried"] += 1; w.state["attempts"] += 1; since_refresh += 1
        if st["tried"] % 500 == 0:
            st["cpu_seconds"] = round(time.process_time(), 1); st["cpu_total"] = round(cpu_prev + time.process_time(), 1); w.save(stats=st)
        if a.fast:
            # screen modulo two primes with float64 matmul (entries < 2^20, n <= 81: sums < 2^53, exact)
            cand = True
            for P, pf in zip(PRIMES, pool_f):
                Mp = np.eye(n)
                for q, i in picks:
                    Mp = np.fmod(Mp @ pf[q][i], P)
                tn = np.fmod(Mp[neg, np.arange(n)].sum(), P); tc = np.fmod(np.trace(Mp), P)
                if a.periods == 1:
                    ok_n = tn in (24.0, P - 24.0); ok_c = tc in (24.0, P - 24.0)
                else:
                    ok_n = tn == 0; ok_c = tc == 0
                if not ok_n and not ok_c:
                    cand = False; break
            if not cand:
                continue
        M = np.eye(n, dtype=object); M3 = None
        for q, i in picks:
            M = M.dot(pools[q][i][0])
        if a.periods == 1:
            # single-period rows (m = 1, 2): the index is Tr(JM) / Tr(M) itself
            vn, vc = trJ(M, neg), int(np.trace(M)); bal_n = bal_c = True; M3 = None
        else:
            bal_n = trJ(M, neg) == 0; bal_c = int(np.trace(M)) == 0
        if bal_n or bal_c:
            if M3 is not None or a.periods != 1:
                M3 = M.dot(M).dot(M)
            for kind, bal, val in (("nega", bal_n, (lambda: vn) if a.periods == 1 else (lambda: trJ(M3, neg))), ("cyc", bal_c, (lambda: vc) if a.periods == 1 else (lambda: int(np.trace(M3))))):
                if not bal:
                    continue
                st["balanced"][kind] += 1; v = val(); av = abs(v)
                if av:
                    if st["min_abs"][kind] is None or av < st["min_abs"][kind]:
                        st["min_abs"][kind] = av
                    if av <= 240:
                        st["small"][kind][str(av)] = st["small"][kind].get(str(av), 0) + 1
                    # near-hit archive (Codex's suggestion): keep the exact descriptors of the smallest
                    # nonzero values so that their neighbourhoods can be explored; bounded to a.near_max entries
                    if av <= a.near_abs and st.get("near_saved", 0) < a.near_max and av != 24:
                        st["near_saved"] = st.get("near_saved", 0) + 1
                        w.append("near_hits.jsonl", {"k": a.k, "kind": kind, "widths": widths, "periods": a.periods, "r": a.r,
                                                     "index": v, "cells": [pools[q][i][2] for q, i in picks],
                                                     "time": time.strftime("%Y-%m-%d %H:%M:%S")})
                if av == 24:
                    cells = [pools[q][i][2] for q, i in picks]
                    if a.hits_max and w.state["hits"] + 1 > a.hits_max:
                        st["hits"][kind] += 1; w.state["hits"] += 1; continue
                    try:
                        if a.periods == 1:
                            G, W = generator_from_cells(cells, a.r, kind == "cyc")
                            rep = {"length": int(G.shape[1]), "rank": int(rank(G)), "self_dual": bool(is_self_dual(G)), "index_banded_DP": banded_index(G, W)}
                            rep["ok"] = bool(rep["self_dual"] and rep["rank"] == G.shape[0] and abs(rep["index_banded_DP"]) == 24)
                        else:
                            rep, G = verify(cells * a.periods, a.r, kind == "cyc", a.periods, m)
                        Gl = G.tolist()
                    except Exception as e:                       # never lose a hit to a verification bug
                        rep = {"ok": False, "error": repr(e)}; Gl = None
                    st["hits"][kind] += 1; st["verified"] += int(bool(rep.get("ok"))); w.state["hits"] += 1
                    if a.hits_max and w.state["hits"] > a.hits_max:
                        continue
                    w.append("hits.jsonl", {"k": a.k, "kind": kind, "widths": widths, "periods": a.periods, "r": a.r, "index": v,
                                            "cells": cells, "generator": Gl, "verification": rep,
                                            "time": time.strftime("%Y-%m-%d %H:%M:%S")})
                    w.save(force=True, stats=st)
        if st["tried"] % 500 == 0:
            w.save(stats=st)
        if time.time() - last_summary > 600:
            w.append("summary.jsonl", {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "stats": st}); last_summary = time.time()
    st["cpu_seconds"] = round(time.process_time(), 1); st["cpu_total"] = round(cpu_prev + time.process_time(), 1); w.save(force=True, stats=st)
    w.finish("stopped by stop.flag" if w.should_stop() else "finished (time limit)")


if __name__ == "__main__":
    main()
