"""Benchmark worker v3 (19 Oct 2026): CPU time per verified, distinct code for the chain-of-cells search,
with the accounting asked for in Codex's publication review (docs/publication_review_2026_10_18.md, item 1
and "worker improvements").

Per lane (folder persistent/<name>/) it keeps, cumulatively over runs:
  stats.draws        all random draws of a cell tuple
  stats.dup          draws rejected because the tuple was already tested in the current pool
  stats.tried        distinct cell tuples tested (= candidates)
  stats.cand_hit     candidates with at least one closure (periodic / anti-periodic) at the target
  stats.hits         {nega, cyc}: exact transfer hits per closure
  stats.verified     closures re-verified from the expanded generator (EVERY hit is verified)
  stats.verify_fail  hits whose verification failed (reported, not counted as codes)
  stats.distinct     distinct RREF rowspace hashes among verified hits, across all pools and runs
  stats.cpu_total, stats.cpu_verify, stats.cpu_search (= total - verify), process CPU seconds
Files:
  runs.jsonl       one line per run: time, pid, seed, argv, source SHA-256s, python / numpy versions
  epochs.jsonl     one line per pool: epoch, run, RNG state before the pool was drawn (pools are reproducible)
  hit_index.jsonl  one line per hit: time, run, epoch, cell-index tuple, closure, index, rowspace hash, verified
  hits.jsonl       full records (cells and generator) for the first --hits-max hits only
  rowspaces.txt    the distinct rowspace hashes (SHA-256 of the RREF over F_3, rows over {0,1,2} joined by "\n",
                   the definition used by the online note)
Distinct rowspaces do not establish inequivalent codes or SCFTs.  A lane stops when its stop.flag or the global
stop file exists (checked every few seconds, not only between runs).
"""
import argparse, hashlib, json, os, platform, sys, time
from math import gcd
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
for p in (HERE.parent, HERE.parent.parent / "kspan", HERE.parent.parent, HERE):
    sys.path.insert(0, str(p))
from common import Worker
from mixed_products import random_cell, trJ
from cells_general import anti_bases
from verify_mixed_hit import generator_from_cells, banded_index
from f3lib import is_self_dual, rank, rref
from run import verify as verify_period, PRIMES
from journal import reconcile

SOURCES = [HERE / "bench_worker.py", HERE / "journal.py", HERE / "run.py", HERE.parent.parent / "kspan" / "mixed_products.py",
           HERE.parent.parent / "kspan" / "cells_general.py", HERE.parent.parent / "kspan" / "verify_mixed_hit.py",
           HERE.parent.parent / "f3lib.py"]


def rowspace_sha256(G):
    R = rref(np.asarray(G) % 3)[0]
    return hashlib.sha256("\n".join("".join(str(int(v)) for v in row) for row in R).encode()).hexdigest()


def signed_check(rep, v):
    """Codex review 19 Oct, item 5: the verifier's signed index must equal the signed transfer value v."""
    d = rep.get("index_banded_DP")
    if isinstance(d, (list, tuple)):
        d = d[0]
    rep["transfer_index"] = int(v)
    rep["signed_match"] = d is not None and int(d) == int(v)
    rep["ok"] = bool(rep.get("ok") and rep["signed_match"])
    return rep


def verify_single_period(cells, r, cyclic):
    G, W = generator_from_cells(cells, r, cyclic)
    rep = {"length": int(G.shape[1]), "rank": int(rank(G)), "self_dual": bool(is_self_dual(G)),
           "index_banded_DP": banded_index(G, W)}
    rep["ok"] = bool(rep["self_dual"] and rep["rank"] == G.shape[0] and abs(rep["index_banded_DP"]) == 24)
    return rep, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True); ap.add_argument("--k", type=int, required=True); ap.add_argument("--widths", required=True)
    ap.add_argument("--periods", type=int, default=3); ap.add_argument("--r", type=int, default=3); ap.add_argument("--pool", type=int, default=40)
    ap.add_argument("--refresh", type=int, default=200000); ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--max-minutes", type=float, default=None); ap.add_argument("--exclude-bases", default="")
    ap.add_argument("--metric2", type=int, default=0); ap.add_argument("--s-refl", type=int, default=0); ap.add_argument("--fast", type=int, default=0)
    ap.add_argument("--hits-max", type=int, default=100, help="full records stored for the first this many hits; all hits are verified")
    ap.add_argument("--global-stop", default=None, help="path of a stop file shared by all lanes")
    a = ap.parse_args()
    widths = [int(w) for w in a.widths.split(",")]
    assert sum(widths) * a.periods == 6 * a.k, "cell widths times periods must equal 6k"
    m = gcd(a.k, 24)
    assert (a.periods == 3 and m in (3, 6)) or (a.periods == 1 and m in (1, 2)), "periods 3 for m in {3,6}, 1 for m in {1,2}"
    folder = HERE.parent / a.name; folder.mkdir(parents=True, exist_ok=True)
    prev = {}
    if (folder / "state.json").exists():
        try:
            prev = json.loads((folder / "state.json").read_text()).get("stats", {}) or {}
        except Exception:
            prev = {}
    w = Worker(folder, a.name)
    run_no = w.state.get("previous_runs", 0)
    st = {"draws": 0, "dup": 0, "tried": 0, "cand_hit": 0, "hits": {"nega": 0, "cyc": 0}, "verified": 0, "verify_fail": 0,
          "distinct": 0, "cpu_total": 0.0, "cpu_verify": 0.0, "cpu_search": 0.0, "min_abs": {"nega": None, "cyc": None}, "small": {"nega": {}, "cyc": {}}}
    st.update(prev)
    seen_rowspaces = set()
    if (folder / "rowspaces.txt").exists():
        seen_rowspaces = set(x.strip() for x in (folder / "rowspaces.txt").read_text().split() if x.strip())
    st["distinct"] = len(seen_rowspaces)
    reconciled, seen_rowspaces = reconcile(folder, st)                 # counters from the durable journals (item 6)
    cpu_prev_total, cpu_prev_verify = float(st["cpu_total"]), float(st["cpu_verify"])
    cpu_verify_run = 0.0
    w.append("runs.jsonl", {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "run": run_no, "pid": os.getpid(), "seed": a.seed, "argv": sys.argv,
                            "python": platform.python_version(), "numpy": np.__version__, "reconciled": reconciled,
                            "sources": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES if p.exists()}})
    rng = np.random.default_rng(a.seed)
    excl = set(b for b in a.exclude_bases.split(",") if b)
    bases = {q: {k: v for k, v in anti_bases(q).items() if k not in excl} for q in set(widths)}
    assert all(bases[q] for q in bases)
    n = 3 ** a.r
    t_end = time.time() + 60 * a.max_minutes if a.max_minutes else None
    epoch = int(prev.get("epochs", 0)); pools = None; since_refresh = 0; refresh_at = 0; seen = set(); full_stored = int(prev.get("full_stored", 0))
    last_check = 0.0; last_save = 0.0

    def stop_requested():
        return w.should_stop() or (a.global_stop is not None and os.path.exists(a.global_stop))

    def save(force=False):
        st["cpu_total"] = round(cpu_prev_total + time.process_time(), 2)
        st["cpu_verify"] = round(cpu_prev_verify + cpu_verify_run, 2)
        st["cpu_search"] = round(st["cpu_total"] - st["cpu_verify"], 2)
        st["epochs"] = epoch; st["full_stored"] = full_stored
        w.save(force=force, stats=st)

    while True:
        now = time.time()
        if now - last_check > 3:
            last_check = now
            if stop_requested() or (t_end is not None and now >= t_end):
                break
            w.wait_for_disk()
        if pools is None or since_refresh >= refresh_at:
            epoch += 1
            w.append("epochs.jsonl", {"epoch": epoch, "run": run_no, "time": time.strftime("%Y-%m-%d %H:%M:%S"), "rng_state": rng.bit_generator.state,
                                      "pool": a.pool, "widths": widths, "r": a.r, "metric2": a.metric2, "s_refl": a.s_refl, "exclude_bases": sorted(excl)})
            pools = {q: [random_cell(rng, q, a.r, bases[q], a.metric2, a.s_refl) for _ in range(a.pool)] for q in sorted(set(widths))}
            neg = pools[widths[0]][0][1]; since_refresh = 0; seen = set()
            refresh_at = max(1, min(a.refresh, (a.pool ** len(widths)) // 2))
            if a.fast:
                pool_f = [{q: [np.asarray(np.mod(T, P), dtype=np.float64) for (T, _, _) in pools[q]] for q in pools} for P in PRIMES]
        picks = [(q, int(rng.integers(len(pools[q])))) for q in widths]
        st["draws"] += 1
        key = tuple(i for _, i in picks)
        if key in seen:
            st["dup"] += 1
            continue
        seen.add(key); st["tried"] += 1; since_refresh += 1; w.state["attempts"] += 1
        if a.fast:
            cand = False
            for P, pf in zip(PRIMES, pool_f):
                Mp = np.eye(n)
                for q, i in picks:
                    Mp = np.fmod(Mp @ pf[q][i], P)
                tn = np.fmod(Mp[neg, np.arange(n)].sum(), P); tc = np.fmod(np.trace(Mp), P)
                if a.periods == 1:
                    ok_n, ok_c = tn in (24.0, P - 24.0), tc in (24.0, P - 24.0)
                else:
                    ok_n, ok_c = tn == 0, tc == 0
                cand = ok_n or ok_c
                if not cand:
                    break
            if not cand:
                if time.time() - last_save > 20:
                    save(); last_save = time.time()
                continue
        M = np.eye(n, dtype=object)
        for q, i in picks:
            M = M.dot(pools[q][i][0])
        if a.periods == 1:
            vals = {"nega": trJ(M, neg), "cyc": int(np.trace(M))}
        else:
            vals = {}
            if trJ(M, neg) == 0 or int(np.trace(M)) == 0:
                M3 = M.dot(M).dot(M)
                if trJ(M, neg) == 0:
                    vals["nega"] = trJ(M3, neg)
                if int(np.trace(M)) == 0:
                    vals["cyc"] = int(np.trace(M3))
        any_hit = False
        for kind, v in vals.items():
            av = abs(v)
            if av:
                if st["min_abs"][kind] is None or av < st["min_abs"][kind]:
                    st["min_abs"][kind] = av
                if av <= 240:
                    st["small"][kind][str(av)] = st["small"][kind].get(str(av), 0) + 1
            if av != 24:
                continue
            any_hit = True; st["hits"][kind] += 1; w.state["hits"] += 1
            cells = [pools[q][i][2] for q, i in picks]
            t0 = time.process_time()
            try:
                if a.periods == 1:
                    rep, G = verify_single_period(cells, a.r, kind == "cyc")
                else:
                    rep, G = verify_period(cells * a.periods, a.r, kind == "cyc", a.periods, m)
                rep = signed_check(rep, v)
                sha = rowspace_sha256(G) if rep.get("ok") else None
            except Exception as e:
                rep, G, sha = {"ok": False, "error": repr(e)}, None, None
            cpu_verify_run += time.process_time() - t0
            if rep.get("ok"):
                st["verified"] += 1
                if sha not in seen_rowspaces:
                    seen_rowspaces.add(sha); st["distinct"] += 1
                    with open(folder / "rowspaces.txt", "a") as f:
                        f.write(sha + "\n")
            else:
                st["verify_fail"] += 1
            w.append("hit_index.jsonl", {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "run": run_no, "epoch": epoch, "tuple": list(key),
                                         "kind": kind, "index": v, "rowspace_sha256": sha, "verified": bool(rep.get("ok"))})
            if full_stored < a.hits_max:
                full_stored += 1
                w.append("hits.jsonl", {"k": a.k, "kind": kind, "widths": widths, "periods": a.periods, "r": a.r, "index": v, "epoch": epoch,
                                        "tuple": list(key), "cells": cells, "generator": None if G is None else G.tolist(), "verification": rep,
                                        "rowspace_sha256": sha, "time": time.strftime("%Y-%m-%d %H:%M:%S")})
        if any_hit:
            st["cand_hit"] += 1
            save(force=True); last_save = time.time()
        elif time.time() - last_save > 20:
            save(); last_save = time.time()
    save(force=True)
    w.finish("stopped" if stop_requested() else "finished (time limit)")


if __name__ == "__main__":
    main()
