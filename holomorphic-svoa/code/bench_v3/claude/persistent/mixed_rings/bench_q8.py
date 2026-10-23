"""Benchmark worker for the Q8-equivariant chains (k = 8, 16), 19 Oct 2026: same accounting and files as
bench_worker.py (benchmark v3), for the search of claude/k8/q8_cells.py (--survey): a random Q8-invariant boundary
space B (r = 4), a pool of random Q8-equivariant anti-orthogonal 12-cells, random L-tuples of cells (L = k/2) with
duplicates rejected within a pool, closures anti-periodic (Tr J M) and periodic (Tr M).

New against the survey: a modular pre-screen (two primes below 2^20, float64 products reduced after each step;
exact, since |index| = 24 implies the residues 24 or P-24), so only screened tuples get the exact integer product;
and every hit is verified inline with the online note's stand-alone verifier (code/verify.py): length, rank and
self-orthogonality, the banded transfer DP index (|index| = 24), and check_q8 with the block-diagonal Q8 action of the
note's certificate of the same length.  Distinct codes = distinct RREF rowspace hashes, written to rowspaces.txt.
"""
import argparse, hashlib, json, os, platform, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
CL = HERE.parent.parent
for p in (HERE.parent, CL / "k8", CL, CL / "online/holomorphic-svoa/code"):
    sys.path.insert(0, str(p))
from common import Worker
import verify as V
from journal import reconcile
from q8_cells import seed_anti, random_anti_equivariant, random_invariant_plane, table, ring_generator, trJ

PRIMES = (1048573.0, 1048571.0)
SOURCES = [HERE / "bench_q8.py", HERE / "journal.py", CL / "k8/q8_cells.py", CL / "online/holomorphic-svoa/code/verify.py"]


def q8_gens(k):
    f = {8: "k08_claude_q8.json", 16: "k16_claude_periodic.json"}[k]
    g = json.load(open(CL / "online/holomorphic-svoa/certificates" / f))["symmetry"]["generators"]
    return {x: (np.array(g[x]["permutation"]), np.array(g[x]["signs"])) for x in ("i", "j")}


def signed_check(rep, v):
    """Codex review 19 Oct, item 5: the verifier's signed index must equal the signed transfer value v."""
    d = rep.get("index_banded_DP")
    if isinstance(d, (list, tuple)):
        d = d[0]
    rep["transfer_index"] = int(v)
    rep["signed_match"] = d is not None and int(d) == int(v)
    rep["ok"] = bool(rep.get("ok") and rep["signed_match"])
    return rep


def verify_hit(B, cells, cyclic, k, gens):
    G = ring_generator(B, cells, cyclic)
    d, n = G.shape
    rep = {"length": int(n), "rank": int(V.rank3(G)), "self_orthogonal": bool(not ((G @ G.T) % 3).any())}
    ci = V.band_index(G, [12] * len(cells))
    rep["index_banded_DP"] = None if ci is None else int(ci[0])
    q = V.check_q8(G, gens)
    rep["q8_ok"] = bool(q["ok"])
    rep["ok"] = bool(n == 12 * k and rep["rank"] == d and 2 * d == n and rep["self_orthogonal"] and ci is not None
                     and abs(ci[0]) == 24 and q["ok"])
    return rep, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True); ap.add_argument("--k", type=int, choices=(8, 16), required=True)
    ap.add_argument("--pool", type=int, default=24); ap.add_argument("--refresh", type=int, default=4000)
    ap.add_argument("--seed", type=int, required=True); ap.add_argument("--max-minutes", type=float, default=None)
    ap.add_argument("--hits-max", type=int, default=100); ap.add_argument("--global-stop", default=None)
    a = ap.parse_args()
    L = a.k // 2
    folder = HERE.parent / a.name; folder.mkdir(parents=True, exist_ok=True)
    prev = {}
    if (folder / "state.json").exists():
        try:
            prev = json.loads((folder / "state.json").read_text()).get("stats", {}) or {}
        except Exception:
            prev = {}
    w = Worker(folder, a.name)
    run_no = w.state.get("previous_runs", 0)
    st = {"draws": 0, "dup": 0, "tried": 0, "screened": 0, "cand_hit": 0, "hits": {"nega": 0, "cyc": 0}, "verified": 0,
          "verify_fail": 0, "distinct": 0, "cpu_total": 0.0, "cpu_verify": 0.0, "cpu_search": 0.0, "epochs": 0, "full_stored": 0}
    st.update(prev)
    seen_rs = set()
    if (folder / "rowspaces.txt").exists():
        seen_rs = set(x for x in (folder / "rowspaces.txt").read_text().split() if x)
    st["distinct"] = len(seen_rs)
    reconciled, seen_rs = reconcile(folder, st)                        # counters from the durable journals (item 6)
    cpu_prev_total, cpu_prev_verify = float(st["cpu_total"]), float(st["cpu_verify"]); cpu_verify_run = 0.0
    w.append("runs.jsonl", {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "run": run_no, "pid": os.getpid(), "seed": a.seed, "argv": sys.argv,
                            "python": platform.python_version(), "numpy": np.__version__, "reconciled": reconciled,
                            "sources": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES if p.exists()}})
    rng = np.random.default_rng(a.seed)
    S0 = seed_anti(); gens = q8_gens(a.k)
    t_end = time.time() + 60 * a.max_minutes if a.max_minutes else None
    epoch = int(st["epochs"]); full_stored = int(st["full_stored"])
    pool = None; seen = set(); since = 0; last_check = 0.0; last_save = 0.0

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
        if pool is None or since >= a.refresh:
            epoch += 1
            w.append("epochs.jsonl", {"epoch": epoch, "run": run_no, "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                                      "rng_state": rng.bit_generator.state, "pool": a.pool, "L": L})
            B = None
            while B is None:
                B = random_invariant_plane(rng)
            cells = pool = [random_anti_equivariant(rng, S0) for _ in range(a.pool)]
            tabs = [table(B, S) for S in cells]
            neg = tabs[0][1]
            exact = [T.astype(object) for T, _ in tabs]
            modp = [[np.mod(T, P).astype(np.float64) for T, _ in tabs] for P in PRIMES]
            n = tabs[0][0].shape[0]; seen = set(); since = 0
        key = tuple(int(x) for x in rng.integers(0, a.pool, L))
        st["draws"] += 1
        if key in seen:
            st["dup"] += 1
            continue
        seen.add(key); since += 1; st["tried"] += 1; w.state["attempts"] += 1
        cand = {"nega": True, "cyc": True}
        for P, pf in zip(PRIMES, modp):
            Mp = pf[key[0]]
            for i in key[1:]:
                Mp = np.fmod(Mp @ pf[i], P)
            tn = np.fmod(Mp[neg, np.arange(n)].sum(), P); tc = np.fmod(np.trace(Mp), P)
            cand["nega"] &= tn in (24.0, P - 24.0); cand["cyc"] &= tc in (24.0, P - 24.0)
            if not (cand["nega"] or cand["cyc"]):
                break
        if not (cand["nega"] or cand["cyc"]):
            if time.time() - last_save > 20:
                save(); last_save = time.time()
            continue
        st["screened"] += 1
        M = exact[key[0]]
        for i in key[1:]:
            M = M.dot(exact[i])
        vals = {"nega": trJ(M, neg), "cyc": int(np.trace(M))}
        any_hit = False
        for kind, v in vals.items():
            if abs(v) != 24:
                continue
            any_hit = True; st["hits"][kind] += 1; w.state["hits"] += 1
            S_cells = [cells[i] for i in key]
            t0 = time.process_time()
            try:
                rep, G = verify_hit(B, S_cells, kind == "cyc", a.k, gens)
                rep = signed_check(rep, v)
                sha = V.rowspace_sha256(G) if rep["ok"] else None
            except Exception as e:
                rep, G, sha = {"ok": False, "error": repr(e)}, None, None
            cpu_verify_run += time.process_time() - t0
            if rep.get("ok"):
                st["verified"] += 1
                if sha not in seen_rs:
                    seen_rs.add(sha); st["distinct"] += 1
                    with open(folder / "rowspaces.txt", "a") as f:
                        f.write(sha + "\n")
            else:
                st["verify_fail"] += 1
            w.append("hit_index.jsonl", {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "run": run_no, "epoch": epoch, "tuple": list(key),
                                         "kind": kind, "index": v, "rowspace_sha256": sha, "verified": bool(rep.get("ok"))})
            if full_stored < a.hits_max:
                full_stored += 1
                w.append("hits.jsonl", {"k": a.k, "kind": kind, "index": v, "epoch": epoch, "tuple": list(key), "B": B.tolist(),
                                        "cells": [S.tolist() for S in S_cells], "verification": rep, "rowspace_sha256": sha,
                                        "time": time.strftime("%Y-%m-%d %H:%M:%S")})
        if any_hit:
            st["cand_hit"] += 1; save(force=True); last_save = time.time()
        elif time.time() - last_save > 20:
            save(); last_save = time.time()
    save(force=True)
    w.finish("stopped" if stop_requested() else "finished (time limit)")


if __name__ == "__main__":
    main()
