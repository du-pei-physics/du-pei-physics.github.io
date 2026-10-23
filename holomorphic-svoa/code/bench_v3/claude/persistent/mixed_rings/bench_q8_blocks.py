"""Benchmark-v3 worker for the Q8-equivariant chains with the batched block engine of the k = 16 worker
(claude/persistent/q8_k16/run.py, claude/k24_transfer/q8_rings.py), 19 Oct 2026.

bench_q8.py multiplies full 81 x 81 transfer matrices one candidate at a time; at k = 16 (chains of eight cells) that
runs about 1,000 candidates per second, a third of the old q8_k16 worker, which (i) splits every transfer table
into the blocks of its Q8-symmetry decomposition (core.blocks) and (ii) screens 20,000 candidates per batch with
stacked modular products (kspan_rings.tr_mod, two primes below 2^20).  This worker uses that engine with the
benchmark-v3 accounting of bench_worker.py / bench_q8.py: tuples are de-duplicated within a pool, each candidate is
tested in both closures, every exact hit is verified with the online note's stand-alone verifier (bench_q8.verify_hit:
rank, self-orthogonality, banded transfer DP, check_q8), distinct RREF rowspaces go to rowspaces.txt, CPU of search
and verification are recorded separately, and the RNG state before each pool is logged in epochs.jsonl.
usage: python bench_q8_blocks.py --name bench3_q8b_k16 --k 16 --seed 1 [--max-minutes 360] [--global-stop FILE]
"""
import argparse, hashlib, json, os, platform, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
CL = HERE.parent.parent
for p in (HERE.parent, CL / "k8", CL, CL / "k24_transfer", CL / "online/holomorphic-svoa/code", HERE):
    sys.path.insert(0, str(p))
from common import Worker
import verify as V
import q8_rings as qr
import kspan_rings as kr
from search_p4 import PR
from bench_q8 import verify_hit, q8_gens, signed_check
from journal import reconcile

SOURCES = [HERE / "bench_q8_blocks.py", HERE / "bench_q8.py", HERE / "journal.py", CL / "k24_transfer/q8_rings.py", CL / "k24_transfer/kspan_rings.py",
           CL / "k8/q8_cells.py", CL / "online/holomorphic-svoa/code/verify.py"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True); ap.add_argument("--k", type=int, choices=(8, 16), required=True)
    ap.add_argument("--pool", type=int, default=150); ap.add_argument("--batch", type=int, default=20000)
    ap.add_argument("--batches-per-pool", type=int, default=6)
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
                            "python": platform.python_version(), "numpy": np.__version__, "engine": "blocks", "reconciled": reconciled,
                            "sources": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES if p.exists()}})
    rng = np.random.default_rng(a.seed)
    gens = q8_gens(a.k)
    t_end = time.time() + 60 * a.max_minutes if a.max_minutes else None
    epoch = int(st["epochs"]); full_stored = int(st["full_stored"])

    def stop_requested():
        return w.should_stop() or (a.global_stop is not None and os.path.exists(a.global_stop))

    def save(force=False):
        st["cpu_total"] = round(cpu_prev_total + time.process_time(), 2)
        st["cpu_verify"] = round(cpu_prev_verify + cpu_verify_run, 2)
        st["cpu_search"] = round(st["cpu_total"] - st["cpu_verify"], 2)
        st["epochs"] = epoch; st["full_stored"] = full_stored
        w.save(force=force, stats=st)

    while not stop_requested() and (t_end is None or time.time() < t_end):
        w.wait_for_disk()
        rng_state = rng.bit_generator.state
        P = qr.q8_pool(rng, a.pool)
        if P is None:
            continue
        epoch += 1
        w.append("epochs.jsonl", {"epoch": epoch, "run": run_no, "time": time.strftime("%Y-%m-%d %H:%M:%S"), "rng_state": rng_state,
                                  "pool": a.pool, "pool_cells": len(P["specs"]), "L": L, "batch": a.batch})
        N = len(P["specs"]); seen = set()
        for _ in range(a.batches_per_pool):
            if stop_requested() or (t_end is not None and time.time() >= t_end):
                break
            idx = rng.integers(0, N, (a.batch, L))
            keep = np.zeros(len(idx), dtype=bool)
            for r, row in enumerate(idx):
                key = row.tobytes()
                if key not in seen:
                    seen.add(key); keep[r] = True
            st["draws"] += len(idx); st["dup"] += int((~keep).sum())
            idx = idx[keep]; st["tried"] += len(idx); w.state["attempts"] += len(idx)
            hit_rows = set()
            for tw in ("nega", "cyc"):
                signs = P["signs"] if tw == "nega" else [1] * len(P["signs"])
                v0 = kr.tr_mod(P["stacks"], signs, idx, 1, PR[0]); cand = idx[(v0 == 24) | (v0 == PR[0] - 24)]
                if len(cand):
                    v1 = kr.tr_mod(P["stacks"], signs, cand, 1, PR[1]); cand = cand[(v1 == 24) | (v1 == PR[1] - 24)]
                st["screened"] += len(cand)
                for row in cand:
                    v = kr.tr_exact(P["stacks"], signs, row, 1)
                    if abs(v) != 24:
                        continue
                    st["hits"][tw] += 1; w.state["hits"] += 1; hit_rows.add(row.tobytes())
                    cells = [P["specs"][int(c)] for c in row]
                    t0 = time.process_time()
                    try:
                        rep, G = verify_hit(P["B"], cells, tw == "cyc", a.k, gens)
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
                    w.append("hit_index.jsonl", {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "run": run_no, "epoch": epoch,
                                                 "tuple": [int(c) for c in row], "kind": tw, "index": int(v), "rowspace_sha256": sha,
                                                 "verified": bool(rep.get("ok"))})
                    if full_stored < a.hits_max:
                        full_stored += 1
                        w.append("hits.jsonl", {"k": a.k, "kind": tw, "index": int(v), "epoch": epoch, "tuple": [int(c) for c in row],
                                                "B": P["B"].tolist(), "cells": [S.tolist() for S in cells], "verification": rep,
                                                "rowspace_sha256": sha, "time": time.strftime("%Y-%m-%d %H:%M:%S")})
                    save(force=True)
            st["cand_hit"] += len(hit_rows)
            save()
    save(force=True)
    w.finish("stopped" if stop_requested() else "finished (time limit)")


if __name__ == "__main__":
    main()
