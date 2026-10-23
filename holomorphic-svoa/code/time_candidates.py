"""Measured CPU cost of the basic search steps (process CPU time, one core, exact integer arithmetic).

    python3 time_candidates.py > ../notes/timings.json        (a few seconds)

Times (1) building the transfer matrix of one random cell (done once per cell of a pool and amortised over all
candidates drawn from that pool) and (2) evaluating one candidate exactly, i.e. multiplying the transfer matrices
of one random choice of cells and taking the traces, for:
  r = 3 boundary spins (27 x 27 matrices): a k = 13 chain of six width-12 cells and one width-6 cell;
  r = 4 boundary spins (81 x 81 matrices, discriminant 2, two reflections on S): a k = 21 period (12,12,12,6),
      cubed for the three periods, with both closures and the balance traces.
Not included: the fast modular screen of mixed_rings.py (used by the k = 18, 21 searches) and the verification of
hits from the expanded generator.
"""
import json, platform, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cells_general import anti_bases          # noqa: E402
from mixed_products import random_cell, trJ   # noqa: E402


def build_pool(rng, q, r, n, metric2=0, s_refl=0):
    bases = anti_bases(q)
    t0 = time.process_time()
    pool = [random_cell(rng, q, r, bases, metric2, s_refl) for _ in range(n)]
    return pool, (time.process_time() - t0) / n


def exact_eval(rng, pools, widths, periods, r, n_cand):
    neg = pools[widths[0]][0][1]
    n = 3 ** r
    t0 = time.process_time()
    for _ in range(n_cand):
        M = np.eye(n, dtype=object)
        for q in widths:
            M = M.dot(pools[q][int(rng.integers(len(pools[q])))][0])
        if periods == 1:
            trJ(M, neg); int(np.trace(M))
        else:
            M3 = M.dot(M).dot(M)
            trJ(M, neg); trJ(M3, neg); int(np.trace(M)); int(np.trace(M3))
    return (time.process_time() - t0) / n_cand


def main():
    out = {"python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine(),
           "system": platform.system(), "unit": "CPU seconds"}
    rng = np.random.default_rng(1)
    p12, c12 = build_pool(rng, 12, 3, 20)
    p6, c6 = build_pool(rng, 6, 3, 20)
    out["r3_build_transfer_matrix_width12"] = round(c12, 5)
    out["r3_build_transfer_matrix_width6"] = round(c6, 5)
    out["r3_k13_evaluate_one_candidate"] = round(exact_eval(rng, {12: p12, 6: p6}, [12] * 6 + [6], 1, 3, 300), 5)
    q12, d12 = build_pool(rng, 12, 4, 12, metric2=1, s_refl=2)
    q6, d6 = build_pool(rng, 6, 4, 12, metric2=1, s_refl=2)
    out["r4_build_transfer_matrix_width12"] = round(d12, 5)
    out["r4_build_transfer_matrix_width6"] = round(d6, 5)
    out["r4_k21_evaluate_one_candidate_exact"] = round(exact_eval(rng, {12: q12, 6: q6}, [12, 12, 12, 6], 3, 4, 40), 5)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
