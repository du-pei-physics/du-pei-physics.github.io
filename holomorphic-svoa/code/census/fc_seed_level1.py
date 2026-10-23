# Copied from claude/paper/fc_seed_level1.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Second invariant for every k = 6 seed: the first excited Ramond level of T = V_{Lambda(C)}//Z_6, and dim T_1.
20 Oct 2026; runs after fc_seed_census.py (same 7,112 seeds, same construction).

The Ramond character of T has integer weights above the ground level (field_content.md 3.2), so the first excited
level is L_0 = c/24 + 1 = 4 for every seed.  It is the sum of three parts:
  untwisted:  (1/6)[ Z^{(0,0)}_R + Z^{(0,2)}_R + Z^{(0,4)}_R ] at q^1.  Z^{(0,0)}_R = Theta_{Lambda+chi/2}/eta^72 gives
              N8 + 72 N, where N8 = #(norm-8 vectors of Lambda + chi/2) = 72 N + 8 A_69 (a full-support word with one
              coordinate lifted to |y| = 5/2, or a weight-69 word with three zero coordinates at y = +-3/2).
              The g^2-fixed terms vanish at q^1 (the fixed coset has no norm-8 vectors and 1/Delta(3 tau) no q^1 term).
              Total (144 N + 8 A_69)/6.
  g^2, g^4-twisted:  each (1/6) sum_b Z^{(2,b)}_R at q^1 = (1/2) Z^{(2,0)}_R[q^1] (b odd: 0; b = 2, 4: T-images,
              phase 1 at integer level).  Z^{(2,0)}_R = Theta_{L* + chi_L/2}/Delta(tau/3), so the two sectors give
              X = [Theta_{L*+chi_L/2} / Delta(tau/3)]_{q^1}, computed from the weight enumerator of Cbar.
  g, g^5-twisted:  built from Delta(3 tau)/Delta(6 tau) alone (no fixed momenta), hence the same for every seed:
              K = 864,288,768 (fc_characters.json, sectors 1 and 5: 432,144,384 each).  g^3-twisted: starts at L_0 = 9/2.
So D_1 = (144 N + 8 A_69)/6 + X(Cbar) + K, split evenly between the two parities (G_0 pairs them above the ground level).
dim T_1 = A_6/6: the 72 Heisenberg currents are projected out (g has no eigenvalue 1), Z_6 acts freely on the A_6
norm-2 vectors (L has no norm-2 vectors and g^3 = -1), and no twisted sector reaches weight 1.

A_69 and A_6 come from truncated transfer DPs over the generator (I_36 | A) (state s_j = a_j Q_R in F_3^r):
codewords with at most 3 zero coordinates, resp. weight at most 6.  The DP also returns A_72 = N (cross-check
with fc_seed_census.py).  Cbar's weight enumerator comes from all 3^12 words.  Exact integer arithmetic.
Spot check: the certificate seed must give A_69 = 637,024,272, A_6 = 0, D_1 = 2,579,935,232 (field_content.md).
usage: python fc_seed_level1.py --part i/m [--limit n]       (i = 1..m; m parallel workers, disjoint parts)
       python fc_seed_level1.py --summarize
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"                                         # single-thread BLAS (one worker = one core)
import argparse, itertools, json, sys, time
from collections import Counter
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fc_seed_census import records, ROOT, h, CONTRACT, sym2_counts   # noqa: E402
from count_k6_seed_codes import code, S                               # noqa: E402
from f3lib import fixed_subcode, rref, all_codewords                 # noqa: E402

K_G = 864288768
OUTDIR = ROOT / "out/census"
ALL12 = np.array(list(itertools.product(range(3), repeat=12)), dtype=np.int64)    # index sum W_j 3^(11-j)
PW = 3 ** np.arange(11, -1, -1)
OUT_IMG = (ALL12 @ S) % 3                                         # c = W S for every W (S = Paley, fixed)
ZEROS_OUT = 12 - np.count_nonzero(OUT_IMG, axis=1)
WT_OUT = np.count_nonzero(OUT_IMG, axis=1)
NZ12 = np.count_nonzero(ALL12, axis=1)
INPUTS = {"zeros": (ALL12[NZ12 >= 9], 12 - NZ12[NZ12 >= 9], ZEROS_OUT, 3),
          "weight": (ALL12[NZ12 <= 6], NZ12[NZ12 <= 6], WT_OUT, 6)}


def graded_trace(Q, R, kind):
    """[x^k] sum over codewords of x^{stat}, k <= cap, stat = #zeros or weight; period j contributes
    stat(a_j) + stat(c_j), c_j = W S with W_K = (a_j Q)_K, W_R = (a_{j-1} Q)_R (negated across the wrap)."""
    V, st_in, st_out, cap = INPUTS[kind]
    R = sorted(R); r = len(R); n = 3 ** r; rest = [j for j in range(12) if j not in R]
    y = (V @ Q) % 3
    sidx = y[:, R] @ (3 ** np.arange(r))
    tail = y[:, rest] @ PW[rest]
    states = np.array(list(itertools.product(range(3), repeat=r)), dtype=np.int64)[:, ::-1]   # idx = sum p_k 3^k
    assert np.array_equal(states @ (3 ** np.arange(r)), np.arange(n))
    off = states @ PW[R]
    neg = ((-states) % 3) @ (3 ** np.arange(r))
    M = np.zeros((cap + 1, n, n), dtype=np.int64)                 # M[k, p, s]
    for p in range(n):
        k = st_in + st_out[tail + off[p]]
        ok = k <= cap
        M[:, p, :] = np.bincount(k[ok] * n + sidx[ok], minlength=(cap + 1) * n).reshape(cap + 1, n)
    Mf = M.astype(np.float64)
    nzk = [k for k in range(cap + 1) if M[k].any()]
    out = [0] * (cap + 1)
    for k0 in nzk:
        for k1 in nzk:
            if k0 + k1 > cap:
                continue
            P = Mf[k0] @ Mf[k1]
            assert P.max() < 2 ** 52                               # exact in float64
            P = np.rint(P).astype(np.int64)[neg, :]
            for k2 in nzk:
                if k0 + k1 + k2 <= cap:
                    out[k0 + k1 + k2] += int((P * M[k2].T).sum())
    return out


def theta_twisted_level1(cbar_enum):
    """[Theta_{L*+chi_L/2} / Delta(tau/3)] at q^1.  Vectors (b + 3/2)/3 per coordinate, b = c mod 3; norm
    sum (2b+3)^2/36; q-exponent norm/2 = sum (2b+3)^2 / 72.  Work in units u = (2b+3)^2 (q^{u/72})."""
    cap = 96                                                     # norm <= 8/3  <->  sum u <= 96
    def coord(res):
        d = Counter()
        for b in range(-20, 21):
            if b % 3 == res:
                u = (2 * b + 3) ** 2
                if u <= cap:
                    d[u] += 1
        return d
    f0, f1 = coord(0), coord(1)
    def mul(a, b):
        out = Counter()
        for i, x in a.items():
            for j, y in b.items():
                if i + j <= cap:
                    out[i + j] += x * y
        return out
    p0 = [Counter({0: 1})]; p1 = [Counter({0: 1})]
    for _ in range(24):
        p0.append(mul(p0[-1], f0)); p1.append(mul(p1[-1], f1))
    th = Counter()
    for w, A in enumerate(cbar_enum):
        if A:
            for u, c in mul(p0[24 - w], p1[w]).items():
                th[u] += A * c
    # 1/Delta(tau/3) = q^{-1/3} sum_n p24(n) q^{n/3};  q^1  <->  u/72 - 1/3 + n/3 = 1  <->  u = 96 - 24 n
    p24 = [1] + [0] * 4
    for m in range(1, 5):
        for _ in range(24):
            for i in range(m, 5):
                p24[i] += p24[i - m]
    assert p24[:4] == [1, 24, 324, 3200]
    return sum(th.get(96 - 24 * n, 0) * p24[n] for n in range(5)), {str(u): th[u] for u in sorted(th)}


def analyse(Q, R):
    z = graded_trace(Q, R, "zeros")
    wt = graded_trace(Q, R, "weight")
    N, A69 = z[0], z[3]
    assert z[1] == z[2] == 0 and wt[0] == 1 and wt[1:6] == [0] * 5
    A6 = wt[6]
    G = code(Q, R)
    Fx = fixed_subcode(G, h); Fc = rref(Fx[:, CONTRACT] % 3)[0]
    cb = np.zeros(25, dtype=np.int64)
    for w in all_codewords(Fc):
        np.add.at(cb, np.count_nonzero(w, axis=1), 1)
    cb = [int(x) for x in cb]
    X, th = theta_twisted_level1(cb)
    assert (144 * N + 8 * A69) % 6 == 0 and A6 % 6 == 0
    D1 = (144 * N + 8 * A69) // 6 + X + K_G
    assert D1 % 2 == 0
    return {"N": N, "A69": A69, "A6": A6, "dim_T1": A6 // 6, "Cbar_enum": cb, "X_twisted": X, "D1": D1,
            "level1_even_odd": [D1 // 2, D1 // 2]}


def run(part, m, limit):
    assert 1 <= part <= m, "parts are 1-based: --part 1/2, --part 2/2"
    recs = records()[part - 1::m]
    out = OUTDIR / f"fc_seed_level1_part{part}of{m}.jsonl"
    done = {tuple(json.loads(l)["key"]) for l in open(out)} if out.exists() else set()
    t0 = time.process_time(); k = 0
    with open(out, "a") as f:
        for key, Q, R, rec in recs:
            if key in done:
                continue
            if limit is not None and k >= limit:
                break
            res = analyse(Q, R)
            assert res["N"] == rec["nmax"], (key, res["N"], rec["nmax"])
            res["key"] = list(key)
            f.write(json.dumps(res) + "\n"); f.flush(); k += 1
    print(f"part {part}/{m}: {k} records, {time.process_time() - t0:.1f} CPU s")


def summarize():
    first = {tuple(json.loads(l)["key"]): json.loads(l) for l in open(OUTDIR / "fc_seed_census.jsonl")}
    sec = {}
    for p in sorted(OUTDIR.glob("fc_seed_level1_part*.jsonl")):
        for l in open(p):
            r = json.loads(l); sec[tuple(r["key"])] = r
    rows = []
    for k, r in first.items():
        if k in sec:
            s = sec[k]
            assert s["N"] == r["N"]
            rows.append({"key": k, "census": r["census"], "class_key": r["class_key"], "n": r["n_total"],
                         "D1": s["D1"], "T1": s["dim_T1"]})
    inv2 = Counter((x["n"], x["D1"]) for x in rows)
    inv3 = Counter((x["n"], x["D1"], x["T1"]) for x in rows)
    classes = {}
    for x in rows:
        classes.setdefault(x["class_key"], set()).add((x["n"], x["D1"], x["T1"]))
    bad = [k for k, v in classes.items() if len(v) > 1]
    hp_classes = {x["class_key"] for x in rows if x["census"] == "hp"}
    # T', T'' Ramond ground counts depend on n alone; their first excited levels too would need more; we report
    # whether the pair (T'' ground counts, D1) separates as (n, D1) does, which holds iff n -> T'' counts is injective.
    tpp = {}
    for v in {x["n"] for x in rows}:
        B, F = (v - 4) // 2, (v + 4) // 2
        Bp, Fp = sym2_counts(B, F); tpp[v] = sym2_counts(Fp, Bp)
    s = {"records_with_both": len(rows), "distinct_n": len({x['n'] for x in rows}),
         "distinct_n_D1": len(inv2), "distinct_n_D1_T1": len(inv3),
         "singletons_n_D1_T1": sum(1 for v in inv3.values() if v == 1),
         "largest_tie_group": max(inv3.values()), "tie_group_size_histogram": {str(k): v for k, v in sorted(Counter(inv3.values()).items())},
         "distinct_D1": len({x["D1"] for x in rows}), "distinct_T1": sorted({x["T1"] for x in rows}),
         "histogram_T1": {str(k): v for k, v in sorted(Counter(x["T1"] for x in rows).items())},
         "classes": len(classes), "hp_classes": len(hp_classes),
         "classes_where_invariant_varies": len(bad), "examples_varying": bad[:3],
         "distinct_invariant_values_per_class_max": max(len(v) for v in classes.values()),
         "classes_merged_by_invariant": len(classes) - len(inv3),
         "n_to_Tpp_injective": len(set(tpp.values())) == len(tpp)}
    (OUTDIR / "fc_seed_level1_summary.json").write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="1/1")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--summarize", action="store_true")
    a = ap.parse_args()
    if a.summarize:
        summarize()
    else:
        i, m = map(int, a.part.split("/"))
        run(i, m, a.limit)
