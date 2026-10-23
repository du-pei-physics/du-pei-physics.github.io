# Copied from claude/paper/fc_seed_level2.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Level 2 of the k = 6 seed census: the complete NS and Ramond characters of T = V_{Lambda(C)}//Z_6 for every seed.
21 Oct 2026; runs after fc_seed_census.py and fc_seed_level1.py (same 7,112 seeds, same construction).

Why this is the last character-level step.  By field_content.md section 6, Z^NS_T = sum_{j=0}^{6} a_j t^j with
t = (theta_3/eta)^12, a_6 = 1, and Z^R_T = sum_j (-1)^j a_j u^j, u = (theta_2/eta)^12.  dim T_{1/2} = 0 forces
a_5 = -144.  The remaining a_4, ..., a_1 are fixed triangularly by dim T_1, T_{3/2}, T_2, T_{5/2}, and a_0 is the number
n of Ramond ground states.  So (n, T_1, T_{3/2}, T_2, T_{5/2}) determines every NS weight space and every Ramond level
of T; higher Ramond levels carry no further information.  The level-1 degeneracy D1 (fc_seed_level1.py), computed
independently from the twisted sectors, must equal -4096 a_1: a per-seed check of the whole orbifold computation.

NS weight spaces below 8/3 come from the untwisted sector only (twisted ground energies 8/3, 19/6, 9/2).  With
N_m = #(norm-m vectors of Lambda), L = Lambda^{g^2} (norms sum b^2, b mod 3 in Cbar; no vectors of norm 1, 2, 4, 5),
and the twining traces of g^j on V (field_content.md section 2; only g^2, g^4 fix momenta, with trivial phase (A1)):
  dim T_1     = A_6 / 6
  dim T_{3/2} = (N_3 + 2 A_3(Cbar)) / 6,                N_3 = 144 + A_9 + 6 A_6
  dim T_2     = (2700 + 2556 + N_4 + 72 N_2) / 6,       N_2 = A_6,  N_4 = A_12 + 9 A_9 + 15 A_6
                (2700 = Heisenberg weight-2 states, 2556 = their g^3-trace; the other g^j traces vanish at q^1, q^2)
  dim T_{5/2} = (N_5 + 72 N_3) / 6,                     N_5 = A_15 + 12 A_12 + 36 A_9 + 152 A_6
A_w (w <= 15) come from a weight-graded transfer trace over the generator (cap 15).  The derived invariants (the
polynomial a_j, the NS dimensions T_3, T_{7/2}, T_4 and the Ramond levels D_2, D_3, D_4) are written out as well.
Certificate seed: a = [13940, -629867, 751404, -129679, 6984, -144, 1] (field_content.md 6.2).
usage: python fc_seed_level2.py --part i/m [--limit n]      (i = 1..m)
       python fc_seed_level2.py --summarize
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import argparse, json, sys, time
from collections import Counter
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fc_seed_level1 as L1                                          # noqa: E402
from fc_seed_census import records, ROOT, sym2_counts                # noqa: E402
import fc_closed_form as CF                                          # noqa: E402

OUTDIR = ROOT / "out/census"
L1.INPUTS["weight15"] = (L1.ALL12, L1.NZ12, L1.WT_OUT, 15)
CERT_A = [13940, -629867, 751404, -129679, 6984, -144, 1]


def polynomial(n, T1, T32, T2, T52):
    """a_6..a_1 from the polar part (exponents -6..-1 in q^{1/2}), a_0 = n."""
    Z = {-6: 1, -5: 0, -4: T1, -3: T32, -2: T2, -1: T52}
    rest = dict(Z); a = [0] * 7
    for j in range(6, 0, -1):
        a[j] = rest.get(-j, 0)
        rest = CF.add(rest, CF.TP[j], coef=[1, -a[j]])
    a[0] = n
    return a


def analyse(Q, R, n, cbar_enum):
    A = L1.graded_trace(Q, R, "weight15")
    assert A[0] == 1 and all(A[w] == 0 for w in range(1, 16) if w % 3) and A[3] == 0
    A6, A9, A12, A15 = A[6], A[9], A[12], A[15]
    N3 = 144 + A9 + 6 * A6; N4 = A12 + 9 * A9 + 15 * A6; N5 = A15 + 12 * A12 + 36 * A9 + 152 * A6
    vals = {"T1": (A6, 6), "T32": (N3 + 2 * cbar_enum[3], 6), "T2": (5256 + N4 + 72 * A6, 6), "T52": (N5 + 72 * N3, 6)}
    for k, (x, d) in vals.items():
        assert x % d == 0, (k, x)
    T1, T32, T2, T52 = (vals[k][0] // 6 for k in ("T1", "T32", "T2", "T52"))
    a = polynomial(n, T1, T32, T2, T52)
    NS = CF.eval_poly(a); RR = CF.eval_poly(a, "R")
    out = {"A": {str(w): A[w] for w in (6, 9, 12, 15)}, "A3_Cbar": cbar_enum[3], "T1": T1, "T32": T32, "T2": T2,
           "T52": T52, "a": a, "D1_pred": -4096 * a[1],
           "T3": NS.get(0, 0), "T72": NS.get(1, 0), "T4": NS.get(2, 0),
           "D2": RR.get(4, 0), "D3": RR.get(6, 0), "D4": RR.get(8, 0)}
    assert all(x >= 0 for x in (out["T3"], out["T72"], out["T4"], out["D2"], out["D3"], out["D4"]))
    return out


def run(part, m, limit):
    assert 1 <= part <= m
    first = {tuple(json.loads(l)["key"]): json.loads(l) for l in open(OUTDIR / "fc_seed_census.jsonl")}
    lev1 = {}
    for p in OUTDIR.glob("fc_seed_level1_part*of2.jsonl"):
        for l in open(p):
            r = json.loads(l); lev1[tuple(r["key"])] = r
    recs = records()[part - 1::m]
    out = OUTDIR / f"fc_seed_level2_part{part}of{m}.jsonl"
    done = {tuple(json.loads(l)["key"]) for l in open(out)} if out.exists() else set()
    t0 = time.process_time(); k = 0
    with open(out, "a") as f:
        for key, Q, R, rec in recs:
            if key in done:
                continue
            if limit is not None and k >= limit:
                break
            res = analyse(Q, R, first[key]["n_total"], lev1[key]["Cbar_enum"])
            assert res["T1"] == lev1[key]["dim_T1"], key
            assert res["D1_pred"] == lev1[key]["D1"], ("D1 check failed", key, res["D1_pred"], lev1[key]["D1"])
            res["key"] = list(key)
            f.write(json.dumps(res) + "\n"); f.flush(); k += 1
    print(f"part {part}/{m}: {k} records, {time.process_time() - t0:.1f} CPU s")


def summarize():
    first = {tuple(json.loads(l)["key"]): json.loads(l) for l in open(OUTDIR / "fc_seed_census.jsonl")}
    lev2 = {}
    for p in OUTDIR.glob("fc_seed_level2_part*of2.jsonl"):
        for l in open(p):
            r = json.loads(l); lev2[tuple(r["key"])] = r
    rows = [(k, first[k], lev2[k]) for k in first if k in lev2]
    inv = {k: tuple(s["a"]) for k, f, s in rows}                  # the full character = (n, T1, T32, T2, T52)
    inv_old = {k: (f["n_total"], s["D1_pred"], s["T1"]) for k, f, s in rows}
    cls = {}
    for k, f, s in rows:
        cls.setdefault(f["class_key"], set()).add(inv[k])
    bad = [c for c, v in cls.items() if len(v) > 1]
    by_val = {}
    for k, f, s in rows:
        by_val.setdefault(inv[k], set()).add(f["class_key"])
    per = Counter(len(v) for v in by_val.values())
    hp = [(k, f) for k, f, s in rows if f["census"] == "hp"]
    s = {"records": len(rows), "distinct_full_character": len(set(inv.values())),
         "distinct_previous_invariant": len(set(inv_old.values())),
         "classes": len(cls), "classes_where_invariant_varies": len(bad), "examples": bad[:3],
         "values_with_one_class": per[1], "values_shared_by_several_classes": sum(v for c, v in per.items() if c > 1),
         "largest_group_by_classes": max(per), "classes_per_value_histogram": {str(c): v for c, v in sorted(per.items())},
         "largest_group_by_seeds": max(Counter(inv.values()).values()),
         "hp_seeds": len(hp), "hp_distinct_values": len({inv[k] for k, f in hp}),
         "hp_classes": len({f["class_key"] for k, f in hp}),
         "refl_distinct_values": len({inv[k] for k, f, s in rows if f["census"] == "refl"}),
         "hp_refl_common_values": len({inv[k] for k, f in hp} & {inv[k] for k, f, s in rows if f["census"] == "refl"}),
         "histogram_T32": {str(x): c for x, c in sorted(Counter(s_["T32"] for _, _, s_ in rows).items())},
         "A3_Cbar_values": dict(Counter(s_["A3_Cbar"] for _, _, s_ in rows))}
    (OUTDIR / "fc_seed_level2_summary.json").write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps({k: v for k, v in s.items() if k != "histogram_T32"}, indent=1))


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
