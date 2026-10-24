#!/usr/bin/env python3
"""Recompute the characters of the 1,636 non-census k = 6 codes of ../../data/k6_union/ from their shipped generators.

Written for the holomorphic-svoa note (not a workspace copy).  For each selected generator in
../../data/k6_union/k6_others_generators.jsonl.gz it checks the systematic form (I_36 | A) and the recorded rowspace
hash, runs k6_general.analyse (self-duality, the symmetry g of the stated closure: nega = signed shift, cyc = -shift,
g^3 = -1; the per-code hypotheses on the g^2-fixed code; the weight counts; the character polynomial; the twisted-sector
check D_1 = -4096 a_1) and compares the result with the shipped record in k6_others_records.jsonl.gz.
Exit status 0 only if every selected code reproduces its shipped record exactly; 1 on any mismatch, on an empty
selection, or if the negative control fails.  --negative-control also compares the first recomputed record with a
deliberately altered copy of its shipped record (a_0 + 1) and requires that comparison to fail.
Output: ../out/census_general/k6_union_rerun.jsonl.  About 2 s per code on one core.
usage: python3 run_k6_union.py [--limit n] [--source bench3_k6|bench2_k6|bench_k6] [--closure nega|cyc] [--negative-control]
"""
import argparse, copy, gzip, hashlib, json, sys, time
from pathlib import Path
import k6_general as KG

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "data" / "k6_union"
OUT = HERE.parent / "out" / "census_general" / "k6_union_rerun.jsonl"


def positive(s):
    v = int(s)
    if v < 1:
        raise argparse.ArgumentTypeError("--limit must be at least 1")
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=positive, default=3)
    ap.add_argument("--source", choices=["bench3_k6", "bench2_k6", "bench_k6"])
    ap.add_argument("--closure", choices=["nega", "cyc"])
    ap.add_argument("--negative-control", action="store_true")
    a = ap.parse_args()
    ship = {json.loads(l)["sha"]: json.loads(l) for l in gzip.open(DATA / "k6_others_records.jsonl.gz", "rt")}
    sel = []
    for l in gzip.open(DATA / "k6_others_generators.jsonl.gz", "rt"):
        g = json.loads(l)
        if (a.source and g["source"] != a.source) or (a.closure and g["closure"] != a.closure):
            continue
        sel.append(g)
        if len(sel) >= a.limit:
            break
    if not sel:
        print("no code matches the selection", file=sys.stderr)
        return 1
    OUT.parent.mkdir(parents=True, exist_ok=True)
    agree, control_ok = 0, None
    with open(OUT, "w") as f:
        for i, g in enumerate(sel):
            if hashlib.sha256("\n".join(g["generator"]).encode()).hexdigest() != g["sha"]:
                print("hash mismatch for", g["sha"], file=sys.stderr)
                return 1
            t0 = time.process_time()
            rec = KG.analyse([[int(c) for c in row] for row in g["generator"]], g["closure"])
            rec.update({"sha": g["sha"], "source": g["source"]})
            same = rec == ship[g["sha"]]
            agree += same
            f.write(json.dumps(rec) + "\n")
            print(g["source"], g["closure"], g["sha"][:12], "a =", rec.get("a"), "identical to shipped record:", same,
                  f"({time.process_time() - t0:.1f} CPU s)")
            if a.negative_control and i == 0:
                altered = copy.deepcopy(ship[g["sha"]])
                altered["a"][0] += 1
                control_ok = rec != altered
                print("negative control (altered a_0 must not match):", "detected" if control_ok else "NOT DETECTED")
    print(f"{agree}/{len(sel)} identical")
    if agree != len(sel) or control_ok is False:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
