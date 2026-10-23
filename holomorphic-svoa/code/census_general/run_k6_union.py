#!/usr/bin/env python3
"""Recompute the characters of the 1,636 non-census k = 6 codes of ../../data/k6_union/ from their shipped generators.

Written for the holomorphic-svoa note (not a workspace copy).  For each generator in
../../data/k6_union/k6_others_generators.jsonl.gz it checks the systematic form (I_36 | A) and the recorded rowspace
hash, runs k6_general.analyse (self-duality, the symmetry g of the stated closure: nega = signed shift, cyc = -shift,
g^3 = -1; the per-code hypotheses on the g^2-fixed code; the weight counts; the character polynomial; the twisted-sector
check D_1 = -4096 a_1) and compares the result with the shipped record in k6_others_records.jsonl.gz.
Outputs: ../out/census_general/k6_union_rerun.jsonl.  About 2 s per code on one core.
usage: python3 run_k6_union.py [--limit n] [--source bench3_k6|bench2_k6|bench_k6] [--closure nega|cyc]
"""
import argparse, gzip, hashlib, json, time
from pathlib import Path
import k6_general as KG

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "data" / "k6_union"
OUT = HERE.parent / "out" / "census_general" / "k6_union_rerun.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=3)
    ap.add_argument("--source"); ap.add_argument("--closure")
    a = ap.parse_args()
    ship = {json.loads(l)["sha"]: json.loads(l) for l in gzip.open(DATA / "k6_others_records.jsonl.gz", "rt")}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    n = agree = 0
    with open(OUT, "w") as f:
        for l in gzip.open(DATA / "k6_others_generators.jsonl.gz", "rt"):
            g = json.loads(l)
            if (a.source and g["source"] != a.source) or (a.closure and g["closure"] != a.closure):
                continue
            if n >= a.limit:
                break
            assert hashlib.sha256("\n".join(g["generator"]).encode()).hexdigest() == g["sha"]
            t0 = time.process_time()
            rec = KG.analyse([[int(c) for c in row] for row in g["generator"]], g["closure"])
            rec.update({"sha": g["sha"], "source": g["source"]})
            same = rec == ship[g["sha"]]
            n += 1; agree += same
            f.write(json.dumps(rec) + "\n")
            print(g["source"], g["closure"], g["sha"][:12], "a =", rec.get("a"), "identical to shipped record:", same,
                  f"({time.process_time() - t0:.1f} CPU s)")
    print(f"{agree}/{n} identical")


if __name__ == "__main__":
    main()
