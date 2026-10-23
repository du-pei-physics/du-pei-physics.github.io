#!/usr/bin/env python3
# Copied from claude/census_general/validate_on_census.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Validate k6_general.analyse on census seeds against the published polynomials (k6_seed_polynomials.csv)."""
import csv, json, random, sys, time
from pathlib import Path
import k6_general as KG
sys.path.insert(0, str(KG.CODE / "census"))
from fc_seed_census import records                       # noqa: E402
from count_k6_seed_codes import code, rowspace_sha256    # noqa: E402

DATA = KG.CODE.parent / "data" / "k6_seed_census" / "k6_seed_polynomials.csv"
rows = [r for r in csv.reader(l for l in open(DATA) if not l.startswith("#"))]
hdr, rows = rows[0], rows[1:]
poly = {r[hdr.index("rowspace_sha256")]: [int(r[hdr.index(f"a{j}")]) for j in range(7)] for r in rows}
recs = records()
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
sample = random.sample(range(len(recs)), int(sys.argv[2]) if len(sys.argv) > 2 else 8)
ok = 0
for i in sample:
    key, Q, R, _ = recs[i]
    G = code(Q, R); sha = rowspace_sha256(G)
    t = time.process_time(); r = KG.analyse(G, "nega"); dt = time.process_time() - t
    good = r["ok"] and r["a"] == poly[sha]
    ok += good
    print(key[0], i, "ok" if good else "MISMATCH", r.get("a"), poly[sha] if not good else "", f"{dt:.1f}s")
print(f"{ok}/{len(sample)} agree")
