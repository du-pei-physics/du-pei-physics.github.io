# Copied from claude/census/count_k6_seed_codes.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Distinct k = 6 seed codes among the 7,112 Hadamard-Paley seed cells of the census (20 Oct 2026).

The paper's yield table counts 7,112 k = 6 seed cells: the distinct records of persistent/hp_census/hits.jsonl
(reflection-free Hadamard-Paley cells, Q = matrix_Q(blocks, parities); 1,081 lines, 1,078 distinct records) and of
persistent/refl_census/hits.jsonl (Q = one or two reflections R_v R_w on the Paley base; 6,034 records), each with a
delayed (boundary) set R, S = the Paley matrix, I1 = 0 and |I3| = 24.  This script builds the length-72 code of every
record exactly as paper/canonical_seed_data.py does (generator (I_36 | A), A = 3x3 block negacyclic matrix of
D = Q[:, K] S[K, :] and E = Q[:, R] S[R, :]), checks that it is self-dual, checks the construction against the k = 6
certificate (same rowspace hash) and against the recorded index on a sample (banded transfer DP over the generator),
and counts distinct codes as distinct rowspaces: (I | A) is already in reduced row echelon form, so the rowspace hash is
the SHA-256 of its rows over {0,1,2} joined by newlines (the online note's definition).
Exact integer / F_3 arithmetic; one process, well under a minute.
usage: python count_k6_seed_codes.py [--sample 60] [--out ../out/census/k6_seed_codes.json]
"""
import argparse, hashlib, json, random, sys
from collections import Counter
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT))
from hpcell import matrix_Q, PALEY                      # noqa: E402
from verify_mixed_hit import banded_index               # noqa: E402

S = PALEY % 3


def reflection(v):                                      # as in persistent/refl_census/run.py
    v = np.asarray(v, dtype=np.int64) % 3
    nn = int(v @ v) % 3
    assert nn
    return (np.eye(12, dtype=np.int64) - 2 * nn * np.outer(v, v)) % 3


def code(Q, R):
    """generator (I_36 | A) of the anti-periodic chain of three copies of the cell (Q, S, R)"""
    R = sorted(R)
    K = [i for i in range(12) if i not in R]
    D = (Q[:, K] @ S[K, :]) % 3
    E = (Q[:, R] @ S[R, :]) % 3
    Z = np.zeros((12, 12), dtype=np.int64)
    A = np.block([[D, E, Z], [Z, D, E], [(-E) % 3, Z, D]]) % 3
    return np.hstack([np.eye(36, dtype=np.int64), A])


def rowspace_sha256(G):                                  # G = (I | A) is in reduced row echelon form already
    return hashlib.sha256("\n".join("".join(str(int(x)) for x in row) for row in G).encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=60, help="records whose index is recomputed from the generator")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    recs = []
    for line in open(ROOT.parent / "data" / "k6_seed_census" / "hp_census_hits.jsonl"):
        h = json.loads(line)
        recs.append(("hp", (json.dumps(h["blocks"]), json.dumps(h["parities"]), tuple(h["R"])),
                     matrix_Q(h["blocks"], h["parities"]) % 3, h["R"], h["I3"]))
    for line in open(ROOT.parent / "data" / "k6_seed_census" / "refl_census_hits.jsonl"):
        h = json.loads(line)
        Q = np.eye(12, dtype=np.int64)
        for v in h["left"]:
            Q = (Q @ reflection(v)) % 3
        recs.append(("refl", (json.dumps(h["left"]), tuple(h["R"])), Q, h["R"], h["I3"]))
    out = {"records": dict(Counter(c for c, *_ in recs))}
    distinct, hashes, self_dual_fail = {}, {}, 0
    for census, key, Q, R, I3 in recs:
        if (census, key) in distinct:
            continue
        G = code(Q, R)
        if ((G @ G.T) % 3).any():
            self_dual_fail += 1
        h = rowspace_sha256(G)
        distinct[(census, key)] = h
        hashes.setdefault(h, []).append((census, I3))
    out["distinct_records"] = dict(Counter(c for c, _ in distinct))
    out["distinct_records_total"] = len(distinct)
    out["self_dual_failures"] = self_dual_fail
    out["distinct_codes_rowspaces"] = len(hashes)
    out["distinct_codes_by_census"] = {c: len({h for (cc, _), h in distinct.items() if cc == c}) for c in ("hp", "refl")}
    out["codes_in_both_censuses"] = sum(1 for v in hashes.values() if {c for c, _ in v} == {"hp", "refl"})
    mult = Counter(len(v) for v in hashes.values())
    out["records_per_code_histogram"] = {str(k): v for k, v in sorted(mult.items())}
    out["codes_with_both_index_signs"] = sum(1 for v in hashes.values() if len({i for _, i in v}) > 1)
    # construction check 1: the paper's canonical seed (certificate k06_claude_canonical.json)
    cert = json.load(open(ROOT.parent / "certificates" / "k06_claude_canonical.json"))
    canon = rowspace_sha256(code(matrix_Q([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]], [0, 0, 0]) % 3, [3, 5, 8]))
    out["canonical_seed_hash_matches_certificate"] = canon == cert["rowspace_sha256"]
    out["canonical_seed_among_census_codes"] = canon in hashes
    # construction check 2: index from the generator alone (banded transfer DP) = recorded I3, on a sample
    rng = random.Random(20261020)
    sample = rng.sample(recs, min(a.sample, len(recs)))
    ok = sum(1 for _, _, Q, R, I3 in sample if banded_index(code(Q, R), [12, 12, 12]) == I3)
    out["index_check_sample"] = {"records": len(sample), "agree": ok}
    print(json.dumps(out, indent=1))
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
