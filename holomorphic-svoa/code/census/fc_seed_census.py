# Copied from claude/paper/fc_seed_census.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Ramond ground states of T = V_{Lambda(C)}//Z_6 for every k = 6 seed code of the census (7,112 codes), 20 Oct 2026.

Seeds: the distinct records of persistent/hp_census/hits.jsonl (Hadamard-Paley cells, 1,078) and of
persistent/refl_census/hits.jsonl (one or two reflections on the Paley base, 6,034), built exactly as in
census/count_k6_seed_codes.py: generator (I_36 | A), A = 3x3 block negacyclic matrix of D = Q_K S_K, E = Q_R S_R.

Invariant (field_content.md 3.2, which shows the derivation for the certificate seed; it uses only
data that every seed shares or that is computed here):
  untwisted R ground states, Z_6-invariant:   (N + 2F)/6,  supertrace I/6
  g^2- and g^4-twisted R ground states:        F/2 each,    supertrace 0
  g, g^5 (L_0 >= 19/6) and g^3 (L_0 >= 9/2):   none
so n_+ + n_- = (N + 8F)/6 and n_+ - n_- = I/6 = +-4, where
  N = number of full-support codewords of C, I = Index(C) = signed count,
  F = number of g^2-fixed full-support codewords (= full-support words of the contracted fixed code Cbar).
This needs, per seed: g^2 fixed-point free of type 3^24 (same g for all seeds), the balance condition (signed fixed
count 0, which GJF's theorem requires anyway), and Cbar self-dual of dimension 12 (twisted-sector defect
multiplicity 3^12/sqrt(det L) = 1).  All three are checked for every seed.  Input A1 (in its corrected form): the
lift of g^2 acts trivially on the g^2-fixed Ramond ground states.
Parity convention: (-1)^F = (-1)^{#2's} on untwisted ground states (field_content.md section 2).  Negating a
coordinate (an Arf twist) swaps n_+ and n_-, so the SCFT invariant is the unordered pair, which is fixed by n_+ + n_-.

Every quantity is computed twice: N, I, F from the one-period cell tables (hpcell.tables: F = tr_J(T_abs),
N = tr_J(T_abs^3), I = tr_J(T^3), I_1 = tr_J(T)), and F, the balance and self-duality of Cbar directly from the
length-72 generator (f3lib.fixed_subcode, enumeration of all 3^12 words of Cbar).  N and I are also compared with the
census record.

Equivalence classes: (Q, S, R) ~ (U Q M_sigma, S, sigma R) for sigma in PGL(2,11) (M_sigma, N_sigma the monomials
with M S = S N, census/symmetry.py) and U a row monomial.  The code changes by the monomial map
diag(U^{-1} per period | N_sigma per period), which commutes with g; so equivalent cells give isomorphic T.
Canonical key: min over sigma of (rows of Q M_sigma normalised to leading entry 1 and sorted, sorted sigma R).

Exact integer / F_3 arithmetic throughout; one process; resumable (appends to the output jsonl).
usage: python fc_seed_census.py [--limit n] [--out ../out/census/fc_seed_census.jsonl]
       python fc_seed_census.py --summarize
"""
import argparse, json, sys, time
from collections import Counter
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for p in (HERE, ROOT):
    sys.path.insert(0, str(p))
import hpcell                                                            # noqa: E402
from count_k6_seed_codes import code, reflection, rowspace_sha256, S    # noqa: E402
from f3lib import monomial_matrix, fixed_subcode, rref, is_self_dual, all_codewords   # noqa: E402
from mixed_rings import shift_monomial                                           # noqa: E402
from nice_reps import pgl_elements                                       # noqa: E402

OUT = ROOT / "out/census/fc_seed_census.jsonl"
SUMMARY = ROOT / "out/census/fc_seed_census_summary.json"
perm_g, signs_g = shift_monomial([12], 3, False, 1)
g = monomial_matrix(perm_g.tolist(), signs_g.tolist()); h = (g @ g) % 3
CONTRACT = list(range(12)) + list(range(36, 48))                        # one coordinate per 3-cycle of g^2
PGL = []
for perm, (de, _) in pgl_elements().items():
    d, e = de
    M = np.zeros((12, 12), dtype=np.int64)
    for i in range(12):
        M[i, perm[i]] = d[i] % 3
    PGL.append((perm, M))


def records():
    recs, seen = [], set()
    for line in open(ROOT.parent / "data" / "k6_seed_census" / "hp_census_hits.jsonl"):
        r = json.loads(line)
        key = ("hp", json.dumps(r["blocks"]), json.dumps(r["parities"]), json.dumps(sorted(r["R"])))
        if key not in seen:
            seen.add(key); recs.append((key, hpcell.matrix_Q(r["blocks"], r["parities"]) % 3, sorted(r["R"]), r))
    for line in open(ROOT.parent / "data" / "k6_seed_census" / "refl_census_hits.jsonl"):
        r = json.loads(line)
        key = ("refl", json.dumps(r["left"]), json.dumps(sorted(r["R"])))
        if key not in seen:
            Q = np.eye(12, dtype=np.int64)
            for v in r["left"]:
                Q = (Q @ reflection(v)) % 3
            seen.add(key); recs.append((key, Q, sorted(r["R"]), r))
    return recs


def rowcanon(Q):
    rows = []
    for row in Q:
        nz = np.nonzero(row)[0][0]
        rows.append(tuple(int(x) for x in ((row * (1 if row[nz] == 1 else 2)) % 3)))
    return tuple(sorted(rows))


def class_key(Q, R):
    return min((rowcanon((Q @ M) % 3), tuple(sorted(perm[k] for k in R))) for perm, M in PGL)


def analyse(Q, R):
    T, Ta, neg = hpcell.tables(Q, R, S, want_abs=True)
    To, Tao = T.astype(object), Ta.astype(object)
    out = {"I1": hpcell.tr_J(T, neg), "F_table": hpcell.tr_J(Ta, neg),
           "N": hpcell.tr_J(Tao.dot(Tao).dot(Tao), neg), "I": hpcell.tr_J(To.dot(To).dot(To), neg)}
    G = code(Q, R)
    assert not ((G @ G.T) % 3).any()
    Fx = fixed_subcode(G, h); Fc = rref(Fx[:, CONTRACT] % 3)[0]
    out["fixed_dim"] = int(Fx.shape[0]); out["Cbar_self_dual"] = bool(Fc.shape[0] == 12 and is_self_dual(Fc))
    e = o = 0
    for w in all_codewords(Fc):
        full = np.all(w != 0, axis=1); par = np.sum(w == 2, axis=1) % 2
        e += int((full & (par == 0)).sum()); o += int((full & (par == 1)).sum())
    out["F_direct_even_odd"] = [e, o]
    F = e + o
    assert out["F_table"] == F and out["I1"] == 0 and e == o and out["Cbar_self_dual"] and out["fixed_dim"] == 12, out
    assert (out["N"] + 8 * F) % 6 == 0 and out["I"] % 6 == 0 and abs(out["I"]) == 24
    n = (out["N"] + 8 * F) // 6; d = out["I"] // 6
    out["n_total"] = n; out["n_plus_minus"] = [(n + d) // 2, (n - d) // 2]
    out["rowspace_sha256"] = rowspace_sha256(G)
    return out


def run(limit, out_path):
    recs = records()
    done = set()
    if out_path.exists():
        done = {tuple(json.loads(l)["key"]) for l in open(out_path)}
    t0 = time.process_time(); k = 0
    with open(out_path, "a") as f:
        for key, Q, R, r in recs:
            if key in done:
                continue
            if limit is not None and k >= limit:
                break
            res = analyse(Q, R)
            assert res["N"] == r["nmax"] and res["I"] == r["I3"], (key, res["N"], r["nmax"], res["I"], r["I3"])
            ck = class_key(Q, R)
            res.update({"key": list(key), "census": key[0], "class_key": json.dumps(ck)})
            f.write(json.dumps(res) + "\n"); f.flush(); k += 1
    print(f"processed {k} records in {time.process_time() - t0:.1f} CPU s; total records {len(recs)}")


def sym2_counts(B, F):
    ue = B * (B + 1) // 2 + F * (F - 1) // 2; uo = B * F
    return ue + B, uo + F


def summarize(out_path):
    rows = [json.loads(l) for l in open(out_path)]
    n = Counter(r["n_total"] for r in rows)
    classes = {}
    for r in rows:
        classes.setdefault(r["class_key"], set()).add(r["n_total"])
    # T' = Sym^2 T and T'' = Sym^2(T' x Arf), normalised so that I(T) = -4 (Arf on T if needed, as in the paper)
    chain = {}
    for v in n:
        B, Fm = (v - 4) // 2, (v + 4) // 2
        Bp, Fp = sym2_counts(B, Fm)
        Bpp, Fpp = sym2_counts(Fp, Bp)
        assert Bp - Fp == 2 and Bpp - Fpp == -1
        chain[v] = {"T": [B, Fm], "Tprime": [Bp, Fp], "Tpp": [Bpp, Fpp]}
    s = {"records": len(rows), "by_census": dict(Counter(r["census"] for r in rows)),
         "distinct_rowspaces": len({r["rowspace_sha256"] for r in rows}),
         "distinct_n_total": len(n),
         "histogram_n_total": {str(k): v for k, v in sorted(n.items())},
         "histogram_F": {str(k): v for k, v in sorted(Counter(r["F_table"] for r in rows).items())},
         "histogram_N": {str(k): v for k, v in sorted(Counter(r["N"] for r in rows).items())},
         "index_signs": dict(Counter(r["I"] for r in rows)),
         "classes": len(classes),
         "classes_with_constant_n": sum(1 for v in classes.values() if len(v) == 1),
         "distinct_Tpp_counts": len({tuple(c["Tpp"]) for c in chain.values()}),
         "Tpp_injective_on_values": len({tuple(c["Tpp"]) for c in chain.values()}) == len(n),
         "Tprime_injective_on_values": len({tuple(c["Tprime"]) for c in chain.values()}) == len(n),
         "chain_counts": {str(k): v for k, v in sorted(chain.items())}}
    cl_census = Counter()
    for k in classes:
        cl_census[frozenset(r["census"] for r in rows if r["class_key"] == k)] += 1
    s["classes_by_census"] = {"+".join(sorted(k)): v for k, v in cl_census.items()}
    s["classes_per_n_total"] = {str(v): sum(1 for c in classes.values() if v in c) for v in sorted(n)}
    SUMMARY.write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps({k: v for k, v in s.items() if k != "chain_counts"}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--summarize", action="store_true")
    a = ap.parse_args()
    if a.summarize:
        summarize(Path(a.out))
    else:
        run(a.limit, Path(a.out))
