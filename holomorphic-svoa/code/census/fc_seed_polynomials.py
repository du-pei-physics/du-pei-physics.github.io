# Copied from claude/paper/fc_seed_polynomials.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""NS character of T = V_{Lambda(C)}//Z_6 as a polynomial in t for all 7,112 k = 6 seeds; data file for the online note.
21 Oct 2026.  Reads fc_seed_census.jsonl and fc_seed_level2_part*of2.jsonl (coefficients a_j from fc_seed_level2.py).

Per-seed checks (all asserted): the q-expansion of sum a_j t^j reproduces dim T_h for h = 0, 1/2, 1, 3/2, 2, 5/2 as
computed from the code, its constant term a_0 equals n_+ + n_- (Ramond ground states, fc_seed_census.py), and
-4096 a_1 equals the first excited Ramond level D1 computed independently from the twisted sectors (fc_seed_level1.py).
Then: number of distinct polynomials (= distinct NS characters), their distribution over seeds and symmetry classes,
and whether distinct Z_T give distinct Z_{T''}: below weight 6 (resp. 12) T' (resp. T'') is the graded-symmetric part
of T x T (resp. T' x T'), so dim T''_h for h <= 3 is dim T_h plus a polynomial in lower dimensions (fc_sym2.sym2_ns);
we compare these tuples.  Output: out/census/k6_seed_polynomials.csv and k6_seed_polynomials_summary.json.
"""
import hashlib, json, sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTDIR = ROOT / "out/census"
sys.path.insert(0, str(ROOT))
import fc_closed_form as CF                     # noqa: E402

HEADER = """# NS characters of the k = 6 seed theories T = V_{Lambda(C)}//Z_6 (c = 72, Witten index +-4), one line per seed code.
# Z_T(tau) = Tr_NS q^(L_0 - 3) (NS sector, no (-1)^F inserted) = a6 t^6 + a5 t^5 + ... + a1 t + a0,
#   t = (theta_3/eta)^12 = eta(tau)^48/(eta(tau/2) eta(2tau))^24 = q^(-1/2) + 24 + 276 q^(1/2) + ...  (NS character of 24 free fermions).
# The same coefficients give Tr_NS (-1)^F q^(L_0-3) = sum (-1)^j a_j (theta_4/eta)^(12j) and the Ramond character
#   Tr_R q^(L_0-3) = sum (-1)^j a_j (theta_2/eta)^(12j); a0 = number of Ramond ground states, -4096*a1 = first excited Ramond level.
# sign_I: sign of the Witten index of T = Index(C)/6, with (-1)^F = (-1)^(number of 2's) on the untwisted Ramond ground states
#   (field_content.md section 2); an Arf twist (negating one coordinate of C) flips it and leaves Z_T unchanged.
# census: hp = Hadamard-Paley cell census, refl = reflection census; key = census record (blocks|parities|R or left|R);
# rowspace_sha256 = SHA-256 of the rows of (I_36 | A) over {0,1,2} joined by newlines; class_id = symmetry class under
#   PGL(2,11) x row monomials (period-uniform monomial code equivalences; equal class => isomorphic T).
# Source: claude/paper/fc_seed_census.py, fc_seed_level1.py, fc_seed_level2.py, fc_seed_polynomials.py (scft-code-search)."""


def sym2_ns(chi, chiAP, hmax):
    """untwisted graded-symmetric square (same formula as fc_sym2.sym2_ns, which runs a script on import)."""
    def mul(a, b):
        out = {}
        for h1, n1 in a.items():
            for h2, n2 in b.items():
                if h1 + h2 <= hmax:
                    out[h1 + h2] = out.get(h1 + h2, 0) + n1 * n2
        return out
    sq, sqAP = mul(chi, chi), mul(chiAP, chiAP)
    dbl = {2 * h: n for h, n in chiAP.items() if 2 * h <= hmax}
    dim, sup = {}, {}
    for h in set(sq) | set(dbl):
        x = sq.get(h, 0) + dbl.get(h, 0); y = sqAP.get(h, 0) + dbl.get(h, 0)
        assert x % 2 == 0 and y % 2 == 0
        dim[h], sup[h] = x // 2, y // 2
    return dim, sup


def main():
    first = {tuple(json.loads(l)["key"]): json.loads(l) for l in open(OUTDIR / "fc_seed_census.jsonl")}
    lev1, lev2 = {}, {}
    for p in OUTDIR.glob("fc_seed_level1_part*of2.jsonl"):
        for l in open(p):
            r = json.loads(l); lev1[tuple(r["key"])] = r
    for p in OUTDIR.glob("fc_seed_level2_part*of2.jsonl"):
        for l in open(p):
            r = json.loads(l); lev2[tuple(r["key"])] = r
    assert set(first) == set(lev1) == set(lev2) and len(first) == 7112
    class_ids = {k: i for i, k in enumerate(sorted({r["class_key"] for r in first.values()}))}
    lines, poly_of = [], {}
    for k in sorted(first, key=lambda k: (first[k]["census"] != "hp", k)):
        f, s1, s2 = first[k], lev1[k], lev2[k]
        a = s2["a"]
        Z = CF.eval_poly(a)
        dims = {-6: 1, -5: 0, -4: s2["T1"], -3: s2["T32"], -2: s2["T2"], -1: s2["T52"]}
        assert all(Z.get(e, 0) == v for e, v in dims.items()), k
        assert a[0] == f["n_total"] and -4096 * a[1] == s1["D1"] and a[6] == 1 and a[5] == -144, k
        assert CF.eval_poly(a, "R").get(0) == f["n_total"]
        poly_of[k] = tuple(a)
        key = "|".join(str(x) for x in k[1:])
        lines.append(",".join([f["census"], '"' + key.replace('"', "'") + '"', f["rowspace_sha256"], str(class_ids[f["class_key"]]),
                               "+" if f["I"] > 0 else "-"] + [str(a[j]) for j in range(6, -1, -1)]))
    out = OUTDIR / "k6_seed_polynomials.csv"
    out.write_text(HEADER + "\ncensus,key,rowspace_sha256,class_id,sign_I,a6,a5,a4,a3,a2,a1,a0\n" + "\n".join(lines) + "\n")
    # certificate seed
    cert = [k for k in first if first[k]["rowspace_sha256"] == json.load(open(ROOT.parent / "certificates" / "k06_claude_canonical.json"))["rowspace_sha256"]]
    assert len(cert) == 1 and list(poly_of[cert[0]]) == [13940, -629867, 751404, -129679, 6984, -144, 1]
    # distribution
    by_poly_seeds = Counter(poly_of.values())
    by_poly_classes = {}
    for k, p in poly_of.items():
        by_poly_classes.setdefault(p, set()).add(first[k]["class_key"])
    per_class = {}
    for k, p in poly_of.items():
        per_class.setdefault(first[k]["class_key"], set()).add(p)
    cls_const = all(len(v) == 1 for v in per_class.values())
    # T'' separation via NS dimensions up to weight 3
    tpp = {}
    for p in by_poly_seeds:
        Z = CF.eval_poly(list(p)); ZAP = CF.eval_poly(list(p), "AP")
        chi = {Fraction(e + 6, 2): v for e, v in Z.items() if e <= 0}
        chiAP = {Fraction(e + 6, 2): v for e, v in ZAP.items() if e <= 0}
        Tp, TpAP = sym2_ns(chi, chiAP, Fraction(3))
        Tpp, _ = sym2_ns(Tp, TpAP, Fraction(3))
        tpp[p] = tuple(Tpp.get(Fraction(h, 2), 0) for h in range(7))
    cert_tpp = tpp[poly_of[cert[0]]]
    assert cert_tpp == (1, 0, 0, 113, 2928, 42001, 975100), cert_tpp
    hp = [k for k in first if first[k]["census"] == "hp"]
    summ = {"seeds": len(first), "distinct_polynomials": len(by_poly_seeds),
            "polynomials_with_one_seed": sum(1 for v in by_poly_seeds.values() if v == 1),
            "polynomials_with_one_class": sum(1 for v in by_poly_classes.values() if len(v) == 1),
            "largest_group_seeds": max(by_poly_seeds.values()),
            "largest_group_classes": max(len(v) for v in by_poly_classes.values()),
            "classes_per_polynomial_histogram": {str(c): n for c, n in sorted(Counter(len(v) for v in by_poly_classes.values()).items())},
            "constant_on_classes": cls_const, "classes": len(class_ids),
            "hp_distinct_polynomials": len({poly_of[k] for k in hp}), "hp_classes": len({first[k]["class_key"] for k in hp}),
            "distinct_Tpp_NS_dims_to_weight_3": len(set(tpp.values())),
            "Z_T_to_Z_Tpp_injective": len(set(tpp.values())) == len(by_poly_seeds),
            "certificate_Tpp_dims_0_to_3": cert_tpp,
            "file": str(out.relative_to(ROOT.parent)), "file_sha256": hashlib.sha256(out.read_bytes()).hexdigest()}
    (OUTDIR / "k6_seed_polynomials_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
