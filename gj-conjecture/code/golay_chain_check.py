"""Stand-alone check of the note on GJF Conjecture 1.4 (index.html), Theorem B: the Golay chain codes C_k have Index(C_k) = (-1)^k 24 at length 12k.

Needs only Python 3 and numpy; nothing else from this folder is imported.  Written from the statement of Section 3 of the note alone.
G is the ternary Golay code with generator (I_6 | A), A the rows 011111, 101221, 110122, 121012, 122101, 112210, and
(positions 0..11, entries mod 3)
    p = e0 + e4,   r = e6 - e7 - e9,   v = e3 - e8 - e10 + e11,   w = -e2 + e3 - e11.
W_k in F_3^(12k) is spanned by P_j = (p on cell j, v on cell j+1) and R_j = (r on cell j, w on cell j+1), 1 <= j < k, and
C_k = (G^k cap W_k^perp) + W_k.  The script checks

  * G: 729 words, self-dual, minimum weight 6, Index(G) = -24;
  * the finite hypotheses of Lemmas 6 and 7: the Gram data of p, r, v, w, the disjoint supports, <p, r> cap G = 0;
  * Lemma 8 and Proposition 9: the 81 x 81 transfer matrix T built from its definition, u = T s, T u = -u exactly, the shapes of T and u,
    s^T T^k s = (-1)^k 24 for k = 1..kmax by direct powers, and the unsigned count s^T U^k s = 3^(6k) = |C_k| (every word of
    C_k is counted once), all in exact Python integers;
  * brute force for k = 1..4: C_k built from the closed form, W_k totally isotropic of dimension 2(k-1), W_k cap G^k = 0,
    C_k self-dual of length 12k, Index(C_k) by enumeration over an information set (2^(6k) words), equal to s^T T^k s;
    and for k <= 3 that negating one coordinate gives a self-dual code of index -Index(C_k).

With --certificate PATH it also writes the sparse certificate of Proposition 9 (state ordering, the non-zero entries of T and of
u = T s, the vector s) as JSON; check_golay_chain_certificate.py verifies T u = -u and s . u = -24 from that file alone.

Usage:  python3 golay_chain_check.py [--kmax 30] [--brute 4] [--json out.json] [--certificate cert.json]   (about 7 s)
Exit status 1 if any check fails.  Licence: MIT (see LICENSE).
"""
import argparse, itertools, json, platform, sys, time
import numpy as np

A = [[0, 1, 1, 1, 1, 1], [1, 0, 1, 2, 2, 1], [1, 1, 0, 1, 2, 2], [1, 2, 1, 0, 1, 2], [1, 2, 2, 1, 0, 1], [1, 1, 2, 2, 1, 0]]
GOLAY = np.hstack([np.eye(6, dtype=np.int64), np.array(A, dtype=np.int64)])


def vec(entries):
    x = np.zeros(12, dtype=np.int64)
    for i, a in entries.items():
        x[i] = a % 3
    return x


P = vec({0: 1, 4: 1})
R = vec({6: 1, 7: -1, 9: -1})
V = vec({3: 1, 8: -1, 10: -1, 11: 1})
W = vec({2: -1, 3: 1, 11: -1})


def rref(M):
    """reduced row echelon form over F_3 (zero rows dropped) and the pivot columns"""
    M = np.array(M, dtype=np.int64) % 3
    if M.ndim == 1:
        M = M.reshape(1, -1)
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            M[[r, p]] = M[[p, r]]
        if M[r, c] == 2:
            M[r] = (2 * M[r]) % 3
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        piv.append(c)
        r += 1
    return M[:r], piv


def rank(M):
    M = np.asarray(M)
    return 0 if M.size == 0 else rref(M)[0].shape[0]


def kernel(M, n):
    """basis of {a in F_3^n : M a = 0}"""
    M = np.asarray(M, dtype=np.int64).reshape(-1, n) % 3
    if not M.any():
        return np.eye(n, dtype=np.int64)
    Rm, piv = rref(M)
    free = [c for c in range(n) if c not in piv]
    B = np.zeros((len(free), n), dtype=np.int64)
    for j, f in enumerate(free):
        B[j, f] = 1
        for i, p in enumerate(piv):
            B[j, p] = (-Rm[i, f]) % 3
    return B


def chi(x):
    x = np.asarray(x) % 3
    return 0 if (x == 0).any() else (1 if (x == 2).sum() % 2 == 0 else -1)


def index_info(G):
    """Index of the row space of G, enumerating the 2^dim words that are nonzero on the pivot columns"""
    Rm, piv = rref(G)
    k, n = Rm.shape
    if k == 0:
        return 0
    B = Rm[:, [c for c in range(n) if c not in piv]]
    shifts = np.arange(k, dtype=np.int64)
    total, N, step = 0, 1 << k, 1 << min(k, 16)
    for start in range(0, N, step):
        idx = np.arange(start, min(N, start + step), dtype=np.int64)
        bits = (idx[:, None] >> shifts) & 1
        sign = 1 - 2 * (bits.sum(axis=1) % 2)
        if B.shape[1]:
            X = ((1 + bits) @ B) % 3
            sign = np.where((X != 0).all(axis=1), sign * (1 - 2 * ((X == 2).sum(axis=1) % 2)), 0)
        total += int(sign.sum())
    return total


def is_self_dual(G, n):
    G = np.asarray(G) % 3
    return rank(G) == n // 2 and not ((G @ G.T) % 3).any()


def golay_words():
    coeff = np.array(list(itertools.product(range(3), repeat=6)), dtype=np.int64)
    return (coeff @ GOLAY) % 3


def junction_vectors(k):
    """the 2(k-1) generators P_j, R_j of W_k"""
    rows = []
    for j in range(k - 1):
        for a, b in ((P, V), (R, W)):
            x = np.zeros(12 * k, dtype=np.int64)
            x[12 * j:12 * j + 12], x[12 * j + 12:12 * j + 24] = a, b
            rows.append(x)
    return np.array(rows, dtype=np.int64).reshape(-1, 12 * k)


def chain_code(k):
    """generator of C_k = (G^k cap W_k^perp) + W_k, and the matrices of G^k and W_k"""
    n = 12 * k
    Gk = np.zeros((6 * k, n), dtype=np.int64)
    for j in range(k):
        Gk[6 * j:6 * j + 6, 12 * j:12 * j + 12] = GOLAY
    Wk = junction_vectors(k)
    if Wk.shape[0] == 0:
        return rref(Gk)[0], Gk, Wk
    D = (kernel(((Gk @ Wk.T) % 3).T, 6 * k) @ Gk) % 3        # coefficient vectors a with (a Gk) . W_k = 0
    return rref(np.vstack([D, Wk]))[0], Gk, Wk


def state(c1, c2, f1, f2):
    return 27 * c1 + 9 * c2 + 3 * f1 + f2


def transfer_matrices(words):
    """T[(c,f),(c',f')] = sum over g in G with psi(g) = -f, phi(g) = f' of chi(g + c'_1 p + c'_2 r + c_1 v + c_2 w);
    U the same with weight 1 (all words)"""
    T = [[0] * 81 for _ in range(81)]
    U = [[0] * 81 for _ in range(81)]
    for g in words:
        phi = (int(g @ P % 3), int(g @ R % 3))
        f = ((-int(g @ V)) % 3, (-int(g @ W)) % 3)            # psi(g) = -f
        for c1, c2, d1, d2 in itertools.product(range(3), repeat=4):
            x = g + d1 * P + d2 * R + c1 * V + c2 * W
            i, j = state(c1, c2, *f), state(d1, d2, *phi)
            T[i][j] += chi(x)
            U[i][j] += 1
    return T, U


def matvec(M, x):
    return [sum(a * b for a, b in zip(row, x)) for row in M]


def write_certificate(path, T, s, u):
    """sparse JSON certificate of Proposition 9, one matrix entry per line"""
    labels = [f"{c1}{c2}{f1}{f2}" for c1, c2, f1, f2 in itertools.product(range(3), repeat=4)]
    assert all(state(*map(int, lab)) == i for i, lab in enumerate(labels))
    header = {
        "description": "Sparse certificate for Proposition 9 of the note on GJF Conjecture 1.4 (Golay chain codes C_k, Index(C_k) = (-1)^k 24): the 81 x 81 "
                       "transfer matrix T of Lemma 8, the vector s and u = T s. Claims: u = T s, T u = -u, s . u = -24. "
                       "Written by golay_chain_check.py; check_golay_chain_certificate.py verifies the claims "
                       "from this file alone.",
        "construction": {"golay_generator_rows": [[int(t) for t in row] for row in GOLAY],
                         "p": [int(t) for t in P], "r": [int(t) for t in R], "v": [int(t) for t in V], "w": [int(t) for t in W],
                         "coordinates": "0..11 in each cell of 12; entries in F_3 = {0, 1, 2}, 2 = -1"},
        "state_ordering": {"index": "27*c1 + 9*c2 + 3*f1 + f2", "labels_c1c2f1f2": labels,
                           "meaning": "the state at the junction after cell j: c = (c1, c2) are the coefficients of P_j and R_j "
                                      "in the word, f = (f1, f2) = (g_j . p, g_j . r) for the Golay word g_j of cell j"},
        "definition": "T[(c,f),(c',f')] = sum over g in G with (g.v, g.w) = -f and (g.p, g.r) = f' of "
                      "chi(g + c'_1 p + c'_2 r + c_1 v + c_2 w), chi(x) = prod_i chi(x_i); row = state before the cell, "
                      "column = state after it; Index(C_k) = s^T T^k s",
        "s": s,
        "T_nonzero_count": sum(1 for row in T for t in row if t),
    }
    entries = [[i, j, T[i][j]] for i in range(81) for j in range(81) if T[i][j]]
    unz = [[i, x] for i, x in enumerate(u) if x]
    lines = ["{"] + [f" {json.dumps(k)}: {json.dumps(v)}," for k, v in header.items()]
    lines.append(' "T_nonzero_entries_row_col_value": [')
    lines += [f"  {json.dumps(e)}" + ("," if n < len(entries) - 1 else "") for n, e in enumerate(entries)]
    lines.append(" ],")
    lines.append(' "u_nonzero_entries_index_value": [')
    lines += [f"  {json.dumps(e)}" + ("," if n < len(unz) - 1 else "") for n, e in enumerate(unz)]
    lines.append(" ],")
    lines.append(f' "claims": {json.dumps({"u_equals_T_s": True, "T_u_equals_minus_u": True, "s_dot_u": sum(a * b for a, b in zip(s, u))})}')
    lines.append("}")
    text = "\n".join(lines) + "\n"
    json.loads(text)                                           # valid JSON
    with open(path, "w") as fh:
        fh.write(text)


class Report:
    def __init__(self):
        self.fail, self.checks = [], {}

    def check(self, name, ok, detail=""):
        self.checks[name] = bool(ok) and self.checks.get(name, True)
        if not ok:
            self.fail.append(f"{name}: {detail}")
            print("FAIL", name, detail, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--kmax", type=int, default=30, help="check s^T T^k s and s^T U^k s for k = 1..kmax")
    ap.add_argument("--brute", type=int, default=4, help="brute-force C_k for k = 1..brute (k = 4: 2^24 words, ~10 s)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--certificate", default=None, help="write the sparse certificate of Proposition 9 to this JSON file")
    a = ap.parse_args()
    rep, t0, out = Report(), time.process_time(), {}

    # the Golay code
    words = golay_words()
    weights = sorted({int((x != 0).sum()) for x in words})
    iG = sum(chi(x) for x in words)
    rep.check("G: 729 distinct words", len({tuple(x) for x in words}) == 729)
    rep.check("G: self-dual", is_self_dual(GOLAY, 12))
    rep.check("G: minimum weight 6", weights[1] == 6, str(weights))
    rep.check("G: Index(G) = -24", iG == -24, str(iG))
    out["golay_weights"], out["golay_index"] = weights, iG

    # finite hypotheses of Lemmas 6 and 7
    gram = {name: int(x @ y % 3) for name, (x, y) in
            {"p.p": (P, P), "v.v": (V, V), "p.r": (P, R), "v.w": (V, W), "r.r": (R, R), "w.w": (W, W),
             "p.v": (P, V), "p.w": (P, W), "r.v": (R, V), "r.w": (R, W)}.items()}
    out["gram"] = gram
    rep.check("Gram: p.p + v.v = p.r + v.w = r.r + w.w = 0",
              (gram["p.p"] + gram["v.v"]) % 3 == 0 and (gram["p.r"] + gram["v.w"]) % 3 == 0 and (gram["r.r"] + gram["w.w"]) % 3 == 0)
    supp_pr = set(np.nonzero(P)[0]) | set(np.nonzero(R)[0])
    supp_vw = set(np.nonzero(V)[0]) | set(np.nonzero(W)[0])
    rep.check("supports of p, r and of v, w are disjoint", not (supp_pr & supp_vw), f"{sorted(supp_pr)} {sorted(supp_vw)}")
    out["supports"] = {"p,r": sorted(int(i) for i in supp_pr), "v,w": sorted(int(i) for i in supp_vw)}
    gset = {tuple(x) for x in words}
    rep.check("<p, r> cap G = 0", all(tuple((s * P + t * R) % 3) not in gset
                                      for s, t in itertools.product(range(3), repeat=2) if (s, t) != (0, 0)))

    # transfer matrix
    T, U = transfer_matrices(words)
    s = [1 if i // 9 == 0 else 0 for i in range(81)]           # indicator of c = 0
    u = matvec(T, s)
    Tu = matvec(T, u)
    rep.check("T u = -u", all(x == -y for x, y in zip(Tu, u)))
    nz_T = [t for row in T for t in row if t]
    out["T_nonzero_entries"], out["T_entry_values"] = len(nz_T), sorted(set(nz_T))
    out["u_nonzero_entries"] = {f"{i // 27}{i // 9 % 3}{i // 3 % 3}{i % 3}": x for i, x in enumerate(u) if x}
    out["u_entry_values"] = sorted({x for x in u if x})
    rep.check("s^T u = Index(G)", sum(x * y for x, y in zip(s, u)) == iG)
    idx, cnt, x, y = [], [], s[:], s[:]
    for k in range(1, a.kmax + 1):
        x, y = matvec(T, x), matvec(U, y)
        idx.append(sum(si * xi for si, xi in zip(s, x)))
        cnt.append(sum(si * yi for si, yi in zip(s, y)))
    rep.check("s^T T^k s = (-1)^k 24 for k = 1..kmax", all(v == (-1) ** k * 24 for k, v in enumerate(idx, 1)), str(idx))
    rep.check("s^T U^k s = 3^(6k) = |C_k| for k = 1..kmax", all(v == 3 ** (6 * k) for k, v in enumerate(cnt, 1)))
    out["transfer_index_k1_to_kmax"] = idx
    if a.certificate:
        write_certificate(a.certificate, T, s, u)

    # brute force
    out["brute_force"] = {}
    for k in range(1, a.brute + 1):
        n = 12 * k
        C, Gk, Wk = chain_code(k)
        rep.check("W_k totally isotropic", not ((Wk @ Wk.T) % 3).any())
        rep.check("dim W_k = 2(k-1) and W_k cap G^k = 0", rank(Wk) == 2 * (k - 1) and rank(np.vstack([Gk, Wk])) == 6 * k + 2 * (k - 1))
        rep.check("C_k self-dual of length 12k", is_self_dual(C, n))
        I = index_info(C)
        rep.check("brute force Index(C_k) = s^T T^k s = (-1)^k 24", I == idx[k - 1] == (-1) ** k * 24, f"k={k} {I}")
        rec = {"length": n, "index": I}
        if k <= 3:
            Cf = C.copy()
            Cf[:, 0] = (-Cf[:, 0]) % 3
            If = index_info(Cf)
            rep.check("one coordinate negated: self-dual, index -Index(C_k)", is_self_dual(Cf, n) and If == -I, f"k={k} {If}")
            rec["index_after_negating_coordinate_0"] = If
        out["brute_force"][str(k)] = rec
        print(f"k = {k}: length {n}, Index(C_k) = {I}", flush=True)

    out.update({"script": "golay_chain_check.py", "python": platform.python_version(), "numpy": np.__version__,
                "kmax": a.kmax, "brute": a.brute, "cpu_seconds": round(time.process_time() - t0, 1),
                "checks": rep.checks, "failures": rep.fail})
    print(json.dumps({k: v for k, v in out.items() if k not in ("u_nonzero_entries", "transfer_index_k1_to_kmax")}, indent=1))
    print(f"{len(rep.checks)} checks, {len(rep.fail)} failures")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1)
            fh.write("\n")
    sys.exit(1 if rep.fail else 0)


if __name__ == "__main__":
    main()
