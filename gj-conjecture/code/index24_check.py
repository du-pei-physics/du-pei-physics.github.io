"""Stand-alone numerical check of the note on GJF Conjecture 1.4 (index.html), Theorem A: 24 | Index(C) for self-dual ternary codes of length 12k.

Needs only Python 3 and numpy; nothing else from this folder is imported.  Index(X) is the signed number of
full-support vectors of a set X of vectors in F_3^n, the sign of x being (-1)^#{i : x_i = 2}.  The script checks

  * the index routine (enumeration over an information set, 2^dim words) against plain enumeration of all 3^dim words
    on random codes, on the ternary Golay code (|Index| = 24) and on Golay + Golay (Index 576);
  * Lemma 1: 3^dim(W) Index(W^perp) = i^n 3^(n/2) Index(W) for random linear W, both sides by plain enumeration
    (n = 3..14: for 4 | n this is the lemma, for n = 2 mod 4 the sign flips, for odd n both sides vanish), and
    Lemma 3, Index(W) even for even n, on the same codes;
  * Lemma 2 along random Kneser-neighbour walks from the sum of n/4 tetracodes: at every step C' = D + <v> with
    v isotropic and not in C, D = C cap v^perp; it checks that C' is self-dual and D = C cap C' of dimension n/2 - 1,
    Index(C) + Index(C') = 4 Index(D) and Index(D) even; for n <= 36 also Index(D^perp) = 3 Index(D) =
    Index(C) + Index(C') - Index(D); for n <= 24 also that the four cosets of D in D^perp outside C and C' contain only
    vectors of nonzero norm and no full-support vector;
  * Lemma 4: the explicit chain C = C_0, ..., C_j = C' built from dual bases of C/E and C'/E, E = C cap C',
    between a code reached by the walk and a signed permutation of the tetracode sum: every C_m is self-dual,
    consecutive ones meet in dimension n/2 - 1, and Index(C_m) + Index(C_m+1) = 0 mod 8;
  * Proposition 5, 8 | Index, for every code met, and, as an observation only (the factor 3 is GJF's Theorem 1.3,
    not re-proved here), 24 | Index.

Usage:  python3 index24_check.py                      # default suite, n = 12, 24, 36 (about a minute)
        python3 index24_check.py --n 48 --steps 4      # one length only
        python3 index24_check.py --json result.json    # also write the summary as JSON
Exit status 1 if any check fails.  Licence: MIT (see LICENSE).
"""
import argparse, itertools, json, platform, sys, time
import numpy as np

TETRA = [[1, 0, 1, 1], [0, 1, 1, 2]]                          # the tetracode; all nonzero words have weight 3


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


def dual(G, n):
    """basis of {y : G y = 0}, the dual of the row space of G"""
    G = np.asarray(G, dtype=np.int64).reshape(-1, n) % 3
    if not G.any():
        return np.eye(n, dtype=np.int64)
    R, piv = rref(G)
    free = [c for c in range(n) if c not in piv]
    B = np.zeros((len(free), n), dtype=np.int64)
    for j, f in enumerate(free):
        B[j, f] = 1
        for i, p in enumerate(piv):
            B[j, p] = (-R[i, f]) % 3
    return B


def inv_mod3(A):
    m = A.shape[0]
    R, piv = rref(np.hstack([np.asarray(A) % 3, np.eye(m, dtype=np.int64)]))
    if piv[:m] != list(range(m)):
        raise ValueError("matrix not invertible mod 3")
    return R[:, m:]


def index_plain(G):
    """Index of the row space of G by enumerating all 3^dim codewords"""
    R, _ = rref(G)
    k = R.shape[0]
    if k == 0:
        return 0                                               # the zero word has no full support
    U = np.array(list(itertools.product(range(3), repeat=k)), dtype=np.int64)
    X = (U @ R) % 3
    full = (X != 0).all(axis=1)
    sign = 1 - 2 * ((X == 2).sum(axis=1) % 2)
    return int(sign[full].sum())


def index_info(G):
    """Index of the row space of G: a full-support word is nonzero on the pivot columns, so only the 2^dim words with
    coordinates in {1, 2} there are enumerated"""
    R, piv = rref(G)
    k, n = R.shape
    if k == 0:
        return 0
    B = R[:, [c for c in range(n) if c not in piv]]
    shifts = np.arange(k, dtype=np.int64)
    total, N, step = 0, 1 << k, 1 << min(k, 16)
    for start in range(0, N, step):
        idx = np.arange(start, min(N, start + step), dtype=np.int64)
        bits = (idx[:, None] >> shifts) & 1                   # 1 means coordinate value 2 on the pivot column
        sign = 1 - 2 * (bits.sum(axis=1) % 2)
        if B.shape[1]:
            X = ((1 + bits) @ B) % 3
            sign = np.where((X != 0).all(axis=1), sign * (1 - 2 * ((X == 2).sum(axis=1) % 2)), 0)
        total += int(sign.sum())
    return total


def tetracode_sum(n):
    G = np.zeros((n // 2, n), dtype=np.int64)
    for b in range(n // 4):
        G[2 * b:2 * b + 2, 4 * b:4 * b + 4] = TETRA
    return G


def golay():
    """ternary Golay code [I | B], B the symmetric conference matrix of order 6 (Legendre symbol mod 5, bordered)"""
    leg = {0: 0, 1: 1, 4: 1, 2: -1, 3: -1}
    B = np.zeros((6, 6), dtype=np.int64)
    B[0, 1:] = B[1:, 0] = 1
    for i in range(5):
        for j in range(5):
            B[1 + i, 1 + j] = leg[(j - i) % 5]
    return np.hstack([np.eye(6, dtype=np.int64), B % 3])


def is_self_dual(G, n):
    G = np.asarray(G) % 3
    return rank(G) == n // 2 and not ((G @ G.T) % 3).any()


def contains(G, H):
    """row space of G contains that of H"""
    return rank(np.vstack([G, H])) == rank(G)


def random_neighbour(G, rng):
    """C' = D + <v> with v isotropic and v not in C, D = {c in C : c.v = 0}; both in reduced echelon form"""
    G = np.asarray(G) % 3
    m, n = G.shape
    while True:
        v = rng.integers(0, 3, n)
        s = (G @ v) % 3
        if (v @ v) % 3 == 0 and s.any():                     # v isotropic and v not in C^perp = C
            break
    D = rref((dual(s.reshape(1, -1), m) @ G) % 3)[0]          # coefficient vectors a with a.s = 0
    return rref(np.vstack([D, v]))[0], D


def signed_permutation(G, rng):
    n = G.shape[1]
    return (G[:, rng.permutation(n)] * rng.choice([1, 2], n)) % 3


class Report:
    def __init__(self):
        self.fail, self.counts, self.indices, self.d_indices = [], {}, {}, {}

    def check(self, name, ok, detail=""):
        self.counts[name] = self.counts.get(name, 0) + 1
        if not ok:
            self.fail.append(f"{name}: {detail}")
            print("FAIL", name, detail, flush=True)

    def note(self, name):
        self.counts[name] = self.counts.get(name, 0) + 1

    def seen(self, n, value, table=None):
        d = (self.indices if table is None else table).setdefault(str(n), {})
        d[str(value)] = d.get(str(value), 0) + 1


def check_index_routine(rep, rng, trials):
    for _ in range(trials):
        n = int(rng.integers(2, 13))
        k = int(rng.integers(1, min(n, 9) + 1))
        G = rng.integers(0, 3, (k, n))
        if rng.random() < 0.5:
            G[0] = rng.choice([1, 2], n)                       # make full-support words likely
        rep.check("index routine: information set = plain enumeration", index_info(G) == index_plain(G), str(G.tolist()))
    g = golay()
    rep.check("Golay code self-dual", is_self_dual(g, 12))
    ig = index_plain(g)
    rep.check("Golay code |Index| = 24", abs(ig) == 24 and index_info(g) == ig, str(ig))
    gg = np.zeros((12, 24), dtype=np.int64)
    gg[:6, :12] = gg[6:, 12:] = g
    rep.check("Golay + Golay Index = 576", index_info(gg) == 576 == ig * ig, str(index_info(gg)))
    rep.check("tetracode self-dual, nonzero weights all 3",
              is_self_dual(TETRA, 4) and {int(w) for w in ((np.array(list(itertools.product(range(3), repeat=2)))[1:]
                                                            @ np.array(TETRA)) % 3 != 0).sum(axis=1)} == {3})


def check_lemmas_1_3(rep, rng, trials):
    for t in range(trials):
        n = int(rng.integers(3, 15))
        d = int(rng.integers(max(1, n - 11), min(n - 1, 11) + 1))
        W = rng.integers(0, 3, (d, n))
        if t % 2 == 0:
            W[0] = rng.choice([1, 2], n)
        R = rref(W)[0]
        d = R.shape[0]
        H = dual(R, n)
        iw, ip = index_plain(R), index_plain(H)
        if n % 2:
            rep.check("Lemma 1, odd n: both indices vanish", iw == 0 == ip, f"n={n} {iw} {ip}")
        else:
            sign = 1 if n % 4 == 0 else -1
            name = "Lemma 1, 4 | n" if n % 4 == 0 else "Lemma 1, n = 2 mod 4 (sign -1)"
            rep.check(name, 3 ** d * ip == sign * 3 ** (n // 2) * iw, f"n={n} dim={d} {iw} {ip}")
            if sign < 0 and iw:
                rep.note("(count) Lemma 1, n = 2 mod 4 with Index(W) != 0: the formula without the sign fails")
            rep.check("Lemma 3: Index even for even n", iw % 2 == 0 and ip % 2 == 0, f"{iw} {ip}")


def coset_check(rep, G, Gp, D, n):
    """the cosets D + a e + b f (a, b nonzero) of D in D^perp: nonzero norm, no full-support vector"""
    e = next(r for r in G if not contains(D, r))
    f = next(r for r in Gp if not contains(D, r))
    rep.check("Lemma 2: D^perp/D is a hyperbolic plane spanned by C/D and C'/D",
              (e @ e) % 3 == 0 and (f @ f) % 3 == 0 and (e @ f) % 3 != 0 and rank(dual(D, n)) == rank(np.vstack([D, e, f])))
    U = np.array(list(itertools.product(range(3), repeat=D.shape[0])), dtype=np.int64)
    words = (U @ D) % 3
    for a in (1, 2):
        for b in (1, 2):
            X = (words + a * e + b * f) % 3
            rep.check("Lemma 2: anisotropic cosets have nonzero norm", ((X * X).sum(axis=1) % 3 != 0).all())
            rep.check("Lemma 2: anisotropic cosets have no full-support vector", not (X != 0).all(axis=1).any())


def walk(rep, n, steps, rng):
    G = rref(tetracode_sum(n))[0]
    I = index_info(G)
    rep.check("base: tetracode sum self-dual with Index 0", is_self_dual(G, n) and I == 0, str(I))
    for _ in range(steps):
        Gp, D = random_neighbour(G, rng)
        rep.check("neighbour: C' self-dual", is_self_dual(Gp, n))
        rep.check("neighbour: D = C cap C' of dimension n/2 - 1",
                  D.shape[0] == n // 2 - 1 and contains(G, D) and contains(Gp, D) and rank(np.vstack([G, Gp])) == n // 2 + 1)
        Ip, ID = index_info(Gp), index_info(D)
        rep.check("Lemma 2: Index(C) + Index(C') = 4 Index(D)", I + Ip == 4 * ID, f"n={n} {I} {Ip} {ID}")
        rep.check("Lemma 3: Index(D) even", ID % 2 == 0, str(ID))
        rep.seen(n, ID, rep.d_indices)
        if n <= 36:
            IDp = index_info(dual(D, n))
            rep.check("Lemma 1 for D: Index(D^perp) = 3 Index(D)", IDp == 3 * ID, f"{IDp} {ID}")
            rep.check("Lemma 2: Index(D^perp) = Index(C) + Index(C') - Index(D)", IDp == I + Ip - ID, f"{IDp} {I} {Ip} {ID}")
        if n <= 24:
            coset_check(rep, G, Gp, D, n)
        rep.check("Theorem: 8 | Index", Ip % 8 == 0, str(Ip))
        rep.check("observed: 24 | Index", Ip % 24 == 0, str(Ip))
        rep.seen(n, Ip)
        G, I = Gp, Ip
    return G


def chain(rep, n, C, Cp):
    """Lemma 4: C_m = E + <f_1..f_m, e_m+1..e_j> with e_a . f_b = delta_ab"""
    E = dual(np.vstack([C, Cp]), n)                            # (C + C')^perp = C cap C'
    E = rref(E)[0] if E.shape[0] and E.any() else np.zeros((0, n), dtype=np.int64)

    def complete(base, rows):
        basis, ext = [r for r in base], []
        for r in rows:
            if rank(np.array(basis + [r])) > len(basis):
                basis.append(r)
                ext.append(r)
        return np.array(ext, dtype=np.int64).reshape(-1, n)
    Ev, Fv = complete(E, C), complete(E, Cp)
    j = n // 2 - E.shape[0]
    rep.check("Lemma 4: C/E and C'/E have dimension n/2 - dim E", Ev.shape[0] == j == Fv.shape[0])
    if j == 0:
        return 0
    F2 = (inv_mod3((Ev @ Fv.T) % 3).T @ Fv) % 3                # dual basis: e_a . f_b = delta_ab
    rep.check("Lemma 4: dual bases", ((Ev @ F2.T) % 3 == np.eye(j, dtype=np.int64)).all())
    codes = [np.vstack([E, F2[:m], Ev[m:]]) % 3 for m in range(j + 1)]
    rep.check("Lemma 4: chain ends at C and C'", contains(codes[0], C) and contains(C, codes[0])
              and contains(codes[-1], Cp) and contains(Cp, codes[-1]))
    idx = []
    for Cm in codes:
        rep.check("Lemma 4: every C_m self-dual", is_self_dual(Cm, n))
        idx.append(index_info(Cm))
        rep.seen(n, idx[-1])
    for a, b, ia, ib in zip(codes, codes[1:], idx, idx[1:]):
        rep.check("Lemma 4: consecutive codes are neighbours", rank(np.vstack([a, b])) == n // 2 + 1)
        rep.check("Theorem: Index(C_m) + Index(C_m+1) = 0 mod 8", (ia + ib) % 8 == 0, f"{ia} {ib}")
    rep.check("Theorem: 8 | Index along the chain", all(i % 8 == 0 for i in idx), str(idx))
    rep.check("observed: 24 | Index along the chain", all(i % 24 == 0 for i in idx), str(idx))
    return j


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n", type=int, default=None, help="one length (multiple of 12); default: 12, 24, 36")
    ap.add_argument("--steps", type=int, default=None, help="neighbour steps per length")
    ap.add_argument("--chains", type=int, default=None, help="Lemma-4 chains per length")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rep, t0 = Report(), time.process_time()
    plan = [(a.n, a.steps or 4, a.chains if a.chains is not None else 1)] if a.n else \
        [(12, a.steps or 60, a.chains if a.chains is not None else 10), (24, a.steps or 60, a.chains if a.chains is not None else 6),
         (36, a.steps or 30, a.chains if a.chains is not None else 2)]
    if any(n % 12 for n, _, _ in plan):
        raise SystemExit("the lengths must be multiples of 12")
    check_index_routine(rep, rng, 200)
    check_lemmas_1_3(rep, rng, 600)
    print(f"index routine, Lemmas 1 and 3 on random codes: {len(rep.fail)} failures", flush=True)
    lengths = {}
    for n, steps, chains in plan:
        end = walk(rep, n, steps, rng)
        js = [chain(rep, n, rref(signed_permutation(end if c % 2 == 0 else tetracode_sum(n), rng))[0],
                    rref(signed_permutation(tetracode_sum(n), rng))[0] if c % 2 == 0 else end) for c in range(chains)]
        lengths[str(n)] = {"neighbour_steps": steps, "chains": chains, "chain_lengths": js}
        print(f"n = {n}: {steps} neighbour steps, {chains} chains of lengths {js}; failures so far {len(rep.fail)}", flush=True)
    out = {"script": "index24_check.py", "seed": a.seed, "python": platform.python_version(), "numpy": np.__version__,
           "cpu_seconds": round(time.process_time() - t0, 1), "lengths": lengths, "checks": rep.counts,
           "failures": rep.fail, "indices_met": rep.indices, "indices_of_D_met": rep.d_indices}
    print(json.dumps({k: v for k, v in out.items() if k != "checks"}, indent=1))
    print(f"{sum(rep.counts.values())} checks, {len(rep.fail)} failures")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1)
            fh.write("\n")
    sys.exit(1 if rep.fail else 0)


if __name__ == "__main__":
    main()
