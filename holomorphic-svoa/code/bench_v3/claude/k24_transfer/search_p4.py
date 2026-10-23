"""Period-p rings of q=12 cells with a common U, done on the tiny (norm class x J) blocks.

All tables T_S with the same U are block diagonal in the same basis (Theorem 1), so a product
M = T_1 ... T_p is computed block by block (blocks are <= 3x3 for r = 2, <= 5x5 for r = 3).
Targets (GJF Thm 5.6, g = shift by p cells to the power L/(3p)):
  P288_p4 : 12 cells, p = 4 :  Tr(J M) = 0, |Tr(J M^3)| = 24   -> c = 288, orbifold index 1
  P288_p2 : 12 cells, p = 2 :  Tr(J M^2) = 0, |Tr(J M^6)| = 24 -> c = 288, orbifold index 1
  W144_p2 :  6 cells, p = 2 :  Tr(J M) = 0, |Tr(J M^3)| = 24   -> c = 144, orbifold index 2
Random quads are screened mod two primes (int64-exact), survivors are recomputed exactly.
usage: python search_p4.py --minutes 10 --seed 1 --r 2 --out results/p4_a.json
"""
import argparse, json, time
from collections import Counter
from pathlib import Path
import numpy as np
import core
from bases import BASES12

PR = (1048573, 1048571)


def pool_blocks(rng, B, N, depth_max=6):
    names = list(BASES12)
    specs, blks = [], []
    shapes = None
    tries = 0
    while len(specs) < N and tries < 20 * N:
        tries += 1
        S = core.random_orth(rng, int(rng.integers(1, depth_max + 1)), BASES12[names[int(rng.integers(4))]])
        T, neg, norm = core.transfer(B, S)
        if not T.any():
            continue
        bl = core.blocks(T, neg, norm)
        shapes = [(nc, j, m.shape[0]) for nc, j, m in bl]
        specs.append(S); blks.append([m for _, _, m in bl])
    signs = [j for _, j, _ in shapes]
    stacks = [np.array([b[k] for b in blks], dtype=np.int64) for k in range(len(shapes))]   # per block: (N, s, s)
    return specs, stacks, signs


def prod_mod(stack, idx, p):
    """Batch product over the index columns of idx (K, m) of matrices stack[(N, s, s)] mod p."""
    P = stack[idx[:, 0]] % p
    for c in range(1, idx.shape[1]):
        P = np.einsum('kab,kbc->kac', P, stack[idx[:, c]] % p) % p
    return P


def trpow_mod(P, L, p):
    X = P.copy()
    for _ in range(L - 1):
        X = np.einsum('kab,kbc->kac', X, P) % p
    return np.trace(X, axis1=1, axis2=2) % p


def indexJ_mod(stacks, signs, idx, L, p):
    tot = np.zeros(idx.shape[0], dtype=np.int64)
    for st, sg in zip(stacks, signs):
        tot = (tot + sg * trpow_mod(prod_mod(st, idx, p), L, p)) % p
    return tot


def indexJ_exact(stacks, signs, cols, L):
    tot = 0
    for st, sg in zip(stacks, signs):
        M = np.eye(st.shape[1], dtype=object)
        for c in cols:
            M = M.dot(st[c].astype(object))
        X = np.eye(st.shape[1], dtype=object)
        for _ in range(L):
            X = X.dot(M)
        tot += sg * int(np.trace(X))
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--minutes", type=float, default=10); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r", type=int, default=2); ap.add_argument("--disc", type=int, default=0)
    ap.add_argument("--pool", type=int, default=400); ap.add_argument("--batch", type=int, default=200000)
    ap.add_argument("--rounds", type=int, default=40); ap.add_argument("--out", required=True)
    ap.add_argument("--dmax", type=int, default=6, help="max reflections in each pool cell")
    ap.add_argument("--cyclic", action="store_true", help="cyclic ring (z^L = +1, g = -z^{L/3}): plain traces instead of J-traces")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed); t_end = time.time() + 60 * a.minutes
    st = Counter(); small = {k: Counter() for k in ("P288_p4", "P288_p2", "W144_p2")}; hits = []; seen = set()
    TARG = [("P288_p4", 4, 1, 3), ("P288_p2", 2, 2, 6), ("W144_p2", 2, 1, 3)]   # name, p, L_balance, L_index
    while time.time() < t_end:
        while True:
            B = core.random_subspace(rng, a.r)
            disc = int(round(np.linalg.det((B @ B.T) % 3))) % 3
            if a.disc in (0, disc):
                break
        specs, stacks, signs = pool_blocks(rng, B, a.pool, depth_max=a.dmax)
        if a.cyclic:
            signs = [1] * len(signs)
        N = len(specs)
        if N < 20:
            continue
        st["pools"] += 1; st["cells"] += N
        for _ in range(a.rounds):
            if time.time() > t_end:
                break
            for name, p, Lb, Li in TARG:
                idx = rng.integers(0, N, (a.batch if p == 4 else a.batch // 4, p))
                st[name + "_tested"] += len(idx)
                b0 = indexJ_mod(stacks, signs, idx, Lb, PR[0])
                keep = b0 == 0
                if not keep.any():
                    continue
                idx = idx[keep]
                i0 = indexJ_mod(stacks, signs, idx, Li, PR[0])
                i1 = indexJ_mod(stacks, signs, idx, Li, PR[1])
                nz = (i0 != 0)
                st[name + "_balanced_nondeg_modp"] += int(nz.sum())
                cand = ((i0 == 24) | (i0 == PR[0] - 24)) & ((i1 == 24) | (i1 == PR[1] - 24))
                # small-value statistics on a sample of balanced nondegenerate tuples
                for row in idx[nz][:5]:
                    if indexJ_exact(stacks, signs, row, Lb) == 0:
                        small[name][abs(indexJ_exact(stacks, signs, row, Li))] += 1
                for row in idx[cand]:
                    if indexJ_exact(stacks, signs, row, Lb) != 0:
                        continue
                    v = indexJ_exact(stacks, signs, row, Li)
                    if abs(v) != 24:
                        continue
                    key = (name, B.tobytes(), tuple(int(x) for x in row))
                    if key in seen:
                        continue
                    seen.add(key)
                    hits.append({"target": name, "r": a.r, "disc": disc, "B": B.tolist(),
                                 "S_period": [specs[int(k)].tolist() for k in row], "index": v})
                    Path(a.out).write_text(json.dumps({"partial": True, "hits": hits}, indent=1))
    summ = {"r": a.r, "twist": "cyclic" if a.cyclic else "negacyclic", "stats": dict(st), "smallest_balanced_abs": {k: sorted(v.items())[:12] for k, v in small.items()},
            "hits": dict(Counter(h["target"] for h in hits))}
    Path(a.out).write_text(json.dumps({"summary": summ, "hits": hits}, indent=1))
    print(json.dumps(summ))


if __name__ == "__main__":
    main()
