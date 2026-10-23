"""Numerical test of the neighbor identity for the index of self-dual ternary codes (19 Oct 2026).

Claim: if C, C' are self-dual codes of length n = 12k meeting in D = C cap C' of codimension 1, then
    Index(C) + Index(C') = 4 Index(D),
where Index(X) = sum over full-support x in X of prod_i chi(x_i) (chi(+1) = 1, chi(-1) = -1).  Since full-support words
of a linear code come in pairs +-x with equal sign (n even), Index(D) is even, so Index(C') = -Index(C) mod 8; any two
self-dual codes of length n are joined by a chain of neighbors, and 3k copies of the tetracode have index 0, so 8 | Index.

The test walks random neighbor chains from the tetracode sum: at each step pick a random isotropic v outside C,
D = {c in C : c.v = 0}, C' = D + <v>; it checks self-duality, computes the three indices by brute force over an
information set (code/verify.py brute_index) and checks the identity and the congruences.
usage: python neighbor_index_test.py --n 24 --steps 200 --seed 1
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))              # code/ of the online note (verify.py)
import verify as V

TETRA = np.array([[1, 0, 1, 1], [0, 1, 1, 2]], dtype=np.int64)


def tetracode_sum(n):
    m = n // 4
    G = np.zeros((2 * m, n), dtype=np.int64)
    for j in range(m):
        G[2 * j:2 * j + 2, 4 * j:4 * j + 4] = TETRA
    return G


def index(G):
    r = V.brute_index(V.rref(np.asarray(G) % 3)[0], max_dim=26)
    assert r is not None
    return r[0]


def neighbor(G, rng):
    """random neighbor: returns (D generator, C' generator)"""
    n = G.shape[1]
    while True:
        v = rng.integers(0, 3, n)
        if (v @ v) % 3 or not ((G @ v) % 3).any():
            continue                                       # need v isotropic and v not in C (= C-perp)
        break
    # D = kernel of c -> c.v on C: coefficients a with a.(G v) = 0
    w = (G @ v) % 3                                        # functional on the row space coordinates
    N = V.left_nullspace(w.reshape(-1, 1))                 # vectors a with a.w = 0
    D = V.rref((N @ G) % 3)[0]
    Cp = V.rref(np.vstack([D, v]) % 3)[0]
    return D, Cp


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=24); ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    G = tetracode_sum(a.n); I = index(G)
    rec = {"n": a.n, "steps": a.steps, "identity_failures": 0, "mod8_failures": 0, "mod24_failures": 0, "indices": {}}
    for s in range(a.steps):
        D, Cp = neighbor(G, rng)
        d, n = Cp.shape
        assert V.rank3(Cp) == n // 2 and not ((Cp @ Cp.T) % 3).any(), "C' not self-dual"
        assert V.rank3(D) == n // 2 - 1
        ID, Ip = index(D), index(Cp)
        if I + Ip != 4 * ID:
            rec["identity_failures"] += 1
            print("identity FAILS", s, I, Ip, ID)
        if Ip % 8:
            rec["mod8_failures"] += 1
        if Ip % 24:
            rec["mod24_failures"] += 1
        rec["indices"][str(Ip)] = rec["indices"].get(str(Ip), 0) + 1
        G, I = Cp, Ip
    print(json.dumps(rec))
    if a.out:
        Path(a.out).write_text(json.dumps(rec, indent=1) + "\n")


if __name__ == "__main__":
    main()
