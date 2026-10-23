"""Independent F_3 code utilities for refereeing the 8 | Index argument.
Nothing imported from the project repo."""
import itertools
import numpy as np

P = 3


def rref(M):
    M = np.array(M, dtype=np.int64) % P
    if M.ndim == 1:
        M = M.reshape(1, -1)
    rows, cols = M.shape
    piv = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            M[[r, p]] = M[[p, r]]
        if M[r, c] == 2:
            M[r] = (2 * M[r]) % P
        for i in range(rows):
            if i != r and M[i, c] != 0:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        piv.append(c)
        r += 1
    return M[:r], piv


def rank(M):
    if len(M) == 0:
        return 0
    return rref(M)[0].shape[0]


def dual(G, n):
    """Basis of {y : G y = 0} = rowspace(G)^perp."""
    G = np.asarray(G)
    if G.size == 0:
        return np.eye(n, dtype=np.int64)
    R, piv = rref(G)
    free = [c for c in range(n) if c not in piv]
    B = np.zeros((len(free), n), dtype=np.int64)
    for j, f in enumerate(free):
        B[j, f] = 1
        for i, p in enumerate(piv):
            B[j, p] = (-R[i, f]) % P
    return B


def inv_mod3(A):
    m = A.shape[0]
    R, piv = rref(np.hstack([A % P, np.eye(m, dtype=np.int64)]))
    assert piv[:m] == list(range(m)), "singular"
    return R[:, m:] % P


def all_words(G):
    G = np.asarray(G) % P
    k = G.shape[0]
    U = np.array(list(itertools.product(range(3), repeat=k)), dtype=np.int64)
    return (U @ G) % P


def chi_prod(X):
    """prod_i chi(x_i) for each row, chi(0)=0, chi(1)=1, chi(2)=-1."""
    X = np.asarray(X)
    full = np.all(X != 0, axis=1)
    sgn = 1 - 2 * (np.sum(X == 2, axis=1) % 2)
    return np.where(full, sgn, 0)


def index_brute(G):
    """Sum over ALL 3^k codewords (no shortcuts); basis via rref, and the
    enumerated words are checked to be distinct."""
    R = rref(G)[0]
    X = all_words(R)
    assert len({tuple(r) for r in X.tolist()}) == X.shape[0]
    return int(chi_prod(X).sum())


def index_fast(G, chunk_bits=20):
    """Only words that are nonzero on an information set can be full support:
    enumerate u in {1,2}^k on the pivot set, 2^k words."""
    R, piv = rref(G)
    k, n = R.shape
    if k == 0:
        return 0 if n > 0 else 1
    nonpiv = [c for c in range(n) if c not in piv]
    B = R[:, nonpiv].astype(np.int64)  # k x (n-k)
    total = 0
    N = 1 << k
    step = 1 << min(chunk_bits, k)
    shifts = np.arange(k, dtype=np.int64)
    for start in range(0, N, step):
        idx = np.arange(start, min(N, start + step), dtype=np.int64)
        bits = (idx[:, None] >> shifts[None, :]) & 1  # 1 -> coordinate value 2
        u = 1 + bits
        su = 1 - 2 * (bits.sum(axis=1) % 2)
        if B.shape[1] == 0:
            total += int(su.sum())
            continue
        X = (u @ B) % P
        full = np.all(X != 0, axis=1)
        sx = 1 - 2 * (np.sum(X == 2, axis=1) % 2)
        total += int((su * sx * full).sum())
    return total


def tetracode():
    return np.array([[1, 0, 1, 1], [0, 1, 1, 2]], dtype=np.int64)


def tetra_sum(n):
    T = tetracode()
    G = np.zeros((n // 2, n), dtype=np.int64)
    for b in range(n // 4):
        G[2 * b:2 * b + 2, 4 * b:4 * b + 4] = T
    return G


def is_self_dual(G, n):
    G = np.asarray(G) % P
    return rank(G) == n // 2 and not ((G @ G.T) % P).any()


def random_neighbor(G, rng):
    """Kneser neighbour: v with v.v=0, v not in C; D = C cap v^perp; C' = D + <v>.
    Returns (G', D)."""
    G = np.asarray(G) % P
    m, n = G.shape
    while True:
        v = rng.integers(0, 3, n)
        if (v @ v) % P != 0:
            continue
        s = (G @ v) % P
        if not s.any():
            continue
        coeff = dual(s.reshape(1, -1), m)  # a with a.s = 0
        D = (coeff @ G) % P
        Gp = np.vstack([D, v.reshape(1, -1)]) % P
        return rref(Gp)[0], rref(D)[0]


def random_monomial(G, rng):
    n = G.shape[1]
    perm = rng.permutation(n)
    signs = rng.choice([1, 2], n)
    return (G[:, perm] * signs[None, :]) % P


def random_orth(m, rng, steps=None):
    """Product of random reflections in O(m, F_3) for the standard form."""
    if steps is None:
        steps = 4 * m + 10
    Q = np.eye(m, dtype=np.int64)
    I = np.eye(m, dtype=np.int64)
    for _ in range(steps):
        while True:
            v = rng.integers(0, 3, m)
            a = int(v @ v) % P
            if a != 0:
                break
        ainv = a  # 1^-1=1, 2^-1=2 mod 3
        S = (I - 2 * ainv * np.outer(v, v)) % P
        Q = (Q @ S) % P
    return Q


def random_selfdual_AAT(n, rng):
    """[I | A0 Q] with A0 A0^T = -I (blocks [[1,1],[1,2]]) and Q orthogonal."""
    m = n // 2
    assert m % 2 == 0
    A0 = np.zeros((m, m), dtype=np.int64)
    for b in range(m // 2):
        A0[2 * b:2 * b + 2, 2 * b:2 * b + 2] = [[1, 1], [1, 2]]
    A = (A0 @ random_orth(m, rng)) % P
    assert not ((A @ A.T + np.eye(m, dtype=np.int64)) % P).any()
    G = np.hstack([np.eye(m, dtype=np.int64), A]) % P
    # random column permutation so the information set is not always the first half
    return G
