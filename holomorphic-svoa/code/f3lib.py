# Copied from claude/f3lib.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Small, self-contained F_3 toolkit used by the Claude audit.

Written from scratch for the audit: it does not import anything from
``scft_code_search`` or ``benchmarks``, so that every check in this folder is
an independent re-implementation of the repository's checks.

Conventions
-----------
* Vectors and matrices are numpy int arrays with entries in {0,1,2}.
* A signed coordinate map ("monomial map") is a pair (perm, signs) meaning
  e_i -> signs[i] * e_{perm[i]}, i.e. (x g)_{perm[i]} = signs[i] * x_i.
  This is the same convention as ``apply_monomial`` in the repository; it was
  re-derived here by checking which convention preserves the saved code.
* A codeword is *maximal* if all coordinates are nonzero.  Its parity is the
  number of coordinates equal to 2, mod 2.  Index(C) = #even - #odd
  (Gaiotto--Johnson-Freyd, Definition 1.1).
"""
from __future__ import annotations

import itertools
import os
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# linear algebra over F_3
# ----------------------------------------------------------------------------
def rref(matrix):
    """Reduced row echelon form over F_3.  Returns (R, pivot_columns)."""
    a = np.array(matrix, dtype=np.int64) % 3
    rows, cols = a.shape
    pivots = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(a[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            a[[r, p]] = a[[p, r]]
        if a[r, c] == 2:          # 2^{-1} = 2 in F_3
            a[r] = (2 * a[r]) % 3
        col = a[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            a[nzr] = (a[nzr] - np.outer(col[nzr], a[r])) % 3
        pivots.append(c)
        r += 1
    return a[:r], pivots


def rank(matrix):
    return len(rref(matrix)[1])


def nullspace(matrix):
    """Basis (as rows) of {x : matrix @ x = 0} over F_3."""
    a = np.array(matrix, dtype=np.int64) % 3
    n = a.shape[1]
    r, piv = rref(a)
    free = [j for j in range(n) if j not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64)
        v[f] = 1
        for i, p in enumerate(piv):
            v[p] = (-r[i, f]) % 3
        basis.append(v)
    return np.array(basis, dtype=np.int64).reshape(len(basis), n)


def is_self_dual(g):
    g = np.asarray(g) % 3
    k, n = g.shape
    return n == 2 * k and rank(g) == k and not np.any((g @ g.T) % 3)


def systematic(g):
    """Return (A, order) with the columns of g permuted by ``order`` row
    equivalent to [I | A]."""
    r, piv = rref(g)
    k, n = r.shape
    rest = [j for j in range(n) if j not in piv]
    order = list(piv) + rest
    return r[:, rest] % 3, order


# ----------------------------------------------------------------------------
# signed coordinate maps
# ----------------------------------------------------------------------------
def monomial_matrix(perm, signs):
    n = len(perm)
    m = np.zeros((n, n), dtype=np.int64)
    for i, (p, s) in enumerate(zip(perm, signs)):
        m[i, p] = s % 3
    return m


def preserves(g, m):
    """Self-dual C = rowspace(g) is preserved by right multiplication with m
    iff (g m) g^T = 0."""
    return not np.any(((np.asarray(g) @ m) % 3) @ np.asarray(g).T % 3)


def order_of(m, limit=1000):
    n = m.shape[0]
    p = np.eye(n, dtype=np.int64)
    for k in range(1, limit + 1):
        p = (p @ m) % 3
        if np.array_equal(p, np.eye(n, dtype=np.int64)):
            return k
    raise ValueError("order exceeds limit")


def perm_and_signs(m):
    """Inverse of monomial_matrix."""
    perm = [int(np.nonzero(row)[0][0]) for row in m]
    signs = [int(row[p]) for row, p in zip(m, perm)]
    return perm, signs


def cycles(perm):
    seen, out = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        c, j = [], i
        while j not in seen:
            seen.add(j)
            c.append(j)
            j = perm[j]
        out.append(c)
    return out


def signed_cycle_data(perm, signs):
    """List of (cycle length, product of signs as +/-1) for each cycle."""
    out = []
    for c in cycles(perm):
        prod = 1
        for i in c:
            prod *= 1 if signs[i] == 1 else -1
        out.append((len(c), prod))
    return out


def frame_shape(perm, signs):
    """Frame shape of the real 72x72 signed permutation matrix:
    a cycle of length l with sign product +1 contributes l^1, with -1
    contributes (2l)^1 l^-1 (char. poly x^l + 1 = (x^{2l}-1)/(x^l-1))."""
    shape = {}
    for l, s in signed_cycle_data(perm, signs):
        if s == 1:
            shape[l] = shape.get(l, 0) + 1
        else:
            shape[2 * l] = shape.get(2 * l, 0) + 1
            shape[l] = shape.get(l, 0) - 1
    return {k: v for k, v in sorted(shape.items()) if v}


# ----------------------------------------------------------------------------
# enumeration of small codes and signed maximal-word counts
# ----------------------------------------------------------------------------
def all_codewords(basis, chunk=1 << 15):
    """Yield arrays of codewords (rows) of the span of ``basis``."""
    basis = np.asarray(basis, dtype=np.int64) % 3
    d = basis.shape[0]
    total = 3 ** d
    powers = 3 ** np.arange(d, dtype=np.int64)
    for start in range(0, total, chunk):
        ids = np.arange(start, min(total, start + chunk), dtype=np.int64)
        coeff = (ids[:, None] // powers) % 3
        yield (coeff @ basis) % 3


def maximal_counts(basis):
    """(#even, #odd) maximal codewords in span(basis)  (exhaustive)."""
    even = odd = 0
    for words in all_codewords(basis):
        alive = np.all(words != 0, axis=1)
        par = np.sum(words == 2, axis=1) % 2
        even += int(np.count_nonzero(alive & (par == 0)))
        odd += int(np.count_nonzero(alive & (par == 1)))
    return even, odd


def index_small(g):
    e, o = maximal_counts(g)
    return e - o


def fixed_subcode(g, m):
    """Basis of {x in C : x m = x} for self-dual C = rowspace(g)."""
    g = np.asarray(g) % 3
    n = g.shape[1]
    # x in C  <=> g x^T = 0 (self-duality);  x m = x <=> (m^T - I) x^T = 0
    cons = np.vstack([g, (m.T - np.eye(n, dtype=np.int64)) % 3]) % 3
    return nullspace(cons)


# ----------------------------------------------------------------------------
# exhaustive full-support enumeration in C  (2^k Gray code)
# ----------------------------------------------------------------------------
_C_SOURCE = HERE / "index_gray.c"
_C_BINARY = HERE / "build" / "index_gray"


def _ensure_binary():
    _C_BINARY.parent.mkdir(exist_ok=True)
    if (not _C_BINARY.exists()
            or _C_BINARY.stat().st_mtime < _C_SOURCE.stat().st_mtime):
        subprocess.run(["cc", "-O3", "-march=native", "-o", str(_C_BINARY),
                        str(_C_SOURCE)], check=True)
    return _C_BINARY


def index_exhaustive(g, workers=1, prefix_bits=6, tag="code"):
    """Exact Index(C) by enumerating all 2^k full-support information words.

    Every maximal codeword has nonzero information coordinates, so the
    2^k words m in {1,2}^k with codeword (m, mA) cover all maximal words.
    ``workers`` independent processes split the leading ``prefix_bits``.
    Returns dict(even, odd, index, seconds).
    """
    import time
    a, _ = systematic(g)
    k, r = a.shape
    assert k <= 64 and r <= 64
    binary = _ensure_binary()
    tmp = HERE / "build" / f"{tag}_A.txt"
    with open(tmp, "w") as fh:
        fh.write(f"{k} {r}\n")
        for row in a:
            fh.write(" ".join(str(int(x)) for x in row) + "\n")
    t0 = time.time()
    procs = [subprocess.Popen(["nice", "-n", "19", str(binary), str(tmp),
                               str(prefix_bits), str(w), str(workers)],
                              stdout=subprocess.PIPE, text=True)
             for w in range(workers)]
    even = odd = 0
    for p in procs:
        out, _ = p.communicate()
        assert p.returncode == 0, out
        e, o = map(int, out.split())
        even += e
        odd += o
    return {"even": even, "odd": odd, "index": even - odd,
            "seconds": time.time() - t0}


# ----------------------------------------------------------------------------
# a few reference codes
# ----------------------------------------------------------------------------
def golay12():
    """Extended ternary Golay code [I | S] with S the Paley conference
    matrix of order 6 (S S^T = 5 I = -I mod 3)."""
    q = [1, 4]                      # quadratic residues mod 5
    chi = lambda x: 0 if x % 5 == 0 else (1 if x % 5 in q else 2)
    s = np.zeros((6, 6), dtype=np.int64)
    for i in range(6):
        for j in range(6):
            if i == 0 and j == 0:
                s[i, j] = 0
            elif i == 0:
                s[i, j] = 1
            elif j == 0:
                s[i, j] = 1          # symmetric conference matrix (5 = 1 mod 4)
            else:
                s[i, j] = chi(j - i)
    g = np.hstack([np.eye(6, dtype=np.int64), s % 3])
    assert is_self_dual(g)
    return g


BALANCED12 = np.array([
    [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
    [0, 1, 2, 0, 1, 2, 0, 2, 1, 0, 0, 0],
    [0, 0, 0, 0, 1, 2, 0, 1, 2, 0, 2, 1],
], dtype=np.int64)
"""The [12,6] code B of the note (same rows as reference_seeds.py).
It is the self-dual code with four (111)-components glued by the tetracode."""


def tetracode():
    g = np.array([[1, 0, 1, 1], [0, 1, 1, 2]], dtype=np.int64)
    assert is_self_dual(g)
    return g


def direct_sum(*gs):
    k = sum(g.shape[0] for g in gs)
    n = sum(g.shape[1] for g in gs)
    out = np.zeros((k, n), dtype=np.int64)
    r = c = 0
    for g in gs:
        out[r:r + g.shape[0], c:c + g.shape[1]] = g
        r += g.shape[0]
        c += g.shape[1]
    return out
