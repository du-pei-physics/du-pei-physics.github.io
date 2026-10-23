"""Base matrices S0 with S0 S0^T = -I over F_3 (q = 12), for structured cells."""
import numpy as np
from core import A0


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def paley12():
    """Paley conference matrix of order 12 (from QR mod 11): P P^T = 11 I = -I mod 3."""
    p = np.zeros((12, 12), dtype=np.int64)
    for i in range(12):
        for j in range(12):
            p[i, j] = 0 if i == j == 0 else (1 if i == 0 else (-1 if j == 0 else legendre(j - i, 11)))
    return p % 3


def conf6x2():
    """Two copies of the symmetric conference matrix of order 6 (QR mod 5): C C^T = 5 I = -I."""
    c6 = np.zeros((6, 6), dtype=np.int64)
    for i in range(6):
        for j in range(6):
            c6[i, j] = 0 if i == j == 0 else (1 if (i == 0 or j == 0) else legendre(j - i, 5))
    g = np.zeros((12, 12), dtype=np.int64)
    g[:6, :6] = c6; g[6:, 6:] = c6
    return g % 3


def bases12():
    a0 = A0(12)
    a0d = a0.copy(); a0d[:, 0] = (-a0d[:, 0]) % 3
    out = {"A0": a0, "A0D": a0d, "Conf6x2": conf6x2(), "Paley12": paley12()}
    for k, m in out.items():
        assert np.array_equal((m @ m.T) % 3, (2 * np.eye(12, dtype=np.int64)) % 3), k
    return out


BASES12 = bases12()
