"""Check the sparse certificate of Proposition 9 of the note on GJF Conjecture 1.4 (index.html) from the JSON file alone (plain Python, no numpy).

The file golay_chain_transfer.json (same folder) lists the non-zero entries of the 81 x 81 transfer matrix T of Lemma 8
(row, column, value; the state ordering is in the file), the vector s and the non-zero entries of u.  This script checks,
in exact integer arithmetic, that u = T s, that T u = -u and that s . u = -24, and as a consequence that
s^T T^k s = (-1)^k 24 for k = 1..20.  That the matrix in the file is the T of Lemma 8, built from the Golay code and the
vectors p, r, v, w, is checked separately by golay_chain_check.py, which writes the file.

Usage:  python3 check_golay_chain_certificate.py [path/to/golay_chain_transfer.json]
Exit status 1 if a check fails.  Licence: MIT (see LICENSE).
"""
import json, sys
from pathlib import Path

path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "golay_chain_transfer.json"   # copy for the standalone note: certificate in the same folder
cert = json.loads(path.read_text())
N = len(cert["state_ordering"]["labels_c1c2f1f2"])
T = {}
for i, j, value in cert["T_nonzero_entries_row_col_value"]:
    assert 0 <= i < N and 0 <= j < N and isinstance(value, int) and value != 0 and (i, j) not in T
    T[(i, j)] = value
s = cert["s"]
u = [0] * N
for i, value in cert["u_nonzero_entries_index_value"]:
    u[i] = value


def apply(x):
    y = [0] * N
    for (i, j), t in T.items():
        y[i] += t * x[j]
    return y


checks = {
    "states": N == 81 and len(s) == 81,
    "number of non-zero entries of T": len(T) == cert["T_nonzero_count"],
    "u = T s": apply(s) == u,
    "T u = -u": apply(u) == [-x for x in u],
    "s . u = -24": sum(a * b for a, b in zip(s, u)) == -24,
}
x, powers = s, []
for k in range(1, 21):
    x = apply(x)
    powers.append(sum(a * b for a, b in zip(s, x)))
checks["s^T T^k s = (-1)^k 24 for k = 1..20"] = powers == [(-1) ** k * 24 for k in range(1, 21)]
for name, ok in checks.items():
    print(("ok    " if ok else "FAIL  ") + name)
print(f"{path.name}: {len(T)} non-zero entries of T, {sum(1 for t in u if t)} of u")
sys.exit(0 if all(checks.values()) else 1)
