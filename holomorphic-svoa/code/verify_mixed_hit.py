# Copied from claude/kspan/verify_mixed_hit.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Independent verification of a mixed-width product-ring hit from mixed_products.py.

Rebuilds the systematic generator G = [I | A] of the ring code directly from the stored cells
(legacy form: Q_j orthogonal, S_j anti-orthogonal, the first r outputs of a_j Q_j delayed into
check block j+1), checks self-duality and rank with f3lib, and recomputes the index with the
generator-only banded transfer DP used for Codex's certificates (copied from
verify_codex_cases.banded_index so that nothing from cells_general.table is reused).

usage: python verify_mixed_hit.py results/kspan/mixed_k13_run1.json [hit_no] [--cert results/kspan/k13_certificate.json]
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from f3lib import is_self_dual, rank


def phi_rows(w):
    return np.where(np.all(w != 0, axis=1), 1 - 2 * (np.sum(w == 2, axis=1) % 2), 0)


def banded_index(G, widths):
    """Tr(T_0 ... T_{L-1}) from the generator alone (block-bidiagonal A, wrap sign inside A)."""
    d = G.shape[0]; A = G[:, d:] % 3; L = len(widths)
    off = np.concatenate([[0], np.cumsum(widths)])
    assert off[-1] == d
    blk = lambda i, j: A[off[i]:off[i + 1], off[j]:off[j + 1]]
    for i in range(L):
        for j in range(L):
            if j not in (i, (i + 1) % L):
                assert not blk(i, j).any(), (i, j)
    ins = [np.array(list(itertools.product((1, 2), repeat=w)), dtype=np.int64) for w in widths]
    phis = [phi_rows(a) for a in ins]
    keys = []
    for j in range(L):
        img = (ins[(j - 1) % L] @ blk((j - 1) % L, j)) % 3
        ks = sorted({tuple(x) for x in img}); keys.append((ks, {k: i for i, k in enumerate(ks)}))
    P = None
    for j in range(L):
        kin, _ = keys[j]; kout, posout = keys[(j + 1) % L]
        out = np.array([posout[tuple(x)] for x in (ins[j] @ blk(j, (j + 1) % L)) % 3])
        base = (ins[j] @ blk(j, j)) % 3
        T = np.zeros((len(kin), len(kout)), dtype=object)
        for i, s in enumerate(kin):
            w = phi_rows((base + np.array(s)) % 3) * phis[j]
            for kk in range(len(kout)):
                T[i, kk] = int(w[out == kk].sum())
        P = T if P is None else P.dot(T)
    return int(np.trace(P))


def generator_from_cells(cells, r, cyclic):
    widths = [len(c["Q"]) for c in cells]
    d = sum(widths); off = np.concatenate([[0], np.cumsum(widths)]); L = len(cells)
    A = np.zeros((d, d), dtype=np.int64)
    for j, c in enumerate(cells):
        Q = np.array(c["Q"]) % 3; S = np.array(c["S"]) % 3
        Sn = np.array(cells[(j + 1) % L]["S"]) % 3
        A[off[j]:off[j + 1], off[j]:off[j + 1]] = (Q[:, r:] @ S[r:, :]) % 3
        nxt = (Q[:, :r] @ Sn[:r, :]) % 3
        if j == L - 1 and not cyclic:
            nxt = (-nxt) % 3
        A[off[j]:off[j + 1], off[(j + 1) % L]:off[(j + 1) % L + 1]] = nxt
    return np.hstack([np.eye(d, dtype=np.int64), A]) % 3, widths


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cert = sys.argv[sys.argv.index("--cert") + 1] if "--cert" in sys.argv else None
    data = json.load(open(args[0])); hit = data["hits"][int(args[1]) if len(args) > 1 else 0]
    r = int(data.get("r", 3)); cyclic = hit["kind"] == "cyc"
    G, widths = generator_from_cells(hit["cells"], r, cyclic)
    rep = {"k": hit["k"], "length": int(G.shape[1]), "widths": widths, "kind": hit["kind"], "claimed_index": hit["index"],
           "rank": int(rank(G)), "self_dual": bool(is_self_dual(G))}
    rep["index_banded_DP"] = banded_index(G, widths)
    rep["ok"] = rep["self_dual"] and rep["rank"] == G.shape[0] and abs(rep["index_banded_DP"]) == 24 and rep["index_banded_DP"] == rep["claimed_index"]
    print(json.dumps(rep))
    if cert and rep["ok"]:
        Path(cert).write_text(json.dumps({
            "what": f"self-dual ternary code of length {G.shape[1]} (c = {G.shape[1]}, k = {hit['k']}) with index {rep['index_banded_DP']}; "
                    f"gcd({hit['k']},24) = 1 so no orbifold is needed: the lattice SVOA itself has the TMF-minimal index 24",
            "construction": "mixed-width product ring: cells of widths " + ",".join(map(str, widths)) + f", common delay count r = {r}, "
                            + ("cyclic" if cyclic else "negacyclic") + " wrap; code = {(a | aA)}, A block-bidiagonal",
            "cells": hit["cells"], "r": r, "cyclic": cyclic, "generator": G.tolist(), "index": rep["index_banded_DP"],
            "verification": rep, "found_by": "kspan/mixed_products.py (seed 2, pool 60, 1,236,959 products)",
            "verified_by": "kspan/verify_mixed_hit.py (generator rebuilt from cells; f3lib self-duality; generator-only banded DP)"}, indent=0))
        print("certificate written:", cert)


if __name__ == "__main__":
    main()
