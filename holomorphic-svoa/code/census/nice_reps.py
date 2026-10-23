# Copied from claude/census/nice_reps.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Find the most structured representative of a seed's PGL(2,11)-orbit.
Score: blocks that are arithmetic progressions (mod 11; a block containing inf scores by its
3 finite points), blocks that are cosets of a subgroup/translates of {0,1,2,3}, etc."""
import itertools, json, sys
import numpy as np
from symmetry import mobius, signs_for, act
from hpcell import evaluate

def pgl_elements():
    seen = {}
    for a, b, c, d in itertools.product(range(11), repeat=4):
        if (a * d - b * c) % 11 == 0:
            continue
        perm = tuple(mobius((a, b, c, d), k) for k in range(12))
        if perm not in seen:
            seen[perm] = (signs_for(perm), (a, b, c, d))
    return seen

def elem(k):
    return None if k == 0 else k % 11

def is_ap(pts):
    f = sorted(elem(k) for k in pts if k != 0)
    if len(f) < 3:
        return True
    for d in range(1, 11):
        for a in f:
            if {(a + i * d) % 11 for i in range(len(f))} == set(f):
                return True
    return False

def describe(cfg, R):
    s = []
    for b, q in cfg:
        pts = ["inf" if k == 0 else str(k % 11) for k in b]
        s.append("{" + ",".join(pts) + "}" + ("-" if q else "+"))
    return " ".join(s) + "  R={" + ",".join("inf" if k == 0 else str(k % 11) for k in R) + "}"

G = pgl_elements()
seeds = [("9-state A", [[0,1,2,3],[4,5,6,9],[7,8,10,11]], [0,1,0], [2,7]),
         ("9-state B", [[0,1,2,3],[4,5,6,9],[7,8,10,11]], [0,1,0], [2,10]),
         ("27-state consecutive", [[0,1,2,3],[4,5,6,7],[8,9,10,11]], [0,0,0], [3,5,8])]
for name, blocks, pars, R in seeds:
    images = {}
    for perm, ((d, e), M) in G.items():
        cfg, R2 = act(perm, d, blocks, pars, R)
        images[(cfg, R2)] = M
    print(f"==== {name}: orbit of (config, R) under PGL(2,11) has {len(images)} elements")
    scored = []
    for (cfg, R2), M in images.items():
        n_ap = sum(is_ap(b) for b, q in cfg)
        n_even = sum(1 for b, q in cfg if q == 0)
        inf_in_R = 0 in R2
        scored.append((-n_ap, -n_even, not inf_in_R, describe(cfg, R2), M))
    scored.sort()
    for s in scored[:6]:
        print("  APs", -s[0], "even blocks", -s[1], "|", s[3], " via", s[4])
