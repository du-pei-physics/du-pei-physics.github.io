# Copied from claude/k8/anomaly/neighbour_chain.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Exact checks for S6 of NOTES_q8_anomaly.md.

(1) Z^96 = Lambda(C') for the Q8-invariant code C' = tetracode^{+24} w.r.t. the frame F = L^24 * u,
    u = 1+i+j (Lipschitz quaternions, left multiplication by Q8 = certificate signed permutations).
(2) An explicit Q8-equivariant neighbour chain from the certificate code C to C'.
All arithmetic exact (integers / F_3).
"""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
from fractions import Fraction

CERT = HERE.parent / "certificates" / "original" / 'q8_certificate.json'
d = json.load(open(CERT))
n = 96

# ---------- signed permutations from the certificate ----------
def mat(sym):
    p = sym['permutation']; sg = [1 if x == 1 else -1 for x in sym['signs']]
    return (tuple(p), tuple(sg))           # e_a -> sg[a] e_{p[a]}
def act(g, v, mod=None):
    p, s = g; w = [0]*n
    for a in range(n):
        w[p[a]] += s[a]*v[a]
    return [x % mod for x in w] if mod else w
I = mat(d['symmetry_i']); J = mat(d['symmetry_j'])

# ---------- F_3 linear algebra ----------
def rref(rows):
    M = [list(r) for r in rows]; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][c] % 3), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = 1 if M[r][c] % 3 == 1 else 2
        M[r] = [(x*inv) % 3 for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % 3:
                f = M[i][c]; M[i] = [(M[i][t] - f*M[r][t]) % 3 for t in range(n)]
        piv.append(c); r += 1
    return M[:r], piv
class Code:
    def __init__(self, rows):
        self.R, self.piv = rref(rows)
    def dim(self): return len(self.R)
    def contains(self, v):
        v = [x % 3 for x in v]
        for row, c in zip(self.R, self.piv):
            if v[c]: f = v[c]; v = [(v[t] - f*row[t]) % 3 for t in range(n)]
        return not any(v)
    def dual(self):
        # kernel of R: free columns
        free = [c for c in range(n) if c not in self.piv]
        rows = []
        for f in free:
            v = [0]*n; v[f] = 1
            for row, c in zip(self.R, self.piv):
                v[c] = (-row[f]) % 3
            rows.append(v)
        return Code(rows)
    def is_self_dual(self):
        return self.dim() == n//2 and all(sum(a*b for a, b in zip(r1, r2)) % 3 == 0 for r1 in self.R for r2 in self.R)
    def is_invariant(self, g):
        return all(self.contains(act(g, r, 3)) for r in self.R)
    def intersect(self, other):
        # C1 ∩ C2 = (C1^⊥ + C2^⊥)^⊥
        return Code(self.dual().R + other.dual().R).dual()
    def plus(self, other):
        return Code(self.R + other.R)
    def meet_perp(self, U):
        # C ∩ U^⊥  = (C^⊥ + U)^⊥
        return Code(self.dual().R + U.R).dual()

C = Code(d['generator'])
assert C.dim() == 48 and C.is_self_dual() and C.is_invariant(I) and C.is_invariant(J)
print("certificate code: self-dual, Q8-invariant, dim 48")

# ---------- (1) the frame code of Z^96 ----------
# Block coordinates (1,i,j,k) = (0,1,2,3). Right multiplication by u = 1+i+j on a block:
# f_1 = 1*u = 1+i+j, f_i = i*u = -1 + i + k, f_j = j*u = -1 + j - k, f_k = k*u = -i + j + k
frame_block = [[1,1,1,0],[-1,1,0,1],[-1,0,1,-1],[0,-1,1,1]]
# check: norms 3, orthogonal, and Q8-equivariance f_{q e_a} = q f_a on block 0
for a in range(4):
    assert sum(x*x for x in frame_block[a]) == 3
    for b in range(a+1, 4):
        assert sum(x*y for x, y in zip(frame_block[a], frame_block[b])) == 0
frame = []
for blk in range(24):
    for a in range(4):
        v = [0]*n
        for b in range(4): v[4*blk + b] = frame_block[a][b]
        frame.append(v)
for g in (I, J):
    for a in range(n):
        gf = act(g, frame[a])
        # g e_a = s e_{p(a)}  =>  g f_a = s f_{p(a)}
        p, s = g
        assert gf == [s[a]*x for x in frame[p[a]]], "frame not Q8-equivariant"
print("frame F = L^24 u: 96 orthogonal norm-3 vectors, Q8-equivariant (g f_a = s_a f_{p(a)})")
# code of Z^96 w.r.t. frame: coordinates of x in frame basis are <x,f_a>/3; x in Z^96 gives <x,f_a> in Z;
# C' = { (<x,f_a>)_a mod 3 : x in Z^96 }  (since Lambda(C') = {y/sqrt3 : y mod 3 in C'} in frame coords)
rows = []
for b in range(n):
    e = [0]*n; e[b] = 1
    rows.append([sum(e[t]*frame[a][t] for t in range(n)) % 3 for a in range(n)])
Cp = Code(rows)
print("frame code C': dim", Cp.dim(), "self-dual", Cp.is_self_dual())
# Q8-invariance in frame coordinates: the Q8 action on frame coordinates is the same signed permutation
print("C' Q8-invariant (same signed permutations):", Cp.is_invariant(I), Cp.is_invariant(J))
# block structure: is C' = direct sum of 24 copies of a length-4 code? and is that the tetracode?
blk = Code([r[:4] + [0]*(n-4) for r in rows[:4]])
print("block-0 code dim:", blk.dim(), "rows:", [r[:4] for r in blk.R])
tet_words = [v for v in blk.R]
# tetracode check: 2-dim, self-dual, all nonzero words weight 3
import itertools
words = set()
for a in range(3):
    for b in range(3):
        w = tuple((a*blk.R[0][t] + b*blk.R[1][t]) % 3 for t in range(4)); words.add(w)
wts = sorted(sum(1 for x in w if x) for w in words if any(w))
print("block code nonzero weights:", wts, "(tetracode iff all 3)")
print("C' = direct sum of blocks:", Cp.dim() == 48 and all(Cp.contains(r) for r in rows))
# C ∩ C'
CC = C.intersect(Cp)
print("dim(C ∩ C') =", CC.dim())
# min weight of C restricted: does C contain any weight-3 word (would be a block-tetracode word)?
has_w3 = any(C.contains(w + [0]*0) for w in [list(r) for r in Cp.R] )  # tetracode words are weight 3
print("C contains some block tetracode generator:", has_w3)

# ---------- (2) equivariant neighbour chain ----------
chain = [C]; steps = []
Ci = C
for b in range(24):
    U = Code([r for r in Cp.R if any(r[4*b:4*b+4]) and not any(r[:4*b]) and not any(r[4*b+4:])])
    # U must be the 2-dim block tetracode
    assert U.dim() == 2 and U.is_invariant(I) and U.is_invariant(J)
    if all(Ci.contains(r) for r in U.R):
        steps.append((b, 'skip (U already in C_i)')); continue
    inter = Ci.intersect(U).dim()
    Cn = Ci.meet_perp(U).plus(U)
    ok = Cn.is_self_dual() and Cn.is_invariant(I) and Cn.is_invariant(J)
    steps.append((b, f"U∩C_i dim {inter}; new code self-dual&Q8-inv: {ok}; dim(C_new ∩ C') = {Cn.intersect(Cp).dim()}"))
    assert ok and inter == 0
    Ci = Cn; chain.append(Cn)
for st in steps: print("block", st[0], ":", st[1])
print("final code == C':", Ci.dim() == 48 and all(Cp.contains(r) for r in Ci.R) and all(Ci.contains(r) for r in Cp.R))
print("number of orbifold steps:", len(chain) - 1)
# record the chain
out = {'frame_block': frame_block, 'steps': steps, 'n_steps': len(chain)-1,
       'intermediate_dims_C_i_cap_Cprime': [c.intersect(Cp).dim() for c in chain]}
json.dump(out, open(HERE / 'out' / 'neighbour_chain_result.json', 'w'), indent=1)
