# Copied from claude/k8/anomaly/psi_cocycle_check.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
import itertools
E=[(a,b) for a in range(3) for b in range(3)]
add=lambda x,y:((x[0]+y[0])%3,(x[1]+y[1])%3)
det=lambda x,y:(x[0]*y[1]-x[1]*y[0])%3
def psi(k):  # exponent of zeta
    return lambda x,y:(k*det(x,y))%3
mats=[((a,b),(c,d)) for a,b,c,d in itertools.product(range(3),repeat=4) if (a*d-b*c)%3]
act=lambda M,x:((M[0][0]*x[0]+M[0][1]*x[1])%3,(M[1][0]*x[0]+M[1][1]*x[1])%3)
ok_cocycle=all((p(b,c)-p(add(a,b),c)+p(a,add(b,c))-p(a,b))%3==0 for k in range(3) for p in [psi(k)] for a in E for b in E for c in E)
print("cocycle:",ok_cocycle, "|GL(2,3)|=",len(mats))
for M in mats:
    d=(M[0][0]*M[1][1]-M[0][1]*M[1][0])%3
    inv=all(psi(1)(act(M,x),act(M,y))==psi(1)(x,y) for x in E for y in E)
    neg=all(psi(1)(act(M,x),act(M,y))==psi(2)(x,y) for x in E for y in E)
    assert (d==1 and inv) or (d==2 and neg)
print("SL(2,3) fixes psi_k on the nose; det=-1 sends k->-k")
# commutator pairing distinguishes classes
for k in range(3):
    print(k,"commutator pairing nondeg:",any((psi(k)(x,y)-psi(k)(y,x))%3 for x in E for y in E))
# right-action associativity: psi(e,qa)+psi(e+qa,qb) == psi(e,q(a+b))+psi(a,b) for all e,a,b  iff psi^q==psi
def rightassoc(p,M):
    return all((p(e,act(M,a))+p(add(e,act(M,a)),act(M,b))-p(e,act(M,add(a,b)))-p(a,b))%3==0 for e in E for a in E for b in E)
print("right assoc, det1:",all(rightassoc(psi(1),M) for M in mats if (M[0][0]*M[1][1]-M[0][1]*M[1][0])%3==1))
print("right assoc fails for every det -1:",all(not rightassoc(psi(1),M) for M in mats if (M[0][0]*M[1][1]-M[0][1]*M[1][0])%3==2))
