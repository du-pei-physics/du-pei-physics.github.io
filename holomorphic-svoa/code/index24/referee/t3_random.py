"""8 | Index for random self-dual codes from [I | A0 Q], A0 A0^T = -I, Q a random
product of reflections (independent of the neighbour walk)."""
import sys
import numpy as np
from collections import Counter
from f3 import *

n = int(sys.argv[1])
num = int(sys.argv[2])
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 3)
vals = Counter()
bad8 = []
for t in range(num):
    G = random_selfdual_AAT(n, rng)
    assert is_self_dual(G, n)
    G = random_monomial(G, rng)
    iC = index_fast(G)
    vals[abs(iC)] += 1
    if iC % 8:
        bad8.append(iC)
print(f"n={n}: {num} codes; |Index| distribution {sorted(vals.items())}")
print("   not divisible by 8:", bad8[:10], " not divisible by 24:",
      sum(c for v, c in vals.items() if v % 24))
