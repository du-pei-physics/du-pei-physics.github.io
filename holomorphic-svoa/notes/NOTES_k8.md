# The k = 8 case (c = 96, target index 3): where Codex stands and how to complement it

*Claude, 17 October 2026 (evening); prepared for publication 20 October (remarks about internal workflow removed or
reworded; the mathematics, the corrections and the data are unchanged). Paths such as `k8/…`, `results/k8/…`,
`anomaly/…` and `persistent/…` refer to the project workspace, and `docs/…`, `data/…` to Codex's files there; they are
not included in this folder unless an included counterpart is named.*

## 1. What Codex has done (read-only survey)

Target: GJF Theorem 5.6 with $k=8$, $m=\gcd(8,24)=8$: a self-dual ternary code $C$ of length
96 with $\operatorname{Index}(C)=24$ and a signed-permutation automorphism $g$ of order 8 with
$g^4=-1$; the $\mathbb Z_8$ orbifold then has index $24/8=3$. Codex's record
(`docs/progress.tex` §§ "Order-eight …", "Audited length-96 index-six orbifold", "Complete
Hermitian invariants", `docs/status_minimal_cases_2026_10_17.md`):

1. **Index-24 codes of length 96 exist** (several hundred, `docs/k8_index24_witness.json` and the
   order-4 campaign), and some carry an order-4 symmetry $h$ with $h^2=-1$: the
   **index-6 orbifold** (`docs/k8_index6_orbifold_witness.json`), the best certified value.
2. **Prescribed order-8 action.** Every $g$-invariant self-dual code has the form
   $C=U\oplus(K_-\cap U^\perp)$ with $U$ an $\mathbb F_9$-subspace of $K_+$ (the ambient space splits
   under $R=\mathbb F_3[g]/(g^4+1)\cong\mathbb F_9\times\mathbb F_9$). Codex samples $U$ uniformly on the
   Grassmannian, with neighbour/maximal-word-conditioned variants. **Every exact index seen is
   divisible by 48**: 14,761 length-96 labels in the continuous collector (minimum nonzero 48),
   plus all sampled lengths 24 and 48. At length 24 Codex **proved** $48\mid$ index for this class
   of actions (pure maximal words pair up in 16s; the mixed-word part by exhaustion of 9104 codes).
3. **Square roots of the order-4 symmetry.** For the index-6 code and 43 further index-24/order-4
   codes, complete low-weight layers (ternary weight 3/6; Hermitian weight 5/6 over $\mathbb F_9$)
   exclude any $g$ with $g^2=h$ for 12 of them; CP-SAT on the rest: UNKNOWN. Twelve codes are
   excluded for *every* signed order-8 action by the orbit argument "$g^4=-1$ acts freely on
   projective lines, so every weight layer has size divisible by 4" (`docs/k8_task2926_global_order8_obstruction.json`).
4. Codex's continuous search at the time (`data/continuous_minimal_cases_v1`): two of three lanes on $k=8$
   (maximal-word-conditioned order-8 neighbours at interface rank 6; broader rank-8 neighbourhoods).

Independent evidence from the Claude team: the transfer-ring search found that σ-equivariant rings
(an internal signed order-8 symmetry, $\sigma^4=-1$) always have $16\mid$ index (359 cells, all
$L\le12$; 70k exact 2-cell products; 118k 4-cell products), while generic cells with the same
delay space have $v_2=3$ (`k24_transfer/NOTE_kspan.md` §4.1).

**Reading.** Codex's main line — generate codes with the prescribed $g$ and hope for index 24 — has
produced a large, consistent body of negative evidence: in every family tried, $g^4=-1$ comes
with an extra factor of 2. If that is a theorem, the line cannot succeed, and the k=8 minimum (if
it is 3) needs a different kind of symmetry. The complementary work is therefore (i) to settle
the divisibility, and (ii) to identify what construction could evade it.

## 2. A reduction valid at every length (proved here)

Let $n=4N$ with $N$ even, $g$ as above, $R=\mathbb F_3[g]/(g^4+1)$, $X\subset R$ the 16 maximal
local words of one 4-block. Facts, all checkable in a few lines (`k8/local_facts.py`):

- $g$ flips the parity of a maximal word by $(-1)^N=+1$, so $\langle g\rangle$ acts freely on the
  maximal words of $C$ with orbits of size 8 of constant sign: **$8\mid\operatorname{Index}$**.
- $X$ consists of two $\langle g\rangle$-orbits, of $(1,1,1,1)$ ("type a") and $(1,1,2,1)$
  ("type b"); every element of $X$ is a unit of $R$. Call a maximal word *pure* if all its $N$
  blocks have the same type, *mixed* otherwise.
- $R^*\cong\mathbb F_9^*\times\mathbb F_9^*$ (order 64) acts on $C$ (multiplication by units
  preserves any $R$-submodule, and $C$ is one), freely on words with unit entries. Averaging the
  index over this action gives
  $$\operatorname{Index}(C)=\sum_{\text{free }R^*\text{-orbits }O}S(O),\qquad
    S(O)=\sum_{\lambda\in R^*}\prod_{i=1}^N\chi(\lambda c_i),$$
  with $\chi=\phi\cdot\mathbf 1_X$ on each block. $S$ is constant on cosets of $\langle g\rangle$,
  so $S=8\sum_{t\in R^*/\langle g\rangle\cong\mathbb Z_8}F_t$ with $F_t\in\{0,\pm1\}$, and $F_t\neq0$ iff
  $t\in\bigcap_i(\bar X-\bar c_i)$, where $\bar X=\{\bar a,\bar b\}\subset\mathbb Z_8$ with
  $\bar b-\bar a=6$. For a pure word this intersection has two elements, for a mixed word one.
- Hence **pure maximal words contribute multiples of 16 to the index, and each mixed
  $\langle g\rangle$-orbit contributes exactly $\pm8$:**
  $$\operatorname{Index}(C)\equiv 8\cdot\#\{\langle g\rangle\text{-orbits of mixed maximal words}\}\pmod{16}.$$

So "$16\mid\operatorname{Index}$ for all $g$-invariant self-dual codes" is equivalent to "the number
of mixed maximal words is divisible by 16", a statement about $C=U\oplus(K_-\cap U^\perp)$ that
does not involve signs at all. (Codex's length-24 proof is exactly the pure-word half of this plus
an exhaustive check of the mixed half.)

## 3. Experiment: is the mixed-orbit count always even?

`k8/order8_div16.py` samples $g$-invariant self-dual codes (random $\mathbb F_9$-subspaces $U$ of
random dimension) and computes exact indices and, for $n\le24$, the pure/mixed split.

| length $n$ ($N=n/4$) | codes sampled | indices | $16\mid$ index | mixed-orbit parity check |
|---|---|---|---|---|
| 16 (4) | 300 | all 0 | ✓ | ✓ (300/300) |
| 24 (6) | 300 | gcd 48, min 48 | ✓ | ✓ (300/300) |
| 32 (8) | 400 | all 0 | ✓ | – |
| 40 (10) | 300 | all 0 | ✓ | – |
| 48 (12) | 150 | gcd 48, min 48 | ✓ | – |
| 96 (24), Codex's collector | 14,761 | gcd 48, min 48 | ✓ | – |

(`results/k8/div16_n*.json`; uniform random $\mathbb F_9$-subspaces $U$ of random dimension.)
Two observations: at $n\not\equiv0\pmod{24}$ every sampled index is $0$, and at
$n\equiv0\pmod{24}$ the index is always a multiple of 48. Not a single exception in 1,450 codes
here plus Codex's 14,761. The pure-word part is a theorem (§2); the mixed-word parity is the
conjecture:

> **Conjecture A.** If a self-dual ternary code has a signed-permutation automorphism $g$ with
> $g^4=-1$, its index is divisible by 16. Equivalently, its mixed maximal words come in multiples
> of 16.

If Conjecture A holds, **GJF Theorem 5.6 cannot produce the minima at $k=8$, $16$ or $24$** through a
code automorphism: every such $g$ forces index $\ge48$, orbifold index $\ge6$. This is consistent
with every negative result Codex has accumulated, and it says the prescribed-order-8 generation
(two of Codex's three current lanes) cannot succeed.

## 4. Consequences and the complementary plan

### 4.1 What Conjecture A rules out

Codex's generation of codes with a prescribed signed $g$, $g^4=-1$ (both of its current $k=8$
lanes, the square-root CP-SAT searches, the Grassmannian sampler) and the σ-rings of the
transfer-ring search all live inside the class covered by Conjecture A. If A holds, none of them
can ever produce index 24, and the best they can give is index 48, i.e. orbifold index 6 — which
is exactly where everything has stalled. The same applies to $k=16$ (order 8) and $k=24$ (order
24 contains a $g$ with $g^4=-1$: $g=z^3$ for an element of order 24 with $z^{12}=-1$).

### 4.2 A way out: non-cyclic groups in which every commuting pair is cyclic

GJF's index computation needs only two things from the symmetry group $G$ of the code: that
every commuting pair $(a,b)$ generates a *cyclic* subgroup (so that $Z(a,b)$ is a modular image of
a twining trace $Z(1,c)$), and that every $c\neq1$ either has $-1$ as a power or (odd order) is
balanced. Then
$$\operatorname{Index}(V_{\Lambda(C)}/\!\!/G)=\operatorname{Index}(C)/|G|,$$
provided the $G$-orbifold is non-anomalous. Groups in which every abelian subgroup is cyclic are
the groups with periodic cohomology: $\mathbb Z_m$, the quaternion group $Q_8$, the dicyclic groups
$\mathrm{Dic}_n$ ($|{\rm Dic}_n|=4n$), $\mathrm{SL}(2,3)$ (order 24), the binary polyhedral groups.
They open exactly the rows closed by Conjecture A:

| target | group | needs | index of code | orbifold |
|---|---|---|---|---|
| $k=8$ | $Q_8$ ($i^2=j^2=k^2=-1$) | nothing beyond the group | 24 | $24/8=3$ |
| $k=8$ | $\mathrm{Dic}_4$ (order 16, contains $\mathbb Z_8$) | nothing | 48 (what Codex's $\mathbb Z_8$ codes have!) | $48/16=3$ |
| $k=16$ | $Q_8$ | nothing | 24 | 3 |
| $k=12$ | $\mathrm{Dic}_3$ (order 12) | the order-3 element balanced | 24 | 2 |
| $k=24$ | $\mathrm{SL}(2,3)$ or $\mathrm{Dic}_6$ (order 24) | order-3 elements balanced | 24 | **1** |

The monomial realisation with $-1\mapsto-I$ is forced: the 4-coordinate module
$\mathrm{Ind}_{\langle-1\rangle}^{Q_8}(\mathrm{sgn})$ with coordinates $(1,i,j,k)$ (left
multiplication), so a length-96 code is 24 copies of it. The three elements $i,j,k$ each act as
48 negative 2-cycles — precisely the order-4 action of Codex's index-6 orbifold — but now there
are three of them, pairwise anticommuting, and no element of order 8.

## 5. Quaternion-equivariant rings: index-24 codes of length 96 with a monomial $Q_8$ (found)

`k8/q8_cells.py` (included as `../code/q8_cells.py`) builds 12-coordinate transfer cells on which $Q_8$ acts internally (three blocks
$(1,i,j,k)$): the ambient cell space is $\mathrm{irr}\otimes\mathbb F_3^6$ with the standard form
$\omega\otimes\Omega$, the equivariant anti-orthogonal $S$ are $1\otimes s$ with
$\Omega(ms,m's)=-\Omega(m,m')$, the equivariant orthogonal maps are $1\otimes\mathrm{Sp}(6,3)$, and
the invariant non-degenerate delay spaces are $\mathrm{irr}\otimes Y$ with $Y$ a symplectic plane
($r=4$, 81 states). Rings of such cells are $Q_8$-invariant self-dual codes whose index is a trace
of a product of tables.

Survey (`results/k8/q8_survey.json`, 8 minutes): identical cells again have $v_2\ge4$ at every
$L$ (the Weyl/Burnside mechanism of Theorem 2 acts on them too), **but products of four different
equivariant cells with a common $U$ have $v_2=3$ in 37% of cases and hit $|{\rm index}|=24$
regularly: 42 hits in 128,000 products** (28 negacyclic, 14 cyclic). All 42 were verified on
their expanded $48\times96$ generators (`k8/verify_q8_hits.py`, `results/k8/q8_hits_verified.json`):
self-dual; invariant under $i,j,k$ acting on all eight blocks; $i^2=j^2=k^2=-1$, $ij=k$, $ji=-k$;
each of $i,j,k$ is 48 negative 2-cycles; no non-identity element fixes a nonzero codeword; index
$\pm24$ by the independent generator-based transfer computation. Certificate with the full
generator and the two generating symmetries: `results/k8/q8_certificate.json` (included as
`../certificates/original/q8_certificate.json`; normalised: `../certificates/k08_claude_q8.json`).

**Consequence.** If the lifted $Q_8$-action on $V_{\Lambda(C)}$ is non-anomalous, then
$$I\bigl(V_{\Lambda(C)}/\!\!/Q_8\bigr)=\tfrac{24}{8}=3\qquad(c=96,\ k=8),$$
the TMF minimum, with a single NS vacuum (every twisted sector has oscillator energy $\ge 9/2$).

## 6. What is not yet verified: the $Q_8$ anomaly

The restriction of the anomaly to each cyclic subgroup $\langle i\rangle,\langle j\rangle,\langle k\rangle\cong\mathbb Z_4$
vanishes by the same data as for the known index-6 orbifold (Frame shape $4^{48}2^{-48}$ = 24
copies of $4^22^{-2}$, Majorana pairing $288\in4\mathbb Z$). But $H^3(Q_8,U(1))\cong\mathbb Z_8$ and the
restriction to the three $\mathbb Z_4$'s has a kernel of order 2 (the class $4e$, $e$ the Euler class of
the quaternionic representation), so the cyclic checks do not decide the full anomaly, and all
genus-one data are blind to it (every commuting pair is cyclic). Two arguments suggest it vanishes:
the action is 24 copies of a single 4-dimensional module (GJF's "24 times" argument for the
anomaly, if it extends to the non-cyclic class), and $24\,p_1$ of the representation vanishes. A
proof needs either a general result on anomalies of lattice-VOA orbifolds by finite groups acting
through a lift of a subgroup of $O(\Lambda)$ (Gemünden–Keller style), or a direct check of the
associativity data of the $i$-, $j$-, $k$-twisted modules. **This is the single remaining step
between the codes found here and a $k=8$ theory of index 3.**

## 7. Recommendations

1. Treat Conjecture A as the working hypothesis; redirect Codex's two $k=8$ lanes from prescribed
   $\mathbb Z_8$ actions to $Q_8$-invariant codes (abundant: this note's family produces one per
   ~3,000 products) and to the anomaly question. A proof of A (parity of the mixed maximal words)
   would be a clean theorem; the pure-word half is proved above.
2. Settle the $Q_8$ anomaly. If it vanishes, $k=8$ is done at index 3 by the certificate code.
3. The same construction gives candidates for $k=16$ ($Q_8$, eight cells) and, with order-3
   elements that are balanced, for $k=24$ ($\mathrm{SL}(2,3)$ or $\mathrm{Dic}_6$ acting internally on
   twelve 12-cells): a *single* orbifold of a lattice SVOA with index 1, i.e. the periodicity
   element without symmetric squares.


## 8. Update (17 October, night): k = 16 joins; k = 24 routes

The transfer-ring search ported the $Q_8$-equivariant cells into its product engine
(`k24_transfer/NOTE_kspan.md` §5): 252 more verified $k=8$ codes, and **two verified length-192
codes with index $-24$ and a monomial $Q_8$ action** (one per twist; `persistent/q8_k16/hits.jsonl`).
I re-verified both independently (`results/k8/q8_k16_*_verified.json`: self-dual, $Q_8$-invariant,
$i,j$ of order 4 with no fixed codewords, index $-24$ by my generator DP). So **$k=16$ is in the same
position as $k=8$: index 3 modulo the $Q_8$ anomaly.**

For $k=24$ that search found that neither $\mathrm{SL}(2,3)$ nor $\mathrm{Dic}_6$ admits an equivariant
cell on the natural modules (no anti-orthogonal $S$, resp. no non-degenerate invariant $U$), while
$Q_8\times\mathbb Z_3$ with the $\mathbb Z_3$ realised by the cell shift does ($p=4$ rings of twelve
cells: $\operatorname{Tr}JM=0$, $|\operatorname{Tr}JM^3|=24$; unobstructed, first values $\sim10^7$;
worker `persistent/q8_k24`). A hit there would be a single-orbifold periodicity element, again modulo
the anomaly of the order-24 group.

## 9. Anomaly update (17 October, late): the $Q_8$ action is non-anomalous (conditionally)

Full trail in NOTES_q8_anomaly.md (this folder). Summary:

* **Lift (S1, proved/cited).** From the certificate, $i^2=j^2=k^2=-I$, $ij=k$, $ji=-k$, the group
  generated has order exactly 8, preserves the 24 blocks and the code; by GJF Theorem 3.1 the
  extension $\hat\Lambda(C)[2].\mathrm{Aut}(C)$ splits, so $Q_8$ lifts to an honest $Q_8\subset
  \mathrm{Aut}_{N=1}(V)$ (unique up to conjugacy since $H^2(Q_8,\Lambda)=0$ by Shapiro), and $-I$ is
  the bosonic lift of $\lambda\mapsto-\lambda$ (GJF Thm 5.6, Ex. 5.7(2)), not $(-1)^F$. Majorana layer
  vanishes for all elements ($\langle\chi,(1-g)\chi\rangle\in\{288,576\}\subset4\mathbb Z$).
* **Bookkeeping (S2, cited).** $H^3(Q_8,U(1))=\mathbb Z_8$, joint kernel of restrictions to the cyclic
  subgroups $=\{0,4e\}$; Yu (2105.02231, Prop. 1.1) proves genus-one data cannot see it. The fermionic
  anomaly group is $\mathrm{Hom}(\Omega_3^{\rm Spin}(BQ_8),U(1))\cong\mathbb Z_8\oplus\mathbb Z_4^2$
  (Bayen–Bruner 1996; Botvinnik–Gilkey 1995; Davighi–Gripaios–Lohitsiri 2207.10700 Table 6):
  **exponent 8**.
* **The argument (S3).** The "24 identical blocks" idea cannot be run on $W=B^{\otimes24}$ (not
  holomorphic). Its correct form: $\Lambda(C)\oplus(-\mathbb Z^{96})\cong\mathbb Z^{96}\oplus(-\mathbb Z^{96})$
  as $Q_8$-lattices (same genus — local isometries at every $p$ via the Lipschitz-order hermitian
  structure, Morita at $p=3$ — plus strong approximation for $\mathrm{Sp}(24,24)$), so the Narain
  CFT $V_{\Lambda(C)}\otimes\overline{V_{\mathbb Z^{96}}}$ deforms $Q_8$-symmetrically to
  $V_{\mathbb Z^{96}}\otimes\overline{V_{\mathbb Z^{96}}}$; deformation invariance of the anomaly gives
  $\alpha(V_{\Lambda(C)})=\alpha(\mathrm{Fer}(192))=24\,\alpha(\mathrm{Fer}(8))=0$ by exponent 8.
  Sanity: restricted to the cyclic subgroups this is GJF's own "24 times" statement.
* **Verdict.** $\mathrm{Index}(V_{\Lambda(C)}/\!/Q_8)=\pm3$ at $c=96$, and $\pm3$ at $c=192$ for the
  length-192 $Q_8$ codes, is **established modulo one physics-standard assumption**: that the 't Hooft
  anomaly is constant along a connected symmetric family of CFTs (not a VOA theorem). All other
  inputs are exact computations or cited theorems. A fully VOA-level proof needs either "anomaly of
  a lifted $G\subset O(\Lambda)$ depends only on the genus of the $G$-lattice" or a direct computation of
  the invisible class (suggested exact check: $\mathrm{Fer}(32)/\!/\langle i\rangle$, whose residual
  $D_8$-or-$Q_8$ symmetry should carry exactly $4e$). Heuristically (η-invariants) the single block has
  anomaly of order 8, so blocks in multiples of 8 — hence $c\in32\mathbb Z$ — is exactly the condition.

### 9.1 Update (18 October): VOA-level proof, no physics assumption

NOTES_q8_anomaly.md (this folder) §S6 replaces the Narain-deformation step by a chiral argument:
(i) $\mathbb Z^{96}=\Lambda(\text{tetracode}^{\oplus24})$ for the $Q_8$-stable 3-frame $L^{24}\cdot(1+i+j)$,
with the *same* monomial $Q_8$ (exact check, `../code/neighbour_chain.py`); (ii) the certificate code
is connected to $\text{tetracode}^{\oplus24}$ by 23 $Q_8$-equivariant neighbour steps
$C_{i+1}=(C_i\cap U_b^\perp)+U_b$, $U_b$ the block tetracodes (exact); (iii) each step is the orbifold of
$V_{\Lambda(C_i)}$ by a $Q_8$-stable $\mathbb Z_3^2\subset C_i^*\subset\mathrm{Aut}_{N=1}$, and a gauging
lemma (anomaly of $\mathbb Z_3^2\rtimes Q_8$ is inflated from $Q_8$ by coprimality; gauging the normal
3-group preserves the $Q_8$ anomaly) gives $\alpha(V_{\Lambda(C_i)})=\alpha(V_{\Lambda(C_{i+1})})$;
(iv) $\alpha(\mathrm{Fer}(192))=24\,\alpha(\mathrm{Fer}(8))=0$ by exponent 8. **Index 3 at $c=96$ and
$c=192$ now rests only on cited theorems** (GJF Thms 3.1/4.2 and the standard twisted-module/algebra
dictionary, additivity under $\otimes$, $\tilde\Omega_3^{\rm Spin}(BQ_8)=\mathbb Z_8\oplus\mathbb Z_4^2$);
the gauging lemma is written as a sketch (S6.3(d)). By-product: all $Q_8$-codes with this block
structure at a given length have the same $Q_8$ anomaly.

### 9.2 Independent audit of §S6 (18 October, morning)

**Re-verified exactly (separate code using `f3lib`, independent of `../code/neighbour_chain.py`):** the frame
$F=L^{24}(1+i+j)$ consists of 96 orthogonal norm-3 vectors mapped to $\pm$ each other by the
certificate's $i,j$ exactly as the standard basis is; the frame code $C'$ of $\mathbb Z^{96}$ is self-dual
of dimension 48, invariant under the same signed permutations, and is the direct sum of 24 tetracodes
(block code of dimension 2, all eight nonzero words of weight 3); $\dim(C\cap C')=2$; the chain
$C_{i+1}=(C_i\cap U_b^\perp)+U_b$ over blocks $b=0,\dots,22$ gives 23 self-dual $Q_8$-invariant codes
with $\dim(C_i\cap C')=4,6,\dots,48$ and ends at $C'$; at every step $U_b$ is totally isotropic,
$C_i\cap U_b=0$, and $C_i\cap C_{i+1}=C_i\cap U_b^\perp$ has codimension 2 in both (so $V\cap V'=V^E$ with
$E\cong\mathbb Z_3^2$); and $i,j$ act on each $U_b$ with determinant $1$ (so $Q_8$ acts trivially on
$H^2(E;U(1))\cong\mathbb Z_3$, as used in S6.3(d)). Also checked by hand: the three layers of
$\widetilde{SH}^3(BQ_8)$ have orders $|H^1(Q_8;\mathbb Z_2)|\cdot|H^2(Q_8;\mathbb Z_2)|\cdot|H^3(Q_8;U(1))|
=4\cdot4\cdot8=128=|\mathbb Z_8\oplus\mathbb Z_4^2|$, so no differential kills $H^3(Q_8;U(1))$, as
S6.3(c) needs; the group $\Omega_3^{\rm Spin}(BQ_8)\cong\mathbb Z_8\oplus\mathbb Z_4^2$ is confirmed by
Davighi–Gripaios–Lohitsiri (arXiv:2207.10700, Table 6) and by a 2026 preprint listing its
generators (arXiv:2607.15349); exponent 8.

**Assessment of the mathematical steps.** (a), (c), (e) are complete proofs (the uniqueness of the lift
via $H^2(Q_8,\Lambda)=0$ by Shapiro is correct: $\Lambda\otimes\mathbb Z_2\cong 24\,L_2$ with $L$ induced
from the centre with sign action, $\hat H^0(\mathbb Z_2,\mathbb Z^{\rm sgn})=0$). (b) is the standard
identification of holomorphic extensions of $V^E$ with Lagrangian algebras of $Z(\mathrm{Vec}_E^{\alpha_E})$;
I checked the count: $N^*/N\cong\mathbb Z_3^4$ carries the split quadratic form, whose 8 Lagrangians
match the 8 pairs $(K\le E,\psi)$, and exactly 3 of them (the three $\psi\in\mathbb Z_3$) meet
$\Lambda(C_i)/N$ trivially, so $V'$ is indeed one of the three $V/\!/_\psi E$. (d), the gauging lemma, is
the standard statement that condensing a non-anomalous normal subgroup whose cocycle is inflated from
the quotient yields the quotient's crossed extension with the same cocycle (bosonic case: Tachikawa,
arXiv:1712.09542; pointed $G$-crossed categories make it an explicit cocycle computation); it is
sketched, not written out, and is the one place where a referee would ask for the full argument,
including the fermionic layers and the independence of the complement.

**Verdict.** I agree with the conclusion of §S6: the $Q_8$ action on the length-96 (and length-192)
code theories is non-anomalous, and the only remaining inputs are cited theorems (GJF 3.1/4.2 and the
twisted-module/algebra dictionary, additivity under $\otimes$, $\Omega_3^{\rm Spin}(BQ_8)$) plus the
standard gauging lemma. The Narain-deformation assumption is gone. Index 3 at $c=96$ and $c=192$,
and the tensor-product rows $k=17,19,20,22,23$, are therefore established at the same level of rigour
as the rest of the programme. What would make it airtight: a written proof of (d) in the pointed
super-crossed setting (two pages), and, as an independent numerical cross-check, the stepwise
computation $\mathrm{Fer}(32)/\!/\langle i\rangle$ suggested in S4.1 (one would like to see the
invisible class actually appear for 4 blocks and vanish for 8).

### 9.3 Gauging lemma written out (18 October)

The gauging lemma S6.3(d) now has a full proof, first written as a working draft (its section labels G0–G13, GP and
GQ are used in §§9.3–9.8); the final version is gauging_lemma.pdf in this folder (source gauging_lemma.tex), and the
correspondence of its numbering is given at the end of §9.8. The proof works in the pointed super-crossed setting:

- $F(V,H)\simeq\pi^*F(V,Q)$, an inflation;
- a $Q$-invariant cocycle exists, and for $\mathbb Z_3^2$ with $Q_8\subset SL(2,3)$ it is $\zeta^{k\det}$ on the nose (exact check `../code/psi_cocycle_check.py`);
- induction $X\mapsto A\otimes X$ is a fully faithful monoidal functor into $A$-bimodules that preserves the Majorana, Gu–Wen and bosonic layers;
- the full symmetry of $V/\!/_\psi E$ is $K=\hat E.Q$, and Schur–Zassenhaus fixes both the complement and the identification of the geometric $Q_8$.

It also corrects the sketch, which said $A\otimes 1_q$ is the only simple bimodule; there are
$|E|$ per coset. One cited input remains: the "twisted modules of the gauged theory $=$
$A$-bimodules" dictionary (D2). It is in the literature in bosonic form. We use it for
super-fusion categories without Majorana objects, for which we know no explicit reference.

*Update (18 October, later):* D2 is now proved (draft §G8; final version: gauging_lemma.pdf, Theorem 5.3). The proof
defines $\Psi(N)=V\boxtimes_{V^E}N$ from $F(V',K)$ to $\mathrm{Bimod}_A F(V,H)$, with
$A=V\boxtimes_{V^E}V'\cong\bigoplus_e M_e$. $\Psi$ is monoidal, using the standard étale-algebra
lemmas. It is fully faithful by a centrality computation: a vector in the charge-$\chi$ sector
commutes with the flux $N_e$ up to the monodromy $\chi(e)$. It is essentially surjective because
the simple bimodules are exactly the $|H|$ objects $B_X^\chi$. The only VOA input is the
"twisted modules $=$ modules of the algebra object in $\mathrm{Rep}(V^H)$" framework on which GJF
Thm 4.2 itself rests. Lemma S6.3(d) therefore adds no input beyond GJF Thms 3.1/4.2, apart from
routine ones such as quantum Galois theory, quantum dimensions and the étale-algebra lemmas.

### 9.4 Independent audit of the gauging lemma (18 October)

I read the whole draft (G0–G8) and checked each step; the exhaustive cocycle check
`../code/psi_cocycle_check.py` reproduces (cocycle identity, SL(2,3)-invariance on the nose, $\det=-1$ sends
$k\mapsto-k$, non-degenerate pairings, and the right-associativity identity G4.2(b) holds for every
determinant-1 matrix and fails for every determinant $-1$ matrix).

**Steps I verified line by line and consider correct.**
* G1.3/G1.4: pull-back of graded categories is natural on the class; with $\tilde\alpha=\pi^*\alpha$ the
  category $F(V,H)$ is an inflation of $F(V,Q)$. Uses only the classification statement attached to
  GJF Thm 4.2 (same class $\Rightarrow$ graded equivalent); G7.3(3) gives the skeletal alternative for $m=0$.
* G2–G3: twisted group algebras on $\bigoplus 1_e$; $|E|$ odd kills the parity homomorphism; the
  averaging trick $\bar\psi=(\prod_q\psi^q)^n$ with $n|Q|\equiv1\ (|E|)$ produces a representative that is
  $Q$-invariant *as a function* with no divisibility of $U(1)$ needed — correct and neat; for
  $\mathbb Z_3^2$ with $Q_8\subset SL(2,3)$ the bilinear $\zeta^{k\det}$ is invariant on the nose.
* G4: I re-derived the right-associativity computation: with $\psi(e,{}^q a)\psi(e\,{}^q a,{}^q b)=
  \psi({}^q a,{}^q b)\psi(e,{}^q(ab))$ from the cocycle identity, the bimodule axiom holds iff
  $\psi^{q}=\psi$ as functions — exactly where the invariant representative is used. Monoidality of Ind
  with no extra scalars (free-module isomorphism; right action of $B_Y$ transported through the
  normalised cocycle) and the associator being $a^G$ are right. The correction G4.6 (there are $|E|$
  simple bimodules per coset, $B_X^\chi$; the sketch's "unique simple" was wrong) is honest and does not
  affect the conclusion, since only the monoidal subcategory $\{B_{X_q}\}$ is used.
* G5.2: $K=\{g\in\mathrm{Aut}(V'):g(V^E)=V^E,\ g|_{V^E}\in H|_{V^E}\}$ is an extension of $Q$ by $\hat E$
  by quantum Galois theory applied twice; Schur–Zassenhaus (abelian kernel, coprime orders) splits it
  and makes complements conjugate. G5.3 then identifies $F(V',\tilde Q)$ with $\mathrm{Ind}(G)$ through
  $\Psi$ and coset compatibility; step 6 (conjugation by $c\in\hat E$) handles the other complements.
* G6: independence of the complement in $H$, of the representative $\psi$ (via $H^1(Q;\hat E)=0$) and
  of $\Phi$. G6.4 reproves S6.3(c).
* G8 (D2): the decisive part. (a) FPdims from Dong–Jiao–Xu; (b) the monodromy lemma
  $c_{N_e,V_\chi}c_{V_\chi,N_e}=\chi(e)^{\pm1}$ from the definition of an $e$-twisted module and Huang's
  braiding, with the image of the restricted vertex operator identified with the invertible
  $V_\chi\boxtimes N_e$; (c) $L\otimes_BN_e\cong M_e$ by adjunction plus equal FPdim ($|H|$ on both
  sides); (d) $A=L\otimes_BL'\cong\bigoplus_eM_e\cong A_{\psi'}$. G8.4 builds $\Psi(N)=L\otimes_BN$ as a
  composite of three monoidal functors (modules in stages, central algebras, transport of bimodules).
  G8.5, Lemma C: I re-did the hexagon/naturality computation; a vector in the charge-$\chi$ summand
  commutes with the flux $N_e$ up to $\chi(e)^{\mp1}$, so only $\chi=1$ vectors are central, giving
  $\mathrm{Hom}_{A\text{-}A}(A,\Psi N)=\mathrm{Hom}_{L'}(L',N)$ and full faithfulness by rigidity. G8.6,
  Lemma B: left modules are free; a simple bimodule is $(A\otimes X_q)\otimes W$ with the right action
  differing from the canonical one by an honest even representation of $E$ on $W$ (the multiplicativity
  uses that the canonical right action commutes with $\mathrm{End}(W)$), hence $P\cong\Pi^pB_{X_q}^\chi$;
  the count $|Q|\,|E|=|H|$ matches the $|K|$ simples of $F(V',K)$, so $\Psi$ is an equivalence.

**Inputs that remain cited (all standard; theorem numbers not re-checked by me either):** (T) twisted
modules $=$ not-necessarily-local modules over the algebra object in $\mathrm{Rep}(W^G)$
(Kirillov, McRae; super case as used by GJF for Thm 4.2); Carnahan–Miyamoto regularity for solvable $G$;
quantum Galois theory (Dong–Mason, Hanaki–Miyamoto–Tambara); quantum dimensions (Dong–Jiao–Xu);
HKL/CKM extension $\leftrightarrow$ algebra; Kirillov–Ostrik, Schauenburg and DMNO on étale algebras;
Huang's braiding. The Majorana-free assumption $m(V,Q_8)=0$ of G8 holds in the application
(GJF Ex. 4.4, S1).

**Verdict.** Lemma S6.3(d) is proved, not sketched, and D2 — the one dictionary the first version
still cited in bosonic form — is now derived from the same input (T) on which GJF Theorem 4.2 rests.
With this, the chain of §S6 is complete: the vanishing of the $Q_8$ anomaly on the length-96 and
length-192 code theories, hence **index 3 at $c=96$ and $c=192$** and the tensor-product rows
$k=17,19,20,22,23$, rests on exact computations, elementary group cohomology written out in these
notes, and cited theorems (GJF 3.1/4.2 with (T), the VOA facts above, additivity of anomalies under
$\otimes$, $\Omega_3^{\rm Spin}(BQ_8)\cong\mathbb Z_8\oplus\mathbb Z_4^2$). What a write-up for
publication would still want: the bibliography re-checked against the papers, and the two-line remark
of G8.7(2) on $q$-type objects if one wants the lemma in full generality.

### 9.5 Codex's review of the gauging lemma (18 October) and my assessment

Codex read the draft independently (Codex's review of 18 October, not included). It
confirms G2–G4 (invariant discrete torsion, the bimodules, preservation of the Gu–Wen and bosonic
layers, the corrected bimodule count) and raises four points. My assessment:

1. **Narrow the statement (G0).** Correct: the proof uses $E$ abelian of odd order and $m(V,Q)=0$, not
   just coprimality. Pure hygiene — the application satisfies both.
2. **"The unit of $G$ placed in degree $s(q)$" (G3.1) does not exist.** Correct: in $\pi^*G$ the unit can
   only be placed in degrees $h\in E$. The repair is conjugation by an invertible $X_q$ of degree $s(q)$;
   because the unit is strict and the associators are inflated, this conjugation is *strictly* monoidal
   on the objects $1_e$ (zigzag + triangle, constant tensorator forced to be 1 by the unit axiom), so it
   sends $A_\psi$ to $A_{\psi^q}$ exactly and the $Q$-action on classes is conjugation in $H$, as used.
   G4 is unaffected (it only places $1_{e'}$ with $e'\in E$).
3. **"$N_e$ is even" (G8.5) is not a reason.** Correct wording issue; the valid reason is that every
   structure morphism in the computation is even and the single test vector is never interchanged past
   another odd morphism, so no Koszul sign can arise.
4. **The substantive one: McRae's theorems need the group to contain $(-1)^F$.** Our $H$, $K$ preserve the
   odd $N=1$ current, so $(-1)^F\notin H,K$ and (T) cannot be cited for $G=H$ as written. The fix is a
   parity adjunction: $\tilde H=H\times\langle(-1)^F\rangle$, common bosonic fixed algebra
   $U_0=(V^H)_0$, $\tilde{\mathcal C}=\mathrm{Rep}(U_0)$ modular, McRae for $(V,\tilde H)$ and $(V',\tilde K)$
   (CKM for the superalgebra objects), then G8.3–G8.6 inside $\tilde{\mathcal C}$ and the dictionary
   between GJF's super-fusion $F(V,H)$ and the $H$-graded part of the bosonic $\tilde H$-crossed category
   (bosonic Hom = full super Hom with the grading forgotten; $\Pi$ from the $(-1)^F$ action). What was
   needed: this appendix, together with fixes 1–3 and the DMNO §3.6 citation. Final status: the parity adjunction is
   §2.1 of gauging_lemma.pdf, and D2 is its Theorem 5.3, proved from CKM and McRae's Thm 4.15 (status table, §9
   there).

Codex also reports that its own workers reached a balanced floor of 48 at both $k=18$ and $k=21$ after
$2\cdot10^7$ products, better than our 240, with its random anti-orthogonal cells of widths 8/10/12 —
worth adopting (§10.3 of the k-span working notes; workspace, not included). None of the review's points suggests the $Q_8$
anomaly could be nonzero; they concern the written chain of citations.

*Revision (18 Oct, after Codex's review of 18 October, not included):* the draft's §G9 narrows Lemma G (E odd abelian, m(V,Q)=0), repairs G3.1 (conjugation by an invertible $X_q$), fixes the parity-sign and DMNO wording, and adds appendix GP deriving D2 from McRae's Main Theorem 1 after adjoining $(-1)^F$ (H̃ = H×⟨σ⟩ over the bosonic $U_0$); the one open point is the definitional match of GJF's $F(V,H)$ with McRae's H-graded part.
Note: arXiv:2207.10700 in §9.2 above is Davighi–Gripaios–Lohitsiri, not Dierigl–Oehlmann–Schimannek.

### 9.6 Re-audit after the revision (18 October, 10:15)

I read the revised draft (G0 narrowed; G3.1 repaired; the Koszul-sign reason
corrected; DMNO cited as the §3.6 remark; new Appendix GP; G9 summary). All five review points are
addressed as proposed, and the revised draft itself now states that its earlier "D2 adds no input beyond
GJF Thm 4.2" was too strong. Checks: the repaired G3.1 is right — with $G$ strict, the conjugation
$c_q=\bar X\otimes-\otimes\bar X^*$ has tensorator $\mathrm{id}\otimes\mathrm{ev}\otimes\mathrm{id}$, and on
the objects $1_e$ both comparison maps become maps $X\otimes X^*\otimes X\otimes X^*\to1$ equal to $1$ after
precomposition with $\mathrm{coev}\otimes\mathrm{coev}$ (zigzag), so $c_q(A_\psi)=A_{\psi^q}$ with no
constant. In Appendix GP, $\sigma\notin H$ and $\sigma'\notin K$ are proved correctly (odd $N=1$ current;
Galois for $(V,\tilde H)$); $U_0=(V^H)_0$ is a regular VOA ($V_0$ a lattice VOA, Carnahan–Miyamoto for
solvable $H$), so McRae's Main Theorem 1 applies to $(V,\tilde H)$, $(V',\tilde K)$, $(V^E,\tilde Q)$ with
its hypothesis "$G\ni$ parity" satisfied; the super ambient category $\tilde{\mathcal C}\boxtimes\mathrm{SVec}$
keeps parity and $\Pi$; the FPdim ratios used in G8.3(c) are unchanged (everything doubles); the
restriction to the $H$- and $K$-graded parts (GP.5) gives exactly the D2 used in G5.3 and, as a
by-product, $m(V',K)=0$.

**What remains open is definitional and bibliographic, not mathematical:** (i) a page-level statement
that GJF's twisted-module category $F(V,H)$, which defines $\alpha$, is the $H$-graded part of McRae's
$\mathrm{Rep}_{S\tilde{\mathcal C}}(V)$ (both are $h$-twisted NS-type modules with intertwining-operator
fusion); (ii) CKM theorem numbers and the HLZ applicability of $\mathrm{Rep}(U_0)$. I consider the
gauging lemma proved in the setting used, and the $Q_8$ anomaly result closed modulo (i)–(ii) and the
inputs already listed in §9.4.

### 9.7 Codex's second review (18 October, 11:15)

Codex's second review (not included): points A–E of the first review are accepted ("no new
counterexample to the categorical argument"). One substantive point remains — the identification of GJF's
$F(V,H)$ with the $H$-graded part of McRae's category must be made at the level of graded monoidal
supercategories (the anomaly is the associator), and McRae's Remark 4.16 leaves the comparison with Huang's
twisted-intertwiner tensor structure open, so "a matter of definitions" is too quick. Codex proposes to
identify GJF's own construction (proof of Thm 4.2) directly with the algebra-object construction over the
common fixed algebra by a universal-property argument, and suggests the status sentence "The categorical
argument proves the lemma under the stated hypotheses and the graded monoidal identification described in
GP.6(1); completing that identification and the precise extension-theory citations closes the application
to the GJF anomaly." I agree with both. What was needed to close GP.6(1): determine how GJF
define the fusion of twisted modules; if via the algebra object, the identification holds by definition and must be
spelled out; if via intertwining operators, a universal-property comparison is needed, citing McRae Thm 4.17 /
Rem. 4.16 for what is proved. Editorial items: exact CKM theorem numbers, an explicit enriched supercategory, lemma
and application separated, inputs kept separate. Final status: gauging_lemma.pdf, Theorem 2.1 (proved in §7)
identifies GJF's twisted-module category, through Kirillov's algebra-object construction, with the graded part of
McRae's category as graded monoidal supercategories; the one interpretive step is Remark 7.2.

*Second revision (18 Oct):* the draft's App. GQ proves GJF's twisted-module category (Kirillov's algebra-object construction, which GJF's Thm 4.2 proof follows) is isomorphic to the graded part of McRae's as graded monoidal supercategories (equal associators), with exact CKM/McRae citations (GQ.5); §G10 records the one interpretive point (McRae Rem. 4.16 matters only if α is defined via Huang's analytic twisted fusion). Final version: gauging_lemma.pdf, Theorem 2.1, §7 and Remarks 7.1–7.2.
Separate inputs unchanged: Ω₃^Spin(BQ8) exponent, lift identification S6.3(e), existence of the lattice extensions.
*Citation audit (18 Oct, later):* the draft's §G11 checks HKL/CKL/CM/DJX numbers against full texts, corrects the Huang reference (Commun. Contemp. Math. 10 (2008) 871–911), and replaces the quantum-dimension input in the application by lattice simple-current fusion; open general-case items: (R) for $V'_0$ and the GJF reading (GQ.1).
*Consolidated proof (18 Oct):* gauging_lemma.tex/.pdf (included in this folder) is the clean write-up of the gauging lemma; (R) is now reduced to $V_0$ regular (Lin Thms 3.10/3.11 + ABD), the one interpretive step being GJF's tensor structure = Kirillov's (Rem. 7.2 there).

### 9.8 Audit of the second revision (18 October, 12:45)

Appendix GQ closes GP.6(1) properly: GJF's proof of Thm 4.2 defers the tensor structure to Kirillov (CMP 2002),
where twisted modules are *defined* as objects of $\mathrm{Rep}\,A$ graded by monodromy (Def. 4.1), the product is
$\otimes_A$ and the anomaly is the associator (eq. 5.3); Theorem GQ gives an explicit isomorphism $\Theta$ of
graded monoidal supercategories between Kirillov's construction over $\mathrm{Rep}^0(V^G)$ and the $G$-graded part of
McRae's $\mathrm{Rep}_{S\tilde{\mathcal C}}(V)$, with tensorator $\tau$ from the universal property of the relative
tensor product (CKM Thm 2.52) and equal associators and units. McRae's Remark 4.16 concerns Huang's *analytic*
twisted fusion, which neither GJF nor the lemma use. I checked the structure of the proof (steps 1–5) and the
citation table GQ.5; the only interpretive point is the reading that GJF's "regularity provides the fusion"
means Kirillov's $\otimes_V$ — the natural reading since Kir02 is the only construction in their cited source.
G11 further removes DJX from the application (lattice coset modules are simple currents) and closes hypothesis
(R) for $V'_0$ in general via Lin 2107.09400 Thms 3.10/3.11. A self-contained LaTeX version (gauging_lemma.tex,
included) accompanies the draft. **Verdict: the gauging lemma and its application to the $Q_8$ anomaly are closed,
modulo the stated reading of GJF's definition and the standard cited theorems; no mathematical gap remains.**
*Fourth pass (18 Oct):* re-reading found and fixed five points in the gauging-lemma proof (draft §G12): the
coprime splitting behind S6.3(c) became a $p$-local isomorphism argument (final version: gauging_lemma.pdf,
Proposition 2.3); the base category of D2; the definition of $V/\!/_\psi E$; the intrinsic form of hypothesis (3);
and a missing non-vanishing citation. The conclusion is unchanged.
*Fifth pass (18 Oct):* the draft's §G13 applies Codex's latest review: the inclusion
$H^3\hookrightarrow r\widetilde{SH}{}^3$ is stated by GJF, so the earlier argument for S6.3(c) is valid too (Remark 2.4 of
the final version); canonical-equivalence wording; cokernel preservation in GQ.3; separate lattice, cocycle and
invariance steps in the application. gauging_lemma.tex/.pdf is the authoritative version.

*Correspondence with the final version (gauging_lemma.pdf).* The statement is Theorem 1.1; the parity adjunction
(draft App. GP) is §2.1; the identification with McRae's category (App. GQ) is Theorem 2.1, proved in §7, with
Remarks 7.1–7.2 (McRae's Remark 4.16, the interpretive step); inflation and coprime splitting are Propositions 2.2–2.3;
invariant cocycles and the case $E=\mathbb F_3^2$, $Q_8\subset SL(2,3)$ are Lemmas 3.2–3.3; the bimodule embedding and
classification are Proposition 4.1 and Lemma 4.2; the gauged theory and the bimodule dictionary D2 are
Proposition 5.1 and Theorem 5.3; independence of choices is Proposition 6.1; the application is §8 and the status
table §9.
