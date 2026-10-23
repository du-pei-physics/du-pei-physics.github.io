> **Annotations for the online note (19 October 2026; prepared for publication 20 October).** This file is the working
> note of 16–17 October. For publication, remarks about internal workflow (task routing, compute allowances,
> internal attributions of checks) were removed or reworded; the mathematics, the corrections and the data
> are unchanged. File names such as `core.py` and paths such as `results/…` refer to the project workspace (folder
> `k24_transfer`); apart from the few archived with the benchmark sources in `../code/bench_v3/claude/k24_transfer/`,
> they are not included. Four blockquotes marked "Annotation" were added (here, after Theorem 1 and twice in the
> proof of Theorem 2). They make three points precise, following Codex's reviews of 19 October: (i) *normalisation:* $F_U$ as defined
> satisfies $F_U^2=3^rJ$, $F_U^4=3^{2r}\cdot1$, so the matrices $F_U$ and $DF_UD^{-1}$ generate an infinite group; the
> $Q_8$ of the proofs is the group of symplectic maps of $W=U\oplus U$ that conjugation by them induces on the Weyl
> operators, which scalars do not affect. (ii) *Scope of Theorem 2:* it concerns the anti-periodic closure
> $\operatorname{Tr}(JM^L)$ only, and needs $M$ to commute with $F_U$ and with $DF_UD^{-1}$, which holds for every
> even width and for products with a common boundary space and form (if $DMD^{-1}=cM$ then
> $(DF_UD^{-1})M=c^{-1}DF_UMD^{-1}=c^{-1}DMF_UD^{-1}=M(DF_UD^{-1})$). (iii) *The shift in step 3 of Theorem 2:* the
> displayed map has order $L$, not $2L$, when $L$ is odd; the anti-periodic shift $T(w)=(w_2,\dots,w_L,-w_1)$ has
> $T^L=-\mathrm{id}$, gives the same group, and comes with an explicit stabiliser lemma. See details.html B.6 and C.3.

# Direct periodicity element from transfer cells?

*Claude, 16–17 October 2026. Follow-up to §2 and §5.3 of the project's earlier notes (workspace, not
included). Everything here is self-contained: numpy, plus `f3lib.py` (included as `../code/f3lib.py`) and
`transfer_cells.py` (workspace) for cross-checks.*

## Answer first

1. **A direct c=288, index-1 orbifold from twelve identical 12-coordinate cells is
   impossible. This is proved.** For *every* cell, $16\mid\operatorname{Tr}(JT^{12})$
   (Theorem 2). This holds for every $r$, every $Q,S$, and also for the new cells whose delay
   space has non-square discriminant. So the 12-cell code never has index $\pm24$, whatever
   symmetry $g$ one uses. The same theorem kills 24 identical cells of size 6, and every
   ring of identical cells with $4\mid L$. The earlier observation "$I_4=0$ never occurred"
   was a side symptom; $I_{12}$ itself is obstructed.
2. **A c=144 family was found, but it has orbifold index 4, not 2** (correction of 17 October,
   after an independent check). A ring of six 12-coordinate cells alternating two cells
   $S_1,S_2$ with a common delay space $U$ (a *period-2 ring*) gives self-dual codes of length 144
   with index $\pm24$, a symmetry $g$ = shift by one period with $g^3=-1$, $g^2$ a fixed-point-free
   product of 3-cycles, and balance 0. **23 such codes were verified** (§5). But $g$ has order 6, so
   GJF Thm 5.6 gives a $\mathbb Z_6$ orbifold of index $24/6=4$. An index-2 theory at $c=144$ needs
   a symmetry of order 12 with $g^6=-1$ and $g^4$ of shape $3^{48}$; the one-cell shift would do,
   but it does not preserve a period-2 code. So these are $c=144$ index-4 theories (like the
   $c=72$ seeds, one level up), and they do **not** shorten the periodicity chain: $\mathrm{Sym}^2$
   sends $\pm4$ to 14 or 2, never to 1. The genuine index-2 target with the shift symmetry is
   stated in §13.

3. **No transfer ring whose GJF symmetry is a power of the cell shift can be a direct c=288
   index-1 orbifold (Proposition, §13).** Index 1 needs $|g|=24$ with $g^{12}=-1$. For a ring of
   $N$ cells of size $q$ ($qN=144$) with period $p$ and $g=z^{p}$ this forces exactly $L=N/p=12$
   periods, and then Theorem 2 gives $16\mid\operatorname{Tr}(JM^{12})$. The period-2 and period-4
   rings searched below have $|g|=12$ and $6$ (orbifold index 2 and 4), so they were never routes
   to index 1; their statistics are kept as data on 4-cell products. Two further obstructions:
   - $r=2$ is 3-adically dead for every ring with $p\ge2$ (§3);
   - 1-dimensional-type eigenvalues are divisible by 3 at $r=3$ (§3).

| c | ring | $g$, order, orbifold index | conditions | status |
|---|---|---|---|---|
| 144 | 6 identical $q{=}12$ cells | $z$, 12, **2** | $I_2=0,\ \lvert I_6\rvert=24$ | allowed; none in 48,728 cells |
| 144 | 12 identical $q{=}6$ cells | $z$, 24 → index 1 impossible; $z^2$, 12 → index 2 needs period 2 | see §13 | open ((6,2) family) |
| 144 | period-2 ring, 6 cells of size 12 | $z^2$, 6, **4** | $\operatorname{Tr}JM=0,\ \lvert\operatorname{Tr}JM^3\rvert=24$ | **23 found and verified** (index-4 theories) |
| 144 | 3 block cells of size 24 | $z$, 6, 4 | same | none in 13,859 cells |
| 288 | 12 identical $q{=}12$ cells | $z$, 24, 1 | $I_4=0,\ \lvert I_{12}\rvert=24$ | **impossible** (Thm 2) |
| 288 | any ring, any $q$, any period, $g$ a power of the shift | forces $L=12$ periods | | **impossible** (Prop. §13 + Thm 2) |
| 288 | period-2 ring of 12 cells | $z^2$, 12, 2 | $\operatorname{Tr}JM^2=0,\ \lvert\operatorname{Tr}JM^6\rvert=24$ | not index 1; none found anyway |
| 288 | period-4 ring of 12 cells | $z^4$, 6, 4 | $\operatorname{Tr}JM=0,\ \lvert\operatorname{Tr}JM^3\rvert=24$ | not index 1; none found anyway |
| 288 | cyclic ring, $g=-z^4$ | 6, 4 | $C_4=0,\ \lvert C_{12}\rvert=24$ | not index 1 (Part II, corrected) |

---

## 1. Setting

### 1.1 General cells

A cell is $(U,S)$: $S\in M_q(\mathbb F_3)$ with $SS^T=-I$, and $U=\operatorname{row}(B)$ a
non-degenerate $r$-dimensional subspace with orthogonal projection $\Pi$. Then
$$A(z)=(I+(z-1)\Pi)\,S\quad\text{over }\mathbb F_3[z]/(z^L\mp1),\qquad AA^*=-I .$$
The legacy form $Q\,\mathrm{diag}(zI_r,I)\,S_0$ is the case $B=(Q_{:,1..r})^T$, $S=QS_0$.
There $U$ has an orthonormal basis, i.e. square discriminant. **Non-square discriminant is
new and equally valid** (`core.py`). The local table is
$$T[c,c']=\sum_{a\in\{1,2\}^q,\ a\Pi=c'B}\phi(a)\,\phi\big((cB+a(I-\Pi))S\big),$$
with $I_L=\operatorname{Tr}(JT^L)$ for $z^L=-1$ and $C_L=\operatorname{Tr}(T^L)$ for the
cyclic twin $z^L=+1$, which is also self-dual of length $24L$.

### 1.2 Periodic rings

Cells with a **common** $U$ but different $S_j$ chain into a self-dual ring:
$A=P(z)\operatorname{diag}(S_0,\dots,S_{L-1})$ with $P$ the conveyor of $U$-components, so
$AA^*=-I$. Its index is $\operatorname{Tr}(J\,T_{S_0}\cdots T_{S_{L-1}})$
(`core.generator_ring`; brute force in `test_composite.py`). With period $p$,
$M=T_{S_1}\cdots T_{S_p}$ is the table of one cell of size $12p$, and the GJF symmetry is the
shift by $p$ cells:
- c=144: $p=2$, 6 cells, $g=$ shift by 2;
- c=288: $p=4$, 12 cells, $g=$ shift by 4; or $p=2$, 12 cells, $g=$ shift by 4.

## 2. Structure theorems (proved)

**Theorem 1 (Weil symmetry).** Let $F_U[c,c']=\omega^{cGc'^T}$ and
$D=\operatorname{diag}(\omega^{cGc^T})$. For every cell ($q$ is necessarily even), $F_UT=TF_U$.
- If $3\mid q$, $T$ is block diagonal by norm, so $DT=TD$ and $T$ commutes with the Weil
  representation of $\mathrm{SL}_2(\mathbb F_3)$.
- If $3\nmid q$, $DTD^{-1}$ is a scalar multiple of $T$.
- In all cases $T$ commutes with $J$ and with $Q_8=\langle F_U,DF_UD^{-1}\rangle$.
- Products of tables with the same $U$ inherit all of this.

> **Annotation (online note).** Read "$T$ commutes with $Q_8=\langle F_U,DF_UD^{-1}\rangle$" as: $T$ commutes with the
> matrices $F_U$ and $DF_UD^{-1}$. Since $F_U^2=3^rJ$, these generate an infinite matrix group; the finite $Q_8$ is
> their induced (conjugation) action on $W=U\oplus U$. When $3\nmid q$, $DTD^{-1}=cT$ and commutation with
> $DF_UD^{-1}$ follows: $(DF_UD^{-1})T=c^{-1}DF_UTD^{-1}=c^{-1}DTF_UD^{-1}=T(DF_UD^{-1})$.

*Proof.*
1. Put $f(x)=\phi(xS)$, $g=\phi$, and kernels $F[u,w]=f(u+w)$, $G[u',w]=g(u'+w)$ for
   $w\in U^\perp$. Then $T=FG^T$.
2. $\hat\chi=i\sqrt3\chi$ and $S^{-1}=-S^T$ make $f$ and $g$ eigenvectors of
   $F_q=F_U\otimes F_{U^\perp}$ with the same eigenvalue $\lambda=(i\sqrt3)^q$, i.e.
   $F_UFF_{U^\perp}=\lambda F$ and $F_UGF_{U^\perp}=\lambda G$.
3. Since $F_U$ and $F_{U^\perp}$ are symmetric, $F_UT=\lambda FF_{U^\perp}^{-1}G^T=TF_U$.
4. The norm statement follows from $\lvert a\rvert^2=\lvert y\rvert^2=q$ and $SS^T=-I$. ∎

Checked in `test_core.py`.

**Consequences.** $T=\bigoplus_\pi\mathrm{Id}_\pi\otimes T_\pi$ over the irreducibles of
$\mathrm{SL}_2(\mathbb F_3)$, of dimensions 1, 1, 1, 2, 2, 2, 3. $J$ acts by $(-1)^r$ on the
1- and 3-dimensional types and by $-(-1)^r$ on the 2-dimensional ones, so
$$I_L=(-1)^r(a_L+3b_L-2c_L),\qquad C_L=a_L+3b_L+2c_L,$$
where $a,b,c$ are the power sums of the 1-, 3- and 2-dimensional types. This is the origin of
the eigenvalue multiplicities 1, 2, 3 in every spectrum. The $K_4$ spectrum of the c=72 seeds
is the permutation representation $1\oplus3$ of $A_4=\mathrm{PSL}_2(\mathbb F_3)$ on
$\mathbb P^1(\mathbb F_3)$ (`weil.py`, `typepolys.py`).

**Theorem 2 (2-adic obstruction).** For every cell, and for every product $M$ of tables with
a common $U$ (period-$p$ ring with $L$ periods), and every $L\ge4$:
$$2^{\,v_2(L)+2}\ \big|\ \operatorname{Tr}(JM^L).$$

*Proof.*
1. **Weyl expansion.** Let $W=U\oplus U$ and
   $W_{(a,b)}e_x=\omega^{\langle x,b\rangle+\langle a,b\rangle/2}e_{x+a}$, a basis of
   $\operatorname{End}\mathbb C[U]$. Conjugation by $F_U$ and by $D$ *permutes* these
   operators, with no phases, via elements of $\mathrm{Sp}(W)$ generating $\mathrm{SL}_2(\mathbb F_3)$.
   So $M=\sum_wt_wW_w$, with $t$ constant on $Q_8$-orbits, $t_{-w}=t_w$, and
   $t_w=3^{-r}\operatorname{Tr}(W_w^\dagger M)\in\mathbb Z_2[\omega]$.

   > **Annotation (online note).** Here $Q_8\subset\mathrm{SL}_2(\mathbb F_3)\subset\mathrm{Sp}(W)$ is generated by the maps
   > $s$ and $usu^{-1}$ induced by conjugation with $F_U$ and with $DF_UD^{-1}$ ($u$ induced by $D$). The coefficients are
   > constant on $Q_8$-orbits because conjugation by $F_U$ and by $DF_UD^{-1}$ fixes $M$; they need not be invariant under
   > $u$ itself when $DMD^{-1}=cM$ with $c\neq1$. The theorem is about the anti-periodic closure $\operatorname{Tr}(JM^L)$ only.
2. **Sum over sequences.**
   $\operatorname{Tr}(JM^L)=\sum_{w\in W^L}\psi(w)$, with
   $\psi(w)=\prod t_{w_i}\operatorname{Tr}(JW_{w_1}\cdots W_{w_L})$.
3. **Group action.** $\psi$ is invariant under the negacyclic shift
   $(w_1,\dots,w_L)\mapsto(-w_2,\dots,-w_L,w_1)$, of order $2L$, and under diagonal $Q_8$.
   Together these form a group of order $8L$.
4. **Stabilisers.** $Q_8$ acts freely on $W\setminus0$. So for $w\neq0$ the only involutions
   that can fix $w$ are (shift by $L/2$)·(element of order 4). Two different ones multiply to
   something fixing only 0, and none of them is a square. Hence the 2-part of
   $\mathrm{Stab}(w)$ is at most 2, and $v_2\lvert\text{orbit}\rvert\ge v_2(L)+2$.

   > **Annotation (online note).** The map displayed in step 3, $R(w)=(-w_2,\dots,-w_L,w_1)$, satisfies
   > $R^L=(-1)^{L-1}\mathrm{id}$: it has order $L$, not $2L$, when $L$ is odd. Use instead the anti-periodic shift
   > $T(w)=(w_2,\dots,w_L,-w_1)=-R(w)$, with $T^L=-\mathrm{id}$ (order $2L$); $\psi$ is $T$-invariant by cyclicity of the trace
   > and $JW_wJ=W_{-w}$. Since $-\mathrm{id}$ is the central element of the diagonal $Q_8$, $\langle T,Q_8\rangle=\langle R,Q_8\rangle$ is
   > the central product $Q_8\circ\mathbb Z_{2L}$ of order $8L$, and step 4 holds for it: the involutions that can fix a
   > nonzero sequence are the $qT^{L/2}$ with $q$ of order 4, at most one of them lies in a stabiliser, and none is a square
   > (details.html C.3, step 3, has the lemma written out).
5. **Zero sequence.** It contributes $t_0^L$ with $t_0=3^{-r}\operatorname{Tr}M$. Here
   $\operatorname{Tr}M$ is the index of a cyclic ring, which is even because maximal words
   come in $\pm$ pairs.
6. Burnside over $\mathbb Z_2[\omega]$ finishes the proof. ∎

`weyl_check.py` verifies every ingredient for $r=1,2$:
- exact permutation of the Weyl operators;
- orbit-constant $t$;
- the free $Q_8$ action;
- the Burnside sum equals $\operatorname{Tr}(JT^L)$;
- invariance of $\psi$ on random sequences.

**Corollaries.**
- $16\mid I_{12}$ for 12 identical cells: no direct c=288 element there.
- Likewise $32\mid I_{24}$ for 24 cells of size 6.
- *Twisted* rings $S_i=P^iS$, with $P$ monomial preserving $U$, reduce to
  $K=T_S\Sigma$, which commutes with $Q_8$, so they are blocked too.
- Only genuinely different cells can work, with at most $v_2(L/p)\le1$.

**Sharpness.** $\min v_2(I_L)=\max(3,v_2(L)+2)$ for all $L\le24$. The data are 3,496 + 202
cells ($q=12$, $r=2..6$), 16,764 ($q=6$) and 148 ($q=18$), in `results/val*.json`.

## 3. Further (partly empirical) constraints

- **$3\nmid q$** (data on 27,126 / 21,118 / 642 cells for $q=8,10,16$): $I_L=0$ unless
  $3\mid L$, and $\min v_2(I_{3m})=m+1$ for $m\ge2$. This blocks c=288 with $(q,L)=(16,9)$
  and $(8,18)$, and c=144 with $(8,9)$. *Not proved.*
- **Type congruences** (600 exact cells): $a_L\equiv b_L\pmod{2^{v_2(L)+2}}$ and
  $c_L\equiv a_L+b_L\pmod{2^{v_2(L)+1}}$. These are 2-modular relations among the
  $\mathrm{SL}_2(\mathbb F_3)$-types, consistent with Theorem 2.
- **3-adic facts** (`v3_products.py`, `results/v3_products.json`).
  - At $r=3$, every eigenvalue of the 1-dimensional types is divisible by 3 (300/300 cells).
  - At $r=2$, $U$ non-square: the 1-type block lies in $3M_2(\mathbb Z)$ and the 2-type
    scalars lie in $3\mathbb Z$. Balance then forces $3^{2p}\mid\cdot$, so products with
    $p\ge2$ have $81\mid I_3$.
  - Observed minimum $v_3$ of balanced products: $r=2$ gives 4 ($p=2$) and 10 ($p=4$);
    $r=3,4$ give 1. **$r=2$ rings are dead; $r\ge3$ are not.**
- **Spectral design** (`spectral_design.py`). Impose all constraints above for $L\le24$ and
  search spectra with at most two monic factors of degree $\le2$ per type.
  - Uniform c=144 ($I_2=0,\lvert I_6\rvert=24$): 484 designs. The smallest is the $K_4$
    core $A=\{-3\},B=\{1\}$ with 2-type sector $\{\pm\sqrt6,\pm i\sqrt3\}$, giving
    $I_2=12-12$ and $I_6=732-756=-24$.
  - c=72-type conditions on a product with 1-type eigenvalues divisible by 9 or 81: designs
    exist, e.g. $B=\{-4\}$, $C=\{-4,-1\pm\sqrt3\}$. **Period-4 is not spectrally obstructed.**
  - Period-2 c=288 ($\operatorname{Tr}JM^2=0$, $\lvert\operatorname{Tr}JM^6\rvert=24$, $9\mid A$):
    no design with at most two factors per type.

## 4. Searches (all exact; 3 workers max, `nice 19`, time-bounded)

| script | family | cells / tuples | balanced | hits | smallest balanced $\lvert I\rvert$ |
|---|---|---|---|---|---|
| `search_shallow.py` ×2 | 6 identical $q{=}12$ cells, ≤6 reflections, both discriminants | 48,728 | 2,245 (1,723 degenerate) | 0 | 1,200 |
| `search_big.py` | block cells $q{=}24$: c=144 ($L{=}3$) / c=288 ($L{=}6$) | 13,859 (4,231 dead) | 752 / 28 | 0 / 0 | 408 / 3,289,248 |
| `search_big.py` | block cells $q{=}48$, c=288 ($L{=}3$) | 206 | 0 | 0 | – |
| `search_composite.py` ×2 | period-2, 6 cells (c=144, $\mathbb Z_6$, index 4) | 135,000 cells, 40.5M pairs | 2.78M | **46** (23 unique, verified) | 24 |
| ″ | period-2, 12 cells (c=288) | 40.5M pairs | 252k (191k degenerate) | 0 | 10,656 |
| ″ | period-4, 12 cells (c=288) | 4.05·10⁹ quads screened | 78.7M mod-$p$ candidates | 0 | 695,952 |
| `search_p4.py` | period-4, $r=3$ | 1.41·10⁸ quads (+3.5·10⁷ pairs) | 155,789 nondegenerate | 0 (c=144 control: 4) | 1,333,296 (period-2 c=288: 1,084,680) |
| `search_p4.py` | period-4, $r=4$ | 1.77·10⁷ quads | 25,972 nondegenerate | 0 (c=144 control: 2) | 827,616 |

Controls: the same runs saw 129 ordinary c=72 seeds, and no violation of 24-divisibility or of
Theorem 2 anywhere.

## 5. Verification (`verify_hit.py`)

Every reported hit is checked on the **expanded generator**:
- self-duality: rank and $GG^T=0$;
- the period shift preserves $C$; $g^3=-1$; $g^2$ is a fixed-point-free product of 3-cycles;
- the index by two independent exact computations: big-integer products of the $(B,S)$
  tables, and per-cell tables read directly off the generator matrix and multiplied around
  the ring;
- balance by **exhaustive** signed count over the $g^2$-fixed subcode ($2^{24}$ words).

Control: Codex's cell gives 24 / balance 0 at $L=3$, and 2,125,032 / balance 312 at $L=6$,
both matching.

**c=144 codes (index-4 orbifolds):** all 23 unique hits pass every check
(`results/c144_seeds_verified.json`, which has full $B,S_1,S_2$). Example:

```
r = 3, discriminant of U = 1, index = -24, Tr(J M^k), k=1..6 = [0, -96, -24, 960, -661080, 21212640], nonzero eigenvalues of M = T1 T2: 15
B (rows):   201102200012  020121200222  021221200220
S1 (rows):  120001021122  002212201220  211101002022  200202021122  112210021100  000222012211  020220122102  121020121010  001000210210  212101011002  122202002202  222002101021
S2 (rows):  122212111110  021201002221  100101202212  211211001100  010012211101  202002011000  222212012111  110201011011  001212000100  221112202111  000011001022  121020202022
ring: cells S1,S2,S1,S2,S1,S2 with common U; g = shift by two cells.
```

The 23 seeds split as: square discriminant 12 (index −24: 10, +24: 2) and non-square discriminant 11 (−24: 6, +24: 5). All have $r=3$; $M=T_1T_2$ has 15–27 nonzero eigenvalues, so none is hand-sized. Six further raw hits from `search_p4.py` (`results/p4_r*_s4*.json`) were not catalogued.

## 6. Files

| file | purpose |
|---|---|
| `core.py` | general cells (any $q$, any discriminant), transfer tables, exact CRT traces, generators of uniform and periodic rings |
| `test_core.py`, `test_composite.py` | agreement with legacy code and certificate; Fourier commutation; periodic-ring index by brute force |
| `weil.py`, `typepolys.py` | Weil group, isotypic decomposition, per-type characteristic polynomials |
| `weyl_check.py` | numerical check of every step of Theorem 2 |
| `survey.py`, `valuations.py`, `valuations_q.py`, `v3_products.py` | divisibility surveys |
| `spectral_design.py` | admissible virtual spectra for the targets |
| `seeds72.py` | type decomposition of the known c=72 seeds |
| `search_l6.py`, `search_shallow.py` | identical-cell $L=6$ searches |
| `bigcells.py`, `search_big.py` | block-diagonal cells of size 24 and 48 |
| `search_composite.py`, `search_p4.py`, `search_p4_descent.py` | periodic rings |
| `verify_hit.py`, `verify_c144_hits.py` | independent verification |
| `results/`, `logs/` | outputs |

## 7. Limitations and next steps

**Retracted consequence.** An earlier version of this note claimed that the 23 codes give
index-$\pm2$ theories at $c=144$ and hence the periodicity element with one symmetric square.
That was wrong: the period-2 ring's symmetry $g=z^2$ has order 6, so the orbifold is
$\mathbb Z_6$ and the index is 4 (independent check, 17 October; see item 2 above).

**Limitations.**
- No c=288 code was found. The period-4 search sampled ~4·10⁹ + 1.6·10⁸ tuples, but the balanced
  $\lvert\operatorname{Tr}JM^3\rvert$ values sit at $10^6$–$10^8$. The density near 24 is tiny, so this is not
  evidence of impossibility.
- Theorem 2 is proved, and every step is checked numerically.
- The $3\nmid q$ law, the type congruences and the $r=2$ 3-adic obstruction are empirical or only
  sketched.
- The orbifold step relies on GJF Thm 5.6, as before.
- Compute: up to three single-thread workers (`nice 19`); nothing is left running.

**Next steps.**
1. (Withdrawn: the c=144 codes are index-4 theories; see item 2 and §13.)
2. For a direct c=288 hit, go after period-4 rings at $r=3,4$ with cells that have small spectra:
   shallow reflection pools (`--dmax 2`), descent (`search_p4_descent.py`), or design-guided pairing
   toward the admissible spectra of §3 (e.g. $B=\{-4\}$, $C=\{-4,-1\pm\sqrt3\}$, $9\mid A$).
3. Prove the $3\nmid q$ law ($v_2(I_{3m})\ge m+1$) and the 3-adic statement for $r=2$. Both are
   probably the same Weyl/Burnside argument with the order-3 unipotent added.

---

# Part II — the cyclic twist (added 17 October, second pass)

*Goal of this pass: a work-around for Theorem 2 at c=288.*

> **Correction (17 October, after an independent check).** The cyclic symmetry
> $g=-z^{L/3}$ has order 6, so every cyclic-ring orbifold below is a $\mathbb Z_6$ orbifold of
> index 4, at $c=72$, 144 **and 288**. A cyclic ring cannot carry a $\mathbb Z_{12}$ or
> $\mathbb Z_{24}$ symmetry of the form $\pm z^k$ at all: $g^{m/2}=-1$ needs $-1\in\langle z\rangle$,
> and in the cyclic ring $z$ is a plain permutation. So Part II does **not** give a route to
> index 1; it is a (negative) study of $c=288$ index-4 code orbifolds with twelve identical cells.
> The structural results (§9: finite-order table, $r\ge4$) are unchanged as statements about that
> target.

## 8. Why the earlier search missed half of the construction

Every cell gives **two** rings: the negacyclic one ($z^L=-1$, index $I_L=\operatorname{Tr}JT^L$)
and the cyclic one ($z^L=+1$, index $C_L=\operatorname{Tr}T^L$). Both are self-dual of length
$24L$. The GJF symmetry exists in both cases:
$$g_{\text{nega}}=z^{L/3},\qquad g_{\text{cyc}}=-\,z^{L/3},\qquad g^3=-1\ \text{in both cases},$$
with $g^2=z^{2L/3}$ a fixed-point-free product of 3-cycles, and balance given by the fixed
words of $z^{L/3}$: $I_{L/3}=0$, resp. $C_{L/3}=0$. `verify_hit.verify_cyclic` checks all
of this on the expanded generator; on Codex's cell the cyclic 3-ring is itself a c=72 seed
($C_1=0$, $C_3=-24$, confirmed by a $2^{36}$ brute-force count: 86,120 − 86,144), and the
cyclic 6-ring has balance 336 = exhaustive count.

**Theorem 2 does not apply to the cyclic index.** Its proof uses the order-$2L$ negacyclic
shift of Weyl sequences; the cyclic shift has order $L$ and the stabiliser argument fails.
The data agree: $\min v_2(C_L)=3$ for every $L\le24$ (all surveys), whereas
$\min v_2(I_L)=\max(3,v_2(L)+2)$.

So a c=288 target with **twelve identical cells** exists in the cyclic ring (orbifold index 4, see the correction above), in the form
$$C_4=\operatorname{Tr}T^4=0,\qquad C_{12}=\operatorname{Tr}T^{12}=\pm24 .$$
(Likewise c=144: $(C_2,C_6)=(0,\pm24)$; c=72: $(C_1,C_3)=(0,\pm24)$.)

## 9. What the cyclic target forces: a finite-order table

`cyclic_design.py` enumerates virtual spectra by $\mathrm{SL}_2(\mathbb F_3)$-type
($C_L=a_L+3b_L+2c_L$, $I_L=\pm(a_L+3b_L-2c_L)$) built from small monic factors (linear,
quadratic, all $\Phi_n$ with $n\le36$, $x^k\pm c$), imposing $24\mid C_L,I_L$ for $L\le24$,
Theorem 2 for $I_L$, the type congruences of §3, and the target. Result (12 solutions at
27 states, 52 at 81 states, all of the same shape):

> **Every admissible spectrum has a finite-order core**: $T^3=\pm1$ on a 24-dimensional
> invariant subspace split as 3 (1-dim types) + 9 (3-dim type) + 12 (2-dim types), or
> $T^9=\pm1$ on 48 dimensions ($\Phi_9,\Phi_{18}$), possibly plus "spectator" factors
> ($x^9\pm8$, $x^{12}+4$) that do not contribute to $C_4,C_{12}$. Nothing else.

Two consequences:
- Such a cell is simultaneously a cyclic c=72 seed ($C_3=\pm24$ or 0) and a c=144 seed
  ($C_6=\pm24$): $C_L=24\,[3\mid L]$ up to signs.
- The 1-dim types must carry unimodular eigenvalues $\{1,\omega,\bar\omega\}$. But at
  $r\le3$ every 1-dim-type eigenvalue is divisible by 3 (215/215 non-nilpotent cells, and at
  $r=2$ also all 2-dim-type eigenvalues); 3-adic units in 1-dim types first appear at $r=4$
  (59/59 cells). **So the cyclic c=288 cell needs $r\ge4$** (81 states; type dimensions
  11 / 40 / 30).

The cyclic twist therefore converts the c=288 problem into a sharp construction problem:
*find a 12-coordinate cell with $r\ge4$ whose local table has finite order on a
24-dimensional range and is nilpotent on the rest.*

## 10. Searches for such a cell (all negative)

| pool | cells | min $\rho(T)$ | $C_4=0$ | smallest balanced $\lvert C_{12}\rvert$ |
|---|---|---|---|---|
| random $q{=}12$, $r=2..5$ (`cyclic_scan.py`) | 14,124 | 3.0 | 5 | 1,326,240 |
| random $q{=}12$, $r=4$ / $r=5$ (ρ survey) | 2,964 / 870 | 4.15 / nilpotent only | – | – |
| random $q{=}6$, $r=1..4$ | 159,103 | – | 0 | – (values quantised: $C_6\in\{840,53760\}$) |
| Paley12 / Conf6x2 / A0 with order-3-orbit $U$ (`sym_cells.py`) | 18,197 | 3.9 | 2 | $4.1\cdot10^9$ |
| free-fermion $S=A_0,A_0D$ with sparse $U$ (`sparse_ff.py`) | 2,800 live of 15,000 | > 3 | 0 | – |

No table with a unimodular (finite-order) spectrum was seen anywhere; $\rho(T)\ge3$ in
every non-nilpotent cell. For comparison the cyclic c=72 target is hit at the same rate as the
negacyclic one (82 vs 87 seeds in the same 14,124 cells), and the cyclic c=144 balanced values
are as small as 1,728, so the method itself works; only the finite-order requirement is out of
reach of sampling.

Non-identical cells with the cyclic twist (period 2 and 4; conditions $(\operatorname{Tr}M^2,
\operatorname{Tr}M^6)=(0,\pm24)$ and $(\operatorname{Tr}M,\operatorname{Tr}M^3)=(0,\pm24)$)
were probed with `search_p4.py --cyclic` ($r=3$: $4.9\cdot10^7$ quads, $1.2\cdot10^7$ pairs; $r=4$: $6.2\cdot10^6$ quads): no hit; smallest balanced values $\lvert\operatorname{Tr}M^3\rvert=2{,}106{,}816$ (period 4) and $\lvert\operatorname{Tr}M^6\rvert=4{,}271{,}544$ (period 2), the same scale as the negacyclic versions. The cyclic period-2 **c=144** control did hit (2 more seeds, not catalogued), so the machinery is fine; the scale of 4-cell products is the problem in both twists.

## 11. The construction problem, made concrete

For $S=A_0$ (six tetracodes) the table is explicit: $\phi(xA_0)\ne0$ iff every coordinate pair
of $x$ has exactly one zero, and then
$$T[u,u+d]=\sum_{\substack{y\in\{1,2\}^{12}\\ (y-d)\ \text{one zero per pair}\\ (y-d)\Pi=u}}
\phi(y)\,(-1)^{\#\{\text{pairs whose first coordinate of }y-d\text{ is }0\}},$$
and a direction $d\in U$ is *live* iff $d$ has no empty pair; each live $d$ contributes exactly
$2^6=64$ terms, pushed forward to $U$. A finite-order table needs the pushforwards of a few live
directions to cancel down to a signed permutation-like structure on 24 states with 3-cycles
$u\to u+d_1\to u+d_1+d_2\to u$ ($d_1+d_2+d_3=0$, all live) of consistent sign. In the
sampled sparse $U$ this never happened (most $U$ have only dead directions, the rest have
$\rho>3$), but the search space of structured $U$ (dimension 4, spanned by one-zero-per-pair
vectors forming hyperbolic planes) is finite and small enough for an exhaustive treatment
with the automorphism group of $A_0$ ($S_6\wr$ pair symmetries) factored out. That, or a
proof that $\rho(T)\ge3$ for every non-nilpotent cell (which would close the cyclic route
too), is the recommended next step.

## 12. Status after the second pass (corrected)

- **Proved:** 2-adic obstruction for identical cells with the negacyclic twist (Thm 2).
- **Cyclic twist:** GJF symmetry verified, but of order 6 only (orbifold index 4, not 1); its
  c=288 cell would need a finite-order table of rank 24 or 48 and $r\ge4$ (classification within
  the factor families, plus the empirical 3-divisibility at $r\le3$). Not a route to the
  periodicity element.
- **Found and verified:** 23 c=144 codes with index $\pm24$ and a balanced $\mathbb Z_6$ symmetry
  (negacyclic period-2 rings): index-4 theories, not index 2 (§7, retracted consequence).
- **Not found:** any c=288 cell or ring, and by the Proposition of §13 no ring whose GJF symmetry
  is a power of the shift can give index 1. What is left for a direct c=288 element are
  symmetries that are not powers of the shift (internal symmetries of the cells, §13.3).

---

# Part III — the corrected targets (17 October, after an independent check)

## 13. Symmetry orders, and what is actually open

GJF Theorem 5.6 at level $k$ needs a symmetry of order $m=\gcd(24,k)$ with $g^{m/2}=-1$,
$g^{m/3}$ of shape $3^{n/3}$ and balanced; the orbifold index is $24/m$. So:
$c=72$: $m=6$ (index 4); $c=144$: $m=12$ (index 2); $c=288$: $m=24$ (index 1).
The earlier parts of this note used $g$ of order 6 throughout, which is only right at $c=72$.

### 13.1 Proposition (no shift-power $\mathbb Z_{24}$ ring at c=288)

Let a ring have $N$ cells of size $q$ ($qN=144$), common $U$, period $p$, and let the GJF
symmetry be a power $g=z^{pj}$ of the cell shift (negacyclic, $z^N=-1$). Then $|g|=24$ and
$g^{12}=-1$ force $z^{12pj}=z^N$, i.e. $12pj\equiv N\pmod{2N}$, and $\operatorname{ord}(z^p)=
2N/\gcd(2N,p)$ divisible by 24. Writing $p=Nu/12$ gives $\gcd(24,u)=1$ and $L=N/p=12/u$, so
$u=1$ and the ring has **exactly 12 periods**. Theorem 2 then gives $16\mid\operatorname{Tr}(JM^{12})$.
*Hence no transfer ring of any cell size or period, with a GJF symmetry that is a power of the
shift, is a direct c=288 index-1 code orbifold.* (Observed in an independent check; it subsumes
the "routes" of Part I.) The cyclic twist cannot help: $\pm z^k$ never has $g^{m/2}=-1$ for
$m\in\{12,24\}$.

### 13.2 The corrected c=144 index-2 target

At $c=144$ the same count forces $L=6$ periods, so $(q,p)\in\{(12,1),(6,2),(4,3),(2,6)\}$,
$g=z^p$ of order 12, balance $\operatorname{Tr}(JM^2)=0$, index $\operatorname{Tr}(JM^6)=\pm24$,
$M=T_{S_1}\cdots T_{S_p}$. Theorem 2 only gives $8\mid\operatorname{Tr}(JM^6)$: the target is open.
- $(2,6)$: **dead** — every 2-coordinate cell is a tetracode with a zero table (no maximal words),
  checked exhaustively over all 8 anti-orthogonal $2\times2$ matrices and all $u$
  (`results/c144_q2_p6_exhaustive.json`).
- $(12,1)$: 48,728 cells before (smallest balanced $\lvert I_6\rvert=1200$); deeper pools were then
  searched by a persistent worker (§14).
- $(6,2)$ and $(4,3)$: new; searched by persistent workers (§14). `verify_hit.verify_ring12` is the verifier (g of order
  exactly 12, $g^6=-1$, $g^4$ fixed-point-free 3-cycles, exhaustive $2^{24}$ balance on the
  $g^4$-fixed subcode, index by two methods); it reproduces Codex's cell as a $(12,1)$ ring
  (index 2,125,032, balance 312 = exhaustive).

Results: see §14.

### 13.3 Internal symmetries: the negacirculant version is dead

The Proposition assumes $g\in\langle z\rangle$. A ring also has symmetries acting *inside* the
cells: $h=(P\text{ on every message block},\,P'_j\text{ on check block }j)$ is an automorphism iff
$P$ preserves $U$ with the same action on $U$ in every cell and $PS_j=S_jP'_j$. A purely internal
$\mathbb Z_{24}$ needs $P^{12}=-1$ with $P^8$ fixed-point-free of order 3, i.e. $P$ a signed
12-cycle, $S_j$ negacirculant and $U$ an ideal of $\mathbb F_3[x]/(x^{12}+1)$. There are 648
anti-orthogonal negacirculants, but **every proper ideal is degenerate** (the dual of the ideal
$(f^ag^b)$ is $(f^{3-b}g^{3-a})$, and they always intersect), so no cell of this kind exists
(`internal_sym.py`, `results/internal_sym.json`).

### 13.4 Internal × shift: an order-8 internal symmetry on a period-2 ring (open, index 1)

Combine an internal symmetry with the shift: $g=h\,z^{p}$ with $h=(\sigma,\sigma)$ on every
block. For a period-$p$ ring of $N=6p$ cells (so $z^{12p}=z^{2N}=1$) and $\sigma^4=-1$,
$$g^{12}=\sigma^{12}z^{12p}=\sigma^4=-1,\qquad g^{8}=z^{8p}=-z^{2p},\qquad g^{24}=1,$$
so $g$ has order 24, $g^8$ is a fixed-point-free product of 3-cycles (it moves every cell), and the
$g^8$-fixed words are the 2-period negacyclic ring. For $q=12$, $N=12$, $p=2$ this is a length-288
code with a $\mathbb Z_{24}$ symmetry, **orbifold index 1**, and the conditions are exactly the
period-2 ones:
$$\operatorname{Tr}(JM^2)=0,\qquad \operatorname{Tr}(JM^6)=\pm24,\qquad M=T_{S_1}T_{S_2}.$$
Theorem 2 gives only $8\mid\operatorname{Tr}(JM^6)$, so nothing forbids $\pm24$. The only new
requirement is that both cells commute with $\sigma$ and $U$ is $\sigma$-invariant:
- $\sigma^4=-1$ forces $\sigma$ to be a product of three signed 4-cycles with sign product $-1$;
  take $\sigma=I_3\otimes P_4$. Over $\mathbb F_3$, $x^4+1=(x^2+x+2)(x^2+2x+2)$ is separable, so the
  centraliser of $\sigma$ is $M_3(\mathbb F_9)\times M_3(\mathbb F_9)$, non-degenerate invariant $U$
  have $r=4$ or 8, and anti-orthogonal equivariant $S$ exist ($S=O\,(I_3\otimes s_4)$ with $s_4$ an
  equivariant anti-orthogonal $4\times4$ matrix and $O$ a Cayley transform of a skew element of
  the centraliser; `search_sigma.py`).
- $T$ then commutes with the order-8 permutation of states induced by $\sigma$, on top of the
  Weil symmetry.

**Closed (17 October, later):** σ-equivariant cells have every ring index divisible by 16 (NOTE_kspan.md §4.1), so $\lvert\operatorname{Tr}JM^6\rvert=24$ is impossible in this family; the route is dead.

`verify_hit.verify_sigma_ring` checks everything on the expanded length-288 generator (both $h$
and $g$ preserve the code, $g^{12}=-1$, order exactly 24, $g^8$ of shape $3^{96}$, index and
balance each by two independent exact methods; the $g^8$-fixed subcode has dimension 48, so no
exhaustive count). Results of the first probe: see §14.

## 14. Addendum: the corrected c=144 index-2 search (state of 17 October)

| family | status | smallest balanced $\lvert\operatorname{Tr}JM^6\rvert$ |
|---|---|---|
| $(2,6)$ | **dead**: all 8 anti-orthogonal $2\times2$ cells have $T=0$ (tetracode, no maximal words) | – |
| $(4,3)$ | **dead**, exhaustive (`c144_q4_exhaustive.py`): 1152 anti-orthogonal $4\times4$ matrices, 138 non-degenerate $U$, at most 6 distinct tables per $U$, all 2,304 ordered triples; every balanced triple has $\operatorname{Tr}JM^6=0$ | none nonzero |
| $(6,2)$ | persistent worker (pools of distinct tables per random $U$, $r=2..4$, all pairs) | 134,784 after 224 pools / 9.5M pairs (values coarsely quantised: 134,784; 245,760; 436,320 recur; most pools have only 16–50 distinct tables) |
| $(12,1)$ | persistent worker (deeper pools: up to 20 reflections, $r=2..6$) | **48** after 382,200 deeper cells (earlier 1,200 in 48,728 shallow cells); worker stopped 17 October 17:05 |

**c=288, index 1, via §13.4 (σ-equivariant period-2 rings):** first probe `search_sigma.py`
(6 min): 72 invariant $U$ ($r=4$), 10,800 equivariant cells, 1.62M ordered pairs, 1,161
balanced, 840 non-degenerate, smallest balanced $\lvert\operatorname{Tr}JM^6\rvert=175{,}370{,}688$
(`results/sigma_a.json`). The single cells have $\rho(T)\approx12$–18 and 40–48 nonzero
eigenvalues, so 2-cell products sit at the $10^8$ scale: the route is open but, like all 4-cell
products before, needs structure rather than sampling. A random equivariant pair was run through
`verify_sigma_ring` as a control: order 24, $g^{12}=-1$, $g^8$ of shape $3^{96}$, index
$1.5\cdot10^{12}$ and balance 14,592, each by two methods.

Any hit is to be verified with `verify_hit.verify_ring12` (order of $g$ exactly 12, $g^6=-1$,
$g^4$ fixed-point-free 3-cycles, exhaustive $2^{24}$ balance, index by two methods).
The workers' state files are in the project workspace (not included).
