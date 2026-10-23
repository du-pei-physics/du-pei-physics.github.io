# Low-weight field content of T (c=72), T′ (c=144) and P (c=288)

Scripts: `../code/fc_weight_enum.py` (weight enumerators, ~4 s), `../code/fc_characters.py` (characters of V and T,
~1 s), `../code/fc_theta.py` (theta series), `../code/fc_sym2.py` (T′ and P), `../code/fc_closed_form.py` (§6). Run them
from `../code/` in this order: `python3 fc_weight_enum.py`, `python3 fc_characters.py`, then `fc_theta.py`, `fc_sym2.py`,
`fc_closed_form.py`, `fc_modular_check.py`. They write the files `fc_*.json` named below into `../code/out/` (these are
not shipped; the input `canonical_seed.json` is read from `../certificates/original/`). All arithmetic is exact: Python
integers, and Z[ω] for the twisted–twined functions. Labels: **[proved]** = proved here or by an elementary argument
given here; **[computed]** = exact machine computation, with the cross-checks listed; **[cited]** = taken from the
literature; **[assumed]** = a convention or input that is not derived here.

## 0. Summary

| theory | c | Witten index | dim V₁/₂ | dim V₁ | dim V₃/₂ | dim V₂ | dim V₅/₂ | dim V₃ | Ramond ground states (even / odd) |
|---|---|---|---|---|---|---|---|---|---|
| V = V_Λ(C) (seed) | 72 | −24 | 0 | 72 | 678 | 15 012 | 252 006 | 5 834 120 | 41 488 / 41 512 |
| T = V//Z₆ | 72 | −4 | 0 | **0** | **113** | **2 928** | 42 001 | 962 444 | **6 968 / 6 972** |
| T′ = Sym²T | 144 | 2 | 0 | **0** | **113** | **2 928** | 42 001 | 968 772 | **48 587 870 / 48 587 868** |
| P = Sym²(T′⊗Arf) | 288 | −1 | 0 | **0** | **113** | **2 928** | 42 001 | 975 100 | **2 360 781 062 549 029 / 2 360 781 062 549 030** |

All rows are **[computed]**. The rows for T, T′ and P use the lift and orbifold inputs A1–A2 of §2 (the standard
construction of GJF; §2 says which parts are proved and which are cited).

* **Vacuum.** The NS vacuum is unique in T, T′ and P **[proved]**. Every twisted sector has positive ground energy:
  8/3, 19/6 and 9/2 for the g², g and g³ sectors of V; c/12 = 6 and 12 for the swap-twisted sectors.
* **No currents, no free fermions.** dim V₁/₂ = dim V₁ = 0 for T, T′ and P **[proved]**, for the paper's seed. The
  weight-1 Lie algebra is zero. Three facts are needed: the seed code has minimum distance 9 (A₃ = A₆ = 0, §1.1), so Λ
  has no vectors of norm 1 or 2 and there are no lattice states of weight 1/2 or 1; the 72 Heisenberg currents of V are
  removed by the projection, because g has no eigenvalue 1 on R⁷²; and every twisted sector starts above weight 1.
  This is particular to this seed: of the 8,748 k = 6 codes, 2,701 give a T with weight-1 fields (details.html §G.3).
* **Weight 3/2.** There are 113 weight-3/2 fields, all of them D₄-symmetric copies of the fields of T. They include
  the supercurrent G.
* **Ramond vacuum is not unique.** P has index −1 but about 4.7·10¹⁵ Ramond ground states. The number of odd
  ground states exceeds the number of even ones by exactly one. T itself has 13 940 Ramond ground states.
* **Weights below 3 are untwisted.** Below weight 3, every field of T lies in the untwisted lattice SVOA sector: the
  Z₆-invariant states of V, including oscillator polynomials (at weight 2, for example, the 876 = (2 700 + 2 556)/6
  invariant Heisenberg bilinears besides the momentum states). The first twisted-sector fields appear at weight 3:
  16 + 16 states from the g²- and g⁴-twisted sectors. The swap-twisted sectors of T′ and P start only at weights 6 and
  12, so below weight 3 the fields of P are untwisted as well (§4); the symmetric squares do add untwisted two-copy
  states, the first ones at weight 3 (C(113, 2) = 6 328 at each step: 962 444 → 968 772 → 975 100).

### What "a description of P in terms of fields" can and cannot mean

P is a holomorphic SVOA with c = 288. It has no weight-1 currents, so there is no Kac–Moody or free-field symmetry
to organise it, and an explicit list of strong generators with their OPEs is not available. Such a list would have to
start with the 113 fields of weight 3/2 and the 2 928 of weight 2, and nothing forces it to stop at a small weight.
What *is* explicit is the following:

1. **A concrete free-field realisation of every field of weight < 3.** Below weight 3, the fields of P are the
   D₄-symmetrised, Z₆⁴-invariant vertex operators :P(∂X)e^{iv·X}: of 288 chiral bosons X = (X₁,…,X₄), Xₐ ∈ ℝ⁷²,
   with momenta in Λ(C)⁴. Here Λ(C) is the Construction-A lattice of the seed code. The counts of these fields
   are exact (§§3–4). At weight 3/2 they come from the 144 vectors ±√3eᵢ and the 534 minimal vectors of Λ(C)
   (one for each weight-9 codeword); Z₆ acts freely on both sets, which gives 24 + 89 = 113 orbit sums.
   The supercurrent G = Σᵢ(e^{√3eᵢ·X} + e^{−√3eᵢ·X}) (with cocycle factors; GJF §2, **[cited]**) is one combination of the 24 orbit sums
   of the first kind. At weight 2 there are 876 Heisenberg bilinears and 2 052 orbit sums of norm-4 vectors.
2. **The first twist fields** appear at weight 3, from the g²- and g⁴-twisted sectors of V: 16 + 16 of them.
   They are the combinations |w⟩ + |−w⟩ of the momentum states w ∈ L* of norm 2/3. These momenta correspond to the
   32 weight-6 words of the contracted fixed code C̄ (§1.2, §3.1). The swap-twisted sectors of T′ and P begin only at
   weights 6 and 12.
3. **A closed form for the NS character** (§6): Z_T = t⁶ − 144t⁵ + 6984t⁴ − 129679t³ + 751404t² − 629867t + 13940
   with t = (ϑ₃/η)¹². Z_{T′} and Z_P follow from it by two explicit Sym² formulas.
4. **The full graded dimensions** of NS and R are given to weight 5 in `fc_sym2.json` / `fc_characters.json`.

So "in terms of fields" means a **lattice-plus-orbifold description**. P is the D₄ ⋉ Z₆⁴ orbifold of the rank-288
lattice SVOA V_{Λ(C)}^{⊗4}, with the twist of NOTES §2.3, and it is described by its untwisted and twisted field
content as above. It does not mean a closed generators-and-relations presentation.

## 1. The seed code C and its lattice

### 1.1 Weight enumerator **[computed]**

C is the self-dual [72,36] code of `canonical_seed.json`, with generator G = [I₃₆ | A] and
A = [[D,E,0],[0,D,E],[−E,0,D]]. Write a codeword as (a₀,a₁,a₂ | c₀,c₁,c₂) with aⱼ ∈ F₃¹² and
cⱼ = aⱼD + aⱼ₋₁E, where the wrap-around term is c₀ = a₀D − a₂E. Since E = Q_R S_R has rank 3, the only information
passed from period j to period j+1 is sⱼ = aⱼQ_R ∈ F₃³ (27 states). The transfer matrix
M_y[t,s](x) = Σ_{a: aQ_R = s} x^{wt(a)+wt(aD + tS_R)} y^{#2's} gives W_y(x) = Σ_s (M_y³)[−s, s]
(`fc_weight_enum.py`). Here y = 1 gives the weight enumerator and y = −1 the signed enumerator; the sign of a word
is (−1)^{number of 2's}.

Checks: the coefficients sum to 3³⁶. The MacWilliams identity W(x) = 3⁻³⁶ Σ A_w(1+2x)^{72−w}(1−x)^w holds
exactly (checked with sympy). A₇₂ = 83 000 with signed count −24, i.e. 41 488 even and 41 512 odd full-support words.
These agree with `canonical_seed.json` and with the banded index DP.

| w | A_w | signed A_w |
|---|---|---|
| 0 | 1 | 1 |
| 9 | 534 | 0 |
| 12 | 7 506 | 62 |
| 15 | 93 894 | 0 |
| 18 | 2 851 406 | 606 |
| 21 | 90 800 514 | 0 |
| 24 | 3 290 190 426 | −51 754 |
| 27 | 127 951 584 574 | 0 |
| 30 | 3 515 424 780 654 | 1 335 234 |
| 33 | 59 322 462 785 676 | 0 |
| 36 | 607 959 503 674 056 | 22 606 072 |
| 39 | 3 798 850 811 948 154 | 0 |
| 42 | 14 442 702 248 266 890 | −63 306 562 |
| 45 | 33 061 834 666 406 330 | 0 |
| 48 | 44 726 944 083 496 866 | 181 059 846 |
| 51 | 34 776 403 818 646 626 | 0 |
| 54 | 14 919 344 347 387 234 | 305 075 206 |
| 57 | 3 327 482 041 224 786 | 0 |
| 60 | 354 229 916 764 050 | 12 784 646 |
| 63 | 15 689 452 318 512 | 0 |
| 66 | 224 546 663 160 | 153 932 |
| 69 | 637 024 272 | 0 |
| 72 | 83 000 | −24 |

The minimum distance is **9** (A₃ = A₆ = 0). The signed counts vanish in odd weight, because c ↦ −c maps
#2 to w − #2.

### 1.2 The symmetry g and the contracted fixed code **[computed]**

g is the negacyclic shift by one period on both halves (`run.shift_monomial([12],3,False,1)`). It has order 6,
g³ = −1 and Frame shape 6²⁴3⁻²⁴. The element g² consists of 24 signed 3-cycles, each with sign pattern (+,−,−) and
sign product +1. Its Frame shape is 3²⁴, and it has the coordinate orbits {i, i+12, i+24} on each half.

C̄ is the contraction of the g²-fixed subcode to one coordinate per 3-cycle. It is a self-dual [24,12] code with
weight enumerator 1 + 32x⁶ + 3 856x⁹ + 62 304x¹² + 242 240x¹⁵ + 198 832x¹⁸ + 24 096x²¹ + 80x²⁴.
Its signed enumerator is 1 − 4x⁶ + 360x¹² + 44x¹⁸. In particular C̄ has minimum distance 6 and **index 0 by
cancellation**: 40 even and 40 odd full-support words. This differs from the witness seed of NOTES §2.1, whose C̄
has weight-3 words.

### 1.3 Theta series **[computed]** (`fc_theta.py`)

Λ = Λ(C) = {x/√3 : x ∈ Z⁷², x mod 3 ∈ C} is odd unimodular with norm v² = Σxᵢ²/3. A codeword of weight w lifts to
vectors whose coordinates are ±1, ±2, ±4, … on the support (one lift of each square 1, 4, 16, 25, … per residue)
and 0, ±3, ±6, … off the support. Hence Θ_Λ = Σ_w A_w θ₀^{72−w} θ₁^w, with θ₀ = Σ_{x∈3Z} q^{x²/6} and
θ₁ = Σ_{x≡1(3)} q^{x²/6}.

| norm | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| #Λ | 1 | 0 | 0 | 678 | 12 312 | 203 190 | 4 877 576 | 155 209 572 | 5 703 710 448 |

The low norms decompose by hand as follows:
* 678 = 144 (the vectors ±√3eᵢ) + 534 (weight-9 words lifted with ±1).
* 12 312 = 7 506 (weight 12) + 9·534 (weight 9 with one entry lifted to ∓2).
* 203 190 = 93 894 + 12·7 506 + C(9,2)·534.

The Ramond coset Λ + χ/2, with χ = √3(1,…,1), has 83 000 vectors of norm 6 (one for each full-support word), none of
norm 7, and 5 102 170 176 of norm 8.

The fixed lattice is L = Λ^{g²} ≅ {b ∈ Z²⁴ : b mod 3 ∈ C̄} with norm Σb², i.e. √3·Λ(C̄). It has det L = 3²⁴,
L* = L/3, and theta series 1 + 32q³ + … (no vectors of norm 1–5, 32 of norm 6).

The dual lattice L* has counts 1 at norm 0, 32 at norm 2/3, 4 096 at norm 1 and 97 488 at norm 4/3. The coset
L* + χ_L/2 has 80 vectors at its minimal norm 2/3.

## 2. Conventions and inputs

* **NS and R sectors of V** **[cited: GJF §2; NOTES §1.3]**
  * NS sector: V_Λ itself, with (−1)^F acting on e^v as (−1)^{v²} = e^{πi⟨v,χ⟩}.
  * R sector: V_{Λ+χ/2}, with (−1)^F acting on e^p as e^{πi⟨p,χ⟩}. The overall sign is fixed so that (−1)^F equals
    (−1)^{#2's} on the ground states e^p, p = (x+3/2·1)/√3, which correspond to the full-support words. Then
    I(V) = Index(C) = −24.
  * Weights are h = v²/2 + oscillator level. Ramond ground states have h = c/24 = 3.
  * The Kitaev/Arf sign of this choice is the only convention that matters (NOTES §2.5).
* **Twisted–twined functions.** Z^{(a,b)}_u = Tr_{g^a-twisted, spin u}(ĝ^b q^{L₀−c/24}), where u ∈ {AA, AP, PA, PP}
  gives the (space, time) periodicities: AA = NS, AP = NS with (−1)^F, PA = R, PP = R with (−1)^F. Encode the spins
  as AA=(1,1), AP=(1,0), PA=(0,1), PP=(0,0). SL₂(Z) acts on the right on both the holonomy pair and the spin
  vector: Z^x_u(γτ) = Z^{xγ}_{uγ}(τ).
* **Untwisted inputs** **[proved, elementary]**. Each Z^{(0,j)} is the momentum sum over the g^j-fixed vectors
  divided by the oscillator trace 1/η_{g^j}.
  * j = 1, 5: Δ(3τ)/Δ(6τ) in NS, both with and without (−1)^F. No momentum is fixed, and the Fock space is bosonic.
    In R the function is 0, since no vector of Λ+χ/2 is fixed.
  * j = 3 (g³ = −1): Δ(τ)³/Δ(2τ)³ in NS and 0 in R.
  * j = 2, 4: Θ_{L,·}/Δ(3τ), with the four spin-structure versions of the fixed-coset theta series
    (L, L with (−1)^{Σb}, L+(3/2)1, L+(3/2)1 with sign).
* **A1 [momentum states: proved, given GJF Thm 3.1; Ramond ground states: cited, GJF proof of Thm 5.6].** The lift ĝ² of the order-3
  element g² acts trivially on e^v for v ∈ L and on the g²-fixed Ramond ground states (the standard lift of an odd-order
  element). These are the only lift phases the computation uses, in Z^{(0,2)} and Z^{(0,4)}; ĝ itself is not trivial
  there, since g³ = −1 gives gv = −v for every g²-fixed v (e.g. v = (x, −x, x) on the three cells). For the momentum
  states: the lift lies in Λ̂[2].Aut(C) (GJF Thm 3.1), so ĝ² multiplies each e^v, v ∈ L, by a sign, and ĝ² has order 3,
  so the sign is +1. For the Ramond ground states: in the proof of their Thm 5.6, GJF note that for an element of odd
  order there is no sign ambiguity when lifting its action, and that its action on the Ramond ground states lifts the
  permutation action on the maximal (full-support) codewords; for h = g² this is triviality on the g²-fixed ground
  states. (Earlier versions of this note called this part an assumed normalisation.) The consistency checks below
  pass. Separately, ĝ³ is the
  standard lift of −1 (GJF Thm 3.1/5.6).
  *Correction (20 October, after Codex's audit):* an earlier wording said ĝ instead of ĝ², which is false as stated;
  the computation always used ĝ² (code/fc_characters.py), and no number changes.
* **A2 [cited: the standard anomaly-free orbifold construction, GJF §4, Thm 5.6].** The action of Z₆ on the twisted sectors
  is the one fixed by modular covariance, with no extra phases. Every function is then determined as follows.
  * Pairs generating Z₆: Z^{(a,b)} = Z^{(0,1)}(γτ) for any γ with second row ≡ (a,b) mod 6. The value uses
    Δ(kγτ) = ((cτ+d)/D_k)¹² Δ((A_kτ+B_k)/D_k), from the Hermite normal form of [[kα,kβ],[c,d]]. Because all
    exponents are multiples of 24, no multiplier phases appear.
  * Pairs in ⟨g²⟩ and ⟨g³⟩: obtained with S and T, using Poisson summation for θ_{L+μ,ν}.

  The resulting g²-twisted sector is Θ_{L*}/η(τ/3)²⁴. Its ground energy is 3 − 1/3 = 8/3, and its defect
  multiplicity is 3¹²/√det L = 1.
* **Spin refinement.** The Z₂ ⊂ Z₆ can be gauged with or without the Kitaev SPT sign (−1)^{Arf(u+ā)+Arf(u)}. That
  sign multiplies only functions that vanish identically here (check `kitaev_sign_only_on_vanishing_functions`).
  The two choices therefore give **identical characters**, so nothing below depends on this choice **[computed]**.

**Consistency checks, all passed [computed]** (`fc_characters.json → checks`):
1. Z^{(0,2)}_{PP} ≡ 0 to the computed order. This is the twined Witten index, 40 − 40.
2. Level matching for the g², g⁴ and g³ sectors (T³ and T² return the same function, with the correct spin).
3. Independence of the order-6 functions from the choice of γ (three different γ per pair).
4. Every projection (1/6)Σ_b Z^{(a,b)} has integer, non-negative coefficients. NS weights lie in ½Z, and the
   super-character equals (−1)^{2h} times the character (spin–statistics). R weights are integers.

## 3. V and T = V//Z₆ (c = 72)

### 3.1 Characters **[computed]**

The table gives graded dimensions (weight h, NS sector). V's super-character is (−1)^{2h} times its character.
The g-twisted columns are the raw twisted sectors before projection; only the Z₆-invariant part enters T.

| h | V | T untwisted (V^{Z₆}) | T from g²-tw. | T from g⁴-tw. | T from g, g⁵, g³-tw. | **T** |
|---|---|---|---|---|---|---|
| 0 | 1 | 1 | 0 | 0 | 0 | **1** |
| 1/2 | 0 | 0 | 0 | 0 | 0 | **0** |
| 1 | 72 | 0 | 0 | 0 | 0 | **0** |
| 3/2 | 678 | 113 | 0 | 0 | 0 | **113** |
| 2 | 15 012 | 2 928 | 0 | 0 | 0 | **2 928** |
| 5/2 | 252 006 | 42 001 | 0 | 0 | 0 | **42 001** |
| 3 | 5 834 120 | 962 412 | 16 | 16 | 0 | **962 444** |
| 7/2 | 171 669 852 | 28 611 642 | 614 400 | 614 400 | 614 400 + 614 400 + 0 | **31 069 242** |
| 4 | 6 089 550 870 | — | — | — | — | **1 879 408 646** |
| 9/2 | | | | | | **136 784 011 042** |
| 5 | | | | | | **6 484 810 801 152** |

Raw twisted sectors (before projection; L₀ values):
* g-twisted: 4 096 at 19/6 (this is 2¹² = √det(1−g)), 98 304 at 10/3, 1 228 800 at 7/2.
* g²-twisted: 1 at 8/3, 56 at 3 (24 oscillator states α₋₁/₃ⁱ|0⟩ plus 32 momentum states |w⟩, w ∈ L*, w² = 2/3),
  4 096 at 19/6.
* g³-twisted: 2³⁶ at 9/2.

**Hand checks [proved].**
* *T₁ = 0.* g has characteristic polynomial (x³+1)²⁴ on ℝ⁷², which has no eigenvalue 1. Λ has no norm-2 vectors.
  No twisted sector reaches weight 1.
* *T₃/₂ = 113.* No g^j with j ≠ 0 fixes a norm-3 vector: L has minimum norm 6, and g³ = −1. So Z₆ acts freely on the
  678 norm-3 vectors, and 678/6 = 113 = 24 + 89.
* *T₂ = 2 928.* The Heisenberg part, Sym²(ℝ⁷²)·α₋₁α₋₁ ⊕ ℝ⁷²·α₋₂, contributes
  (1/6)[2 700 + 2 556] = 876. Only j = 0 and j = 3 contribute: tr g³ = −72 gives (5 184+72)/2 − 72 = 2 556. The
  12 312 norm-4 vectors are permuted freely and give 2 052. Total 876 + 2 052 = 2 928, as computed.

**The first twisted fields [computed].** In the g²-twisted sector at L₀ = 3, the computed twining traces are
Tr ĝ² = 56 and Tr ĝ = Tr ĝ³ = −24. On the ground state |0⟩ (L₀ = 8/3), ĝ acts as ω, so the ground state is
projected out. On the weight-3 states, ĝ² = 1 and ĝ³ is −1 on the 24 oscillator states and swaps |w⟩ ↔ |−w⟩. So
exactly the 16 combinations |w⟩ + |−w⟩ survive. The same happens in the g⁴-twisted sector. These 32 fields, all
bosonic, are the lowest twist fields of T.

### 3.2 Ramond ground states of T **[computed]**

| sector | R ground states (L₀ = 3) | Z₆-invariant | even / odd |
|---|---|---|---|
| untwisted | 83 000 (full-support words) | (83 000 + 80 + 80)/6 = 13 860 | 6 928 / 6 932 |
| g²-twisted | 80 (minimal vectors of L* + χ_L/2, from the 80 full-support words of C̄) | 40 | 20 / 20 |
| g⁴-twisted | 80 | 40 | 20 / 20 |
| g, g³, g⁵ | none (ground energies 19/6 and 9/2 > 3) | 0 | — |
| **total** | | **13 940** | **6 968 / 6 972** (net −4 ✓) |

How the parities are fixed:
* Untwisted: the traces of ĝ^j on the 83 000 ground states are 0 for odd j (nothing is fixed). For j = 2, 4 the trace
  is 80 (the fixed words; A1), with supertrace 40 − 40 = 0. Hence the invariant part has supertrace −24/6 = −4.
* g²-twisted: Z^{(2,b)}_{PP} ≡ 0 for all b, and Z^{(2,odd)}_{PA} = 0 (these are images of the vanishing R-twining
  function of g). So the 40 invariant states split 20/20.

The Witten index −4 = −24/6 agrees with NOTES §2.1. Here every nontrivial pair contributes 0, because the balanced
fixed words give str ĝ² = 0.

The full R character of T is 13 940 q⁰ + 2 579 935 232 q¹ + 12 668 385 656 832 q² + … (powers of q^{L₀−3}).
The super-character Σ_a Z_{PP} is the constant −4.

## 4. T′ = Sym²T (c = 144) and P = Sym²(T′ ⊗ Arf) (c = 288)

### 4.1 Formulas **[proved, standard; NOTES §1.5, §2.2]**

Let X be a holomorphic N=1 SVOA with central charge c_X.
* **NS sector.** The swap-twisted NS sector of Sym²X is the Ramond sector of X on a doubled circle, with
  L₀ = L₀^R/2 + c_X/16 ≥ c_X/12. This bound is 6 for X = T and 12 for X = T′. Below it, Sym²X is the
  graded-symmetric part of X ⊗ X. With the Koszul sign, Tr(σq^{L₀}) = Tr((−1)^Fσq^{L₀}) = χ^{AP}_X(q²), so
  χ_{Sym²X} = ½[χ_X(q)² + χ^{AP}_X(q²)] and χ^{AP}_{Sym²X} = ½[χ^{AP}_X(q)² + χ^{AP}_X(q²)].
  Arf ⊗ (·) does not change the NS sector, and the Kitaev sign does not affect untwisted NS functions.
* **Ramond ground states.** Write B/F for the even/odd ground states of X.
  * Untwisted part: graded-symmetric, so even = S²B ⊕ Λ²F and odd = B ⊗ F.
  * Swap-twisted part: one copy of the ground states of X with the same parities. Here Z^{PA}(σ,1) = χ^R_X(q^{1/2})
    and Z^{PA}(σ,σ) = its T-image, which equals +1 on the ground states. Also Z^{PP}(σ,1) = Z^{PP}(σ,σ) = I(X).

  The index is then (I² + I)/2 + I = (I² + 3I)/2. X ⊗ Arf exchanges B and F.

### 4.2 Results **[computed]** (`fc_sym2.py`)

| h | T | T′ | P |
|---|---|---|---|
| 0 | 1 | 1 | 1 |
| 1/2, 1 | 0 | 0 | 0 |
| 3/2 | 113 | 113 | 113 |
| 2 | 2 928 | 2 928 | 2 928 |
| 5/2 | 42 001 | 42 001 | 42 001 |
| 3 | 962 444 | 968 772 = 962 444 + Λ²(113) | 975 100 = 968 772 + Λ²(113) |
| 7/2 | 31 069 242 | 31 400 106 | 31 730 970 |
| 4 | 1 879 408 646 | 1 888 442 815 | 1 897 476 984 |
| 9/2 | 136 784 011 042 | 137 015 746 142 | 137 248 196 306 |
| 5 | 6 484 810 801 152 | 6 492 021 682 530 | 6 499 288 479 924 |

In words: up to weight 5/2, T′ and P have exactly the fields of a single copy of T, symmetrised over the copies.
The first "multi-copy" fields are the bilinears ψ⁽ⁱ⁾ψ⁽ʲ⁾ of two weight-3/2 fields at weight 3, giving
Λ²(113) = 6 328 new fields at each step.

| Ramond ground states | untwisted (even / odd) | twisted (even / odd) | total (even / odd) | index |
|---|---|---|---|---|
| T′ | 48 580 902 / 48 580 896 | 6 968 / 6 972 | **48 587 870 / 48 587 868** | 2 ✓ |
| P | 2 360 781 013 961 161 / 2 360 781 013 961 160 | 48 587 868 / 48 587 870 | **2 360 781 062 549 029 / 2 360 781 062 549 030** | −1 ✓ |

The indices agree with (I²+3I)/2: −4 ↦ 2, and then (via Arf, I = −2) ↦ −1. As a single polynomial in I = I(T),
P has index (I⁴ + 6I³ + 3I² − 18I)/8. This equals the form (I⁴ − 6I³ + 3I² + 18I)/8 used in the original statement of the question, evaluated at −I(T) = +4;
NOTES §2.3 is written for the index +4 seed U. Both give −1.

**Answer to "is the Ramond vacuum literally unique?"** No **[computed]**. P has 4 721 562 125 098 059 Ramond ground
states. The Witten index −1 is the net count (one more odd state than even). The NS vacuum *is* unique **[proved]**.
The 113 weight-3/2 fields contain the supercurrent and no other distinguished field; the weight-1 space is zero.

## 5. Limitations

* A1 and A2 are the standard inputs of the orbifold construction, not re-derived from scratch here. The momentum part
  of A1 is proved in §2 (given GJF Thm 3.1); its Ramond part is GJF's treatment of odd-order lifts (proof of their
  Thm 5.6). A2 is the standard modular-covariant definition of the twisted-sector Z₆ action of the anomaly-free orbifold
  (GJF §4, Thm 5.6); its spin refinement changes no character (§2). The integrality and spin–statistics checks of §2
  all pass.
* The orbifold characters are exact through L₀ = 5 (c = 72). The weights 4–5 are included only as consistency data.
* An explicit operator basis (actual vertex operators with cocycles) was not written out. The counts and the
  orbit descriptions above determine it only up to the choice of basis within each Z₆-orbit or D₄-orbit.
* Corrections to the original statement of the question: g³ = −1 has Frame shape 2⁷²1⁻⁷² (not 1⁻¹²2¹²). The g-twisted ground energy is 19/6 with
  multiplicity 2¹², the g²-twisted ground energy is 8/3 with multiplicity 1, and the g³-twisted ground energy is 9/2
  with multiplicity 2³⁶.

## 6. Closed forms for the NS characters (`fc_closed_form.py` → `fc_closed_form.json`)

### 6.1 Structure **[proved, given modular invariance]**

Let X be a holomorphic N=1 SVOA with 24 | c (as for T, T′, P: c = 72, 144, 288) and Ramond weights ≥ c/24 (unitarity,
G₀² = L₀ − c/24).
Write Z_X(τ) = Tr_NS q^{L₀−c/24}. Then Z_X is invariant under Γ_θ = ⟨S, T²⟩: S-invariance is the modular
invariance of the NS–NS torus **[cited: standard, e.g. Höhn's thesis]**, and T² acts trivially because 2·c/24 ∈ ℤ.
At the cusp ∞, Z_X has a pole of order c/12 in q^{1/2}. At the other cusp (τ = 1), Z_X turns into the Ramond
character Z^R_X, which starts at q⁰, so Z_X is regular there. Γ_θ has genus 0 with Hauptmodul

  **t(τ) = (ϑ₃/η)¹² = η(τ)⁴⁸ / (η(τ/2) η(2τ))²⁴ = q^{−1/2} + 24 + 276 q^{1/2} + 2048 q + 11202 q^{3/2} + 49152 q² + …**

This is the NS character of 24 free fermions; it has a simple pole at ∞ and a zero at the cusp 1. Therefore

  **Z_X = Σ_{j=0}^{c/12} a_j t^j,  with a_{c/12} = 1.**

The same coefficients a_j give the other two spin structures:
* **Z̃_X(τ) = Z_X(τ+1) = Σ_j (−1)^j a_j (ϑ₄/η)^{12j}.**
* **Z^R_X(τ) = Σ_j (−1)^j a_j (ϑ₂/η)^{12j} = Σ_j a_j (−4096)^j (η(2τ)/η(τ))^{24j}.**

In particular, **a₀ is the number of Ramond ground states.** The Witten index is not visible in Z_X; it needs the
Z^{PP} sector. (With 12 | c only, the vacuum phase e^{−2πi c/24} gives an extra sign (−1)^{c/12} in the first two
formulas, and the Ramond character is (−1)^{c/12} Σ_j (−1)^j a_j (ϑ₂/η)^{12j}, with ground-state count (−1)^{c/12} a₀; for
c = 12, P(t) = t − 24 and the Ramond character is u + 24. For 24 | c the formulas are as stated.)

There is an equivalent form that shows the "simplification". Let t₀ = t − 24, and let F_k be the Faber polynomials:
monic of degree k in t₀, with F_k(t₀) = q^{−k/2} + O(q^{1/2}) and no constant term. Then

  **Z_X = Σ_{0 ≤ h < c/24} (dim X_h) · F_{c/12 − 2h}(t₀) + dim X_{c/24}.**

So the whole NS character is fixed by the dimensions of the weight spaces X_h with h ≤ c/24. The first Faber
polynomials are:
* F₁ = t₀
* F₂ = t₀² − 552
* F₃ = t₀³ − 828t₀ − 6144
* F₄ = t₀⁴ − 1104t₀² − 8192t₀ + 107544
* F₅ = t₀⁵ − 1380t₀³ − 10240t₀² + 324870t₀ + 2580480
* F₆ = t₀⁶ − 1656t₀⁴ − 12288t₀³ + 618372t₀² + 6488064t₀ − 12019872

### 6.2 T (c = 72) **[computed]**

  **Z_T = t⁶ − 144 t⁵ + 6 984 t⁴ − 129 679 t³ + 751 404 t² − 629 867 t + 13 940**

  **  = F₆(t₀) + 113 F₃(t₀) + 2 928 F₂(t₀) + 42 001 F₁(t₀) + 962 444.**

The polynomial was fitted from the weights h ≤ 3 alone. It then reproduces the independently computed orbifold
dimensions at h = 7/2, 4, 9/2 and 5, and the whole Ramond character
13 940 + 2 579 935 232 q + 12 668 385 656 832 q² + … from §3.2. Its constant term a₀ = 13 940 is the
Ramond ground-state count. The coefficient −144 of t⁵ just records dim T₁/₂ = 0.

### 6.3 T′ and P from Z_T **[computed; formula proved as in §4.1]**

The swap orbifold has a closed form in terms of the characters of X. With I = I(X) and ε = +1 for Sym²X,
ε = −1 for Sym²(X⊗Arf):

  **Z_{Sym²X}(τ) = ½ [ Z_X(τ)² + Z̃_X(2τ) + Z^R_X(τ/2) + ε I ],**
  **Z^R_{Sym²X}(τ) = ½ [ Z^R_X(τ)² + ε I + Z^R_X(τ/2) + Z^R_X((τ+1)/2) ].**

The four terms come from the pairs (1,1), (1,σ), (σ,1) and (σ,σ). The constant ε I is the (σ,σ) term,
Z^{AP}(σ,1) = I, shifted by T. The Kitaev sign falls only on it.

Composing the two steps gives Z_P explicitly in terms of t:
* T′: ε = +1 with X = T, I = −4.
* P: ε = −1 with X = T′, I = 2.

In both cases the result is again a polynomial in t, which checks the formula; the script checks it through q^{12}
(q^{15} for T′).
* **Z_{T′}** = t¹² − 288t¹¹ + 34 704t¹⁰ − 2 270 863t⁹ + 87 648 096t⁸ − 2 030 562 659t⁷ + 27 545 157 180t⁶
  − 204 469 234 212t⁵ + 731 620 354 115t⁴ − 973 544 209 008t³ + 365 122 692 738t² − 11 873 213 572t
  + **97 175 738**. The constant term equals 48 587 870 + 48 587 868 Ramond ground states ✓.
* **Z_P** = Σ_{j=0}^{24} a_j t^j. The coefficients are listed in `fc_closed_form.json`. They run from a₂₄ = 1,
  a₂₃ = −576, a₂₂ = 152 352 down to a₀ = **4 721 562 125 098 059**, which is the Ramond ground-state count of §4.2 ✓.

In Faber form:

  **Z_P = F₂₄ + 113 F₂₁ + 2 928 F₂₀ + 42 001 F₁₉ + 975 100 F₁₈ + 31 730 970 F₁₇ + … + dim P₁₂**

The weights 0 to 12 contribute as follows. The dimensions from h = 6 on include the twisted sectors of T′ (h ≥ 6)
and of P itself (h ≥ 12):

| h | dim P_h | h | dim P_h |
|---|---|---|---|
| 0 | 1 | 6½ | 89 638 685 414 344 486 |
| 3/2 | 113 | 7 | 1 412 443 362 317 858 730 |
| 2 | 2 928 | 7½ | 20 062 277 606 109 789 280 |
| 5/2 | 42 001 | 8 | 280 295 415 365 778 999 446 |
| 3 | 975 100 | 8½ | 4 621 082 593 072 333 903 898 |
| 7/2 | 31 730 970 | 9 | 118 843 932 996 730 698 362 774 |
| 4 | 1 897 476 984 | 9½ | 4 523 772 280 273 612 426 910 310 |
| 9/2 | 137 248 196 306 | 10 | 181 863 088 276 262 876 538 953 233 |
| 5 | 6 499 288 479 924 | 10½ | 6 508 374 198 522 608 847 025 630 113 |
| 11/2 | 206 120 857 647 702 | 11 | 200 002 421 246 142 519 615 740 165 006 |
| 6 | 4 818 700 471 201 380 | 11½ | 5 303 360 359 410 139 763 326 259 501 900 |
| | | 12 | 123 101 670 625 668 498 528 898 165 575 027 |

(dim P₁/₂ = dim P₁ = 0.)

**The most compact statement.** Everything is generated by the seven integers of Z_T and two Sym² steps:

  Z_T = t⁶ − 144t⁵ + 6984t⁴ − 129679t³ + 751404t² − 629867t + 13940,  t = (ϑ₃/η)¹²,
  Z_{T′} = ½[Z_T² + Z̃_T(2τ) + Z^R_T(τ/2) − 4],
  Z_P = ½[Z_{T′}² + Z̃_{T′}(2τ) + Z^R_{T′}(τ/2) − 2],

where Z̃(τ) = Z(τ+1) and Z^R = the same polynomial evaluated at −(ϑ₂/η)¹². The fitted polynomial for T is
exact as an identity of q-series to the computed order (L₀ ≤ 5). Given §6.1, it is exact for all orders, since a
polynomial of degree 6 is fixed by its 7 leading coefficients.

### 6.4 Ramond partition function (bosons and fermions counted together) **[computed; derivation proved]**

Let Z^R_X(τ) = Tr_R q^{L₀−c/24}, with no (−1)^F inserted. It is the image of the NS character under T and then S:

  Z_X(τ) --T--> Z̃_X(τ) = Z_X(τ+1) = Tr_NS (−1)^F q^{L₀−c/24} --S--> Z^R_X(τ) = Z̃_X(−1/τ) = Z_X(1 − 1/τ).

For the Hauptmodul, t(τ+1) = −(ϑ₄/η)¹² because η(τ+1)¹² = −η(τ)¹², and ϑ₄(−1/τ)/η(−1/τ) = ϑ₂(τ)/η(τ). Hence

  **t(1 − 1/τ) = −u(τ),  u := (ϑ₂/η)¹² = 4096 η(2τ)²⁴/η(τ)²⁴ = 4096 q ∏_{n≥1}(1+qⁿ)²⁴ = 4096q + 98 304q² + 1 228 800q³ + …**

u is the Ramond character of 24 free fermions. So **if Z_X = P_X(t), then Z^R_X = P_X(−u)**. The coefficients
b_j = (−1)^j a_j turn out to be **all positive**:

* **Z^R_T** = 13 940 + 629 867 u + 751 404 u² + 129 679 u³ + 6 984 u⁴ + 144 u⁵ + u⁶
  = 13 940 + 2 579 935 232 q + 12 668 385 656 832 q² + 9 517 357 430 358 016 q³ + …
* **Z^R_{T′}** = 97 175 738 + 11 873 213 572 u + 365 122 692 738 u² + 973 544 209 008 u³ + 731 620 354 115 u⁴
  + 204 469 234 212 u⁵ + 27 545 157 180 u⁶ + 2 030 562 659 u⁷ + 87 648 096 u⁸ + 2 270 863 u⁹ + 34 704 u¹⁰ + 288 u¹¹ + u¹²
* **Z^R_P** = Σ_{j=0}^{24} b_j u^j, where b₀ = 4 721 562 125 098 059, b₁ = 1 155 284 118 797 296 712, …,
  b₂₂ = 152 352, b₂₃ = 576, b₂₄ = 1. The full list is `R_poly_in_u` in `fc_closed_form.json`.

The constant term is the number of Ramond ground states, and u^j starts at q^j, so the level-n coefficient of
Z^R only involves b₀,…,b_n. The top coefficient b_{c/12} = 1 is the image of the NS vacuum.

**Bosons and fermions separately.** Z^{PP}_X is the constant I(X), so

  **Z^R_{X,even} = ½(Z^R_X + I),  Z^R_{X,odd} = ½(Z^R_X − I).**

Above the ground level, bosons and fermions are paired by G₀. The first levels are:

| level q^n | T (even = odd for n ≥ 1) | T′ | P |
|---|---|---|---|
| q⁰ | 6 968 / 6 972 | 48 587 870 / 48 587 868 | 2 360 781 062 549 029 / 2 360 781 062 549 030 |
| q¹ | 1 289 967 616 each | 24 316 341 395 456 each | 2 366 021 875 296 863 666 176 each |
| q² | 6 334 192 828 416 each | 3 063 454 733 477 019 648 each | 994 356 398 634 044 869 265 719 296 each |
| q³ | 4 758 678 715 179 008 each | 33 597 749 420 880 416 833 536 each | 231 913 873 738 540 900 846 571 492 671 488 each |

Checks:
* *T.* P_T(−u) reproduces the Ramond character computed directly from the orbifold (§3.2).
* *T′ and P.* P(−u) agrees with the independent Sym² Ramond formula of §6.3,
  ½[Z^R_X² + εI + Z^R_X(τ/2) + Z^R_X((τ+1)/2)], through q^{12} (`fc_closed_form.py` asserts both).
* *Numerical modular check (`fc_modular_check.py`, mpmath, 60 digits).* t(1 − 1/τ) = −(ϑ₂/η)¹² holds at three test points to relative
  error ~10⁻⁶⁰. At τ = 0.1+1.6i and −0.3+2.2i, the transformed NS polynomial Σ a_j t(1−1/τ)^j agrees with the
  q-series of Z^R for T, T′ and P to between 10⁻³⁹ and 10⁻⁶¹. The worst case, 10⁻³⁹, is P at Im τ = 1.6, where it
  is limited by the truncation of the q-series at q²⁵.
