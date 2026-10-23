# Is the monomial Q8 action on the length-96, index-24 code theory non-anomalous?

Notes of 17 October 2026 (Claude), extended step by step; prepared for publication on 20 October (remarks about
internal workflow removed or reworded; the mathematics, the corrections and the data are unchanged).
Companion files in this folder: the certificate `../certificates/k08_claude_q8.json` (original:
`../certificates/original/q8_certificate.json`) and NOTES_k8.md (§4.2, §5, §6).

## Goal

V = lattice SVOA of Construction-A of a self-dual ternary code C of length 96
(c = 96, k = 8), Index(V) = ±24. Q8 acts monomially on C (24 blocks of 4
coordinates (1,i,j,k), left multiplication by unit quaternions as signed
permutations). Decide whether the induced Q8 action on V is non-anomalous,
i.e. whether V//Q8 exists as a holomorphic SVOA with Index 24/8 = 3, the TMF
minimum at c = 96. Every cyclic subgroup is already known to be non-anomalous
(NOTES_k8.md §6); the open bit is the part of H^3(Q8,U(1)) = Z8 invisible to
cyclic subgroups.

## Status table

| Step | Content | Status |
|------|---------|--------|
| S1 | Lift of Q8 ⊂ 2^96:S_96 to Aut(V); action of −I | done (proved/cited) |
| S2 | H^3(Q8,U(1)), restrictions, 2105.02231, Ω_3^Spin(BQ8) exponent | done (cited; exponent 8) |
| S3 | 24-identical-blocks argument (Narain/genus form) | done: α = 0, conditional on items (A)–(E) listed in S3.4 |
| S4 | Direct computation for the c=4 block B | partial: not needed for the verdict; heuristic block anomaly recorded; stepwise-orbifold check outlined |
| S5 | Verdict and §9 of NOTES_k8.md | done (superseded by S6 if S6 closes) |
| S6 | VOA-level proof without the Narain assumption (equivariant neighbour chain to Z^96) | done: α = 0 with no physics assumption; inputs listed in S6.4 |

## S1. Lift

**Independent check of the certificate (exact, 17 October).** The check script is not included; the same group
relations are checked by `../code/verify.py` (for `../certificates/k08_claude_q8.json`: i and j preserve the code,
i² = j² = (ij)² = −I, ji = −(ij), the group generated has order 8, commuting pairs are cyclic), and
`../code/neighbour_chain.py` checks that i and j are compatible with the frame $L^{24}(1+i+j)$ built on the 4-blocks
($g f_a=\pm f_{p(a)}$); the block-by-block form of the action is not re-checked by an included script. The original
check reads `symmetry_i`, `symmetry_j` (convention e_a → signs[a]·e_{perm[a]}, signs ∈ {1,2} ⊂ F_3) and
checks as signed permutation matrices: i² = j² = (ij)² = −I, ji = −(ij), i⁴ = I, the group generated
has order exactly 8, every element preserves the 24 blocks {4b,…,4b+3}, and on block 0 i acts as
e_1→e_i, e_i→−e_1, e_j→e_k, e_k→−e_j and j as e_1→e_j, e_i→−e_k, e_j→−e_1, e_k→e_i, i.e. left
multiplication by the unit quaternions in the basis (1,i,j,k). Also re-verified: the generator matrix
has rank 48, G·Gᵀ ≡ 0 mod 3 (self-dual), and i, j map every generator row into the code. So
Q8 = ⟨i,j⟩ ⊂ 2^96:S_96 with centre ⟨−I⟩ is an honest quaternion group of monomial code
automorphisms. **Status: proved (computation).**

**How a monomial automorphism acts on V = V_{Λ(C)} (GJF Theorem 3.1, p. 9).** The signed
permutation group 2^c:S_c acts on V_{√3Z^c} = V_{√3Z}^{⊗c}: S_c permutes tensor factors, and each
reflection acts on its factor by the order-2 automorphism Γ_{+√3} ↔ Γ_{−√3} of V_{√3Z} (the lift of
n√3 ↦ −n√3), which fixes the N=1 generator G = Γ_{+√3} + Γ_{−√3}. This is a genuine group action
(no cocycle), so the extension (√3Z^c)^[2].(2^c:S_c) ⊂ Aut(V_{√3Z^c}) splits. Because √3Z^c has odd
index 3^{c/2} in Λ(C), the restriction map Λ̂(C)[2] → (√3Z^c)^[2] is an isomorphism, hence the
extension Λ̂(C)[2].Aut(C) ⊂ Aut(V_{Λ(C)}) splits too, and the splitting is the unique one compatible
with the explicit action on the subalgebra V_{√3Z^c}. (GJF's own words: "the extension
Λ̂(C)[2].Aut(C) splits", and the lift "manifestly preserves the N=1 structure".)

Consequences for us:
* Each g ∈ Aut(C) has a canonical lift ĝ ∈ Aut_{N=1}(V), and g ↦ ĝ is a **group homomorphism**
  (it is the pull-back of a splitting). Therefore ⟨î, ĵ⟩ ⊂ Aut_{N=1}(V) is a quaternion group of order
  exactly 8, with î² = ĵ² = (îĵ)² = (−I)^ and îĵ = k̂. **No double cover arises.** (A double cover would
  need a non-split extension by signs; Theorem 3.1 rules it out.) **Status: cited (GJF Thm 3.1) +
  verified relations.**
* The action is faithful (it is faithful on the weight-1 Heisenberg currents h_a, a = 1..96, which
  transform in the 96-dimensional signed-permutation representation). Q8 is solvable, so V^{Q8}
  is regular [CM16] and GJF Theorem 4.2 applies: there is a well-defined anomaly α ∈ S̃H³(BQ8),
  and each trivialisation gives an orbifold V//Q8. **Status: cited.**
* GJF's Example 4.4 ("Majorana layer vanishes iff ⟨χ,(1−g)χ⟩ ∈ 4Z", χ = √3(1,…,1)) applied here:
  for g = −I, ⟨χ,(1−g)χ⟩ = 2χ² = 6c = 576 ∈ 4Z; for g = i, j, k (48 negative 2-cycles each, so
  ⟨χ,gχ⟩ = 3·(48−48) = 0), ⟨χ,(1−g)χ⟩ = χ² = 288 ∈ 4Z. So the Majorana layer α⁽¹⁾ ∈ H¹(Q8;Z2)
  vanishes: every element of Q8 acts **evenly** on the Ramond sector. **Status: proved
  (GJF Ex. 4.4 + arithmetic).**

**How −I acts; it is a bosonic Z2, not (−1)^F.** GJF treat g = −1 ∈ Aut(C) exactly as an
automorphism of the SVOA (Theorem 5.6 with m = 2: "If m is even, then g^{m/2} = −1 ∈ Aut(C) is the
central element"; Example 5.7(2), p. 20: "We set g = −1 ∈ Aut(C)" for the length-24, index-24 code,
giving index 12). In the proof of Theorem 5.6 they lift g "via Theorem 3.1", check the Majorana
layer ⟨χ,(1−g)χ⟩ ∈ 4Z, and compute its anomaly from the Frame shape 1^{−12k}2^{12k} as a lattice
automorphism λ ↦ −λ (Example 4.4; the same computation appears for OddLeech, p. 16, "the Frame
shape of g is 1^{−24}2^{24}" and "the action of g on V_{OddLeech} is nonanomalous"). Concretely, on
V_{√3Z}^{⊗96} the lift of −I is ⊗ of the reflection Γ_{+√3} ↔ Γ_{−√3}; it sends every Heisenberg
current h_a to −h_a, whereas (−1)^F fixes all currents. Hence (−I)^ ≠ (−1)^F, and (−I)^ is an
even (bosonic) order-2 automorphism of the NS sector commuting with G. It is also distinct from
(−1)^F·(−I)^. **Status: cited (GJF Thm 5.6 proof, Ex. 5.7(2), Ex. 4.4) + direct observation.**

**Remark on lift choices.** Other lifts of Q8 (twist by a character Q8 → C* ⊂ Aut_{N=1}(V), or by
(−1)^F on a quotient) are possible in principle, but the index argument in the certificate
(`orbifold_argument`) and NOTES_k8 §4.2 use *this* lift, so this is the one whose anomaly matters.


## S2. Bookkeeping: anomaly groups

### S2.1 Bosonic anomalies: H³(Q8,U(1)) ≅ Z8 and its restrictions (standard; verified by hand)

Q8 ⊂ SU(2) acts freely on S³, so it has periodic cohomology of period 4 and
H³(Q8;U(1)) ≅ H⁴(Q8;Z) ≅ Z/|Q8| = Z8, generated by e = c₂(τ), the second Chern class (Euler class)
of the 2-dimensional irreducible representation τ: Q8 → SU(2) (Cartan–Eilenberg, Ch. XII;
Yu 2105.02231 eq. (13)–(14) also lists H^i(BG;U(1)) for finite G ⊂ S³).

Restrictions, computed from characters:
* To Z4 = ⟨i⟩: τ|_{Z4} = χ ⊕ χ⁻¹ with χ a faithful character, so e|_{Z4} = c₁(χ)c₁(χ⁻¹) = −x²,
  where x = c₁(χ) generates H²(Z4;Z) = Z4 and x² generates H⁴(Z4;Z) = Z4. Hence
  res: Z8 → Z4 is **surjective**, kernel {0,4}. Same for ⟨j⟩, ⟨k⟩ (all conjugate under the outer
  automorphism permuting i, j, k).
* To the centre Z2 = ⟨−1⟩: τ|_{Z2} = sgn ⊕ sgn, e|_{Z2} = c₁(sgn)² = y² ≠ 0 ∈ H⁴(Z2;Z) = Z2, so res is
  reduction mod 2, kernel {0,2,4,6}.
* **Joint kernel of all restrictions to cyclic subgroups: {0, 4} ⊂ Z8**, i.e. exactly the class 4e.
  **Status: proved (hand computation above).**

### S2.2 What arXiv:2105.02231 (M. Yu, "Genus-one data and anomaly detection") proves

* **Proposition 1.1** (p. 2): "The genus-one data applied to detect anomalies for the symmetry given
  by the dicyclic group of order 4N, Dic_N, has an undetectable Z2 kernel." Q8 = Dic_2 (N = 2).
* **Proposition 3.3** (p. 7): "Let G be a subgroup of SU(2), then G has no joint CSD" (CSD = the
  cohomology classes detectable by commuting pairs / genus-one data). Proof via S³/G → BG.
* **Proposition 3.4** (p. 9): for Dic_N with anomaly a 4N-th root of unity, the most informative
  genus-one data is (g = rotation generator, h = reflection), and even that cannot distinguish the
  elements of the Z2 kernel (p. 9–10: "we are unable to distinguish the elements in the kernel").
* §4 works the example SU(2)_k WZW with Q8 ⊂ SU(2)_L: the Z8 anomaly k mod 8 is only seen mod 4.

Genus-one data = all twisted–twined torus partition functions Z(g,h), [g,h] = 1, with their
SL(2,Z) action. For Q8 every commuting pair lies in a cyclic subgroup, so genus-one data is
exactly the data of the cyclic subgroups; Proposition 1.1 says precisely that the class 4e
of S2.1 is invisible to it. **Consequence: no computation of (twisted, twined) characters of V —
NS or R, with or without (−1)^F — can decide whether α ∈ {0, 4e}.** This rules out the
"Frame shape / eta-product multiplier" method of GJF Example 4.4 and the modular-data method of
Gemünden–Keller 1909.09626 for the remaining bit (GK say so themselves, p. 5: "For the cases we
consider however, it turns out that ω is determined by its restriction to various cyclic
subgroups"). **Status: cited.**

### S2.3 The correct anomaly group for a fermionic theory with bosonic Q8 symmetry

V is a holomorphic SVOA (fermionic) and Q8 acts by even automorphisms commuting with (−1)^F, so
the full symmetry type is Spin × Q8 (G_f = Z2^f × Q8). GJF Theorem 4.2: the anomaly lives in
extended supercohomology S̃H³_gp(Q8), filtered by the layers
H³(Q8;U(1)) = Z8 (bosonic) ⊂ … → H²(Q8;Z2) = Z2² (Gu–Wen) → H¹(Q8;Z2) = Z2² (Majorana),
so |SH³(BQ8)| ≤ 8·4·4 = 128. (H*(BQ8;Z2) = Z2[x,y,w]/(x²+xy+y², x²y+xy²), |x|=|y|=1, |w|=4; DGL
eq. (B.1).) By Gaiotto–Johnson-Freyd 1712.07950 (and Freed–Hopkins), SH is the Anderson dual of
spin bordism truncated to the relevant homotopy, and in degree 3 for a finite group this gives

    S̃H³(BQ8) ≅ Hom(Ω̃₃^Spin(BQ8), U(1)).

(The AHSS for the Anderson dual in degree 3 only involves π₀ = U(1), π₋₁ = Z2, π₋₂ = Z2, the same
three layers, since Ω₃^Spin(pt) = 0.) **Status: cited.**

### S2.4 Ω̃₃^Spin(BQ8) and its exponent

Three independent sources agree:

| Source | Statement |
|---|---|
| Bayen–Bruner, *Real connective K-theory and the quaternion group*, Trans. AMS 348 (1996) 2201–2216 | ko_*(BQ8) via the stable splitting BQ8 ≃ BSL(2,F3) ∨ Σ⁻¹BS³/BN ∨ Σ⁻¹BS³/BN (2-locally); degree 3: Z/8 ⊕ Z/4 ⊕ Z/4 (quoted as (A.9)–(A.10) with n=0 in DGL) |
| Botvinnik–Gilkey, *An analytic computation of ko_{4ν−1}(BQ8)*, Topol. Methods Nonlinear Anal. 6 (1995) 127–135 | same groups by η-invariants |
| Barrera-Yanez–Gilkey, math/0205084, Thm 5.1 | k̃o_{4k−1}(BQ_ℓ) = Z_{2^{k+1}} ⊕ Z_{2^{k+1}} ⊕ B_{k+1} for k odd, of order 4^{k+1}ℓ^k; k = 1, ℓ = 8: Z4 ⊕ Z4 ⊕ B₂ with |B₂| = 128/16 = 8 |
| Davighi–Gripaios–Lohitsiri, 2207.10700 (JHEP 09 (2022) 147), Table 6 and the table on p. 9 | Ω̃₃^Spin(BQ8) = Z/8 × (Z/4)² (also Ω̃₁ = Ω̃₂ = (Z/2)², Ω̃₄ = Z/2, Ω̃₅ = (Z/2)², Ω̃₆ = Z/2) |
| Fraiman–Larotonda–Tartaglia, 2607.15349, §4.1 eq. (4.17)–(4.26) | Ω₃^Spin(BQ8) ≃ Z8 ⊕ Z4^{⊕2}, "strictly bigger than" Ω₃^Spin(BZ4) ≃ Z8 ⊕ Z2; order-8 generator X_{Q8} = S³/τ(Q8) with η(ρ_i,X_{Q8}) = η(ρ_j,·) = η(ρ_k,·) = 1/2, η(ρ_τ,X_{Q8}) = 7/8; the Z4 ⊕ Z4 is generated by any two of X_i = S³/τ(⟨i⟩), X_j, X_k, with η(ρ_τ, X_g) = 3/4 and η(ρ_g, X_h) = (1−δ_{gh})/2 |

(Ω₃^Spin(pt) = 0, so reduced = unreduced in degree 3; and Ω̃₃^Spin ≅ k̃o₃ since bordism and ko agree
below degree 8.) The order is 128, saturating the layer bound of S2.3, so all three layers are
fully realised and the identification S̃H³(BQ8) ≅ Hom(Z8 ⊕ Z4², U(1)) ≅ **Z8 ⊕ Z4 ⊕ Z4** holds.

**Exponent of the anomaly group: 8.** In particular **24·α = 0 for every α ∈ SH³(BQ8)**; in fact
already 8·α = 0. **Status: cited (five sources above).**

Remarks. (i) The η-invariant table shows the relation X_i + X_j + X_k = −2X_{Q8} (check: η(ρ_τ) = 9/4 ≡
1/4 ≡ −7/4 mod 2; η(ρ_i) = 1 ≡ −1 mod 2), so the images of the cyclic subgroups span an index-2
subgroup of Ω₃^Spin(BQ8); dually the anomalies invisible to all cyclic subgroups form exactly one Z2,
which must be the bosonic class 4e of S2.1 (consistent with Yu's Proposition 1.1). (ii) For our
action the Majorana layer vanishes (S1), so α ∈ rSH³(BQ8), of order 32 and exponent 8.


## S3. The "24 identical blocks" argument

### S3.1 Why the argument cannot be run on W = B^{⊗24} as literally stated

B = V_{√3Z⁴} (c = 4) and W = B^{⊗24} = V_{√3Z^{96}} are lattice SVOAs of the **non-unimodular**
lattice √3Z^{96} (determinant 3^{96}), so W is not holomorphic and V = V_{Λ(C)} is a simple-current
extension of W of index 3^{48} (by the code C ≅ Λ(C)/√3Z^{96}). For a non-holomorphic W the
G-twisted modules form a G-crossed braided extension of the (super-)MTC Rep(W) ≅ pointed category
of F_3^{96}; this extension always exists (W^G is regular), and its "associator" is not a class in
SH³(BG). So "the anomaly of G on W" is not defined in the same group as α(V), and statement (a) of
the question as first posed is not well-posed as written. One can define the anomaly of V from the pair (twisted
W-modules, the G-stable algebra V ∈ Rep(W)), but then it depends on the code C and additivity over
blocks is not automatic. **Status: the literal version of (a) is withdrawn; replaced by S3.2–S3.4.**

### S3.2 The correct form: anomalies of lattice theories depend only on the genus of the G-lattice

The physical content of GJF's phrase "the 't Hooft anomaly of the action of g is 24 times that of an
automorphism of Frame shape (m/2)^{−k/m}m^{k/m}" (proof of Thm 5.6) is that the anomaly depends only
on the Q8-representation Λ ⊗ Q, not on the integral lattice. For cyclic groups this follows from
the eta-product formula (genus-one data); for Q8 it does not, by S2.2. The following replaces it.

**Setup.** Let Λ = Λ(C) (odd unimodular, rank 96) with its Q8 action, and let Γ := Λ ⊕ (−Z^{96})
with Q8 acting diagonally (on Z^{96} by the same signed permutations). Γ is an odd unimodular lattice
of signature (96,96). The Narain CFT T(Γ, P) of Γ at the decomposition P: Γ⊗R = Λ⊗R ⊕ (−Z^{96})⊗R
is the non-chiral fermionic CFT V_Λ ⊗ \overline{V_{Z^{96}}}, and its Q8 anomaly is

    α(T) = α(V_Λ, Q8) − α(V_{Z^{96}}, Q8) ∈ Hom(Ω₃^Spin(BQ8), U(1))

(antiholomorphic sector contributes with orientation reversed, i.e. inverted).

**Claim 1 (lattice statement).** There is a Q8-equivariant isometry
Λ(C) ⊕ (−Z^{96}) ≅ Z^{96} ⊕ (−Z^{96}).

**Claim 2 (connectedness).** The locus of Narain moduli at which a fixed Q8 ⊂ O(Γ) is a symmetry is
connected.

**Claim 3 (physics).** The anomaly of a symmetry is constant along a connected family of CFTs
preserving it.

Given Claims 1–3: at the far end of the family the CFT is V_{Z^{96}} ⊗ \overline{V_{Z^{96}}} with
Q8 acting by the same lift on both sides, whose anomaly is α(V_{Z^{96}}) − α(V_{Z^{96}}) = 0. Hence
**α(V_Λ, Q8) = α(V_{Z^{96}}, Q8).** The right-hand side is a free-fermion anomaly: V_{Z^{96}} = Fer(192)
= Fer(8)^{⊗24} and the Q8 lift is block-diagonal, so by additivity of anomalies under tensor products
with the diagonal action (S3.3) α(V_{Z^{96}}) = 24·α(Fer(8), Q8), and 24·(anything) = 0 in a group of
exponent 8 (S2.4). **Therefore α(V_Λ, Q8) = 0.**

### S3.3 Proofs / justifications of the claims

**Claim 1.** Write L = Z⟨1,i,j,k⟩ ⊂ H_Q for the Lipschitz order; Z^{96} = L^{24} with Q8 acting by
left multiplication, and Λ(C) = (1/√3)M with M = {x ∈ Z^{96} : x mod 3 ∈ C}, 3Z^{96} ⊂ M ⊂ Z^{96},
[Z^{96}:M] = 3^{48}, M a left L-submodule (C is Q8-invariant, S1). A Q8-invariant Z-bilinear form b
on a left L-module determines the L-hermitian form h(x,y) = b(x,y) + b(x,iy)i + b(x,jy)j + b(x,ky)k
(h(lx,y) = l h(x,y), h(x,ly) = h(x,y) l̄, h(y,x) = \overline{h(x,y)}, Re h = b); on L^{24} with the
Euclidean form this is h(x,y) = Σ_a x_a ȳ_a. Because L is self-dual for Re(x ȳ), the h-dual of a
lattice equals its b-dual, so (Λ(C), h) is a **unimodular hermitian L-lattice**, as is (Z^{96}, h).
Equivariant isometries of Q8-lattices = isometries of hermitian L-lattices.

*Local isometry at every prime p* (Λ(C)_p ≅ Z_p^{96}):
* p ≠ 3: M_p = Z_p^{96} = L_p^{24}, so Λ(C)_p ≅ (L_p^{24}, (1/3)·std). Put u = (1−i−j)/3 ∈ L_p (3 ∈ Z_p^×);
  N(u) = 3/9 = 1/3 ∈ Z_p^×, so u ∈ L_p^× and right multiplication R_u: x ↦ xu is a left-L_p-linear
  bijection with h(xu, yu) = h(x,y) u ū = (1/3) h(x,y). So R_u: (L_p^{24}, std) → (L_p^{24}, (1/3)std) is a
  Q8-equivariant isometry. (This includes p = 2, where H_Q is ramified and L_2 is non-maximal: the
  only thing used is that u is a unit of L_2.)
* p = 3: H_Q splits at 3 and L has 2-power discriminant, so L_3 ≅ M_2(Z_3) with the canonical
  involution (= adjugate, of symplectic type). Both Λ(C)_3 and Z_3^{96} are unimodular hermitian
  L_3-lattices of L_3-rank 24. Morita equivalence (Knus, *Quadratic and Hermitian Forms over Rings*,
  Ch. I §9) takes unimodular hermitian lattices over (M_2(Z_3), adj) to unimodular **alternating**
  Z_3-lattices of rank 48, and unimodular alternating lattices over a PID are unique up to isometry
  (symplectic basis). Hence Λ(C)_3 ≅ Z_3^{96} as hermitian L_3-lattices.
* p = ∞: both Λ(C)⊗R and Z^{96}⊗R are positive definite of H-rank 24.

So Λ(C) and Z^{96} lie in the same **genus** of hermitian L-lattices, and therefore so do the
indefinite lattices Γ = Λ(C) ⊕ (−Z^{96}) and Γ₀ = Z^{96} ⊕ (−Z^{96}) (hermitian signature (24,24)).

*Genus ⇒ class.* The unitary group U(h_Γ) of a hermitian form over a quaternion algebra with its
canonical involution is an algebraic group of type C (over R it is Sp(24,24)), **simply connected**
and **non-compact at ∞**, so Kneser's strong approximation theorem applies: U(Q) is dense in U(A_f),
hence the class set U(Q)\U(A_f)/∏_p Stab(Γ_p) of lattices in the genus of Γ is a single point,
i.e. **Γ ≅ Γ₀ as hermitian L-lattices**, which is Claim 1. (Kneser, *Starke Approximation in
algebraischen Gruppen I*, J. reine angew. Math. 218 (1965); the same argument gives the classical
fact that indefinite unimodular lattices are unique in their genus.) **Status: proved here modulo
two cited theorems (Morita equivalence for hermitian forms; strong approximation).**

*Side remark (needed for Claim 3's lift bookkeeping).* All lifts of the Q8-action from Γ to the
Narain CFT are conjugate: the lifts of a split extension form a torsor under H¹(Q8, Γ̂),
Γ̂ = Hom(Γ,U(1)), and 0 → Γ* → Hom(Γ,R) → Γ̂ → 0 with Hom(Γ,R) uniquely divisible gives
H¹(Q8, Γ̂) ≅ H²(Q8, Γ*) ≅ H²(Q8, Γ). Now Γ ⊗ Z₂ ≅ 48·L₂ (M₂ = Z₂^{96}) and
L = Z[Q8] ⊗_{Z[⟨−1⟩]} Z^{sgn} is induced from the centre, so by Shapiro's lemma
H²(Q8, L) ≅ H²(Z2, Z^{sgn}) = Ĥ⁰(Z2, Z^{sgn}) = 0; H²(Q8,Γ) is a 2-group, so H²(Q8,Γ) = 0.
The same computation with Λ in place of Γ shows that **the GJF lift of Q8 to V_{Λ(C)} is the unique
lift up to conjugacy** (strengthening S1). **Status: proved.**

**Claim 2.** Γ⊗R ≅ 48·H as a real Q8-representation (H = the 4-dimensional irreducible, quaternionic
type), so the centralizer of Q8 in O(96,96) is the quaternionic unitary group U(24,24;H) = Sp(24,24),
and the set of Q8-stable decompositions Γ⊗R = V₊ ⊕ V₋ is the Riemannian symmetric space
Sp(24,24)/(Sp(24)×Sp(24)), which is connected. By Claim 1 the two points (Λ⊗R, −Z^{96}⊗R) and
(Z^{96}⊗R, −Z^{96}⊗R) (transported by the isometry) lie in this one fixed locus. The Q8 lift to the
CFT is defined uniformly along the locus (it is a cocycle on Γ, independent of the modulus), and by
the side remark any two lifts are conjugate anyway. **Status: proved (standard).**

**Claim 3.** The anomaly of a finite symmetry of a 2d fermionic QFT is the class of a 3d invertible
spin-TQFT, a deformation invariant; continuous families of CFTs preserving the symmetry have constant
anomaly. (For cyclic groups this is visible on the genus-one data, which varies continuously in
Narain moduli but takes values in roots of unity.) This is standard physics but **not a theorem about
VOAs**; it is the one genuinely non-mathematical input. **Status: assumed (physics-standard).**

**Additivity under tensor products.** For holomorphic V₁, V₂ with G acting diagonally,
Rep_g(V₁⊗V₂) = Rep_g(V₁) ⊠ Rep_g(V₂) and the G-graded super-fusion categories multiply, so the
associator classes add: α(V₁⊗V₂) = α(V₁) + α(V₂) in SH³(BG). GJF use exactly this ("24 times")
in the proof of Thm 5.6. For free fermions this is the linearity of the η-invariant in the
representation (Grigoletto 2112.01485 §4, the map ch₂: RO(G) → Hom(Ω₃^Spin(BG),U(1)); Fraiman et al.
§4.1). **Status: cited/standard.**

### S3.4 Conclusion of S3

Under (A) Claim 3 [assumed, physics], (B) strong approximation for Sp(24,24) [cited theorem],
(C) the local isometries of S3.3 [proved here], (D) Ω̃₃^Spin(BQ8) ≅ Z8 ⊕ Z4² of exponent 8 and
SH³ ≅ its Pontryagin dual [cited, S2.3–S2.4], (E) additivity under tensor products [standard]:

    α(V_{Λ(C)}, Q8) = α(Fer(192), Q8) = 24·α(Fer(8), Q8) = 0.

**The monomial Q8 action on V_{Λ(C)} is non-anomalous** for every Q8-invariant self-dual code of
length 96 with the block structure of the certificate — and in fact for every length c = 4·(number of
blocks) with the number of blocks divisible by 8 (c ∈ 32Z), since already 8·α(Fer(8)) = 0. In
particular the same holds for the length-192 codes (48 blocks, k = 16).

### S3.5 Sanity checks against known cyclic results

* Restricting to Z4 = ⟨i⟩: GJF's Frame-shape formula for one block B (shape 2^{−2}4²) gives
  Σ_k d_k (m/k)(k²−1) = −2·2·3 + 2·1·15 = 18, i.e. eigenvalue exp(2πi·18/24) = exp(2πi·6/8): the
  single-block Z4 anomaly is 6 ∈ Z8 = rSH³(Z4), of order 4; for 24 blocks, 24·6 ≡ 0 mod 8 ✓ (and
  Σ = 432 ≡ 0 mod 24, GJF's own check). For 4 blocks the Z4 anomaly already vanishes (24 ≡ 0 mod 8).
* Restricting to the centre Z2 = ⟨−1⟩: Frame shape 1^{−96}2^{96}, Σ = −96·2·0 + 96·1·3 = 288 ≡ 0 mod 24 ✓.
* The free-fermion model Fer(8) with the same signed-permutation lift has the same Frame shapes and
  the same Majorana-layer values (⟨χ,(1−g)χ⟩ = 4 ∈ 4Z for g = i, 8 for g = −1), hence the same cyclic
  anomalies — consistent with Claim 1 (which says more: the full Q8 anomalies agree).
* The argument predicts nothing inconsistent: on cyclic subgroups it reproduces GJF's "24 times"
  statement, and it only adds the invisible class 4e, which it kills by 8 | 24.


## S4. Direct computation on the block B

S3 closes the question (conditionally on Claim 3), so S4 is a cross-check, not a load-bearing step.

### S4.1 What a direct computation on B = V_{√3Z⁴} could and could not do

* B is not holomorphic (det √3Z⁴ = 81), so "ω_B ∈ SH³(BQ8)" is not defined (S3.1); the natural
  holomorphic c = 4 block with the same Q8 action is Fer(8) = V_{Z⁴} (Lipschitz quaternions, Q8 =
  unit quaternions by left multiplication), and by S3 the two give the same Q8 anomaly once 24
  copies are glued into a unimodular lattice.
* Gemünden–Keller's non-abelian machinery (1909.09626) decides triviality of the 3-cocycle ω from
  modular data of V^G, hence only "when ω is determined by its restriction to various cyclic
  subgroups" (their p. 5) — exactly what fails for the class 4e. **Not usable for the invisible bit.**
* The stepwise route proposed in the original question needs a holomorphic theory on which ⟨i⟩ is non-anomalous: the
  single block has Z4 anomaly 6 ∈ Z8 (S3.5), so the smallest candidate is 4 blocks, Fer(32) = V_{Z^{16}}
  (Z4 anomaly 24 ≡ 0). Then Fer(32)//⟨i⟩ is a holomorphic c = 16 SVOA with residual symmetry
  K = ⟨Ẑ4 (quantum symmetry), j⟩, an extension 1 → Z4 → K → Z2 → 1 with j acting on Z4 by inversion
  (j i j⁻¹ = i⁻¹), so K ≅ D8 or Q8 according to whether j² acts on the twisted sectors as the trivial
  or the order-2 element of Ẑ4. The Q8 anomaly of Fer(32) (predicted below to be exactly 4e ≠ 0)
  would reappear as the K-anomaly of Fer(32)//⟨i⟩, detectable by genus-one data if K ≅ D8. This is
  a feasible exact computation (c = 16 free fermions, 2¹⁶-type lattice bookkeeping) but was **not
  executed** here (bounded effort; not needed for the verdict). **Status: open (optional check).**

### S4.2 Heuristic value of the single-block anomaly from η-invariants (convention-dependent)

On Fer(8) the GJF-type lift of the signed permutations acts on the eight Majorana–Weyl fermions
ψ_a = Γ_{e_a} + Γ_{−e_a} (fixed by reflections) and ψ̄_a ∝ Γ_{e_a} − Γ_{−e_a} (sign-flipped) through the
real representation ρ_block = (permutation rep of Q8 on 4 points) ⊕ H = 1 ⊕ ρ_i ⊕ ρ_j ⊕ ρ_k ⊕ H,
H ⊗ C = 2τ. If the anomaly of MW fermions in a real representation ρ evaluated on (X, spin,
Q8-bundle) is exp(iπ·η_Dirac(ρ⊗C, X)) with the reduced η-invariants of Fraiman et al. (4.24)–(4.26):

| generator X | ½·[η(ρ_i)+η(ρ_j)+η(ρ_k)+2η(ρ_τ)] mod 1 | order |
|---|---|---|
| X_{Q8} = S³/Q8 | ½·(3/2 + 7/4) = 13/8 ≡ 5/8 | 8 |
| X_i = S³/⟨i⟩ | ½·(0 + 1/2 + 1/2 + 3/2) = 5/4 ≡ 1/4 | 4 |

So in this normalisation α(Fer(8), Q8) has order exactly 8, its restriction to each Z4 has order 4
(matching 6 ∈ Z8 from S3.5), and 4·α(Fer(8)) is a non-zero class killed by every cyclic subgroup —
necessarily the invisible class 4e. Predictions (heuristic): Fer(32) with 4 blocks is Q8-anomalous
with anomaly exactly 4e (a genus-one-invisible anomaly!), Fer(64) with 8 blocks is not, and
V_{Λ(C)} with 24 blocks is not. **Status: numerical/heuristic (depends on the η normalisation; the
order-8 statement is robust to a factor of 2 only in one direction).** Nothing in S3 depends on it.


## S5. Verdict

### Status table (what rests on what)

| Item | Statement | Status |
|---|---|---|
| S1 | Q8 ⊂ 2^96:S_96 is an honest Q8 (relations, order 8, blocks, code invariance) | **proved** (exact computation) |
| S1 | It lifts to an honest Q8 ⊂ Aut_{N=1}(V) (no double cover); the lift is unique up to conjugacy | **cited** (GJF Thm 3.1) + **proved** (H²(Q8,Λ)=0, S3.3) |
| S1 | −I acts as a bosonic (even) Z2, ≠ (−1)^F; Majorana layer of every g ∈ Q8 vanishes | **cited** (GJF Thm 5.6, Ex. 4.4, Ex. 5.7(2)) + arithmetic |
| S2 | H³(Q8,U(1)) = Z8; joint kernel of restrictions to cyclic subgroups = {0,4e} | **proved** (hand) |
| S2 | Genus-one data (all twisted–twined characters) cannot see 4e | **cited** (Yu 2105.02231 Prop. 1.1, 3.3, 3.4) |
| S2 | Anomaly group = Hom(Ω̃₃^Spin(BQ8),U(1)) ≅ Z8 ⊕ Z4 ⊕ Z4, **exponent 8** | **cited** (Bayen–Bruner 96; Botvinnik–Gilkey 95; Barrera-Yanez–Gilkey Thm 5.1; DGL 2207.10700 Table 6; Fraiman et al. (4.17)) |
| S3 | Λ(C) ⊕ (−Z^96) ≅ Z^96 ⊕ (−Z^96) as Q8-lattices (local isometries + strong approximation) | **proved here modulo cited theorems** (Morita; Kneser) |
| S3 | Q8-fixed locus in Narain moduli is connected (Sp(24,24)/Sp(24)×Sp(24)) | **proved** (standard) |
| S3 | Anomaly constant along the connected symmetric family | **assumed** (physics-standard; not a VOA theorem) |
| S3 | Anomalies add under ⊗ with diagonal action; 24·α(Fer(8)) = 0 | **standard** + exponent 8 |
| S3 | **α(V_{Λ(C)}, Q8) = 0** | follows from the above |
| S4 | Single-block anomaly has order 8; 4 blocks anomalous (=4e), 8 blocks not | **heuristic** (η normalisation) |
| S4 | Independent VOA-level check via Fer(32)//⟨i⟩ and its K ∈ {D8,Q8} symmetry | **open** (optional) |

### Verdict (one paragraph) — as of 17 Oct; superseded by the S6 update below

The monomial Q8 action on the index-24, length-96 code theory is **non-anomalous, conditionally**:
the only non-mathematical input is the deformation invariance of the 't Hooft anomaly along the
connected Q8-symmetric Narain family joining V_{Λ(C)} ⊗ \overline{V_{Z^96}} to V_{Z^96} ⊗ \overline{V_{Z^96}}
(Claim 3 of S3), everything else being either an exact computation done here, a lattice argument
(same genus + strong approximation for Sp(24,24)), or a cited theorem (Ω̃₃^Spin(BQ8) ≅ Z8 ⊕ Z4² of
exponent 8, SH³ = Pontryagin dual of spin bordism, GJF Thms 3.1/4.2/5.6). The argument is the
rigorous form of the "24 identical blocks" intuition: it reduces the anomaly to 24 copies of a
free-fermion block in a group of exponent 8, and in fact shows that **any** number of blocks divisible
by 8 gives a non-anomalous Q8. Consequently V_{Λ(C)}//Q8 exists (for each of the SH²(BQ8)-many
trivialisations), is a holomorphic N=1 SVOA at c = 96, and by the vanishing of all twisted–twined RR
traces (every commuting pair of Q8 is cyclic and every non-identity element has −1 as a power, GJF
Thm 5.6 mechanism; discrete torsion only multiplies vanishing traces by phases) has
**Index = ±24/8 = ±3, the TMF minimum at k = 8.** The same applies verbatim to the two verified
length-192 Q8-codes (48 blocks), giving **index 3 at k = 16**, and via tensor products with the
known index-1 theories to k = 17, 19, 20, 22. What remains for a VOA-level proof: either a theorem
that the GJF anomaly of a lattice (S)VOA with a lifted finite G ⊂ O(Λ) depends only on the genus of
the G-lattice, or a direct exact computation of the invisible class (e.g. the Fer(32)//⟨i⟩ check of
S4.1, or the associator of the i-, j-, k-twisted modules). Nothing in S1–S4 points towards an
anomaly; the single-block heuristic (S4.2) even explains why GJF's "24" is exactly what is needed.



### Verdict update (18 Oct, after S6)

The physics assumption (Claim 3 of S3) is no longer needed. S6 gives a VOA-level proof:
Z^96 = Λ(tetracode^{⊕24}) with the same monomial Q8 (exact), the certificate code is joined to
tetracode^{⊕24} by 23 Q8-equivariant neighbour steps (exact), each step is an orbifold by a
Q8-stable Z_3^2 inside GJF's torus C* whose Q8 anomaly is preserved (gauging lemma S6.3, standard
categorical argument + coprimality), and α(Fer(192), Q8) = 24·α(Fer(8), Q8) = 0 by exponent 8.
**Status: Index 3 at c = 96 (and at c = 192 for the length-192 Q8 codes) is established modulo
cited theorems only** (GJF Thms 3.1/4.2 with the standard extension↔algebra dictionary for twisted
modules; additivity under ⊗; Ω̃₃^Spin(BQ8) ≅ Z8 ⊕ Z4²). The one step written as a sketch rather
than a full proof is the gauging lemma S6.3(d); it is standard tensor-category material and does not
depend on any computation specific to this code.

## S6. Towards a VOA-level proof (2026-10-18)

Goal: remove Claim 3 of S3 (deformation invariance of the anomaly along a Narain family), the only
non-mathematical input. Plan:

1. Z^96 is itself a Construction-A lattice: with u = 1+i+j ∈ L (N(u) = 3), right multiplication
   R_u: L^24 → L^24 is Q8-equivariant (commutes with left multiplication) and scales norms by 3, so
   F := L^24·u ≅ √3Z^96 is a Q8-stable 3-frame inside Z^96 and Z^96 = Λ(C') for a Q8-invariant
   self-dual code C' ⊂ F_3^96 with the SAME block Q8-module as the certificate code C. (Expect
   C' = tetracode^{⊕24}: per block, Z^4/L·u has order 9, and the unique self-dual ternary code of
   length 4 is the tetracode, Aut = GL(2,3) ⊃ Q8.)  → exact check below.
2. Connect C to C' by Q8-equivariant neighbour steps C_{i+1} = (C_i ∩ U^⊥) + U with U ⊂ C' an
   irreducible Q8-submodule (a block tetracode) not contained in C_i. Each step is self-dual and
   Q8-invariant; dim(C_i ∩ C') grows by ≥ 2 per step, so ≤ 24 steps.  → explicit chain below.
3. Each step is an orbifold of V_{Λ(C_i)} by E = (C_i + U)/C_i ≅ Z_3^2 ⊂ C_i^* = F_3^96/C_i, i.e. by a
   Q8-stable 3-subgroup of GJF's torus C* ⊂ Aut_{N=1}(V). Show α_{Q8}(V//E) = α_{Q8}(V): anomaly of
   E ⋊ Q8 is inflated from Q8 (coprimality + existence of the E-orbifold), and gauging a normal
   subgroup with inflated anomaly preserves the quotient anomaly (categorical lemma).
4. α(V_{Z^96}) = α(Fer(192)) = 24·α(Fer(8)) = 0 (tensor additivity + exponent 8).


### S6.1 Z^96 is the Construction-A lattice of tetracode^{⊕24} (proved, exact computation)

Script: `neighbour_chain.py` (output `neighbour_chain_result.json`). In block coordinates (1,i,j,k),
right multiplication by u = 1+i+j gives the frame vectors
f_1 = (1,1,1,0), f_i = (−1,1,0,1), f_j = (−1,0,1,−1), f_k = (0,−1,1,1):
norms 3, pairwise orthogonal, and **g·f_a = s_a f_{p(a)} for the certificate's i and j** (the frame
transforms by the same signed permutations as the standard basis, because R_u commutes with left
multiplication). So F = L^24·u ≅ √3Z^96 is a Q8-stable 3-frame in Z^96 and
Z^96 = Λ(C') with C' := {(⟨x,f_a⟩ mod 3)_a : x ∈ Z^96} ⊂ F_3^96.
Verified: C' is self-dual of dimension 48, Q8-invariant under the certificate's signed permutations,
block-diagonal, and each block is the **tetracode** (2-dimensional, all 8 nonzero words of weight 3;
Aut(tetracode) = GL(2,3) ⊃ Q8). Hence **V_{Z^96} = Fer(192) = V_{Λ(tetracode^{⊕24})}, with the
monomial Q8 of the certificate acting on it exactly as on the certificate code theory.** Also
dim(C ∩ C') = 2 and C contains no weight-3 word.

### S6.2 The Q8-equivariant neighbour chain (proved, exact computation)

For a Q8-stable totally isotropic U ⊂ F_3^96 and a Q8-invariant self-dual C_i, the code
C_{i+1} := (C_i ∩ U^⊥) + U is Q8-invariant and self-dual: (C_i ∩ U^⊥)^⊥ = C_i + U, so
C_{i+1}^⊥ = (C_i + U) ∩ U^⊥ = (C_i ∩ U^⊥) + U by the modular law (U ⊂ U^⊥). Taking U = U_b, the
tetracode of block b (an irreducible Q8-module ≅ σ, the unique 2-dimensional irreducible of
F_3[Q8]), one has U_b ∩ C_i = 0 unless U_b ⊂ C_i (irreducibility), and
dim(C_{i+1} ∩ C') ≥ dim((C_i ∩ C') + U_b) = dim(C_i ∩ C') + 2. Since C' = ⊕_b U_b, the chain stops
only at C'. Computed chain for the certificate code: blocks b = 0,…,22 each give a valid step
(U_b ∩ C_i = 0, new code self-dual and Q8-invariant), dim(C_i ∩ C') = 2,4,6,…,48, and block 23 is
already contained; **23 steps from C to C' = tetracode^{⊕24}.**

### S6.3 Each step preserves the Q8 anomaly (VOA/categorical argument)

Fix one step V := V_{Λ(C_i)}, U := U_b, V' := V_{Λ(C_{i+1})}.

(a) **The step is an orbifold by a 3-group inside GJF's torus.** E := (C_i + U)/C_i ≅ U ≅ Z_3^2 acts on
V through C_i^* = F_3^96/C_i ⊂ Aut_{N=1}(V) (GJF Thm 3.1: Γ_λ ↦ ζ^{χ(λ)}Γ_λ for a character χ of
Λ(C_i)/√3Z^96 = C_i). V^E = V_{Λ(C_i ∩ U^⊥)} ⊂ V' (as C_i ∩ U^⊥ ⊂ C_{i+1}), and V' ∩ V = V^E
(C_i ∩ C_{i+1} = C_i ∩ U^⊥). The GJF lift Q̂8 normalises E (q̂ e q̂⁻¹ is the character χ∘q⁻¹, and U is
Q8-stable) and Q̂8 ∩ E = 1, so H := E ⋊ Q8 ⊂ Aut(V) is an honest group of order 72. **Proved.**

(b) **E is non-anomalous and V' = V//_ψ E.** By GJF Thm 4.2 the E-twisted modules form an E-graded
super-fusion category F_E = ⊕_{e∈E} Rep_e(V) with each component ≃ SVec and associator
α_E ∈ S̃H³(BE) = H³(E;U(1)) (E has odd order, so the Z2-layers vanish and E acts evenly). V' is a
holomorphic extension of V^E whose intersection with V is V^E, i.e. it contains exactly one simple
module from each e-twisted sector with the e = 1 component V^E only; its vertex operators make
⊕_e 1_e an associative algebra in F_E, which exists iff α_E = 0, and then the algebra structures are
the twisted group algebras C^ψ[E], ψ ∈ H²(E;U(1)) ≅ Z3. So α_E = 0 and V' = V//_ψ E for some ψ.
(Equivalently, Lagrangian algebras of Rep(V^E) ≅ Z(Vec_E^{α_E}) are the pairs (K ≤ E, ψ) with
α_E|_K = 0, and only K = E has untwisted component V^E alone.) **Status: standard dictionary
(Kirillov; GJF Thm 4.2 and its proof); applied here.**

(c) **The anomaly of H on V is inflated from Q8.** S̃H³(BH) → S̃H³(BQ8) ⊕ H³(E;U(1))^{Q8},
x ↦ (res_{Q8}x, res_E x), is injective: H^1(H;Z2) ≅ H^1(Q8;Z2) and H²(H;Z2) ≅ H²(Q8;Z2) (E has odd
order, Lyndon–Hochschild–Serre), so an x with res_{Q8}x = 0 has vanishing Majorana and Gu–Wen
layers and comes from y ∈ H³(H;U(1)); H³(H;U(1)) ≅ H³(Q8;U(1)) ⊕ H³(E;U(1))^{Q8} via the two
restrictions (LHS with H^p(Q8; H^q(E;U(1))) = 0 for p,q ≥ 1 by coprimality, split since the pieces
have coprime orders), and H³ ↪ rSH³ in degree 3 (GJF p. 11), so y = 0. *[18 Oct: GJF state this inclusion explicitly just before Thm 4.2, because the connecting maps are stable operations vanishing in low degrees, so the step is valid as written. An alternative direct proof, a p-local comparison showing that (res_{Q8}, res_E) is an isomorphism, is in `gauging_lemma.md` §G12.1 and tex Prop. 2.3.]* The anomaly α̃ of H has
res_{Q8} α̃ = α := α(V,Q8) and res_E α̃ = α_E = 0 (by (b)), as does inf(α) := π*α for
π: H → Q8; hence **α̃ = inf(α)**. **Proved.**

(d) **Gauging lemma.** Let H = E ⋊ Q act on a holomorphic SVOA V with |E| coprime to |Q|, E
non-anomalous, anomaly α̃ = inf(α) with α ∈ S̃H³(BQ), and ψ ∈ H²(E;U(1)) Q-invariant. Then Q acts on
V//_ψ E, and this action has anomaly α. Sketch: by GJF Thm 4.2, F := ⊕_{h∈H} Rep_h(V) is an H-graded
extension of SVec with pointed components; choose cocycle data (σ, ω) for α on Q and inflate: F ≃
SVec_H^{inf(σ,ω)}, whose associator and Gu–Wen signs are literally independent of the E-coordinate.
Choose a Q-invariant cocycle representative of ψ (possible: the class is Q-invariant and |Q| is
prime to |E|, U(1)-cochains being divisible; here Q8 ⊂ SL(2,3) acts trivially on H²(Z_3²;U(1)) ≅ det).
Then A := C^ψ[E] = ⊕_{e∈E} 1_e is an algebra in F_E ≅ Vec_E on which Q acts by algebra
automorphisms (conjugation). The objects A⊗1_q (q ∈ Q) are A-bimodules; since the cocycle is
inflated, the bimodule maps 1_e ⊗ 1_q → 1_{eq} carry no ω- or σ-factors, (A⊗1_{q1}) ⊗_A (A⊗1_{q2}) ≅
A⊗1_{q1q2} with the parity sign σ(q1,q2) and the associator ω(q1,q2,q3). By the standard dictionary
(twisted modules of the extension V//_ψE ⊂-restricted to V are the A-modules in the coset components
F_{qE}, and A⊗1_q is the unique simple one up to odd isomorphism), ⊕_q Rep_q(V//_ψE) ≃ SVec_Q^{(σ,ω)}:
the anomaly of Q on V//_ψE is α. The complement Q ⊂ H is unique up to E-conjugacy (Schur–Zassenhaus),
and conjugate complements give isomorphic bimodules (1_e ∈ A), so the result does not depend on it.
**Status: standard categorical argument, sketched here; the only VOA input is the dictionary
"twisted modules of a simple-current extension = modules of the corresponding algebra object",
the same one underlying GJF Thm 4.2.**
*Full proof (18 Oct): `gauging_lemma.md`. It also corrects "unique simple one": there are |E| simple bimodules per coset, and the extra ones are the dual Ê sectors.*

(e) **The induced Q8 on V' is the GJF lift.** The complement Q8 acts on V' = V_{Λ(C_{i+1})} as a
group of automorphisms acting on the Heisenberg currents by the signed-permutation representation,
hence as a lift of the monomial Q8 ⊂ O(Λ(C_{i+1})). Any two homomorphic lifts differ by a cocycle
in Z¹(Q8, Λ̂(C_{i+1})), and H¹(Q8, Λ̂) ≅ H²(Q8, Λ) = 0 for every lattice in this family
(Λ ⊗ Z_2 = Z_2^96 = 24·L_2 and Shapiro, S3.3), so the induced action is conjugate to the GJF lift of
Theorem 3.1, and conjugate actions have the same anomaly. **Proved.**

Combining (a)–(e): **α(V_{Λ(C_i)}, Q8) = α(V_{Λ(C_{i+1})}, Q8) for every step of the chain.**

### S6.4 Conclusion: VOA-level proof

α(V_{Λ(C)}, Q8) = α(V_{Λ(C')}, Q8) = α(Fer(192), Q8) = 24·α(Fer(8), Q8) = 0,

using S6.1 (Λ(C') = Z^96 with the same Q8; V_{Z^96} = V_Z^{⊗96} = Fer(8)^{⊗24} with the blockwise lift,
conjugate to any other lift by S6.3(e)), S6.2 (23 steps), S6.3 (each step), additivity of anomalies
under tensor products with the diagonal action, and exponent 8 of S̃H³(BQ8) ≅ Hom(Ω₃^Spin(BQ8),U(1))
≅ Z8 ⊕ Z4 ⊕ Z4 (S2.4). **No deformation/Narain assumption is used.** The inputs now are:
(i) exact computations [S1, S6.1, S6.2]; (ii) GJF Theorems 3.1 and 4.2 and the standard
"extension ↔ algebra object" dictionary for twisted modules [cited/standard]; (iii) additivity under
⊗ [standard]; (iv) the spin-bordism group of BQ8 [cited, three independent computations];
(v) elementary group cohomology (coprimality splitting, Shapiro) [proved here].

Remarks. (1) The chain never leaves the N=1 world: each E ⊂ C_i^* ⊂ Aut_{N=1}, so all V_{Λ(C_i)} are
N=1 SCFTs with Q8 ⊂ Aut_{N=1}; this is not needed for the anomaly. (2) The same argument works for
every Q8-invariant self-dual code with this block structure, at every length 4m: α = m·α(Fer(8)),
so the Q8 orbifold is non-anomalous iff 8 | m (granted the heuristic order-8 value of S4.2; the
"if" direction is a theorem). (3) The argument is the chiral counterpart of the Narain argument of S3:
"same genus" is replaced by "connected by equivariant neighbours", and deformation invariance by the
gauging lemma. (4) It also shows directly that two Q8-codes with the same block module have equal
anomalies, so the anomaly of the certificate code equals that of any of the other 41 (or 294) codes.
