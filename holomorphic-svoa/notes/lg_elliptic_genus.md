# The elliptic genus of the LG orbifold (Σ Φᵢ³)/C and the lattice double

Supplementary note to the paper's §6.3–6.4 and Appendix A (Claim 1, the gapping test, refined elliptic genera, and the
comparison of all four spin-structure partition functions). Status labels: **proved**, **cited**, **computed**,
**assumed**. Checks: `../code/lg/` (numpy and the standard library; times below are on one core). Written by GPT-6 Astra,
Fable 5.1 and Opus 5.5 (AI agents), supervised by Du Pei; CC BY 4.0.

What is proved here are identities between the functions written down (sections 2–4) and between the four torus
partition functions of two specified orbifold data (section 5). The physical statement of Claim 1, that the LG orbifold
flows to the lattice double, additionally uses the cited LG/minimal-model correspondence and the standard assumption that
gauging a symmetry preserved along an RG flow commutes with the flow (section 6). Equality of characters or partition
functions is evidence for, not a proof of, an identification of the full CFTs including their operator products.

## 1. Setup and conventions

- C is a linear ternary code of length n; from section 4 on it is self-dual (so 4 | n). For w in F₃ⁿ with all wᵢ ≠ 0,
  write sign(w) = (−1)^#{i : wᵢ = 2} and ŵ ∈ {+1, −1}ⁿ for the lift 1 → +1, 2 → −1. Index(C) = Σ over full-support
  w ∈ C of sign(w) (the paper's §2, with − = 2). c ≤ w means cᵢ ∈ {0, wᵢ} for every i; wt(c) is the number of nonzero
  entries.
- q = e^{2πiτ}, y = e^{2πiz}, ω = e^{2πi/3};
  θ₁(τ, u) = −i q^{1/8} e^{iπu} ∏_{n≥1} (1 − qⁿ)(1 − e^{2πiu} qⁿ)(1 − e^{−2πiu} q^{n−1}),
  with θ₁(τ, u + aτ + b) = (−1)^{a+b} q^{−a²/2} e^{−2πiau} θ₁(τ, u) for integers a, b. Throughout, a and b denote fixed
  integer representatives in {0, 1, 2}.
- **R-charge convention.** Φ has the left-moving R-charge 1/3 (vector R-charge 2/3, the usual normalisation of W = Φ³ of
  charge 2); the fugacity y couples to the left-moving U(1) current J₀. In this section the lattice SVOA V and its
  conjugate are placed as in the paper's Appendix A: the U(1)_L refinement acts on V̄ = V̄_Λ(C), the other factor gives
  the Witten index. (Elsewhere in the paper and in these notes the holomorphic theory is written as V; the chirality
  labels are a convention.)
- **Lattice side.** V = V_Λ(C), Λ = { x/√3 : x ∈ ℤⁿ, x mod 3 ∈ C }, with the frame N = 2 structure (one copy of the
  ℓ = 1 N = 2 algebra per coordinate: J = i∂X/√3, G± = e^{±i√3X}, c = 1). Its Ramond sector is Λ + δ, δ = (3/2)(1, …, 1)/√3
  (2δ is characteristic and lies in Λ). For one coordinate and a label u ∈ ℤ₃ (u = 0, 1, 2 for 0, +, −) put

      ψ_u(τ, z) = η(τ)⁻¹ Σ_{x ∈ ℤ, x ≡ u mod 3} (−1)^x y^{(x + 3/2)/3} q^{(x + 3/2)²/6},

  and ι_w = ψ_w(τ, 0) (exactly a constant): ι₀ = 0, ι₁ = 1, ι₂ = −1 (the per-coordinate Witten index of a label).
  The per-coordinate fermion-parity factor is e^{iπ(x + 3/2)} = (−1)^x (−i); the constants (−i) multiply to (−i)ⁿ = 1
  (4 | n) and are dropped. Then

      χ^R_V(τ, z) = Σ_{u ∈ C} ∏ᵢ ψ_{uᵢ}(τ, z)

  is the Ramond-sector **supertrace** Tr_R (−1)^F q^{L₀ − n/24} y^{J₀} (in the rest of these notes "Ramond character"
  means the unsigned trace), and the elliptic genus of the double (fugacity on one side, Witten index on the other) is
  Index(C) χ^R_V(τ, z).
- **LG side.** For one field with W = Φ³ and the ℤ₃ phase twisted by ωᵃ in space and ωᵇ in time,

      F_{a,b}(τ, z) = y^{−a/3} θ₁(τ, −2z/3 + (aτ + b)/3) / θ₁(τ, z/3 + (aτ + b)/3),

  and G_{a,k} = (1/3) Σ_{b ∈ ℤ₃} ω^{−bk} F_{a,b}, so that F_{a,b} = Σ_k ω^{bk} G_{a,k}. The prefactor y^{−a/3} is part of
  the formula (without it the identities below fail). The elliptic genus of the orbifold (Σᵢ Φᵢ³)/C, C acting by
  Φᵢ ↦ ω^{cᵢ} Φᵢ, without discrete torsion, is

      Z_{LG/C}(τ, z) = (1/|C|) Σ_{c, c′ ∈ C} ∏ᵢ F_{cᵢ, c′ᵢ}(τ, z).

- **Where F_{a,b} comes from.** The untwisted genus F_{0,0} is Witten's (Int. J. Mod. Phys. A 9 (1994) 4783,
  hep-th/9304026, §2.2.2, an infinite product without flavor fugacity; the θ₁ form is the same product; **cited**). The
  time twists F_{0,b} follow from the flavor action Φ ↦ ωᵇ Φ, which multiplies each field by e^{2πibr} with its
  R-charge r, i.e. shifts z ↦ z + b. The space twists follow by spectral flow z ↦ z + aτ, the orbifold formula of
  Kawai, Yamada, Yang (Nucl. Phys. B 414 (1994) 191, hep-th/9306096; **cited**), applied per field, with prefactor
  q^{a²/6} y^{a/3} for ĉ = 1/3. The overall phase per sector is not fixed by that formula for fractional ĉ; it is fixed
  here (**proved** by direct computation): the displayed F_{a,b} does not depend on the integer representatives of a, b;
  F_{a,b}(τ, 0) = −2 for (a, b) = (0, 0) (a limit of a 0/0 expression) and 1 otherwise, the signed fixed-point counts of the
  Morse function of section 7; and it is covariant without extra phases,

      F_{a,b}(τ + 1, z) = F_{a, a+b}(τ, z),        F_{a,b}(−1/τ, z/τ) = e^{πiz²/(3τ)} F_{b, −a}(τ, z).

  Given F_{0,b}, covariance fixes every sector (c, c′) with c, c′ linearly dependent. The remaining freedom is a discrete
  torsion, an alternating bicharacter on C × C; "no discrete torsion" means the plain product above (**assumed**, as the
  definition of the orbifold, as in the paper). A nontrivial torsion changes the Witten index (for the Golay code from
  576 to 288), so the torsion-free orbifold is the one whose index is Index(C)².

## 2. Theorem A (pointwise identity; proved)

For every linear ternary code C of length n and every full-support codeword w of C,

    N(w) := Σ_{c ∈ C, c ≤ w} (−1)^{n − wt(c)} = (−1)ⁿ sign(w) Index(C).

In particular N(w) = sign(w) Index(C) for even n (e.g. self-dual codes).

*Proof.* Let φ_w multiply coordinate i by wᵢ ∈ {±1}. Then C′ = φ_w(C) is a linear code containing 1 = φ_w(w), and φ_w
maps {c ∈ C : c ≤ w} bijectively onto C′ ∩ {0,1}ⁿ, preserving weights. A full-support u has sign(φ_w(u)) = ∏ᵢ wᵢuᵢ =
sign(w) sign(u) (entries read as ±1), so Index(C′) = sign(w) Index(C). In C′, b ↦ 1 + b is a bijection from
C′ ∩ {0,1}ⁿ onto the full-support words (inverse u ↦ u − 1), with sign(1 + b) = (−1)^{wt(b)}. Hence
N(w) = Σ_b (−1)^{n − wt(b)} = (−1)ⁿ Index(C′) = (−1)ⁿ sign(w) Index(C). ∎
Self-duality is not used. For n = 0 the empty word has full support and both sides are 1. The statement is vacuous
when C has no full-support word.

*Corollary A (the Witten index of LG/C; the paper's gapping test).* For self-dual C,
(1/|C|) Σ_{c, c′} ∏ᵢ f(cᵢ, c′ᵢ) = Index(C)², where f(0,0) = −2 and f = 1 otherwise (the values F_{a,b}(τ, 0)), i.e.
(1/|C|) Σ_{c, c′ ∈ C} (−2)^{#{i : cᵢ = c′ᵢ = 0}} = Index(C)².
*Proof.* Write f(a, b) = Σ_k ω^{bk} g(a, k): g(0,0) = 0, g(0, k) = −1 for k ≠ 0, g(a, k) = [k = 0] for a ≠ 0. Then
(1/|C|) Σ_{c′} ∏ᵢ f(cᵢ, c′ᵢ) = Σ_{k ∈ C⊥} ∏ᵢ g(cᵢ, kᵢ) (the sum over c′ ∈ C is |C| times the projection onto C⊥), whose
nonzero terms have supp k = complement of supp c and value (−1)^{n − wt(c)}. With C⊥ = C, w = c + k runs over the
full-support codewords with c ≤ w, so the total is Σ_w N(w) = Index(C) Σ_w sign(w) = Index(C)² by Theorem A. ∎

## 3. Theorem B (per-field identity; proved)

For all a, k ∈ ℤ₃: G_{a,k} = − ψ_{2a+k} ι_{a+k}. (This holds for the functions as defined; the inputs of section 1 are
needed only to read F_{a,b} as the LG twisted-twined genera.)

*Proof.* Branches: X^{1/2} := e^{iπw}, as in the factor e^{iπu} of θ₁, and fixed integer representatives a, b.
(i) Put w = (z + aτ + b)/3. Then z/3 + (aτ + b)/3 = w and −2z/3 + (aτ + b)/3 = −2w + aτ + b. By quasi-periodicity,
θ₁(τ, −2w + aτ + b) = −(−1)^{a+b} q^{−a²/2} e^{4πiaw} θ₁(τ, 2w), and e^{4πiaw} = y^{2a/3} q^{2a²/3} ω^{2ab}, so

    F_{a,b} = −(−1)^{a+b} y^{a/3} q^{a²/6} ω^{2ab} R(w),     R(w) = θ₁(τ, 2w)/θ₁(τ, w).

(ii) With X = e^{2πiw}, the product formula gives, after splitting (1 − X²qⁿ)(1 − X⁻²q^{n−1}) into even and odd n,

    η(τ) R(w) = q^{1/24} X^{1/2} ∏_{n≥1} (1 − qⁿ)(1 + Xqⁿ)(1 + X⁻¹q^{n−1})(1 − X²q^{2n−1})(1 − X⁻²q^{2n−1}).

The quintuple product identity (**cited**: Watson 1929; survey S. Cooper, Int. J. Number Theory 2 (2006) 115–161)

    Σ_{n ∈ ℤ} q^{n(3n+1)/2} (x^{3n} − x^{−3n−1}) = ∏_{n≥1} (1 − qⁿ)(1 − xqⁿ)(1 − x⁻¹q^{n−1})(1 − x²q^{2n−1})(1 − x⁻²q^{2n−1})

(the form used is **computed** exactly to q¹⁴ by `../code/lg/quintuple_check.py`) at x = −X gives

    η R(w) = Σ_{n ∈ ℤ} (−1)ⁿ q^{(6n+1)²/24} (X^{3n+1/2} + X^{−3n−1/2}) = Σ_{x ∈ ℤ} ι_x (−1)^x X^{x+3/2} q^{(x+3/2)²/6}:

after n ↦ −n in the second term, the first term is x = 3n − 2 (x ≡ 1 mod 3, ι = +1) and the second x = 3n − 1
(x ≡ 2 mod 3, ι = −1).
(iii) Insert X = y^{1/3} q^{a/3} ωᵇ into (ii) and multiply by the prefactor in (i). Completing the square,
(x + 3/2)²/6 + a(x + 3/2)/3 + a²/6 = (x + a + 3/2)²/6, and y^{a/3} y^{(x+3/2)/3} = y^{(x+a+3/2)/3}. Shift x ↦ x − a.
The signs combine as (−1)^{a+b} (−1)^{x−a} ω^{3b/2} with ω^{3b/2} = e^{iπb} = (−1)ᵇ (not 3-periodic in b; it is
compensated by (−1)^{a+b}, which is why one representative is used throughout), giving (−1)^x, and
ω^{2ab} ω^{b(x−a)} = ω^{b(x+a)}. Grouping x by its residue u:

    F_{a,b} = − Σ_{u ∈ ℤ₃} ι_{u−a} ω^{b(u+a)} ψ_u.

Then G_{a,k} = (1/3) Σ_b ω^{−bk} F_{a,b} = − Σ_u ι_{u−a} ψ_u [u + a = k] = − ι_{k−2a} ψ_{k−a}, and k − a = 2a + k,
k − 2a = a + k in ℤ₃. ∎
Checks (**computed**): `../code/lg/per_field_genus.py` (exact in ℚ(ω), to q², all nine (a, k), searching GL₂(F₃) for the
relabeling; it finds two, among them the one above); an independent audit (23 October 2026) compared direct
theta-series division with the lattice sums exactly to q⁶ at two rational values of y^{1/6} (18 identities), with
negative controls without the prefactor y^{−a/3}.

## 4. Corollary (the refined identity of the paper's Appendix A; proved for the functions as defined)

For every self-dual ternary code C: Z_{LG/C}(τ, z) = Index(C) χ^R_V(τ, z).

*Proof.* Insert F_{a,b} = Σ_k ω^{bk} G_{a,k}: (1/|C|) Σ_{c′} ∏ᵢ F_{cᵢ, c′ᵢ} = Σ_{k ∈ C⊥} ∏ᵢ G_{cᵢ, kᵢ}. With C = C⊥
and Theorem B, Z_{LG/C} = Σ_{(c, k) ∈ C × C} ∏ᵢ (−ψ_{2cᵢ+kᵢ} ι_{cᵢ+kᵢ}). The map (c, k) ↦ (u, w) = (2c + k, c + k) has
determinant 2 − 1 = 1 and maps C × C onto C × C, so Z_{LG/C} = (−1)ⁿ (Σ_{u ∈ C} ∏ᵢ ψ_{uᵢ})(Σ_{w ∈ C} ∏ᵢ ι_{wᵢ})
= χ^R_V Index(C), since ∏ᵢ ι_{wᵢ} vanishes unless w has full support, where it is sign(w), and n is even. ∎
At y = 1 this is Corollary A; the q⁰ term is the charge-resolved form of Theorem A (`../code/lg/pointwise_identity.py`,
23 codes: the Golay code, the Pless code P24, the extended QR code of length 24 and 20 codes from random neighbour walks;
about 5 minutes). `../code/lg/lg_count_dense.py` checks Corollary A for P24 and QR24 (index 48, 2,304 = 48²).

## 5. Proposition C (all four spin structures, for the IR orbifold Mⁿ/C; proved for the stated data)

Let M be the ℓ = 1 N = 2 minimal model, the IR limit of W = Φ³ (**cited**: Vafa–Warner, Phys. Lett. B 218 (1989) 51;
Martinec, Phys. Lett. B 217 (1989) 431), written as the diagonal theory of V_√3ℤ: one field has labels x_L, x_R ∈ ℤ
(momenta X = x for NS, X = x + 3/2 for R; weights X²/6 plus oscillators), with x_L ≡ x_R mod 3 in the untwisted sector.
Data per field: the a-twisted spatial sector x_L − x_R ≡ a mod 3; the ℤ₃ phase of the temporal twist b is ω^{b x_L} in
the NS sector and ω^{b(x_L + a)} in the R sector (the phase appearing in F_{a,b} = −Σ_u ι_{u−a} ω^{b(u+a)} ψ_u with
u = x_L mod 3); fermion parity (−1)^{x_L + x_R}; no discrete torsion. The identification of Φ ↦ ωΦ with ω^{x_L} on the
untwisted sector of M is an **input**. Then the four torus partition functions of Mⁿ/C (NS and R sectors, each with and
without (−1)^F) equal those of V_Λ(C) ⊗ V̄_Λ(C).

*Proof.* The η factors and the fermion-parity signs are common to both sides (the insertion of (−1)^F multiplies both by
the same sign (−1)^{x_L + x_R}, which gives the other two spin structures), so it suffices to compare the weight of a
lattice point (x_L, x_R) ∈ (ℤⁿ)². On the orbifold side it is (1/|C|) Σ_{c ∈ C} [x_L − x_R ≡ c] Σ_{c′ ∈ C} ω^{c′·v} with
v = x_L (NS) or v = x_L + c (R); the inner sum is |C| [v mod 3 ∈ C⊥]. With C = C⊥ and c ∈ C, both cases give
[x_L mod 3 ∈ C][x_R mod 3 ∈ C], the weight of the double. ∎
*Remarks.* (1) A per-field bilinear phase ω^{ab} multiplies to ω^{c·c′} = 1 on C × C and drops out. This is also what
makes the NS phase ω^{b x_L} legitimate: per field it is not level-matched (for a = 1 the point (x_L, x_R) = (0, −1)
survives with spin h − h̄ = −1/6), but over C × C the offending phases cancel. An isolated coordinate with this chiral
labeling is therefore not a consistent orbifold by itself; nothing is inferred from that about LG anomalies.
(2) A constant extra per-field R phase ω^{κb} multiplies the sector (c, c′) by ω^{κ(1·c′)}. If 1 ∈ C (possible only when
12 | n) it is trivial. If 1 ∉ C the R projection changes, although the elliptic genus need not change (for the
tetracode it is 0 either way), so a vanishing or unchanged genus does not fix the R phases. Such a family is excluded by
T-covariance: under τ ↦ τ + 1 the sector (c, c′) goes to (c, c + c′) and picks up ω^{κ(1·c)}, nontrivial when
1 ∉ C⊥ = C. Checks (**computed**, exact in ℤ[ω]): `../code/lg/spin_structure_check.py` (tetracode: all points of the
box [−2, 2]⁸, 390,625 points, no mismatch in NS or R, 9,125 mismatches for each κ = 1, 2; Golay: 20,000 random points,
no mismatch, 9,981 for each κ; about 30 s). The independent audit of 23 October 2026 checked the tetracode projector on
all 6,561 residue pairs in each of the NS and R prescriptions (13,122 checks).

## 6. Scope and what is assumed

- Theorems A and B do not use self-duality. Self-duality enters in Corollary A, the corollary of section 4 and
  Proposition C. For self-orthogonal codes that are not self-dual nothing beyond Theorems A and B is claimed here.
- The UV-to-IR step "LG/C flows to Mⁿ/C" uses that W = Φ³ flows to M (cited) and that orbifolding by a geometric symmetry
  commutes with the flow (standard; **assumed**). Sections 2–4 are statements about the protected elliptic genus (UV and
  IR); section 5 is about the IR orbifold. Together they support Claim 1 of the paper: Mⁿ/C has the sectors and the four
  partition functions of V_Λ(C) ⊗ V̄_Λ(C). Identifying the full CFTs (the simple-current extension by the code, with its
  operator products) is the further step the paper takes as part of the claim.
- The phases of the twisted sectors are fixed by the flavor action and modular covariance (section 1), not by section 5;
  "no discrete torsion" is the choice of the plain product, the one with Witten index Index(C)².
- Self-orthogonality is the locality condition for the chiral code extension (each added word must be orthogonal to the
  code and to itself). It is not a condition for gauging the geometric, non-chiral phase symmetry of the LG model, which
  acts in the same way on both chiralities: for example W = Φ³ with its full ℤ₃ (the length-1 code F₃, not
  self-orthogonal) has an ordinary orbifold, whose index by the formula of Corollary A is (−2 + 8)/3 = 2; it just has no
  description as a diagonal code double.

## 7. The gapping deformation (the paper's §6.4)

The real superpotential h = Re Σᵢ Φᵢ³ + ε Σᵢ |Φᵢ|² preserves N = (1,1) supersymmetry and the phase symmetry of C. Per
field, with Φ = x + iy, h = x³ − 3xy² + ε(x² + y²), so ∂ₓh = 3x² − 3y² + 2εx and ∂ᵧh = 2y(ε − 3x). The critical points
are x = y = 0, (x, y) = (−2ε/3, 0) and (ε/3, ±ε/√3), i.e. Φ ∈ {0, −2ε/3, −2εω/3, −2εω²/3}: four non-degenerate critical
points **per field**, 4ⁿ before orbifolding. The Hessian is [[6x + 2ε, −6y], [−6y, −6x + 2ε]]: at the origin it is 2ε·1,
a minimum for ε > 0 (bosonic vacuum, fixed by ℤ₃); at (−2ε/3, 0) it is diag(−2ε, 6ε), Morse index 1 (fermionic), and
the other two are its ℤ₃ images. In the gapped theory a coordinate with both holonomies trivial contributes 1 − 3 = −2 to the twisted-twined Witten
index and otherwise 1 (only the fixed vacuum at 0), which gives Corollary A as the Witten index of the gapped orbifold.
This gaps the diagonal double of a single seed (index Index(C)², 576 for the paper's seed); it does not gap P̄ ⊗ P with
the chiral symmetries G_L × G_R that build T, T′, P = T″, which act on one chirality only and are invisible in the LG
description. Flows that keep the (2,2) superconformal algebra end at theories with U(1) R-currents, so the current-free
T, T′, P of the paper's seed are not endpoints of such flows; this says nothing about (1,1) or (0,1) UV theories.
