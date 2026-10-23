# Power operations on TMF and the symmetric square: does "P_2 then gauge" reproduce the index arithmetic?

Note of 18 October 2026 (Claude); prepared for publication on 20 October (remarks about internal workflow removed or
reworded; the mathematics, the corrections and the data are unchanged). Each claim is labelled
**proved** / **cited** (with theorem number where visible) / **computed** / **assumed** / **open**.

Question (from the paper's last remark): the physics index arithmetic of the
symmetric square of a holomorphic N=1 SCFT is I -> I(I+3)/2 (plain) and I -> I(I-3)/2 (Arf-twisted). Does the
homotopy-theoretic composite "total power operation P_2 on TMF, then the Sigma_2 transfer" send
-4 Delta^6 (pi_144) to 2 Delta^12 (pi_288), and (Arf-twisted) -2 Delta^12 to Delta^24 (pi_576), at least modulo torsion?

Sections: S1 literature; S2 formulation; S3 sympy checks; S4 integrality/torsion; S5 D_4 interpretation; S6 verdict.

## S1. Literature (bounded; statements read from the arXiv PDFs, text extracted with ghostscript)

### S1.1 Ganter, "Stringy power operations in Tate K-theory", arXiv:math/0701565 (v4, 2013) — **cited**

Setting: Devoto's equivariant Tate K-theory K_{Dev,G}(X) (twisted sectors indexed by [g], coefficients in
R(C_g)[[q^{1/|g|}]]), with a level-2 Hopkins–Kuhn–Ravenel character theory: an element is a class function on
commuting pairs (g,h) with values in q-series (Def. 1.1 and §3.1).

* **Def. 1.1** (orbifold genus = "gauging"): for χ ∈ Ell_G(pt) with HKR theory, ε_G(χ) := (1/|G|) Σ_{gh=hg} χ(g,h).
  This is literally the physics orbifold formula (sum over commuting pairs = G-bundles on the torus).
* **Def. 5.1**: H_∞-structure on a compatible family of equivariant theories = power operations
  P_n : E_G(X) → E_{G≀Σ_n}(X^n) with (a) P_1 = id, P_0 = 1; (b) P_j(x)P_k(x) = res^{Σ_{j+k}}_{Σ_j×Σ_k} P_{j+k}(x);
  **(c) P_j(P_k(x)) = res^{Σ_{jk}}_{Σ_k≀Σ_j} P_{jk}(x)**; (d) multiplicativity. (c) is the composition law used in S5.
* **Def. 5.2 / eq. (12)**: the "Atiyah power operations" on K_G(X)[[q]] (untwisted sector):
  P^top_n(Σ_j x_j q^j) = Σ_{j∈N^n} x_{j_1}q^{j_1}⋯x_{j_n}q^{j_n}, with Σ_n permuting the factors (with sign);
  in particular x q^j ↦ x^{⊗n} q^{jn}. **Prop. 5.3**: these satisfy Def. 5.1.
* **Def. 5.9–5.10, Thm. 5.11**: the stringy power operations P^string_n on all of K_{Dev,G}(X), built from P^top
  on each cycle type via the substitution k(−): q^{j/(r|g|)} ↦ q^{j/(rk|g|)} (the "long circle" of length k);
  Thm. 5.11: they satisfy the axioms of Def. 5.1.
* **Thm. 5.13**: the stringy Todd-genus map Td^string_G : MU_G → K_{Dev} is **H^2_∞** (i.e. power operations
  on classes of even degree, via Thom classes; H^d_∞ in the sense of Bruner–May–McClure–Steinberger).
  **Prop. 5.5**: the topological Witten genus φ^top_G : MSpin_G → K_G[[q]] is H^4_∞ (H^2_∞ in the complex case).
* **Def. 5.14, Remark 5.15 (Hecke operators)**: T_n(x) := (1/n) Σ_{α∈A} ψ_α(x), ψ_α(x) = eval_α(res_δ P_n(x)),
  α running over conjugacy classes of commuting pairs (σ,τ) in Σ_n generating a transitive subgroup.
  On q-expansions (G trivial): (T_n x)(z) = (1/n) Σ_{Nk=n, 0≤m<k} x((Nz+m)/k) — the classical Hecke operator;
  "If G is the trivial group, we recover the operators of [And00]" (Ando, Trans. AMS 2000).
* **§6, Prop. 6.1, Cor. 6.2 (DMVV)**: sym^n_string(x) := (1/n!) Σ_{α commuting pairs in Σ_n} ψ_α(x) and
  Sym_t(x) = Σ_n sym^n t^n = exp(Σ_{m≥1} T_m(x) t^m); for X = pt this is the DMVV formula.
  **Consequence (proved, trivial algebra):** for n = 2 the four commuting pairs (1,1),(1,s),(s,1),(s,s) give
  sym^2_string(x)(z) = (1/2)[x(z)^2 + x(2z) + x(z/2) + x((z+1)/2)], i.e. exactly the physics formula for the
  symmetric square of a weight-0 Ramond partition function (untwisted / swap-twined / twisted / twisted-twined).
* **Thom classes and the c/24 factor (Def. 3.4, §3.2, Prop. 3.9)**: for a G-equivariant U²(d)-bundle the
  stringy Thom class is u^string_g(V) = q^{d/12 − age(g)/2} · u^Td(V_1) · (√det V)^{−1} · ∏_{j≥1} Λ_{−q^{j/l}}(V_C)_{ζ^j};
  the q^{d/12} comes from ζ-regularising Σ_{j>0} j = −1/12 (§3.2). Prop. 3.9: the orbifold genus takes values in
  Z[[q]][q^{±1/12}]. So in Ganter's framework the fractional power q^{d/12} is **part of the Thom class**, and the
  H_∞ property P_n(u(V)) = u(V^n) forces P_n(q^{d/12}·…) = q^{nd/12}·… in the untwisted sector, consistent with
  Def. 5.2 (x q^j ↦ x^{⊗n} q^{jn}); in the twisted sectors (k-cycles) the substitution q ↦ q^{1/k} is applied to the
  whole series including the fractional prefactor (Def. 5.9), i.e. q^{d/12} ↦ q^{d/(12k)}. Note Def. 5.12/Thm 5.13
  prove the H^2_∞ property for the *Td* variant, in which "the terms q^{−d/12+age(g)/2} e^{iπ age(h)} ∏ e^{x_i}"
  are dropped; the paper remarks that for SU-bundles these are "just a renormalization constant".

### S1.2 Ganter, "Orbifold genera, product formulas and power operations", arXiv:math/0407021, Adv. Math. 205 (2006) 84–133 — **cited**

* §4.1 Def. (H_∞ ring spectra and families of equivariant power operations); **Prop. 4.5** ("compare [BMMS86, VIII.1.1]")
  lists the axioms of the P_n (sums to products, composition, restriction); **Prop. 4.6/Cor. 4.7**: Thom classes
  compatible with power operations ⇒ the equivariant genus φ_G preserves power operations.
* §7.4–7.5 generalized orbifold genera and symmetric powers; **Prop. 7.9**: over Q the augmentation/"orbifold"
  map η_*: E^0(BG) → E^0 is the HKR average over commuting h-tuples; **Cor. 7.11**: orbifold genus is integral.
* **Thm. 9.2 (= Thm. 1.5)**: for an H_∞-orientation φ of a Morava E-theory E_h,
  Σ_n φ_orb(M^n//Σ_n) t^n = exp(Σ_m T_m φ(M) t^m), the p-typical DMVV formula; §9.3 identifies T_m with
  Atiyah–Tall–Grothendieck-type Hecke operators (averages over isogenies).

### S1.3 Other references (abstract pages fetched) — **cited (abstracts)**

* Ando–Hopkins–Strickland, arXiv:math/0204053: "if E is a K(2)-local H_∞ elliptic spectrum, then the sigma
  orientation is a map of H_∞ spectra" (abstract). Ganter's stringy paper notes that restricted to MU⟨6⟩ the Witten
  genus becomes the σ-orientation of K_Tate (§1, citing [AHS01]).
* Ando–Hopkins–Rezk, "Multiplicative orientations of KO-theory and of tmf" (rezk.web.illinois.edu/koandtmf.pdf,
  v0.5 2010): the Witten genus is the value on homotopy of an E_∞ map MString → tmf (abstract). Theorem numbers
  for the E_∞ map tmf → KO((q)) extracted below if the text is readable.
* Ganter, "Power operations in orbifold Tate K-theory", arXiv:1301.2754: axioms of an orbifold theory with power
  operations; orbifold Tate K-theory with power operations; "symmetric powers, exterior powers and Hecke operators"
  (abstract). Not read further.
* Barthel–Stapleton, arXiv:1502.01987 (HTML v3): **Thm. 9.1 (Thm. A)**: the HKR character map intertwines the total
  power operation P_m with P^φ_m on class functions, where (§5, before Def. 5.6) P^φ_m(f) at a class
  ⊕_i (H_i,[α_i]) of G≀Σ_m is the product over summands of φ*_{H_i} f([α_i ∘ ψ*_{H_i}]): a product over orbits of
  f transported along the isogeny φ_{H} attached to the subgroup H ⊂ Z_p^n (height n). Thm. 10.1: after inverting p
  the formula is independent of choices. This is the height-n analogue of Ganter's Remark 5.15.


## S2. Formulation: Sym²(f) := ε_{Σ_2}(P_2(f)) on π_{2w}

### S2.1 Degree zero — **cited + proved (trivial algebra)**

For E = Tate K-theory with Ganter's stringy power operations (S1.1) and G = Σ_2, X = pt, define
Sym²(x) := ε_{Σ_2}(P_2^string(x)) = (1/2) Σ_{(g,h)∈Σ_2², gh=hg} P_2^string(x)(g,h)  (Ganter Def. 1.1 + §6).
By Remark 5.15 / Prop. 6.1 (n = 2: the commuting pairs are (1,1), (1,s), (s,1), (s,s); the pair (1,s) is the
Adams operation ψ², k=1,N=2; the pairs (s,1),(s,s) are the long circle k=2,N=1 with m=0,1):

    Sym²(x)(τ) = (1/2)[ x(τ)² + x(2τ) + x(τ/2) + x((τ+1)/2) ].                                   (2.1)

This is exactly the DMVV/physics formula for the symmetric square of a weight-0 Ramond partition function,
and ε_{Σ_2} ("gauging") is literally the orbifold sum over the four Σ_2-bundles on the torus.

**Status of "gauging = transfer".** (i) In Ganter's framework ε_G is a *definition* (Def. 1.1) made to match the
physics. (ii) It is a *theorem* that ε_G agrees with the K(h)-local transfer/Strickland inner product for Morava
E-theory: Ganter (math/0407021) Thm. 7.6 (citing Strickland, "K(n)-local duality for finite groups and groupoids",
Topology 39 (2000)) and Prop. 7.9: "Over Q the augmentation map η_*: E^0(BG) → E^0 is" the HKR average over
commuting h-tuples. So for h = 2 (K(2)-local TMF, i.e. Morava E_2 at a supersingular curve) "P_2 then the K(2)-local
transfer" is, rationally, the sum over commuting pairs. (iii) For integral tmf (not K(2)-localised) a transfer
tmf^0(BΣ_2) → π_0 tmf that reproduces the commuting-pair average is an **expectation** (it should be the
pushforward along M_{Σ_2}(E) → M_{Ell} in genuine Σ_2-equivariant elliptic cohomology, Lurie/Gepner–Meier); the
Borel transfer BΣ_2+ → S^0 goes the other way and is *not* it. We label (iii) **assumed** below.

### S2.2 Nonzero degree: Thom classes — **proved (bookkeeping), cited, computed**

A class f ∈ π_{2w}E is a map S^{2w} → E; its power operation lands on the extended power
D_2(S^{2w}) = (BΣ_2)^{V_w}, the Thom spectrum of V_w := C^w ⊗ ρ_2 = C^w ⊕ (C^w ⊗ sgn) over BΣ_2
(ρ_2 = permutation representation). To turn P_2(f) ∈ E^0((BΣ_2)^{V_w}) into an element of π_{4w}E by a transfer
one needs three things:
 (i) an E-Thom class u(V_w) for V_w over BΣ_2, giving E^0((BΣ_2)^{V_w}) ≅ E^{−4w}(BΣ_2);
 (ii) a transfer E^{*}(BΣ_2) → π_{−*}E (S2.1: K(2)-locally Strickland's transfer = HKR average; for Devoto/Ganter
     Tate K-theory the orbifold sum ε_{Σ_2} of Def. 1.1);
 (iii) the nonequivariant degree-4w Thom class u(C^{2w}) to land in π_{4w}.
(This is the H^d_∞ formalism of Bruner–May–McClure–Steinberger, used by Ganter math/0407021 end of §4.3:
"we can follow [GM97] to extend the power operations to P_n: E^{2d}(X) → E^{2nd}(X^n)" via Thom isomorphisms.)
So Sym²(f) := ε_{Σ_2}( P_2(f) / u(V_w) ) · u(C^{2w}).

**The Thom classes are those of the String (σ/Witten) orientation** (Ganter math/0701565 Def. 3.4; for K_Tate
the σ-orientation of MU⟨6⟩ is the Witten genus, AHR §1 / AHS01; AHR Thm 12.3: tmf has E_∞ String orientations
refining the σ-orientation). Characters of these Thom classes (zero modes give a Bott/Â factor, nonzero modes give
Λ_{−q^{j/l}}, zero point by ζ-regularisation, Ganter §3.2):
 * nonequivariant: u(C^w) = q^{w/12} ∏_{n≥1}(1−q^n)^{2w} · β^w = η(τ)^{2w} β^w = **Δ^{w/12} β^w** (w = 12k: Δ^k β^w);
 * u(V_w) in the four sectors (computed in S3 from mode products):
   (1,1): Δ(τ)^{2k};  (1,s): 2^w Δ(2τ)^k;  (s,1): Δ(τ/2)^k;  (s,s): Δ((τ+1)/2)^k,
   where in the twisted sectors the C^w⊗sgn modes are antiperiodic, ∏_{j odd}(1 ∓ q^{j/2})^{2w}, with zero point
   q^{−1/24} per complex dimension (ζ(−1,1/2) = +1/24, versus q^{+1/12} from ζ(−1) = −1/12 for periodic modes) and,
   in sector (s,s), the phase (−1)^k from the s-action on the twisted ground state (√det N)^{−1}.
   **(A1, assumed/conventional):** the ζ-regularised zero point q^{−1/24} per antiperiodic complex dimension. I could
   not read Ganter's literal age exponent from the extracted PDF; A1 is the value for which (H) below holds.

**Key point (proved + computed):** the String Thom classes are H_∞-compatible, P_2(u(C^w)) = u(V_w). For the
untwisted sector this is Ganter's Prop. 5.5 (φ^top is H_∞); for all four sectors it is the content of her
Thm. 5.13 (stated for the Td-normalised variant) and is checked here directly on q-expansions: S3 (H) verifies
P_2^string(Δ^k) = u(V_w) in all four sectors for k = 1, 2, 6. Writing f = Z · u(C^w) (Z = a/Δ^k of degree 0,
the physics partition function), multiplicativity gives P_2(f) = P_2(Z) · u(V_w), hence

    Sym²(f) = ε_{Σ_2}(P_2(Z)) · u(C^{2w}) = (1/2)[ Z(τ)² + Z(2τ) + Z(τ/2) + Z((τ+1)/2) ] · Δ(τ)^{2k}.          (2.2)

**This is exactly the physics prediction** (the original question's formula with weight factors (Δ(τ)²/Δ(Λ_H))^{w/12} is (2.2)
rewritten in terms of f). In particular **no extra Euler-class factor survives**: the Euler classes of C^w⊗sgn
(the characters 2^wΔ(2τ)^k, Δ(τ/2)^k, Δ((τ+1)/2)^k above) appear in both P_2(f) and u(V_w) and cancel. They are
units in Z((q^{1/2}))[1/2] — unlike in K-theory, where the sign-representation Euler class 1 − sgn is a zero
divisor. That invertibility (the twisted sector has no zero modes) is what makes a degree-shifted Σ_2 operation
possible at all.

### S2.3 What goes wrong with other choices — **computed (S3)**

* **Borel / K(1)-local transfer** (only the untwisted sectors (1,1), (1,s); this is what one gets from tmf → K((q))
  followed by the Borel transfer in K((q)), height 1): Sym²(IΔ^k) = (I²+I)/2 Δ^{2k}; −4Δ^6 ↦ 6Δ^{12}, not 2Δ^{12}.
  So the twisted sectors (height-2 information: commuting *pairs*) are essential.
* **No Thom isomorphism** (averaging P_2(a) without dividing by u(V_w)): (1/2)[a² + Σ_H a(τ_H)] mixes weights 2w and
  w (leading term q^{3} = Δ(τ/2)^6 for k = 6): not a class of degree 4w.
* **Other Thom class:** the only freedom is the choice in S2.4, which produces the Arf twist.

### S2.3' Structure of the three twisted sectors — **interpretation (cited at abstract level)**

The three nontrivial Σ_2-bundles on T² are permuted transitively by SL(2,Z) with stabiliser Γ_0(2) (for (1,s)), so
Σ_H Z(τ_H) = Tr_{Γ_0(2)}^{SL(2,Z)} Z(2τ). In genuine Σ_2-equivariant elliptic cohomology (Lurie; Gepner–Meier,
"On equivariant topological modular forms", arXiv:2004.10254, Compositio, abstract: equivariant TMF for compact
abelian groups via Lurie's construction) the Σ_2-fixed points live over the 2-torsion E[2] = {0} ⊔ E[2]^×, and
E[2]^× ≅ the Γ_0(2)-moduli; the transfer to TMF is the pushforward along E[2] → M_Ell, i.e. identity on the
{0}-summand plus the trace from Γ_0(2) on the other. That is the algebraic shape of (2.2). Identifying this
pushforward, composed with tmf's E_∞ power operation P_2, with (2.2) on the nose is **open** (not in the
literature I could find); rationally it follows from the K(2)-local statement plus the q-expansion principle.

### S2.4 The two spin structures / the Arf twist — **proved (counting) + interpretation**

What has to be oriented is the virtual Σ_2-bundle V_w − C^{2w} = C^w ⊗ (sgn − 1) = w(L − 1) over BΣ_2, with
L = sgn ⊗ C, c_1(L) = y ∈ H²(BZ/2;Z) = Z/2. For w = 12k: w_1 = 0; w_2 = w·y mod 2 = 0; 2(L−1) is spin with
p_1/2 = ±y² ∈ H⁴(BZ/2;Z) = Z/2, so w(L−1) = 6k·2(L−1) has p_1/2 = 6k·y² = 0. Spin structures on it form a torsor
over H¹(BΣ_2;Z/2) = Z/2; string structures over a fixed spin structure form a torsor over H³(BΣ_2;Z) = 0.
**Hence there are exactly two Σ_2-equivariant string
structures on V_w − C^{2w} (proved, standard characteristic-class computation), i.e. exactly two Thom isomorphisms
in step (i) of S2.2 built from the String orientation.** One of them is the H_∞-canonical one, u(V_w) = P_2(u(C^w)),
and gives the plain formula (2.2). Their
Thom classes differ, sector by sector, by a sign ε(g,h) ∈ {±1} with ε(1,1) = 1; in K-theory the change of spin structure by the nonzero
class of H¹(BΣ_2;Z/2) multiplies the Thom class by the sign representation, so ε(1,s) = −1; modular covariance
(SL(2,Z) permutes the three nontrivial Σ_2-bundles on T²) then forces ε = (+1, −1, −1, −1) on
((1,1),(1,s),(s,1),(s,s)) (computed in S3 (M'): the two tested patterns with ε(1,s) = −1 but not all three twisted signs −1
do not give modular forms).
That is exactly the Arf-twisted symmetric square: (I²+3I)/2 ↦ (I²−3I)/2. **Interpretation (not proved here):** the
two equivariant spin structures on the swap correspond to Sym²(T) and Sym²(T⊗Arf); the sign ε(g,h) is
(−1)^{Arf(ρ_{RR} + α_{g,h})+1} for the Z/2-bundle α_{g,h} on T².


## S3. Exact q-expansion checks — **computed**

Script: `power_operations_check.py`, included as `../code/power_operations_check.py` (pure Python, exact `Fraction` arithmetic, series in t = q^{1/2}
with tracked precision; ≈2 min single thread at precision q^60). Output: `power_operations_check.out` (workspace, not included; the table below copies it).
Implementation: P_2^string sector characters per Ganter Def. 5.2/5.9/Rem. 5.15 ((1,1): a², (1,s): 2^w a(q²),
(s,1): a(q^{1/2}), (s,s): coefficient of q^{n/2} times (−1)^n); u(V_w) built **from mode products** (not from Δ
formulas) with the zero points of S2.2; Sym²_ε(f) := ½ Σ ε(g,h) P_2(f)(g,h)/u(V_w)(g,h) · Δ^{2k}.

| check | result |
|---|---|
| (H) u(C^w) = Δ^k and P_2^string(u(C^w)) = u(V_w) in all four sectors, k = 1, 2, 6 | True |
| (A) Sym²_plain(−4Δ^6) = 2Δ^{12} (to q^41) | True |
| (A) Sym²_Arf(2Δ^{12}) = −Δ^{24} (to q^53) | True |
| (A) Sym²_plain(−2Δ^{12}) = −Δ^{24} | True |
| (A) Sym²_Arf(−2Δ^{12}) = 5Δ^{24} | True |
| (A) untwisted-only (Borel/K(1)) on −4Δ^6 gives 6Δ^{12} (not 2Δ^{12}) | True |
| (W) averaging without the Thom isomorphism: leading term q^3 (wrong weight) | as stated |
| (M) f = IΔ^k + cE_4^{3k} (k = 1, 2, 6; both ε): Sym²(f) ∈ M_{24k}, integral q-coefficients | True (all six) |
| (M') sign patterns (+,−,+,+), (+,−,−,+): Sym² **not** in M_24 | True (both fail, as they should) |

What this establishes:
* **Index arithmetic (computed; also proved by one line from (2.2)):** for f = IΔ^k, Z = I is constant and (2.2)
  gives (I² + 3I)/2 · Δ^{2k} (plain) and (I² − 3I)/2 · Δ^{2k} (other string structure). −4Δ^6 ↦ 2Δ^{12};
  2Δ^{12} ↦ −Δ^{24} with the Arf choice.
* **Sign remark on the original question.** "(with the Arf twist) −2Δ^{12} ↦ Δ^{24}" should read: the Arf-twisted square
  of **+2Δ^{12}** (= the plain square of −2Δ^{12}, since Sym²_Arf(f) = Sym²_plain(−f) identically by (2.2)) is
  **−Δ^{24}**. The Arf-twisted square of −2Δ^{12} is 5Δ^{24}. The physics chain (index −4 → 2 → −1) gives the class
  −Δ^{24}, i.e. "±Δ^{24}" as in the paper. The identity Sym²_ε(f) = Sym²_plain(−f) is the TMF shadow of
  Sym²(T ⊗ Arf), since the class of T ⊗ Arf is −[T] (Arf has index −1).
* **SL(2,Z)-covariance and integrality (computed):** (2.2) maps M_w (tested on non-pure-Δ inputs) to M_{2w} with
  integral coefficients, for both string structures; mixed sign patterns break modularity, so the Arf pattern is
  the only alternative. The trivial identity Σ_H f(Λ_H)(Δ²/Δ(Λ_H))^{w/12} = 3Δ^{2k} for f = Δ^k is (H) in disguise.
* **Not tested here:** that (2.2) is *literally* the composite of tmf's E_∞ power operation with a transfer in some
  model of tmf (see S2.1(iii), S2.3'). What is verified is the formula in Ganter's equivariant Tate K-theory with
  String Thom classes, plus its consistency (modularity, integrality, H_∞-compatibility of the Thom classes).

## S4. Integrality and torsion — **cited**

* **Image of π_* tmf in MF_*** — Hopkins, "Algebraic topology and modular forms", ICM 2002, arXiv:math/0212397,
  **Prop. 4.6** (Hopkins–Mahowald): the image of π_{2*}tmf → MF_* has basis a_{i,j,k} c_4^i c_6^j Δ^k with
  a = 1 (i > 0, j = 0), 2 (j = 1), and **24/gcd(24,k) (i = j = 0)**. Hence:
  −4Δ^6: 24/gcd(24,6) = 4 → in the image, minimal. 2Δ^{12}: 24/gcd(24,12) = 2 → in the image, minimal.
  Δ^{24}: 24/24 = 1 → in the image. So all three classes of the chain are realised by tmf classes and are the
  minimal positive multiples in their degrees (the "TMF-minimal" property).
* **Torsion in the relevant degrees** — Bruner–Rognes, *The Adams spectral sequence for topological modular forms*,
  AMS Surveys 253 (2021) (preliminary version, Rognes's home page):
  - Introduction: the kernel and cokernel of the edge homomorphism are torsion groups annihilated by 24 (stated for
    TMF), so torsion is only 2- and 3-primary.
  - p = 2: **Thm. 9.27** π_*(tmf) ≅ (Γ_B N_* ⊕ ⊕_{k=0}^{7} ko[k]) ⊗ Z[M], |M| = 192; **Thm. 9.26 / Table 9.4**:
    Γ_B N_* is finite and concentrated in degrees 3 ≤ * ≤ 164, with **no entry in degrees 96 or 144** (rows run …,
    90, 99, …, 142, 147, …); the remaining 2-torsion is the η-family Z/2[B]{η_k, ηη_k, …} of the ko[k], in degrees
    ≡ 1, 2 mod 8. Degrees 144, 288 = 96 + 192 and 576 = 3·192 are ≡ 0 mod 8 and meet no row of Table 9.4.
    **So π_{144}, π_{288}, π_{576} tmf have no 2-torsion.**
  - p = 3: π_*(tmf)_3 ≅ N_* ⊗ Z[H], |H| = 72 (Ch. 13 introduction, Thm. 13.18); the B-power torsion
    (= all 3-torsion) Γ_B N_* "is finite and concentrated in degrees 3 ≤ * ≤ 40" (§13.4, proof of the duality
    statement before Prop. 13.21). Degrees 144, 288, 576 are ≡ 0 mod 72, so **no 3-torsion** either.
  - **Conclusion (cited + one-line deduction):** π_{144}tmf, π_{288}tmf and π_{576}tmf are torsion-free, and the edge
    map to MF is injective in these degrees. "Modulo torsion" is vacuous for this question. The earlier project
    note "torsion in degrees 288, 576" is not supported by Bruner–Rognes. The torsion
    used by Tachikawa–Yamashita (heterotic anomalies) sits in other degrees.
  - Caveat: this is for connective tmf. In periodic TMF the degrees are the same mod 576, and the 2- and
    3-primary periodicities 192 | 576 and 72 | 576 give the same answer.


## S5. The D_4 = Σ_2 ≀ Σ_2 interpretation — **cited (axiom) + computed + interpretation**

* **Composition law (cited):** for an H_∞ structure, P_j(P_k(x)) = res^{Σ_{jk}}_{Σ_k ≀ Σ_j} P_{jk}(x)
  (Ganter math/0701565 **Def. 5.1(c)**; Ganter math/0407021 **Prop. 4.5**, "compare [BMMS86, VIII.1.1]"; the
  underlying map D_j D_k E → D_{jk} E of extended powers is in Bruner–May–McClure–Steinberger, *H_∞ ring spectra
  and their applications*, LNM 1176 (1986), Ch. I; I did not re-read the chapter or its statement numbers).
  Ganter's stringy operations satisfy it (**Thm. 5.11**, via Prop. 4.11 / Cor. 4.12 on iterated symmetric powers).
* **Sym² ∘ Sym² as a D_4 orbifold (proved modulo one axiom):** Sym²(Sym²(x)) = ε_{Σ_2} P_2(ε_{Σ_2} P_2 x)
  = ε_{D_4}(res^{Σ_4}_{D_4} P_4(x)). This uses (c) and the compatibility of power operations with transfers
  (P_2 ∘ ε_H = ε_{H≀Σ_2} ∘ P_2), which is part of being a global power functor. Barthel–Stapleton
  arXiv:1502.01987 **Thm. C (Thm. 10.4, Cor. 10.5)** proves that the HKR character map is a map of global power
  functors (rationally, Morava E-theory).
* **Computed** (`../code/power_operations_d4.py`): for constant Z = I the sum over all 40 commuting pairs in D_4,
  each weighted by I^{#orbits of ⟨a,b⟩ on 4 points}, gives (1/8)Σ = Sym²(Sym²(I)) for I = −6..6. With the sign −1
  on pairs whose outer Z_2-bundle is nontrivial it gives Sym²_Arf(Sym²(I)). For I = −4: plain/plain gives 5 and
  plain/Arf gives −1.
* **The spin twist (interpretation):** the outer Arf sign is the change of D_4-equivariant spin structure on the
  virtual representation being oriented. The change is by the class α ∈ H¹(BD_4; Z/2) of the outer homomorphism
  D_4 → Z_2. On the torus sector of a commuting pair (a,b) it contributes (−1)^{Arf(ρ_RR + α(a,b)) − Arf(ρ_RR)},
  which is −1 exactly when α(a,b) ≠ 0. That matches S2.4: on the TMF side the "spin twist" of the D_4 orbifold is
  the choice of one of the two equivariant string structures, i.e. the choice of Thom class in the H^d_∞
  structure. It is not an extra invertible phase in the transfer. Equivalently, it multiplies the outer factor by
  [Arf] = −1 ∈ π_0 TMF.


## S6. Verdict

**Answer.** Yes, at the level where power operations on elliptic cohomology are explicitly known, and the
answer reproduces the physics exactly. The topological operation is defined as follows:
"P_2, then Thom isomorphism with the String Thom class of C^w ⊗ ρ_2, then the Σ_2 transfer over commuting pairs".
In Ganter's equivariant Tate K-theory (and K(2)-locally, rationally, in Morava E_2 via Strickland's transfer
and Barthel–Stapleton) it sends a degree-2w class ZΔ^k to ½[Z(τ)² + Z(2τ) + Z(τ/2) + Z((τ+1)/2)]·Δ^{2k}.
That is precisely the DMVV/physics formula. The Euler classes of the sign representation cancel against the
H_∞-compatible Thom class and leave no extra factor. On f = IΔ^k this is (I² + 3I)/2·Δ^{2k}. The other
Σ_2-equivariant string structure gives (I² − 3I)/2·Δ^{2k}, and there are exactly two such structures.
So the chain is −4Δ^6 ↦ 2Δ^{12} ↦ −Δ^{24}, with the Arf choice at the second step (S3, computed). All three
classes are the minimal tmf classes in their degrees (Hopkins Prop. 4.6). π_{144}, π_{288} and π_{576} tmf are
torsion-free (Bruner–Rognes), so no "modulo torsion" qualifier is needed for the *values*.

**What is not established (assumptions / open):**
1. **(open)** That the composite is induced by an operation on *integral* tmf. That would need P_2 on tmf (E_∞, so
   it exists), followed by a transfer that sees the twisted sectors. The Borel transfer does not see them (S2.3:
   it gives 6Δ^{12}). The natural candidate is the pushforward in genuine Σ_2-equivariant TMF (Lurie;
   Gepner–Meier arXiv:2004.10254), but I found no reference that computes it on q-expansions. What is proved is the
   formula after mapping to Tate K-theory / K(2)-local E-theory, where the operations are Ganter's
   (Def. 5.10, Thm. 5.11) resp. the HKR characters of the E_∞ power operations (Barthel–Stapleton Thm. 9.1).
   The comparison with tmf's own E_∞ structure is cited only via AHS math/0204053 (σ-orientation H_∞ for
   K(2)-local elliptic spectra) and AHR Thm. 12.3 (E_∞ String orientation of tmf refining σ).
2. **(assumed, A1)** The ζ-regularised zero point q^{−1/24} per antiperiodic complex dimension in the twisted
   Thom class. It is forced by H_∞-compatibility, as check (H) shows. I could not read Ganter's literal exponent
   from the PDF text.
3. **(definition)** "Gauging = transfer" is Ganter's Def. 1.1 in Tate K-theory. It is a theorem K(2)-locally over Q
   (Ganter math/0407021 Prop. 7.9 with Strickland's transfer).
4. **(interpretation)** The second string structure is the Arf-twisted symmetric square. This is supported by the
   identity Sym²_Arf(f) = Sym²_plain(−f), which is proved, and by [Arf] = −1.
5. **Caveat on content.** For f = IΔ^k the partition function is constant. The whole computation then reduces to
   counting the 1 + 3 Σ_2-bundles on T² with signs, and the nontrivial input is that the Thom-class bookkeeping
   introduces no further factor. So this is a consistency check of the physics with the E_∞ structure. It does not
   give independent evidence that the physical SCFTs exist.

**Suggested wording for the paper's last remark** (replacing "… is a statement about power operations on TMF that
we have not verified"):

> Our construction is compatible with the E_∞-ring structure of TMF. The symmetric square is the field-theoretic
> shadow of the total power operation P_2 followed by the Σ_2-transfer, the sum over commuting pairs, i.e. over
> Σ_2-bundles on the torus. The degree shift is absorbed by the String Thom class of the Σ_2-representation C^w ⊗ ρ_2.
> In equivariant Tate K-theory, where these operations were constructed by Ganter [Ganter07, Def. 5.10,
> Thm. 5.11, Rem. 5.15], the composite sends Z·Δ^k to ½[Z(τ)² + Z(2τ) + Z(τ/2) + Z((τ+1)/2)]·Δ^{2k}, hence I·Δ^k to
> ½(I² + 3I)·Δ^{2k}. The other Σ_2-equivariant string structure gives ½(I² − 3I)·Δ^{2k}. Thus −4Δ^6 ↦ 2Δ^{12} ↦ −Δ^{24},
> in agreement with the index arithmetic above. Since π_{288}tmf and π_{576}tmf are torsion-free [BrunerRognes,
> Thm. 9.27, Table 9.4, Ch. 13], these classes are determined by their q-expansions. The two steps combine into the
> restriction of P_4 along Σ_2 ≀ Σ_2 ⊂ Σ_4. We have not verified that this composite is induced by a transfer in
> genuine Σ_2-equivariant tmf; this would make the statement one about π_* tmf itself.

(Bib keys to add: Ganter07 = arXiv:math/0701565; Ganter06 = Adv. Math. 205 (2006) 84–133, arXiv:math/0407021;
BrunerRognes = AMS Math. Surveys Monogr. 253 (2021); Hopkins ICM 2002 = arXiv:math/0212397; optionally
BarthelStapleton = arXiv:1502.01987, GepnerMeier = arXiv:2004.10254.)
