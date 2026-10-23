# Landau–Ginzburg orbifolds and refined elliptic genera: checks (details.html §L)

Checks of the identities proved in `../../notes/lg_elliptic_genus.md` (the paper's §6.3–6.4 and Appendix A). Plain Python 3
with numpy; exact arithmetic (integers, fractions, ℚ(ω), ℤ[ω]). Run from this folder, e.g. `python3 -B quintuple_check.py`.

| script | what it checks | time (one core) |
|---|---|---|
| `quintuple_check.py` | the quintuple product identity in the form used in Theorem B, exactly to q¹⁴ | < 1 s |
| `per_field_genus.py` | Theorem B: the charge projections G_{a,k} of the one-field genus F_{a,b} against −ψ_{2a+k} ι_{a+k}, exact series in ℚ(ω) to q², searching GL₂(F₃) for the relabeling (finds two, among them the one used) | < 1 s |
| `spin_structure_check.py` | Proposition C: the weights of the IR orbifold and of the code double agree in all four spin structures (tetracode: all 390,625 points of a box; Golay: 20,000 random points); a constant extra R phase changes the projection; a bilinear phase drops out | about 30 s |
| `lg_count_dense.py` | Corollary A, (1/\|C\|) Σ (−2)^{#common zeros} = Index(C)², for the Pless code P24 and the extended QR code of length 24 (2,304 = 48²) | 1–2 s |
| `pointwise_identity.py` | Theorem A and the charge-resolved q⁰ term of the refined identity on 23 codes (Golay, P24, QR24, 20 codes from random neighbour walks); importing `lg_identity_random.py` also runs its 30 walk checks of Corollary A | about 5 min |

These are copies of scripts of the companion project GappedCFT (`code/`); only pointers to its working note were replaced
by pointers to the public note (by `tools/package_lg.py` in the workspace). Licence: MIT (see `../LICENSE`).
