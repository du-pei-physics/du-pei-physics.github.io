# Changes

## Version of 23 October 2026 (revised paper: Claim 1, §6.3–6.4, Appendix A)

Added
- `data/k6_union/`: all 8,748 k = 6 codes of the paper's Table 1 (the 7,112 census codes plus 1,636 benchmark codes, with
  the generators of the added ones) and the character polynomial of each; 4,579 distinct polynomials, 4,563 already
  through weight 5/2. `code/census_general/` recomputes them. The census dataset `data/k6_seed_census/` is unchanged.
- details.html §L and `notes/lg_elliptic_genus.md`: the Landau–Ginzburg material of the paper (Claim 1 and its inputs,
  the gapping test, the index identity, the refined elliptic genus, the four spin-structure comparison), with checks in
  `code/lg/`.
- details.html §E.4: the precise scope of the extended permutation-orbifold search (family, finished k ≤ 23 closure,
  runs still in progress, exactness of the screening, what is not covered).

Changed
- G.3 and the overview: inequivalence counts for all 8,748 codes (at least 4,579 theories; 4,563 from the low-weight NS
  spectra); the census subset keeps its counts (2,961, 2,958).
- H.1 and `notes/field_content.md`: the Ramond part of input A1 is now cited from GJF (odd-order lifts, proof of Thm 5.6)
  instead of being called an assumption; A2 is the standard anomaly-free orbifold construction.
- `notes/field_content.md`: a stray q^j in a formula of §6.1 removed; §6.1 restricted to 24 | c (with the sign for
  12 | c); the summary now states the minimum-distance argument for the absence of weight-1 fields, the oscillator states
  below weight 3, and the new untwisted states of T′ and P at weight 3; script paths.
- Smaller scope statements (k = 15 certificate, signs in the k = 12, 24 rows), datelines, references.
- Release hygiene: bytecode caches removed; the manifest build refuses caches and unexpected file types.

Unchanged: the canonical construction (T, T′, T″ and their numbers), the 17 certificates and `verify.py`, the frozen
benchmark-v3 table and its audit bundle, and all earlier proofs.
