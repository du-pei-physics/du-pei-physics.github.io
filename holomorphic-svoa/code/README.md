# Code for the holomorphic-SVOA note

Everything here runs from this folder with Python 3. `verify.py` needs only **numpy**; the reproduction scripts also
use **sympy** (and `fc_modular_check.py` **mpmath**). Inputs are read from `../certificates/original/`; outputs go to
`out/`. All arithmetic is exact (integers, F_3, fractions, Z[ω]).

Tested with: Python 3.12.15, numpy 2.5.4, sympy 1.14.0, mpmath 1.3.0 on macOS (Darwin, arm64), single core.
Licence: MIT (see `LICENSE`).

## Verify the certificates

```
python3 verify.py --selftest     # banded transfer DP vs brute force on random small self-dual chains
python3 verify.py                # every ../certificates/k*.json, about 30 s; prints one row per code
python3 verify.py --json report.json ../certificates/k21_codex.json   # one file, full report
```

`verify.py` is self-contained and independent of the search code below: it rebuilds each code from its generator
matrix (R in the paper), checks that its rowspace hash matches the certificate (distinct hashes = distinct codes),
checks self-duality, computes the index and the number of full-support codewords (Ramond ground states) by a
banded transfer DP over the generator alone (the cell widths stored in the certificate are only a hint; the zero
pattern they imply is checked), checks the stated symmetry (preservation, order, cycle type, `g^(m/2) = -1`,
3-cycles of `h = g^(m/3)` and the balance of its fixed full-support words via the contracted fixed subcode; for Q8 the
quaternion relations, group order, cyclic commuting pairs) and prints the orbifold index next to the TMF minimum
24/gcd(24,k). The report separates `finite_checks_pass` (exact facts about the code and its symmetry) from
`orbifold_index_status`: GJF Thm 5.6 (cited, its finite hypotheses checked) for cyclic groups,
`orbifold_index_conditional_on_anomaly_lemma` for Q8 (the non-anomaly argument of details.html §F.3, which no finite check
replaces), or no orbifold (m = 1). The modular int64 matrix products assert inner dimension × (p−1)² < 2⁶³ and reject
inputs that would overflow. It ends with the index arithmetic of the rows without a code and an independent sum over the 40
commuting pairs of D4 for the twisted D4 orbifold of four copies of the k = 6 theory (the paper's figure). Our output
is `../verify_output.txt`. Notation follows the paper: chains of cells (periodic / anti-periodic), boundary spins, transfer
matrix M; a *candidate* is one random choice of the cells of a chain. The code itself still says "ring", "cyc"/"nega"
and "delayed" in places (older names for chain, periodic/anti-periodic and boundary).

## Reproduce the constructions

| script | what it does | time |
|---|---|---|
| `canonical_seed_data.py` | rebuilds the k = 6 seed of the paper (Hadamard–Paley cell, R = {3,5,8}) with its checks; writes `out/canonical_seed.json` (identical generator to the certificate) and the LaTeX matrices | 1 s |
| `hpcell.py` | fast exact transfer matrices for Hadamard–Paley cells with any set of boundary coordinates (the census engine) | library |
| `cells_general.py` | cell bases (Paley, conference, Hadamard blocks, random anti-orthogonal), transfer matrices, traces | library |
| `mixed_products.py` | random chains of cells of mixed widths, one or several periods (k = 13, 15 searches) | `--minutes` |
| `mixed_rings.py`, `worker_common.py` | the checkpointed period-chain worker used for k = 15, 18, 21 and the benchmark (3 or 4 boundary spins, discriminant 2, fast modular screen, inline verification of hits) | `--max-minutes` |
| `verify_mixed_hit.py`, `verify_period_hits.py` | rebuild and verify hits of the two searches above | seconds |
| `q8_cells.py` | Q8-equivariant cells and chains (k = 8, 16); `--selftest` | 1 s |
| `neighbour_chain.py` | Z^96 = Λ(tetracode^24) with a Q8-stable frame, and the 23-step Q8-equivariant neighbour chain from the k = 8 code | 1 s |
| `psi_cocycle_check.py` | the SL(2,3)-invariance of the Z_3^2 cocycles used in the gauging lemma | < 1 s |
| `perm_orbifolds.py`, `check_wreath.py` | index polynomials of all permutation orbifolds of degree ≤ 6 with all spin twists; the D4 twist table | 40 s, 1 s |
| `fc_weight_enum.py` → `fc_characters.py` → `fc_theta.py`, `fc_sym2.py`, `fc_closed_form.py`, `fc_modular_check.py` | field content of T, T′, T″ (called P in their output; run in this order) | ≤ 5 s each |
| `power_operations_check.py [M]`, `power_operations_d4.py` | exact q-series checks of the TMF power-operation computation (`M` = precision, default 40 ≈ 2 min; 30 takes 25 s; very small values such as 12 fail) | |
| `time_candidates.py` | measured CPU cost of building one cell's transfer matrix and of evaluating one candidate exactly (3 and 4 boundary spins); our output is `../notes/timings.json` | 3 s |
| `f3lib.py`, `index_gray.c` | F_3 linear algebra used by the search code; `index_gray.c` is compiled on demand by `f3lib.index_exhaustive` (needs a C compiler) | |

## GJF's Conjecture 1.4 (details.html §A.7, §A.8)

`index24/` checks the two proofs: that 24 | Index(C) for every self-dual ternary code of length divisible by 12
(`python3 index24/index24_check.py`, numpy only, stand-alone, about 30 s, tests each lemma of §A.7 separately;
`index24/referee/` holds the independent referee's scripts with a re-run, `index24/walk_records/` the first workspace
test), and that the Golay chain codes of §A.8 have index (-1)^k 24 at every length 12k (`python3
index24/golay_chain_check.py`, numpy only, stand-alone, about 7 s; `index24/check_golay_chain_certificate.py` checks the
sparse certificate `../certificates/golay_chain_transfer.json` of §A.8.5 from the file alone). See `index24/README.md`.

## The k = 6 seed census (details.html §C.2, §G.3)

`census/` builds the codes of the 7,112 k = 6 seed cells of the two censuses (data: `../data/k6_seed_census/`) and
computes invariants of the seed theories T: the number of Ramond ground states (`census/fc_seed_census.py`) and the
number of states at the first excited Ramond level (`census/fc_seed_level1.py`), which together take 2,961 values, and
the complete NS character as a polynomial in t (`census/fc_seed_level2.py`, `census/fc_seed_polynomials.py`; 2,961
distinct polynomials); `census/count_k6_seed_codes.py` checks that the 7,112 codes are distinct. A few seconds on a
handful of seeds; the full runs took about 1.6, 0.8 and 1.8 CPU-hours. Outputs go to `out/census/`. See `census/README.md`.

## All 8,748 k = 6 codes (details.html §G.3)

`census_general/` computes the same characters for any k = 6 generator of a chain of three width-12 cells (anti-periodic
or periodic closure) and recomputes the 1,636 non-census codes of `../data/k6_union/` from their shipped generators. With
the census, the 8,748 codes give 4,579 distinct character polynomials (4,563 already through weight 5/2). See
`census_general/README.md`.

## Landau–Ginzburg orbifolds and refined elliptic genera (details.html §L)

`lg/` checks the identities of the paper's §6.3–6.4 and Appendix A (proofs in `../notes/lg_elliptic_genus.md`): the
quintuple-product form, the per-field refined genus, the index identity and the four spin-structure comparison. Seconds
each, except `lg/pointwise_identity.py` (about 5 min). See `lg/README.md`.

## Benchmark-v3 worker sources

`bench_v3/` archives the exact worker sources behind the benchmark table of the summary page (`bench_worker.py`,
`bench_q8.py`, `bench_q8_blocks.py`, `run.py` and their imports, in their workspace layout), recovered by the SHA-256 the
runs recorded, together with each lane's run records (seed, argv, versions). See `bench_v3/README.md` and
`bench_v3/SOURCES.json`; `mixed_rings.py` above is the older, path-patched worker kept for the examples below.

Examples (a few seconds each):

```
python3 mixed_products.py --k 13 --widths 12,12,12,12,12,12,6 --minutes 0.1 --seed 1 --out out/k13_probe.json
python3 mixed_rings.py --name probe_k21 --k 21 --widths 12,12,12,6 --r 4 --metric2 1 --s-refl 2 --exclude-bases A0 --fast 1 --max-minutes 0.1 --seed 7
python3 q8_cells.py --selftest
python3 neighbour_chain.py
```

Hits are rare (see the cost table on the summary page): short runs will normally find none. These scripts are copies
of the project's working code; only import and output paths were changed (by `tools/build_code_folder.py` in the
workspace).
