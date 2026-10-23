# The k = 6 seed census: codes and invariants of the seed theories (details.html §C.2, §G.3)

For each of the 7,112 k = 6 seed cells of the two censuses these scripts build the length-72 code (generator (I₃₆ | A),
the anti-periodic chain of three copies of the cell) and compute invariants of the seed theory T = V_{Λ(C)}//Z₆ that
separate theories (details.html §G.3), up to the complete NS and Ramond characters (a polynomial in t per seed).
Inputs: `../../data/k6_seed_census/` (the census hit files); outputs: `../out/census/`. They use `hpcell.py`, `f3lib.py`,
`mixed_rings.py`, `verify_mixed_hit.py` and `fc_closed_form.py` from `code/`. Exact
integer and F₃ arithmetic; one process each. Run from `code/` with the same Python as the rest of the folder
(numpy, sympy).

| script | what it does | time |
|---|---|---|
| `count_k6_seed_codes.py` | builds all 7,112 codes, checks self-duality, counts distinct rowspace hashes (7,112), checks the construction against the k = 6 certificate (same hash) and the index on a random sample (banded transfer DP over the generator) | 6 s with `--sample 5`; the default 60 takes longer |
| `fc_seed_census.py` | Ramond ground states of T: n₊ + n₋ = (N + 8F)/6, n₊ − n₋ = I/6, with N, F, I computed twice (from the cell's transfer tables and from the length-72 generator) and the hypotheses checked per seed (C̄ self-dual of dimension 12, balanced fixed words); the symmetry class of the cell (canonical form under PGL(2,11) and row monomials); `--limit n`, resumable; `--summarize` writes the summary | about 0.8 s per seed (1.6 CPU-hours for all) |
| `fc_seed_level1.py` | states at the first excited Ramond level, D₁ = (144N + 8A₆₉)/6 + X(C̄) + 864,288,768, and dim T₁ = A₆/6, from truncated transfer DPs over the generator and the weight enumerator of C̄; `--part i/m` splits the seeds into m disjoint parts (we ran 2), `--summarize` joins them with the output of `fc_seed_census.py` | about 0.4 s per seed |
| `fc_seed_level2.py` | the complete NS and Ramond characters of T: A₆, A₉, A₁₂, A₁₅ from a weight-graded transfer trace over the generator give dim T_h for h ≤ 5/2, which together with n fix the degree-6 polynomial in t (details.html §G.3); asserts −4096 a₁ = D₁ of `fc_seed_level1.py` for every seed; `--part i/m`, `--summarize` (needs the outputs of the two scripts above) | about 0.9 s per seed |
| `fc_seed_polynomials.py` | writes `out/census/k6_seed_polynomials.csv` (the polynomial of every seed) and its summary, after per-seed checks: the q-expansion reproduces dim T_h for h ≤ 5/2, a₀ = n, −4096 a₁ = D₁; the paper's seed gives its polynomial; and distinct polynomials give distinct NS dimensions of T″ up to weight 3 | 12 s |
| `nice_reps.py`, `symmetry.py` | the 1,320 monomial automorphisms of the Paley matrix, one per element of PGL(2,11), used for the symmetry classes; importing `nice_reps.py` prints its demonstration (the PGL(2,11)-orbits of three seeds), as in our runs | library |

The docstrings of `fc_seed_census.py` and `fc_seed_level1.py` give the sector-by-sector derivation of the two formulas
(given the inputs A1–A2 of details.html §H.1).

Examples (a few seconds each):

```
python3 census/count_k6_seed_codes.py --sample 5
python3 census/fc_seed_census.py --limit 3                # the first three Hadamard–Paley seeds
python3 census/fc_seed_level1.py --part 1/2 --limit 2
```

Each record written should equal the corresponding record in `../../data/k6_seed_census/`: from this folder we reran
three Hadamard–Paley and two reflection seeds with `fc_seed_census.py` and two of each with `fc_seed_level1.py`, and
every record came out identical; likewise `fc_seed_level2.py` on two seeds of each kind. To rerun the summaries on our
full records, unpack `../../data/k6_seed_census/*.jsonl.gz` into `out/census/` and run `python3 census/fc_seed_census.py
--summarize`, `python3 census/fc_seed_level1.py --summarize`, `python3 census/fc_seed_level2.py --summarize` and `python3
census/fc_seed_polynomials.py`: the three summaries and `k6_seed_polynomials.csv` come out byte for byte as the ones in the
data folder (the polynomial summary differs only in the path it records for the CSV).

These are copies of the workspace scripts (`claude/paper/fc_seed_census.py`, `fc_seed_level1.py`, `fc_seed_level2.py`,
`fc_seed_polynomials.py`, `claude/census/count_k6_seed_codes.py`, `nice_reps.py`, `symmetry.py`); only import, input and output paths were
changed (by `tools/build_code_folder.py` in the workspace). Licence: MIT (see `../LICENSE`).
