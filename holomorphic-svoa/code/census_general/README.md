# Characters of arbitrary k = 6 generators (details.html §G.3; data in ../../data/k6_union/)

| script | what it does | time |
|---|---|---|
| `k6_general.py` | `analyse(G, closure)`: for a generator (I₃₆ \| A) of a chain of three width-12 cells, anti-periodic (`nega`, g = signed shift) or periodic (`cyc`, g = −shift), checks self-duality, that g preserves C and g³ = −1, the hypotheses on the g²-fixed code C̄ (dimension 12, self-dual, balanced full-support words, no words of weight 1, 2, 4, 5), computes N, F, I, the weight counts A₆ … A₁₅, A₆₉ by transfer traces over the generator, the dimensions T₁ … T₅/₂, the character polynomial in t and D₁ from the twisted sectors, and checks D₁ = −4096 a₁. Both closures give a g of the same type (24 signed 3-cycles with sign product −1, g³ = −1), which is all the formulas of `../census/` use | about 2 s per code |
| `run_k6_union.py` | recomputes codes of `../../data/k6_union/` from the shipped generators and compares with the shipped records; `--limit`, `--source`, `--closure`; output `../out/census_general/` (written for this folder) | seconds |
| `validate_on_census.py [seed] [count]` | runs the general engine on random census seeds and compares with `../../data/k6_seed_census/k6_seed_polynomials.csv` | seconds |

It imports `f3lib.py`, `mixed_rings.py`, `verify.py`, `fc_closed_form.py` from `code/` and `fc_seed_level1.py`,
`fc_seed_level2.py` from `code/census/`; numpy only. Run from this folder, e.g. `python3 -B run_k6_union.py --limit 2
--closure cyc`. The full run over the 1,636 codes took about 3,100 CPU s (one core). The generators of the 1,535
benchmark-v3 hits were rebuilt in the workspace from the recorded pools (the worker stored full records only for its
first 100 hits); that replay script needs the private run state and is not included, so the rebuilt generators are
shipped instead, each checked against its recorded hash.

`k6_general.py` and `validate_on_census.py` are copies of workspace scripts (`claude/census_general/`); only import
paths were changed (by `tools/build_code_folder.py` in the workspace). Licence: MIT (see `../LICENSE`).
