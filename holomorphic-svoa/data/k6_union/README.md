# All 8,748 k = 6 codes: union dataset (details.html §C.2, §G.3)

The paper's Table 1 counts 8,748 distinct k = 6 codes (distinct systematic-generator rowspaces, not inequivalent codes or
theories). This folder extends the frozen census dataset `../k6_seed_census/` (7,112 codes, unchanged) by the other
1,636: 1,535 benchmark-v3 hits (lane `bench3_k6`, generators rebuilt from the recorded random pools and accepted only if
the rowspace hash matches the recorded one) and 101 codes from the stored records of the earlier benchmark lanes
(`bench2_k6`: 100, `bench_k6`: 1). None of them is a census code. `tools/package_k6_union.py` (workspace) wrote this
folder; `FILES.json` lists sources and SHA-256.

| file | what it is |
|---|---|
| `k6_union_polynomials.csv` | **one row per code, 8,748 rows**: `source` (hp, refl = census; bench3_k6, bench2_k6, bench_k6), `closure` (nega = anti-periodic chain, symmetry g = signed shift; cyc = periodic chain, g = −shift; both g³ = −1), `census_key` and `census_class_id` (census rows only: the class ids classify the 7,112 census seeds, not the union), `rowspace_sha256`, `sign_I`, `dim_T1`, `D1`, and `a6`, …, `a0`: Z_T = Tr_NS q^(L₀−3) = Σ a_j t^j, t = (ϑ₃/η)¹² (conventions in the comment header and in `../k6_seed_census/README.md`) |
| `k6_others_generators.jsonl.gz` | the 1,636 added generators (I₃₆ \| A), rows over {0,1,2}, with `sha`, `source`, `closure` |
| `k6_others_records.jsonl.gz` | their per-code records from `code/census_general/k6_general.py`: per-code hypotheses (fixed subcode of dimension 12, C̄ self-dual, balanced fixed words, no words of weight 1, 2, 4, 5 in C̄, \|index\| = 24), N, F, I, weight counts A₆, A₉, A₁₂, A₁₅, A₆₉, dimensions T₁ … T₅/₂, polynomial `a` (a₀ first), D₁ and its check D₁ = −4096 a₁ |
| `k6_union_summary.json` | counts: 4,579 distinct polynomials (2,961 among the census codes, 1,618 new, none shared), 4,563 distinct polar parts (NS spectra through weight 5/2), the 16 pairs differing only in a₀, 4,579 distinct (n, D₁), distinct values of a₀ … a₆, 2,701 codes with dim T₁ > 0 (dim T₁ up to 16; not every value in between occurs), the three values of X(C̄) |

Checks (all pass; asserted by the packaging tool or by the scripts): every union hash occurs once; each added generator is
systematic and its hash is the recorded one; for every added code `k6_general.analyse` verifies self-duality, that g
preserves C with g³ = −1, the hypotheses above, D₁ = −4096 a₁ and non-negativity of the derived dimensions.
`code/census_general/run_k6_union.py` recomputes any of them from the generators (identical records for samples of all
three sources and both closures, from a clean copy of this folder), and `validate_on_census.py` reproduces census
polynomials with the general engine. The census subset keeps its own statistics in `../k6_seed_census/`.

Licence: CC BY 4.0, like the rest of `data/` (see `../../LICENSE`).
