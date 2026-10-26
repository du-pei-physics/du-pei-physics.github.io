# GJF's Conjecture 1.4 (index 24 for self-dual ternary codes of length 12k): numerical checks

The proofs are in `details.html`. Section A.7 (`../../details.html#A7`): 8 | Index(C) for every self-dual ternary code
of length divisible by 12, by an elementary neighbour argument, and with GJF's Theorem 1.3 (12 | Index) the
divisibility half of GJF's Conjecture 1.4, 24 | Index(C). Section A.8 (`../../details.html#A8`): an explicit chain C_k
of k Golay codes with index (-1)^k 24 at length 12k, the existence half. Everything in this folder is a numerical check
of these proofs; Python 3 and numpy are enough.

| file | what it is |
|---|---|
| `index24_check.py` | stand-alone check of every lemma of A.7; it imports nothing else from this folder. `python3 index24_check.py` runs the default suite (n = 12, 24, 36; about 30 s); our n = 48 run was `python3 index24_check.py --n 48 --steps 6 --chains 0 --seed 2` (about 75 s) |
| `index24_check_output.json`, `index24_check_output_n48.json` | our two runs (Python 3.12.15, numpy 2.5.4): 4,022 and 1,217 checks, no failure; the indices met (all divisible by 24) and the values of Index(D) (6, 18, 30, ... occur, which is why Lemma 3 is needed) |
| `referee/` | the scripts of the independent referee, a separate agent that wrote its own F_3 code (`f3.py`) and tested each lemma: `t1_L1.py` (Lemmas 1 and 3 on 1,200 random codes, n = 3-16), `t2_L2.py n` (Lemma 2 on neighbour pairs, with a brute-force split of D^perp into its nine cosets), `t3_random.py n count` (8 \| Index on random self-dual codes), `t4_chain.py n count` (explicit Lemma-4 chains); copied unchanged. `referee_rerun.txt` is a re-run of all of them for this page (13 runs, all pass) |
| `golay_chain_check.py` | stand-alone check of A.8, written from its statement (imports nothing else from this folder): the Golay code, the finite facts used, the 81 x 81 transfer matrix T built from its definition, T u = -u, s^T T^k s = (-1)^k 24 and the word count 3^(6k) for k <= 30 in exact integers, and C_1, ..., C_4 built from the closed form with their indices by enumeration (-24, 24, -24, 24). `python3 golay_chain_check.py`, about 7 s; with `--certificate ../../certificates/golay_chain_transfer.json` it also writes the sparse certificate |
| `golay_chain_check_output.json` | our run: 16 checks, no failure; T has 423 non-zero entries, u = T s has 33 (listed) |
| `check_golay_chain_certificate.py` | plain Python (no numpy): reads `../../certificates/golay_chain_transfer.json` (state ordering, the non-zero entries of T and u, the vector s) and checks u = T s, T u = -u, s . u = -24 and s^T T^k s = (-1)^k 24 for k <= 20 from that file alone |
| `walk_records/` | the workspace test that first checked the proof of A.7: `neighbor_index_test.py` (it uses `../../verify.py`; only its import path was changed) and its records `neighbor_index_n*.json`, 492 neighbour steps at n = 12, 24, 36, 48 with no failure; the seeds were not recorded, so other seeds give other walks |

For A.7 three independent implementations of the index are involved: `index24_check.py` (enumeration over an information
set, cross-checked against plain enumeration), the referee's `f3.index_fast` / `f3.index_brute`, and `verify.py`'s
`brute_index` (used by `walk_records/neighbor_index_test.py`).

Citable standalone version of both proofs, with copies of these scripts and the certificate: `../../../gj-conjecture/`
(to be published at https://du-pei-physics.github.io/gj-conjecture/); the paper cites that note.

Status: both halves proved (A.7, A.8), checked numerically here, by the referee (A.7) and by independent
recomputations (A.8). Licence: MIT (see `../LICENSE`).
