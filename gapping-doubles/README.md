# Gapping the non-chiral doubles of holomorphic SCFTs

`index.html` is a research note on the supersymmetric gapping problem for the non-chiral doubles Vbar (x) V of
holomorphic N = 1 SVOAs, with the index-1 theory P at c = 288 as the target: can one show that P is a periodic element by gapping
Pbar (x) P supersymmetrically? It covers:

- what is established;
- where the obstruction sits;
- the relation to the 576-fold periodicity of the spectrum SQFT;
- the open lanes.

Every statement carries a status label. The note accompanies the companion notes `../holomorphic-svoa/` and
`../gj-conjecture/` but is not cited by D. Pei, *Construction of holomorphic SCFTs with minimal index*,
[arXiv:2610.33165](https://arxiv.org/abs/2610.33165).

Version of 28 October 2026 (first version: 26 October 2026; see the update note in the page). Authors: Fable 5.1, GPT-6 Astra and Opus 5.5 (AI agents), supervised by Du Pei.

| file | what it is |
|---|---|
| `index.html` | the note (MathJax 3 from cdnjs; otherwise self-contained) |
| `code/tci_kw_parity.py` | sympy. Kramers-Wannier parities of the tricritical Ising primaries from the modular S-matrix |
| `code/t_poly_check.py` | plain Python. NS characters at c = 24 and of the k = 6 N = 2 intermediate as polynomials in t = (theta_3/eta)^12 |
| `code/mirror_tqft.py` | numpy. Vacuum bookkeeping of the Landau-Ginzburg Morse phase after the code orbifold (tetracode, Golay code) |
| `code/u13_duality.py` | numpy. Neveu-Schwarz Z3 self-duality of the c = 1 N = 2 minimal model (the IR of Phi^3): exact over Z[omega], and numerically from the lattice sums |
| `code/*_output.txt` | our runs |
| `code/LICENSE` | MIT, for the code |

```
cd code
python3 tci_kw_parity.py
python3 t_poly_check.py
python3 mirror_tqft.py
python3 u13_duality.py
```

Text: CC BY 4.0.
