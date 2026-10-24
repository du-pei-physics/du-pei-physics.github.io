# The index of self-dual ternary codes: a proof of a conjecture of Gaiotto and Johnson-Freyd

`index.html` is a standalone note with complete proofs of both halves of GJF Conjecture 1.4 (arXiv:1811.00589).
Theorem A: the index of every self-dual ternary code of length 12k is divisible by 24.
Theorem B: an explicit code of index (-1)^k 24 exists at every length 12k.

Version of 23 October 2026, cited by D. Pei,
*Construction of holomorphic SCFTs with minimal index* (arXiv preprint, 2026). Not refereed. Authors: GPT-6 Astra, Fable 5.1 and Opus 5.5 (AI
agents), supervised by Du Pei.

| file | what it is |
|---|---|
| `index.html` | the note (MathJax 3 from cdnjs; otherwise self-contained) |
| `code/check_golay_chain_certificate.py` | plain Python. Checks Proposition 9 (u = T s, T u = -u, s . u = -24) from the certificate alone |
| `code/golay_chain_transfer.json` | the certificate: the nonzero entries of the 81 x 81 matrix T, the vector u = T s, and the vector s |
| `code/golay_chain_check.py` | numpy. Rebuilds T from its definition, checks Section 3, and writes the certificate (`--certificate`) |
| `code/index24_check.py` | numpy. Checks the lemmas of Section 2 lemma by lemma |
| `code/*_output.txt` | our runs (Python 3.12.15, numpy 2.5.4; certificate check with Python 3.9.6) |
| `code/LICENSE` | MIT, for the code |
| `figures/` | diagrams of the matrix A of the Golay code and of the chains C_k (PNG, exact JSON; Codex team), see `figures/README.md`; CC BY 4.0 |

```
cd code
python3 check_golay_chain_certificate.py
python3 golay_chain_check.py
python3 index24_check.py
```

The scripts are copies of those in `../holomorphic-svoa/code/index24/`, with section references renumbered for this note.
The certificate's data are identical to `../holomorphic-svoa/certificates/golay_chain_transfer.json`; only its
description line differs.
