# Regular complex Hadamard matrices of order 7,8

**GPT and Teo**, November 2026

A complex Hadamard matrix is called *regular* when each scalar product between distinct rows decomposes into rotated cycles of roots of unity. This paper reviews the general theory of regular complex Hadamard matrices, including the known classification at order 6, and studies the classification at orders 7 and 8.

## Main results

- **Order 7:** Classification of regular complex Hadamard matrices, up to equivalence, in terms of the Fourier matrix and the families \(P_7^q\) and \(R_7^q\).
- **Order 8:** Classification results organized by the possible prime-cycle decompositions of row quotients, including continuous families and exceptional matrices.
- **Deformations:** Results on affine deformations of Butson matrices at orders 7 and 8.
- **Exact computations:** Enumeration and verification procedures for row configurations, completions, equivalences, and exceptional cases.

## Contents

- `002.pdf` — the paper (64 pages).
- `002.tex` — LaTeX source.
- `code002/` — Python and C++ programs, computational inputs, certificates, and verification procedures.
- `code002/README.md` — detailed code documentation and reproduction instructions.

## Running the code

Python 3 is required; some computations use SymPy, and some order-8 routines use C++.

From the `code002` directory:

```
python3 -m venv ../code002_env
source ../code002_env/bin/activate
python -m pip install -r requirements.txt
python run.py list
python run.py check
```

The saved archives `N7_data.zip` and `N8_data.zip` contain computational inputs and results. See `code002/README.md` and `code002/SOURCE_INDEX.md` for further details. Full reproduction of some enumerations can require substantial computation.
