Regular complex Hadamard matrices of order 7,8

GPT and Teo, November 2026

A complex Hadamard matrix is called regular when the scalar products between its rows decompose into trivially vanishing sums of roots of unity, rotated by scalars. This paper reviews the general theory of regular complex Hadamard matrices, including the known classification at order 6, and investigates the classification at orders 7 and 8.

The main results are:

General theory: every regular complex Hadamard matrix is an affine deformation of a Butson matrix.

Order 7: up to equivalence, the regular matrices are the Fourier matrix and the 2 families (P_7^q) and (R_7^q), related by conjugation and monomial normalization.

Order 8: classification into 9 continuous families and 4 isolated matrices, organized by the possible prime-cycle decompositions of row quotients.

Computational classification: exact enumeration and verification procedures for the relevant row configurations, completions, equivalences and exceptional cases.

Contents

002.pdf — the paper (64 pages).

002.tex — its LaTeX source, compiled with pdflatex.

code002/ — Python and C++ programs, input data, saved certificates and verification procedures supporting the classification results.

code002/README.md — detailed documentation and instructions for reproducing the computations.

Running the code

Python 3.10 or later is required. Some order-8 checks use SymPy; the more extensive enumeration programs also use C++17.

From the code002 directory:

python3 -m venv ../code002_env
source ../code002_env/bin/activate
python -m pip install -r requirements.txt
python run.py list
python run.py check

The default verification procedure runs 16 named checks against the saved certificates. Full enumeration and reproduction of the classification are separate procedures, described in code002/README.md, and some can require substantial computation.

The archives N7_data.zip and N8_data.zip contain the retained computational inputs and results. The file SOURCE_INDEX.md maps the original source paths to their reorganized locations.
