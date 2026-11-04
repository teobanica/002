# Code for the regular classification at orders 7 and 8

Start here. The mathematical programs are unchanged; the folders and supporting
files have been reorganized. This package accompanies `002.tex` and does
not contain another copy of the manuscript.

## What is where

| Folder | Contents | Source files |
| --- | --- | ---: |
| `N7` | Order 7: triples, fourth rows, completions, square families | 42 |
| `N8_AB` | Order 8: A/B completions and the BBB exclusion | 37 |
| `N8_AAA` | Generic AAA completions, family comparisons, catalogue checks | 34 |
| `N8_C` | C-pair extensions and the three level-30 classes | 11 |
| `N8_supports` | Exceptional AAA classes, exact solvers and supporting computations | 62 |
| `N8_checks` | Row and column cycle censuses | 1 |

All 187 original source files are directly inside these folders. Identical copies
have been consolidated; different solvers remain separate. `SOURCE_INDEX.md`
maps every old program path to its new location.

`N7_data.zip` and `N8_data.zip` contain the inputs and saved results, including
all 1,286 JSON files from the previous package. These are ordinary ZIP files:
you can inspect or extract them without special software. Their contents have
not been replaced by summaries. The two TeX matrix files read by the cycle
checker are also retained there.

## Checking the saved certificates

Use Python 3.10 or later. Order-7 checks use only the standard library. Some
order-8 checks need SymPy 1.14.0. From this folder:

```sh
python3 -m venv ../code002_env
source ../code002_env/bin/activate
python -m pip install -r requirements.txt
python run.py list
python run.py check
```

`check` runs the 16 named checks listed by `list`. These include direct family
identities, independent order-7 verifications, C-sector certificates, AAA family
equivalences and the final A/B completion checks. It checks existing enumeration
records; it does not repeat the entire exhaustive classification. Some final
checkers also use previously saved verification records. See `VALIDATION.json`
for the checks actually run during this cleanup.

To run just one check:

```sh
python run.py check --only n7-squares
python run.py check --only n8-ab-final
```

The launcher first checks the packaged files and creates a separate sibling
folder, `code002_work`, with the relative paths expected by the original programs.
Results and compiler outputs go there. The browsing folder `code002` stays small.
The execution folder necessarily has more files and the older directory layout;
this is what lets the mathematical source remain byte-for-byte unchanged.
It is not included in the publication ZIP.

To unpack without running calculations, or to choose a different work folder:

```sh
python run.py prepare
python run.py --work ../another_code002_work prepare
```

Use a new work-folder name for a fresh run. The launcher refuses to overwrite
an unrelated folder. Never use Python's `-O` option: these programs use assertions
as mathematical checks.

## Running other retained programs

Supply the original path, as shown in `SOURCE_INDEX.md`, to the launcher:

```sh
python run.py run N7/run_all.py --list
python run.py run N8/n8_AB_completion/code/reproduce.py --final-only
```

For the more extensive existing A/B certificate checks, omit `--final-only`.
Running the files directly from the shallow source folders will not in general
work: they still use their original imports and relative paths. After `prepare`,
you can also run them directly from the corresponding `code002_work` directories.
Paths cited in the paper refer to that restored layout.

The original full-reproduction entry points are retained:

| Original path, passed after `python run.py run` | Purpose |
| --- | --- |
| `N7/run_all.py` | All eight order-7 stages; `--stage` selects a stage |
| `N8/n8_aaa_complete/code/rebuild_generic.py` | Generic AAA enumeration |
| `N8/n8_C_sector/code/rebuild.py` | C-sector enumeration and checks |
| `N8/n8_AB_completion/code/reproduce.py --full` | Full A/B pipeline |
| `N8/proofs/AAA_class21/reproduce.sh` | Exceptional class 21 |
| `N8/proofs/AAA_class22/reproduce.sh` | Exceptional class 22 |
| `N8/proofs/AAA_classes11_19/reproduce.sh` | Exceptional classes 11 and 19 |
| `N8/proofs/consolidation/reproduce.sh` | Supporting triple classification |

Full searches can take many hours. The C++ programs require GNU C++ with C++17
support, available under the name `g++` used by the original scripts. Several
use GNU headers, so Apple's default Clang alias named `g++` is not sufficient.
This cleanup was checked with GNU g++ 13.3.0 on Linux; a complete Mac build has
not been tested. The default 16 certificate checks do not compile C++.

The removed `expanded_report/reproduce_report.sh` was an old PDF wrapper. Its
three mathematical sources remain under `N8_supports/expanded_report_*`.
Its seed check can still be run after preparing the work folder:

```sh
cd ../code002_work/N8/proofs/expanded_report
g++ -O3 -std=c++17 code/n8_seed_extend_bruteforce.cpp -o checks/seed_check
python code/verify_seed.py
```

## What changed in this cleanup

- Removed 12 obsolete document assembly, PDF-checking and packaging scripts.
- Consolidated 30 byte-identical source copies into their retained originals.
- Removed old logs, restart fingerprints, old manuscript TeX and a nested handoff.
- Preserved every original JSON record and all other retained input bytes.
- Restored the missing `N7/sources/BH-7-6.txt` from Teo's earlier uploaded file.
  The previous ZIP omitted this file even though the square checker needs it.
- Added this guide, the launcher, a source index and one manifest. The manifest
  records every retained path, checksum, duplicate mapping and removal reason.

The old document-building dependencies (`pypdf`, Pillow, PyMuPDF and LaTeX) are
no longer needed for this code package. Source attribution is in
`SOURCE_ATTRIBUTION.txt`. No new license for the authors' code is imposed here.
