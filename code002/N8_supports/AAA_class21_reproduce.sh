#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
mkdir -p checks outputs
"$PYTHON" code/parameterize.py
g++ -O2 -std=c++17 code/archived_canonical_aaa.cpp -o checks/canonical
"$PYTHON" code/run_exact.py --log checks/canonical.log outputs/canonical_rerun.json -- checks/canonical
"$PYTHON" - <<'PY'
import json
a=json.load(open('inputs/canonical.json'))['triples']['AAA']
b=json.load(open('outputs/canonical_rerun.json'))['triples']['AAA']
assert a==b
print('PASS: archived AAA canonical census reproduced exactly.')
PY
"$PYTHON" code/archived_square_verify.py
g++ -O3 -std=c++17 code/aaa_scan.cpp -o checks/aaa_scan
"$PYTHON" code/run_exact.py --records 38 --log checks/generic_AAA.log outputs/generic_AAA.json -- checks/aaa_scan inputs/components.txt
"$PYTHON" code/verify_witnesses.py
g++ -O3 -std=c++17 code/class21_events.cpp -o checks/events
"$PYTHON" code/run_exact.py --records 2 --log checks/class21_A_events.log outputs/class21_A_events.json -- checks/events inputs/components.txt
"$PYTHON" code/class21_lines.py
"$PYTHON" code/run_exact.py --records 32 --log checks/class21_line_scan.log outputs/class21_line_scan.json -- checks/aaa_scan inputs/class21_lines.txt "$(cat outputs/line_modulus.txt)"
"$PYTHON" code/class21_A_residue.py
"$PYTHON" - <<'PY'
import json
d=json.load(open('outputs/class21_A_residue.json'))
with open('inputs/extra_points.txt','w') as f:
 for p in d[0]['off_good_points']:f.write(' '.join(str(v) for q in p for v in q)+'\n')
PY
g++ -O3 -std=c++17 code/class21_full.cpp -o checks/full
for branch in 0 1; do
 "$PYTHON" code/run_exact.py --log "checks/class21_full${branch}.log" "outputs/class21_full${branch}.json" -- checks/full inputs/components.txt "$branch" inputs/extra_points.txt
done
"$PYTHON" code/prepare_pairs.py
g++ -O3 -std=c++17 code/pair_events.cpp -o checks/pairs
"$PYTHON" code/run_exact.py --records 2 outputs/global_regular_pair_events.json -- checks/pairs inputs/pairs_generic.txt 30
"$PYTHON" code/run_exact.py --records 16 outputs/line_regular_pair_events.json -- checks/pairs inputs/pairs_lines.txt 60
"$PYTHON" code/final_residue.py
g++ -O3 -std=c++17 code/finite_scan.cpp -o checks/finite
for branch in 0 1; do
 "$PYTHON" code/run_exact.py --records 64 --log "checks/finite${branch}.log" "outputs/finite${branch}.json" -- checks/finite inputs/final_points.txt "$(cat outputs/final_point_modulus.txt)" "$branch"
done
"$PYTHON" code/verify_class21.py
"$PYTHON" code/export_templates.py
printf '%s\n' 'COMPLETE: AAA foundation audit and class-21 exact completion locus.'
