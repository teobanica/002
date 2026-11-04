#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
mkdir -p checks outputs
"$PYTHON" code/parameterize.py
"$PYTHON" code/prepare_class22.py components
g++ -O3 -std=c++17 code/aaa_scan.cpp -o checks/aaa_scan
"$PYTHON" code/run_exact.py --records 2 --log checks/generic_AAA.log outputs/generic_AAA.json -- checks/aaa_scan inputs/class22_components.txt
g++ -O3 -std=c++17 code/class22_events.cpp -o checks/events
"$PYTHON" code/run_exact.py --records 2 --log checks/class22_A_events.log outputs/class22_A_events.json -- checks/events inputs/class22_components.txt
"$PYTHON" code/class22_lines.py
"$PYTHON" code/run_exact.py --records 24 --log checks/class22_line_scan.log outputs/class22_line_scan.json -- checks/aaa_scan inputs/class22_lines.txt "$(cat outputs/line_modulus.txt)"
"$PYTHON" code/class22_A_residue.py
"$PYTHON" code/prepare_class22.py strata
g++ -O3 -std=c++17 code/class22_full.cpp -o checks/full
for branch in 0 1; do
 "$PYTHON" code/run_exact.py --log "checks/class22_full${branch}.log" "outputs/class22_full${branch}.json" -- checks/full inputs/class22_components.txt "$branch" inputs/extra_points.txt inputs/good_lines.txt
done
"$PYTHON" code/prepare_pairs.py
g++ -O3 -std=c++17 code/pair_events.cpp -o checks/pairs
"$PYTHON" code/run_exact.py --records 2 outputs/global_regular_pair_events.json -- checks/pairs inputs/pairs_generic.txt 30
"$PYTHON" code/run_exact.py --records 6 outputs/line_regular_pair_events.json -- checks/pairs inputs/pairs_lines.txt "$(cat outputs/line_modulus.txt)"
"$PYTHON" code/final_residue.py
g++ -O3 -std=c++17 code/finite_scan.cpp -o checks/finite
for branch in 0 1; do
 "$PYTHON" code/run_exact.py --records 50 --log "checks/finite${branch}.log" "outputs/finite${branch}.json" -- checks/finite inputs/final_points.txt "$(cat outputs/final_point_modulus.txt)" "$branch"
done
"$PYTHON" code/verify_class22.py
"$PYTHON" code/export_templates.py
"$PYTHON" code/verify_displayed_data.py
printf '%s\n' 'COMPLETE: exact nine-circle class-22 completion theorem, both components.'
