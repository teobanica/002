#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
mkdir -p checks outputs
"$PYTHON" code/parameterize.py
"$PYTHON" code/prepare_components.py
compile() {
 g++ -O3 -std=c++17 "code/$1.cpp" -o "checks/$1.pending"
 test -s "checks/$1.pending"
 mv "checks/$1.pending" "checks/$1"
}
compile lattice_check
"$PYTHON" code/run_exact.py --records 1200 outputs/lattice_tests.json -- checks/lattice_check
"$PYTHON" code/verify_lattice.py
compile events
"$PYTHON" code/run_exact.py --records 2 --log checks/events.log outputs/A_events.json -- checks/events inputs/components.txt
"$PYTHON" code/planes.py
compile aaa_scan
"$PYTHON" code/run_exact.py --records 42 --log checks/planes_scan.log outputs/planes_scan.json -- checks/aaa_scan inputs/planes.txt "$(cat outputs/planes_modulus.txt)"
"$PYTHON" code/prepare_good.py
for cid in 11 19; do
 compile full
 "$PYTHON" code/run_exact.py --log "checks/full${cid}.log" "outputs/full${cid}.json" -- checks/full inputs/components.txt "$cid" "inputs/good${cid}.txt"
done
"$PYTHON" code/prepare_pairs.py planes
compile pairs
"$PYTHON" code/run_exact.py --records 26 outputs/planes_pairs.json -- checks/pairs inputs/planes_pairs.txt "$(cat outputs/planes_pairs_modulus.txt)"
"$PYTHON" code/strata.py
"$PYTHON" code/prepare_pairs.py lines
compile pairs
"$PYTHON" code/run_exact.py --records 294 outputs/lines_pairs.json -- checks/pairs inputs/lines_pairs.txt "$(cat outputs/lines_pairs_modulus.txt)"
"$PYTHON" code/final_points.py
"$PYTHON" code/finite_cases.py
compile finite_scan
"$PYTHON" code/run_exact.py --records 132 --log checks/finite_scan.log outputs/finite_scan.json -- checks/finite_scan inputs/finite_cases.txt "$(cat outputs/finite_modulus.txt)"
"$PYTHON" code/verify_final.py
printf '%s\n' 'COMPLETE: exact completion-existence criteria for AAA classes 11 and 19.'
