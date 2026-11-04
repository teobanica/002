#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
mkdir -p checks outputs report
compile() {
 g++ -O3 -std=c++17 "code/$1.cpp" -o "checks/$1.pending"
 test -s "checks/$1.pending"
 mv "checks/$1.pending" "checks/$1"
}
compile n8_pilot
"$PYTHON" code/run_exact.py --log checks/pilot.log outputs/pilot.json -- checks/n8_pilot
compile n8_canonical
"$PYTHON" code/run_exact.py --log checks/canonical.log outputs/canonical.json -- checks/n8_canonical
"$PYTHON" code/affine.py > checks/affine.log
"$PYTHON" code/pure_normal_forms.py
for name in BBB CCC; do
 if [ "$name" = BBB ]; then cid=100; else cid=101; fi
 compile pure_extensions
 "$PYTHON" code/run_exact.py --log "checks/${name}_extensions.log" "outputs/${name}_extensions.json" -- checks/pure_extensions inputs/pure_components.txt "$cid" inputs/no_good_loci.txt
done
compile canonical_squares
"$PYTHON" code/run_exact.py --records 60 --log checks/ccc_canonical.log outputs/ccc_canonical.json -- checks/canonical_squares inputs/ccc_squares.txt
"$PYTHON" code/verify_pure.py
"$PYTHON" code/finite_cases.py
compile finite_scan
"$PYTHON" code/run_exact.py --records 14 --log checks/finite_scan.log outputs/finite_scan.json -- checks/finite_scan inputs/finite_cases.txt "$(cat outputs/finite_modulus.txt)"
"$PYTHON" code/verify_final.py
"$PYTHON" code/maximal_squares.py
"$PYTHON" code/build_tables.py
printf '%s\n' 'COMPLETE: consolidated census, BBB/CCC extension proofs and corrected CCC equivalence.'
