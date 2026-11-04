"""Rerun the original exact CCC seed search and check its nonextension result."""
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
result = subprocess.run([str(root / 'checks/seed_check')], capture_output=True, text=True)
assert result.returncode == 1, (result.returncode, result.stdout, result.stderr)
assert not result.stdout.strip(), result.stdout
assert result.stderr.splitlines() == ['signatures=3870', 'tested=3923640 found=0'], result.stderr
message = 'PASS: expected nonextension exit status 1.\n' + result.stderr
(root / 'checks/ccc_seed_check.txt').write_text(message)
print(message, end='')
