"""Run both standard-library exact checks. No external dependencies."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
output = []
for name in ('classify_ab.py', 'independent_check.py'):
    result = subprocess.run([sys.executable, str(root/name)], cwd=root,
                            text=True, capture_output=True)
    print(result.stdout, end='')
    output.append(result.stdout)
    if result.returncode:
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)
(root/'run.log').write_text(''.join(output))
print('All exact checks passed.')
