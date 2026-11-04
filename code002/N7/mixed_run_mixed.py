"""Run the complete exact classification and both verification passes."""
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
lines=[]
for name in ('enumerate_mixed.py','verify_mixed.py','full_support_check.py'):
    result=subprocess.run([sys.executable,'-u',str(root/name)],cwd=root,text=True,capture_output=True)
    print(result.stdout,end='',flush=True)
    lines.append(result.stdout)
    if result.returncode:
        print(result.stderr,file=sys.stderr)
        raise SystemExit(result.returncode)
(root/'run.log').write_text(''.join(lines))
print('All mixed-type checks passed.',flush=True)
