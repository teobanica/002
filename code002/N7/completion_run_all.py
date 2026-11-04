"""Reproduce the final completion step and its independent verification."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for script,log in [('f2.py','f2_run.log'),('f3.py','f3_run.log'),('check.py','check_run.log')]:
    print('Running',script,flush=True)
    with (root/log).open('w') as f:
        result=subprocess.run([sys.executable,'-u',str(root/script)],cwd=root,stdout=f,stderr=subprocess.STDOUT)
    if result.returncode:
        print((root/log).read_text());raise SystemExit(result.returncode)
    print('Passed:',script,flush=True)
print('Both remaining sectors are completely covered by the F4 completion theorem.')
