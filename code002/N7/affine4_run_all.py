from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
jobs=[('census.py','run.log'),('independent.py','independent_run.log'),('analyze.py','analysis_run.log'),
      ('f4_special.py','f4_special_run.log'),('verify_export.py','verification_run.log')]
for script,log in jobs:
    print('Running',script,flush=True)
    with (root/log).open('w') as f:
        p=subprocess.run([sys.executable,'-u',str(root/script)],cwd=root,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:print((root/log).read_text());raise SystemExit(p.returncode)
    print('Passed:',script,flush=True)
print('Complete six-family fourth-row census and all reduction checks passed.')
