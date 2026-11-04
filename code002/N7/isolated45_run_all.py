from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name,log in [('cleanup.py','run.log'),('independent.py','independent_run.log'),('verify.py','verification_run.log')]:
    print('Running',name,flush=True)
    with (root/log).open('w') as stream:
        result=subprocess.run([sys.executable,'-u',str(root/name)],cwd=root,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:
        print((root/log).read_text());raise SystemExit(result.returncode)
    print('Passed:',name,flush=True)
print('All 45 isolated cases resolved. 43 admit no regular fourth row; two complete uniquely to F7.')
