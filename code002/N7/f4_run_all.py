from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name in ['pilot.py','independent.py','verify_candidates.py','compatibility.py']:
    print('Running',name,flush=True)
    subprocess.run([sys.executable,'-u',str(root/name)],cwd=root,check=True)
print('All exact pilot checks passed.')
