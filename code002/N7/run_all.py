"""Reproduce every exact certificate in the consolidated first draft."""
from pathlib import Path
import argparse,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
STAGES=[('ccc','run_ccc.py'),('ab','run_ab.py'),('mixed','run_mixed.py'),
        ('square','verify_square_families.py'),('f4','run_all.py'),
        ('isolated45','run_all.py'),('affine4','run_all.py'),('completion','run_all.py')]
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',action='append',choices=[s for s,e in STAGES],help='Run only the selected stage; repeat to select several.')
    parser.add_argument('--list',action='store_true',help='List the stages without running them.')
    args=parser.parse_args()
    if args.list:
        for stage,entry in STAGES:print(stage,'->','code/'+stage+'/'+entry)
        return
    if not __debug__:raise SystemExit('Assertions must be enabled; do not use -O.')
    (ROOT/'logs').mkdir(exist_ok=True)
    for stage,entry in STAGES:
        if args.stage and stage not in args.stage:continue
        start=time.monotonic();print('Running',stage,flush=True)
        log=ROOT/'logs'/(stage+'.log')
        with log.open('w') as stream:
            result=subprocess.run([sys.executable,'-u',entry],cwd=ROOT/'code'/stage,stdout=stream,stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.read_text());raise SystemExit(result.returncode)
        print('Passed:',stage,'in',round(time.monotonic()-start,2),'seconds',flush=True)
    print('All selected exact-certificate stages passed.',flush=True)
if __name__=='__main__':main()
