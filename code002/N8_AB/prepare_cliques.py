from pathlib import Path
import sys,json
HOME_ROOT=Path(__file__).resolve().parents[1]
R=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HOME_ROOT
strata=json.loads((R/'outputs/strata.json').read_text())['strata']
pairs=json.loads((R/'outputs/pairs_strata_maximal.json').read_text())
with (R/'inputs/cliques_strata.txt').open('w') as f:
 f.write(str(len(strata))+'\n')
 for s,p in zip(strata,pairs):
  assert s['id']==p['class_id']
  fs=s['families'];f.write(f"{s['id']} {s['base']} {len(fs)}\n"+' '.join(str(x['free']) for x in fs)+'\n'+str(len(p['pairs']))+'\n')
  for q in p['pairs']:
   f.write(f"{q['i']} {q['j']} {len(q['relations'])}\n")
   for r in q['relations']:
    f.write(str(len(r['eq']))+'\n')
    for e in r['eq']:f.write(' '.join(map(str,e))+'\n')
