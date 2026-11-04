from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
ds=json.loads((R/'outputs/fixed_pairs.json').read_text())
ps=json.loads((R/'outputs/fixed_pair_compatibility.json').read_text())
with (R/'inputs/fixed_pair_cliques.txt').open('w') as f:
 f.write('2\n')
 for ci,(d,p) in enumerate(zip(ds,ps)):
  fs=[d['families'][i] for i in d['keep']]
  f.write(f'{ci} 0 {len(fs)}\n'+' '.join(str(x['free']) for x in fs)+'\n'+str(len(p['pairs']))+'\n')
  for x in p['pairs']:
   f.write(f"{x['i']} {x['j']} {len(x['relations'])}\n")
   for rel in x['relations']:
    f.write(str(len(rel['eq']))+'\n')
    for row in rel['eq']:f.write(' '.join(map(str,row))+'\n')
print('Prepared both fixed C-pair clique inputs.')
