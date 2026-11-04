"""Deterministic text adapters between the independent exact engines."""
from phase_algebra import *
import argparse
p=argparse.ArgumentParser();p.add_argument('stage',choices=['maximal','cliques']);a=p.parse_args()
if a.stage=='maximal':
    ds=read('generic_all.json')[:2]
    with (R/'inputs/maximal_families.txt').open('w') as f:
        f.write('2\n')
        for d in ds:
            fs=d['families'];f.write(f"{d['class_id']} {d['dim']} {len(fs)}\n")
            for x in fs:f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n')
            (R/f"outputs/generic{d['class_id']}.json").write_text(json.dumps([d],separators=(',',':'))+'\n')
else:
    ds=read('generic_reduced.json');ps=read('reduced_pairs.json')
    with (R/'inputs/reduced_cliques.txt').open('w') as f:
        f.write(str(len(ds))+'\n')
        for ci,(d,pair) in enumerate(zip(ds,ps)):
            fs=[d['families'][i] for i in d['maximal_family_indices']]
            f.write(f"{ci} {d['dim']} {len(fs)}\n"+' '.join(str(x['free']) for x in fs)+'\n'+str(len(pair['pairs']))+'\n')
            for x in pair['pairs']:
                f.write(f"{x['i']} {x['j']} {len(x['relations'])}\n")
                for rel in x['relations']:
                    f.write(str(len(rel['eq']))+'\n')
                    for row in rel['eq']:f.write(' '.join(map(str,row))+'\n')
print('Prepared',a.stage,'inputs')
