"""Include regular (A/B/C) pair events, not only antipodal pair events."""
import json
from pathlib import Path
from torus import *
from fractions import Fraction as F
from math import lcm
R=Path(__file__).resolve().parents[1]
GP=json.loads((R/'outputs/global_regular_pair_events.json').read_text())
LP=json.loads((R/'outputs/line_regular_pair_events.json').read_text())
L=json.loads((R/'outputs/class21_line_parameters.json').read_text())
C=[c for c in json.loads((R/'outputs/components.json').read_text()) if c['class_id']==21]
out=[];inputs=[];mod=30
for br in range(2):
    D=json.loads((R/f'outputs/class21_full{br}.json').read_text())
    AR=json.loads((R/'outputs/class21_A_residue.json').read_text())[br]
    good=[decode_line(l) for l in AR['good_lines']]
    points={decode_point(f['point']) for f in D['fibers']}
    for rr in GP[br]['profiles'].values():
        points.update(decode_point(p) for p in rr['points'])
        assert all(decode_line(l) in good for l in rr['lines'])
    for lp in LP:
        if lp['branch']!=br:continue
        l=next(l for l in L if l['class_id']==lp['id'] and l['branch']==br)
        p,v=line_param(decode_line(l['line']))
        for rr in lp['profiles'].values():
            assert not rr['points']
            for eq in rr['lines']:
                a,b,c=decode_line(eq);assert (a,b)==(1,0)
                points.add(((p[0]+v[0]*c)%1,(p[1]+v[1]*c)%1))
    off=sorted(p for p in points if not any(on_line(p,l) for l in good))
    old={decode_point(f['point']) for f in D['fibers']}
    print('branch',br,'final off-good points',len(off),'additional',len(set(off)-old))
    out.append(dict(branch=br,points=off,additional=sorted(set(off)-old)))
    for i,p in enumerate(off):
        rows=[[(F(z,30)+a*p[0]+b*p[1])%1 for z,a,b in r] for r in C[br]['rows']]
        for r in rows:
            for q in r:mod=lcm(mod,q.denominator)
        inputs.append((i,br,rows))
(R/'outputs/final_residue.json').write_text(json.dumps(encode(out),indent=2))
with (R/'inputs/final_points.txt').open('w') as f:
    f.write(str(len(inputs))+'\n')
    for i,br,rows in inputs:
        f.write(f'{i} {br} 0\n')
        for r in rows:f.write(' '.join(str(int(q*mod)) for q in r)+'\n')
(R/'outputs/final_point_modulus.txt').write_text(str(mod))
print('modulus',mod)
