import json
from pathlib import Path
from fractions import Fraction as F
from math import lcm
from itertools import combinations
from torus import *
R=Path(__file__).resolve().parents[1]
C=[c for c in json.loads((R/'outputs/components.json').read_text()) if c['class_id']==22]
G=[c for c in json.loads((R/'outputs/generic_AAA.json').read_text()) if c['class_id']==22]
EV=json.loads((R/'outputs/class22_A_events.json').read_text())
out=[];mod=30
for c,g,ev in zip(C,G,EV):
    lines=set(decode_line(l) for l in ev['lines']);points=set(decode_point(p) for p in ev['points']);pairlines=set();pairpoints=set()
    for r,s in combinations(g['rows'],2):
        diff=[[F(x[0]-y[0],30)]+[x[k]-y[k] for k in (1,2)] for x,y in zip(r,s)]
        for m in MATCH:
            eq=[(diff[a][1]-diff[b][1],diff[a][2]-diff[b][2],F(1,2)-diff[a][0]+diff[b][0]) for a,b in m]
            full,ll,pp=solve(eq);pairlines.update(ll);pairpoints.update(pp)
    lines.update(pairlines);points.update(pairpoints)
    for i,line in enumerate(sorted(lines)):
        p,v=line_param(line);rows=[]
        for r in c['rows']:
            rr=[]
            for z,a,b in r:
                q=(F(z,30)+a*p[0]+b*p[1])%1;d=a*v[0]+b*v[1];rr.append([q,d]);mod=lcm(mod,q.denominator)
            rows.append(rr)
        out.append(dict(class_id=i,branch=c['branch'],dim=1,line=line,rows=rows))
    print('branch',c['branch'],'row-event lines',len(ev['lines']),'pair-event lines',len(pairlines),'combined',len(lines),'isolated events',len(points))
(R/'outputs/class22_line_parameters.json').write_text(json.dumps(encode(out),indent=2))
with (R/'inputs/class22_lines.txt').open('w') as f:
    f.write(str(len(out))+'\n')
    for c in out:
        f.write(f"{c['class_id']} {c['branch']} 1\n")
        for r in c['rows']:f.write(' '.join(f'{int(q*mod)} {a}' for q,a in r)+'\n')
(R/'outputs/line_modulus.txt').write_text(str(mod))
print('modulus',mod)
