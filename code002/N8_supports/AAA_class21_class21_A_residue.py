"""Reduce the A-only converse to a finite set, with exact event equations."""
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from torus import *
R=Path(__file__).resolve().parents[1]
L=json.loads((R/'outputs/class21_line_parameters.json').read_text())
S=json.loads((R/'outputs/class21_line_scan.json').read_text())
C=[c for c in json.loads((R/'outputs/components.json').read_text()) if c['class_id']==21]
EV=json.loads((R/'outputs/class21_A_events.json').read_text())
mod=int((R/'outputs/line_modulus.txt').read_text());result=[]
for br in range(2):
    pairs=[(l,s) for l,s in zip(L,S) if l['branch']==br]
    good=[decode_line(l['line']) for l,s in pairs if s['k5']]
    bad=[(l,s) for l,s in pairs if not s['k5']]
    points=set(decode_point(p) for p in EV[br]['points'])
    alllines=[decode_line(l['line']) for l,s in pairs]
    for a,b in combinations(alllines,2):points.update(solve([a,b])[2])
    jump=[]
    for l,s in bad:
        assert s['disconnected_supports']==0
        line=decode_line(l['line']);p,v=line_param(line);ts=set()
        for rr,ss in combinations(s['rows'],2):
            diff=[(F(a[0]-b[0],mod),a[1]-b[1]) for a,b in zip(rr,ss)]
            for m in MATCH:
                eq=[(diff[a][1]-diff[b][1],(F(1,2)-diff[a][0]+diff[b][0])%1) for a,b in m]
                if any(k==0 and c for k,c in eq):continue
                nonzero=next(((k,c) for k,c in eq if k),None)
                if nonzero is None:continue
                k,c=nonzero
                for j in range(abs(k)):
                    t=((c+j)/k)%1
                    if all((kk*t-cc)%1==0 for kk,cc in eq):ts.add(t)
        pp={((p[0]+v[0]*t)%1,(p[1]+v[1]*t)%1) for t in ts}
        points.update(pp);jump.append(dict(line=line,parameters=sorted(ts),points=sorted(pp)))
    off=sorted(p for p in points if not any(on_line(p,l) for l in good))
    result.append(dict(branch=br,good_lines=good,bad_lines=[decode_line(l['line']) for l,s in bad],pair_jumps=jump,all_event_points=sorted(points),off_good_points=off))
    print('branch',br,'good circles',len(good),'bad circles',len(bad),'event points',len(points),'off good',len(off))
(R/'outputs/class21_A_residue.json').write_text(json.dumps(encode(result),indent=2))
out=[];mm=30
for res,c in zip(result,C):
    for i,p in enumerate(res['off_good_points']):
        rows=[[(F(z,30)+a*p[0]+b*p[1])%1 for z,a,b in r] for r in c['rows']]
        for r in rows:
            for q in r:mm=lcm(mm,q.denominator)
        out.append((i,c['branch'],rows))
with (R/'inputs/class21_A_points.txt').open('w') as f:
    f.write(str(len(out))+'\n')
    for i,br,rows in out:
        f.write(f'{i} {br} 0\n')
        for r in rows:f.write(' '.join(str(int(q*mm)) for q in r)+'\n')
(R/'outputs/A_point_modulus.txt').write_text(str(mm))
print('point modulus',mm)
