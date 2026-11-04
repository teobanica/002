from pathlib import Path
import json,sys
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from exact_rows import profiles,quotient
D=json.loads((R/'outputs/c_pair_rows.json').read_text())
H1=[[0],[0],[0],[0],[0],[0],[0],[0]]
base=[0,10,20,0,6,12,18,24]
def prof(row):return profiles([(F(e[0],30),)+tuple(e[1:]) for e in row])
def contained(a,b):
    na,nb=a['free'],b['free'];ra,rb=a['row'],b['row']
    if na>nb:return None
    vals=[]
    for k in range(nb):
        j=next(j for j,e in enumerate(rb) if e[k+1])
        assert rb[j][k+1]==1 and sum(rb[j][1:])==1
        vals.append([(ra[j][0]-rb[j][0])%30]+ra[j][1:])
    for x,y in zip(ra,rb):
        z=[y[0]]+[0]*na
        for k,v in enumerate(vals):
            for t,u in enumerate(v):z[t]+=y[k+1]*u
        z[0]%=30
        if z!=x:return None
    return vals
out=[]
for phase in [0,15]:
    allrows={}
    for ix,f in enumerate(D['families']):
        if any((eq[0]*phase+eq[1])%30 for eq in f['base_equations']):continue
        row=[[(e[0]+e[1]*phase)%30]+e[2:] for e in f['row']]
        key=tuple(tuple(e) for e in row)
        if key not in allrows:allrows[key]={'row':row,'free':f['free'],'sources':[]}
        allrows[key]['sources'].append(ix)
    fs=list(allrows.values());keep=[];maps=[]
    for i in sorted(range(len(fs)),key=lambda i:(-fs[i]['free'],i)):
        for j in keep:
            vals=contained(fs[i],fs[j])
            if vals is not None:maps.append({'i':i,'j':j,'parameters':vals});break
        else:keep.append(i)
    h=[(base[j]+(phase if j>=3 else 0))%30 for j in range(8)]
    for f in fs:
        f['profiles']=[prof(f['row']),prof([[(e[0]-h[j])%30]+e[1:] for j,e in enumerate(f['row'])])]
        assert all(f['profiles'])
    out.append({'base_phase_mod30':phase,'base_row':h,'families':fs,'keep':keep,'containment_maps':maps})
    print('z =',1 if phase==0 else -1,'all rows',len(fs),'free counts',dict(Counter(f['free'] for f in fs)))
    print('retained',len(keep),'profiles',dict(Counter((fs[i]['free'],''.join(fs[i]['profiles'])) for i in keep)))
(R/'outputs/fixed_pairs.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
with (R/'inputs/fixed_pair_families.txt').open('w') as f:
    f.write('2\n')
    for ci,d in enumerate(out):
        fs=[d['families'][i] for i in d['keep']]
        f.write(f'{ci} 0 {len(fs)}\n')
        for x in fs:f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n')
