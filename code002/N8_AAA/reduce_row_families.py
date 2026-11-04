"""Remove affine row families contained in other row families, exactly."""
import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
data=json.loads((R/'outputs/generic_all.json').read_text())
def contained(a,b,dim):
    # Does a lie in b? Each free parameter of b is a free DSU component
    # rotation, and has an anchor entry with coefficient one.
    na,nb=a['free'],b['free']
    if na>nb:return None
    ra,rb=a['row'],b['row'];values=[]
    for k in range(nb):
        j=next(j for j,e in enumerate(rb) if e[dim+1+k])
        assert rb[j][dim+1+k]==1 and sum(rb[j][dim+1:])==1
        values.append([(ra[j][0]-rb[j][0])%30]+[ra[j][t]-rb[j][t] for t in range(1,dim+1)]+ra[j][dim+1:])
    for x,y in zip(ra,rb):
        z=y[:dim+1]+[0]*na
        for k,v in enumerate(values):
            for t,u in enumerate(v):z[t]+=y[dim+1+k]*u
        z[0]%=30
        if z!=x:return None
    return values
out=[]
for d in data:
    fs=d['families'];keep=[];maps=[]
    order=sorted(range(len(fs)),key=lambda i:(-fs[i]['free'],i))
    for i in order:
        for j in keep:
            v=contained(fs[i],fs[j],d['dim'])
            if v is not None:maps.append({'i':i,'j':j,'parameters':v});break
        else:keep.append(i)
    out.append({**d,'maximal_family_indices':sorted(keep),'containment_maps':maps})
    print(d['class_id'],d['branch'],'families',len(fs),'->',len(keep),'free-counts',{f:sum(fs[i]['free']==f for i in keep) for f in range(3)})
(R/'outputs/generic_reduced.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
with (R/'inputs/reduced_families.txt').open('w') as f:
    f.write(str(len(out))+'\n')
    for ci,d in enumerate(out):
        fs=[d['families'][i] for i in d['maximal_family_indices']]
        f.write(f"{ci} {d['dim']} {len(fs)}\n")
        for x in fs:f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n')
