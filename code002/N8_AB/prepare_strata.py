"""Split the complete fourth-row lists into generic circles and exceptional points.

Every point with a non-generic fourth row is listed. Generic row lists remain
valid on the entire circle, so their later square search also retains any
additional special points imposed by mutual compatibility. Equivalent point
seeds are identified by an explicit ordinary monomial transformation.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
from math import lcm
import json
R=Path(__file__).resolve().parents[1]
data=json.loads((R/'outputs/reduced_rows_all.json').read_text())
seeds=json.loads((R/'inputs/circle_components.json').read_text())
MOD=60

def contains(a,b,base):
    if a['free']>b['free']:return None
    w=1+base+a['free'];values=[]
    for k in range(b['free']):
        j=next(j for j,e in enumerate(b['row']) if e[1+base+k])
        assert b['row'][j][1+base+k]==1 and sum(b['row'][j][1+base:])==1
        values.append([a['row'][j][s]-b['row'][j][s] for s in range(1+base)]+a['row'][j][1+base:])
    for x,y in zip(a['row'],b['row']):
        z=y[:1+base]+[0]*a['free']
        for k,v in enumerate(values):
            for s in range(w):z[s]+=y[1+base+k]*v[s]
        if (z[0]-x[0])%MOD or z[1:]!=x[1:]:return None
    return values

def reduce(fs,base):
    keep=[];maps=[]
    for i in sorted(range(len(fs)),key=lambda i:(-fs[i]['free'],i)):
        for j in keep:
            v=contains(fs[i],fs[j],base)
            if v is not None:maps.append(dict(source=i,target=j,substitution=v));break
        else:keep.append(i)
    return [fs[i] for i in keep],dict(keep=keep,maps=maps)

def canonical_seed(H):
    best=None;cert=None
    for rp in permutations(range(3)):
        for p in range(8):
            cols=[tuple((H[r][j]-H[r][p]-H[rp[0]][j]+H[rp[0]][p])%MOD for r in rp) for j in range(8)]
            cp=sorted(range(8),key=lambda j:cols[j]);key=tuple(cols[j] for j in cp)
            if best is None or key<best:best=key;cert=dict(rows=rp,pivot_column=p,columns=cp)
    return best,cert

strata=[];point_index={};point_sources=[]
for d,s in zip(data,seeds):
    assert d['id']==s['id'];cid=d['id'];fs=d['families']
    generic=[dict(f,row=[[e[0]*2]+e[1:] for e in f['row']],raw_index=i) for i,f in enumerate(fs) if not f['eq']]
    generic,reduction=reduce(generic,1)
    strata.append(dict(id=len(strata),circle=cid,kind='circle',base=1,modulus=MOD,rows=[[[e[0]*2]+e[1:] for e in row] for row in s['rows']],families=generic,reduction=reduction))
    points=sorted({(F(j)-F(c,30))/k%1 for f in fs for k,c in f['eq'] for j in range(abs(k))})
    for t in points:
        assert (t*MOD).denominator==1
        H=[[(e[0]*2+e[1]*int(t*MOD))%MOD for e in row] for row in s['rows']]
        key,cert=canonical_seed(H)
        origin=dict(circle=cid,point=[t.numerator,t.denominator],canonical_transform=cert)
        if key in point_index:
            origin['representative']=point_index[key];point_sources.append(origin);continue
        good=[]
        for i,f in enumerate(fs):
            if any((k*t+F(c,30))%1 for k,c in f['eq']):continue
            good.append(dict(f,eq=[],row=[[(e[0]*2+e[1]*int(t*MOD))%MOD]+e[2:] for e in f['row']],raw_index=i))
        reduced,reduction=reduce(good,0)
        sid=len(strata);point_index[key]=sid;origin['representative']=sid;point_sources.append(origin)
        strata.append(dict(id=sid,circle=cid,kind='point',point=origin['point'],base=0,modulus=MOD,rows=[[[e] for e in row] for row in H],families=reduced,reduction=reduction,canonical_transform=cert))

(R/'outputs/strata.json').write_text(json.dumps(dict(modulus=MOD,strata=strata,point_sources=point_sources),separators=(',',':'))+'\n')
with (R/'inputs/pairs_strata.txt').open('w') as f:
    f.write(str(len(strata))+'\n')
    for s in strata:
        f.write(f"{s['id']} {s['base']} {len(s['families'])}\n")
        for x in s['families']:
            f.write(str(x['free'])+'\n'+' '.join(str(v) for e in x['row'] for v in e)+'\n0\n')
print('Strata',len(strata),'circles',sum(s['kind']=='circle' for s in strata),'distinct point seeds',len(point_index),'point presentations',len(point_sources))
for s in strata:print(s['id'],s['circle'],s['kind'],s.get('point',''),len(s['families']),[f['free'] for f in s['families']])
