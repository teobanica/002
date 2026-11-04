"""Build the complete curve arrangement outside certified completion planes."""
from torus3 import *
from itertools import combinations
R=Path(__file__).resolve().parents[1]
out=[];extra=[]
goods=json.load(open(R/'outputs/good_planes.json'))
planes=json.load(open(R/'outputs/planes.json'))
plane_pairs=json.load(open(R/'outputs/planes_pairs.json'))
for cid in (11,19):
    full=json.load(open(R/f'outputs/full{cid}.json'))
    good=[decode(x['torus']) for x in goods if x['class_id']==cid]
    H=next(x['rows'] for x in json.load(open(R/'inputs/components.json')) if x['class_id']==cid)
    ts={};ps={}
    def add(t):
        if any(subset(t,g) for g in good):return
        if t['dim']==1:ts[key(t)]=t
        elif t['dim']==0:ps[key(t)]=t
    for fam in full['row_families']:
        if len(fam['eq'])>=2:
            for t in split(fam['eq']):add(t)
    for ev in plane_pairs:
        if ev['class_id']==cid:
            for eq in ev['loci']:
                if len(eq)>=2:
                    for t in split(eq,int((R/'outputs/planes_pairs_modulus.txt').read_text())):add(t)
    pls=[decode(x['torus']) for x in planes if x['class_id']==cid and not any(subset(decode(x['torus']),g) for g in good)]
    for a,b in combinations(pls,2):
        for t in intersection(a,b):add(t)
    for k,t in sorted(ts.items()):out.append({'class_id':cid,'torus':t,'H':H})
    for t in ps.values():extra.append({'class_id':cid,'torus':t})
    print(cid,len(ts),'curves',len(ps),'pair-event points',flush=True)
mod=write_scan(R/'inputs/lines.txt',out)
(R/'outputs/lines.json').write_text(json.dumps(enc(out),indent=2))
(R/'outputs/lines_modulus.txt').write_text(str(mod))
(R/'outputs/plane_pair_points.json').write_text(json.dumps(enc(extra),indent=2))
# Collect every full regular row family containing each curve. This includes B.
scans=[]
for ix,it in enumerate(out):
    cid=it['class_id'];t=it['torus'];full=json.load(open(R/f'outputs/full{cid}.json'));rows=set()
    for fam in full['row_families']:
        ok=True
        for e in fam['eq']:
            if (F(e[3],30)+sum(F(e[i])*t['o'][i] for i in range(3)))%1:ok=False;break
            if any(sum(e[i]*t['V'][i][j] for i in range(3)) for j in range(t['dim'])):ok=False;break
        if ok:
            for r in specialize(fam['rows'],t,mod):rows.add(tuple(map(tuple,r)))
    scans.append({'class_id':cid,'id':ix,'dim':1,'rows':sorted(rows),'k5':[]})
(R/'outputs/lines_scan.json').write_text(json.dumps(scans,indent=2))
print('Complete curve candidates saved, modulus',mod)
