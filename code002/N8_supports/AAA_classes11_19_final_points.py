from torus3 import *
from exact_rows import *
from itertools import combinations
R=Path(__file__).resolve().parents[1]
planes=json.load(open(R/'outputs/planes.json'));goods=json.load(open(R/'outputs/good_planes.json'))
lines=json.load(open(R/'outputs/lines.json'));scans=json.load(open(R/'outputs/lines_scan.json'));pairs=json.load(open(R/'outputs/lines_pairs.json'))
line_mod=int((R/'outputs/lines_modulus.txt').read_text());pair_mod=int((R/'outputs/lines_pairs_modulus.txt').read_text())
good_lines=[]
for it,sc,pa in zip(lines,scans,pairs):
    t=decode(it['torus']);rr=[[(F(e[0],line_mod),)+tuple(e[1:]) for e in row] for row in sc['rows']]
    H=[[(F(e[0],line_mod),)+tuple(e[1:]) for e in row] for row in specialize(it['H'],t,line_mod)]
    for row in rr:
        for h in H:assert profiles(quotient(row,h))
    k,cl,ed=graph(rr);assert bool(k)==bool(pa['k5'])
    if k:
        square=H+[rr[i] for i in k]
        for i,j in combinations(range(8),2):assert profile(quotient(square[i],square[j]),[2]*4)
        good_lines.append({**it,'square':square,'clique':k})
(R/'outputs/good_lines.json').write_text(json.dumps(enc(good_lines),indent=2))
print('Verified',len(good_lines),'additional A-only completion circles.',flush=True)
def on(p,t):return all((sum(x*y for x,y in zip(a,p))-c)%1==0 for a,c in zip(t['A'],t['c']))
def cross(line,t):
    o=line['o'];v=[x[0] for x in line['V']];eq=[]
    for a,c in zip(t['A'],t['c']):
        h=sum(x*y for x,y in zip(a,v));b=(c-sum(x*y for x,y in zip(a,o)))%1
        if not h:
            if b:return []
        else:eq.append((h,b))
    if not eq:return []
    h,b=min(eq,key=lambda x:abs(x[0]));out=[]
    for k in range(abs(h)):
        q=(b+k)/h
        if all((aa*q-bb)%1==0 for aa,bb in eq):out.append(tuple((x+q*y)%1 for x,y in zip(o,v)))
    return out
out=[];summary=[]
for cid in (11,19):
    full=json.load(open(R/f'outputs/full{cid}.json'))
    good=[decode(x['torus']) for x in goods+good_lines if x['class_id']==cid]
    pls=[decode(x['torus']) for x in planes if x['class_id']==cid and not any(subset(decode(x['torus']),g) for g in good)]
    ls=[decode(x['torus']) for x in lines if x['class_id']==cid and not any(subset(decode(x['torus']),g) for g in good)]
    pts=set()
    def add(p):
        p=tuple(p)
        if not any(on(p,g) for g in good):pts.add(p)
    for f in full['fibers']:add([F(*q) for q in f['point']])
    for pa in pairs:
        if pa['class_id']==cid:
            for eq in pa['loci']:
                if len(eq)==3:
                    for t in split(eq,pair_mod):add(t['o'])
    for pp in json.load(open(R/'outputs/plane_pair_points.json')):
        if pp['class_id']==cid:add(decode(pp['torus'])['o'])
    for ln in ls:
        for pl in pls:
            for p in cross(ln,pl):add(p)
    for a,b in combinations(ls,2):
        for p in cross(a,b):add(p)
    print(cid,'critical isolated points:',len(pts),flush=True)
    pointrows={tuple(F(*q) for q in f['point']):{tuple(F(*e) for e in row) for row in f['rows']} for f in full['fibers']}
    cs=next(x['rows'] for x in json.load(open(R/'inputs/components.json')) if x['class_id']==cid)
    def value(row,p):return tuple((F(e[0],30)+sum(F(e[i+1])*p[i] for i in range(3)))%1 for e in row)
    stats=Counter()
    for ix,p in enumerate(sorted(pts)):
        rr=set(pointrows.get(p,set()))
        for fam in full['row_families']:
            if all((F(e[3],30)+sum(F(e[i])*p[i] for i in range(3)))%1==0 for e in fam['eq']):
                for row in fam['rows']:rr.add(value(row,p))
        if p in pointrows:assert rr==pointrows[p]
        rows=sorted(rr);H=[value(row,p) for row in cs]
        symbolic=[[(e,) for e in row] for row in rows]
        for row in symbolic:
            for h in H:assert profiles(quotient(row,[(e,) for e in h]))
        k,cl,edges=graph(symbolic);stats[(len(rows),cl)]+=1
        out.append({'class_id':cid,'point':p,'rows':rows,'clique_number_up_to_5':cl,'edges':edges,'k5':k})
        if (ix+1)%200==0:print(cid,'checked',ix+1,'/',len(pts),flush=True)
    summary.append({'class_id':cid,'points':len(pts),'histogram':[[a,b,n] for (a,b),n in sorted(stats.items())]})
    print(cid,dict(stats),flush=True)
(R/'outputs/final_points.json').write_text(json.dumps(enc(out),indent=2))
(R/'outputs/final_points_summary.json').write_text(json.dumps(summary,indent=2))
