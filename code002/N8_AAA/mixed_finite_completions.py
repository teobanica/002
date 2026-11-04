"""Enumerate every square over the 32 exceptional class-11 triples."""
from phase_algebra import *
P=R
pts=json.loads((P/'inputs/prior_final_points.json').read_text())
cs={c['class_id']:c['rows'] for c in json.loads((P/'inputs/prior_class11_components.json').read_text())}
out=[];summary=[]
for x in pts:
    if not x['k5']:continue
    assert x['class_id']==11
    p=[F(*v) for v in x['point']]
    H=[tuple(((F(e[0],30)+sum(F(e[i+1])*p[i] for i in range(3)))%1,) for e in row) for row in cs[11]]
    rr=[tuple((F(*e),) for e in row) for row in x['rows']]
    adj=[[i!=j and bool(profiles(quotient(a,b))) for j,b in enumerate(rr)] for i,a in enumerate(rr)]
    ks=cliques(adj,5);summary.append({'point':x['point'],'rows':len(rr),'squares':len(ks)})
    for k in ks:
        sq=H+[rr[i] for i in k]
        typ=Counter(profiles(quotient(sq[i],sq[j])) for i,j in combinations(range(8),2))
        assert '' not in typ
        assert all((e[0]*12).denominator==1 for row in sq for e in row)
        out.append({'point':x['point'],'clique':k,'profiles':dict(typ),'matrix_mod12':[[int(e[0]*12) for e in row] for row in sq]})
assert len(summary)==32
(R/'outputs/mixed_finite_squares.json').write_text(json.dumps(out,indent=2))
(R/'outputs/mixed_finite_summary.json').write_text(json.dumps(summary,indent=2))
(R/'inputs/mixed_squares.txt').write_text(f'{len(out)} 12\n'+'\n'.join(' '.join(str(e) for row in x['matrix_mod12'] for e in row) for x in out)+'\n')
print('All completions:',len(out),'cliques per point:',Counter(x['squares'] for x in summary),'profiles:',Counter(str(x['profiles']) for x in out))
