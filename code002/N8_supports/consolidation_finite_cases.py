"""Group exact point problems under explicit row/column Hadamard equivalences."""
from torus3 import *
from itertools import permutations
R=Path(__file__).resolve().parents[1]

cases={};maps=[]
for idx,x in enumerate(json.load(open(R/'outputs/pure_points.json'))):
    cid=x['class_id'];p=[F(*q) for q in x['point']]
    H=[tuple(F(*e) for e in row) for row in x['H']]
    best=None;choice=None
    for a,b,c in permutations(range(3)):
        u=[(H[b][j]-H[a][j])%1 for j in range(8)];v=[(H[c][j]-H[a][j])%1 for j in range(8)]
        for col in range(8):
            order=sorted(range(8),key=lambda j:((u[j]-u[col])%1,(v[j]-v[col])%1))
            key0=tuple(((u[j]-u[col])%1,(v[j]-v[col])%1) for j in order)
            if best is None or key0<best:best=key0;choice=a,b,c,col,order
    a,b,c,col,order=choice
    rows=set()
    for row in x['rows']:
        w=[F(*e) for e in row];rows.add(tuple((w[j]-H[a][j]-w[col]+H[a][col])%1 for j in order))
    # Use type-independent canonical triples; C is included when needed by modulus.
    if best not in cases:cases[best]={'id':len(cases),'H':[[F(0)]*8,[q[0] for q in best],[q[1] for q in best]],'rows':sorted(rows),'members':[]}
    case=cases[best];assert set(case['rows'])==rows
    case['members'].append(idx);maps.append({'point_index':idx,'case_id':case['id'],'row_order':[a,b,c],'pivot':col,'column_order':order})
out=sorted(cases.values(),key=lambda x:x['id'])
M=lcm(30,*(e.denominator for c in out for row in c['H'] for e in row))
with (R/'inputs/finite_cases.txt').open('w') as f:
    f.write(str(len(out))+'\n')
    for c in out:
        f.write(f"{c['id']} 0 0\n")
        for row in c['H']:f.write(' '.join(str(int(e*M)) for e in row)+'\n')
(R/'outputs/finite_cases.json').write_text(json.dumps(enc(out),indent=2))
(R/'outputs/finite_case_maps.json').write_text(json.dumps(maps,indent=2))
(R/'outputs/finite_modulus.txt').write_text(str(M))
print(len(maps),'finite triples reduce to',len(out),'explicit canonical triples; modulus',M)
