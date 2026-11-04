"""Independent entrywise checks of finite reductions and normal-form certificates."""
import json
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from exact_rows import profiles, quotient

R=Path(__file__).resolve().parents[1]
def read(name): return json.loads((R/'outputs'/name).read_text())
def frac(rows): return [tuple(F(*x) for x in row) for row in rows]
for old,new in [('n8_pilot_results.json','pilot.json'),('n8_canonical_results_verified.json','canonical.json'),('n8_affine_results.json','affine.json')]:
    assert json.loads((R/'inputs'/old).read_text())==read(new),new

points=read('pure_points.json');cases=read('finite_cases.json');maps=read('finite_case_maps.json');scans=read('finite_scan.json')
assert len(points)==len(maps)==165 and len(cases)==len(scans)==14
assert [c['id'] for c in cases]==[s['id'] for s in scans]==list(range(14))
for c,s in zip(cases,scans):
    assert set(frac(c['rows']))=={tuple(F(x,s['modulus']) for x in row) for row in s['rows']}
assert sorted(m['point_index'] for m in maps)==list(range(165))
for m in maps:
    p=points[m['point_index']];c=cases[m['case_id']]
    H=frac(p['H']);rr=m['row_order'];cc=m['column_order'];pivot=m['pivot'];a=rr[0]
    assert sorted(rr)==list(range(3)) and sorted(cc)==list(range(8))
    def transform(row):return tuple((row[j]-H[a][j]-row[pivot]+H[a][pivot])%1 for j in cc)
    assert [transform(H[i]) for i in rr]==frac(c['H'])
    assert {transform(row) for row in frac(p['rows'])}==set(frac(c['rows']))
    assert m['point_index'] in c['members']

# Check the graph itself, its two components, and each component's profiles.
square_inventory=set()
for p in points:
    if p['class_id']!=101:continue
    rows=[[(e,) for e in row] for row in frac(p['rows'])]
    H=[[(e,) for e in row] for row in frac(p['H'])]
    adj=[set() for _ in rows]
    for i,j in combinations(range(len(rows)),2):
        if profiles(quotient(rows[i],rows[j])):adj[i].add(j);adj[j].add(i)
    todo=set(range(10));components=[]
    while todo:
        component={min(todo)};front=list(component)
        while front:
            for j in adj[front.pop()]-component:component.add(j);front.append(j)
        todo-=component;components.append(component)
    assert sorted(map(len,components))==[5,5]
    for component in components:
        assert all(adj[i]==component-{i} for i in component)
        types=Counter(''.join(profiles(quotient(rows[i],h)) for h in H) for i in component)
        assert types=={'AAA':4,'CCC':1}
        square_inventory.add((tuple(frac(p['H'])),tuple(sorted(frac(p['rows'])[i] for i in component))))
old=json.loads((R/'inputs/n8_ccc_square_results.json').read_text())
old_inventory=set()
for square in old['square_completions']:
    H=[tuple(F(x,30) for x in row) for row in square['rows_exp30']]
    old_inventory.add((tuple(H[:3]),tuple(sorted(H[3:]))))
assert len(square_inventory)==len(old_inventory)==len(old['square_completions'])==60
assert square_inventory==old_inventory

nf=read('pure_normal_forms.json')
target=next(c['rows'] for c in json.loads((R/'inputs/pure_components.json').read_text()) if c['name']=='BBB')
signs=Counter()
for c in nf['BBB_equivalence_certificates']:
    H=c['component'];rr=c['row_order'];cc=c['columns'];a,b=c['a_phase'],c['b_phase'];sg=c['target_root_sign'];signs[sg]+=1
    assert sorted(rr)==list(range(3)) and sorted(cc)==list(range(8))
    assert abs(a[1]*b[2]-a[2]*b[1])==1
    for i in range(3):
        for j in range(8):
            got=[H[rr[i]][cc[j]][k]-H[rr[0]][cc[j]][k]-H[rr[i]][cc[0]][k]+H[rr[0]][cc[0]][k] for k in range(3)]
            e=target[i][j]
            want=[sg*e[0]+e[1]*a[0]+e[2]*b[0],e[1]*a[1]+e[2]*b[1],e[1]*a[2]+e[2]*b[2]]
            got[0]%=60;want[0]%=60
            assert got==want
assert signs=={1:4,-1:4}
summary={'census_outputs_matching_archive':3,'finite_enumerations_matching':14,'finite_triple_maps_verified':165,'CCC_graphs_explicitly_two_K5':30,'CCC_square_inventory_exactly_reproduced':60,'BBB_symbolic_normal_form_certificates_verified':8,'BBB_sign_counts':dict(signs)}
(R/'outputs/final_verification.json').write_text(json.dumps(summary,indent=2))
print('PASS:',json.dumps(summary))
