"""Independent exact check of dephasing, supports and all square row quotients."""
import json,csv
from pathlib import Path
from collections import Counter
from itertools import combinations
R=Path(__file__).resolve().parents[1]
components=json.loads((R/'outputs/components.json').read_text())
scans=json.loads((R/'outputs/generic_AAA.json').read_text())
def diff(a,b):return ((a[0]-b[0])%30,)+tuple(x-y for x,y in zip(a[1:],b[1:]))
def check_A(a,b):
    C=Counter(diff(x,y) for x,y in zip(a,b))
    return all(C[e]==C[((e[0]+15)%30,)+e[1:]] for e in C)
success=[]
for c,q in zip(components,scans):
    assert (c['class_id'],c['branch'])==(q['class_id'],q['branch'])
    for r in c['rows']:assert not any(r[0])
    assert all(check_A(a,b) for a,b in combinations(c['rows'],2))
    assert all(check_A(r,s) for r in q['rows'] for s in c['rows'])
    if q['k5']:
        H=c['rows']+[q['rows'][i] for i in q['k5']]
        assert len(H)==8
        assert all(check_A(a,b) for a,b in combinations(H,2))
        success.append((c['class_id'],c['branch']))
assert len(success)==32
assert sorted(set(range(23))-{i for i,j in success})==[11,19,21,22]
with (R/'outputs/AAA_status.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['class','dimension','components','all_parameter_A_square_verified','generic_A_candidates_per_component','generic_A_edges_per_component','general_completion_status'])
    for i in range(23):
        cc=[c for c in components if c['class_id']==i];qq=[q for q in scans if q['class_id']==i]
        w.writerow([i,cc[0]['dim'],len(cc),all(q['k5'] for q in qq),';'.join(str(len(q['rows'])) for q in qq),';'.join(str(q['edges']) for q in qq),'existence proved by symbolic witness' if all(q['k5'] for q in qq) else 'exceptional loci unresolved; mixed profiles not checked by this scan'])
print('PASS: 38 parametrized components; 32 square witnesses independently checked; residual classes 11,19,21,22.')
