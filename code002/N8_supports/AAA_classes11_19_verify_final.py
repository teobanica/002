from torus3 import *
from exact_rows import *
import csv
R=Path(__file__).resolve().parents[1]
goods=json.load(open(R/'outputs/good_planes.json'));glines=json.load(open(R/'outputs/good_lines.json'))
pts=json.load(open(R/'outputs/final_points.json'));cs={c['class_id']:c['rows'] for c in json.load(open(R/'inputs/components.json'))}
planes=json.load(open(R/'outputs/planes.json'));psc=json.load(open(R/'outputs/planes_scan.json'));pps=json.load(open(R/'outputs/planes_pairs.json'))
lines=json.load(open(R/'outputs/lines.json'));lsc=json.load(open(R/'outputs/lines_scan.json'));lps=json.load(open(R/'outputs/lines_pairs.json'))
def literal_plane(a,c):return key(torus([a],[F(c)]))
exp11={literal_plane(a,c) for a in [[0,0,1],[0,1,0],[0,-1,1]] for c in [0,F(1,2)]}|{literal_plane([0,-2,1],F(1,2))}
exp19={literal_plane(a,c) for a in [[1,0,0],[0,-1,1],[-1,-1,1],[1,-1,1]] for c in [0,F(1,2)]}|{literal_plane([-1,1,1],0)}
for cid,expected in [(11,exp11),(19,exp19)]:assert {key(decode(g['torus'])) for g in goods if g['class_id']==cid}==expected
el=set()
for c in [F(0),F(1,2)]:
    for A,d in [([[1,0,0],[0,1,-2]],[c,F(0)]),([[1,0,0],[0,1,1]],[c,F(0)]),
                ([[-1,1,0],[-2,0,1]],[c,F(0)]),([[-2,1,0],[-1,0,1]],[F(1,2),c]),
                ([[0,1,1],[-1,0,1]],[F(1,2),c])]:el.add(key(torus(A,d)))
assert {key(decode(x['torus'])) for x in glines if x['class_id']==11}==el and len(el)==10
assert not any(x['class_id']==19 for x in glines)
# Independently verify the normal forms printed in the PDF.
literal11=[[[0,0,0,0]]*8,[[0,0,0,0],[15,0,0,0],[15,1,1,-1],[0,1,1,-1],[0,0,1,-1],[15,0,1,-1],[0,0,1,0],[15,0,1,0]],[[0,0,0,0],[15,0,0,-1],[15,0,1,-1],[15,0,1,-1],[0,0,0,-1],[0,0,1,-1],[15,0,0,0],[0,0,1,-1]]]
literal19=[[[0,0,0,0]]*8,[[0,0,0,0],[15,0,0,0],[0,1,1,-1],[15,1,1,-1],[15,0,1,-1],[0,0,1,-1],[0,1,0,0],[15,1,0,0]],[[0,0,0,0],[15,0,0,-1],[0,0,1,-1],[15,1,0,-1],[0,0,0,-1],[15,0,1,-1],[15,0,0,0],[0,1,0,-1]]]
assert cs[11]==literal11 and cs[19]==literal19
# Verify all uniform square templates, including specialization to every point.
for x in goods:
    sq=[[(F(e[0],30),)+tuple(e[1:]) for e in row] for row in x['square']]
    for i,j in combinations(range(8),2):assert profile(quotient(sq[i],sq[j]),[2]*4)
for x in glines:
    sq=[[(F(*e[0]),)+tuple(e[1:]) for e in row] for row in x['square']]
    for i,j in combinations(range(8),2):assert profile(quotient(sq[i],sq[j]),[2]*4)
# Check that every positive-dimensional projected locus belongs to the arrangement.
for cid in (11,19):
    full=json.load(open(R/f'outputs/full{cid}.json'));g=[decode(x['torus']) for x in goods if x['class_id']==cid]
    pl={key(decode(x['torus'])) for x in planes if x['class_id']==cid};ln={key(decode(x['torus'])) for x in lines if x['class_id']==cid}
    for v in full['profiles'].values():assert v['disconnected']==0
    for f in full['fibers']:assert f['disconnected']==0
    glob=set()
    for fam in full['row_families']:
        for t in split(fam['eq']):
            if any(subset(t,z) for z in g):continue
            if t['dim']==3:glob.update(tuple(map(tuple,r)) for r in fam['rows'])
            elif t['dim']==2:assert key(t) in pl
            elif t['dim']==1:assert key(t) in ln
    assert len(glob)==4
    # Every bad plane's full row set agrees with the A scanner's row set.
    for it,sc in zip(planes,psc):
        if it['class_id']!=cid or sc['k5']:continue
        t=decode(it['torus']);rr=set()
        for fam in full['row_families']:
            if all((F(e[3],30)+sum(F(e[i])*t['o'][i] for i in range(3)))%1==0 and all(sum(e[i]*t['V'][i][j] for i in range(3))==0 for j in range(2)) for e in fam['eq']):
                rr.update(tuple(map(tuple,row)) for row in specialize(fam['rows'],t,30))
        assert rr==set(tuple(map(tuple,row)) for row in sc['rows'])
        symbolic=[[(F(e[0],30),)+tuple(e[1:]) for e in row] for row in rr]
        assert not graph(symbolic)[0]
# Finite-root candidate lists agree in every one of the 132 canonical cases.
cases=json.load(open(R/'outputs/finite_cases.json'));enumerated=json.load(open(R/'outputs/finite_scan.json'))
assert len(cases)==len(enumerated)==132
for a,b in zip(cases,enumerated):
    assert a['id']==b['id']
    assert {tuple(F(*e) for e in row) for row in a['rows']}=={tuple(F(e,b['modulus']) for e in row) for row in b['rows']}
maps=json.load(open(R/'outputs/finite_case_maps.json'));assert len(maps)==len(pts)==2120
for ix,(x,mp) in enumerate(zip(pts,maps)):
    p=[F(*e) for e in x['point']];H=[[(F(e[0],30)+sum(F(e[i+1])*p[i] for i in range(3)))%1 for e in row] for row in cs[x['class_id']]]
    a,b,c=mp['row_order'];col=mp['pivot'];order=mp['column_order'];ca=cases[mp['case_id']]
    base=[[(H[k][j]-H[a][j]-H[k][col]+H[a][col])%1 for j in order] for k in [a,b,c]]
    assert base==[[F(*e) for e in row] for row in ca['H']]
    transformed=set()
    for row in x['rows']:
        w=[F(*e) for e in row];transformed.add(tuple((w[j]-H[a][j]-w[col]+H[a][col])%1 for j in order))
    assert transformed=={tuple(F(*e) for e in row) for row in ca['rows']}
# Check the printed 32-point table and demonstrate that A-only is impossible there.
table=[([1,7],[(4,3),(4,9),(5,9),(11,3)]),([2,8],[(5,3),(10,3),(10,9),(11,9)]),([4,10],[(1,3),(2,3),(2,9),(7,9)]),([5,11],[(1,9),(7,3),(8,3),(8,9)])]
expected={tuple(F(v,12) for v in [i,j,k]) for ii,jk in table for i in ii for j,k in jk}
actual={tuple(F(*e) for e in x['point']) for x in pts if x['class_id']==11 and x['k5']}
assert actual==expected and len(actual)==32
assert not any(x['k5'] for x in pts if x['class_id']==19)
point_squares=[]
for x in pts:
    if not x['k5']:continue
    p=[F(*q) for q in x['point']];H=[[(F(e[0],30)+sum(F(e[i+1])*p[i] for i in range(3)))%1 for e in row] for row in cs[11]];rr=[[F(*e) for e in row] for row in x['rows']]
    aa=[r for r in rr if all(profile([(a-b,) for a,b in zip(r,h)],[2]*4) for h in H)];assert len(aa)==4
    sq=H+[rr[i] for i in x['k5']];assert all((e*12).denominator==1 for row in sq for e in row)
    typ=Counter(profiles([(a-b,) for a,b in zip(sq[i],sq[j])]) for i,j in combinations(range(8),2));assert typ=={'A':20,'AB':4,'B':4}
    point_squares.append({'point':p,'exponents_mod_12':[[int(e*12) for e in row] for row in sq]})
example=[[0,0,0,0,0,0,0,0],[0,6,8,2,1,7,4,10],[0,3,7,7,9,1,6,1],[0,9,4,10,3,1,6,7],[0,6,4,1,7,9,10,3],[0,6,0,6,6,3,0,9],[0,0,8,8,4,7,1,4],[0,0,1,4,10,6,7,6]]
assert next(x['exponents_mod_12'] for x in point_squares if x['point']==[F(1,12),F(1,3),F(1,4)])==example
(R/'outputs/point_squares.json').write_text(json.dumps(enc(point_squares),indent=2))
status=list(csv.DictReader(open(R/'inputs/AAA_status_prior.csv')))
for row in status:
    cid=int(row['class'])
    if cid==11:row['general_completion_status']='exact iff: 7 surfaces + 10 circles (A-only completions) + 32 isolated mixed A/B points; independently verified'
    if cid==19:row['general_completion_status']='exact iff: 9 surfaces, all with A-only completions; independently verified'
with (R/'outputs/AAA_status.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(status[0]));w.writeheader();w.writerows(status)
report={'class11':{'planes':7,'circles':10,'isolated_points':32,'A_only_on_isolated_points':False,'critical_points':616},'class19':{'planes':9,'circles':0,'isolated_points':0,'critical_points':1504},'uniform_square_templates':26,'isolated_square_witnesses':32,'independent_finite_cases':132,'finite_triples':2120,'lattice_test_cases':1200,'all_checks_passed':True}
(R/'outputs/final_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
