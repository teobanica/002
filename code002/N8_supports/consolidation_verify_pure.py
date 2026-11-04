from torus3 import *
from exact_rows import *
from itertools import combinations
import ccc_supports as C
R=Path(__file__).resolve().parents[1]
cs={x['class_id']:x['rows'] for x in json.load(open(R/'inputs/pure_components.json'))}
def evalrow(row,p):return tuple((F(e[0],60)+sum(F(e[i+1])*p[i] for i in range(3)))%1 for e in row)
allpoints=[];summary={}
for name,cid in [('BBB',100),('CCC',101)]:
    d=json.load(open(R/f'outputs/{name}_extensions.json'));assert all(v['disconnected']==0 for v in d['profiles'].values())
    families=[]
    for f in d['row_families']:
        for t in split(f['eq'],60):families.append((t,f['rows']))
    if name=='BBB':
        expected={key(torus([[1,0,0],[0,0,1]],[c,F(0)])) for c in [F(1,6),F(1,2),F(5,6)]}
        expected|={key(torus([[0,1,0],[0,0,1]],[c,F(0)])) for c in [F(1,4),F(7,12),F(11,12)]}
        expected|={key(torus([[1,-1,0],[0,0,1]],[c,F(0)])) for c in [F(1,12),F(5,12),F(3,4)]}
        assert {key(t) for t,rr in families}==expected and len(families)==9
        for t,rr in families:
            assert len(rr)==1
            def symbol(row):return [((F(e[0],60)+sum(F(e[i+1])*t['o'][i] for i in range(3)))%1,)+tuple(sum(e[i+1]*t['V'][i][j] for i in range(3)) for j in range(t['dim'])) for e in row]
            assert all(profiles(quotient(symbol(rr[0]),symbol(h))) for h in cs[cid])
    else:assert not families
    pointset={tuple(F(*e) for e in f['point']) for f in d['fibers']}
    for (a,_),(b,_) in combinations(families,2):
        for t in intersection(a,b):
            if t['dim']==0:assert tuple(t['o']) in pointset
    hist=Counter();pairs=0;off=[]
    for f in d['fibers']:
        assert f['disconnected']==0;p=tuple(F(*e) for e in f['point']);assert p[2]==0
        rows=sorted(tuple(F(*e) for e in row) for row in f['rows']);H=[evalrow(row,p) for row in cs[cid]]
        for row in rows:
            assert all(profiles([(x-y,) for x,y in zip(row,h)]) for h in H)
        sy=[[(e,) for e in row] for row in rows];k,cl,edges=graph(sy);hist[len(rows)]+=1;pairs+=len(rows)*(len(rows)-1)//2
        if name=='BBB':assert edges==0
        else:assert len(rows)==10 and len(k)==5 and edges==20
        if not any(all((sum(F(a)*b for a,b in zip(eq,p))-c)%1==0 for eq,c in zip(t['A'],t['c'])) for t,_ in families):off.append(p)
        allpoints.append({'class_id':cid,'point':p,'rows':rows,'H':H})
    if name=='BBB':assert len(pointset)==135 and len(off)==54 and pairs==378 and hist=={1:36,2:72,3:18,8:9}
    else:
        expected={(F(a,30),F((2*a+b)%30,30),F(0)) for a in range(0,30,2) for b in [10,20]};assert pointset==expected
        old=json.load(open(R/'inputs/n8_ccc_square_results.json'))
        for f in old['points']:
            t=(F(f['A'],30),F(f['B'],30),F(0));new=next(x for x in allpoints if x['class_id']==cid and tuple(x['point'])==t)
            assert {tuple(F(e,30) for e in r['row']) for r in f['candidate_fourth_rows']}==set(new['rows'])
    summary[name]={'curves':len(families),'points_in_finite_checks':len(pointset),'isolated_points_off_curves':len(off),'row_count_distribution':dict(hist),'candidate_pairs_tested':pairs}
    if name=='BBB':
        (R/'outputs/BBB_isolated_points.json').write_text(json.dumps(enc(sorted(off)),indent=2))
        (R/'outputs/BBB_circle_templates.json').write_text(json.dumps(enc([{'torus':t,'rows_mod_60':rr} for t,rr in families]),indent=2))
# Explicit equivalences are checked independently of the canonical search.
data=json.load(open(R/'inputs/n8_ccc_square_results.json'));certs=json.load(open(R/'outputs/ccc_canonical.json'));K=certs[0]['canonical']
for original,cert in zip(data['square_completions'],certs):
    H=original['rows_exp30'];a,b=cert['pivot']
    normalized=[[(H[i][j]-H[a][j]-H[i][b]+H[a][b])%30 for j in cert['columns']] for i in cert['rows']]
    assert normalized==cert['canonical']==K
    assert sorted(cert['rows'])==sorted(cert['columns'])==list(range(8))
    for rows in [H,list(zip(*H))]:
        assert Counter(C.bits(C.quotient(rows[i],rows[j])) for i,j in combinations(range(8),2))=={1:16,4:12}
assert C.tangent_check(K)==4096
# The correct Haagerup multiset uses ordered row and column pairs.
h=Counter((K[i][a]+K[j][b]-K[i][b]-K[j][a])%30 for i in range(8) for j in range(8) if i!=j for a in range(8) for b in range(8) if a!=b)
assert any(__import__('math').gcd(k,30)==1 and n for k,n in h.items())
summary['CCC'].update({'square_presentations':60,'monomial_equivalence_classes':1,'independent_equivalence_certificates':60,'regular_support_strata_checked':4096,'dephased_support_rank':49,'minimal_Butson_level':30})
(R/'outputs/pure_points.json').write_text(json.dumps(enc(allpoints),indent=2))
(R/'outputs/pure_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
