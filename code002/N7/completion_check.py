"""Independent direct verification of the residual completion certificates.

Uses only JSON, fractions, and finite combinatorics. Does not import either
congruence solver, chart canonicalizer, recognition search, or primary code.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from collections import defaultdict
import json
ROOT=Path(__file__).resolve().parent

def read(p):return json.loads((ROOT/p).read_text())
def exponent(x):
    v=6*(Q(x)%1);assert v.denominator==1
    return int(v)
def reduce_polynomial(c):
    c=list(c)
    for j in range(len(c)-1,1,-1):
        c[j-1]+=c[j];c[j-2]-=c[j];c[j]=0
    return tuple(c[:2])
def product_polynomial(a,b):
    out=defaultdict(lambda:[0]*6)
    for (e,n),(f,m) in zip(a,b):out[n-m][(e-f)%6]+=1
    return {n:reduce_polynomial(c) for n,c in out.items() if reduce_polynomial(c)!=(0,0)}
def cycles(a,b):
    q=[((e-f)%6,n-m) for (e,n),(f,m) in zip(a,b)]
    for tri in combinations(range(7),3):
        e,n=q[tri[0]]
        if any(q[j][1]!=n for j in tri) or {(q[j][0]-e)%6 for j in tri}!={0,2,4}:continue
        rest=[j for j in range(7) if j not in tri]
        for j in rest[1:]:
            other=[k for k in rest if k not in (rest[0],j)]
            if all(q[u][1]==q[v][1] and (q[u][0]-q[v][0])%6==3 for u,v in [(rest[0],j),tuple(other)]):return True
    return False
def matrix(record):
    return [[(exponent(c),b[0] if b else 0) for c,b in zip(rc,rb)]
            for rc,rb in zip(record['constants'],record['basis'])]
def verify_square(M):
    assert len(M)==7 and all(len(row)==7 for row in M)
    for i,j in combinations(range(7),2):
        assert not product_polynomial(M[i],M[j]),(i,j)
        assert cycles(M[i],M[j]),(i,j)
def witness(M,w):
    rows=[r-1 for r in w['rows']];cols=[j-1 for j in w['columns']]
    assert len(set(rows))==3 and all(0<=i<len(M) for i in rows)
    assert sorted(cols)==list(range(7))
    ac=exponent(w['parameter_constant']);ab=w['parameter_basis'];slope=ab[0] if ab else 0
    E=[[0]*7,[0,2,4,2,5,0,3],[0,2,4,0,2,3,5]]
    D=[[0]*7,[0,0,0,1,1,1,1],[0,0,0,1,1,1,1]]
    for h,r in enumerate(rows):
        for k,j in enumerate(cols):
            ee=(M[r][j][0]-M[rows[0]][j][0]-M[r][cols[0]][0]+M[rows[0]][cols[0]][0])%6
            nn=M[r][j][1]-M[rows[0]][j][1]-M[r][cols[0]][1]+M[rows[0]][cols[0]][1]
            assert ee==(E[h][k]+D[h][k]*ac)%6 and nn==D[h][k]*slope

def point_rows(raw,theta):
    rows=set()
    for rec in raw['maximal_components']:
        c=[Q(x) for x in rec['constant']];B=rec['basis']
        if rec['dimension']==0:
            if tuple(c[:2])==theta:rows.add(tuple(c[2:]))
        else:
            assert rec['dimension']==1
            vary=[i for i in range(2) if B[i][0]];assert vary
            h=vary[0];assert abs(B[h][0])==1
            t=((theta[h]-c[h])/B[h][0])%1
            if all((c[i]+B[i][0]*t-theta[i])%1==0 for i in range(2)):
                rows.add(tuple((c[i]+B[i][0]*t)%1 for i in range(2,9)))
    return rows

def seed_rows(raw,theta):
    s=raw['family']
    return [[(exponent(Q(e,6)+sum(x*y for x,y in zip(b,theta))),0) for e,b in zip(E,B)]
            for E,B in zip(s['constants'],s['basis'])]

def check_f2():
    raw=read('input/F2.json');out=read('results/F2_completions.json')
    circles=out['circles'];assert len(circles)==8
    consumed=[]
    for circle in circles:
        M=matrix(circle['matrix']);verify_square(M);witness(M,circle['f4_witness'])
        base=circle['base'];bc=[Q(x) for x in base['constant']];BB=base['basis']
        # Confirm every appended row and original seed is exactly the census input.
        for h,(E,B) in enumerate(zip(raw['family']['constants'],raw['family']['basis'])):
            expected=[(exponent(Q(e,6)+sum(x*y for x,y in zip(b,bc))),sum(b[k]*BB[k][0] for k in range(2))) for e,b in zip(E,B)]
            assert M[h]==expected
        for h,idx in enumerate(circle['fourth_row_ids'],3):
            rec=raw['maximal_components'][idx-1];assert rec['dimension']==1
            assert rec['constant'][:2]==base['constant'] and rec['basis'][:2]==base['basis']
            assert M[h]==[(exponent(c),b[0]) for c,b in zip(rec['constant'][2:],rec['basis'][2:])]
            consumed.append(idx)
    assert sorted(consumed)==[i for i,r in enumerate(raw['maximal_components'],1) if r['dimension']==1]
    assert all(r['dimension'] in [0,1] for r in raw['maximal_components'])
    # Exceptional seed values: every point component plus every circle crossing.
    exceptional={tuple(Q(c) for c in r['constant'][:2]) for r in raw['maximal_components'] if r['dimension']==0}
    for a,b in combinations(circles,2):
        aa=a['base'];bb=b['base'];ia=next(i for i in range(2) if aa['basis'][i][0]==0);ib=next(i for i in range(2) if bb['basis'][i][0]==0)
        if ia!=ib:
            p=[None,None];p[ia]=Q(aa['constant'][ia]);p[ib]=Q(bb['constant'][ib]);exceptional.add(tuple(p))
        else:assert aa['constant'][ia]!=bb['constant'][ib]
    assert exceptional=={tuple(Q(x) for x in r['theta']) for r in out['exceptional_seeds']}
    assert len(exceptional)==32
    total=0
    for r in out['exceptional_seeds']:
        theta=tuple(Q(t) for t in r['theta']);rows=[tuple(Q(v) for v in row) for row in r['rows']]
        assert len(rows)==len(set(rows)) and set(rows)==point_rows(raw,theta)
        M=[[(exponent(v),0) for v in row] for row in rows]
        edges={p for p in combinations(range(len(M)),2) if not product_polynomial(M[p[0]],M[p[1]])}
        assert edges=={tuple(i-1 for i in p) for p in r['edges']}
        actual={p for p in combinations(range(len(M)),4) if all(q in edges for q in combinations(p,2))}
        assert actual=={tuple(i-1 for i in c['selected_rows']) for c in r['completions']}
        assert edges=={q for p in actual for q in combinations(p,2)}
        for comp in r['completions']:
            selected=[i-1 for i in comp['selected_rows']];full=seed_rows(raw,theta)+[M[i] for i in selected]
            verify_square(full);witness(full,comp['f4_witness'])
            chosen={tuple(v[0] for v in M[i]) for i in selected}
            for label in comp['circle_specializations']:
                circ=circles[label-1];base=circ['base'];free=base['pivot'][0];t=theta[free]
                assert all((Q(c)+b[0]*t-theta[h])%1==0 for h,(c,b) in enumerate(zip(base['constant'],base['basis'])))
                cm=matrix(circ['matrix']);sp={tuple(exponent(Q(e,6)+n*t) for e,n in row) for row in cm[3:]}
                assert sp==chosen
        total+=len(actual)
    assert total==40
    print('Independent F2 check: eight circle completions; all 32 exceptional fibres; 40 finite completions; all are circle specializations.')

def on_circle(point,rec):
    c=[Q(x) for x in rec['constant']];B=rec['basis'];h=next(i for i in range(2,9) if B[i][0])
    assert abs(B[h][0])==1
    t=((point[h]-c[h])/B[h][0])%1
    return all((c[i]+B[i][0]*t-point[i])%1==0 for i in range(9))

def check_f3():
    raw=read('input/F3.json');out=read('results/F3_completions.json');assert len(out)==2
    for r in out:
        q=Q(r['theta']);theta=(q,q);vertical=[];points=[]
        for rec in raw['maximal_components']:
            c=[Q(x) for x in rec['constant']];B=rec['basis']
            if rec['dimension']==0:
                if tuple(c[:2])==theta:points.append(c)
            else:
                assert rec['dimension']==1
                vary=[i for i in range(2) if B[i][0]]
                if not vary:
                    if tuple(c[:2])==theta:vertical.append(rec)
                else:
                    h=vary[0];assert abs(B[h][0])==1;t=((theta[h]-c[h])/B[h][0])%1
                    p=[(c[i]+B[i][0]*t)%1 for i in range(9)]
                    if tuple(p[:2])==theta:points.append(p)
        assert len(vertical)==4
        key=lambda rec:(tuple(rec['constant']),tuple(tuple(b) for b in rec['basis']))
        assert {key(c) for c in vertical}=={key(c) for c in r['row_circles']}
        assert all(any(on_circle(p,c) for c in vertical) for p in points)
        # The first three pair equations force t2=-t1, t3=-alpha/t1,
        # t4=alpha/t1. Their coefficients are nonzero, and no division by
        # a phase variable loses any unit-circle solution.
        assert len(r['pair_grams'])==6
        eq={tuple(p['rows']):(p['equation'],Q(p['rhs'])) for p in r['pair_grams']}
        A,b=eq[(1,2)];assert A==[-1,1,0,0] and b==Q(1,2)
        A,b=eq[(1,3)];assert A==[-1,0,-1,0] and b==(-q-Q(1,2))%1
        A,b=eq[(1,4)];assert A==[-1,0,0,-1] and b==(-q)%1
        expected_c=(0,Q(1,2),(q+Q(1,2))%1,q);expected_B=(1,1,-1,-1)
        sol=r['solutions'];assert len(sol)==1
        assert tuple(Q(x) for x in sol[0]['constant'])==expected_c and tuple(b[0] for b in sol[0]['basis'])==expected_B
        # Recompute all Gram binomials without importing the primary arithmetic.
        for pair in r['pair_grams']:
            i,j=[h-1 for h in pair['rows']];left=r['row_circles'][i];right=r['row_circles'][j]
            terms=defaultdict(lambda:[0]*6)
            for col in range(2,9):
                degree=[0]*4;degree[i]=left['basis'][col][0];degree[j]-=right['basis'][col][0]
                e=exponent(Q(left['constant'][col])-Q(right['constant'][col]));terms[tuple(degree)][e]+=1
            reduced={k:reduce_polynomial(c) for k,c in terms.items() if reduce_polynomial(c)!=(0,0)}
            assert reduced=={tuple(t['degree']):tuple(t['coefficient']) for t in pair['terms']}
            assert len(pair['terms'])==2
            u,v=pair['terms'];assert pair['equation']==[a-b for a,b in zip(u['degree'],v['degree'])]
            assert tuple(u['coefficient'])!=(0,0)
            k=exponent(pair['rhs']);p=[0]*6
            for h,c in enumerate(u['coefficient']):p[(h+k)%6]+=c
            for h,c in enumerate(v['coefficient']):p[h]+=c
            assert reduce_polynomial(p)==(0,0)
        for rec in r['row_circles']:
            assert sorted(b[0] for b in rec['basis'][2:])==[0,0,0,0,0,1,1]
        M=matrix(r['matrix']);verify_square(M);witness(M,r['f4_witness'])
        assert M[:3]==seed_rows(raw,theta)
        for i,circ in enumerate(r['row_circles']):
            expected=[(exponent(Q(c)+b[0]*expected_c[i]),b[0]*expected_B[i]) for c,b in zip(circ['constant'][2:],circ['basis'][2:])]
            assert M[3+i]==expected
    print('Independent F3 check: both complete fibres, all Gram equations, the unique completion circles, and F4 certificates.')

if __name__=='__main__':
    check_f2();check_f3()
    (ROOT/'results/independent_verification.json').write_text(json.dumps(dict(F2=True,F3=True,imports_no_primary_solver=True,all_passed=True),indent=2)+'\n')
