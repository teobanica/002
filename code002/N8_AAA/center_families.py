"""Choose simpler roots-of-unity centers and analyze the mixed circle exactly."""
from phase_algebra import *
from functools import lru_cache
for ti in [7,9]:
    t=read('ten_representatives.json')[ti];H=[[[F(*x) for x in e] for e in row] for row in t['square']];d=len(H[0][0])-1
    A=S.Matrix([e[1:] for row in H for e in row]);c=S.Matrix([S.Rational(e[0]) for row in H for e in row]);sd,ud,vd=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ));D,U,V=[x.to_Matrix() for x in [sd,ud,vd]]
    for m in [1,2,3,4,6,12]:
        b=-U*c*m
        if all(x.q==1 for x in b[d:,0]):
            vals=V*S.Matrix([b[i]/D[i,i]/m for i in range(d)])
            K=[[((e[0]+sum(F(vals[k])*e[k+1] for k in range(d)))%1,)+tuple(e[1:]) for e in row] for row in H]
            assert all(profiles(quotient(a,b)) for a,b in combinations(K,2))
            out={**t,'square':enc(K),'origin_parameter_shift':enc([F(x) for x in vals]),'center_level':m}
            (R/f'outputs/centered_family_{ti}.json').write_text(json.dumps(out,indent=2));break

@lru_cache(None)
def matchings(ix):
    if not ix:return ((),)
    a=ix[0];out=[]
    for b in ix[1:]:
        rest=tuple(i for i in ix if i not in [a,b])
        out.extend((((a,b),)+p) for p in matchings(rest))
    return tuple(out)

H=[[tuple(F(*x) for x in e) for e in row] for row in read('centered_family_9.json')['square']]
special=[]
for i,j in combinations(range(8),2):
    q=quotient(H[i],H[j])
    if 'A' in profiles(q):continue
    points=set()
    for p in matchings(tuple(range(8))):
        eq=[(q[b][1]-q[a][1],(F(1,2)-q[b][0]+q[a][0])%1) for a,b in p]
        if any(not k and c for k,c in eq):continue
        if not any(k for k,c in eq):raise AssertionError('Generic A was missed')
        k,c=next((k,c) for k,c in eq if k)
        assert k.denominator==1
        for n in range(abs(int(k))):
            t=((c+n)/k)%1
            if all((a*t-b)%1==0 for a,b in eq):points.add(t)
    special.append({'pair':[i,j],'A_parameters':enc(sorted(points))})
allpoints=set.intersection(*[{F(*p) for p in x['A_parameters']} for x in special])
assert allpoints=={F(i,12) for i in range(12)}
t=F(1,18);sq=[[(e[0]+e[1]*t)%1 for e in row] for row in H]
assert all((e*18).denominator==1 for row in sq for e in row)
prof=Counter(profiles([(a-b,) for a,b in zip(u,v)]) for u,v in combinations(sq,2));assert prof=={'A':24,'B':4}
assert any(e.denominator==18 for row in sq for e in row)
(R/'outputs/mixed_circle_analysis.json').write_text(json.dumps({'exceptional_pairs':special,'all_A_parameters':enc(sorted(allpoints)),'BH18_parameter':enc(t),'BH18_matrix':[[int(e*18) for e in row] for row in sq],'BH18_profiles':dict(prof)},indent=2))
print('Mixed circle A-only exactly at 12th roots; exact BH(8,18) mixed example certified.')
