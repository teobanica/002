"""Independent exact checks and completion coverage for maximal generic AAA.

Rows are additive phases (constant in Q/Z, followed by formal exponents).
No specialization to numerical roots or floating point is used.
"""
from pathlib import Path
import json,sys
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter
import sympy as S
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from exact_rows import profiles,quotient

R=Path(__file__).resolve().parents[1]
def read(n):return json.loads((R/'outputs'/n).read_text())
def decode(row):return tuple((F(e[0],30),)+tuple(map(F,e[1:])) for e in row)
def pad(row,n):return tuple(e+(F(0),)*n for e in row)
def enc(o):
    if isinstance(o,F):return [o.numerator,o.denominator]
    if isinstance(o,(list,tuple,set)):return [enc(x) for x in sorted(o)] if isinstance(o,set) else [enc(x) for x in o]
    if isinstance(o,dict):return {k:enc(v) for k,v in o.items()}
    return o
def substitute(row,base,values,extra=0):
    result=[]
    for e in row:
        v=list(e[:base+1]) + [F(0)]*extra
        for k,value in enumerate(values):
            for j,x in enumerate(value):v[j]+=e[base+1+k]*x
        v[0]%=1;result.append(tuple(v))
    return tuple(result)
def parameter_for(row,target,base):
    j=next(j for j,e in enumerate(row) if e[base+1])
    k=row[j][base+1]
    assert k==1
    value=tuple((target[j][i]-row[j][i]) if i else (target[j][0]-row[j][0])%1 for i in range(base+1))
    return value if substitute(row,base,[value])==target else None
def solutions(eq,f,base):
    if not eq:
        return [[(F(0),)*(base+1)+tuple(F(i==j) for j in range(f)) for i in range(f)]]
    A=S.Matrix([e[:f] for e in eq]);B=S.Matrix([e[f:f+base]+[S.Rational(e[-1],30)] for e in eq])
    sd,ud,vd=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ))
    D,U,V=[x.to_Matrix() for x in (sd,ud,vd)];assert U*A*V==D
    rank=A.rank();C=-U*B
    assert all(x==0 for x in C[rank:,:])
    out=[]
    for ks in product(*(range(abs(int(D[i,i]))) for i in range(rank))):
        Y=S.zeros(f,base+1+f-rank)
        for i in range(rank):
            for j in range(base+1):Y[i,j]=(C[i,j]+(ks[i] if j==base else 0))/D[i,i]
        for i in range(rank,f):Y[i,base+1+i-rank]=1
        X=V*Y
        out.append([tuple([F(X[i,base])%1]+[F(X[i,j]) for j in range(base)]+[F(X[i,j]) for j in range(base+1,X.cols)]) for i in range(f)])
    return out
def cliques(adj,need):
    out=[]
    def rec(xs,sol):
        if len(sol)==need:out.append(sol);return
        while len(xs)+len(sol)>=need:
            i=xs[0];xs=xs[1:];rec([j for j in xs if adj[i][j]],sol+[i])
    rec(list(range(len(adj))),[]);return out

components=json.loads((R/'inputs/components.json').read_text())
pairs=read('maximal_pairs.json');summaries=[];templates=[]
for d in pairs:
    cid=d['class_id'];full=read(f'generic{cid}.json')[0];base=full['dim'];assert base==4
    fs=full['families'];rows=[decode(x['row']) for x in fs]
    fixed=[i for i,x in enumerate(fs) if x['free']==0];free=[i for i,x in enumerate(fs) if x['free']==1]
    fixedset={rows[i] for i in fixed};assert len(fixedset)==len(fixed)
    H=[decode(r) for r in next(c['rows'] for c in components if c['class_id']==cid)]
    for i,x in enumerate(fs):
        assert x['profiles']==['AAA']
        assert all('A' in profiles(quotient(rows[i],pad(h,x['free']))) for h in H)
    adj=[[False]*len(fixed) for _ in fixed]
    for i,j in combinations(range(len(fixed)),2):
        p=profiles(quotient(rows[fixed[i]],rows[fixed[j]]));assert p in ('','A')
        adj[i][j]=adj[j][i]=bool(p)
    kk=[[fixed[i] for i in x] for x in cliques(adj,5)]
    assert len(kk)==[21,8][cid]
    assert sum(map(sum,adj))//2==[150,65][cid]
    bypair={(p['i'],p['j']):p for p in d['pairs']}
    cplusplus_edges={(p['i'],p['j']) for p in d['pairs'] if p['free']==0}
    assert cplusplus_edges=={(fixed[i],fixed[j]) for i,j in combinations(range(len(fixed)),2) if adj[i][j]}
    family_info=[]
    for i in free:
        universal=[j for j in fixed if 'A' in profiles(quotient(rows[i],pad(rows[j],1)))]
        assert len(universal)==3
        assert all('A' in profiles(quotient(rows[j],rows[k])) for j,k in combinations(universal,2))
        half=(F(1,2),)+(F(0),)*base+(F(1),)
        negrow=substitute(rows[i],base,[half],1)
        square=[pad(h,1) for h in H]+[pad(rows[j],1) for j in universal]+[rows[i],negrow]
        assert all('A' in profiles(quotient(a,b)) for a,b in combinations(square,2))
        columns=list(zip(*square))
        assert all('A' in profiles(quotient(a,b)) for a,b in combinations(columns,2))
        covered=[]
        for ki,k5 in enumerate(kk):
            for j in k5:
                value=parameter_for(rows[i],rows[j],base)
                if value is None:continue
                value2=((value[0]+F(1,2))%1,)+value[1:]
                if set([rows[z] for z in universal]+[substitute(rows[i],base,[value]),substitute(rows[i],base,[value2])])=={rows[z] for z in k5}:
                    covered.append({'clique':ki,'parameter':value});break
        family_info.append({'row_family':i,'triangle':universal,'connected_clique_coverage':covered})
        templates.append({'class_id':cid,'row_family':i,'base_dimension':base,'square':square})
    assert set(c['clique'] for f in family_info for c in f['connected_clique_coverage'])==set(range(len(kk)))
    relation_components=0
    critical_parameters={i:set() for i in free}
    for p in d['pairs']:
        if not p['free']:continue
        i,j=p['i'],p['j'];nf=p['free'];fi=fs[i]['free'];fj=fs[j]['free']
        for rel in p['relations']:
            assert rel['profiles']=='A'
            for vals in solutions(rel['eq'],nf,base):
                extra=len(vals[0])-base-1;relation_components+=1
                a=substitute(rows[i],base,vals[:fi],extra)
                b=substitute(rows[j],base,vals[fi:],extra)
                assert 'A' in profiles(quotient(a,b))
                if i==j:
                    assert fi==fj==1 and len(rel['eq'])==1
                    diff=tuple(x-y for x,y in zip(vals[0],vals[1]));assert diff[0]%1==F(1,2) and all(x==0 for x in diff[1:])
                elif nf==1:
                    if any(not z['eq'] for z in p['relations']):continue
                    varied=a if fi else b
                    assert extra==0
                    critical_parameters[i if fi else j].add(vals[0])
                else:
                    # Distinct free families: at least one row becomes one of
                    # the connected candidates, even on the torsion branches.
                    def fixed_candidate(r):return all(all(x==0 for x in e[base+1:]) for e in r) and tuple(e[:base+1] for e in r) in fixedset
                    assert fixed_candidate(a) or fixed_candidate(b),(cid,p,rel,vals)
                    for index,value in zip([i,j],vals):
                        if all(v==0 for v in value[base+1:]):critical_parameters[index].add(value[:base+1])
    augmented=set(fixedset)
    for i,values in critical_parameters.items():
        for value in list(values):values.add(((value[0]+F(1,2))%1,)+value[1:])
        for value in values:augmented.add(substitute(rows[i],base,[value]))
    augmented=sorted(augmented)
    aa=[[False]*len(augmented) for _ in augmented]
    for i,j in combinations(range(len(augmented)),2):aa[i][j]=aa[j][i]=bool(profiles(quotient(augmented[i],augmented[j])))
    augmented_cliques=cliques(aa,5)
    uncovered=[]
    for ki,k5 in enumerate(augmented_cliques):
        target={augmented[j] for j in k5};covered=False
        for info in family_info:
            i=info['row_family'];triangle=[rows[j] for j in info['triangle']]
            for j in k5:
                value=parameter_for(rows[i],augmented[j],base)
                if value is None:continue
                value2=((value[0]+F(1,2))%1,)+value[1:]
                if set(triangle+[substitute(rows[i],base,[value]),substitute(rows[i],base,[value2])])==target:covered=True
        if not covered:uncovered.append(k5)
    print('Critical-row enlargement:',len(fixed),'->',len(augmented),'rows;',len(augmented_cliques),'K5;',len(uncovered),'uncovered.')
    assert len(uncovered)==cid
    for k5 in uncovered:
        square=H+[augmented[j] for j in k5]
        assert all('A' in profiles(quotient(a,b)) for a,b in combinations(square,2))
        templates.append({'class_id':cid,'row_family':'additional_fixed_completion','base_dimension':base,'square':square})
    summaries.append({'class_id':cid,'connected_candidates':len(fixed),'free_row_families':len(free),'connected_edges':sum(map(sum,adj))//2,'connected_K5':len(kk),'cliques':kk,'completion_families':family_info,'checked_pair_components':relation_components,'augmented_candidates':len(augmented),'augmented_K5':len(augmented_cliques),'augmented_rows':augmented,'additional_cliques':uncovered,'critical_parameters':critical_parameters})
    print('PASS class',cid,': complete generic coverage by',len(free),'five-parameter families and',len(uncovered),'additional fixed completions; all pair relations independently verified.')
(R/'outputs/maximal_classification.json').write_text(json.dumps(enc(summaries),indent=2)+'\n')
(R/'outputs/maximal_square_templates.json').write_text(json.dumps(enc(templates),indent=2)+'\n')
