"""Exact affine phases, cycle witnesses, and positive F4 certificates."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations, permutations, product
from collections import defaultdict
import sys,json
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'engine'))
from families import family
from census import PAT,canonical,deserialize,serialize,contains,solve_cached
from independent import canon
from analyze import POWERS,valid_patterns

def maximal(cs):
    out=[]
    for dim in sorted({len(c[2]) for c in cs},reverse=True):
        out.extend(c for c in cs if len(c[2])==dim and not any(contains(b,c) for b in out))
    return sorted(out,key=lambda c:(-len(c[2]),c))

def restrict(comp,indices,values):
    c,B,p=comp;k=len(p)
    A=tuple(B[i] for i in indices);rhs=tuple((v-c[i])%1 for i,v in zip(indices,values))
    for u,V in solve_cached(A,rhs,k):
        m=len(V[0]) if V else 0
        cc=tuple((v+sum(a*b for a,b in zip(row,u)))%1 for v,row in zip(c,B))
        BB=tuple(tuple(sum(row[h]*V[h][j] for h in range(k)) for j in range(m)) for row in B)
        yield canon(cc,BB)

def data(name):
    return json.loads((ROOT/'input'/(name+'.json')).read_text())

def seed_rows(name,pc,PB):
    seed=family(name);k=len(PB[0])
    cs=[];Bs=[]
    for rc,rb in zip(seed['constants'],seed['basis']):
        cs.append(tuple((Q(e,6)+sum(a*b for a,b in zip(v,pc)))%1 for e,v in zip(rc,rb)))
        Bs.append(tuple(tuple(sum(v[h]*PB[h][j] for h in range(len(pc))) for j in range(k)) for v in rb))
    return cs,Bs

def fibre(name,theta):
    d=len(theta)
    return maximal({f for rec in data(name)['maximal_components']
                    for f in restrict(deserialize(rec),range(d),theta)})

def gram(c,B,e,D):
    out=defaultdict(lambda:[0,0])
    for a,b,u,v in zip(c,B,e,D):
        z=6*((a-u)%1);assert z.denominator==1
        p=POWERS[int(z)]
        deg=tuple(x-y for x,y in zip(b,v))
        out[deg][0]+=p[0];out[deg][1]+=p[1]
    return {k:tuple(v) for k,v in sorted(out.items()) if v!=[0,0]}

def verify_matrix(cs,Bs):
    witnesses=[]
    for i,j in combinations(range(len(cs)),2):
        assert not gram(cs[i],Bs[i],cs[j],Bs[j]),(i,j)
        c=tuple((a-b)%1 for a,b in zip(cs[i],cs[j]));B=tuple(tuple(x-y for x,y in zip(a,b)) for a,b in zip(Bs[i],Bs[j]))
        pat=valid_patterns(c,B);assert pat,(i,j)
        witnesses.append(dict(rows=[i+1,j+1],cycle_edges=pat[0]))
    return witnesses

def verify_f4(cs,Bs,w):
    target=family('F4');rows=[r-1 for r in w['rows']];cols=[c-1 for c in w['columns']]
    ac=Q(w['parameter_constant']);ab=w['parameter_basis'];r0=rows[0];j0=cols[0]
    for h,r in enumerate(rows):
        for k,j in enumerate(cols):
            c=(cs[r][j]-cs[r0][j]-cs[r][j0]+cs[r0][j0])%1
            B=tuple(Bs[r][j][v]-Bs[r0][j][v]-Bs[r][j0][v]+Bs[r0][j0][v] for v in range(len(ab)))
            e=Q(target['constants'][h][k],6);s=target['basis'][h][k][0]
            assert (c-e-s*ac)%1==0 and B==tuple(s*b for b in ab)

def f4_match(cs,Bs):
    target=family('F4');tc=tuple(Q(v,6) for row in target['constants'][1:] for v in row)
    tb=tuple(row[0] for rr in target['basis'][1:] for row in rr)
    for base,u,v in permutations(range(len(cs)),3):
        cc=[tuple((cs[r][j]-cs[base][j])%1 for j in range(7)) for r in (u,v)]
        bb=[tuple(tuple(x-y for x,y in zip(Bs[r][j],Bs[base][j])) for j in range(7)) for r in (u,v)]
        for pat in valid_patterns(cc[0],bb[0]):
            tri=(pat[0][1],pat[0][0],pat[1][0])
            if pat[0][2]==4:tri=(tri[0],tri[2],tri[1])
            P=(pat[2][1],pat[2][0]);S=(pat[3][1],pat[3][0])
            for shift,swap,rp,rs in product(range(3),range(2),range(2),range(2)):
                TT=tri[shift:]+tri[:shift];PP=P[::-1] if rp else P;SS=S[::-1] if rs else S
                if swap:PP,SS=SS,PP
                cols=TT+PP+SS;a=cols[0]
                nc=tuple((cc[r][j]-cc[r][a])%1 for r in range(2) for j in cols)
                nb=tuple(tuple(x-y for x,y in zip(bb[r][j],bb[r][a])) for r in range(2) for j in cols)
                ac=nc[5];ab=nb[5]
                if all((nc[i]-tc[i]-tb[i]*ac)%1==0 and nb[i]==tuple(tb[i]*x for x in ab) for i in range(14)):
                    w=dict(rows=[base+1,u+1,v+1],columns=[j+1 for j in cols],parameter_constant=str(ac),parameter_basis=ab)
                    verify_f4(cs,Bs,w);return w
    return None

def serialize_matrix(cs,Bs):
    return dict(constants=[[str(c) for c in row] for row in cs],basis=Bs)
