import json, itertools, math
from fractions import Fraction
from collections import Counter, defaultdict
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'outputs/canonical.json'
OUT=ROOT/'outputs/affine.json'
D0=json.load(open(IN))
D=D0['triples']
Pstd={
 'A':[[0,1],[2,3],[4,5],[6,7]],
 'B':[[0,1,2],[3,4,5],[6,7]],
 'C':[[3,4,5,6,7],[0,1,2]],
}

def build(P,Q,R):
    parts=[P,Q,R]; signs=[1,-1,-1]; offs=[]; n=0
    for S in parts: offs.append(n); n+=len(S)
    A=[[0]*n for _ in range(8)]
    for S,off,sg in zip(parts,offs,signs):
        for bi,b in enumerate(S):
            for j in b: A[j][off+bi]=sg
    return A

def label_options(P,p):
    blocks=[b for b in P if len(b)==p]
    if not blocks: return [[0]*8]
    per=[]
    for b in blocks:
        opts=[]
        for rest in itertools.permutations(range(1,p)):
            v=[0]*8
            v[b[0]]=0
            for i in range(1,p): v[b[i]]=rest[i-1]
            opts.append(v)
        per.append(opts)
    ans=[[0]*8]
    for opts in per:
        ans2=[]
        for a in ans:
            for o in opts: ans2.append([(a[j]+o[j])%p for j in range(8)])
        ans=ans2
    return ans

def mat_vec(M,v):
    return [sum(int(M[i][j])*v[j] for j in range(len(v))) for i in range(len(M))]

def frac_mod1(q):
    # canonical in [0,1)
    n=q.numerator % q.denominator
    return Fraction(n,q.denominator)

def lcm(a,b): return abs(a//math.gcd(a,b)*b) if a and b else 0

def root_order(phases):
    m=1
    for q in phases:
        q=frac_mod1(q); m=lcm(m,q.denominator)
    return m

def dm_to_lists(dm):
    return [[int(x) for x in row] for row in dm.to_list()]

def compatible_prime(U, rank, P,Q,R,p):
    OP,OQ,OR=label_options(P,p),label_options(Q,p),label_options(R,p)
    out=[]
    for ia,a in enumerate(OP):
      for ib,b in enumerate(OQ):
       for ic,c in enumerate(OR):
        rhs=[(-a[j]+b[j]+c[j])%p for j in range(8)]
        ur=mat_vec(U,rhs)
        if all(z%p==0 for z in ur[rank:]): out.append((a,b,c))
    return out

summary={}
all_components=0
all_lifts=0
for T,tx in D.items():
    P=Pstd[T[0]]
    type_orders=Counter(); type_comp_orders=Counter(); type_snf=Counter(); type_tors=Counter(); type_dims=Counter()
    max_order=1; max_example=None; total_lifts=0; total_components=0
    for cl in tx['classes']:
        Q,R=cl['Q'],cl['R']; A=build(P,Q,R); n=len(A[0])
        dm=DomainMatrix([[ZZ(z) for z in row] for row in A],(8,n),ZZ)
        S,U,V=smith_normal_decomp(dm)
        Sm=dm_to_lists(S); Um=dm_to_lists(U); Vm=dm_to_lists(V)
        diag=[]
        for i in range(min(8,n)):
            z=abs(Sm[i][i])
            if z: diag.append(z)
        rank=len(diag); tors=math.prod(diag) if diag else 1
        # Sanity dimension
        dim=n-rank-2
        assert dim==cl['dim'], (T,cl['id'],dim,cl['dim'])
        # verify U A V = S using sympy cheaply for each class
        assert U*dm*V == S
        primes=[]
        for ch in T:
            for p in {'A':[2],'B':[2,3],'C':[3,5]}[ch]:
                if p not in primes: primes.append(p)
        cps={p:compatible_prime(Um,rank,P,Q,R,p) for p in primes}
        lifts=math.prod(len(cps[p]) for p in primes)
        assert lifts==cl['normalized_label_lifts'], (T,cl['id'],lifts,cl['normalized_label_lifts'])
        comp_count=lifts*tors
        total_lifts+=lifts; total_components+=comp_count
        type_snf[tuple(diag)]+=1; type_tors[tors]+=1; type_dims[dim]+=1
        # enumerate full lifts and SNF torsion branches; choose free SNF coordinates = 0,
        # then normalize the two matrix rows by their first entries.
        class_orders=Counter()
        for combo in itertools.product(*(cps[p] for p in primes)):
            la=[Fraction(0) for _ in range(8)]; lb=[Fraction(0) for _ in range(8)]; lc=[Fraction(0) for _ in range(8)]
            for p,(a,b,c) in zip(primes,combo):
                for j in range(8):
                    la[j]+=Fraction(a[j],p); lb[j]+=Fraction(b[j],p); lc[j]+=Fraction(c[j],p)
            rhs=[-la[j]+lb[j]+lc[j] for j in range(8)]
            bp=[sum(Fraction(Um[i][j])*rhs[j] for j in range(8)) for i in range(8)]
            assert all(frac_mod1(bp[i])==0 for i in range(rank,8))
            ranges=[range(diag[i]) for i in range(rank)]
            for ks in itertools.product(*ranges):
                y=[Fraction(0) for _ in range(n)]
                for i in range(rank): y[i]=(bp[i]+ks[i])/Sm[i][i]
                x=[sum(Fraction(Vm[i][j])*y[j] for j in range(n)) for i in range(n)]
                assert all(frac_mod1(sum(Fraction(A[i][j])*x[j] for j in range(n))-rhs[i])==0 for i in range(8))
                # rows u and v from P,Q phases + labels
                bpidx={j:i for i,b in enumerate(P) for j in b}
                bqidx={j:i for i,b in enumerate(Q) for j in b}
                u=[x[bpidx[j]]+la[j] for j in range(8)]
                v=[x[len(P)+bqidx[j]]+lb[j] for j in range(8)]
                u=[q-u[0] for q in u]; v=[q-v[0] for q in v]
                ro=root_order(u+v)
                class_orders[ro]+=1; type_comp_orders[ro]+=1
                if ro>max_order:
                    max_order=ro; max_example={'class_id':cl['id'],'dim':dim,'smith':diag,'torsion_branch':list(ks),'root_order':ro,'Q':Q,'R':R}
        assert sum(class_orders.values())==comp_count
        cl['smith']=diag; cl['torsion_components_per_lift']=tors; cl['affine_components']=comp_count
        cl['seed_orders_zero_free_snf']={str(k):v for k,v in sorted(class_orders.items())}
        type_orders[min(class_orders) if class_orders else 1]+=1
    summary[T]={
      'feasible_classes':tx['feasible_classes'],
      'dimension_distribution':{str(k):v for k,v in sorted(type_dims.items())},
      'snf_distribution':{str(k):v for k,v in sorted(type_snf.items(),key=lambda kv:str(kv[0]))},
      'torsion_distribution':{str(k):v for k,v in sorted(type_tors.items())},
      'normalized_label_lifts':total_lifts,
      'affine_components_before_support_automorphisms':total_components,
      'zero_free_snf_seed_order_distribution_over_components':{str(k):v for k,v in sorted(type_comp_orders.items())},
      'max_zero_free_snf_seed_order':max_order,
      'max_example':max_example,
    }
    all_lifts+=total_lifts; all_components+=total_components
    print(T,'lifts',total_lifts,'components',total_components,'maxseed',max_order,'orders',dict(sorted(type_comp_orders.items())))

out={'partition_counts':D0['partition_counts'],'triples':D,'summary':summary,'global':{'feasible_support_classes':sum(x['feasible_classes'] for x in D.values()),'normalized_label_lifts':all_lifts,'affine_components_before_support_automorphisms':all_components}}
json.dump(out,open(OUT,'w'),indent=2)
print('GLOBAL',out['global'])
