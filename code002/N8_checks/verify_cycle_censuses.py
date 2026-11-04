"""Exact, independently organized checks for the two added cycle censuses."""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
from functools import lru_cache
from fractions import Fraction
import json, re, sys
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'n8_aaa_complete/code'))
from exact_rows import profiles

def read_exponents(name):
    t=(ROOT/'report'/name).read_text()
    t=t.split(r'\begin{pmatrix}')[1].split(r'\end{pmatrix}')[0]
    a=[[int(x) for x in row.strip().split('&')] for row in t.split(r'\\') if row.strip()]
    assert len(a)==8 and all(len(row)==8 for row in a)
    return a

def admissible(q,vals,parts):
    @lru_cache(None)
    def rec(v,ps):
        if not ps:return not v
        for p in set(ps):
            if q%p:continue
            # The cycle containing the first occurrence can be rotated to start there.
            cycle=Counter((v[0]+k*(q//p))%q for k in range(p))
            cc=Counter(v)
            if any(cc[x]<n for x,n in cycle.items()):continue
            cc.subtract(cycle)
            rest=tuple(sorted(x for x,n in cc.items() for _ in range(n)))
            pp=list(ps);pp.remove(p)
            if rec(rest,tuple(pp)):return True
        return False
    return rec(tuple(sorted(vals)),tuple(parts))

def census(E,q):
    pairs={};checks=0
    x=S.symbols('x');phi=S.Poly(S.cyclotomic_poly(q,x),x)
    for i,j in combinations(range(8),2):
        d=[(a-b)%q for a,b in zip(E[i],E[j])]
        assert S.rem(S.Poly(sum(x**v for v in d),x),phi).is_zero
        p=''.join(t for t,ps in [('A',[2]*4),('B',[2,3,3]),('C',[3,5])] if admissible(q,d,ps))
        assert p==profiles([(Fraction(v,q),) for v in d])
        pairs[i+1,j+1]=p;checks+=1
    exclusive=Counter();inclusive=Counter()
    for triple in combinations(range(1,9),3):
        labs=[pairs[ij] for ij in combinations(triple,2)]
        types=tuple(sorted(set(''.join(sorted(p)) for p in product(*labs))))
        exclusive[types]+=1;inclusive.update(types)
    assert sum(exclusive.values())==56
    return {'pairs':{f'{i},{j}':p for (i,j),p in pairs.items()},
            'pair_counts':dict(Counter(pairs.values())),
            'disjoint_triples':{' / '.join(k):v for k,v in sorted(exclusive.items())},
            'overlapping_triples':dict(sorted(inclusive.items())),
            'exact_orthogonality_checks':checks}

def main():
    out={}
    for name,file,q in [('I12','isolated_matrix.tex',12),('K','ccc_matrix.tex',30)]:
        E=read_exponents(file);assert all(0<=v<q for r in E for v in r)
        c=census(E,q);ct=census(list(zip(*E)),q)
        assert c['pair_counts']==ct['pair_counts']
        assert c['disjoint_triples']==ct['disjoint_triples']
        out[name]={'q':q,'exponents':E,'rows':c,'columns':ct}
    assert out['I12']['rows']['disjoint_triples']=={'AAA':16,'AAA / AAB':16,'AAB':16,'AAB / ABB':8}
    assert out['I12']['rows']['overlapping_triples']=={'AAA':32,'AAB':40,'ABB':8}
    assert out['K']['rows']['disjoint_triples']=={'AAC':48,'CCC':8}
    assert out['I12']['rows']['pair_counts']=={'A':20,'B':4,'AB':4}
    assert out['K']['rows']['pair_counts']=={'A':16,'C':12}
    Kpairs=out['K']['rows']['pairs']
    for group in ([1,2,4,7],[3,5,6,8]):
        assert all(Kpairs[f'{i},{j}']=='C' for i,j in combinations(group,2))
    raw=(ROOT/'n8_aaa_complete/inputs/aalto/BH-8-12.txt').read_bytes().splitlines()[739]
    assert [v for r in out['I12']['exponents'] for v in r]==[x-48 for x in raw[6:]]
    out['all_checks_passed']=True
    (ROOT/'checks/cycle_censuses.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: both matrices, 112 exact row/column orthogonality and profile checks, 224 triples; catalogue line 740 matches I12.')

if __name__=='__main__':main()
