"""Exact symbolic coverage certificates, without parameter sampling."""
from itertools import permutations, product
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,time
from independent_enumeration import decode,ROOT
from classify_components import (canonical,contains,dimension,oriented_cycles,regular_family,serialize)
from families import templates

def transformed_raw(f,rows,cols):
    c,B,_=f; d=len(B[0])
    RC=[(0,)*7,c[:7],c[7:]]
    RB=[((0,)*d,)*7,B[:7],B[7:]]
    base,first,second=rows; anchor=cols[0]
    cc=tuple((RC[r][j]-RC[base][j]-RC[r][anchor]+RC[base][anchor])%24
             for r in (first,second) for j in cols)
    BB=tuple(tuple(RB[r][j][h]-RB[base][j][h]-RB[r][anchor][h]+RB[base][anchor][h]
                   for h in range(d)) for r in (first,second) for j in cols)
    return cc,BB

def transforms(f):
    """All equivalent images normalized to u=(1,w,w^2,a,-a,b,-b)."""
    c,B,_=f;d=len(B[0])
    RC=[(0,)*7,c[:7],c[7:]]; RB=[((0,)*d,)*7,B[:7],B[7:]]
    for rows in permutations(range(3)):
        base,first,second=rows
        uc=tuple((RC[first][j]-RC[base][j])%24 for j in range(7))
        UB=tuple(tuple(RB[first][j][h]-RB[base][j][h] for h in range(d)) for j in range(7))
        for (T,P,Q),sgn in oriented_cycles(uc,UB):
            T=T if sgn==1 else (T[0],T[2],T[1])
            for shift,swap,fp,fq in product(range(3),range(2),range(2),range(2)):
                TT=T[shift:]+T[:shift]
                PP=P[::-1] if fp else P;QQ=Q[::-1] if fq else Q
                cols=TT+(QQ+PP if swap else PP+QQ)
                cc,BB=transformed_raw(f,rows,cols)
                yield canonical(cc,BB),{'rows':rows,'columns':cols},(cc,BB)

def as24(f):
    c,B,P=f
    v=[x*24 for x in c]
    assert all(x.denominator==1 for x in v),'Constant outside 24th roots.'
    return canonical([int(x) for x in v],B)

def fits_raw(target,raw):
    bc,B,P=target;sc,S=raw
    return all((sc[i]-bc[i]-sum(B[i][h]*sc[p] for h,p in enumerate(P)))%24==0
       and all(S[i][k]==sum(B[i][h]*S[p][k] for h,p in enumerate(P)) for k in range(len(S[0])))
       for i in range(14))

def run():
    data=json.loads((ROOT/'enumeration.json').read_text())
    source=[as24(decode(f)) for f in data['components']]
    targets=templates()
    assert len(targets)==17 and all(regular_family(f) for _,f in targets)
    normalized=[]
    for name,f in targets:
        g,trace,raw=min(transforms(f),key=lambda r:(r[0],r[1]['rows'],r[1]['columns']))
        normalized.append((name,g,trace))
    certificates=[];tally=Counter()
    for i,f in enumerate(source):
        found=None
        for _,trace,raw in transforms(f):
            for name,g,_ in normalized:
                if dimension(g)>=dimension(f) and fits_raw(g,raw):
                    found={'source':i,'target':name,'transformation':trace}
                    break
            if found:break
        assert found is not None,('UNCOVERED',i)
        certificates.append(found);tally[found['target']]+=1
    # Verify each certificate from the source and the specified finite transformations.
    for cert in certificates:
        f=source[cert['source']]
        trace=cert['transformation']
        raw=transformed_raw(f,trace['rows'],trace['columns'])
        target=next(g for name,g,_ in normalized if name==cert['target'])
        assert fits_raw(target,raw)
    result={'targets':[{'name':name,'family':serialize(f)} for name,f in targets],
            'normalized_targets':[{'name':name,'family':serialize(f),'transformation':trace}
                                  for name,f,trace in normalized],
            'coverage_counts':dict(tally),'certificates':certificates}
    (ROOT/'coverage.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print('Covered all',len(source),'normalized solution tori:',dict(tally),flush=True)
    # Exhaustive symbolic test: no representative family is contained in any other.
    # Source cycle enumeration is exhaustive for equivalences into a normalized target.
    for i,(name,f) in enumerate(targets):
        for _,trace,raw in transforms(f):
            for j,(other,g,_) in enumerate(normalized):
                if i!=j and dimension(g)>=dimension(f):
                    assert not fits_raw(g,raw),(name,'contained in',other)
    print('All six positive-dimensional families are necessary; eleven point classes are distinct and outside them.',flush=True)
    return result

if __name__=='__main__':
    start=time.perf_counter();run();print('Elapsed',round(time.perf_counter()-start,2),'seconds.')

