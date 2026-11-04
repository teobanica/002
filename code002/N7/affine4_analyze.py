"""Projection loci, direct validity, F4-triple witnesses, and pilot regression."""
from fractions import Fraction as Q
from itertools import permutations,product,combinations
from collections import Counter,defaultdict
from math import gcd
from functools import reduce
from pathlib import Path
import json,time
from families import family
from census import PAT,deserialize,canonical
from independent import kernel
ROOT=Path(__file__).resolve().parent
POWERS=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]

def four_rows(seed,comp):
    c,B,p=comp;d=seed['parameters'];k=len(p);cs=[];Bs=[]
    for oldc,oldB in zip(seed['constants'],seed['basis']):
        cs.append(tuple((Q(v,6)+sum(a*c[h] for h,a in enumerate(row)))%1 for v,row in zip(oldc,oldB)))
        Bs.append(tuple(tuple(sum(row[h]*B[h][j] for h in range(d)) for j in range(k)) for row in oldB))
    cs.append(c[d:]);Bs.append(B[d:]);return cs,Bs

def valid_patterns(c,B):
    return [p for p in PAT if all(B[i]==B[j] and (c[i]-c[j]-Q(e,6))%1==0 for i,j,e in p)]

def f4_match(cs,Bs):
    target=family('F4');tc=tuple(Q(v,6) for row in target['constants'][1:] for v in row)
    tb=tuple(row[0] for rr in target['basis'][1:] for row in rr)
    for base,u,v in permutations(range(4),3):
        cc=[tuple((cs[r][j]-cs[base][j])%1 for j in range(7)) for r in (u,v)]
        bb=[tuple(tuple(x-y for x,y in zip(Bs[r][j],Bs[base][j])) for j in range(7)) for r in (u,v)]
        for pat in valid_patterns(cc[0],bb[0]):
            tri=(pat[0][1],pat[0][0],pat[1][0])
            if pat[0][2]==4:tri=(tri[0],tri[2],tri[1])
            P=(pat[2][1],pat[2][0]);S=(pat[3][1],pat[3][0])
            for shift,swap,revp,revs in product(range(3),range(2),range(2),range(2)):
                TT=tri[shift:]+tri[:shift];PP=P[::-1] if revp else P;SS=S[::-1] if revs else S
                if swap:PP,SS=SS,PP
                cols=TT+PP+SS;a=cols[0]
                nc=tuple((cc[r][j]-cc[r][a])%1 for r in range(2) for j in cols)
                nb=tuple(tuple(x-y for x,y in zip(bb[r][j],bb[r][a])) for r in range(2) for j in cols)
                ac=nc[5];ab=nb[5]
                if all((nc[i]-tc[i]-tb[i]*ac)%1==0 and nb[i]==tuple(tb[i]*x for x in ab) for i in range(14)):
                    return dict(rows=[base+1,u+1,v+1],columns=[j+1 for j in cols],parameter_constant=str(ac),parameter_basis=ab)
    return None

def base_locus(comp,d):
    c,B,p=comp;k=len(p)
    assert k<=1
    vec=[B[i][0] for i in range(d)] if k else [0]*d
    g=reduce(gcd,(abs(x) for x in vec),0)
    if g==0:
        eq=tuple(tuple(int(i==j) for j in range(d)) for i in range(d));rhs=c[:d]
        return (eq,rhs),dict(projection_dimension=0,extra_row_parameters=k,fibre_components=1,equations=eq,rhs=[str(v) for v in rhs])
    primitive=[x//g for x in vec]
    if next(x for x in primitive if x)<0:primitive=[-x for x in primitive]
    L=kernel((tuple(primitive),),d);eq=tuple(zip(*L))
    rhs=tuple(sum(a*b for a,b in zip(row,c[:d]))%1 for row in eq)
    return (eq,rhs),dict(projection_dimension=1,extra_row_parameters=0,fibre_components=g,equations=eq,rhs=[str(v) for v in rhs])

def main():
    summary=[]
    for name in ['F4','F5','F6','F2','F3','F1']:
        data=json.loads((ROOT/'results'/(name+'.json')).read_text());seed=family(name);d=seed['parameters']
        records=[];groups={}
        for idx,r in enumerate(data['maximal_components'],1):
            comp=deserialize(r);cs,Bs=four_rows(seed,comp);witness=[]
            for i,j in combinations(range(4),2):
                c=tuple((x-y)%1 for x,y in zip(cs[i],cs[j]));B=tuple(tuple(x-y for x,y in zip(a,b)) for a,b in zip(Bs[i],Bs[j]))
                ps=valid_patterns(c,B);assert ps,(name,idx,i,j)
                witness.append(dict(rows=[i+1,j+1],cycle_edges=ps[0]))
            key,info=base_locus(comp,d);match=f4_match(cs,Bs)
            if key not in groups:groups[key]=dict(label=len(groups)+1,**info,branches=[])
            groups[key]['branches'].append(idx)
            records.append(dict(id=idx,base_locus=groups[key]['label'],**info,f4_witness=match,cycle_certificates=witness))
        nfixed=sum(r['dimension']==0 for r in data['maximal_components'])
        stats=dict(name=name,seed_parameters=d,components=len(records),circle_components=len(records)-nfixed,points=nfixed,
                   base_loci=len(groups),base_curves=sum(g['projection_dimension']==1 for g in groups.values()),
                   base_points=sum(g['projection_dimension']==0 for g in groups.values()),
                   extra_phase_branches=sum(r['extra_row_parameters']>0 for r in records),
                   wholly_f4=sum(r['f4_witness'] is not None for r in records),
                   remaining=sum(r['f4_witness'] is None for r in records))
        out=dict(summary=stats,base_loci=list(groups.values()),branches=records)
        (ROOT/'results'/(name+'_analysis.json')).write_text(json.dumps(out,indent=2)+'\n')
        print(stats,flush=True);summary.append(stats)
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
