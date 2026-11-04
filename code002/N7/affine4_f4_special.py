"""Exact F4 containment at EVERY point of each residual one-parameter branch."""
from fractions import Fraction as Q
from itertools import permutations,product
from pathlib import Path
import json,time
from families import family
from census import PAT,deserialize
from analyze import four_rows
ROOT=Path(__file__).resolve().parent

def restrict(allowed,c,n):
    c%=6
    if n==0:return allowed if c==0 else frozenset()
    if allowed is not None:return frozenset(t for t in allowed if (Q(c,6)+n*t)%1==0)
    if n<0:c,n=-c,-n
    return frozenset(((Q(-c,6)+j)/n)%1 for j in range(n))

def special(cs,Bs):
    cs=[[int(6*v)%6 for v in row] for row in cs]
    Bs=[[row[0] for row in rr] for rr in Bs]
    target=family('F4');tc=sum(target['constants'][1:],[]);tb=[v[0] for rr in target['basis'][1:] for v in rr]
    found={}
    for base,u,v in permutations(range(4),3):
        cc=[[(cs[r][j]-cs[base][j])%6 for j in range(7)] for r in (u,v)]
        bb=[[Bs[r][j]-Bs[base][j] for j in range(7)] for r in (u,v)]
        for pat in PAT:
            domain=None
            for i,j,e in pat:
                domain=restrict(domain,cc[0][i]-cc[0][j]-e,bb[0][i]-bb[0][j])
                if domain==frozenset():break
            if domain==frozenset() or (domain is not None and domain.issubset(found)):continue
            tri=(pat[0][1],pat[0][0],pat[1][0])
            if pat[0][2]==4:tri=(tri[0],tri[2],tri[1])
            P=(pat[2][1],pat[2][0]);S=(pat[3][1],pat[3][0])
            for shift,swap,rp,rs in product(range(3),range(2),range(2),range(2)):
                TT=tri[shift:]+tri[:shift];PP=P[::-1] if rp else P;SS=S[::-1] if rs else S
                if swap:PP,SS=SS,PP
                cols=TT+PP+SS;anchor=cols[0]
                nc=[(cc[r][j]-cc[r][anchor])%6 for r in range(2) for j in cols]
                nb=[bb[r][j]-bb[r][anchor] for r in range(2) for j in cols]
                ac,ab=nc[5],nb[5];current=domain
                # First impose v's ternary block, which rejects most column orders.
                for i in [8,9]+list(range(14)):
                    current=restrict(current,nc[i]-tc[i]-tb[i]*ac,nb[i]-tb[i]*ab)
                    if current==frozenset():break
                if current==frozenset():continue
                assert current is not None,'A supposedly residual branch is wholly F4.'
                for t in current:
                    if t not in found:
                        found[t]=dict(t=str(t),rows=[base+1,u+1,v+1],columns=[j+1 for j in cols],f4_parameter=str((Q(ac,6)+ab*t)%1))
    return [found[t] for t in sorted(found)]

def main():
    summary=[]
    for name in ['F4','F5','F6','F2','F3','F1']:
        start=time.time();data=json.loads((ROOT/'results'/(name+'.json')).read_text());an=json.loads((ROOT/'results'/(name+'_analysis.json')).read_text())
        records=[]
        for r,a in zip(data['maximal_components'],an['branches']):
            if a['f4_witness'] is not None:
                records.append(dict(id=a['id'],status='all',special=[]));continue
            if r['dimension']==0:
                records.append(dict(id=a['id'],status='none',special=[]));continue
            cs,Bs=four_rows(family(name),deserialize(r));s=special(cs,Bs)
            records.append(dict(id=a['id'],status='finite' if s else 'none',special=s))
        result=dict(name=name,branches=records,seconds=time.time()-start)
        (ROOT/'results'/(name+'_f4_loci.json')).write_text(json.dumps(result,indent=2)+'\n')
        stat=dict(name=name,all=sum(r['status']=='all' for r in records),finite=sum(r['status']=='finite' for r in records),none=sum(r['status']=='none' for r in records),special_points=sum(len(r['special']) for r in records),seconds=result['seconds'])
        summary.append(stat);print(stat,flush=True)
    (ROOT/'f4_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
