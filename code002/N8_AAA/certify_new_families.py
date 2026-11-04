from local_regular import *
from math import factorial

def matching(q):
    c=Counter(tuple(e) for e in q);ans=1
    while c:
        e=min(c);n=c.pop(e);o=((e[0]+F(1,2))%1,)+e[1:]
        if c.pop(o,0)!=n:return 0
        ans*=factorial(n)
    return ans

targets=read('ten_representatives.json');info=[];local=[]
for i,t in enumerate(targets):
    H=[[tuple(F(*x) for x in e) for e in row] for row in t['square']]
    inv=[];pro=[]
    for K in (H,list(zip(*H))):
        inv.append(dict(sorted(Counter(matching(quotient(a,b)) for a,b in combinations(K,2)).items())))
        pro.append(dict(Counter(profiles(quotient(a,b)) for a,b in combinations(K,2))))
    d=len(H[0][0])-1;assert all('' not in x for x in pro)
    if i in [7,8,9]:
        a=local_locus(H);a['representative']=i;local.append(a)
        assert len(a['families'])==1 and a['families'][0]['dimension']==d
    info.append({'representative':i,'dimension':d,'matching_invariants':inv,'profiles':pro})
    print(i,d,inv,pro,flush=True)
assert len({str(x['matching_invariants']) for x in info})==10
(R/'outputs/ten_family_invariants.json').write_text(json.dumps(info,indent=2))
(R/'outputs/new_generic_local_loci.json').write_text(json.dumps(local,indent=2))
