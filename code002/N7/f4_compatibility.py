"""All regular extensions of F4(a), and exact identification of its completions."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,permutations
from collections import defaultdict
import json
BASE=Path(__file__).parent
DATA=json.loads((BASE/'F4_results.json').read_text())
G=[d for d in DATA['loci'] if d['parameter'] is None]
X=[[int(6*Q(e['constant'])) for e in d['row']] for d in G]
XD=[[e['slope'] for e in d['row']] for d in G]
SE=[[0]*7,[0,2,4,2,5,0,3],[0,2,4,0,2,3,5]]
SD=[[0]*7,[0,0,0,1,1,1,1],[0,0,0,1,1,1,1]]
P=[[3,0,2,0,2,0,2],[0,3,2,0,0,2,2],[2,2,5,0,2,2,0],
 [0,0,0,3,2,2,2],[2,0,2,2,5,2,0],[0,2,2,2,2,5,0],[2,2,0,2,0,0,3]]
PD=[[int(i<2 and j<2)-int(4<=i<=5 and 4<=j<=5) for j in range(7)] for i in range(7)]
R=[[0,0,0,0,0,0,0],[0,0,0,2,3,3,4],[0,1,4,4,0,3,2],
 [0,2,3,0,4,1,4],[0,3,5,3,3,0,1],[0,4,2,2,0,5,3],[0,4,2,5,2,3,0]]
RD=[[int((2<=i<=4 and 4<=j<=5) or (5<=i<=6 and 3<=j<=6)) for j in range(7)] for i in range(7)]
POWERS=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]

def gram(e,d,f,g):
    out=defaultdict(lambda:[0,0])
    for a,b,c,h in zip(e,d,f,g):
        x,y=POWERS[(a-c)%6];out[b-h][0]+=x;out[b-h][1]+=y
    return {k:v for k,v in sorted(out.items()) if v!=[0,0]}

def cycle(e,d,f,g):
    phases=[(a-b)%6 for a,b in zip(e,f)];deg=[a-b for a,b in zip(d,g)]
    for tri in combinations(range(7),3):
        rest=[j for j in range(7) if j not in tri]
        for k in rest[1:]:
            pair=(rest[0],k);last=tuple(j for j in rest if j not in pair)
            blocks=(tri,pair,last)
            if all(len({deg[j] for j in b})==1 and
                   {(phases[j]-phases[b[0]])%6 for j in b}==set(range(0,6,6//len(b))) for b in blocks):
                return [[j+1 for j in b] for b in blocks]
    return None

def dephase(E,rows,cols,mod=None):
    out=[[E[r][c]-E[rows[0]][c]-E[r][cols[0]]+E[rows[0]][cols[0]] for c in cols] for r in rows]
    return out if mod is None else [[v%mod for v in r] for r in out]

def match(E,D,T,TD):
    # Find D((source(omega^shift*a^sign))[rows,columns]) = target(a).
    target=dephase(T,range(7),range(7),6);td=dephase(TD,range(7),range(7))
    targetcols=[tuple((target[i][j],td[i][j]) for i in range(1,7)) for j in range(7)]
    targetdegrees={tuple(td[i][j] for i in range(1,7)) for j in range(7)}
    for sign in [1,-1]:
        DD=[[sign*v for v in row] for row in D]
        for rows in permutations(range(7)):
            diff=[[DD[r][j]-DD[rows[0]][j] for j in range(7)] for r in rows[1:]]
            for anchor in range(7):
                degrees=[tuple(diff[i][j]-diff[i][anchor] for i in range(6)) for j in range(7)]
                if set(degrees)!=targetdegrees:continue
                for shift in range(6):
                    EE=[[(c+shift*n)%6 for c,n in zip(e,d)] for e,d in zip(E,D)]
                    ec=dephase(EE,rows,range(7),6)
                    ec=[[(v-row[anchor])%6 for v in row] for row in ec]
                    lookup={tuple((ec[i][j],degrees[j][i-1]) for i in range(1,7)):j for j in range(7)}
                    if set(lookup)!=set(targetcols):continue
                    cols=[lookup[k] for k in targetcols]
                    assert dephase(EE,rows,cols,6)==target
                    assert dephase(DD,rows,cols)==td
                    return dict(shift=shift,sign=sign,rows=[j+1 for j in rows],columns=[j+1 for j in cols])
    return None

def main():
    edges=set();pairs=[]
    for i,j in combinations(range(7),2):
        coeff=gram(X[i],XD[i],X[j],XD[j])
        assert not coeff or list(coeff)==[0] # all Gram entries independent of a
        cyc=cycle(X[i],XD[i],X[j],XD[j])
        if not coeff:
            assert cyc is not None;edges.add((i+1,j+1))
        else:assert cyc is None
        pairs.append(dict(pair=[i+1,j+1],gram_coefficients=coeff,cycle=cyc))
    cliques={}
    for k in range(1,8):
        cliques[k]=[s for s in combinations(range(1,8),k) if all(p in edges for p in combinations(s,2))]
    maximal=[s for k in cliques for s in cliques[k] if not any(set(s)<set(t) for kk in cliques for t in cliques[kk])]
    assert maximal==[(1,3,5,6),(2,4,5,7)]
    assert {k:len(v) for k,v in cliques.items()}=={1:7,2:12,3:8,4:2,5:0,6:0,7:0}
    # Check the cross Gram matrix printed in the proof, in order U=(1,3,6), V=(2,4,7).
    cross=[[[1,-3],[1,4],[5,-1]],[[1,4],[1,-3],[5,-1]],[[5,-1],[5,-1],[-3,2]]]
    for r,i in enumerate([1,3,6]):
        for c,j in enumerate([2,4,7]):
            assert gram(X[i-1],XD[i-1],X[j-1],XD[j-1])=={0:cross[r][c]}
    out=dict(pairs=pairs,edges=sorted(edges),cliques=cliques,maximal_cliques=maximal,completions=[])
    for name,selection in zip(['H_plus','H_minus'],maximal):
        E=SE+[X[j-1] for j in selection];D=SD+[XD[j-1] for j in selection]
        for i,j in combinations(range(7),2):
            assert not gram(E[i],D[i],E[j],D[j]);assert cycle(E[i],D[i],E[j],D[j]) is not None
        mapping=None
        for family,e,d in [('P',P,PD),('R',R,RD)]:
            mapping=match(e,d,E,D)
            if mapping:
                mapping['family']=family;break
        assert mapping is not None,name
        out['completions'].append(dict(name=name,selected_rows=selection,phase_exponents=E,parameter_exponents=D,equivalence=mapping))
        print(name,mapping,flush=True)
    (BASE/'compatibility_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Edges:',sorted(edges));print('Maximal cliques:',maximal)
    print('Number of fixed-seed extensions at 4,5,6,7 rows: 7,12,8,2, for every a.')

if __name__=='__main__':main()
