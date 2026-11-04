"""Exact Laurent-polynomial checks for P7(q), R7(q), and their conjugacy.

Entries are omega^E q^D, omega=exp(2*pi*i/6). P uses w=omega^2.
All computations are standard-library integer/rational arithmetic.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
P=[[3,0,2,0,2,0,2],[0,3,2,0,0,2,2],[2,2,5,0,2,2,0],
   [0,0,0,3,2,2,2],[2,0,2,2,5,2,0],[0,2,2,2,2,5,0],[2,2,0,2,0,0,3]]
PD=[[int(i<2 and j<2)-int(4<=i<=5 and 4<=j<=5) for j in range(7)] for i in range(7)]
R=[[0,0,0,0,0,0,0],[0,0,0,2,3,3,4],[0,1,4,4,0,3,2],
   [0,2,3,0,4,1,4],[0,3,5,3,3,0,1],[0,4,2,2,0,5,3],[0,4,2,5,2,3,0]]
RD=[[int((2<=i<=4 and 4<=j<=5) or (5<=i<=6 and 3<=j<=6)) for j in range(7)] for i in range(7)]
POWERS=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))


def check_cycles(E,D):
    certificates=[]
    for i,j in combinations(range(7),2):
        phases=[(E[i][k]-E[j][k])%6 for k in range(7)]
        degrees=[D[i][k]-D[j][k] for k in range(7)]
        coeff=defaultdict(lambda:[0,0])
        for e,d in zip(phases,degrees):
            for h in range(2):coeff[d][h]+=POWERS[e][h]
        assert all(value==[0,0] for value in coeff.values())
        found=None
        for T in combinations(range(7),3):
            rest=[k for k in range(7) if k not in T]
            for A in combinations(rest,2):
                B=tuple(k for k in rest if k not in A)
                blocks=(T,A,B)
                if all(len({degrees[k] for k in block})==1 and
                       {(phases[k]-phases[block[0]])%6 for k in block}
                       ==set(range(0,6,6//len(block))) for block in blocks):
                    found=blocks;break
            if found:break
        assert found is not None
        certificates.append({'rows':[i+1,j+1],'cycle_supports':[[k+1 for k in b] for b in found],
                             'phases':phases,'degrees':degrees})
    return certificates


def dephase(E,rows,cols,modulus=None):
    out=[[E[r][c]-E[rows[0]][c]-E[r][cols[0]]+E[rows[0]][cols[0]]
          for c in cols] for r in rows]
    return out if modulus is None else [[x%modulus for x in row] for row in out]


def equivalence(E,target):
    target=dephase(target,tuple(range(7)),tuple(range(7)),6)
    target_cols=[tuple(target[i][j] for i in range(1,7)) for j in range(7)]
    for rows in permutations(range(7)):
        B=[[(E[r][j]-E[rows[0]][j])%6 for j in range(7)] for r in rows[1:]]
        for anchor in range(7):
            lookup={tuple((B[i][j]-B[i][anchor])%6 for i in range(6)):j for j in range(7)}
            if set(lookup)==set(target_cols):
                cols=[lookup[c] for c in target_cols]
                assert dephase(E,rows,cols,6)==target
                return {'rows':[i+1 for i in rows],'columns':[j+1 for j in cols]}
    return None


def rank(A):
    A=[[F(x) for x in r] for r in A];r=0
    for j in range(len(A[0])):
        k=next((k for k in range(r,len(A)) if A[k][j]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r]
        p=A[r][j];A[r]=[x/p for x in A[r]]
        for k in range(r+1,len(A)):
            if A[k][j]:
                p=A[k][j];A[k]=[x-p*y for x,y in zip(A[k],A[r])]
        r+=1
        if r==len(A):break
    return r


def tangent(E,dephased):
    variables=[(i,j) for i in range(7) for j in range(7)
               if not dephased or (i and j)]
    re=(2,1,-1,-2,-1,1);im=(0,1,1,0,-1,-1)
    equations=[]
    for i,k in combinations(range(7),2):
        for values in (re,im):
            equations.append([((a==i)-(a==k))*values[(E[i][j]-E[k][j])%6]
                              for a,j in variables])
    r=rank(equations)
    return {'variables':len(variables),'rank':r,'dimension':len(variables)-r}


def main():
    result={'P':{'phase_exponents':P,'parameter_exponents':PD,'cycles':check_cycles(P,PD)},
            'R':{'phase_exponents':R,'parameter_exponents':RD,'cycles':check_cycles(R,RD)}}
    rows=(0,1,6,2,3,5,4);cols=(3,6,2,4,1,0,5)
    assert dephase([[-x for x in r] for r in P],rows,cols,6)==R
    assert dephase([[-x for x in r] for r in PD],rows,cols)==RD
    result['conjugacy']={'identity':'R(q) = dephase(conjugate(P(q))[rows, columns])',
                         'rows':[x+1 for x in rows],'columns':[x+1 for x in cols]}
    seeds=[]
    for line in (ROOT.parents[1]/'sources/BH-7-6.txt').read_text().splitlines():
        assert line.startswith('BH1776') and len(line[6:])==49
        seeds.append([list(map(int,line[6+7*i:6+7*(i+1)])) for i in range(7)])
    assert R==seeds[1]
    mapping=equivalence(P,seeds[0]);assert mapping is not None
    assert equivalence(P,R) is None
    result['seed1_equivalence']=mapping
    result['seed_classes_distinct']=True
    for name,E in (('P',P),('R',R)):
        result[name]['undephased_tangent']=tangent(E,False)
        result[name]['dephased_tangent']=tangent(E,True)
        assert result[name]['undephased_tangent']['dimension']==16
        assert result[name]['dephased_tangent']['dimension']==3
    (ROOT/'square_family_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    text=('P and R: all 21 row pairs verified coefficient by coefficient in q.\n'
          'Both families have symbolic 3+2+2 cycle decompositions for every row pair.\n'
          'R(q) = dephase(conjugate(P(q))[1,2,7,3,4,6,5; 4,7,3,5,2,1,6]) exactly.\n'
          'P(1) matches catalogue seed 1; R(1) equals seed 2; those seeds are inequivalent.\n'
          'Both seeds: undephased tangent dimension 16; dephased defect 3.\n')
    (ROOT/'run.log').write_text(text);print(text,end='')


if __name__=='__main__':
    main()
