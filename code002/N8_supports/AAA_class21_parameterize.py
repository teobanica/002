"""Reconstruct dephased AAA components by exact Smith decomposition."""
import json,itertools
from pathlib import Path
import sympy as s
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
ROOT=Path(__file__).resolve().parents[1]
D=json.loads((ROOT/'inputs/canonical.json').read_text())['triples']['AAA']
out=[]
for cl in D['classes']:
    P=[[0,1],[2,3],[4,5],[6,7]]; Q=cl['Q']; R=cl['R']
    parts=[P,Q,R]; labels=[[int(j==b[1]) for j in range(8) for b in S if j in b] for S in parts]
    A=s.zeros(8,10)
    for k,S in enumerate(parts):
        for i,b in enumerate(S):
            if k<2 and i==0: continue
            ix=([0,3,6][k]+i-(k<2))
            for j in b:A[j,ix]=[1,-1,-1][k]
    dm=DomainMatrix.from_Matrix(A).convert_to(ZZ)
    sd,ud,vd=smith_normal_decomp(dm)
    S,U,V=[z.to_Matrix() for z in (sd,ud,vd)]
    assert U*A*V==S
    rank=A.rank(); dim=10-rank
    assert dim==cl['dim']
    rhs=s.Matrix([s.Rational(-labels[0][j]+labels[1][j]+labels[2][j],2) for j in range(8)])
    b=U*rhs
    assert all(z.q==1 for z in b[rank:])
    for branch,ks in enumerate(itertools.product(*(range(abs(int(S[i,i]))) for i in range(rank)))):
        y=s.zeros(10,1)
        for i in range(rank):y[i]=(b[i]+ks[i])/S[i,i]
        c=V*y; F=V[:,rank:]
        rows=[[[0]*(dim+1) for j in range(8)]]
        for k in range(2):
            row=[]
            for j in range(8):
                bi=next(i for i,bb in enumerate(parts[k]) if j in bb)
                ix=3*k+bi-1
                cc=s.Rational(labels[k][j],2)+(c[ix] if bi else 0)
                assert (30*cc).q==1
                row.append([int(30*cc)%30]+([int(z) for z in F[ix,:]] if bi else [0]*dim))
            rows.append(row)
        out.append(dict(class_id=cl['id'],branch=branch,dim=dim,rows=rows,smith=[int(S[i,i]) for i in range(rank)],Q=Q,R=R))
assert len(out)==38
(ROOT/'outputs/components.json').write_text(json.dumps(out,indent=2))
with (ROOT/'inputs/components.txt').open('w') as f:
    f.write(str(len(out))+'\n')
    for c in out:
        f.write(f"{c['class_id']} {c['branch']} {c['dim']}\n")
        for r in c['rows']:f.write(' '.join(str(v) for e in r for v in e)+'\n')
print('Reconstructed',len(out),'components from',len(D['classes']),'support classes.')
for c in out:
    if c['class_id']==21:print(json.dumps(c))
