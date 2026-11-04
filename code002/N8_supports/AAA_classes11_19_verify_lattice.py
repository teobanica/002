import json
from pathlib import Path
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form as hnf
R=Path(__file__).resolve().parents[1]
for x in json.load(open(R/'outputs/lattice_tests.json')):
    d=x['dim'];m=x['mod'];base=[0]*d+[m];M=S.Matrix(x['input']+[base]);H=hnf(M[:,::-1].T).T
    B=H[::-1,::-1];pure=[r[-1] for r in B.tolist() if not any(r[:-1])];assert len(pure)==1
    assert x['empty']==(abs(pure[0])!=m),(x,B)
    if not x['empty']:
        out=S.Matrix(x['basis']+[base]);assert hnf(out.T)==hnf(M.T),(x,B)
        assert out==B,(x,B)
print('PASS: 1200 integer lattice systems independently checked against SymPy Hermite normal form.')
