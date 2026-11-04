"""Check the mathematical data transcribed into the readable PDF."""
import json
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
R=Path(__file__).resolve().parents[1]
def parse(x):
    sign=x.startswith('-');x=x[1:] if sign else x
    assert x in ('1','y')
    return (int(sign),int(x=='y'))
rows=['1 1 1 1 1 1 1 1','1 -1 y -y -y y -y y','1 -1 -1 y 1 -y -1 1',
      '1 -1 1 -1 y -1 1 -y','1 -1 -y 1 -1 1 y -1','1 1 -1 -y y y -1 -y',
      '1 1 y y -1 -y -y -1','1 1 -y -1 -y -1 y y']
H=[[parse(x) for x in row.split()] for row in rows]
for r,s in combinations(H,2):
    c=Counter(((a[0]-b[0])%2,a[1]-b[1]) for a,b in zip(r,s))
    assert all(c[(1-h,k)]==n for (h,k),n in c.items())
table={1:[3,4,9,10],2:[3,4,5,9,11],3:[1,2,4,5,7,8,10,11],4:[1,3,7,9],5:[2,3,8,9],7:[3,4,9,10],8:[3,5,9,11],9:[1,2,4,5,7,8,10,11],10:[1,3,7,8,9],11:[2,3,8,9]}
expected={(F(i,12),F(j,12)) for i,js in table.items() for j in js}
D=json.loads((R/'outputs/final_residue.json').read_text())
for d in D:
    actual={(F(*p[0]),(F(*p[1])+F(1,2))%1) for p in d['points']}
    assert actual==expected and len(actual)==50
print('PASS: displayed square and complete 50-point PDF table.')
