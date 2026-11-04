"""Exact ordinary equivalence of affine circles of three prescribed rows.

Normalize row/column phases, parameter orientation and common slope gcd. A
canonical nonzero slope entry fixes the parameter translation, with every
finite choice of its zero retained. This gives a complete canonical key under
these transformations and explicit certificates of every positive match.
"""
from fractions import Fraction as F
from itertools import permutations
from math import gcd
from functools import reduce

def canonical(H):
 g=reduce(gcd,(abs(e[1]) for row in H for e in row),0);assert g
 shape=None;configs=[]
 for rp in permutations(range(3)):
  for pivot in range(8):
   slope=[tuple((H[r][j][1]-H[r][pivot][1]-H[rp[0]][j][1]+H[rp[0]][pivot][1])//g for r in rp) for j in range(8)]
   for sign in [1,-1]:
    v=[tuple(sign*x for x in col) for col in slope];key=tuple(sorted(v))
    if shape is None or key<shape:shape=key;configs=[]
    if key==shape:configs.append((rp,pivot,sign,v))
 anchor=next(v for v in shape if any(v));row=next(i for i,x in enumerate(anchor) if x)
 best=None;certificate=None
 for rp,pivot,sign,slopes in configs:
  const=[tuple((H[r][j][0]-H[r][pivot][0]-H[rp[0]][j][0]+H[rp[0]][pivot][0])%1 for r in rp) for j in range(8)]
  for j in range(8):
   if slopes[j]!=anchor:continue
   k=anchor[row]
   for branch in range(abs(k)):
    shift=(F(branch)-const[j][row])/k
    cols=[tuple((slopes[x][r],(const[x][r]+slopes[x][r]*shift)%1) for r in range(3)) for x in range(8)]
    cp=sorted(range(8),key=lambda x:cols[x]);key=tuple(cols[x] for x in cp)
    if best is None or key<best:
     best=key;certificate=dict(rows=list(rp),pivot_column=pivot,columns=cp,slope_gcd=g,parameter_sign=sign,parameter_shift=[shift.numerator,shift.denominator])
 return best,certificate

def verify(H,key,w):
 p=w['rows'];cp=w['columns'];pivot=w['pivot_column'];g=w['slope_gcd'];sign=w['parameter_sign'];shift=F(*w['parameter_shift'])
 assert sorted(p)==list(range(3)) and sorted(cp)==list(range(8)) and sign in [-1,1] and g>0
 cols=[]
 for j in cp:
  col=[]
  for r in p:
   c=H[r][j][0]-H[r][pivot][0]-H[p[0]][j][0]+H[p[0]][pivot][0]
   k=H[r][j][1]-H[r][pivot][1]-H[p[0]][j][1]+H[p[0]][pivot][1]
   assert k%g==0;k=k//g*sign
   col.append((k,(c+k*shift)%1))
  cols.append(tuple(col))
 assert tuple(cols)==key
 return True
