"""Direct symbolic validity and redundancy checks for the seven F4 branches."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from itertools import combinations
import json
from pilot import seed, MOD
from independent import PAT
BASE=Path(__file__).parent
D=json.loads((BASE/'F4_results.json').read_text())
G=[d for d in D['loci'] if d['parameter'] is None]
assert len(G)==7 and all(not d['free_components'] for d in G)
rows=[[(Q(t['constant']),t['slope']) for t in d['row']] for d in G]
# Independently transcribed sixth-root exponents of the seven printed rows.
printed_constants=[[0,4,3,4,2,0,1],[0,4,1,0,4,2,3],
 [0,5,2,2,4,2,5],[0,3,2,4,0,4,1],[0,0,0,3,2,4,3],
 [0,3,1,5,5,3,2],[0,5,3,3,3,1,0]]
for i,e in enumerate(printed_constants):
    slopes=[0]*7 if i==4 else [0,0,0,1,0,0,1]
    assert rows[i]==[(Q(c,6),n) for c,n in zip(e,slopes)]
assert len({tuple(r[:3]) for r in rows})==7 # distinct for every parameter
for d in D['loci']:
    if d['parameter'] is None:continue
    q=Q(d['parameter']);r=tuple(Q(t['constant']) for t in d['row'])
    assert r in [tuple((c+n*q)%1 for c,n in s) for s in rows]
S=[[(Q(c,MOD),n) for c,n in r] for r in seed('F4')]
phi6=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
cert=[]
for i,row in enumerate(rows):
    pairs=[]
    for k,old in enumerate(S):
        dif=[((c-d)%1,n-m) for (c,n),(d,m) in zip(row,old)]
        coeff=defaultdict(lambda:[0,0])
        for c,n in dif:
            assert (6*c).denominator==1
            a,b=phi6[int(6*c)%6];coeff[n][0]+=a;coeff[n][1]+=b
        assert all(v==[0,0] for v in coeff.values())
        matches=[p for p in PAT['C'] if all(dif[a][1]==dif[b][1] and (dif[a][0]-dif[b][0]-Q(e,MOD))%1==0 for a,b,e in p)]
        assert matches
        pairs.append(dict(seed_row=k+1,cycle_edges=matches[0]))
    cert.append(dict(branch=i+1,quotients=pairs))
(BASE/'candidate_checks.json').write_text(json.dumps(dict(seven_distinct_for_all_a=True,special_loci_redundant=42,certificates=cert),indent=2)+'\n')
print('Seven distinct symbolic rows checked directly against all three seed rows.')
print('All 42 fixed-parameter loci are already contained in the seven branches.')
