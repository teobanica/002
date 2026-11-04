"""Compare two complete enumerations and certify the generic reduction."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json
R=Path(__file__).resolve().parents[1]
a=json.loads((R/'outputs/c_pair_rows.json').read_text())
b=json.loads((R/'outputs/independent_c_rows.json').read_text())
assert a['attempted']==b['attempted']==3913**2
assert a['profile_systems']==b['profile_systems']
assert b['non_mu30_roots']==0
def key(phase,free,row):return phase,free,tuple(tuple(e) for e in row)
expected=Counter();restricted_phases=set()
for f in a['families']:
    eq=f['base_equations']
    if not eq:
        expected[key(-1,f['free'],f['row'])]+=f['multiplicity'];continue
    k,c=eq[0];assert k>0
    for r in range(k):
        t=(F(-c,30)+r)/k%1
        if any((kk*t+F(cc,30))%1 for kk,cc in eq):continue
        assert (t*30).denominator==1
        phase=int(t*30);restricted_phases.add(phase)
        row=[[(e[0]+e[1]*phase)%30]+e[2:] for e in f['row']]
        expected[key(phase,f['free'],row)]+=f['multiplicity']
observed=Counter({key(f['phase'],f['free'],f['row']):f['multiplicity'] for f in b['records']})
assert expected==observed, (len(expected),len(observed),list((expected-observed).items())[:3],list((observed-expected).items())[:3])
generic=[f for f in a['families'] if not f['base_equations']]
circles=[f for f in generic if f['free']]
assert len(generic)==93 and len(circles)==3
assert all(f['profiles']==['CC'] and f['free']==1 for f in circles)
slopes=set();containment=[]
for f in circles:
    row=f['row'];assert row[:3]==[[0,0,0],[20,0,0],[10,0,0]]
    assert all(e[1:]==[0,1] for e in row[3:]) and row[3][0]==0
    s=row[4][0]//6;assert [e[0] for e in row[3:]]==[(6*s*j)%30 for j in range(5)]
    slopes.add(s)
assert slopes=={2,3,4}
for i,f in enumerate(generic):
    if f['free']:continue
    for j,c in enumerate(circles):
        v=f['row'][3]
        row=[[(e[0]+e[2]*v[0])%30,e[1]+e[2]*v[1]] for e in c['row']]
        if row==f['row']:
            containment.append({'generic_fixed_index':i,'circle':j,'parameter':v});break
    else:raise AssertionError(('uncovered generic row',f))
assert len(containment)==90
# The pair phase changes by 6b-10a in Z/30 under independent cyclic
# permutations of the pentagon/triangle blocks and rephasing the second row.
orbit=sorted({(6*b-10*a)%30 for a in range(3) for b in range(5)})
assert orbit==list(range(0,30,2))
summary={'passed':True,'support_pairs':a['attempted'],'profile_systems':a['profile_systems'],
 'independent_split_records':len(expected),'independent_multiplicities_agree':True,
 'exceptional_base_phases_mod30':sorted(restricted_phases),'generic_records':93,
 'generic_circles':3,'generic_fixed_rows_contained':90,'generic_slopes':sorted(slopes),
 'even_phase_orbit':orbit,'containment':containment}
(R/'checks/reduction_verified.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS: independent support enumeration, all base torsion, all multiplicities, 93 generic rows reduced to three circles; phases reduce to +1/-1.')
