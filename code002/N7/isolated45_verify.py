"""Input crosswalk, regularity certificates and Fourier identification."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json
from seeds import seeds,MOD
from cleanup import direct_type,quotient
BASE=Path(__file__).resolve().parent
ss={s['label']:s for s in seeds()}
ref=json.loads((BASE/'input_reference.json').read_text())['representatives']
assert len(ref)==45 and {r['cleanup_label'] for r in ref}==set(ss)
for r in ref:
    expected=[[(x-row[0])*(MOD//r['modulus'])%MOD for x in row] for row in r['exponents']]
    assert ss[r['cleanup_label']]['rows']==expected
out=[]
for label,seed in ss.items():
    d=json.loads((BASE/'results'/(label+'.json')).read_text());assert d['seed']==seed
    oldtypes=[]
    for i,j in combinations(range(3),2):
        t,p=direct_type(quotient(seed['rows'][i],seed['rows'][j]));assert t and p;oldtypes.append(t)
    assert ''.join(sorted(oldtypes))==seed['kind']
    cert=[]
    for locus in d['loci']:
        assert locus['dimension']==0 and locus['components']==[0]*7
        r=locus['constants'];quotients=[]
        for old in seed['rows']:
            t,p=direct_type(quotient(r,old));assert t and p
            quotients.append(dict(type=t,cycle_supports=p))
        cert.append(dict(row=r,quotients=quotients))
    if label in ('A1','A2'):
        initial=[0,1,2] if label=='A1' else [0,1,3]
        missing=sorted(set(range(7))-set(initial))
        expected=[[j*k*(MOD//7)%MOD for j in range(7)] for k in missing]
        assert sorted(expected)==sorted(p['constants'] for p in d['loci'])
        assert d['compatibility']['extension_counts']=={'4':4,'5':6,'6':4,'7':1}
        allrows=seed['rows']+expected
        row_order=initial+missing
        assert allrows==[[j*k*(MOD//7)%MOD for j in range(7)] for k in row_order]
        fourier=dict(missing_characters=missing,square_row_order=row_order)
    else:
        assert not d['loci'];fourier=None
    out.append(dict(label=label,kind=seed['kind'],additional_rows=len(cert),certificates=cert,fourier_identification=fourier))
ind=json.loads((BASE/'independent_results.json').read_text())
assert len(ind)==45 and all(r['matches_primary'] for r in ind)
assert sum(r['solution_loci'] for r in ind)==8
summary=dict(seeds=45,blocked_at_four_rows=43,fourier_survivors=['A1','A2'],
 no_continuous_loci=True,
 independent_totals={key:sum(r[key] for r in ind) for key in ['a_candidates','a_survivors','bc_pair_systems','consistent_bc_pairs','bc_pair_loci','bc_third_systems','consistent_bc_thirds']},results=out)
(BASE/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print('45 inputs cross-checked. Both complete enumerations agree.')
print('43 seeds have no regular fourth row. A1,A2 each have exactly four remaining Fourier rows.')
print('Each survivor has complete compatibility graph K4 and a unique Fourier square completion.')
