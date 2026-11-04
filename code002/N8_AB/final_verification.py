"""Final exact certificate audit for the complete A/B family cover.

Checks the saved finite outputs. Exhaustive enumeration is reproduced by
reproduce.py --full; a positive witness check alone is not an exhaustion.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from collections import Counter
import json,sys,hashlib
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from exact_rows import profiles,quotient
from verify_affine import verify

def read(p):return json.loads(p.read_text())
def fixed(rows,mod=None):
 return [[F(e[0],mod) if mod else F(*e[0]) for e in row] for row in rows]
def image(H,w):
 rp,cp,p=w['rows'],w['columns'],w['pivot_column']
 assert sorted(rp)==list(range(3)) and sorted(cp)==list(range(8)) and 0<=p<8
 return tuple(tuple((H[r][j]-H[r][p]-H[rp[0]][j]+H[rp[0]][p])%1 for r in rp) for j in cp)
# Independent verification of all fixed-triple equivalences in the first stages.
a=read(R/'outputs/strata.json');seeds=read(R/'inputs/circle_components.json');aliases=0
for w in a['point_sources']:
 s=seeds[w['circle']];t=F(*w['point']);H=[[(F(e[0],30)+e[1]*t)%1 for e in row] for row in s['rows']]
 target=a['strata'][w['representative']]
 assert image(H,w['canonical_transform'])==image(fixed(target['rows'],a['modulus']),target['canonical_transform'])
 aliases+=1
b=read(R/'stage11_19/outputs/strata.json');raw=read(R/'outputs/strata_11_19.json')
for w in b['point_equivalences']:
 s=raw['strata'][w['source']];q=w['target'];st=a if q['stage']=='21_22' else b;t=st['strata'][q['id']]
 assert image(fixed(s['rows']),w['transform'])==image(fixed(t['rows'],st['modulus']),q['transform'])
 aliases+=1
bbb=read(R/'stage_bbb/inputs/components.json');br=read(R/'stage_bbb/inputs/reduced_components.json')
for w in br['point_equivalences']:
 s=bbb['components'][w['source']];q=w['target'];t=br['components'][q['id']]
 assert image(fixed(s['rows'],bbb['modulus']),w['transform'])==image(fixed(t['rows'],br['modulus']),q['transform'])
 aliases+=1
assert sum(c['label_lifts'] for c in bbb['counts'])==652
assert sum(c['smith_branches'] for c in bbb['counts'])==1040
assert bbb['uniformly_AAA_skipped']==392 and len(br['components'])==102
print('Verified',aliases,'additional fixed-triple certificates.',flush=True)

@lru_cache(maxsize=250000)
def checked_profile(q):
 p=profiles(q);assert 'A' in p or 'B' in p
 return p

targets=read(R.parent/'n8_aaa_complete/outputs/ten_representatives.json')
summary={'additional_fixed_triple_certificates':aliases,'stages':{}}
for stage,want,searches,loci,branches in [('',396,55,28846,32534),('stage11_19',1557,177,106633,119397),('stage_uniform',10048,400,845666,896184)]:
 W=R/stage;ss=read(W/'outputs/strata.json');ex=read(W/'outputs/square_extraction_summary.json');sq=read(W/'outputs/completion_squares.json');emb=read(W/'outputs/completion_embeddings.json')
 assert len(ss['strata'])==len(ex['counts'])==searches and ex['all_strata_finished']
 assert [s['id'] for s in ss['strata']]==list(range(searches))
 assert sum(c['loci'] for c in ex['counts'])==loci and sum(c['smith_branches'] for c in ex['counts'])==branches
 assert len(sq)==want==ex['unique'];assert len(emb)==want and all(w['target']>=0 for w in emb)
 assert {str(k):v for k,v in Counter(s['dimension'] for s in sq).items()}==ex['dimensions']
 # Every complete per-stratum output is present; totals are not read from logs.
 for s,c in zip(ss['strata'],ex['counts']):
  d=read(W/f"outputs/cliques_stratum_{s['id']:02d}.json");assert len(d)==1 and d[0]['component_index']==s['id']==c['stratum']
  assert len(d[0]['squares'])==c['loci']==d[0]['loci']
 # All Smith branches have an origin, after exact duplicate-leaf removal.
 assert sum(len(s['origins']) for s in sq)==sum(c.get('distinct_leaf_smith_branches',c['smith_branches']) for c in ex['counts'])
 assert verify(targets,sq,emb,ex['modulus'])==want
 for s in sq:
  H=tuple(tuple(tuple(F(*x) for x in e) for e in row) for row in s['square']);dim=s['dimension']
  assert len(H)==8 and all(len(row)==8 and all(len(e)==dim+1 for e in row) for row in H)
  assert all(not any(e) for e in H[0]) and all(not any(row[0]) for row in H)
  pp=[checked_profile(tuple(tuple(e) for e in quotient(H[i],H[j]))) for i,j in combinations(range(8),2)]
  assert pp==s['profiles']
 tc=dict(Counter(w['target'] for w in emb))
 check=dict(all_covered=True,verified=want,total=want,target_counts=tc)
 (W/'checks/completion_embedding_verification.json').write_text(json.dumps(check,indent=2)+'\n')
 summary['stages'][stage or 'classes21_22']=dict(squares=want,searches=searches,lattice_loci=loci,smith_branches=branches,all_embeddings_verified=True,all_symbolic_cycles_verified=True,target_counts=tc)
 print(stage or 'classes21_22',want,'complete formulas verified.',flush=True)
# Previously executed independent checks must all be present and successful.
requirements=[('checks/row_witness_verification.json','all_circle_witnesses_valid'),('stage_uniform/checks/row_witness_verification.json','all_witnesses_valid'),('stage_uniform/checks/alias_verification.json','all_equivalences_verified'),('stage_uniform/checks/frame_pruning_verification.json','all_certificates_verified'),('stage_uniform/checks/smith_verification.json','all_checks_passed'),('stage_bbb/checks/final_verification.json','all_checks_passed')]
for p,k in requirements:assert read(R/p)[k],(p,k)
assert read(R/'stage_uniform/checks/row_witness_verification.json')['total']==548310
assert read(R/'checks/row_witness_verification.json')['classes11_19_row_witnesses']==34006
assert read(R/'stage_uniform/checks/alias_verification.json')['total']==15415
assert read(R/'stage_uniform/checks/frame_pruning_verification.json')['searched_nodes']==400
assert read(R/'stage_bbb/checks/final_verification.json')['pure_BBB_has_no_AB_fifth_row']
assert read(R.parent/'n8_aaa_complete/outputs/final_verification.json')['all_checks_passed']
assert read(R.parent/'n8_C_sector/checks/final_verification.json')['all_checks_passed']
summary.update(all_checks_passed=True,new_square_formulas=12001,exceptional_AAA_complete=True,uniform_AAA_all_special_strata_complete=True,pure_B_fifth_row_obstruction=True,AB_squares_always_contain_AAA=True,regular_order_eight_family_cover_complete=True,continuous_family_types=9,isolated_ordinary_classes=['I12','K','L30','L30_transpose'],scope='Complete family cover; no claim of a disjoint parameter atlas or classification of nonregular matrices.')
(R/'checks/final_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
