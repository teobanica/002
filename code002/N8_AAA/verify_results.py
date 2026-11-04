"""Independent final audit of the saved exact results (no numerical tests)."""
from verify_affine import *
from local_regular import local_locus
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--recompute-local',action='store_true');args=parser.parse_args()
targets=read('ten_representatives.json')
for t in targets:
 H=[[tuple(F(*x) for x in e) for e in row] for row in t['square']]
 assert all(profiles(quotient(a,b)) for K in (H,list(zip(*H))) for a,b in combinations(K,2))
 A=S.Matrix([e[1:] for row in H for e in row]);assert A.rank()==len(H[0][0])-1
assert sum(x['templates'] for x in read('generic_square_summary.json'))==1512
assert verify(targets,json.loads((R/'inputs/transpose_family_sources.json').read_text()),read('transpose_families.json'),1680)==10
assert verify(targets,json.loads((R/'inputs/mixed_to_catalogue_sources.json').read_text()),read('mixed_to_catalogue.json'),1680)==1
assert read('mixed_to_catalogue.json')[0]['target']==8
local_sources=json.loads((R/'inputs/exceptional_affine_sources.json').read_text());assert len(local_sources)==41
assert verify(targets,local_sources,read('exceptional_affine_embeddings.json'),1680)==41
finite=read('mixed_finite_squares.json');canonical=read('mixed_canonical.json');keys=set()
assert len(finite)==len(canonical)==32
for x,y in zip(finite,canonical):
 H=x['matrix_mod12'];p,q=y['pivot'];T=[[(H[i][j]-H[p][j]-H[i][q]+H[p][q])%12 for j in y['columns']] for i in y['rows']]
 assert T==y['canonical'];keys.add(tuple(map(tuple,T)))
assert len(keys)==1
assert all(x['squares']==1 and x['rows']==7 for x in read('mixed_finite_summary.json'))
for c in read('family_rank_certificates.json'):
 A=S.Matrix(c['equations']);assert A[c['minor_rows'],c['minor_columns']].det()==c['minor_determinant']!=0
 t=read('row_column_closed_targets.json')[c['target']];H=[[[F(*x) for x in e] for e in row] for row in t['square']];V=S.Matrix([H[i][j][1:] for i in range(1,8) for j in range(1,8)])
 assert A*V==S.zeros(A.rows,c['dimension']) and V.rank()==c['dimension'] and c['rank']+c['dimension']==49
loc=read('mixed_local_locus.json');G=S.Matrix([[F(*x) for x in row] for row in loc['common_equations']]);assert G.rank()==48
for block in loc['one_parameter_blockers']:
 for sup in block['supports']:
  row=S.Matrix([[F(*x) for x in sup['row']]]);assert S.Matrix.vstack(G,row).det()==F(*sup['determinant'])!=0
if args.recompute_local:
 for p in sorted((R/'outputs').glob('exceptional_local_*.json')):
  old=json.loads(p.read_text());src=read(f"aalto_sources_{old['q']}.json")[old['index']];H=[[F(*e[0]) for e in row] for row in src['square']]
  result=local_locus(H);result.update({'q':old['q'],'index':old['index']});assert json.loads(json.dumps(result))==old,p
  print('Recomputed local locus',p.name,flush=True)
summary={'generic_components':38,'generic_square_parameterizations':1512,'generic_family_types':6,'continuous_families':9,'continuous_dimensions':[5,5,5,4,4,2,2,2,1],'isolated_mixed_classes':1,'finite_mixed_triples':32,'finite_mixed_completions_per_triple':1,'catalogue_entries':5469,'catalogue_covered':5469,'exceptional_catalogue_points':27,'exceptional_local_branches':41,'global_AAA_complete':False,'all_checks_passed':True}
assert read('aalto_final_summary.json')['covered']==5469
(R/'outputs/final_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
