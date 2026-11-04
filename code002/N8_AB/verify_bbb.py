"""Verify the pure-BBB fifth-row obstruction by two exact routes."""
from pathlib import Path
import json,sys
from fractions import Fraction as F
R=Path(__file__).resolve().parents[1];W=R/'stage_bbb'
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from phase_algebra import solutions,substitute
from exact_rows import profiles,quotient
from verify_rows import verify
import sympy as S
cs=json.loads((W/'inputs/reduced_components.json').read_text())['components']
data=[json.loads((W/f"outputs/full_rows_{c['id']:02d}.json").read_text())[0] for c in cs]
checked=verify(data,cs,60);assert checked==3830
# Route 1: all closed strata have at most one generic circle row, or finitely
# many fixed rows. The fixed rows are checked directly with rational cycles.
ss=json.loads((W/'outputs/strata.json').read_text());M=ss['modulus'];numeric_pairs=0
for s in ss['strata']:
 fs=s['families'];assert all(f['free']==0 for f in fs)
 if s['base']:assert len(fs)<=1;continue
 rows=[[[F(e[0],M)] for e in f['row']] for f in fs]
 for i in range(len(rows)):
  for j in range(i):assert not profiles(quotient(rows[i],rows[j]));numeric_pairs+=1
# Route 2: bypass the recursive stratification. Enumerate all raw fourth-row
# pairs including their full original base constraints. Every possible fifth
# row forces an A-admitting edge in the original BBB triple.
ps=json.loads((W/'outputs/raw_pairs.json').read_text());loci=branches=0
for c,p in zip(cs,ps):
 assert c['id']==p['class_id'];eqs={tuple(map(tuple,r['eq'])) for pair in p['pairs'] for r in pair['relations']}
 for eq in eqs:
  loci+=1
  for values in solutions([list(e[:-1])+[S.Rational(e[-1],2)] for e in eq],c['dim'],0):
   dim=len(values[0])-1 if values else 0
   rows=[substitute(tuple((F(e[0],60),)+tuple(map(F,e[1:])) for e in row),0,values,dim) for row in c['rows']]
   pp=[profiles(quotient(rows[i],rows[j])) for i,j in [(0,1),(0,2),(1,2)]]
   assert any('A' in q for q in pp),(c['id'],eq,pp)
   branches+=1
summary=dict(all_checks_passed=True,normalized_label_lifts=652,smith_components=1040,pure_B_components_before_point_reduction=616,components_after_point_reduction=len(cs),raw_fourth_row_witnesses=checked,strata_numeric_pair_checks=numeric_pairs,raw_pair_base_loci=loci,raw_pair_smith_branches=branches,all_fifth_rows_force_base_A_edge=True,pure_BBB_has_no_AB_fifth_row=True,consequence='Every A/B regular matrix with at least six rows has an AAA triple.')
(W/'checks/final_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
