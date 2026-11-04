"""Independent graph agreement, equivalence witnesses, and representative checks."""
from pathlib import Path
import sys,json
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from math import gcd,lcm
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from exact_rows import profiles,quotient
import sympy as S
read=lambda n:json.loads((R/n).read_text())
ss=read('outputs/c_squares.json');ind=read('outputs/independent_fixed_squares.json');cs=read('outputs/c_canonical.json')
assert read('checks/reduction_verified.json')['passed']
for ci in range(2):
 expected={tuple(sorted(tuple(2*x for x in row) for row in d['matrix_mod30'][2:])) for d in ss if d['pair_index']==ci}
 observed={tuple(tuple(row) for row in sq) for sq in ind[ci]['square_rowsets_mod60']}
 assert expected==observed and len(observed)==(16 if ci==0 else 0)
 assert ind[ci]['rank_zero_loci']==ind[ci]['free_line_triples']==0
 # Every clique of three, four, or five added rows lies in one square.
 sqs=[list(map(tuple,s)) for s in ind[ci]['square_rowsets_mod60']]
 for k in [3,4,5,6]:
  contained=Counter(tuple(z) for s in sqs for z in combinations(s,k))
  assert len(contained)==ind[ci]['clique_counts'][k]
  assert all(n==1 for n in contained.values())
 ds=read('outputs/fixed_pairs.json')[ci];fs=[ds['families'][i] for i in ds['keep']]
 linepoints=set()
 for i,j,c in ind[ci]['free_lines']:
  for t in range(60):
   row=lambda f,u:tuple((2*e[0]+e[1]*u)%60 for e in f['row'])
   linepoints.add(tuple(sorted((row(fs[i],t),row(fs[j],t+2*c)))))
 rows=list(map(tuple,ind[ci]['rows_mod60']))
 edges={tuple((rows[i],rows[j])) for i,j in ind[ci]['edge_indices']}
 assert linepoints<=edges
 discrete=edges-linepoints
 assert len(discrete)==(772 if ci==0 else 112)
 ind[ci]['four_row_circle_count']=len(ind[ci]['free_lines'])
 ind[ci]['four_row_isolated_points']=len(discrete)
 (R/f'outputs/four_row_isolated_{ci}.json').write_text(json.dumps(sorted(discrete),separators=(',',':'))+'\n')
# Recover exact matrices consumed by the canonicalizer.
data=list(map(int,(R/'inputs/canonical_squares.txt').read_text().split()));n,mod=data[:2];assert mod==30
mats=[[[data[2+64*i+8*r+c] for c in range(8)] for r in range(8)] for i in range(n)]
groups={}
for rec,H in zip(cs,mats):
 a,b=rec['pivot'];rp=rec['rows'];cp=rec['columns']
 assert sorted(rp)==sorted(cp)==list(range(8))
 out=[[(H[i][j]-H[a][j]-H[i][b]+H[a][b])%30 for j in cp] for i in rp]
 assert out==rec['canonical']
 groups.setdefault(tuple(map(tuple,out)),[]).append(rec['id'])
assert sorted(map(len,groups.values()))==[2,2,13]
K=cs[-1]['canonical'];L=cs[12]['canonical'];LT=cs[13]['canonical']
assert K==list(map(list,zip(*K))) and LT==list(map(list,zip(*L)))
assert len({tuple(map(tuple,h)) for h in [K,L,LT]})==3
z=S.symbols('z');phi=S.Poly(S.cyclotomic_poly(30,z),z)
def cycle_certificate(row,parts):
 def rec(unused,ps):
  if not ps:return []
  a=min(unused)
  for p in sorted(set(ps)):
   chosen=[]
   for k in range(p):
    want=(row[a]+30*k//p)%30
    js=[j for j in sorted(unused) if row[j]==want and j not in chosen]
    if not js:break
    chosen.append(js[0])
   if len(chosen)==p:
    pp=list(ps);pp.remove(p);tail=rec(unused-set(chosen),pp)
    if tail is not None:return [chosen]+tail
  return None
 return rec(set(range(8)),parts)
reps=[]
for name,H in [('K',K),('L30',L),('L30T',LT)]:
 result={'name':name,'matrix_mod30':H,'minimal_level':30}
 assert lcm(*(30//gcd(30,x) for row in H for x in row))==30
 for direction,hh in [('rows',H),('columns',list(map(list,zip(*H))))]:
  ps={};cert=[]
  for i,j in combinations(range(8),2):
   q=[(a-b)%30 for a,b in zip(hh[i],hh[j])]
   assert S.rem(S.Poly(sum(z**x for x in q),z),phi).is_zero
   p=profiles([(F(x,30),) for x in q]);assert p in ['A','C'];ps[i,j]=p
   cycles=cycle_certificate(q,[2]*4 if p=='A' else [3,5]);assert cycles is not None
   cert.append({'pair':[i+1,j+1],'profile':p,'cycles':[[k+1 for k in block] for block in cycles]})
  triples=Counter(''.join(sorted(ps[min(i,j),max(i,j)] for i,j in combinations(t,2))) for t in combinations(range(8),3))
  todo=set(range(8));cliques=[]
  while todo:
   i=min(todo);block={i}|{j for j in todo if j!=i and ps[min(i,j),max(i,j)]=='C'}
   assert all(ps[min(a,b),max(a,b)]=='C' for a,b in combinations(block,2));todo-=block;cliques.append([j+1 for j in sorted(block)])
  result[direction]={'profiles':dict(Counter(ps.values())),'triple_census':dict(triples),'C_cliques':cliques,'cycle_certificates':cert}
 reps.append(result)
summary={'all_checks_passed':True,'global_C_sector_complete':True,'global_AAA_complete':False,
 'square_presentations':16,'classes':3,'class_members':list(groups.values()),
 'pair_extension_counts':[{k:v for k,v in d.items() if k not in ['free_lines','square_rowsets_mod60','rows_mod60','edge_indices']} for d in ind],
 'every_five_to_seven_row_C_partial_has_unique_square_completion':True,
 'canonicalization_witnesses_verified':len(cs),'representatives':reps}
(R/'outputs/c_representatives.json').write_text(json.dumps(reps,indent=2)+'\n')
(R/'checks/final_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS: independent square sets agree exactly, all 17 equivalence witnesses verified, three ordinary classes K/L30/L30T, exact row/column orthogonality and all cycle certificates.')
