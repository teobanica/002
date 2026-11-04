"""Split every C-pair square locus exactly and verify every pair profile."""
from pathlib import Path
import sys,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from phase_algebra import solutions,decode,substitute,enc
from exact_rows import profiles,quotient
ds=json.loads((R/'outputs/fixed_pairs.json').read_text())
ls=json.loads((R/'outputs/c_square_loci.json').read_text())
out=[];unique={}
for ci,(d,l) in enumerate(zip(ds,ls)):
 fs=[d['families'][i] for i in d['keep']];count=0
 for li,x in enumerate(l['squares']):
  selected=[fs[i] for i in x['families']];nf=sum(f['free'] for f in selected)
  for branch,values in enumerate(solutions(x['eq'],nf,0)):
   extra=len(values[0])-1 if nf else 0
   H=[tuple((F(0),)+(F(0),)*extra for j in range(8)),tuple((F(c,30),)+(F(0),)*extra for c in d['base_row'])];pos=0
   for fam in selected:
    f=fam['free'];H.append(substitute(decode(fam['row']),0,values[pos:pos+f],extra));pos+=f
   ps=[profiles(quotient(a,b)) for a,b in combinations(H,2)]
   assert all(ps) and len(H)==8
   assert extra==0,('continuous survivor',ci,li,branch,extra)
   assert all((e[0]*30).denominator==1 for row in H for e in row)
   matrix=[[int(e[0]*30)%30 for e in row] for row in H]
   cp=[profiles(quotient(a,b)) for a,b in combinations(tuple(zip(*H)),2)]
   assert all(cp)
   rowkey=tuple(sorted(tuple(row) for row in matrix));ix=unique.setdefault(rowkey,len(unique))
   out.append({'pair_index':ci,'locus_index':li,'torsion_branch':branch,'unique_rowset_index':ix,'matrix_mod30':matrix,'row_profiles':dict(Counter(ps)),'column_profiles':dict(Counter(cp))})
   count+=1
 print('pair',ci,'loci',len(l['squares']),'verified square branches',count,flush=True)
(R/'outputs/c_squares.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
# Add the known K as the final canonicalization input.
K=[[0,0,0,0,0,0,0,0],[0,0,6,10,12,18,20,24],[0,6,13,28,27,21,12,15],[0,10,28,20,16,4,10,22],[0,12,27,16,1,15,24,9],[0,18,21,4,15,19,6,3],[0,20,12,10,24,6,0,18],[0,24,15,22,9,3,18,7]]
matrices=[list(map(list,k)) for k,v in sorted(unique.items(),key=lambda kv:kv[1])]+[K]
with (R/'inputs/canonical_squares.txt').open('w') as f:
 f.write(f'{len(matrices)} 30\n')
 for H in matrices:f.write(' '.join(str(x) for row in H for x in row)+'\n')
print('Unique row sets:',len(unique),'and known K added for equivalence checking.')
