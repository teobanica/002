"""Verify exceptional-frame certificates and discharge covered subtrees."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import json
from verify_rows import matchings
R=Path(__file__).resolve().parents[1];W=R/'stage_uniform'
doc=json.loads((W/'outputs/strata_equivalent.json').read_text());ss=doc['strata'];ws=json.loads((W/'outputs/exceptional_frames.json').read_text());certs={w['stratum']:w for w in ws if w['exceptional_class']>=0}
classes={c['id']:c for c in json.loads((R.parent/'proofs/consolidation/outputs/canonical.json').read_text())['triples']['AAA']['classes']};ms=list(matchings(list(range(8))))
for sid,w in certs.items():
 s=ss[sid];H=s['rows'];part=[]
 for mid,(a,b) in zip(w['matchings'],[(0,1),(0,2),(1,2)]):
  blocks=[0]*8
  for k,(i,j) in enumerate(ms[mid]):
   assert (F(*H[a][i][0])-F(*H[b][i][0])-F(*H[a][j][0])+F(*H[b][j][0]))%1==F(1,2)
   assert all(H[a][i][q]-H[b][i][q]-H[a][j][q]+H[b][j][q]==0 for q in range(1,s['base']+1))
   blocks[i]=blocks[j]=k
  part.append(blocks)
 actual=[0]*64
 for j in range(8):actual[16*part[0][j]+4*part[1][j]+part[2][j]]+=1
 assert sorted(w['axis'])==[0,1,2] and all(sorted(p)==[0,1,2,3] for p in w['blocks'])
 image=[0]*64
 for k,count in enumerate(classes[w['exceptional_class']]['tensor']):
  a=[k//16,(k//4)%4,k%4];v=[w['blocks'][i][a[w['axis'][i]]] for i in range(3)];image[16*v[0]+4*v[1]+v[2]]+=count
 assert actual==image
adj=defaultdict(list)
for a,b in doc['edges']:adj[a].append(b)
seen=set();search=[];covered=[]
def visit(i):
 if i in seen:return
 seen.add(i);s=ss[i]
 if i in certs:covered.append(i);return
 if 'alias'in s:
  a=s['alias']
  if a['stage']=='uniform':visit(a['id'])
  return
 search.append(i)
 for j in adj[i]:visit(j)
for r in doc['roots']:visit(r['root'])
pruned=dict(strata=[ss[i] for i in sorted(search)],source='strata_equivalent.json',reachable_nodes=sorted(seen),exceptional_discharge_nodes=sorted(covered),all_frame_certificates_verified=True)
(W/'outputs/strata_pruned.json').write_text(json.dumps(pruned,separators=(',',':'))+'\n')
summary=dict(all_certificates_verified=True,certificates=len(certs),searched_nodes=len(search),dimensions=dict(Counter(ss[i]['base'] for i in search)),reachable_nodes=len(seen),exceptional_discharge_nodes=len(covered))
(W/'checks/frame_pruning_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
