from pathlib import Path
from fractions import Fraction as F
import json
from canonical_circles import canonical,verify
R=Path(__file__).resolve().parents[1];W=R/'stage_uniform';seen={};certificates=[]
for stage,path in [('21_22',R),('11_19',R/'stage11_19')]:
 doc=json.loads((path/'outputs/strata.json').read_text());M=doc['modulus']
 for s in doc['strata']:
  if s['base']!=1:continue
  H=[[[F(e[0],M),e[1]] for e in row] for row in s['rows']];key,cert=canonical(H);verify(H,key,cert)
  if key not in seen:seen[key]=dict(stage=stage,id=s['id'],transform=cert)
print('Known distinct circles',len(seen),flush=True)
doc=json.loads((W/'outputs/strata_uniform.json').read_text());aliases=0
for s in doc['strata']:
 if s['base']!=1 or 'alias'in s:continue
 H=[[[F(*e[0]),e[1]] for e in row] for row in s['rows']];key,cert=canonical(H);verify(H,key,cert)
 if key in seen:
  s['alias']=seen[key];s['transform']=cert;s['alias_kind']='affine_circle';aliases+=1
 else:seen[key]=dict(stage='uniform',id=s['id'],transform=cert)
 if s['id']%500==0:print('stratum',s['id'],'circle aliases',aliases,'distinct',len(seen),flush=True)
(W/'outputs/strata_equivalent.json').write_text(json.dumps(doc,separators=(',',':'))+'\n')
summary=dict(all_positive_certificates_verified=True,circle_aliases=aliases,distinct_circles_including_prior_stages=len(seen))
(W/'checks/circle_equivalence_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)
