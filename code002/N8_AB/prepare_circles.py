from pathlib import Path
import json
R=Path(__file__).resolve().parents[1];out=[]
for cid in [21,22]:
 p=R.parent/f'proofs/AAA_class{cid}/outputs/class{cid}_square_templates.json'
 for t in json.loads(p.read_text()):
  M=t['constant_modulus'];rows=[]
  for row in t['matrix'][:3]:
   assert all(e[0]*30%M==0 for e in row)
   rows.append([[e[0]*30//M]+e[1:] for e in row])
  out.append(dict(id=len(out),class_id=cid,branch=t['branch'],line=t['line'],dim=1,rows=rows))
assert len(out)==34
(R/'inputs/circle_components.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
with (R/'inputs/circle_components.txt').open('w') as f:
 f.write(str(len(out))+'\n')
 for s in out:
  f.write(f"{s['id']} {s['branch']} 1\n")
  for row in s['rows']:f.write(' '.join(str(x) for e in row for x in e)+'\n')
print('Prepared all 34 successful circles from the preserved exact templates.')
