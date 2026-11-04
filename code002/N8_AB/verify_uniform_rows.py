from pathlib import Path
import json
from verify_rows import verify
R=Path(__file__).resolve().parents[1];W=R/'stage_uniform'
cs=json.loads((W/'inputs/components.json').read_text());total=0;summary=[]
for s in cs:
 i=s['component_index'];d=json.loads((W/f'outputs/full_rows_{i:02d}.json').read_text())[0]
 count=verify([d],[dict(s,id=i)],30);total+=count;summary.append(dict(component=i,class_id=s['class_id'],dim=s['dim'],row_witnesses=count));print(i,count,flush=True)
(W/'checks/row_witness_verification.json').write_text(json.dumps(dict(all_witnesses_valid=True,total=total,components=summary),indent=2)+'\n');print('Verified',total,'row families.',flush=True)
