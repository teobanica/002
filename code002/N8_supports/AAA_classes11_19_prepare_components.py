import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
cs=[c for c in json.load(open(R/'outputs/components.json')) if c['class_id'] in (11,19)]
assert len(cs)==2 and all(c['dim']==3 and c['branch']==0 for c in cs)
(R/'inputs/components.json').write_text(json.dumps(cs,indent=2))
with (R/'inputs/components.txt').open('w') as f:
    f.write('2\n')
    for c in cs:
        f.write(f"{c['class_id']} 0 3\n")
        for row in c['rows']:f.write(' '.join(str(x) for e in row for x in e)+'\n')
print('Selected the unique three-parameter components of classes 11 and 19.')
