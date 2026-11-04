import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
for name,source,mod in [('generic','generic_AAA.json',30),('lines','class22_line_scan.json',60)]:
    D=json.loads((R/'outputs'/source).read_text())
    if name=='generic':D=[c for c in D if c['class_id']==22]
    else:D=[c for c in D if not c['k5']]
    with (R/f'inputs/pairs_{name}.txt').open('w') as f:
        f.write(str(len(D))+'\n')
        for c in D:
            f.write(f"{c['class_id']} {c['branch']} {c['dim']} {len(c['rows'])}\n")
            for r in c['rows']:f.write(' '.join(str(v) for e in r for v in e)+'\n')
