"""Produce focused input files, without hard-coding the discovered locus."""
import json,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
mode=sys.argv[1]
if mode=='components':
    D=[c for c in json.loads((R/'outputs/components.json').read_text()) if c['class_id']==22]
    assert len(D)==2
    with (R/'inputs/class22_components.txt').open('w') as f:
        f.write('2\n')
        for c in D:
            f.write(f"22 {c['branch']} 2\n")
            for r in c['rows']:f.write(' '.join(str(v) for e in r for v in e)+'\n')
elif mode=='strata':
    D=json.loads((R/'outputs/class22_A_residue.json').read_text())
    assert D[0]['good_lines']==D[1]['good_lines']
    with (R/'inputs/good_lines.txt').open('w') as f:
        for a,b,c in D[0]['good_lines']:f.write(f'{a} {b} {c[0]} {c[1]}\n')
    with (R/'inputs/extra_points.txt').open('w') as f:
        for p in D[0]['off_good_points']:f.write(' '.join(str(v) for q in p for v in q)+'\n')
else:raise ValueError(mode)
