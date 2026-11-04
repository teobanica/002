import sys
from torus3 import *
R=Path(__file__).resolve().parents[1]
stage=sys.argv[1]
items=json.load(open(R/f'outputs/{stage}.json'));scans=json.load(open(R/f'outputs/{stage}_scan.json'))
mod=int((R/f'outputs/{stage}_modulus.txt').read_text())
out=[]
for ix,(a,b) in enumerate(zip(items,scans)):
    if b['k5']:continue
    t=decode(a['torus']);rows=[lift(row,t,mod) for row in b['rows']]
    out.append((a['class_id'],ix,t,rows))
M=lcm(30,*(q[0].denominator for _,_,_,rs in out for row in rs for q in row),*(q.denominator for _,_,t,_ in out for q in t['c']))
with (R/f'inputs/{stage}_pairs.txt').open('w') as f:
    f.write(str(len(out))+'\n')
    for cid,ix,t,rows in out:
        eq=equations(t,M);f.write(f'{cid} {ix} {len(eq)} {len(rows)}\n')
        for e in eq:f.write(' '.join(map(str,e))+'\n')
        for row in rows:f.write(' '.join(str(x) for e in row for x in [int(e[0]*M)]+e[1:])+'\n')
(R/f'outputs/{stage}_pairs_modulus.txt').write_text(str(M))
print(stage,len(out),'pair scans, modulus',M)
