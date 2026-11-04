from torus3 import *
R=Path(__file__).resolve().parents[1]
cs={c['class_id']:c['rows'] for c in json.load(open(R/'inputs/components.json'))}
out=[]
for c in json.load(open(R/'outputs/A_events.json')):
    seen={}
    for ev in c['loci']:
        if len(ev['eq'])==1:
            for t in split(ev['eq']):seen[key(t)]=t
    for k,t in sorted(seen.items()):out.append({'class_id':c['class_id'],'torus':t,'H':cs[c['class_id']]})
mod=write_scan(R/'inputs/planes.txt',out)
(R/'outputs/planes.json').write_text(json.dumps(enc(out),indent=2))
(R/'outputs/planes_modulus.txt').write_text(str(mod))
print('planes',len(out),'mod',mod)
