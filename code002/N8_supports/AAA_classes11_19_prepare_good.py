from torus3 import *
R=Path(__file__).resolve().parents[1]
p=json.load(open(R/'outputs/planes.json'));s=json.load(open(R/'outputs/planes_scan.json'))
good=[]
def isA(row,m):
    from collections import Counter
    c=Counter(map(tuple,row))
    return all(c[((e[0]+m//2)%m,)+e[1:]]==n for e,n in c.items())
for a,b in zip(p,s):
    if not b['k5']:continue
    t=decode(a['torus']);H=specialize(a['H'],t,30);square=H+[b['rows'][i] for i in b['k5']]
    assert len(square)==8
    for i in range(8):
        for j in range(i):
            q=[[((x[0]-y[0])%30)]+[xx-yy for xx,yy in zip(x[1:],y[1:])] for x,y in zip(square[i],square[j])]
            assert isA(q,30)
    good.append({**a,'square':square,'modulus':30})
(R/'outputs/good_planes.json').write_text(json.dumps(good,indent=2))
for cid in [11,19]:
    ts=[decode(x['torus']) for x in good if x['class_id']==cid]
    with (R/f'inputs/good{cid}.txt').open('w') as f:
        f.write(str(len(ts))+'\n')
        for t in ts:
            assert len(t['A'])==1
            f.write(' '.join(map(str,equations(t,30)[0]))+'\n')
    print(cid,':',len(ts),'independently checked square templates')
