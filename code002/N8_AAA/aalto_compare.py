"""Parse the public Aalto catalogue, check regularity, prepare exact embeddings."""
from phase_algebra import *
import hashlib,urllib.request,concurrent.futures,argparse
parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');args=parser.parse_args()
root=R/'inputs/aalto';root.mkdir(exist_ok=True)
counts={2:1,4:15,6:36,8:143,10:299,12:756,14:1412,16:2807}
if args.download:
    def fetch(q):
        path=root/f'BH-8-{q}.txt'
        if path.exists():return
        url=f'https://wiki.aalto.fi/download/attachments/120482261/BH-8-{q}.txt?api=v2'
        data=urllib.request.urlopen(url,timeout=40).read();assert len(data.splitlines())==counts[q]
        path.write_bytes(data)
        print('Downloaded',q,len(data),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(fetch,counts))
targets=read('row_column_closed_targets.json');sources=[];meta=[]
for q,count in counts.items():
    path=root/f'BH-8-{q}.txt'
    if not path.exists():continue
    data=path.read_bytes();lines=data.splitlines();assert len(lines)==count
    rows=[];hist=Counter()
    for ix,line in enumerate(lines):
        assert len(line)==70 and line[:6]==b'BH188'+bytes([48+q])
        es=[v-48 for v in line[6:]];assert all(0<=v<q for v in es)
        H=[[(F(es[8*i+j],q),) for j in range(8)] for i in range(8)]
        assert all(H[0][j]==(0,) and H[j][0]==(0,) for j in range(8))
        ps=[profiles(quotient(H[i],H[j])) for i,j in combinations(range(8),2)]
        assert all(ps),(q,ix,ps)
        hist[tuple(sorted(Counter(ps).items()))]+=1
        pair_profiles=dict(zip(combinations(range(8),2),ps))
        aaa=sum(all('A' in pair_profiles[i,j] for i,j in combinations(t,2)) for t in combinations(range(8),3))
        rows.append({'q':q,'index':ix,'square':enc(H),'profiles':dict(Counter(ps)),'AAA_triples':aaa})
    # A common modulus per job; target signs are q/2, exponents unchanged.
    def writeH(f,t):
        H=[[[F(*x) for x in e] for e in row] for row in t['square']];d=len(H[0][0])-1
        f.write(str(d)+'\n')
        for row in H:
            assert all((e[0]*q).denominator==1 and all(x.denominator==1 for x in e[1:]) for e in row)
            f.write(' '.join(str(int(x)) for e in row for x in [e[0]*q]+list(e[1:]))+'\n')
        return H,d
    with (R/f'inputs/aalto_equivalence_{q}.txt').open('w') as f:
        f.write(str(len(targets))+'\n')
        for t in targets:
            H,d=writeH(f,t);A=S.Matrix([e[1:] for row in H for e in row]);piv=list(A.T.rref()[1]);inv=A[piv,:].inv()
            assert all(x.q==1 for x in inv)
            f.write(' '.join(map(str,piv))+'\n')
            for row in inv.tolist():f.write(' '.join(map(str,row))+'\n')
        f.write(str(len(rows))+'\n')
        for t in rows:writeH(f,t)
    (R/f'outputs/aalto_sources_{q}.json').write_text(json.dumps(rows,separators=(',',':'))+'\n')
    meta.append({'q':q,'matrices':count,'sha256':hashlib.sha256(data).hexdigest(),'row_profiles':[{'counts':dict(k),'matrices':v} for k,v in hist.items()],'AAA_triple_counts':dict(Counter(x['AAA_triples'] for x in rows))})
    print('Verified catalogue',q,count,'profiles',dict(hist),flush=True)
(R/'outputs/aalto_inventory.json').write_text(json.dumps(meta,indent=2))
