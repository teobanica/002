"""Exact local regular families at catalogue matrices outside the current cover."""
from local_regular import *
for q in [6,8,10,12,14,16]:
    path=R/f'outputs/aalto_embeddings_{q}.json'
    if not path.exists():continue
    embeddings=json.loads(path.read_text());sources=read(f'aalto_sources_{q}.json')
    for em in embeddings:
        if em['target']>=0:continue
        i=em['source'];dest=R/f'outputs/exceptional_local_{q}_{i}.json'
        if dest.exists():continue
        print('Catalogue point',q,i,flush=True)
        H=[[F(*e[0]) for e in row] for row in sources[i]['square']]
        result=local_locus(H);result.update({'q':q,'index':i})
        dest.write_text(json.dumps(result,indent=2))
        print('Local dimensions:',[x['dimension'] for x in result['families']],flush=True)
