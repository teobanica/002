"""Independent verification of the complete finite Aalto comparison."""
from verify_affine import *
targets=read('ten_representatives.json');sources=json.loads((R/'inputs/aalto_all_sources.json').read_text());es=read('aalto_all_embeddings.json')
covered=verify(targets,sources,es,1680);summary=[]
for q in [2,4,6,8,10,12,14,16]:
    pairs=[(s,w) for s,w in zip(sources,es) if s['q']==q]
    summary.append({'q':q,'catalogue_entries':len(pairs),'covered':sum(w['target']>=0 for s,w in pairs),'uncovered_indices':[s['index'] for s,w in pairs if w['target']<0],'target_histogram':dict(Counter(w['target'] for s,w in pairs)),'regular':len(pairs),'AAA':sum(s['AAA_triples']>0 for s,w in pairs),'admits_A_on_every_row_pair':sum('B' not in s['profiles'] and 'C' not in s['profiles'] for s,w in pairs)})
(R/'outputs/aalto_final_summary.json').write_text(json.dumps({'checked':len(es),'covered':covered,'by_q':summary},indent=2))
print(json.dumps(summary,indent=2))
