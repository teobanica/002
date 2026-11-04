from phase_algebra import *
import sys
expanded='--expanded' in sys.argv
targets=read('generic_cover_targets.json' if expanded else 'maximal_square_templates.json')
sources=[{'origin':'maximal','index':i,**t} for i,t in enumerate(targets)]
for ci in range(38):
    ds=read(f'generic_squares_{ci:02d}.json')
    sources += [{'origin':'generic','component_index':ci,**t} for t in ds['templates']]
def ints(H):
    out=[]
    for row in H:
        r=[]
        for e in row:
            v=[F(*x) for x in e];assert (v[0]*2).denominator==1 and all(x.denominator==1 for x in v[1:])
            r.append([int(v[0]*2)%2]+[int(x) for x in v[1:]])
        out.append(r)
    return out
def write_matrix(f,t):
    H=ints(t['square']);d=len(H[0][0])-1;f.write(str(d)+'\n')
    for row in H:f.write(' '.join(str(x) for e in row for x in e)+'\n')
    return H,d
with (R/('inputs/equivalence_expanded.txt' if expanded else 'inputs/equivalence.txt')).open('w') as f:
    f.write(str(len(targets))+'\n')
    for t in targets:
        H,d=write_matrix(f,t);A=S.Matrix([e[1:] for r in H for e in r]);piv=list(A.T.rref()[1]);inv=A[piv,:].inv()
        assert len(piv)==d and all(x.q==1 for x in inv)
        f.write(' '.join(map(str,piv))+'\n')
        for row in inv.tolist():f.write(' '.join(map(str,row))+'\n')
    f.write(str(len(sources))+'\n')
    for t in sources:write_matrix(f,t)
(R/('inputs/equivalence_expanded_sources.json' if expanded else 'inputs/equivalence_sources.json')).write_text(json.dumps(sources,separators=(',',':'))+'\n')
print('Prepared',len(targets),'targets and',len(sources),'sources.')
