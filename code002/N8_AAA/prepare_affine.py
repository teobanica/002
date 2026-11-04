"""Prepare affine-family comparisons with rational inverses and torsion branches."""
from phase_algebra import *

def write_input(path,targets,sources,M):
    def writeH(f,H):
        d=len(H[0][0])-1;f.write(str(d)+'\n')
        for row in H:
            assert all((e[0]*M).denominator==1 and all(x.denominator==1 for x in e[1:]) for e in row)
            f.write(' '.join(str(int(x)) for e in row for x in [e[0]*M]+list(e[1:]))+'\n')
        return d
    with path.open('w') as f:
        f.write(str(len(targets))+'\n')
        for t in targets:
            H=[[[F(*x) for x in e] for e in row] for row in t['square']];d=writeH(f,H)
            if d:
                A=S.Matrix([e[1:] for row in H for e in row]);piv=list(A.T.rref()[1]);inv=A[piv,:].inv();den=int(S.ilcm(1,*[x.q for x in inv]))
            else:piv=[];inv=S.zeros(0,0);den=1
            assert len(piv)==d
            f.write(' '.join(map(str,piv))+'\n'+str(den)+'\n')
            for row in inv.tolist():f.write(' '.join(str(int(x*den)) for x in row)+'\n')
        f.write(str(len(sources))+'\n')
        for t in sources:writeH(f,[[[F(*x) for x in e] for e in row] for row in t['square']])

if __name__=='__main__':
    targets=read('row_column_closed_targets.json')
    for p in sorted((R/'outputs').glob('exceptional_local_*.json')):
        x=json.loads(p.read_text())
        for i,f in enumerate(x['families']):targets.append({'origin':p.name,'local_branch':i,'square':f['matrix']})
    (R/'outputs/expanded_affine_targets.json').write_text(json.dumps(targets,separators=(',',':'))+'\n')
    write_input(R/'inputs/affine_candidates.txt',targets,targets,1680)
    print('Prepared',len(targets),'targets and sources.')
