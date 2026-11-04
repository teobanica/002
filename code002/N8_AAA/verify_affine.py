"""Independent Fraction verification of C++ affine equivalence witnesses."""
from phase_algebra import *
def verify(targets,sources,embeddings,M):
    assert len(sources)==len(embeddings)
    for si,(s,w) in enumerate(zip(sources,embeddings)):
        assert w['source']==si
        if w['target']<0:continue
        t=targets[w['target']];H=[[[F(*x) for x in e] for e in row] for row in s['square']];K=[[[F(*x) for x in e] for e in row] for row in t['square']]
        p,q=w['rows'],w['columns'];assert sorted(p)==sorted(q)==list(range(8));d=len(H[0][0])-1
        den=w.get('parameter_denominator',1)
        X=[tuple([F(v[0],M*den)]+[F(x,den) for x in v[1:]]) for v in w['parameter_map']]
        for i in range(8):
            for j in range(8):
                lhs=[H[p[i]][q[j]][k]-H[p[0]][q[j]][k]-H[p[i]][q[0]][k]+H[p[0]][q[0]][k] for k in range(d+1)]
                rhs=[K[i][j][0]]+[F(0)]*d
                for a,x in enumerate(X):
                    for k in range(d+1):rhs[k]+=K[i][j][a+1]*x[k]
                assert (lhs[0]-rhs[0])%1==0 and lhs[1:]==rhs[1:],(si,i,j)
    return sum(x['target']>=0 for x in embeddings)
if __name__=='__main__':
    targets=read('expanded_affine_targets.json');es=read('affine_candidates.json')
    assert verify(targets,targets,es,1680)==len(targets)
    keep=sorted({x['target'] for x in es});out=[targets[i] for i in keep]
    (R/'outputs/ten_representatives.json').write_text(json.dumps(out,indent=2))
    print('Verified',len(es),'embeddings into',len(keep),'representatives:',keep)
